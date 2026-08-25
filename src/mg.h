// mg.h — geometric multigrid V-cycles used as PCG preconditioners.
//
// Two operator families appear in the projection method:
//   MGCell : cell-centred  a x - div(sigma grad x)   (5-point; MAC projection with a=0, viscous solve with a=rho/dt)
//   MGNode : node-centred  sum_cells sigma_c int grad N_a . grad N_b   (Q1 finite-element nodal Laplacian)
// Smoother: damped Jacobi (symmetric, so the V-cycle is a valid SPD preconditioner);
// restriction: 4-cell average (cells) / 9-point full weighting (nodes); prolongation: bilinear.
// Coarse coefficients are re-discretised (face sigma averaged, cell sigma averaged), not Galerkin.
// The number of levels is limited by the largest power of two dividing nx and ny.
#pragma once
#include "grid.h"
#include <vector>
#include <functional>

// Boundary treatment of the unknown at each side: periodic, even (homogeneous Neumann), odd (homogeneous Dirichlet on the face).
enum MGBC { MG_PERIODIC = 0, MG_EVEN = 1, MG_ODD = 2 };
struct MGBCSet { MGBC xlo = MG_PERIODIC, xhi = MG_PERIODIC, ylo = MG_PERIODIC, yhi = MG_PERIODIC; };

class MGCell {
public:
    // sx(i,j): sigma on x-face i-1/2 (i = 0..nx), sy(i,j): sigma on y-face j-1/2 (j = 0..ny); a(i,j) >= 0 cell coefficient.
    MGCell(int nx, int ny, double dx, double dy, const Field2D& sx, const Field2D& sy, const Field2D& a, const MGBCSet& bc,
           int nu_pre = 2, int nu_post = 2, double omega = 0.8);
    void apply_operator(int lev, Field2D& x, Field2D& y) const;      // y = A x (fills ghosts of x)
    void vcycle(const Field2D& b, Field2D& x) const;                  // one V-cycle from x (interior), result in x
    int levels() const { return (int)nx_.size(); }
    void fill(int lev, Field2D& x) const;
private:
    struct Level { int nx, ny; double dx, dy; Field2D sx, sy, a, diag; };
    std::vector<int> nx_; std::vector<Level> L_; MGBCSet bc_; int nu1_, nu2_; double om_;
    void smooth(int lev, const Field2D& b, Field2D& x, int n) const;
    void residual(int lev, const Field2D& b, Field2D& x, Field2D& r) const;
    void restrict_(int lev, const Field2D& rf, Field2D& rc) const;
    void prolong_add(int lev, const Field2D& ec, Field2D& xf) const;
    void cycle(int lev, const Field2D& b, Field2D& x) const;
};

class MGNode {
public:
    // sig(i,j): cell coefficient (i = 0..nx-1). Nodes are 0..nx (or 0..nx-1 if periodic). dirichlet(i,j) marks p=0 nodes (outflow).
    MGNode(int nx, int ny, double dx, double dy, bool per_x, bool per_y, const Field2D& sig,
           const std::function<bool(int, int, int, int)>& dirichlet,   // (i, j, nx_lev, ny_lev) -> is Dirichlet node at this level
           int nu_pre = 2, int nu_post = 2, double omega = 0.8);
    void apply_operator(int lev, Field2D& x, Field2D& y) const;
    void vcycle(const Field2D& b, Field2D& x) const;
    int levels() const { return (int)L_.size(); }
    void fill(int lev, Field2D& x) const;
    int ni(int lev) const { return perx_ ? L_[lev].nx : L_[lev].nx + 1; }
    int nj(int lev) const { return pery_ ? L_[lev].ny : L_[lev].ny + 1; }
private:
    struct Level { int nx, ny; double dx, dy; Field2D sig, diag; };
    std::vector<Level> L_; bool perx_, pery_; std::function<bool(int, int, int, int)> dir_; int nu1_, nu2_; double om_;
    void smooth(int lev, const Field2D& b, Field2D& x, int n) const;
    void residual(int lev, const Field2D& b, Field2D& x, Field2D& r) const;
    void restrict_(int lev, const Field2D& rf, Field2D& rc) const;
    void prolong_add(int lev, const Field2D& ec, Field2D& xf) const;
    void cycle(int lev, const Field2D& b, Field2D& x) const;
};
