// flow.h — constant-density incompressible Navier–Stokes, single grid:
// the approximate projection method of Almgren et al. 1998 §3 (A98) with the
// unsplit Godunov predictor and the equation numbering of Sussman et al. 1999 §3 (S99),
// which is the paper the routine comments cite (its §3 reduces to the present
// scheme when the two densities and viscosities are equal and sigma = 0).
//
// Solved equations (rho, mu constant, ns.rho / ns.mu):
//   U_t + (U.grad)U = -(1/rho) grad p + (1/rho) div( mu (grad U + grad U^T) ) + g e_y
//   div U = 0
//
// Unknowns and centring (S99 §3):
//   u, v : cell centres                   umac, vmac : faces (t^{n+1/2}, discretely divergence-free)
//   p    : nodes (p^{n-1/2})              gpx, gpy   : cell-centred lagged pressure gradient  Gp^{n-1/2}, eq. (33)
//
// One step (S99 §3.1):
//   1. predict face velocities U^{n+1/2}, eq. (20)-(25), and MAC-project them, eq. (27)-(28)
//   2. Godunov edge states of u, v, upwinded with U^ADV, eq. (29)-(30); advective derivatives eq. (31)-(32)
//   3. Crank–Nicolson viscous solve for U*, eq. (16)
//   4. approximate nodal projection, eq. (17), S99 §3.4 / A98 §3.3
#pragma once
#include "grid.h"
#include "bc.h"
#include "godunov.h"
#include "linsolve.h"
#include "mg.h"
#include <string>

struct FlowParams {
    double rho = 1.0;                      // constant density (ns.rho)
    double mu = 0.0;                       // constant dynamic viscosity (ns.mu); 0 = inviscid
    double gravity = 0.0;                  // acceleration in y (ns.gravity, negative means downward)
    double cfl = 0.5;
    double fixed_dt = 0.0;
    int    init_iters = 3;                 // initial pressure iterations, S99 §3.6
    double lin_rtol = 1e-10; int lin_maxiter = 5000;
    std::string lin_precond = "mg";        // "mg" (multigrid V-cycle) | "jacobi"
    int mg_nu = 2; double mg_omega = 0.8;  // smoother sweeps (pre = post) and Jacobi damping
    int    verbose = 1;
    GodunovParams godunov;
};

// Velocity boundary codes per side (ns.lo_bc / ns.hi_bc): 0 periodic, 2 outflow, 3 symmetry, 4 slip, 5 no-slip.
// vt_* is the tangential velocity of a no-slip wall (ns.wall_vel_lo / ns.wall_vel_hi): 0 for a fixed
// wall, non-zero for a moving wall such as the lid of the driven cavity.
struct VelBC {
    int xlo = 0, xhi = 0, ylo = 0, yhi = 0;
    double vt_xlo = 0, vt_xhi = 0, vt_ylo = 0, vt_yhi = 0;
};

class FlowSolver {
public:
    FlowSolver(const Geometry& g, const VelBC& vbc, const FlowParams& fp);

    // --- state (public for initialisation and diagnostics) ---
    Field2D u, v, gpx, gpy;                 // cell-centred, ng ghosts
    Field2D p;                              // node-centred, indices 0..nx, 0..ny (+ghosts)
    Field2D umac, vmac;                     // umac(i,j): x-face i-1/2 of cell i; vmac(i,j): y-face j-1/2 of cell j
    double t = 0; int step = 0;

    void initialise();                       // project U^0 and iterate for p^{1/2} (S99 §3.6)
    double estimate_dt() const;              // S99 §3.1.1
    void advance(double dt);                 // one full step

    // diagnostics
    double kinetic_energy() const;
    double max_divergence_mac() const;
    double max_divergence_nodal() const;     // eq. (33)-transpose divergence of the cell-centred velocity
    void   vorticity(Field2D& w) const;
    SolveStats last_mac, last_visc, last_nodal;

    const Geometry& geom() const { return g_; }
    // Fill the velocity ghost cells from the boundary conditions.  `homogeneous` drops the
    // wall velocities vt_*, which is what the linear operator of the viscous solve needs.
    void fill_velocity_ghosts(Field2D& uu, Field2D& vv, bool homogeneous = false) const;
    void fill_scalar(Field2D& s) const { fill_scalar_ghosts(s, sbc_); }

private:
    Geometry g_; VelBC vbc_; BCSet sbc_; FlowParams fp_; int ng_ = 4;

    void viscous_operator(const Field2D& uu, const Field2D& vv, Field2D& Lu, Field2D& Lv, bool include_transpose) const;
    void forcing_for_predictor(Field2D& su, Field2D& sv) const;                 // -Gp/rho + L(U^n)/rho + F
    void predict_mac_velocity(double dt, const Field2D& su, const Field2D& sv); // eq. (20)-(25) -> umac, vmac
    void mac_project();                                                         // eq. (27)-(28)
    void viscous_solve(double dt, const Field2D& aofs_u, const Field2D& aofs_v, Field2D& us, Field2D& vs);
    void nodal_project(double dt, const Field2D& us, const Field2D& vs);
    void pressure_gradient_from_nodes(const Field2D& pn, Field2D& gx, Field2D& gy) const;   // eq. (33)
    void fill_node_ghosts(Field2D& pn) const;
    int  n_nodes_x() const { return g_.periodic_x ? g_.nx : g_.nx + 1; }
    int  n_nodes_y() const { return g_.periodic_y ? g_.ny : g_.ny + 1; }
};
