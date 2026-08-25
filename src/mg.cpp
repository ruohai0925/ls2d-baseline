// mg.cpp — geometric multigrid V-cycles for the cell-centred (MGCell) and nodal (MGNode) operators of mg.h.
#include "mg.h"
#include <cmath>

static const int NG = 2;   // ghost layers used on every level (the stencils need one)

// Unpreconditioned CG on the coarsest level (the coarse problems are small; Jacobi alone
// would need O(n^2) sweeps on elongated domains such as 8 x 32).
template <class Op>
static void coarse_cg(const Op& A, const Field2D& bin, Field2D& x, int ni, int nj, int maxit, double rtol, bool singular) {
    Field2D b = bin, r(x), p(x), Ap(x);
    if (singular) {   // constants are in the null space: make b consistent and keep x mean-free
        double mb = 0, mx = 0; for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) { mb += b(i, j); mx += x(i, j); }
        mb /= double(ni) * nj; mx /= double(ni) * nj;
        for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) { b(i, j) -= mb; x(i, j) -= mx; }
    }
    A(x, Ap);
    double rr = 0, bb = 0;
    for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) { r(i, j) = b(i, j) - Ap(i, j); rr += r(i, j) * r(i, j); bb += b(i, j) * b(i, j); }
    if (bb == 0.0) return;
    for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) p(i, j) = r(i, j);
    for (int it = 0; it < maxit && rr > rtol * rtol * bb; ++it) {
        A(p, Ap); double pAp = 0; for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) pAp += p(i, j) * Ap(i, j);
        if (pAp <= 0) break;
        const double alpha = rr / pAp; double rr2 = 0;
        for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) { x(i, j) += alpha * p(i, j); r(i, j) -= alpha * Ap(i, j); rr2 += r(i, j) * r(i, j); }
        const double beta = rr2 / rr; rr = rr2;
        for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) p(i, j) = r(i, j) + beta * p(i, j);
    }
}

// ===========================================================================
// MGCell
// ===========================================================================
// Build the level hierarchy: copy the fine coefficients, then coarsen by 2 (face sigma: average of the two
// fine faces; cell a: average of the four fine cells) while both dimensions stay even and >= 4.
// Also precomputes the operator diagonal on every level for the Jacobi smoother.
MGCell::MGCell(int nx, int ny, double dx, double dy, const Field2D& sx, const Field2D& sy, const Field2D& a, const MGBCSet& bc,
               int nu_pre, int nu_post, double omega) : bc_(bc), nu1_(nu_pre), nu2_(nu_post), om_(omega) {
    Level l0; l0.nx = nx; l0.ny = ny; l0.dx = dx; l0.dy = dy;
    l0.sx.define(nx, ny, NG); l0.sy.define(nx, ny, NG); l0.a.define(nx, ny, NG); l0.diag.define(nx, ny, NG);
    for (int j = 0; j < ny; ++j) for (int i = 0; i <= nx; ++i) l0.sx(i, j) = sx(i, j);
    for (int j = 0; j <= ny; ++j) for (int i = 0; i < nx; ++i) l0.sy(i, j) = sy(i, j);
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) l0.a(i, j) = a(i, j);
    L_.push_back(l0); nx_.push_back(nx);
    while (L_.back().nx % 2 == 0 && L_.back().ny % 2 == 0 && L_.back().nx >= 4 && L_.back().ny >= 4) {
        const Level& f = L_.back(); Level c; c.nx = f.nx / 2; c.ny = f.ny / 2; c.dx = 2 * f.dx; c.dy = 2 * f.dy;
        c.sx.define(c.nx, c.ny, NG); c.sy.define(c.nx, c.ny, NG); c.a.define(c.nx, c.ny, NG); c.diag.define(c.nx, c.ny, NG);
        for (int j = 0; j < c.ny; ++j) for (int i = 0; i <= c.nx; ++i) c.sx(i, j) = 0.5 * (f.sx(2 * i, 2 * j) + f.sx(2 * i, 2 * j + 1));
        for (int j = 0; j <= c.ny; ++j) for (int i = 0; i < c.nx; ++i) c.sy(i, j) = 0.5 * (f.sy(2 * i, 2 * j) + f.sy(2 * i + 1, 2 * j));
        for (int j = 0; j < c.ny; ++j) for (int i = 0; i < c.nx; ++i) c.a(i, j) = 0.25 * (f.a(2 * i, 2 * j) + f.a(2 * i + 1, 2 * j) + f.a(2 * i, 2 * j + 1) + f.a(2 * i + 1, 2 * j + 1));
        L_.push_back(c); nx_.push_back(c.nx);
    }
    for (Level& l : L_) for (int j = 0; j < l.ny; ++j) for (int i = 0; i < l.nx; ++i)
        l.diag(i, j) = l.a(i, j) + (l.sx(i + 1, j) + l.sx(i, j)) / (l.dx * l.dx) + (l.sy(i, j + 1) + l.sy(i, j)) / (l.dy * l.dy);
}

// Fill one ghost layer of a cell-centred unknown on level lev: periodic wrap, even (Neumann) or odd (Dirichlet) reflection.
void MGCell::fill(int lev, Field2D& x) const {
    const int nx = L_[lev].nx, ny = L_[lev].ny;
    auto val = [](MGBC b, double in) { return b == MG_ODD ? -in : in; };
    for (int j = 0; j < ny; ++j) {
        x(-1, j) = bc_.xlo == MG_PERIODIC ? x(nx - 1, j) : val(bc_.xlo, x(0, j));
        x(nx, j) = bc_.xhi == MG_PERIODIC ? x(0, j) : val(bc_.xhi, x(nx - 1, j));
    }
    for (int i = -1; i <= nx; ++i) {
        x(i, -1) = bc_.ylo == MG_PERIODIC ? x(i, ny - 1) : val(bc_.ylo, x(i, 0));
        x(i, ny) = bc_.yhi == MG_PERIODIC ? x(i, 0) : val(bc_.yhi, x(i, ny - 1));
    }
}
// y = a x - div(sigma grad x)  (5-point) on level lev; fills ghosts of x first.
void MGCell::apply_operator(int lev, Field2D& x, Field2D& y) const {
    const Level& l = L_[lev]; fill(lev, x);
    for (int j = 0; j < l.ny; ++j) for (int i = 0; i < l.nx; ++i)
        y(i, j) = l.a(i, j) * x(i, j) - ((l.sx(i + 1, j) * (x(i + 1, j) - x(i, j)) - l.sx(i, j) * (x(i, j) - x(i - 1, j))) / (l.dx * l.dx)
                                       + (l.sy(i, j + 1) * (x(i, j + 1) - x(i, j)) - l.sy(i, j) * (x(i, j) - x(i, j - 1))) / (l.dy * l.dy));
}
// r = b - A x on level lev.
void MGCell::residual(int lev, const Field2D& b, Field2D& x, Field2D& r) const {
    apply_operator(lev, x, r); const Level& l = L_[lev];
    for (int j = 0; j < l.ny; ++j) for (int i = 0; i < l.nx; ++i) r(i, j) = b(i, j) - r(i, j);
}
// n sweeps of damped Jacobi:  x += omega * r / diag.
void MGCell::smooth(int lev, const Field2D& b, Field2D& x, int n) const {
    const Level& l = L_[lev]; Field2D r(l.nx, l.ny, NG);
    for (int s = 0; s < n; ++s) { residual(lev, b, x, r); for (int j = 0; j < l.ny; ++j) for (int i = 0; i < l.nx; ++i) x(i, j) += om_ * r(i, j) / l.diag(i, j); }
}
// Restrict a fine residual to level lev+1 by averaging the four fine cells of each coarse cell.
void MGCell::restrict_(int lev, const Field2D& rf, Field2D& rc) const {
    const Level& c = L_[lev + 1];
    for (int j = 0; j < c.ny; ++j) for (int i = 0; i < c.nx; ++i) rc(i, j) = 0.25 * (rf(2 * i, 2 * j) + rf(2 * i + 1, 2 * j) + rf(2 * i, 2 * j + 1) + rf(2 * i + 1, 2 * j + 1));
}
// Add the bilinear interpolation of the coarse correction ec (level lev+1) to xf (level lev).
void MGCell::prolong_add(int lev, const Field2D& ec, Field2D& xf) const {
    const Level& f = L_[lev]; Field2D e = ec; fill(lev + 1, e);
    for (int j = 0; j < f.ny; ++j) for (int i = 0; i < f.nx; ++i) {          // bilinear from coarse cell centres
        const int ic = i / 2, jc = j / 2; const int di = (i % 2 == 0) ? -1 : 1, dj = (j % 2 == 0) ? -1 : 1;
        xf(i, j) += 0.5625 * e(ic, jc) + 0.1875 * (e(ic + di, jc) + e(ic, jc + dj)) + 0.0625 * e(ic + di, jc + dj);
    }
}
// Recursive V-cycle: pre-smooth, restrict residual, recurse, prolong-add, post-smooth.
// The coarsest level is solved with CG (singular when all-periodic/Neumann and a == 0).
void MGCell::cycle(int lev, const Field2D& b, Field2D& x) const {
    const Level& l = L_[lev];
    if (lev == (int)L_.size() - 1) {
        auto Aop = [&](Field2D& q, Field2D& y) { apply_operator(lev, q, y); };
        const bool singular = bc_.xlo != MG_ODD && bc_.xhi != MG_ODD && bc_.ylo != MG_ODD && bc_.yhi != MG_ODD && l.a.max_abs_interior() == 0.0;
        coarse_cg(Aop, b, x, l.nx, l.ny, 500, 1e-6, singular); return;
    }
    smooth(lev, b, x, nu1_);
    Field2D r(l.nx, l.ny, NG), rc(L_[lev + 1].nx, L_[lev + 1].ny, NG), ec(L_[lev + 1].nx, L_[lev + 1].ny, NG);
    residual(lev, b, x, r); restrict_(lev, r, rc);
    cycle(lev + 1, rc, ec);
    prolong_add(lev, ec, x);
    smooth(lev, b, x, nu2_);
}
// One V-cycle from the finest level; x is both the initial guess and the result.
void MGCell::vcycle(const Field2D& b, Field2D& x) const { cycle(0, b, x); }

// ===========================================================================
// MGNode  (Q1 stiffness; nodes ordered 00,10,01,11 within a cell)
// ===========================================================================
// Element stiffness of the Q1 (bilinear) element split into x- and y-derivative parts:
//   K = sigma * (dy/dx * KX + dx/dy * KY) / 6   (the same stencil as AMReX MLNodeLaplacian).
static const double KX[4][4] = {{2,-2,1,-1},{-2,2,-1,1},{1,-1,2,-2},{-1,1,-2,2}};
static const double KY[4][4] = {{2,1,-2,-1},{1,2,-1,-2},{-2,-1,2,1},{-1,-2,1,2}};

// Build the nodal level hierarchy: cell coefficients are coarsened by 4-cell averaging; the diagonal of the
// assembled stiffness matrix is accumulated per level, and Dirichlet nodes get diag = 1 (identity rows).
MGNode::MGNode(int nx, int ny, double dx, double dy, bool per_x, bool per_y, const Field2D& sig,
               const std::function<bool(int, int, int, int)>& dirichlet, int nu_pre, int nu_post, double omega)
    : perx_(per_x), pery_(per_y), dir_(dirichlet), nu1_(nu_pre), nu2_(nu_post), om_(omega) {
    Level l0; l0.nx = nx; l0.ny = ny; l0.dx = dx; l0.dy = dy; l0.sig.define(nx, ny, NG); l0.diag.define(nx + 1, ny + 1, NG);
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) l0.sig(i, j) = sig(i, j);
    L_.push_back(l0);
    while (L_.back().nx % 2 == 0 && L_.back().ny % 2 == 0 && L_.back().nx >= 4 && L_.back().ny >= 4) {
        const Level& f = L_.back(); Level c; c.nx = f.nx / 2; c.ny = f.ny / 2; c.dx = 2 * f.dx; c.dy = 2 * f.dy;
        c.sig.define(c.nx, c.ny, NG); c.diag.define(c.nx + 1, c.ny + 1, NG);
        for (int j = 0; j < c.ny; ++j) for (int i = 0; i < c.nx; ++i) c.sig(i, j) = 0.25 * (f.sig(2 * i, 2 * j) + f.sig(2 * i + 1, 2 * j) + f.sig(2 * i, 2 * j + 1) + f.sig(2 * i + 1, 2 * j + 1));
        L_.push_back(c);
    }
    for (int lev = 0; lev < (int)L_.size(); ++lev) {
        Level& l = L_[lev]; const int ni_ = ni(lev), nj_ = nj(lev);
        for (int j = 0; j < nj_; ++j) for (int i = 0; i < ni_; ++i) l.diag(i, j) = 0;
        for (int jc = 0; jc < l.ny; ++jc) for (int ic = 0; ic < l.nx; ++ic) {
            const int ii[4] = {ic, ic + 1, ic, ic + 1}, jj[4] = {jc, jc, jc + 1, jc + 1};
            for (int a = 0; a < 4; ++a) {
                const int ia = perx_ ? ii[a] % l.nx : ii[a], ja = pery_ ? jj[a] % l.ny : jj[a];
                l.diag(ia, ja) += l.sig(ic, jc) * (l.dy / l.dx * KX[a][a] + l.dx / l.dy * KY[a][a]) / 6.0;
            }
        }
        for (int j = 0; j < nj_; ++j) for (int i = 0; i < ni_; ++i) if (dir_(i, j, l.nx, l.ny)) l.diag(i, j) = 1.0;
    }
}
// Periodic ghost images of a nodal unknown (nothing to do on non-periodic sides: natural boundary).
void MGNode::fill(int lev, Field2D& x) const {
    const Level& l = L_[lev];
    if (perx_) for (int j = -1; j <= l.ny + 1; ++j) { x(l.nx, j) = x(0, j); x(-1, j) = x(l.nx - 1, j); x(l.nx + 1, j) = x(1, j); }
    if (pery_) for (int i = -1; i <= l.nx + 1; ++i) { x(i, l.ny) = x(i, 0); x(i, -1) = x(i, l.ny - 1); x(i, l.ny + 1) = x(i, 1); }
}
// y = A x: scatter the Q1 element stiffness over every cell of level lev; Dirichlet rows become y = x.
void MGNode::apply_operator(int lev, Field2D& x, Field2D& y) const {
    const Level& l = L_[lev]; const int ni_ = ni(lev), nj_ = nj(lev); fill(lev, x);
    for (int j = 0; j < nj_; ++j) for (int i = 0; i < ni_; ++i) y(i, j) = 0;
    for (int jc = 0; jc < l.ny; ++jc) for (int ic = 0; ic < l.nx; ++ic) {
        const int ii[4] = {ic, ic + 1, ic, ic + 1}, jj[4] = {jc, jc, jc + 1, jc + 1};
        double xv[4]; for (int b = 0; b < 4; ++b) xv[b] = x(perx_ ? ii[b] % l.nx : ii[b], pery_ ? jj[b] % l.ny : jj[b]);
        for (int a = 0; a < 4; ++a) {
            double acc = 0; for (int b = 0; b < 4; ++b) acc += (l.dy / l.dx * KX[a][b] + l.dx / l.dy * KY[a][b]) / 6.0 * xv[b];
            y(perx_ ? ii[a] % l.nx : ii[a], pery_ ? jj[a] % l.ny : jj[a]) += l.sig(ic, jc) * acc;
        }
    }
    for (int j = 0; j < nj_; ++j) for (int i = 0; i < ni_; ++i) if (dir_(i, j, l.nx, l.ny)) y(i, j) = x(i, j);
}
// r = b - A x on level lev.
void MGNode::residual(int lev, const Field2D& b, Field2D& x, Field2D& r) const {
    apply_operator(lev, x, r); for (int j = 0; j < nj(lev); ++j) for (int i = 0; i < ni(lev); ++i) r(i, j) = b(i, j) - r(i, j);
}
// n sweeps of damped Jacobi on the nodal unknown.
void MGNode::smooth(int lev, const Field2D& b, Field2D& x, int n) const {
    const Level& l = L_[lev]; Field2D r(l.nx + 1, l.ny + 1, NG);
    for (int s = 0; s < n; ++s) { residual(lev, b, x, r); for (int j = 0; j < nj(lev); ++j) for (int i = 0; i < ni(lev); ++i) x(i, j) += om_ * r(i, j) / l.diag(i, j); }
}
// Full-weighting restriction of a nodal residual (transpose of the bilinear prolongation), see the note inside.
void MGNode::restrict_(int lev, const Field2D& rf, Field2D& rc) const {
    Field2D r = rf; fill(lev, r);
    const Level& f = L_[lev];
    auto val = [&](int i, int j) {   // zero outside the node range for non-periodic sides (natural boundary)
        if (!perx_ && (i < 0 || i > f.nx)) return 0.0;
        if (!pery_ && (j < 0 || j > f.ny)) return 0.0;
        return r(i, j);
    };
    for (int j = 0; j < nj(lev + 1); ++j) for (int i = 0; i < ni(lev + 1); ++i) {
        const int I = 2 * i, J = 2 * j; double s = 0;
        for (int dj = -1; dj <= 1; ++dj) for (int di = -1; di <= 1; ++di) s += (di == 0 ? 2 : 1) * (dj == 0 ? 2 : 1) * val(I + di, J + dj);
        // Transpose of the bilinear prolongation: weights [1 2 1; 2 4 2; 1 2 1]/4 (sum 4). The Q1 stiffness is
        // scale-free in h, so the coarse residual must carry the 4x larger integration area — 1/16 (the
        // cell-centred "average" restriction) would make every coarse correction 4x too small.
        rc(i, j) = s / 4.0;
        // Restricted residual at a fine Dirichlet node is meaningless; coarse Dirichlet nodes carry 0.
        if (dir_(i, j, L_[lev + 1].nx, L_[lev + 1].ny)) rc(i, j) = 0;
    }
}
// Add the bilinear interpolation of the coarse nodal correction to xf; Dirichlet nodes are left untouched.
void MGNode::prolong_add(int lev, const Field2D& ec, Field2D& xf) const {
    Field2D e = ec; fill(lev + 1, e);
    const Level& f = L_[lev], &c = L_[lev + 1];
    auto val = [&](int i, int j) {
        if (!perx_ && (i < 0 || i > c.nx)) i = std::min(std::max(i, 0), c.nx);
        if (!pery_ && (j < 0 || j > c.ny)) j = std::min(std::max(j, 0), c.ny);
        return e(i, j);
    };
    for (int j = 0; j < nj(lev); ++j) for (int i = 0; i < ni(lev); ++i) {
        const int ic = i / 2, jc = j / 2; const bool ox = i % 2, oy = j % 2;
        double v;
        if (!ox && !oy) v = val(ic, jc);
        else if (ox && !oy) v = 0.5 * (val(ic, jc) + val(ic + 1, jc));
        else if (!ox && oy) v = 0.5 * (val(ic, jc) + val(ic, jc + 1));
        else v = 0.25 * (val(ic, jc) + val(ic + 1, jc) + val(ic, jc + 1) + val(ic + 1, jc + 1));
        if (!dir_(i, j, f.nx, f.ny)) xf(i, j) += v;
    }
}
// Recursive V-cycle for the nodal operator; coarsest level solved by CG (singular if no Dirichlet node).
void MGNode::cycle(int lev, const Field2D& b, Field2D& x) const {
    const Level& l = L_[lev];
    if (lev == (int)L_.size() - 1) {
        auto Aop = [&](Field2D& q, Field2D& y) { apply_operator(lev, q, y); };
        bool has_dir = false; for (int j = 0; j < nj(lev) && !has_dir; ++j) for (int i = 0; i < ni(lev); ++i) if (dir_(i, j, l.nx, l.ny)) { has_dir = true; break; }
        coarse_cg(Aop, b, x, ni(lev), nj(lev), 500, 1e-6, !has_dir); return;
    }
    smooth(lev, b, x, nu1_);
    Field2D r(l.nx + 1, l.ny + 1, NG), rc(L_[lev + 1].nx + 1, L_[lev + 1].ny + 1, NG), ec(L_[lev + 1].nx + 1, L_[lev + 1].ny + 1, NG);
    residual(lev, b, x, r); restrict_(lev, r, rc);
    cycle(lev + 1, rc, ec);
    prolong_add(lev, ec, x);
    smooth(lev, b, x, nu2_);
}
// One V-cycle from the finest level; x is both the initial guess and the result.
void MGNode::vcycle(const Field2D& b, Field2D& x) const { cycle(0, b, x); }
