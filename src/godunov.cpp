// godunov.cpp — unsplit second-order Godunov advection of a cell-centred scalar
// with the MAC velocities (Sussman 1999 §3.2, Almgren 1998 §3.2).
#include "godunov.h"
#include <cmath>

// Colella (1985) monotonicity-limited slopes as implemented in AMReX-Hydro PLM
// (hydro_godunov_plm.H): second-order MC limiter, optionally sharpened to the
// fourth-order estimate  4/3 dcen - 1/6 (df_{i+1} + df_{i-1}), re-limited.
static double mc_slope(double sm, double s0, double sp) {
    const double dlft = 2.0 * (s0 - sm), drgt = 2.0 * (sp - s0), dcen = 0.5 * (sp - sm);
    const double dsgn = (dcen >= 0) ? 1.0 : -1.0;
    const double dlim = (dlft * drgt >= 0.0) ? std::min(std::fabs(dlft), std::fabs(drgt)) : 0.0;
    return dsgn * std::min(dlim, std::fabs(dcen));
}
// Limited slope of s in direction dir (0 = x, 1 = y) at (i,j), returned in units of s (already times h).
// fourth_order = false: plain MC limiter; true: Colella's 4th-order estimate re-limited (A98).
double plm_slope(const Field2D& s, int i, int j, int dir, bool fourth_order) {
    auto v = [&](int k) { return dir == 0 ? s(i + k, j) : s(i, j + k); };
    const double df0 = mc_slope(v(-1), v(0), v(1));
    if (!fourth_order) return df0;
    const double dfm = mc_slope(v(-2), v(-1), v(0)), dfp = mc_slope(v(0), v(1), v(2));
    const double dlft = 2.0 * (v(0) - v(-1)), drgt = 2.0 * (v(1) - v(0)), dcen = 0.5 * (v(1) - v(-1));
    const double dsgn = (dcen >= 0) ? 1.0 : -1.0;
    const double dlim = (dlft * drgt >= 0.0) ? std::min(std::fabs(dlft), std::fabs(drgt)) : 0.0;
    const double dtemp = 4.0 / 3.0 * dcen - (dfp + dfm) / 6.0;
    return dsgn * std::min(dlim, std::fabs(dtemp));
}

// Time-centred face states of s (Sussman 1999 eq. 20-25 / Almgren 1998 §3.2):
//   1. normal predictor from each neighbouring cell (eq. 21, 23), optional +dt/2 src (eq. 20, 22)
//   2. "hat" states = upwind of the normal-only predictions with the face velocity
//   3. transverse correction -(dt/2) v ds/dy from the hat states, then upwind (eq. 30)
// xed(i,j) is the state on the x-face between cells i-1 and i; yed(i,j) on the y-face between j-1 and j.
void godunov_edge_states(const Field2D& s, const Geometry& g, double dt,
                         const Field2D& ucc, const Field2D& vcc, const Field2D& umac, const Field2D& vmac,
                         const GodunovParams& gp, Field2D& xed, Field2D& yed, const Field2D* src) {
    const int nx = g.nx, ny = g.ny, ng = s.ng();
    const double dx = g.dx, dy = g.dy;
    // Normal predictor only ("hat" states), eq. (21)/(23) without the transverse term,
    // stored at faces:  xhat(i,j) = state at face between cells i-1 and i.
    Field2D xlo(nx, ny, ng), xhi(nx, ny, ng), ylo(nx, ny, ng), yhi(nx, ny, ng);
    // lo = extrapolated from the cell on the low side (the "L" state), hi = from the high side ("R").
    for (int j = -1; j < ny + 1; ++j) for (int i = -1; i < nx + 2; ++i) {
        const double sx_m = plm_slope(s, i - 1, j, 0, gp.fourth_order_slopes);   // slope in cell i-1
        const double sx_p = plm_slope(s, i, j, 0, gp.fourth_order_slopes);       // slope in cell i
        xlo(i, j) = s(i - 1, j) + 0.5 * (1.0 - dt * ucc(i - 1, j) / dx) * sx_m;   // eq. (21): from cell i-1 to its right face
        xhi(i, j) = s(i, j)     - 0.5 * (1.0 + dt * ucc(i, j) / dx) * sx_p;       // eq. (23): from cell i to its left face
        if (src) { xlo(i, j) += 0.5 * dt * (*src)(i - 1, j); xhi(i, j) += 0.5 * dt * (*src)(i, j); }
    }
    for (int j = -1; j < ny + 2; ++j) for (int i = -1; i < nx + 1; ++i) {
        const double sy_m = plm_slope(s, i, j - 1, 1, gp.fourth_order_slopes);
        const double sy_p = plm_slope(s, i, j, 1, gp.fourth_order_slopes);
        ylo(i, j) = s(i, j - 1) + 0.5 * (1.0 - dt * vcc(i, j - 1) / dy) * sy_m;
        yhi(i, j) = s(i, j)     - 0.5 * (1.0 + dt * vcc(i, j) / dy) * sy_p;
        if (src) { ylo(i, j) += 0.5 * dt * (*src)(i, j - 1); yhi(i, j) += 0.5 * dt * (*src)(i, j); }
    }
    // Upwind the hat states with the face velocity to get single-valued transverse states (A98 §3.2).
    auto upw = [](double lo, double hi, double uface) { return uface > 0 ? lo : (uface < 0 ? hi : 0.5 * (lo + hi)); };
    Field2D xhat(nx, ny, ng), yhat(nx, ny, ng);
    for (int j = -1; j < ny + 1; ++j) for (int i = -1; i < nx + 2; ++i) xhat(i, j) = upw(xlo(i, j), xhi(i, j), umac(i, j));
    for (int j = -1; j < ny + 2; ++j) for (int i = -1; i < nx + 1; ++i) yhat(i, j) = upw(ylo(i, j), yhi(i, j), vmac(i, j));
    // Full left/right states: add the transverse term  -(dt/2) v ds/dy  (non-conservative form) then upwind.
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx + 1; ++i) {
        double L = xlo(i, j), R = xhi(i, j);
        if (gp.use_transverse) {
            const double vbar_m = 0.5 * (vmac(i - 1, j) + vmac(i - 1, j + 1));   // cell i-1
            const double vbar_p = 0.5 * (vmac(i, j) + vmac(i, j + 1));           // cell i
            L -= 0.5 * dt * vbar_m * (yhat(i - 1, j + 1) - yhat(i - 1, j)) / dy;
            R -= 0.5 * dt * vbar_p * (yhat(i, j + 1) - yhat(i, j)) / dy;
        }
        xed(i, j) = upw(L, R, umac(i, j));                                       // eq. (30)
    }
    for (int j = 0; j < ny + 1; ++j) for (int i = 0; i < nx; ++i) {
        double L = ylo(i, j), R = yhi(i, j);
        if (gp.use_transverse) {
            const double ubar_m = 0.5 * (umac(i, j - 1) + umac(i + 1, j - 1));
            const double ubar_p = 0.5 * (umac(i, j) + umac(i + 1, j));
            L -= 0.5 * dt * ubar_m * (xhat(i + 1, j - 1) - xhat(i, j - 1)) / dx;
            R -= 0.5 * dt * ubar_p * (xhat(i + 1, j) - xhat(i, j)) / dx;
        }
        yed(i, j) = upw(L, R, vmac(i, j));
    }
}

// Non-conservative advective derivative  [U.grad s]^{n+1/2}  at cell centres from face states
// and face velocities, Sussman 1999 eq. (32).
void godunov_advective_derivative(const Field2D& xed, const Field2D& yed, const Field2D& umac, const Field2D& vmac,
                                  const Geometry& g, Field2D& aofs) {
    for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)                 // eq. (32)
        aofs(i, j) = 0.5 * (umac(i + 1, j) + umac(i, j)) * (xed(i + 1, j) - xed(i, j)) / g.dx
                   + 0.5 * (vmac(i, j + 1) + vmac(i, j)) * (yed(i, j + 1) - yed(i, j)) / g.dy;
}

