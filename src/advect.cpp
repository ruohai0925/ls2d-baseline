// advect.cpp — prescribed-velocity level-set advection with upwind ENO2 + Heun (simple reference scheme, see advect.h).
#include "advect.h"
#include "levelset.h"
#include <cmath>

// rhs = -(u phi_x + v phi_y) at interior cells; phi_x, phi_y are the one-sided ENO2
// derivatives of Sussman 1999 eq. (58)-(59), chosen on the upwind side of the velocity.
void advect_rhs(const Field2D& phi, const Geometry& g, const VelocityFn& vel, double t, Field2D& rhs) {
    for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
        double u, v; vel(g.xc(i), g.yc(j), t, u, v);
        double DLx, DRx, DLy, DRy;
        eno2_onesided(phi, i, j, 0, g.dx, DLx, DRx);
        eno2_onesided(phi, i, j, 1, g.dy, DLy, DRy);
        const double px = (u > 0) ? DLx : DRx;      // upwind: information comes from the side the flow comes from
        const double py = (v > 0) ? DLy : DRy;
        rhs(i, j) = -(u * px + v * py);
    }
}

// Heun / RK2 update of phi from t to t+dt using advect_rhs; ghosts are refilled before each RHS evaluation.
void advect_step_rk2(Field2D& phi, const Geometry& g, const BCSet& bc, const VelocityFn& vel, double t, double dt) {
    Field2D k1(g.nx, g.ny, phi.ng()), p1(g.nx, g.ny, phi.ng()), k2(g.nx, g.ny, phi.ng());
    fill_scalar_ghosts(phi, bc);
    advect_rhs(phi, g, vel, t, k1);
    for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) p1(i, j) = phi(i, j) + dt * k1(i, j);
    fill_scalar_ghosts(p1, bc);
    advect_rhs(p1, g, vel, t + dt, k2);
    for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) phi(i, j) += 0.5 * dt * (k1(i, j) + k2(i, j));
    fill_scalar_ghosts(phi, bc);
}

// Largest |u| or |v| over cell centres at time t (CFL estimate).
double max_velocity(const Geometry& g, const VelocityFn& vel, double t) {
    double m = 0;
    for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) { double u, v; vel(g.xc(i), g.yc(j), t, u, v); m = std::max(m, std::max(std::fabs(u), std::fabs(v))); }
    return m;
}
