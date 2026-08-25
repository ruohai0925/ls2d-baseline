// advect.h — level-set advection with a *prescribed* cell-centred velocity
// (prescribed-velocity tests only).  Non-conservative form  phi_t + u phi_x + v phi_y = 0, i.e.
// Sussman 1999 eq. (12) with the Godunov edge-state construction replaced by
// second-order ENO one-sided derivatives (same operator as the redistance
// step, eq. 58-59) and Heun time stepping.  The unsplit Godunov predictor of
// §3.2 (`godunov.cpp`) is the scheme used for all results; this simpler one is kept for a side-by-side.
#pragma once
#include "grid.h"
#include "bc.h"
#include <functional>

// vel(x, y, t, u, v): analytic velocity at a point.
using VelocityFn = std::function<void(double, double, double, double&, double&)>;

// Right-hand side  -(u phi_x + v phi_y) at interior cells, upwind-ENO2.
void advect_rhs(const Field2D& phi, const Geometry& g, const VelocityFn& vel, double t, Field2D& rhs);

// One Heun (RK2) step from t to t+dt; refills ghosts with bc.
void advect_step_rk2(Field2D& phi, const Geometry& g, const BCSet& bc, const VelocityFn& vel, double t, double dt);

// Max |u|,|v| over cell centres at time t (for the CFL time step).
double max_velocity(const Geometry& g, const VelocityFn& vel, double t);
