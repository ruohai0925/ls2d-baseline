// godunov.h — unsplit second-order Godunov (predictor–corrector) advection of a
// cell-centred scalar, Sussman 1999 §3.2 / Almgren 1998 §3.2 (Colella 1990).
//
// The cell-centred velocity at t^n feeds the normal predictor (eq. 21, 23) and the
// MAC-projected face velocity at t^{n+1/2} plays the role of U^ADV (eq. 28).
#pragma once
#include "grid.h"
#include "bc.h"

struct GodunovParams {
    bool fourth_order_slopes = true;   // Colella's 4th-order limited slope (A98); false -> 2nd-order MC-limited
    bool use_transverse = true;        // include the transverse derivative terms of eq. (21)/(23)
};

// Monotonicity-limited slope of s in direction dir at (i,j), in units of s (already multiplied by h).
double plm_slope(const Field2D& s, int i, int j, int dir, bool fourth_order);

// Face states s_{i+1/2,j} (stored at x-face (i+1,j), i.e. xed(i,j) is the face between cells i-1 and i)
// and s_{i,j+1/2} (yed(i,j) between cells j-1 and j), upwinded with the face velocities.
void godunov_edge_states(const Field2D& s, const Geometry& g, double dt,
                         const Field2D& ucc, const Field2D& vcc,          // cell-centred velocity at t^n
                         const Field2D& umac, const Field2D& vmac,        // face velocity at t^{n+1/2}
                         const GodunovParams& gp, Field2D& xed, Field2D& yed,
                         const Field2D* src = nullptr);                        // optional +dt/2 * src in the predictor (eq. 20, 22)

// Non-conservative advective derivative, Sussman 1999 eq. (32).
void godunov_advective_derivative(const Field2D& xed, const Field2D& yed, const Field2D& umac, const Field2D& vmac,
                                  const Geometry& g, Field2D& aofs);
