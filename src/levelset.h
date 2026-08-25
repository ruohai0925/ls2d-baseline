// levelset.h — level-set utilities: smoothed Heaviside, material properties,
// redistancing (Sussman 1999 §3.5 with the Sussman–Fatemi volume constraint)
// and diagnostics.  Equation numbers refer to Sussman et al., JCP 148 (1999).
#pragma once
#include "grid.h"
#include "bc.h"
#include <string>

struct ReinitParams {
    // Interface half-thickness eps = alpha * dx used by H_eps and delta_eps (eq. 50, 68).
    double alpha_heaviside = 2.0;        // eps = alpha * dx for H_eps (S99 eq. 50), default 2
    // eps used for the smoothed sign function (eq. 57).  Paper: dx (alpha=1);
    // (Some codes use the Heaviside width, alpha = 2, here as well.)
    double alpha_sign = 2.0;
    // eps used by H'_eps in the volume constraint (eq. 68); paper: same as alpha_heaviside. <=0 -> same; 1 keeps the sharp area much better on thin filaments (results.md §1.2).
    double alpha_delta = -1.0;
    int    n_iter = 4;                   // number of pseudo-time steps (paper: tau = 0..eps  ->  2*alpha)
    double dtau_factor = 0.5;            // dtau = factor * dxmin  (paper: dx/2)
    std::string upwind = "godunov";      // "godunov" (eq. 62: 0 in the expansion case) | "average" ((DL+DR)/2 in the expansion case, an alternative found in other codes)
    std::string sign   = "frozen";       // "frozen" (eq. 57, from d^0) | "peng" (Peng et al. 1999: recomputed from the current d every iteration)
    bool   volume_fix = true;            // Sussman–Fatemi constraint, eq. (66)-(71)
    double vf_weight_centre = 16.0;      // 9-point quadrature weights of Sussman–Fatemi 1999 eq. (4.8): h^2/24 (16 g_ij + sum of 8 neighbours)
    double vf_weight_neigh  = 1.0;
    bool   sign_guard = true;            // guard: undo sign flips inside the band (not in the paper; never triggers in the validation cases)
};

// ---- pointwise kernels -----------------------------------------------------
double heaviside_eps(double phi, double eps);   // eq. (50)
double delta_eps(double phi, double eps);       // eq. (68)  = dH_eps/dphi
double sign_eps(double phi, double eps);        // eq. (57)  = 2 (H_eps - 1/2)

// ---- field operations ------------------------------------------------------
// rho = rho2 + (rho1 - rho2) H_eps(phi); fluid 1 is where phi > 0 (input keys ns.rho_w / ns.mu_w for phi > 0, ns.rho_a / ns.mu_a for phi < 0).
void material_from_phi(const Field2D& phi, double val_phi_pos, double val_phi_neg, double eps, Field2D& out);

// One-sided second-order ENO derivatives, eq. (58)-(59), (63)-(65).
// Returns (D^L, D^R) of `d` along x (dir=0) or y (dir=1) at cell (i,j).
void eno2_onesided(const Field2D& d, int i, int j, int dir, double h, double& DL, double& DR);

// Godunov-type upwind selection for the redistance equation, eq. (60)-(62).
double upwind_gradient(double DL, double DR, double S, const std::string& mode);

// Redistance `phi` in place:  d_tau = S(phi0)(1 - |grad d|), eq. (51)-(57),
// RK2 (Heun) in pseudo-time, eq. (54)-(55), then the volume constraint after
// every iteration, eq. (66)-(71).  Ghost cells of `phi` are refilled with `bc`.
struct ReinitStats { int sign_guard_hits = 0; double max_grad_err_band = 0; double vol_before = 0, vol_after = 0; };
ReinitStats redistance(Field2D& phi, const Geometry& g, const BCSet& bc, const ReinitParams& p);

// ---- diagnostics -----------------------------------------------------------
double volume_heaviside(const Field2D& phi, const Geometry& g, double eps);   // eq. (81), smoothed
double volume_sharp(const Field2D& phi, const Geometry& g, int nsub);          // sub-cell bilinear sampling of H(phi)
double perimeter_heaviside(const Field2D& phi, const Geometry& g, double eps); // sum delta_eps(phi) |grad phi| dA
// max and mean of | |grad phi| - 1 | over cells with |phi| < band (central differences).
void grad_norm_error(const Field2D& phi, const Geometry& g, double band, double& maxerr, double& meanerr);
