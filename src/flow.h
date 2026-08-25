// flow.h — variable-density incompressible Navier–Stokes with a level-set
// interface: the single-grid algorithm of Sussman et al. 1999 §3 (S99) /
// Almgren et al. 1998 §3 (A98)
// (NavierStokes::advance_semistaggered_twophase_ls).
//
// Unknowns and centring (S99 §3):
//   u, v, phi, rho, mu : cell centres        umac, vmac : faces (t^{n+1/2}, divergence-free)
//   p                  : nodes (p^{n-1/2})   gpx, gpy   : cell-centred lagged pressure gradient  Gp^{n-1/2}, eq. (33)
//
// One step (S99 §3.1; redistance is done right after the phi update, before the projection):
//   1. predict face velocities U^{n+1/2}, eq. (20)-(25), and MAC-project them, eq. (27)-(28)
//   2. Godunov edge states of u, v, phi, upwinded with U^ADV, eq. (29)-(30); advective derivatives eq. (31)-(32)
//   3. phi^{n+1} = phi^n - dt [U.grad phi]^{n+1/2}, eq. (12); redistance; rho, mu from phi^{n+1/2}, eq. (13)-(15)
//   4. Crank–Nicolson viscous solve for U*, eq. (16)
//   5. approximate nodal projection, eq. (17), S99 §3.4 / A98 §3.3
#pragma once
#include "grid.h"
#include "bc.h"
#include "levelset.h"
#include "godunov.h"
#include "linsolve.h"
#include "mg.h"
#include <string>

struct FlowParams {
    double rho_pos = 1.0, rho_neg = 1.0;   // density where phi > 0 (input ns.rho_w) and phi < 0 (ns.rho_a)
    double mu_pos = 0.0, mu_neg = 0.0;     // viscosity (ns.mu_w, ns.mu_a)
    double gravity = 0.0;                  // acceleration in y (ns.gravity, negative means downward)
    double sigma = 0.0;                    // surface tension coefficient (S99 §3.3, term M); 0 = off
    double cfl = 0.5;
    double fixed_dt = 0.0;
    int    do_phi = 1;                     // 0: single-phase (rho, mu constant = *_pos)
    int    do_reinit = 1, reinit_interval = 1;
    int    reinit_after_projection = 0;    // 1: paper order (redistance last), 0: redistance right after the phi update (default)
    int    init_iters = 3;                 // initial pressure iterations, S99 §3.6
    double lin_rtol = 1e-10; int lin_maxiter = 5000;
    std::string lin_precond = "mg";       // "mg" (multigrid V-cycle) | "jacobi"
    int mg_nu = 2; double mg_omega = 0.8;  // smoother sweeps (pre = post) and Jacobi damping
    int    verbose = 1;
    GodunovParams godunov;
    ReinitParams reinit;
};

// Velocity boundary codes per side (ns.lo_bc / ns.hi_bc): 0 periodic, 2 outflow, 3 symmetry, 4 slip, 5 no-slip.
struct VelBC { int xlo = 0, xhi = 0, ylo = 0, yhi = 0; };

class FlowSolver {
public:
    FlowSolver(const Geometry& g, const VelBC& vbc, const FlowParams& fp);

    // --- state (public for initialisation and diagnostics) ---
    Field2D u, v, phi, rho, mu, gpx, gpy;   // cell-centred, ng ghosts
    Field2D p;                              // node-centred, indices 0..nx, 0..ny (+ghosts)
    Field2D umac, vmac;                     // umac(i,j): x-face i-1/2 of cell i; vmac(i,j): y-face j-1/2 of cell j
    double t = 0; int step = 0;

    void set_material_from_phi();            // rho, mu from phi (S99 eq. 5-6 with H_eps)
    void initialise();                       // project U^0 and iterate for p^{1/2} (S99 §3.6)
    double estimate_dt() const;              // S99 §3.1.1
    void advance(double dt);                 // one full step

    // diagnostics
    double kinetic_energy() const;
    double max_divergence_mac() const;
    double max_divergence_nodal() const;     // eq. (33)-transpose divergence of the cell-centred velocity
    void   vorticity(Field2D& w) const;
    double curvature_max() const;             // diagnostic: max |kappa| on the band
    SolveStats last_mac, last_visc, last_nodal;

    const Geometry& geom() const { return g_; }
    void fill_velocity_ghosts(Field2D& uu, Field2D& vv) const;
    void fill_scalar(Field2D& s) const { fill_scalar_ghosts(s, sbc_); }

private:
    Geometry g_; VelBC vbc_; BCSet sbc_; FlowParams fp_; int ng_ = 4;
    Field2D phi_half_;                       // phi^{n+1/2} of the current step (eq. 13), used by the viscous RHS
    double eps_h() const { return fp_.reinit.alpha_heaviside * g_.dxmin(); }

    void viscous_operator(const Field2D& uu, const Field2D& vv, const Field2D& muf, Field2D& Lu, Field2D& Lv, bool include_transpose) const;
    void forcing_for_predictor(Field2D& su, Field2D& sv) const;                // -Gp/rho + L(U^n)/rho + F - M/rho
    void surface_tension(const Field2D& ph, Field2D& Mx, Field2D& My) const;   // M = sigma * kappa * G H^node, S99 eq. (34)-(39)
    void predict_mac_velocity(double dt, const Field2D& su, const Field2D& sv); // eq. (20)-(25) -> umac, vmac
    void mac_project(const Field2D& rho_face_src);                             // eq. (27)-(28)
    void viscous_solve(double dt, const Field2D& rho_h, const Field2D& mu_h, const Field2D& aofs_u, const Field2D& aofs_v, Field2D& us, Field2D& vs);
    void nodal_project(double dt, const Field2D& rho_h, const Field2D& us, const Field2D& vs);
    void pressure_gradient_from_nodes(const Field2D& pn, Field2D& gx, Field2D& gy) const;   // eq. (33)
    void fill_node_ghosts(Field2D& pn) const;
    int  n_nodes_x() const { return g_.periodic_x ? g_.nx : g_.nx + 1; }
    int  n_nodes_y() const { return g_.periodic_y ? g_.ny : g_.ny + 1; }
};
