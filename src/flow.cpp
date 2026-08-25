// flow.cpp — the two-phase projection step of Sussman 1999 §3 / Almgren 1998 §3 (algorithm outline in flow.h).
#include "flow.h"
#include <cmath>
#include <cstdio>

static const double PI_ = 3.141592653589793238462643383279502884197;

// Allocate all state fields; scalar ghost BCs are periodic where the velocity is periodic, else zero-gradient.
FlowSolver::FlowSolver(const Geometry& g, const VelBC& vbc, const FlowParams& fp) : g_(g), vbc_(vbc), fp_(fp) {
    sbc_.xlo = (vbc.xlo == BC_PERIODIC) ? BC_PERIODIC : BC_OUTFLOW; sbc_.xhi = (vbc.xhi == BC_PERIODIC) ? BC_PERIODIC : BC_OUTFLOW;
    sbc_.ylo = (vbc.ylo == BC_PERIODIC) ? BC_PERIODIC : BC_OUTFLOW; sbc_.yhi = (vbc.yhi == BC_PERIODIC) ? BC_PERIODIC : BC_OUTFLOW;
    const int nx = g.nx, ny = g.ny;
    for (Field2D* f : {&u, &v, &phi, &rho, &mu, &gpx, &gpy, &umac, &vmac}) f->define(nx, ny, ng_);
    p.define(nx + 1, ny + 1, ng_);
    rho.fill(fp.rho_pos); mu.fill(fp.mu_pos);
}

// ---------------------------------------------------------------------------
// ghost cells
// ---------------------------------------------------------------------------
// Velocity: periodic wrap; walls reflect (normal component odd; tangential odd
// for no-slip, even for slip/symmetry); outflow zero-gradient.
void FlowSolver::fill_velocity_ghosts(Field2D& uu, Field2D& vv) const {
    const int nx = g_.nx, ny = g_.ny, ng = ng_;
    auto side = [&](int code, Field2D& f, bool is_normal, int i_in, int i_gh, int j, bool xdir) {
        double val;
        const double in = xdir ? f(i_in, j) : f(j, i_in);
        switch (code) {
            case BC_PERIODIC: { const int n = xdir ? nx : ny; const int src = (i_gh < 0) ? i_gh + n : i_gh - n; val = xdir ? f(src, j) : f(j, src); break; }
            case BC_OUTFLOW:  val = in; break;
            case BC_NOSLIPWALL: val = -in; break;
            case BC_SLIPWALL: case BC_SYMMETRY: val = is_normal ? -in : in; break;
            default: val = in;
        }
        if (xdir) f(i_gh, j) = val; else f(j, i_gh) = val;
    };
    for (Field2D* f : {&uu, &vv}) {
        const bool isu = (f == &uu);
        for (int j = 0; j < ny; ++j) for (int g = 1; g <= ng; ++g) {
            side(vbc_.xlo, *f, isu, g - 1, -g, j, true);                 // reflect about the face x = xlo: ghost -g <-> interior g-1
            side(vbc_.xhi, *f, isu, nx - g, nx - 1 + g, j, true);
        }
        for (int i = -ng; i < nx + ng; ++i) for (int g = 1; g <= ng; ++g) {
            side(vbc_.ylo, *f, !isu, g - 1, -g, i, false);
            side(vbc_.yhi, *f, !isu, ny - g, ny - 1 + g, i, false);
        }
    }
}

// Periodic images of a node-centred field (node nx == node 0); nothing is done on non-periodic sides.
void FlowSolver::fill_node_ghosts(Field2D& pn) const {
    const int nx = g_.nx, ny = g_.ny, ng = ng_;
    if (g_.periodic_x) { for (int j = -ng; j < ny + 1 + ng; ++j) for (int g = 0; g <= ng; ++g) { pn(nx + g, j) = pn(g, j); if (g > 0) pn(-g, j) = pn(nx - g, j); } }
    if (g_.periodic_y) { for (int i = -ng; i < nx + 1 + ng; ++i) for (int g = 0; g <= ng; ++g) { pn(i, ny + g) = pn(i, g); if (g > 0) pn(i, -g) = pn(i, ny - g); } }
}

// rho, mu = smoothed-Heaviside blend of the two phases from the current phi (S99 eq. 5-6); constant if do_phi = 0.
void FlowSolver::set_material_from_phi() {
    if (fp_.do_phi) { material_from_phi(phi, fp_.rho_pos, fp_.rho_neg, eps_h(), rho); material_from_phi(phi, fp_.mu_pos, fp_.mu_neg, eps_h(), mu); }
    else { rho.fill(fp_.rho_pos); mu.fill(fp_.mu_pos); }
}

// ---------------------------------------------------------------------------
// viscous term  L = div( mu (grad U + grad U^T) ),  S99 §3.3 (first component written out there)
// include_transpose=false gives only div(mu grad U), the part treated implicitly.
// ---------------------------------------------------------------------------
void FlowSolver::viscous_operator(const Field2D& uu, const Field2D& vv, const Field2D& muf, Field2D& Lu, Field2D& Lv, bool include_transpose) const {
    const double dx = g_.dx, dy = g_.dy;
    auto mux = [&](int i, int j) { return 0.5 * (muf(i - 1, j) + muf(i, j)); };   // mu at x-face i-1/2
    auto muy = [&](int i, int j) { return 0.5 * (muf(i, j - 1) + muf(i, j)); };   // mu at y-face j-1/2
    for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i) {
        // div(mu grad u)
        double lu = (mux(i + 1, j) * (uu(i + 1, j) - uu(i, j)) - mux(i, j) * (uu(i, j) - uu(i - 1, j))) / (dx * dx)
                  + (muy(i, j + 1) * (uu(i, j + 1) - uu(i, j)) - muy(i, j) * (uu(i, j) - uu(i, j - 1))) / (dy * dy);
        double lv = (mux(i + 1, j) * (vv(i + 1, j) - vv(i, j)) - mux(i, j) * (vv(i, j) - vv(i - 1, j))) / (dx * dx)
                  + (muy(i, j + 1) * (vv(i, j + 1) - vv(i, j)) - muy(i, j) * (vv(i, j) - vv(i, j - 1))) / (dy * dy);
        if (include_transpose) {
            // div(mu grad U^T): x-component = d/dx(mu u_x) + d/dy(mu v_x); y-component = d/dx(mu u_y) + d/dy(mu v_y)
            lu += (mux(i + 1, j) * (uu(i + 1, j) - uu(i, j)) - mux(i, j) * (uu(i, j) - uu(i - 1, j))) / (dx * dx)
                + (muy(i, j + 1) * (vv(i + 1, j + 1) - vv(i - 1, j + 1) + vv(i + 1, j) - vv(i - 1, j))
                 - muy(i, j)     * (vv(i + 1, j) - vv(i - 1, j) + vv(i + 1, j - 1) - vv(i - 1, j - 1))) / (4 * dx * dy);
            lv += (muy(i, j + 1) * (vv(i, j + 1) - vv(i, j)) - muy(i, j) * (vv(i, j) - vv(i, j - 1))) / (dy * dy)
                + (mux(i + 1, j) * (uu(i + 1, j + 1) - uu(i + 1, j - 1) + uu(i, j + 1) - uu(i, j - 1))
                 - mux(i, j)     * (uu(i, j + 1) - uu(i, j - 1) + uu(i - 1, j + 1) - uu(i - 1, j - 1))) / (4 * dx * dy);
        }
        Lu(i, j) = lu; Lv(i, j) = lv;
    }
}

// ---------------------------------------------------------------------------
// Surface tension, S99 §3.3 eq. (34)-(39):  M = sigma * kappa * G(H^node)
//   N = G phi / |G phi| at nodes (eq. 35-36, cell -> node gradient), kappa = D N at cells (eq. 37-38, node -> cell
//   divergence), H^node = average of the four cell Heavisides (eq. 39), G H^node at cells (same G as eq. 33).
// The momentum equation carries -M/rho (eq. 16), i.e. the continuum-surface-force  -sigma kappa grad H  with
// kappa = div(n), n = grad phi/|grad phi| pointing into the phi > 0 fluid.
// ---------------------------------------------------------------------------
void FlowSolver::surface_tension(const Field2D& ph, Field2D& Mx, Field2D& My) const {
    const int nx = g_.nx, ny = g_.ny, ng = ng_; const double dx = g_.dx, dy = g_.dy, eps = eps_h();
    Mx.fill(0); My.fill(0);
    if (fp_.sigma == 0.0) return;
    Field2D Nx(nx + 1, ny + 1, ng), Ny(nx + 1, ny + 1, ng), Hn(nx + 1, ny + 1, ng);   // node-centred
    for (int j = -ng + 1; j < ny + ng; ++j) for (int i = -ng + 1; i < nx + ng; ++i) {   // node (i,j) = corner shared by cells i-1,i / j-1,j
        const double gx = (ph(i, j) + ph(i, j - 1) - ph(i - 1, j) - ph(i - 1, j - 1)) / (2 * dx);   // eq. (36)
        const double gy = (ph(i, j) + ph(i - 1, j) - ph(i, j - 1) - ph(i - 1, j - 1)) / (2 * dy);
        const double m = std::sqrt(gx * gx + gy * gy);
        Nx(i, j) = m > 1e-12 ? gx / m : 0.0; Ny(i, j) = m > 1e-12 ? gy / m : 0.0;                    // eq. (35)
        Hn(i, j) = 0.25 * (heaviside_eps(ph(i, j), eps) + heaviside_eps(ph(i - 1, j), eps)
                         + heaviside_eps(ph(i, j - 1), eps) + heaviside_eps(ph(i - 1, j - 1), eps)); // eq. (39)
    }
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        // nodes of cell (i,j): (i,j) (i+1,j) (i,j+1) (i+1,j+1)
        const double kappa = (Nx(i + 1, j) + Nx(i + 1, j + 1) - Nx(i, j) - Nx(i, j + 1)) / (2 * dx)          // eq. (37)-(38)
                           + (Ny(i, j + 1) + Ny(i + 1, j + 1) - Ny(i, j) - Ny(i + 1, j)) / (2 * dy);
        const double gHx = (Hn(i + 1, j) + Hn(i + 1, j + 1) - Hn(i, j) - Hn(i, j + 1)) / (2 * dx);          // G H^node, eq. (33) form
        const double gHy = (Hn(i, j + 1) + Hn(i + 1, j + 1) - Hn(i, j) - Hn(i + 1, j)) / (2 * dy);
        Mx(i, j) = fp_.sigma * kappa * gHx; My(i, j) = fp_.sigma * kappa * gHy;                               // eq. (34)
    }
    fill_scalar_ghosts(Mx, sbc_); fill_scalar_ghosts(My, sbc_);
}

// Diagnostic: max |kappa grad H| over the domain (0 without surface tension).
double FlowSolver::curvature_max() const {
    if (fp_.sigma == 0.0) return 0.0;
    Field2D Mx(g_.nx, g_.ny, ng_), My(g_.nx, g_.ny, ng_); surface_tension(phi, Mx, My);
    double m = 0; for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i) m = std::max(m, std::hypot(Mx(i, j), My(i, j)));
    return m / fp_.sigma;   // |kappa grad H| max
}

// Source term of the normal predictor, eq. (20)/(22):  -Gp^{n-1/2}/rho^n + L^n/rho^n - M^n/rho^n + F
void FlowSolver::forcing_for_predictor(Field2D& su, Field2D& sv) const {
    Field2D Lu(g_.nx, g_.ny, ng_), Lv(g_.nx, g_.ny, ng_), Mx(g_.nx, g_.ny, ng_), My(g_.nx, g_.ny, ng_);
    viscous_operator(u, v, mu, Lu, Lv, true);
    surface_tension(phi, Mx, My);
    for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i) {
        su(i, j) = (-gpx(i, j) + Lu(i, j) - Mx(i, j)) / rho(i, j);
        sv(i, j) = (-gpy(i, j) + Lv(i, j) - My(i, j)) / rho(i, j) + fp_.gravity;
    }
    fill_scalar_ghosts(su, sbc_); fill_scalar_ghosts(sv, sbc_);
}

// ---------------------------------------------------------------------------
// Face-velocity predictor, S99 eq. (20)-(25): extrapolate u to x-faces and v to
// y-faces (normal PLM predictor + transverse correction with the "hat" states),
// then choose the face value with the upwind rule (25) applied to the face's own
// left/right states.  Mirrors AMReX-Hydro Godunov::ExtrapVelToFaces.
// ---------------------------------------------------------------------------
void FlowSolver::predict_mac_velocity(double dt, const Field2D& su, const Field2D& sv) {
    const int nx = g_.nx, ny = g_.ny, ng = ng_; const double dx = g_.dx, dy = g_.dy;
    const bool fo = fp_.godunov.fourth_order_slopes;
    Field2D ulo(nx, ny, ng), uhi(nx, ny, ng), vlo(nx, ny, ng), vhi(nx, ny, ng);   // u at x-faces, v at y-faces (normal-only)
    Field2D uylo(nx, ny, ng), uyhi(nx, ny, ng), vxlo(nx, ny, ng), vxhi(nx, ny, ng); // u at y-faces, v at x-faces (transverse hats)
    for (int j = -1; j < ny + 2; ++j) for (int i = -1; i < nx + 2; ++i) {
        // u to x-face i-1/2
        ulo(i, j) = u(i - 1, j) + 0.5 * (1.0 - dt * u(i - 1, j) / dx) * plm_slope(u, i - 1, j, 0, fo) + 0.5 * dt * su(i - 1, j);
        uhi(i, j) = u(i, j)     - 0.5 * (1.0 + dt * u(i, j) / dx)     * plm_slope(u, i, j, 0, fo)     + 0.5 * dt * su(i, j);
        // v to y-face j-1/2
        vlo(i, j) = v(i, j - 1) + 0.5 * (1.0 - dt * v(i, j - 1) / dy) * plm_slope(v, i, j - 1, 1, fo) + 0.5 * dt * sv(i, j - 1);
        vhi(i, j) = v(i, j)     - 0.5 * (1.0 + dt * v(i, j) / dy)     * plm_slope(v, i, j, 1, fo)     + 0.5 * dt * sv(i, j);
        // transverse hats: u to y-faces, v to x-faces (normal predictor only)
        uylo(i, j) = u(i, j - 1) + 0.5 * (1.0 - dt * v(i, j - 1) / dy) * plm_slope(u, i, j - 1, 1, fo);
        uyhi(i, j) = u(i, j)     - 0.5 * (1.0 + dt * v(i, j) / dy)     * plm_slope(u, i, j, 1, fo);
        vxlo(i, j) = v(i - 1, j) + 0.5 * (1.0 - dt * u(i - 1, j) / dx) * plm_slope(v, i - 1, j, 0, fo);
        vxhi(i, j) = v(i, j)     - 0.5 * (1.0 + dt * u(i, j) / dx)     * plm_slope(v, i, j, 0, fo);
    }
    // Riemann-type selection, eq. (25), for a face whose normal velocity is the unknown itself.
    auto riemann = [](double L, double R) { if (L > 0 && L + R > 0) return L; if (R < 0 && L + R < 0) return R; return 0.0; };
    auto upw = [](double L, double R, double w) { return w > 0 ? L : (w < 0 ? R : 0.5 * (L + R)); };
    Field2D uhat_x(nx, ny, ng), vhat_y(nx, ny, ng), uhat_y(nx, ny, ng), vhat_x(nx, ny, ng);
    for (int j = -1; j < ny + 2; ++j) for (int i = -1; i < nx + 2; ++i) {
        uhat_x(i, j) = riemann(ulo(i, j), uhi(i, j));          // normal velocity at x-faces (no transverse yet)
        vhat_y(i, j) = riemann(vlo(i, j), vhi(i, j));          // normal velocity at y-faces
        uhat_y(i, j) = upw(uylo(i, j), uyhi(i, j), vhat_y(i, j));   // u transported across y-faces
        vhat_x(i, j) = upw(vxlo(i, j), vxhi(i, j), uhat_x(i, j));   // v transported across x-faces
    }
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx + 1; ++i) {
        double L = ulo(i, j), R = uhi(i, j);
        if (fp_.godunov.use_transverse) {
            L -= 0.5 * dt * 0.5 * (vhat_y(i - 1, j) + vhat_y(i - 1, j + 1)) * (uhat_y(i - 1, j + 1) - uhat_y(i - 1, j)) / dy;
            R -= 0.5 * dt * 0.5 * (vhat_y(i, j) + vhat_y(i, j + 1))         * (uhat_y(i, j + 1) - uhat_y(i, j)) / dy;
        }
        umac(i, j) = riemann(L, R);
    }
    for (int j = 0; j < ny + 1; ++j) for (int i = 0; i < nx; ++i) {
        double L = vlo(i, j), R = vhi(i, j);
        if (fp_.godunov.use_transverse) {
            L -= 0.5 * dt * 0.5 * (uhat_x(i, j - 1) + uhat_x(i + 1, j - 1)) * (vhat_x(i + 1, j - 1) - vhat_x(i, j - 1)) / dx;
            R -= 0.5 * dt * 0.5 * (uhat_x(i, j) + uhat_x(i + 1, j))         * (vhat_x(i + 1, j) - vhat_x(i, j)) / dx;
        }
        vmac(i, j) = riemann(L, R);
    }
    // Physical boundaries: no normal flow through walls; periodic wrap is handled by the ghost fill of the callers.
    if (vbc_.xlo != BC_PERIODIC && vbc_.xlo != BC_OUTFLOW) for (int j = 0; j < ny; ++j) umac(0, j) = 0.0;
    if (vbc_.xhi != BC_PERIODIC && vbc_.xhi != BC_OUTFLOW) for (int j = 0; j < ny; ++j) umac(nx, j) = 0.0;
    if (vbc_.ylo != BC_PERIODIC && vbc_.ylo != BC_OUTFLOW) for (int i = 0; i < nx; ++i) vmac(i, 0) = 0.0;
    if (vbc_.yhi != BC_PERIODIC && vbc_.yhi != BC_OUTFLOW) for (int i = 0; i < nx; ++i) vmac(i, ny) = 0.0;
    if (g_.periodic_x) for (int j = 0; j < ny; ++j) umac(nx, j) = umac(0, j);
    if (g_.periodic_y) for (int i = 0; i < nx; ++i) vmac(i, ny) = vmac(i, 0);
}

// ---------------------------------------------------------------------------
// MAC projection, eq. (27)-(28): cell-centred p_mac, sigma = 1/rho on faces.
// ---------------------------------------------------------------------------
void FlowSolver::mac_project(const Field2D& rho_src) {
    const int nx = g_.nx, ny = g_.ny, ng = ng_; const double dx = g_.dx, dy = g_.dy;
    Field2D sx(nx, ny, ng), sy(nx, ny, ng);                       // sigma on x-faces / y-faces
    for (int j = -ng + 1; j < ny + ng; ++j) for (int i = -ng + 1; i < nx + ng; ++i) {
        sx(i, j) = 2.0 / (rho_src(i - 1, j) + rho_src(i, j));
        sy(i, j) = 2.0 / (rho_src(i, j - 1) + rho_src(i, j));
    }
    Field2D rhs(nx, ny, ng), pm(nx, ny, ng), diag(nx, ny, ng);
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        rhs(i, j) = (umac(i + 1, j) - umac(i, j)) / dx + (vmac(i, j + 1) - vmac(i, j)) / dy;
        diag(i, j) = (sx(i + 1, j) + sx(i, j)) / (dx * dx) + (sy(i, j + 1) + sy(i, j)) / (dy * dy);
    }
    // Wall faces carry no flux: drop them from the operator (homogeneous Neumann on p_mac). Outflow: p_mac = 0 (Dirichlet).
    const bool wx_lo = vbc_.xlo != BC_PERIODIC && vbc_.xlo != BC_OUTFLOW, wx_hi = vbc_.xhi != BC_PERIODIC && vbc_.xhi != BC_OUTFLOW;
    const bool wy_lo = vbc_.ylo != BC_PERIODIC && vbc_.ylo != BC_OUTFLOW, wy_hi = vbc_.yhi != BC_PERIODIC && vbc_.yhi != BC_OUTFLOW;
    if (wx_lo) for (int j = 0; j < ny; ++j) sx(0, j) = 0.0;
    if (wx_hi) for (int j = 0; j < ny; ++j) sx(nx, j) = 0.0;
    if (wy_lo) for (int i = 0; i < nx; ++i) sy(i, 0) = 0.0;
    if (wy_hi) for (int i = 0; i < nx; ++i) sy(i, ny) = 0.0;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i)
        diag(i, j) = (sx(i + 1, j) + sx(i, j)) / (dx * dx) + (sy(i, j + 1) + sy(i, j)) / (dy * dy);
    auto fill_p = [&](Field2D& q) {  // periodic wrap or zero-gradient (Neumann); outflow -> odd reflection gives p=0 on the boundary face
        BCSet b = sbc_; fill_scalar_ghosts(q, b);
        auto odd = [&](int code) { return code == BC_OUTFLOW; };
        if (odd(vbc_.xlo)) for (int j = -ng; j < ny + ng; ++j) q(-1, j) = -q(0, j);
        if (odd(vbc_.xhi)) for (int j = -ng; j < ny + ng; ++j) q(nx, j) = -q(nx - 1, j);
        if (odd(vbc_.ylo)) for (int i = -ng; i < nx + ng; ++i) q(i, -1) = -q(i, 0);
        if (odd(vbc_.yhi)) for (int i = -ng; i < nx + ng; ++i) q(i, ny) = -q(i, ny - 1);
    };
    auto A = [&](Field2D& x, Field2D& y) {
        fill_p(x);
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i)
            y(i, j) = -((sx(i + 1, j) * (x(i + 1, j) - x(i, j)) - sx(i, j) * (x(i, j) - x(i - 1, j))) / (dx * dx)
                      + (sy(i, j + 1) * (x(i, j + 1) - x(i, j)) - sy(i, j) * (x(i, j) - x(i, j - 1))) / (dy * dy));
    };
    const bool has_outflow = vbc_.xlo == BC_OUTFLOW || vbc_.xhi == BC_OUTFLOW || vbc_.ylo == BC_OUTFLOW || vbc_.yhi == BC_OUTFLOW;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) rhs(i, j) = -rhs(i, j);     // A p = -div U  (A is -div sigma grad)
    if (fp_.lin_precond == "mg") {
        auto mgbc = [&](int code) { return code == BC_PERIODIC ? MG_PERIODIC : (code == BC_OUTFLOW ? MG_ODD : MG_EVEN); };
        MGBCSet mb; mb.xlo = mgbc(vbc_.xlo); mb.xhi = mgbc(vbc_.xhi); mb.ylo = mgbc(vbc_.ylo); mb.yhi = mgbc(vbc_.yhi);
        Field2D zero(nx, ny, ng); MGCell mg(nx, ny, dx, dy, sx, sy, zero, mb, fp_.mg_nu, fp_.mg_nu, fp_.mg_omega);
        Precond M = [&](const Field2D& r, Field2D& z) { mg.vcycle(r, z); };
        last_mac = pcg(A, M, rhs, pm, nx, ny, fp_.lin_rtol, fp_.lin_maxiter, !has_outflow);
    } else
        last_mac = pcg(A, diag, rhs, pm, nx, ny, fp_.lin_rtol, fp_.lin_maxiter, !has_outflow);
    fill_p(pm);
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx + 1; ++i) umac(i, j) -= sx(i, j) * (pm(i, j) - pm(i - 1, j)) / dx;   // eq. (28)
    for (int j = 0; j < ny + 1; ++j) for (int i = 0; i < nx; ++i) vmac(i, j) -= sy(i, j) * (pm(i, j) - pm(i, j - 1)) / dy;
}

// ---------------------------------------------------------------------------
// Crank–Nicolson viscous solve, eq. (16):
//   rho_h/dt U* - 1/2 div(mu_h grad U*) = rho_h [U^n/dt - A(U) + F] - Gp^{n-1/2} + 1/2 div(mu_h grad U^n) + div(mu_h grad U^n^T)
// The transpose part is taken explicitly (it vanishes for constant mu since div U = 0).
// ---------------------------------------------------------------------------
void FlowSolver::viscous_solve(double dt, const Field2D& rho_h, const Field2D& mu_h, const Field2D& au, const Field2D& av, Field2D& us, Field2D& vs) {
    const int nx = g_.nx, ny = g_.ny, ng = ng_; const double dx = g_.dx, dy = g_.dy;
    Field2D Mx(nx, ny, ng), My(nx, ny, ng); surface_tension(phi_half_, Mx, My);   // M^{n+1/2}, eq. (16)
    Field2D Lu(nx, ny, ng), Lv(nx, ny, ng), Tu(nx, ny, ng), Tv(nx, ny, ng), full_u(nx, ny, ng), full_v(nx, ny, ng);
    viscous_operator(u, v, mu_h, Lu, Lv, false);
    viscous_operator(u, v, mu_h, full_u, full_v, true);
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) { Tu(i, j) = full_u(i, j) - Lu(i, j); Tv(i, j) = full_v(i, j) - Lv(i, j); }
    Field2D rhs_u(nx, ny, ng), rhs_v(nx, ny, ng), diag(nx, ny, ng);
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        rhs_u(i, j) = rho_h(i, j) * (u(i, j) / dt - au(i, j)) - gpx(i, j) + 0.5 * Lu(i, j) + Tu(i, j) - Mx(i, j);
        rhs_v(i, j) = rho_h(i, j) * (v(i, j) / dt - av(i, j) + fp_.gravity) - gpy(i, j) + 0.5 * Lv(i, j) + Tv(i, j) - My(i, j);
        const double mxl = 0.5 * (mu_h(i - 1, j) + mu_h(i, j)), mxr = 0.5 * (mu_h(i, j) + mu_h(i + 1, j));
        const double myl = 0.5 * (mu_h(i, j - 1) + mu_h(i, j)), myr = 0.5 * (mu_h(i, j) + mu_h(i, j + 1));
        diag(i, j) = rho_h(i, j) / dt + 0.5 * ((mxl + mxr) / (dx * dx) + (myl + myr) / (dy * dy));
    }
    const bool inviscid = (fp_.mu_pos == 0.0 && fp_.mu_neg == 0.0);
    for (int comp = 0; comp < 2; ++comp) {
        Field2D& x = comp == 0 ? us : vs; const Field2D& b = comp == 0 ? rhs_u : rhs_v;
        x = comp == 0 ? u : v;
        if (inviscid) { for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) x(i, j) = b(i, j) / (rho_h(i, j) / dt); continue; }
        Field2D other(nx, ny, ng);          // dummy second component for the ghost fill of the unknown
        auto A = [&](Field2D& q, Field2D& y) {
            Field2D Lq(nx, ny, ng), Ld(nx, ny, ng), zero(nx, ny, ng);
            if (comp == 0) fill_velocity_ghosts(q, other); else fill_velocity_ghosts(other, q);
            if (comp == 0) viscous_operator(q, zero, mu_h, Lq, Ld, false); else viscous_operator(zero, q, mu_h, Ld, Lq, false);
            for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) y(i, j) = rho_h(i, j) / dt * q(i, j) - 0.5 * Lq(i, j);
        };
        if (fp_.lin_precond == "mg") {
            // BC of this velocity component: periodic / even (Neumann) / odd (Dirichlet 0), same rule as fill_velocity_ghosts
            auto mgbc = [&](int code, bool is_normal) {
                if (code == BC_PERIODIC) return MG_PERIODIC;
                if (code == BC_NOSLIPWALL) return MG_ODD;
                if (code == BC_SLIPWALL || code == BC_SYMMETRY) return is_normal ? MG_ODD : MG_EVEN;
                return MG_EVEN;
            };
            MGBCSet mb; mb.xlo = mgbc(vbc_.xlo, comp == 0); mb.xhi = mgbc(vbc_.xhi, comp == 0); mb.ylo = mgbc(vbc_.ylo, comp == 1); mb.yhi = mgbc(vbc_.yhi, comp == 1);
            Field2D sxm(nx, ny, ng), sym(nx, ny, ng), am(nx, ny, ng);
            for (int j = 0; j < ny; ++j) for (int i = 0; i <= nx; ++i) sxm(i, j) = 0.25 * (mu_h(i - 1, j) + mu_h(i, j));
            for (int j = 0; j <= ny; ++j) for (int i = 0; i < nx; ++i) sym(i, j) = 0.25 * (mu_h(i, j - 1) + mu_h(i, j));
            for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) am(i, j) = rho_h(i, j) / dt;
            MGCell mg(nx, ny, dx, dy, sxm, sym, am, mb, fp_.mg_nu, fp_.mg_nu, fp_.mg_omega);
            Precond M = [&](const Field2D& r, Field2D& z) { mg.vcycle(r, z); };
            last_visc = pcg(A, M, b, x, nx, ny, fp_.lin_rtol, fp_.lin_maxiter, false);
        } else
            last_visc = pcg(A, diag, b, x, nx, ny, fp_.lin_rtol, fp_.lin_maxiter, false);
    }
    fill_velocity_ghosts(us, vs);
}

// Cell-centred gradient of the nodal pressure: average of the two node differences on each side, S99 eq. (33).
void FlowSolver::pressure_gradient_from_nodes(const Field2D& pn, Field2D& gx, Field2D& gy) const {
    for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i) {
        gx(i, j) = (pn(i + 1, j) + pn(i + 1, j + 1) - pn(i, j) - pn(i, j + 1)) / (2 * g_.dx);
        gy(i, j) = (pn(i, j + 1) + pn(i + 1, j + 1) - pn(i, j) - pn(i + 1, j)) / (2 * g_.dy);
    }
    fill_scalar_ghosts(gx, sbc_); fill_scalar_ghosts(gy, sbc_);
}

// ---------------------------------------------------------------------------
// Nodal approximate projection, eq. (17), S99 §3.4 / A98 §3.3.
//   V = U*/dt + Gp^{n-1/2}/rho_h ;  solve  L_sigma p = D V  with sigma = 1/rho_h ;  U^{n+1} = dt (V - sigma G p)
// D: node-centred divergence of cell-centred data (transpose of eq. 33, scaled by the cell area);
// L: bilinear finite-element stiffness  sum_cells sigma_c int grad N_a . grad N_b  (the AMReX MLNodeLaplacian stencil).
// Updates p, gpx, gpy, u, v.
// ---------------------------------------------------------------------------
void FlowSolver::nodal_project(double dt, const Field2D& rho_h, const Field2D& us, const Field2D& vs) {
    const int nx = g_.nx, ny = g_.ny, ng = ng_; const double dx = g_.dx, dy = g_.dy;
    const int ni = n_nodes_x(), nj = n_nodes_y();
    Field2D Vx(nx, ny, ng), Vy(nx, ny, ng), sig(nx, ny, ng);
    for (int j = -ng; j < ny + ng; ++j) for (int i = -ng; i < nx + ng; ++i) sig(i, j) = 1.0 / rho_h(i, j);
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) { Vx(i, j) = us(i, j) / dt + gpx(i, j) / rho_h(i, j); Vy(i, j) = vs(i, j) / dt + gpy(i, j) / rho_h(i, j); }
    fill_velocity_ghosts(Vx, Vy);     // wall reflection makes the boundary-node divergence see zero normal flux (A98)
    // Local Q1 stiffness (nodes ordered 00,10,01,11):  K = sigma * (dy/dx * Kx + dx/dy * Ky)
    static const double Kx[4][4] = {{2,-2,1,-1},{-2,2,-1,1},{1,-1,2,-2},{-1,1,-2,2}};
    static const double Ky[4][4] = {{2,1,-2,-1},{1,2,-1,-2},{-2,-1,2,1},{-1,-2,1,2}};
    auto wrap = [&](int a, int n, bool per) { return per ? ((a % n) + n) % n : a; };
    auto is_dirichlet = [&](int i, int j) {   // outflow boundary nodes: p = 0
        return (vbc_.xlo == BC_OUTFLOW && i == 0) || (vbc_.xhi == BC_OUTFLOW && i == nx) || (vbc_.ylo == BC_OUTFLOW && j == 0) || (vbc_.yhi == BC_OUTFLOW && j == ny);
    };
    // Operator: y = A x by scattering over cells.
    auto A = [&](Field2D& x, Field2D& y) {
        fill_node_ghosts(x);
        for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) y(i, j) = 0.0;
        for (int jc = 0; jc < ny; ++jc) for (int ic = 0; ic < nx; ++ic) {   // every real cell once; periodic images enter through the node wrap
            const double s = sig(ic, jc);
            const int ii[4] = {ic, ic + 1, ic, ic + 1}, jj[4] = {jc, jc, jc + 1, jc + 1};
            double xv[4]; for (int b = 0; b < 4; ++b) xv[b] = x(wrap(ii[b], nx, g_.periodic_x), wrap(jj[b], ny, g_.periodic_y));
            for (int a = 0; a < 4; ++a) {
                const int ia = wrap(ii[a], nx, g_.periodic_x), ja = wrap(jj[a], ny, g_.periodic_y);
                if (ia < 0 || ia >= ni || ja < 0 || ja >= nj) continue;
                double acc = 0; for (int b = 0; b < 4; ++b) acc += s * (dy / dx * Kx[a][b] + dx / dy * Ky[a][b]) / 6.0 * xv[b];
                y(ia, ja) += acc;
            }
        }
        for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) if (is_dirichlet(i, j)) y(i, j) = x(i, j);
    };
    Field2D b(nx + 1, ny + 1, ng), diag(nx + 1, ny + 1, ng), pn(nx + 1, ny + 1, ng);
    for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) { b(i, j) = 0; diag(i, j) = 0; }
    for (int jc = 0; jc < ny; ++jc) for (int ic = 0; ic < nx; ++ic) {
        const double s = sig(ic, jc), uc = Vx(ic, jc), vc = Vy(ic, jc);
        const int ii[4] = {ic, ic + 1, ic, ic + 1}, jj[4] = {jc, jc, jc + 1, jc + 1};
        const double sx[4] = {-1, 1, -1, 1}, sy[4] = {-1, -1, 1, 1};
        for (int a = 0; a < 4; ++a) {
            const int ia = wrap(ii[a], nx, g_.periodic_x), ja = wrap(jj[a], ny, g_.periodic_y);
            if (ia < 0 || ia >= ni || ja < 0 || ja >= nj) continue;
            // (A p)_a = sum_c V_c . int_c grad N_a  ⇔  div(sigma grad p) = div V   (int_c grad N_a = (sx dy/2, sy dx/2))
            b(ia, ja) += sx[a] * 0.5 * dy * uc + sy[a] * 0.5 * dx * vc;
            diag(ia, ja) += s * (dy / dx * Kx[a][a] + dx / dy * Ky[a][a]) / 6.0;
        }
    }
    bool has_dirichlet = false;
    for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) if (is_dirichlet(i, j)) { has_dirichlet = true; b(i, j) = 0; diag(i, j) = 1; }
    for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) pn(i, j) = p(i, j);
    if (fp_.lin_precond == "mg") {
        auto dirlev = [&](int i, int j, int nxl, int nyl) {
            return (vbc_.xlo == BC_OUTFLOW && i == 0) || (vbc_.xhi == BC_OUTFLOW && i == nxl) || (vbc_.ylo == BC_OUTFLOW && j == 0) || (vbc_.yhi == BC_OUTFLOW && j == nyl);
        };
        MGNode mg(nx, ny, dx, dy, g_.periodic_x, g_.periodic_y, sig, dirlev, fp_.mg_nu, fp_.mg_nu, fp_.mg_omega);
        Precond M = [&](const Field2D& r, Field2D& z) { mg.vcycle(r, z); };
        last_nodal = pcg(A, M, b, pn, ni, nj, fp_.lin_rtol, fp_.lin_maxiter, !has_dirichlet);
    } else
        last_nodal = pcg(A, diag, b, pn, ni, nj, fp_.lin_rtol, fp_.lin_maxiter, !has_dirichlet);
    // periodic images of the solved nodes so that eq. (33) can address node nx / ny
    if (g_.periodic_x) for (int j = 0; j < nj; ++j) pn(nx, j) = pn(0, j);
    if (g_.periodic_y) for (int i = 0; i < nx + 1; ++i) pn(i, ny) = pn(i, 0);
    fill_node_ghosts(pn);
    p = pn;
    pressure_gradient_from_nodes(p, gpx, gpy);
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) { u(i, j) = dt * (Vx(i, j) - gpx(i, j) / rho_h(i, j)); v(i, j) = dt * (Vy(i, j) - gpy(i, j) / rho_h(i, j)); }
    fill_velocity_ghosts(u, v);
}

// ---------------------------------------------------------------------------
// One full time step t^n -> t^{n+1} following the outline in flow.h (S99 §3.1; redistance after the phi update):
// predictor + MAC projection, Godunov advective derivatives, phi update (+ redistance), Crank–Nicolson
// viscous solve, nodal projection, then rho^{n+1}, mu^{n+1} from phi^{n+1}.
// ---------------------------------------------------------------------------
void FlowSolver::advance(double dt) {
    const int nx = g_.nx, ny = g_.ny, ng = ng_;
    fill_velocity_ghosts(u, v); fill_scalar(phi); fill_scalar(rho); fill_scalar(mu);
    // 1. predictor + MAC projection (rho^n on faces, eq. 28)
    Field2D su(nx, ny, ng), sv(nx, ny, ng);
    forcing_for_predictor(su, sv);
    predict_mac_velocity(dt, su, sv);
    mac_project(rho);
    // 2. Godunov edge states with U^ADV; advective derivatives eq. (31)-(32)
    Field2D xed(nx, ny, ng), yed(nx, ny, ng), au(nx, ny, ng), av(nx, ny, ng), aphi(nx, ny, ng);
    godunov_edge_states(u, g_, dt, u, v, umac, vmac, fp_.godunov, xed, yed, &su);  godunov_advective_derivative(xed, yed, umac, vmac, g_, au);
    godunov_edge_states(v, g_, dt, u, v, umac, vmac, fp_.godunov, xed, yed, &sv);  godunov_advective_derivative(xed, yed, umac, vmac, g_, av);
    Field2D phi_old = phi, rho_h(nx, ny, ng), mu_h(nx, ny, ng);
    if (fp_.do_phi) {
        godunov_edge_states(phi, g_, dt, u, v, umac, vmac, fp_.godunov, xed, yed); godunov_advective_derivative(xed, yed, umac, vmac, g_, aphi);
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) phi(i, j) -= dt * aphi(i, j);   // eq. (12)
        fill_scalar(phi);
        if (fp_.do_reinit && !fp_.reinit_after_projection && ((step + 1) % fp_.reinit_interval == 0)) redistance(phi, g_, sbc_, fp_.reinit);
        Field2D& phi_h = phi_half_; phi_h.define(nx, ny, ng);
        for (int j = -ng; j < ny + ng; ++j) for (int i = -ng; i < nx + ng; ++i) phi_h(i, j) = 0.5 * (phi_old(i, j) + phi(i, j));   // eq. (13)
        material_from_phi(phi_h, fp_.rho_pos, fp_.rho_neg, eps_h(), rho_h);                                                   // eq. (14)
        material_from_phi(phi_h, fp_.mu_pos, fp_.mu_neg, eps_h(), mu_h);                                                      // eq. (15)
    } else { rho_h = rho; mu_h = mu; phi_half_ = phi; }
    // 3. viscous solve, 4. nodal projection
    Field2D us(nx, ny, ng), vs(nx, ny, ng);
    viscous_solve(dt, rho_h, mu_h, au, av, us, vs);
    nodal_project(dt, rho_h, us, vs);
    if (fp_.do_phi) {
        if (fp_.do_reinit && fp_.reinit_after_projection && ((step + 1) % fp_.reinit_interval == 0)) redistance(phi, g_, sbc_, fp_.reinit);
        set_material_from_phi();                                   // rho^{n+1}, mu^{n+1} from phi^{n+1}
    }
    t += dt; ++step;
}

// Start-up (S99 §3.6): project U^0 to be discretely divergence-free, then run init_iters dummy steps
// keeping only the pressure gradient Gp^{-1/2} so the first real step has a consistent lagged pressure.
void FlowSolver::initialise() {
    fill_velocity_ghosts(u, v); fill_scalar(phi); set_material_from_phi(); fill_scalar(rho); fill_scalar(mu);
    gpx.fill(0); gpy.fill(0); p.fill(0);
    // Project the initial velocity (S99 §3.6): nodal projection with dt = 1 and U* = U^0, then discard the pressure.
    { Field2D us = u, vs = v; nodal_project(1.0, rho, us, vs); gpx.fill(0); gpy.fill(0); p.fill(0); }
    // Iterate for p^{1/2}: take a step, keep Gp, restore U and phi.
    Field2D u0 = u, v0 = v, phi0 = phi, rho0 = rho, mu0 = mu;
    for (int it = 0; it < fp_.init_iters; ++it) {
        const double dt = estimate_dt();
        advance(dt);
        u = u0; v = v0; phi = phi0; rho = rho0; mu = mu0; t = 0; step = 0;
        if (fp_.verbose) std::printf("  init iteration %d: dt=%.3e nodal iters=%d\n", it + 1, dt, last_nodal.iterations);
    }
}

// Time-step estimate, S99 §3.1.1: CFL on |U|, gravity, capillary and (two-viscosity) explicit-viscous
// constraints, all multiplied by cfl; ns.fixed_dt overrides everything.
double FlowSolver::estimate_dt() const {
    double umax = 1e-300, dt = 1e30;
    for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i) umax = std::max(umax, std::max(std::fabs(u(i, j)), std::fabs(v(i, j))));
    dt = std::min(dt, g_.dxmin() / umax);
    if (fp_.gravity != 0.0) dt = std::min(dt, std::sqrt(2.0 * g_.dxmin() / std::fabs(fp_.gravity)));
    if (fp_.sigma > 0.0) dt = std::min(dt, std::sqrt((fp_.rho_pos + fp_.rho_neg) * std::pow(g_.dxmin(), 3) / (8.0 * PI_ * fp_.sigma)));   // capillary, S99 §3.1.1
    // Viscous constraint of S99 §3.1.1, 3/14 rho dx^2/mu, evaluated per phase (rho and mu of the same fluid);
    // needed because the transpose part of the stress is explicit (it vanishes for constant mu).
    if (fp_.mu_pos != fp_.mu_neg) {
        double r = 1e300;
        if (fp_.mu_pos > 0) r = std::min(r, fp_.rho_pos / fp_.mu_pos);
        if (fp_.mu_neg > 0) r = std::min(r, fp_.rho_neg / fp_.mu_neg);
        dt = std::min(dt, 3.0 / 14.0 * r * g_.dxmin() * g_.dxmin());
    }
    dt *= fp_.cfl;
    if (fp_.fixed_dt > 0) dt = fp_.fixed_dt;
    return dt;
}

// Diagnostic: integral of 1/2 rho |U|^2 over the domain.
double FlowSolver::kinetic_energy() const {
    double e = 0; for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i) e += 0.5 * rho(i, j) * (u(i, j) * u(i, j) + v(i, j) * v(i, j));
    return e * g_.cell_area();
}
// Diagnostic: max |div U^ADV| of the face velocities (should be ~ solver tolerance after the MAC projection).
double FlowSolver::max_divergence_mac() const {
    double m = 0; for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i) m = std::max(m, std::fabs((umac(i + 1, j) - umac(i, j)) / g_.dx + (vmac(i, j + 1) - vmac(i, j)) / g_.dy));
    return m;
}
// Diagnostic: max of the node-centred divergence (transpose of eq. 33) of the cell-centred velocity, interior nodes.
double FlowSolver::max_divergence_nodal() const {
    double m = 0;
    for (int j = 1; j < g_.ny; ++j) for (int i = 1; i < g_.nx; ++i)
        m = std::max(m, std::fabs((u(i, j) + u(i, j - 1) - u(i - 1, j) - u(i - 1, j - 1)) / (2 * g_.dx) + (v(i, j) + v(i - 1, j) - v(i, j - 1) - v(i - 1, j - 1)) / (2 * g_.dy)));
    return m;
}
// Diagnostic: w = v_x - u_y by central differences at cell centres.
void FlowSolver::vorticity(Field2D& w) const {
    for (int j = 0; j < g_.ny; ++j) for (int i = 0; i < g_.nx; ++i)
        w(i, j) = (v(i + 1, j) - v(i - 1, j)) / (2 * g_.dx) - (u(i, j + 1) - u(i, j - 1)) / (2 * g_.dy);
}
