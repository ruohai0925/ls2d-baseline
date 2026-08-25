// problems_flow.cpp — set-up of the Navier–Stokes test cases listed in problems_flow.h.
#include "problems_flow.h"
#include <cmath>
#include <stdexcept>

static const double PI = 3.141592653589793238462643383279502884197;

// Read geometry, velocity BCs, material / solver / reinit parameters from the input file (`ns.*` key names),
// then attach the initial velocity, initial phi and (Taylor–Green) exact solution for prob.type.
FlowProblem make_flow_problem(const Params& p) {
    FlowProblem pr; pr.name = p.get("prob.type", "taylorgreen");
    auto lo = p.get_vec<double>("geometry.prob_lo", {0, 0}), hi = p.get_vec<double>("geometry.prob_hi", {1, 1});
    auto n  = p.get_vec<int>("amr.n_cell", {64, 64});
    pr.geom.define(n[0], n[1], lo[0], hi[0], lo[1], hi[1]);
    auto per = p.get_vec<int>("geometry.is_periodic", {0, 0});
    pr.geom.periodic_x = per[0]; pr.geom.periodic_y = per[1];
    auto lobc = p.get_vec<int>("ns.lo_bc", {0, 0}), hibc = p.get_vec<int>("ns.hi_bc", {0, 0});
    pr.vbc.xlo = per[0] ? 0 : lobc[0]; pr.vbc.xhi = per[0] ? 0 : hibc[0]; pr.vbc.ylo = per[1] ? 0 : lobc[1]; pr.vbc.yhi = per[1] ? 0 : hibc[1];
    pr.t_final = p.get("stop_time", 0.0);
    FlowParams& fp = pr.fp;
    fp.do_phi = p.get("ns.do_phi", 0);
    fp.rho_pos = p.get("ns.rho_w", 1.0); fp.rho_neg = p.get("ns.rho_a", fp.rho_pos);
    fp.mu_pos = p.get("ns.mu_w", p.get("ns.vel_visc_coef", 0.0)); fp.mu_neg = p.get("ns.mu_a", fp.mu_pos);
    fp.gravity = p.get("ns.gravity", 0.0); fp.sigma = p.get("ns.sigma", 0.0);
    fp.cfl = p.get("ns.cfl", 0.5); fp.fixed_dt = p.get("ns.fixed_dt", 0.0);
    fp.do_reinit = p.get("ns.do_reinit", 1); fp.reinit_interval = p.get("ns.lev0step_of_reinit", 1);
    fp.reinit_after_projection = p.get("reinit.after_projection", 0);
    fp.init_iters = p.get("ns.init_iter", 3);
    fp.lin_rtol = p.get("lin.rtol", 1e-10); fp.lin_maxiter = p.get("lin.maxiter", 5000); fp.lin_precond = p.get("lin.precond", "mg"); fp.mg_nu = p.get("lin.mg_nu", 2); fp.mg_omega = p.get("lin.mg_omega", 0.8);
    fp.verbose = p.get("ns.v", 1);
    fp.godunov.fourth_order_slopes = p.get("godunov.fourth_order_slopes", 1) != 0; fp.godunov.use_transverse = p.get("godunov.use_transverse", 1) != 0;
    ReinitParams& rp = fp.reinit;
    rp.alpha_heaviside = p.get("ns.epsilon", 2.0); rp.alpha_sign = p.get("reinit.alpha_sign", rp.alpha_heaviside); rp.alpha_delta = p.get("reinit.alpha_delta", -1.0);
    rp.n_iter = p.get("ns.number_of_reinit", 4); rp.dtau_factor = p.get("reinit.dtau_factor", 0.5);
    rp.upwind = p.get("reinit.upwind", "godunov"); rp.sign = p.get("reinit.sign", "frozen");
    rp.volume_fix = p.get("reinit.volume_fix", 1) != 0; rp.vf_weight_centre = p.get("reinit.vf_weight_centre", 16.0); rp.vf_weight_neigh = p.get("reinit.vf_weight_neigh", 1.0);
    rp.sign_guard = p.get("reinit.sign_guard", 1) != 0;

    if (pr.name == "taylorgreen") {
        // u =  sin(2πx) cos(2πy) e^{-8π²νt},  v = -cos(2πx) sin(2πy) e^{-8π²νt},  p = -¼ (cos 4πx + cos 4πy) e^{-16π²νt}
        const double nu = fp.mu_pos / fp.rho_pos;
        pr.vel0 = [](double x, double y, double& u, double& v) { u = std::sin(2 * PI * x) * std::cos(2 * PI * y); v = -std::cos(2 * PI * x) * std::sin(2 * PI * y); };
        pr.phi0 = [](double, double) { return 1.0; };
        pr.exact = [=](double x, double y, double t, double& u, double& v, double& pp) {
            const double d = std::exp(-8 * PI * PI * nu * t);
            u = std::sin(2 * PI * x) * std::cos(2 * PI * y) * d; v = -std::cos(2 * PI * x) * std::sin(2 * PI * y) * d;
            pp = 0.25 * (std::cos(4 * PI * x) + std::cos(4 * PI * y)) * d * d;   // sign verified from Bernoulli: -p_x = u u_x + v u_y = pi sin(4 pi x)
        };
        pr.has_exact = true;
    } else if (pr.name == "hydrostatic") {
        const double yi = p.get("prob.interface_y", 0.5);
        pr.vel0 = [](double, double, double& u, double& v) { u = 0; v = 0; };
        pr.phi0 = [=](double, double y) { return y - yi; };          // phi>0 above the interface (rho_pos = upper fluid)
    } else if (pr.name == "rt") {
        // Guermond–Salgado 2009 §5.2 / Tryggvason 1988: heavy fluid above, eta(x) = -amp*d*cos(2πx/d), d = Lx.
        const double d = pr.geom.xhi - pr.geom.xlo, amp = p.get("prob.perturbation_amplitude", 0.1), y0 = p.get("prob.interface_y", 0.5 * (pr.geom.ylo + pr.geom.yhi));
        pr.vel0 = [](double, double, double& u, double& v) { u = 0; v = 0; };
        // Interface y0 + amp*d*cos(2πx/d) in absolute x: on [0,1] the heavy
        // fluid dips at x = 1/2 (spike at the centre); on (-1/2,1/2) the bubble sits at x = 0, as in Guermond–Salgado's frames.
        pr.phi0 = [=](double x, double y) { return y - (y0 + amp * d * std::cos(2 * PI * x / d)); };   // phi>0 = heavy (upper) fluid
    } else if (pr.name == "bubble") {
        // Circular bubble/drop of the phi > 0 fluid, radius R at (cx, cy) — static drop (Laplace pressure, spurious
        // currents) or Hysing et al. 2009 rising bubble (rho_w/mu_w = bubble, rho_a/mu_a = surrounding liquid).
        const double cx = p.get("prob.blob_x", 0.5), cy = p.get("prob.blob_y", 0.5), R = p.get("prob.blob_radius", 0.25);
        pr.vel0 = [](double, double, double& u, double& v) { u = 0; v = 0; };
        pr.phi0 = [=](double x, double y) { return R - std::hypot(x - cx, y - cy); };
    } else throw std::runtime_error("unknown flow problem " + pr.name);
    return pr;
}
