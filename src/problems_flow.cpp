// problems_flow.cpp — set-up of the Navier–Stokes test cases listed in problems_flow.h.
#include "problems_flow.h"
#include <cmath>
#include <stdexcept>

static const double PI = 3.141592653589793238462643383279502884197;

// Read geometry, velocity BCs, material and solver parameters from the input file (`ns.*` key names),
// then attach the initial velocity and (Taylor–Green) exact solution for prob.type.
FlowProblem make_flow_problem(const Params& p) {
    FlowProblem pr; pr.name = p.get("prob.type", "taylorgreen");
    auto lo = p.get_vec<double>("geometry.prob_lo", {0, 0}), hi = p.get_vec<double>("geometry.prob_hi", {1, 1});
    auto n  = p.get_vec<int>("amr.n_cell", {64, 64});
    pr.geom.define(n[0], n[1], lo[0], hi[0], lo[1], hi[1]);
    auto per = p.get_vec<int>("geometry.is_periodic", {0, 0});
    pr.geom.periodic_x = per[0]; pr.geom.periodic_y = per[1];
    auto lobc = p.get_vec<int>("ns.lo_bc", {0, 0}), hibc = p.get_vec<int>("ns.hi_bc", {0, 0});
    pr.vbc.xlo = per[0] ? 0 : lobc[0]; pr.vbc.xhi = per[0] ? 0 : hibc[0]; pr.vbc.ylo = per[1] ? 0 : lobc[1]; pr.vbc.yhi = per[1] ? 0 : hibc[1];
    // Tangential velocity of a moving no-slip wall, one value per side (see VelBC in flow.h).
    auto wlo = p.get_vec<double>("ns.wall_vel_lo", {0, 0}), whi = p.get_vec<double>("ns.wall_vel_hi", {0, 0});
    pr.vbc.vt_xlo = wlo[0]; pr.vbc.vt_ylo = wlo[1]; pr.vbc.vt_xhi = whi[0]; pr.vbc.vt_yhi = whi[1];
    pr.t_final = p.get("stop_time", 0.0);
    pr.steady_tol = p.get("ns.steady_tol", 0.0);
    FlowParams& fp = pr.fp;
    fp.rho = p.get("ns.rho", 1.0);
    fp.mu  = p.get("ns.mu", 0.0);
    fp.gravity = p.get("ns.gravity", 0.0);
    fp.cfl = p.get("ns.cfl", 0.5); fp.fixed_dt = p.get("ns.fixed_dt", 0.0);
    fp.init_iters = p.get("ns.init_iter", 3);
    fp.lin_rtol = p.get("lin.rtol", 1e-10); fp.lin_maxiter = p.get("lin.maxiter", 5000); fp.lin_precond = p.get("lin.precond", "mg"); fp.mg_nu = p.get("lin.mg_nu", 2); fp.mg_omega = p.get("lin.mg_omega", 0.8);
    fp.verbose = p.get("ns.v", 1);
    fp.godunov.fourth_order_slopes = p.get("godunov.fourth_order_slopes", 1) != 0; fp.godunov.use_transverse = p.get("godunov.use_transverse", 1) != 0;

    if (pr.name == "taylorgreen") {
        // u =  sin(2πx) cos(2πy) e^{-8π²νt},  v = -cos(2πx) sin(2πy) e^{-8π²νt},  p = -¼ (cos 4πx + cos 4πy) e^{-16π²νt}
        const double nu = fp.mu / fp.rho;
        pr.vel0 = [](double x, double y, double& u, double& v) { u = std::sin(2 * PI * x) * std::cos(2 * PI * y); v = -std::cos(2 * PI * x) * std::sin(2 * PI * y); };
        pr.exact = [=](double x, double y, double t, double& u, double& v, double& pp) {
            const double d = std::exp(-8 * PI * PI * nu * t);
            u = std::sin(2 * PI * x) * std::cos(2 * PI * y) * d; v = -std::cos(2 * PI * x) * std::sin(2 * PI * y) * d;
            pp = 0.25 * (std::cos(4 * PI * x) + std::cos(4 * PI * y)) * d * d;   // sign verified from Bernoulli: -p_x = u u_x + v u_y = pi sin(4 pi x)
        };
        pr.has_exact = true;
    } else if (pr.name == "cavity") {
        // Lid-driven cavity (Ghia, Ghia & Shin, JCP 48 (1982) 387): fluid at rest in a square box with
        // no-slip walls, the top wall sliding at u = prob.lid_velocity (1 in the benchmark).
        // Re = rho * U * L / mu with L the cavity side.
        const double U = p.get("prob.lid_velocity", 1.0);
        if (!p.has("ns.lo_bc") && !p.has("ns.hi_bc")) pr.vbc.xlo = pr.vbc.xhi = pr.vbc.ylo = pr.vbc.yhi = BC_NOSLIPWALL;
        if (!p.has("ns.wall_vel_hi")) pr.vbc.vt_yhi = U;
        pr.vel0 = [](double, double, double& u, double& v) { u = 0; v = 0; };
    } else throw std::runtime_error("unknown flow problem " + pr.name);
    return pr;
}
