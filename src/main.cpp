// main.cpp — ls2d-baseline driver for the interface-transport tests (prescribed
// velocity + Sussman redistancing).  Navier–Stokes problems (prob.flow = 1) are
// dispatched to run_flow() in main_flow.cpp.
//
//   ./ls2d inputs.zalesak
//
// Outputs (in `out_dir`): diag.csv (time series), phi_XXXX.vtk snapshots,
// summary.txt (final errors) — the numbers compared in docs/results.md.
#include "grid.h"
#include "params.h"
#include "bc.h"
#include "levelset.h"
#include "advect.h"
#include "godunov.h"
#include "problems.h"
#include "io.h"
#include "main_flow.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

// Heaviside-difference error, Sussman 1999 eq. (80) / Enright 2002 eq. (14):
// sum over cells of the sub-cell integral of |H(phi_exact) - H(phi)|, with phi
// bilinearly interpolated to nsub x nsub sub-cell midpoints.
static double interface_error(const Field2D& phi, const Geometry& g, const std::function<int(double, double)>& exact, int nsub) {
    double err = 0;
    for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
        for (int sj = 0; sj < nsub; ++sj) for (int si = 0; si < nsub; ++si) {
            const double fx = (si + 0.5) / nsub - 0.5, fy = (sj + 0.5) / nsub - 0.5;
            const int i0 = fx < 0 ? i - 1 : i, j0 = fy < 0 ? j - 1 : j;
            const double tx = fx < 0 ? fx + 1.0 : fx, ty = fy < 0 ? fy + 1.0 : fy;
            const double v = (1 - tx) * (1 - ty) * phi(i0, j0) + tx * (1 - ty) * phi(i0 + 1, j0)
                           + (1 - tx) * ty * phi(i0, j0 + 1) + tx * ty * phi(i0 + 1, j0 + 1);
            const int Hc = v > 0 ? 1 : 0;
            const int He = exact(g.xc(i) + fx * g.dx, g.yc(j) + fy * g.dy);
            err += std::abs(Hc - He);
        }
    return err * g.cell_area() / (nsub * nsub);
}

// Driver: read inputs, build the problem, (optionally) redistance the initial phi, then loop
// advect -> redistance -> log until stop_time.  Delegates to run_flow() when prob.flow = 1.
int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: ls2d <inputs>\n"; return 1; }
    Params p; p.read(argv[1]);
    if (p.get("prob.flow", 0)) return run_flow(p);      // Navier–Stokes problems (results.md §2–3); otherwise prescribed-velocity interface tests (§1)
    Problem pr = make_problem(p);
    const Geometry& g = pr.geom;
    const std::string out = p.get("out_dir", "out");
    mkdir(out.c_str(), 0755);

    ReinitParams rp;
    rp.alpha_heaviside = p.get("ns.epsilon", 2.0);
    rp.alpha_sign      = p.get("reinit.alpha_sign", rp.alpha_heaviside);
    rp.alpha_delta     = p.get("reinit.alpha_delta", -1.0);
    rp.n_iter          = p.get("ns.number_of_reinit", 4);
    rp.dtau_factor     = p.get("reinit.dtau_factor", 0.5);
    rp.upwind          = p.get("reinit.upwind", "godunov");
    rp.sign            = p.get("reinit.sign", "frozen");
    rp.volume_fix      = p.get("reinit.volume_fix", 1) != 0;
    rp.vf_weight_centre = p.get("reinit.vf_weight_centre", 16.0);
    rp.vf_weight_neigh  = p.get("reinit.vf_weight_neigh", 1.0);
    rp.sign_guard       = p.get("reinit.sign_guard", 1) != 0;
    const int do_reinit = p.get("ns.do_reinit", 1);
    const int reinit_interval = p.get("ns.lev0step_of_reinit", 1);
    const int init_reinit_iters = p.get("reinit.init_iters", 0);    // extra redistancing of the initial data (e.g. CSG shapes)
    const double cfl = p.get("ns.cfl", 0.5);
    const double fixed_dt = p.get("ns.fixed_dt", 0.0);
    const int plot_int = p.get("amr.plot_int", 0), max_step = p.get("max_step", 1000000);
    const int nsub = p.get("diag.nsub", 10);
    const std::string advection = p.get("ns.advection_scheme", "Godunov_PLM");   // Godunov_PLM (paper, S99 §3.2) | ENO2_RK2 (simpler reference scheme, advect.cpp)
    GodunovParams gp; gp.fourth_order_slopes = p.get("godunov.fourth_order_slopes", 1) != 0; gp.use_transverse = p.get("godunov.use_transverse", 1) != 0;
    const double eps_h = rp.alpha_heaviside * g.dxmin();

    std::cout << "ls2d-baseline  problem=" << pr.name << "  grid=" << g.nx << "x" << g.ny << "  dx=" << g.dx << "\n";
    p.print(std::cout);

    const int ng = 4;   // 4th-order PLM slopes on the predictor stencil need 4 ghost layers
    Field2D phi(g.nx, g.ny, ng);
    for (int j = -ng; j < g.ny + ng; ++j) for (int i = -ng; i < g.nx + ng; ++i) phi(i, j) = pr.phi0(g.xc(i), g.yc(j));
    fill_scalar_ghosts(phi, pr.bc);
    if (init_reinit_iters > 0) { ReinitParams r0 = rp; r0.n_iter = init_reinit_iters; redistance(phi, g, pr.bc, r0); }

    const double area0_h = volume_heaviside(phi, g, eps_h), area0_s = volume_sharp(phi, g, nsub);
    const double L0 = perimeter_heaviside(phi, g, eps_h);   // interface length, used to normalise the error as in Enright 2002 eq. (14)
    const double err0 = interface_error(phi, g, pr.indicator, nsub);
    std::ofstream diag(out + "/diag.csv");
    diag << "step,time,dt,area_heaviside,area_sharp,area_loss_pct,interface_error,grad_err_max_band,grad_err_mean_band,sign_guard_hits\n";
    auto log = [&](int step, double t, double dt, const ReinitStats* st) {
        double gmax, gmean; grad_norm_error(phi, g, eps_h, gmax, gmean);
        const double as = volume_sharp(phi, g, nsub);
        diag << step << "," << t << "," << dt << "," << volume_heaviside(phi, g, eps_h) << "," << as << ","
             << 100.0 * (pr.exact_area - as) / pr.exact_area << "," << interface_error(phi, g, pr.indicator, nsub) << ","
             << gmax << "," << gmean << "," << (st ? st->sign_guard_hits : 0) << "\n";
    };
    log(0, 0.0, 0.0, nullptr);
    write_vtk(out + "/phi_0000.vtk", g, {{"phi", &phi}});

    double t = 0; int step = 0;
    ReinitStats st;
    while (t < pr.t_final - 1e-12 && step < max_step) {
        double dt = fixed_dt > 0 ? fixed_dt : cfl * g.dxmin() / std::max(max_velocity(g, pr.vel, t), 1e-300);
        if (t + dt > pr.t_final) dt = pr.t_final - t;
        if (pr.name != "disk") {
            if (advection == "Godunov_PLM") godunov_step(phi, g, pr.bc, pr.vel, t, dt, gp);
            else advect_step_rk2(phi, g, pr.bc, pr.vel, t, dt);
        }
        t += dt; ++step;
        if (do_reinit && step % reinit_interval == 0) st = redistance(phi, g, pr.bc, rp);
        if (step % p.get("diag.interval", 10) == 0 || t >= pr.t_final - 1e-12) log(step, t, dt, &st);
        if (plot_int > 0 && step % plot_int == 0) { char b[64]; std::snprintf(b, 64, "/phi_%04d.vtk", step); write_vtk(out + b, g, {{"phi", &phi}}); }
        if (step % 100 == 0) std::cout << "step " << step << " t=" << t << "\n";
    }
    write_vtk(out + "/phi_final.vtk", g, {{"phi", &phi}});
    write_ascii(out + "/phi_final.txt", g, phi);

    const double area_s = volume_sharp(phi, g, nsub), err = interface_error(phi, g, pr.indicator, nsub);
    double gmax, gmean; grad_norm_error(phi, g, eps_h, gmax, gmean);
    std::ofstream sum(out + "/summary.txt"); sum.precision(6);
    sum << "problem " << pr.name << "\ngrid " << g.nx << " " << g.ny << "\nsteps " << step << "\nt_final " << t << "\n"
        << "exact_area " << pr.exact_area << "\narea_sharp_initial " << area0_s << "\narea_sharp_final " << area_s
        << "\narea_loss_pct " << 100.0 * (pr.exact_area - area_s) / pr.exact_area
        << "\narea_heaviside_initial " << area0_h << "\narea_heaviside_final " << volume_heaviside(phi, g, eps_h)
        << "\nperimeter_initial " << L0
        << "\ninterface_error_initial " << err0 << "\ninterface_error_final " << err
        << "\ninterface_error_final_over_L " << err / L0
        << "\ngrad_err_max_band " << gmax << "\ngrad_err_mean_band " << gmean << "\nsign_guard_hits_last " << st.sign_guard_hits << "\n";
    std::cout << "done: area_loss_pct=" << 100.0 * (pr.exact_area - area_s) / pr.exact_area << "  interface_error=" << err << "  error/L=" << err / L0
              << "  grad_err_max_band=" << gmax << "\n";
    return 0;
}
