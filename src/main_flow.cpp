// main_flow.cpp — driver for the Navier–Stokes problems (V4–V7).
#include "main_flow.h"
#include "problems_flow.h"
#include "io.h"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sys/stat.h>

// Build the FlowSolver from the inputs, set U^0 and phi^0, run the initial pressure iterations, then
// time-step to stop_time writing diag.csv / VTK snapshots / summary.txt (with Taylor–Green errors if exact).
int run_flow(const Params& p) {
    FlowProblem pr = make_flow_problem(p);
    const Geometry& g = pr.geom;
    const std::string out = p.get("out_dir", "out"); mkdir(out.c_str(), 0755);
    const int plot_int = p.get("amr.plot_int", 0), max_step = p.get("max_step", 1000000), diag_int = p.get("diag.interval", 10);
    FlowSolver S(g, pr.vbc, pr.fp);
    const int ng = S.u.ng();
    for (int j = -ng; j < g.ny + ng; ++j) for (int i = -ng; i < g.nx + ng; ++i) {
        pr.vel0(g.xc(i), g.yc(j), S.u(i, j), S.v(i, j)); S.phi(i, j) = pr.phi0(g.xc(i), g.yc(j));
    }
    std::cout << "ls2d-baseline flow  problem=" << pr.name << " grid=" << g.nx << "x" << g.ny << " rho=" << pr.fp.rho_pos << "/" << pr.fp.rho_neg
              << " mu=" << pr.fp.mu_pos << "/" << pr.fp.mu_neg << " g=" << pr.fp.gravity << "\n";
    S.initialise();
    auto write_plot = [&](const std::string& fn) {
        Field2D w(g.nx, g.ny, ng), pc(g.nx, g.ny, ng); S.vorticity(w);
        for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) pc(i, j) = 0.25 * (S.p(i, j) + S.p(i + 1, j) + S.p(i, j + 1) + S.p(i + 1, j + 1));
        write_vtk(fn, g, {{"u", &S.u}, {"v", &S.v}, {"p", &pc}, {"phi", &S.phi}, {"rho", &S.rho}, {"vorticity", &w}});
    };
    write_plot(out + "/flow_0000.vtk");
    std::ofstream diag(out + "/diag.csv");
    diag << "step,time,dt,kinetic_energy,max_div_mac,max_div_nodal,area_pos,interface_ymin,interface_ymax,mac_iters,visc_iters,nodal_iters,umax,yc,vc,circularity,pmax_minus_pmin\n";
    auto interface_extent = [&](double& ymin, double& ymax) {   // y-range of the zero contour (RT spike/bubble tips)
        ymin = 1e30; ymax = -1e30;
        for (int j = 0; j < g.ny - 1; ++j) for (int i = 0; i < g.nx; ++i) {
            const double a = S.phi(i, j), b = S.phi(i, j + 1);
            if (a * b <= 0 && a != b) { const double y = g.yc(j) + g.dy * a / (a - b); ymin = std::min(ymin, y); ymax = std::max(ymax, y); }
        }
    };
    // Hysing 2009 benchmark quantities for the phi > 0 phase: centre of mass y_c, rise velocity v_c, circularity
    // (perimeter of the area-equivalent circle / actual perimeter, from the smoothed Heaviside and delta).
    auto bubble_quantities = [&](double& yc, double& vc, double& circ) {
        const double eps = pr.fp.reinit.alpha_heaviside * g.dxmin(); double A = 0, Y = 0, V = 0, L = 0;
        for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
            const double H = heaviside_eps(S.phi(i, j), eps); A += H; Y += H * g.yc(j); V += H * S.v(i, j);
            const double gx = (S.phi(i + 1, j) - S.phi(i - 1, j)) / (2 * g.dx), gy = (S.phi(i, j + 1) - S.phi(i, j - 1)) / (2 * g.dy);
            L += delta_eps(S.phi(i, j), eps) * std::sqrt(gx * gx + gy * gy);
        }
        yc = Y / A; vc = V / A; circ = 2 * std::sqrt(3.141592653589793 * A * g.cell_area()) / (L * g.cell_area());
    };
    auto log = [&](double dt) {
        double ymin, ymax; interface_extent(ymin, ymax);
        double yc, vc, circ; bubble_quantities(yc, vc, circ);
        double pmin = 1e300, pmax = -1e300; for (int j = 0; j <= g.ny; ++j) for (int i = 0; i <= g.nx; ++i) { pmin = std::min(pmin, S.p(i, j)); pmax = std::max(pmax, S.p(i, j)); }
        double umax = 0; for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) umax = std::max(umax, std::max(std::fabs(S.u(i, j)), std::fabs(S.v(i, j))));
        diag << S.step << "," << S.t << "," << dt << "," << S.kinetic_energy() << "," << S.max_divergence_mac() << "," << S.max_divergence_nodal() << ","
             << volume_heaviside(S.phi, g, pr.fp.reinit.alpha_heaviside * g.dxmin()) << "," << ymin << "," << ymax << ","
             << S.last_mac.iterations << "," << S.last_visc.iterations << "," << S.last_nodal.iterations << "," << umax << ","
             << yc << "," << vc << "," << circ << "," << pmax - pmin << "\n";
    };
    log(0);
    while (S.t < pr.t_final - 1e-12 && S.step < max_step) {
        double dt = S.estimate_dt(); if (S.t + dt > pr.t_final) dt = pr.t_final - S.t;
        S.advance(dt);
        { double um = 0; for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) um = std::max(um, std::fabs(S.u(i, j)) + std::fabs(S.v(i, j)));
          if (!(um < 1e6) || dt < 1e-12) { std::cerr << "BLOW-UP at step " << S.step << " t=" << S.t << " umax=" << um << " dt=" << dt << "\n"; return 2; } }
        if (S.step % diag_int == 0 || S.t >= pr.t_final - 1e-12) log(dt);
        if (plot_int > 0 && S.step % plot_int == 0) { char b[64]; std::snprintf(b, 64, "/flow_%04d.vtk", S.step); write_plot(out + b); }
        if (pr.fp.verbose && S.step % 50 == 0) std::printf("step %d t=%.5f dt=%.3e KE=%.6e divMAC=%.2e divN=%.2e iters mac/visc/nodal=%d/%d/%d\n",
            S.step, S.t, dt, S.kinetic_energy(), S.max_divergence_mac(), S.max_divergence_nodal(), S.last_mac.iterations, S.last_visc.iterations, S.last_nodal.iterations);
    }
    write_plot(out + "/flow_final.vtk");
    std::ofstream sum(out + "/summary.txt"); sum.precision(8);
    sum << "problem " << pr.name << "\ngrid " << g.nx << " " << g.ny << "\nsteps " << S.step << "\nt_final " << S.t << "\nkinetic_energy " << S.kinetic_energy()
        << "\nmax_div_mac " << S.max_divergence_mac() << "\nmax_div_nodal " << S.max_divergence_nodal() << "\n";
    if (pr.has_exact) {
        double e2u = 0, e2v = 0, emax = 0, e2p = 0, pm_num = 0, pm_ex = 0;
        for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
            double ue, ve, pe; pr.exact(g.xc(i), g.yc(j), S.t, ue, ve, pe);
            e2u += std::pow(S.u(i, j) - ue, 2); e2v += std::pow(S.v(i, j) - ve, 2); emax = std::max(emax, std::fabs(S.u(i, j) - ue));
            double pn = 0.25 * (S.p(i, j) + S.p(i + 1, j) + S.p(i, j + 1) + S.p(i + 1, j + 1)); pm_num += pn; pm_ex += pe;
        }
        pm_num /= g.nx * g.ny; pm_ex /= g.nx * g.ny;
        for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
            double ue, ve, pe; pr.exact(g.xc(i), g.yc(j), S.t, ue, ve, pe);
            double pn = 0.25 * (S.p(i, j) + S.p(i + 1, j) + S.p(i, j + 1) + S.p(i + 1, j + 1));
            e2p += std::pow((pn - pm_num) - (pe - pm_ex), 2);
        }
        sum << "L2_error_u " << std::sqrt(e2u * g.cell_area()) << "\nL2_error_v " << std::sqrt(e2v * g.cell_area()) << "\nLinf_error_u " << emax
            << "\nL2_error_p " << std::sqrt(e2p * g.cell_area()) << "\n";
        std::cout << "L2_error_u " << std::sqrt(e2u * g.cell_area()) << "  Linf_u " << emax << "  L2_p " << std::sqrt(e2p * g.cell_area()) << "\n";
    }
    double ymin, ymax; interface_extent(ymin, ymax); sum << "interface_ymin " << ymin << "\ninterface_ymax " << ymax << "\n";
    return 0;
}
