// main_flow.cpp — ls2d-baseline (single-phase) driver.
//
//   ./ls2d tests/inputs.taylorgreen
//
// Outputs (in `out_dir`): diag.csv (time series), flow_XXXX.vtk snapshots, summary.txt (final
// numbers), centreline_u.txt / centreline_v.txt (profiles through the domain centre) — the
// quantities compared with the literature in docs/results.md.
#include "main_flow.h"
#include "problems_flow.h"
#include "io.h"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sys/stat.h>

// Linear interpolation of a cell-centred field along x to the abscissa xq, for every row j.
// (For an even cell count the domain centre falls on a face, i.e. halfway between two cell centres.)
static void column_profile(const Field2D& f, const Geometry& g, double xq, std::vector<double>& y, std::vector<double>& val) {
    const double s = (xq - g.xlo) / g.dx - 0.5;                 // continuous i index of xq
    const int i0 = std::max(0, std::min(g.nx - 2, static_cast<int>(std::floor(s))));
    const double w = (xq - g.xc(i0)) / g.dx;
    y.clear(); val.clear();
    for (int j = 0; j < g.ny; ++j) { y.push_back(g.yc(j)); val.push_back((1 - w) * f(i0, j) + w * f(i0 + 1, j)); }
}
// Same along y, to the ordinate yq, for every column i.
static void row_profile(const Field2D& f, const Geometry& g, double yq, std::vector<double>& x, std::vector<double>& val) {
    const double s = (yq - g.ylo) / g.dy - 0.5;
    const int j0 = std::max(0, std::min(g.ny - 2, static_cast<int>(std::floor(s))));
    const double w = (yq - g.yc(j0)) / g.dy;
    x.clear(); val.clear();
    for (int i = 0; i < g.nx; ++i) { x.push_back(g.xc(i)); val.push_back((1 - w) * f(i, j0) + w * f(i, j0 + 1)); }
}
static void write_profile(const std::string& fn, const std::string& head, const std::vector<double>& a, const std::vector<double>& b) {
    std::ofstream os(fn); os.precision(10); os << "# " << head << "\n";
    for (size_t k = 0; k < a.size(); ++k) os << a[k] << " " << b[k] << "\n";
}

// Build the FlowSolver from the inputs, set U^0, run the initial pressure iterations, then time-step to
// stop_time (or to the steady state) writing diag.csv / VTK snapshots / summary.txt / centreline profiles.
int run_flow(const Params& p) {
    FlowProblem pr = make_flow_problem(p);
    const Geometry& g = pr.geom;
    const std::string out = p.get("out_dir", "out"); mkdir(out.c_str(), 0755);
    const int plot_int = p.get("amr.plot_int", 0), max_step = p.get("max_step", 1000000), diag_int = p.get("diag.interval", 10);
    FlowSolver S(g, pr.vbc, pr.fp);
    const int ng = S.u.ng();
    for (int j = -ng; j < g.ny + ng; ++j) for (int i = -ng; i < g.nx + ng; ++i) pr.vel0(g.xc(i), g.yc(j), S.u(i, j), S.v(i, j));
    const double L = std::min(g.xhi - g.xlo, g.yhi - g.ylo);
    std::cout << "ls2d-baseline (single phase)  problem=" << pr.name << " grid=" << g.nx << "x" << g.ny
              << " rho=" << pr.fp.rho << " mu=" << pr.fp.mu << " g=" << pr.fp.gravity;
    if (pr.fp.mu > 0) std::cout << "  Re=" << pr.fp.rho * 1.0 * L / pr.fp.mu << " (per unit velocity)";
    std::cout << "\n";
    S.initialise();
    auto write_plot = [&](const std::string& fn) {
        Field2D w(g.nx, g.ny, ng), pc(g.nx, g.ny, ng); S.vorticity(w);
        for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) pc(i, j) = 0.25 * (S.p(i, j) + S.p(i + 1, j) + S.p(i, j + 1) + S.p(i + 1, j + 1));
        write_vtk(fn, g, {{"u", &S.u}, {"v", &S.v}, {"p", &pc}, {"vorticity", &w}});
    };
    write_plot(out + "/flow_0000.vtk");
    std::ofstream diag(out + "/diag.csv");
    diag << "step,time,dt,kinetic_energy,max_div_mac,max_div_nodal,mac_iters,visc_iters,nodal_iters,umax,max_dudt,rms_dudt\n";
    auto log = [&](double dt, double dudt, double rms) {
        double umax = 0; for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) umax = std::max(umax, std::max(std::fabs(S.u(i, j)), std::fabs(S.v(i, j))));
        diag << S.step << "," << S.t << "," << dt << "," << S.kinetic_energy() << "," << S.max_divergence_mac() << "," << S.max_divergence_nodal() << ","
             << S.last_mac.iterations << "," << S.last_visc.iterations << "," << S.last_nodal.iterations << "," << umax << "," << dudt << "," << rms << "\n";
    };
    log(0, 0, 0);
    const auto wall0 = std::chrono::steady_clock::now();
    bool steady = false; double dudt = 0, rms = 0; int dudt_i = -1, dudt_j = -1;
    while (S.t < pr.t_final - 1e-12 && S.step < max_step && !steady) {
        double dt = S.estimate_dt(); if (S.t + dt > pr.t_final) dt = pr.t_final - S.t;
        Field2D uold = S.u, vold = S.v;
        S.advance(dt);
        // |dU/dt|: measures the approach to a steady state (lid-driven cavity) and detects blow-up.
        // Both the max (with the cell it occurs in, to tell a localised ringing from a global transient)
        // and the domain rms are reported.
        dudt = 0; rms = 0; double um = 0;
        for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
            const double du = std::fabs(S.u(i, j) - uold(i, j)) / dt, dv = std::fabs(S.v(i, j) - vold(i, j)) / dt;
            if (std::max(du, dv) > dudt) { dudt = std::max(du, dv); dudt_i = i; dudt_j = j; }
            rms += du * du + dv * dv;
            um = std::max(um, std::fabs(S.u(i, j)) + std::fabs(S.v(i, j)));
        }
        rms = std::sqrt(rms * g.cell_area());
        if (!(um < 1e6) || dt < 1e-12) { std::cerr << "BLOW-UP at step " << S.step << " t=" << S.t << " umax=" << um << " dt=" << dt << "\n"; return 2; }
        if (pr.steady_tol > 0 && dudt < pr.steady_tol) steady = true;
        if (S.step % diag_int == 0 || S.t >= pr.t_final - 1e-12 || steady) log(dt, dudt, rms);
        if (plot_int > 0 && S.step % plot_int == 0) { char b[64]; std::snprintf(b, 64, "/flow_%04d.vtk", S.step); write_plot(out + b); }
        if (pr.fp.verbose && S.step % 50 == 0) std::printf("step %d t=%.5f dt=%.3e KE=%.6e divMAC=%.2e divN=%.2e max|dU/dt|=%.3e iters mac/visc/nodal=%d/%d/%d\n",
            S.step, S.t, dt, S.kinetic_energy(), S.max_divergence_mac(), S.max_divergence_nodal(), dudt, S.last_mac.iterations, S.last_visc.iterations, S.last_nodal.iterations);
    }
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
    write_plot(out + "/flow_final.vtk");
    // Centreline profiles: u(y) at the vertical centreline and v(x) at the horizontal one
    // (Ghia, Ghia & Shin 1982, Tables I and II, are tabulated in exactly these two cuts).
    { std::vector<double> a, b;
      column_profile(S.u, g, 0.5 * (g.xlo + g.xhi), a, b); write_profile(out + "/centreline_u.txt", "y  u  at x = " + std::to_string(0.5 * (g.xlo + g.xhi)), a, b);
      row_profile(S.v, g, 0.5 * (g.ylo + g.yhi), a, b);    write_profile(out + "/centreline_v.txt", "x  v  at y = " + std::to_string(0.5 * (g.ylo + g.yhi)), a, b); }
    std::ofstream sum(out + "/summary.txt"); sum.precision(8);
    sum << "problem " << pr.name << "\ngrid " << g.nx << " " << g.ny << "\nsteps " << S.step << "\nt_final " << S.t
        << "\nsteady " << (steady ? 1 : 0) << "\nmax_dudt " << dudt << "\nmax_dudt_cell " << dudt_i << " " << dudt_j
        << "\nrms_dudt " << rms << "\nwall_time_s " << wall
        << "\nkinetic_energy " << S.kinetic_energy()
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
    std::printf("done: %d steps, t=%g, max|dU/dt|=%.3e at cell (%d,%d), rms %.3e%s, wall %.1f s\n",
                S.step, S.t, dudt, dudt_i, dudt_j, rms, steady ? " (steady)" : "", wall);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: ls2d <inputs>\n"; return 1; }
    Params p; p.read(argv[1]);
    p.print(std::cout);
    return run_flow(p);
}
