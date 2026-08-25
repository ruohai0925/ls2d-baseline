// problems_flow.h — initial conditions for the Navier–Stokes tests (docs/results.md §2–3).
//   taylorgreen  : decaying vortex, analytic solution (single phase); p = +1/4(cos4πx+cos4πy)e^{-16π²νt}   §2.1
//   hydrostatic  : flat interface at rest under gravity, density ratio rho_pos/rho_neg                V5
//   rt           : Rayleigh–Taylor, Tryggvason 1988 setup as in Guermond–Salgado 2009 §5.2             V6
//   bubble       : circular drop/bubble — static drop (Laplace) and Hysing 2009 rising bubble           V8
#pragma once
#include "flow.h"
#include "params.h"
#include <functional>
#include <string>

struct FlowProblem {
    std::string name;
    Geometry geom; VelBC vbc; FlowParams fp;
    double t_final = 0;
    std::function<void(double, double, double&, double&)> vel0;    // initial velocity at a point
    std::function<double(double, double)> phi0;                   // initial level set (+ = fluid "pos")
    std::function<void(double, double, double, double&, double&, double&)> exact;   // exact u, v, p (Taylor–Green only)
    bool has_exact = false;
};

FlowProblem make_flow_problem(const Params& p);
