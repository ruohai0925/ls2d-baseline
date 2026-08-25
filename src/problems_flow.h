// problems_flow.h — initial conditions for the Navier–Stokes tests (docs/results.md).
//   taylorgreen : decaying vortex, analytic solution; p = 1/4 (cos 4πx + cos 4πy) e^{-16π²νt}   §1
//   cavity      : lid-driven cavity, Ghia, Ghia & Shin, JCP 48 (1982)                          §2
#pragma once
#include "flow.h"
#include "params.h"
#include <functional>
#include <string>

struct FlowProblem {
    std::string name;
    Geometry geom; VelBC vbc; FlowParams fp;
    double t_final = 0;
    double steady_tol = 0;                                         // stop when max|dU/dt| falls below this (0 = off)
    std::function<void(double, double, double&, double&)> vel0;    // initial velocity at a point
    std::function<void(double, double, double, double&, double&, double&)> exact;   // exact u, v, p (Taylor–Green only)
    bool has_exact = false;
};

FlowProblem make_flow_problem(const Params& p);
