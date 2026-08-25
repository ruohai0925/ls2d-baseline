// problems.h — initial conditions and exact solutions for the interface tests.
//
//  disk      static circle: signed distance, no flow            (V1)
//  zalesak   slotted disk in solid-body rotation, Enright 2002 §3.1 (Zalesak 1979)     (V2)
//  vortex    single vortex (Bell–Colella–Glaz) with LeVeque time reversal cos(pi t/T)   (V3)
// `indicator(x,y)` returns 1 inside fluid 1 (phi>0) for the *exact* initial shape;
// used for the Heaviside-difference error, Sussman 1999 eq. (80).
#pragma once
#include "grid.h"
#include "advect.h"
#include "params.h"
#include <string>
#include <functional>

struct Problem {
    std::string name;
    Geometry geom;
    BCSet bc;
    double t_final = 0;
    double exact_area = 0;                                    // area of the phi>0 region
    std::function<double(double, double)> phi0;               // initial level set (signed distance, + inside)
    std::function<int(double, double)> indicator;             // exact indicator of the initial shape
    VelocityFn vel;                                           // prescribed velocity
};

Problem make_problem(const Params& p);
