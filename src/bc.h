// bc.h — ghost-cell filling for cell-centred scalar fields.
//
// Boundary codes are read from `ns.lo_bc / ns.hi_bc` (one integer per side):
//   0 periodic   1 inflow   2 outflow   3 symmetry   4 slip wall   5 no-slip wall
// For a scalar such as the level-set function every non-periodic boundary is
// treated as zero-gradient (first-order extrapolation), which is what a typical
// FillPatch with FOEXTRAP does for `phi`.  Velocity components get their own
// treatment in the flow solver (flow.cpp).
#pragma once
#include "grid.h"

// Boundary-condition codes.
enum BCType { BC_PERIODIC = 0, BC_INFLOW = 1, BC_OUTFLOW = 2, BC_SYMMETRY = 3, BC_SLIPWALL = 4, BC_NOSLIPWALL = 5 };

// One BCType code per domain side.
struct BCSet { int xlo = 0, xhi = 0, ylo = 0, yhi = 0; };

// Fill all ghost layers of a cell-centred scalar.  Corners are filled by
// applying x then y, which is correct for periodic and extrapolation BCs.
inline void fill_scalar_ghosts(Field2D& f, const BCSet& bc) {
    const int nx = f.nx(), ny = f.ny(), ng = f.ng();
    for (int j = -ng; j < ny + ng; ++j) {
        for (int g = 1; g <= ng; ++g) {
            if (j >= 0 && j < ny) {
                f(-g, j)       = (bc.xlo == BC_PERIODIC) ? f(nx - g, j) : f(0, j);
                f(nx - 1 + g, j) = (bc.xhi == BC_PERIODIC) ? f(g - 1, j)  : f(nx - 1, j);
            }
        }
    }
    for (int i = -ng; i < nx + ng; ++i) {
        for (int g = 1; g <= ng; ++g) {
            f(i, -g)         = (bc.ylo == BC_PERIODIC) ? f(i, ny - g) : f(i, 0);
            f(i, ny - 1 + g) = (bc.yhi == BC_PERIODIC) ? f(i, g - 1)  : f(i, ny - 1);
        }
    }
    // x-ghosts of the y-ghost rows (corners) — needed by the 9-point stencil.
    for (int j = -ng; j < ny + ng; ++j) {
        if (j >= 0 && j < ny) continue;
        for (int g = 1; g <= ng; ++g) {
            f(-g, j)         = (bc.xlo == BC_PERIODIC) ? f(nx - g, j) : f(0, j);
            f(nx - 1 + g, j) = (bc.xhi == BC_PERIODIC) ? f(g - 1, j)  : f(nx - 1, j);
        }
    }
}
