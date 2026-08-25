// io.h — field output: legacy VTK (ParaView).  The time series (diag.csv) and the
// centreline profiles are written directly by the driver, src/main_flow.cpp.
#pragma once
#include "grid.h"
#include <string>
#include <vector>
#include <utility>

// Write cell-centred scalar fields as a STRUCTURED_POINTS legacy VTK file (CELL_DATA).
void write_vtk(const std::string& filename, const Geometry& g, const std::vector<std::pair<std::string, const Field2D*>>& fields);
