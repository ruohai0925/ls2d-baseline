// io.h — output: legacy VTK (ParaView) for fields, CSV for time series.
#pragma once
#include "grid.h"
#include <string>
#include <vector>
#include <utility>

// Write cell-centred scalar fields as a STRUCTURED_POINTS legacy VTK file (CELL_DATA).
void write_vtk(const std::string& filename, const Geometry& g, const std::vector<std::pair<std::string, const Field2D*>>& fields);
// Write a cell-centred scalar as a plain text matrix (row j, column i) — for python post-processing.
void write_ascii(const std::string& filename, const Geometry& g, const Field2D& f);
