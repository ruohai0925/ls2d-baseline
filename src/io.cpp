// io.cpp — legacy VTK and plain-text field writers (see io.h).
#include "io.h"
#include <fstream>
#include <stdexcept>

// Legacy VTK STRUCTURED_POINTS file with one CELL_DATA scalar per (name, field) pair; ghosts are not written.
void write_vtk(const std::string& filename, const Geometry& g, const std::vector<std::pair<std::string, const Field2D*>>& fields) {
    std::ofstream os(filename); if (!os) throw std::runtime_error("cannot write " + filename);
    os << "# vtk DataFile Version 3.0\nls2d\nASCII\nDATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << g.nx + 1 << " " << g.ny + 1 << " 1\n";
    os << "ORIGIN " << g.xlo << " " << g.ylo << " 0\nSPACING " << g.dx << " " << g.dy << " 1\n";
    os << "CELL_DATA " << static_cast<long>(g.nx) * g.ny << "\n";
    for (auto& f : fields) {
        os << "SCALARS " << f.first << " double 1\nLOOKUP_TABLE default\n";
        for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) os << (*f.second)(i, j) << "\n";
    }
}

// Plain-text matrix of the interior cells, one row per j, 15 significant digits.
void write_ascii(const std::string& filename, const Geometry& g, const Field2D& f) {
    std::ofstream os(filename); os.precision(15);
    for (int j = 0; j < g.ny; ++j) { for (int i = 0; i < g.nx; ++i) os << f(i, j) << (i + 1 < g.nx ? " " : ""); os << "\n"; }
}
