// grid.h — uniform 2D Cartesian grid and a cell-centred field with ghost cells.
//
// Index convention (AMReX-style):
//   cell (i,j), i = 0..nx-1, j = 0..ny-1, centre at x_i = xlo + (i+0.5)*dx.
//   Ghost cells: i = -ng..-1 and nx..nx+ng-1 (same in j).
//   Node (i,j) sits at the lower-left corner of cell (i,j), x = xlo + i*dx.
//   x-face (i,j) is the left face of cell (i,j); y-face (i,j) the bottom face.
#pragma once
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>

// Uniform grid geometry: cell counts, physical extent, spacings and periodicity flags.
struct Geometry {
    int nx = 0, ny = 0;          // number of cells
    double xlo = 0, ylo = 0;     // physical lower corner
    double xhi = 1, yhi = 1;     // physical upper corner
    double dx = 0, dy = 0;
    bool periodic_x = false, periodic_y = false;

    // Set the extent and cell counts; dx, dy follow.
    void define(int nx_, int ny_, double xlo_, double xhi_, double ylo_, double yhi_) {
        nx = nx_; ny = ny_; xlo = xlo_; xhi = xhi_; ylo = ylo_; yhi = yhi_;
        dx = (xhi - xlo) / nx; dy = (yhi - ylo) / ny;
    }
    double xc(int i) const { return xlo + (i + 0.5) * dx; }   // cell centre
    double yc(int j) const { return ylo + (j + 0.5) * dy; }
    double xn(int i) const { return xlo + i * dx; }           // node / x-face
    double yn(int j) const { return ylo + j * dy; }
    double dxmin() const { return std::min(dx, dy); }
    double cell_area() const { return dx * dy; }
};

// A scalar field stored with `ng` ghost layers on each side.  Data are stored
// row-major with j slowest, so the innermost loop should run over i.
// The same class is used for cell-, face- and node-centred data; the caller
// keeps track of the centring (see comments at each use).
class Field2D {
public:
    Field2D() = default;
    // Allocate nx x ny interior cells plus ng ghost layers, all set to `init`.
    Field2D(int nx, int ny, int ng, double init = 0.0) { define(nx, ny, ng, init); }
    void define(int nx, int ny, int ng, double init = 0.0) {
        nx_ = nx; ny_ = ny; ng_ = ng;
        sx_ = nx + 2 * ng; sy_ = ny + 2 * ng;
        data_.assign(static_cast<size_t>(sx_) * sy_, init);
    }
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int ng() const { return ng_; }
    // Element access; ghost indices (-ng..-1, nx..nx+ng-1) are allowed.
    double& operator()(int i, int j) {
        assert(i >= -ng_ && i < nx_ + ng_ && j >= -ng_ && j < ny_ + ng_);
        return data_[static_cast<size_t>(j + ng_) * sx_ + (i + ng_)];
    }
    double operator()(int i, int j) const {
        assert(i >= -ng_ && i < nx_ + ng_ && j >= -ng_ && j < ny_ + ng_);
        return data_[static_cast<size_t>(j + ng_) * sx_ + (i + ng_)];
    }
    // Set every entry (interior and ghosts) to v.
    void fill(double v) { std::fill(data_.begin(), data_.end(), v); }
    // Copy only the interior (valid) cells.
    void copy_interior_from(const Field2D& o) {
        for (int j = 0; j < ny_; ++j) for (int i = 0; i < nx_; ++i) (*this)(i, j) = o(i, j);
    }
    // Max |f| over the interior cells only (ghosts excluded).
    double max_abs_interior() const {
        double m = 0; for (int j = 0; j < ny_; ++j) for (int i = 0; i < nx_; ++i) m = std::max(m, std::fabs((*this)(i, j)));
        return m;
    }
private:
    int nx_ = 0, ny_ = 0, ng_ = 0, sx_ = 0, sy_ = 0;
    std::vector<double> data_;
};

// minmod limiter: the argument of smaller magnitude.  Sussman 1999 eq. (63).
inline double minmod(double a, double b) {
    return (std::fabs(a) <= std::fabs(b)) ? a : b;
}
