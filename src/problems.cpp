// problems.cpp — set-up of the prescribed-velocity interface tests (disk, zalesak, vortex), see problems.h.
#include "problems.h"
#include <cmath>
#include <stdexcept>

static const double PI = 3.141592653589793238462643383279502884197;

// Signed distance to a solid axis-aligned rectangle, positive inside.
static double sdf_rect_inside(double x, double y, double cx, double cy, double hx, double hy) {
    const double qx = std::fabs(x - cx) - hx, qy = std::fabs(y - cy) - hy;
    const double outside = std::sqrt(std::pow(std::max(qx, 0.0), 2) + std::pow(std::max(qy, 0.0), 2));
    const double inside = std::min(std::max(qx, qy), 0.0);
    return -(outside + inside);
}

// Read geometry / BCs / stop time from the input file and attach phi0, indicator, velocity and exact area
// for the requested prob.type.  Throws on an unknown type.
Problem make_problem(const Params& p) {
    Problem pr; pr.name = p.get("prob.type", "disk");
    auto lo = p.get_vec<double>("geometry.prob_lo", {0, 0}), hi = p.get_vec<double>("geometry.prob_hi", {1, 1});
    auto n  = p.get_vec<int>("amr.n_cell", {64, 64});
    pr.geom.define(n[0], n[1], lo[0], hi[0], lo[1], hi[1]);
    auto per = p.get_vec<int>("geometry.is_periodic", {0, 0});
    pr.geom.periodic_x = per[0]; pr.geom.periodic_y = per[1];
    pr.bc.xlo = pr.bc.xhi = per[0] ? BC_PERIODIC : BC_OUTFLOW;
    pr.bc.ylo = pr.bc.yhi = per[1] ? BC_PERIODIC : BC_OUTFLOW;
    pr.t_final = p.get("stop_time", 0.0);

    if (pr.name == "disk") {
        const double cx = p.get("prob.blob_x", 0.5), cy = p.get("prob.blob_y", 0.5), R = p.get("prob.blob_radius", 0.15);
        // Optional distortion so that redistancing has something to do: phi0 = f(dist) with f(0)=0, f'(0)=slope.
        const double slope = p.get("prob.init_slope", 1.0), cubic = p.get("prob.init_cubic", 0.0);
        pr.phi0 = [=](double x, double y) { const double d = R - std::hypot(x - cx, y - cy); return slope * d + cubic * d * d * d; };
        pr.indicator = [=](double x, double y) { return std::hypot(x - cx, y - cy) < R ? 1 : 0; };
        pr.vel = [](double, double, double, double& u, double& v) { u = 0; v = 0; };
        pr.exact_area = PI * R * R;
    } else if (pr.name == "zalesak") {
        // Enright 2002 §3.1: domain [0,100]^2, disk centre (50,75), R=15, slot width 5, slot length 25,
        // u = (pi/314)(50-y), v = (pi/314)(x-50); one revolution = 628.
        const double cx = 50, cy = 75, R = 15, sw = 5, sl = 25;
        const double rcx = cx, rcy = (cy - R) + 0.5 * sl, rhx = 0.5 * sw, rhy = 0.5 * sl;    // slot rectangle
        pr.phi0 = [=](double x, double y) {
            const double phi_disk = R - std::hypot(x - cx, y - cy);
            const double phi_slot = sdf_rect_inside(x, y, rcx, rcy, rhx, rhy);
            return std::min(phi_disk, -phi_slot);        // CSG difference; exact distance except near re-entrant corners
        };
        pr.indicator = [=](double x, double y) {
            const bool in_disk = std::hypot(x - cx, y - cy) < R;
            const bool in_slot = std::fabs(x - rcx) < rhx && std::fabs(y - rcy) < rhy;
            return (in_disk && !in_slot) ? 1 : 0;
        };
        pr.vel = [](double x, double y, double, double& u, double& v) { u = (PI / 314.0) * (50.0 - y); v = (PI / 314.0) * (x - 50.0); };
        // Exact area: disk minus the part of the slot inside the disk.  The slot bottom is tangent to the
        // circle, so the removed part is  sw*sl - sliver, sliver = int_{-sw/2}^{sw/2} (R - sqrt(R^2 - x^2)) dx.
        auto F = [&](double x) { return 0.5 * x * std::sqrt(R * R - x * x) + 0.5 * R * R * std::asin(x / R); };
        const double sliver = sw * R - (F(0.5 * sw) - F(-0.5 * sw));
        pr.exact_area = PI * R * R - (sw * sl - sliver);                                   // = 582.26 (Enright 2002: 582.2)
    } else if (pr.name == "vortex") {
        // Enright 2002 §3.2: unit box, circle R=0.15 at (0.5,0.75), psi = sin^2(pi x) sin^2(pi y)/pi, reversal cos(pi t/T).
        const double cx = p.get("prob.blob_x", 0.5), cy = p.get("prob.blob_y", 0.75), R = p.get("prob.blob_radius", 0.15);
        const double T = p.get("prob.reverse_period", 8.0);
        pr.phi0 = [=](double x, double y) { return R - std::hypot(x - cx, y - cy); };
        pr.indicator = [=](double x, double y) { return std::hypot(x - cx, y - cy) < R ? 1 : 0; };
        pr.vel = [=](double x, double y, double t, double& u, double& v) {
            const double c = (T > 0) ? std::cos(PI * t / T) : 1.0;
            u = -2.0 * std::pow(std::sin(PI * x), 2) * std::sin(PI * y) * std::cos(PI * y) * c;   // u =  psi_y
            v =  2.0 * std::sin(PI * x) * std::cos(PI * x) * std::pow(std::sin(PI * y), 2) * c;   // v = -psi_x
        };
        pr.exact_area = PI * R * R;
    } else throw std::runtime_error("unknown prob.type " + pr.name);
    return pr;
}
