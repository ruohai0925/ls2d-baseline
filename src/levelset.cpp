// levelset.cpp — see levelset.h for the equation map.
#include "levelset.h"
#include <cmath>

static const double PI = 3.141592653589793238462643383279502884197;

// Smoothed Heaviside H_eps(phi), Sussman 1999 eq. (50): 0 for phi < -eps, 1 for phi > eps, smooth between.
double heaviside_eps(double phi, double eps) {
    if (phi >  eps) return 1.0;
    if (phi < -eps) return 0.0;
    return 0.5 * (1.0 + phi / eps + std::sin(PI * phi / eps) / PI);
}
// Smoothed delta function dH_eps/dphi, eq. (68); zero outside |phi| <= eps.
double delta_eps(double phi, double eps) {
    if (std::fabs(phi) > eps) return 0.0;
    return 0.5 * (1.0 + std::cos(PI * phi / eps)) / eps;
}
// Smoothed sign function S_eps(phi) = 2 (H_eps - 1/2), eq. (57).
double sign_eps(double phi, double eps) {
    return 2.0 * (heaviside_eps(phi, eps) - 0.5);
}

// out = vneg + (vpos - vneg) H_eps(phi) on interior and ghost cells, eq. (5)-(6) / (13)-(15) (rho or mu from phi).
void material_from_phi(const Field2D& phi, double vpos, double vneg, double eps, Field2D& out) {
    const int ng = phi.ng();
    for (int j = -ng; j < phi.ny() + ng; ++j)
        for (int i = -ng; i < phi.nx() + ng; ++i)
            out(i, j) = vneg + (vpos - vneg) * heaviside_eps(phi(i, j), eps);
}

// One-sided second-order ENO derivatives of d at (i,j): DL from the left/below, DR from the right/above,
// eq. (58)-(59) with the minmod-limited second difference of eq. (63)-(65).  Uses a 5-point stencil.
void eno2_onesided(const Field2D& d, int i, int j, int dir, double h, double& DL, double& DR) {
    // Central second differences D+D- d at i-1, i, i+1 (eq. 64-65 composed).
    auto v = [&](int s) { return dir == 0 ? d(i + s, j) : d(i, j + s); };
    const double dm  = (v(0)  - v(-1)) / h;                       // D^- d_i
    const double dp  = (v(1)  - v(0))  / h;                       // D^+ d_i
    const double d2m = (v(0)  - 2 * v(-1) + v(-2)) / (h * h);     // D+D- d_{i-1}
    const double d2c = (v(1)  - 2 * v(0)  + v(-1)) / (h * h);     // D+D- d_i
    const double d2p = (v(2)  - 2 * v(1)  + v(0))  / (h * h);     // D+D- d_{i+1}
    DL = dm + 0.5 * h * minmod(d2c, d2m);                          // eq. (58)
    DR = dp - 0.5 * h * minmod(d2c, d2p);                          // eq. (59)
}

// Godunov upwind choice between DL and DR for the redistance equation given the sign S, eq. (60)-(62).
// In the "expansion" case (neither side selected) the paper uses 0; mode == "average" uses (DL+DR)/2 (an alternative used by other codes).
double upwind_gradient(double DL, double DR, double S, const std::string& mode) {
    const double wL = DL * S, wR = DR * S;                         // eq. (60)-(61)
    if (wL > 0.0 && wL + wR > 0.0) return DL;                      // eq. (62), case 1
    if (wR < 0.0 && wL + wR < 0.0) return DR;                      //           case 2
    return (mode == "average") ? 0.5 * (DL + DR) : 0.0;            //           case 3: paper 0 (expansion), or the average
}

// |grad d| at every interior cell, using the frozen-sign upwind choice.
static void grad_magnitude(const Field2D& d, const Field2D& S, const Geometry& g, const std::string& mode, Field2D& G) {
    for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
        double DLx, DRx, DLy, DRy;
        eno2_onesided(d, i, j, 0, g.dx, DLx, DRx);
        eno2_onesided(d, i, j, 1, g.dy, DLy, DRy);
        const double dx_ = upwind_gradient(DLx, DRx, S(i, j), mode);
        const double dy_ = upwind_gradient(DLy, DRy, S(i, j), mode);
        G(i, j) = std::sqrt(dx_ * dx_ + dy_ * dy_);
    }
}

// Coefficient multiplying (1 - |grad d|): frozen S(d0) (paper) or Peng's
// S(d) = d / sqrt(d^2 + (|grad d| * 2 dtau)^2)  (Peng et al. 1999; with dtau = dx/2 this is d/sqrt(d^2 + |grad d|^2 dx^2)).
static double sign_coeff(const ReinitParams& p, double S0, double d, double G, double dtau) {
    if (p.sign == "peng") return d / std::sqrt(d * d + std::pow(G * 2.0 * dtau, 2) + 1e-300);
    return S0;
}

// Redistance phi in place by n_iter Heun steps of  d_tau = S(d0)(1 - |grad d|), eq. (51)-(57), (54)-(55),
// each followed by the Sussman–Fatemi volume constraint eq. (66)-(71).  Returns band statistics.
ReinitStats redistance(Field2D& phi, const Geometry& g, const BCSet& bc, const ReinitParams& p) {
    ReinitStats st;
    const int nx = g.nx, ny = g.ny, ng = phi.ng();
    const double eps_h = p.alpha_heaviside * g.dxmin();
    const double eps_s = p.alpha_sign * g.dxmin();
    const double eps_d = (p.alpha_delta > 0 ? p.alpha_delta : p.alpha_heaviside) * g.dxmin();
    const double dtau = p.dtau_factor * g.dxmin();

    Field2D d0(nx, ny, ng), S(nx, ny, ng), G0(nx, ny, ng), G1(nx, ny, ng), d1(nx, ny, ng), dtil(nx, ny, ng);
    fill_scalar_ghosts(phi, bc);
    d0 = phi;                                                   // step 1: d(x,0) = phi
    for (int j = -ng; j < ny + ng; ++j) for (int i = -ng; i < nx + ng; ++i) S(i, j) = sign_eps(d0(i, j), eps_s);   // eq. (57)
    st.vol_before = volume_heaviside(phi, g, eps_h);

    for (int k = 1; k <= p.n_iter; ++k) {
        // ---- RK2 / Heun, eq. (54)-(55) ------------------------------------
        grad_magnitude(phi, S, g, p.upwind, G0);                // L(d^k)   = S (1 - G0)
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            const double s = sign_coeff(p, S(i, j), phi(i, j), G0(i, j), dtau);
            d1(i, j) = phi(i, j) + dtau * s * (1.0 - G0(i, j));   // eq. (54)
        }
        fill_scalar_ghosts(d1, bc);
        grad_magnitude(d1, S, g, p.upwind, G1);                 // L(d^(1)) = S (1 - G1)
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            const double s0 = sign_coeff(p, S(i, j), phi(i, j), G0(i, j), dtau);
            const double s1 = sign_coeff(p, S(i, j), d1(i, j),  G1(i, j), dtau);
            dtil(i, j) = phi(i, j) + 0.5 * dtau * (s0 * (1.0 - G0(i, j)) + s1 * (1.0 - G1(i, j)));   // eq. (55), Heun form
        }
        // Optional guard (not in the paper): a sign change of d inside the band is undone.
        if (p.sign_guard)
            for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i)
                if (dtil(i, j) * phi(i, j) < 0.0 && std::fabs(phi(i, j)) <= eps_h) { dtil(i, j) = 0.1 * phi(i, j); ++st.sign_guard_hits; }
        fill_scalar_ghosts(dtil, bc);

        // ---- volume constraint, eq. (66)-(71) -----------------------------
        if (p.volume_fix) {
            const double tau = k * dtau;                         // tau^k - tau^0
            Field2D num(nx, ny, ng), den(nx, ny, ng), lam(nx, ny, 0);
            for (int j = -ng; j < ny + ng; ++j) for (int i = -ng; i < nx + ng; ++i) {
                // H'_eps(d^0): the support test uses the *frozen* d^0, as in eq. (68).
                const double test = d0(i, j);
                const double del = (std::fabs(test) <= eps_d) ? 0.5 * (1.0 + std::cos(PI * d0(i, j) / eps_d)) / eps_d : 0.0;
                const double ld = (dtil(i, j) - d0(i, j)) / tau;  // (d~^k - d^0)/(tau^k - tau^0)
                num(i, j) = -del * ld;                             // integrand of the numerator of eq. (71)
                den(i, j) = del * del;                             // integrand of the denominator
            }
            for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
                double n = 0, dd = 0;                              // 9-point quadrature, Sussman–Fatemi eq. (4.8)
                for (int jj = -1; jj <= 1; ++jj) for (int ii = -1; ii <= 1; ++ii) {
                    const double w = (ii == 0 && jj == 0) ? p.vf_weight_centre : p.vf_weight_neigh;
                    n += w * num(i + ii, j + jj); dd += w * den(i + ii, j + jj);
                }
                lam(i, j) = (dd > 0.0) ? n / dd : 0.0;             // eq. (71)
            }
            for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
                const double del0 = delta_eps(d0(i, j), eps_d);
                dtil(i, j) += lam(i, j) * tau * del0;              // eq. (70)
            }
        }
        phi.copy_interior_from(dtil);
        fill_scalar_ghosts(phi, bc);
    }
    st.vol_after = volume_heaviside(phi, g, eps_h);
    double mean; grad_norm_error(phi, g, eps_h, st.max_grad_err_band, mean);
    return st;
}

// Area of the phi > 0 region measured with the smoothed Heaviside, eq. (81).
double volume_heaviside(const Field2D& phi, const Geometry& g, double eps) {
    double v = 0; for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) v += heaviside_eps(phi(i, j), eps);
    return v * g.cell_area();
}

// Interface length  sum delta_eps(phi) |grad phi| dA  with central-difference gradients.
double perimeter_heaviside(const Field2D& phi, const Geometry& g, double eps) {
    double L = 0;
    for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
        const double gx = (phi(i + 1, j) - phi(i - 1, j)) / (2 * g.dx), gy = (phi(i, j + 1) - phi(i, j - 1)) / (2 * g.dy);
        L += delta_eps(phi(i, j), eps) * std::sqrt(gx * gx + gy * gy);
    }
    return L * g.cell_area();
}

// Area of the phi > 0 region by sampling the sign of the bilinearly interpolated phi at nsub x nsub
// sub-cell points per cell (the sharp counterpart of volume_heaviside).
double volume_sharp(const Field2D& phi, const Geometry& g, int nsub) {
    // Bilinear interpolation of phi from cell centres to nsub x nsub sub-cell midpoints, count phi > 0.
    double count = 0;
    for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
        for (int sj = 0; sj < nsub; ++sj) for (int si = 0; si < nsub; ++si) {
            const double fx = (si + 0.5) / nsub - 0.5, fy = (sj + 0.5) / nsub - 0.5;   // offset from cell centre in cells
            const int i0 = fx < 0 ? i - 1 : i, j0 = fy < 0 ? j - 1 : j;
            const double tx = fx < 0 ? fx + 1.0 : fx, ty = fy < 0 ? fy + 1.0 : fy;
            const double v = (1 - tx) * (1 - ty) * phi(i0, j0) + tx * (1 - ty) * phi(i0 + 1, j0)
                           + (1 - tx) * ty * phi(i0, j0 + 1) + tx * ty * phi(i0 + 1, j0 + 1);
            if (v > 0) count += 1.0;
        }
    return count * g.cell_area() / (nsub * nsub);
}

// Max and mean of | |grad phi| - 1 | over the band |phi| < band (central differences); measures how
// far phi is from a signed-distance function after redistancing.
void grad_norm_error(const Field2D& phi, const Geometry& g, double band, double& maxerr, double& meanerr) {
    maxerr = 0; meanerr = 0; int n = 0;
    for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
        if (std::fabs(phi(i, j)) >= band) continue;
        const double gx = (phi(i + 1, j) - phi(i - 1, j)) / (2 * g.dx), gy = (phi(i, j + 1) - phi(i, j - 1)) / (2 * g.dy);
        const double e = std::fabs(std::sqrt(gx * gx + gy * gy) - 1.0);
        maxerr = std::max(maxerr, e); meanerr += e; ++n;
    }
    if (n) meanerr /= n;
}
