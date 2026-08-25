// linsolve.cpp — preconditioned conjugate gradient (see linsolve.h for the operators it serves).
#include "linsolve.h"
#include <cmath>

// Interior inner product over ni x nj unknowns.
static double dot(const Field2D& a, const Field2D& b, int ni, int nj) {
    double s = 0; for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) s += a(i, j) * b(i, j); return s;
}
// Subtract the interior mean (projects out the constant null vector).
static void remove_mean(Field2D& f, int ni, int nj) {
    double m = 0; for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) m += f(i, j);
    m /= double(ni) * nj; for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) f(i, j) -= m;
}

// Jacobi-preconditioned variant: wraps diag into a preconditioner and calls the general pcg below.
SolveStats pcg(const std::function<void(Field2D&, Field2D&)>& A, const Field2D& diag, const Field2D& bin, Field2D& x,
               int ni, int nj, double rtol, int maxiter, bool nullspace_constant) {
    Precond M = [&](const Field2D& r, Field2D& z) { for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) z(i, j) = r(i, j) / diag(i, j); };
    return pcg(A, M, bin, x, ni, nj, rtol, maxiter, nullspace_constant);
}

// Standard PCG (Saad, Alg. 9.1) on the interior of x.  Converges when ||r||_2 <= rtol ||b||_2.
// With nullspace_constant the RHS, the initial guess and every preconditioned residual are made mean-free.
SolveStats pcg(const std::function<void(Field2D&, Field2D&)>& A, const Precond& M, const Field2D& bin, Field2D& x,
               int ni, int nj, double rtol, int maxiter, bool nullspace_constant) {
    SolveStats st;
    Field2D b = bin, r(x), z(x), p(x), Ap(x);
    if (nullspace_constant) { remove_mean(b, ni, nj); remove_mean(x, ni, nj); }
    A(x, Ap);
    for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) r(i, j) = b(i, j) - Ap(i, j);
    const double bnorm = std::sqrt(dot(b, b, ni, nj));
    double rnorm = std::sqrt(dot(r, r, ni, nj));
    if (bnorm == 0.0 || rnorm <= rtol * bnorm) { st.converged = true; st.residual = rnorm; return st; }
    z.fill(0); M(r, z); if (nullspace_constant) remove_mean(z, ni, nj);
    p = z;
    double rz = dot(r, z, ni, nj);
    for (int it = 1; it <= maxiter; ++it) {
        A(p, Ap);
        const double alpha = rz / dot(p, Ap, ni, nj);
        for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) { x(i, j) += alpha * p(i, j); r(i, j) -= alpha * Ap(i, j); }
        rnorm = std::sqrt(dot(r, r, ni, nj));
        st.iterations = it; st.residual = rnorm;
        if (rnorm <= rtol * bnorm) { st.converged = true; break; }
        z.fill(0); M(r, z); if (nullspace_constant) remove_mean(z, ni, nj);
        const double rz_new = dot(r, z, ni, nj);
        const double beta = rz_new / rz; rz = rz_new;
        for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) p(i, j) = z(i, j) + beta * p(i, j);
    }
    if (nullspace_constant) remove_mean(x, ni, nj);
    return st;
}
