// linsolve.h — preconditioned conjugate gradient for the symmetric positive
// (semi-)definite systems of the projection method:
//   * MAC projection      (cell-centred p,  5-point  div(sigma grad))
//   * nodal projection    (node-centred p,  9-point bilinear-FE  div(sigma grad))
//   * Crank–Nicolson viscous solve (cell-centred, 5-point  a I - b div(mu grad))
// The operator is supplied as a functor A(x, y) computing y = A x on the
// interior; the caller is responsible for filling ghost cells of x inside A.
// Jacobi preconditioning; the diagonal is supplied by the caller.
#pragma once
#include "grid.h"
#include <functional>

struct SolveStats { int iterations = 0; double residual = 0; bool converged = false; };

// Solve A x = b to relative tolerance `rtol` (on the 2-norm of the residual), at most `maxiter` iterations.
// `nullspace_constant`: the operator annihilates constants (all-periodic / all-Neumann); the RHS is made
// mean-free and the solution mean is pinned to zero.
SolveStats pcg(const std::function<void(Field2D&, Field2D&)>& A, const Field2D& diag, const Field2D& b, Field2D& x,
               int ni, int nj, double rtol, int maxiter, bool nullspace_constant);
// Same, with an arbitrary SPD preconditioner z = M^{-1} r (e.g. one multigrid V-cycle).
using Precond = std::function<void(const Field2D&, Field2D&)>;
SolveStats pcg(const std::function<void(Field2D&, Field2D&)>& A, const Precond& M, const Field2D& b, Field2D& x,
               int ni, int nj, double rtol, int maxiter, bool nullspace_constant);
