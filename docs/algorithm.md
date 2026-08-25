# Algorithm ↔ code correspondence (2D, single grid, single phase)

Sources: **S99** = Sussman, Almgren, Bell, Colella, Howell, Welcome, JCP 148 (1999); **A98** = Almgren, Bell,
Colella, Howell, Welcome, JCP 142 (1998). Equation numbers refer to the printed papers; every routine named
below carries the same reference in its header comment. S99 is the two-phase paper this branch was stripped
down from — its §3 reduces to exactly the scheme below when the two fluids are identical and σ = 0, so its
equation numbering is kept as the reference for the discretisation.

Equations solved (ρ, μ constant, `ns.rho` / `ns.mu`):

    U_t + (U·∇)U = -(1/ρ) ∇p + (1/ρ) ∇·( μ (∇U + ∇Uᵀ) ) + g e_y ,     ∇·U = 0

Staggering (S99 §3): `u, v` cell-centred; `u^ADV` face-centred (MAC); `p` node-centred.

## A. One time step (S99 §3.1) — `FlowSolver::advance` in `src/flow.cpp`
| Step | Form | Source | ls2d |
|---|---|---|---|
| 1a. Predictor to faces | Taylor extrapolation of U to faces at t^{n+½}; 4th-order limited normal slopes; upwinded transverse terms; forcing term | S99 (20)–(24); A98 §3.2 | `plm_slope` (`src/godunov.cpp`), `FlowSolver::predict_mac_velocity`, `FlowSolver::forcing_for_predictor` |
| 1b. Upwind face velocity | choose left/right state by the sign rule | S99 (25) | `predict_mac_velocity` (`riemann` lambda) |
| 1c. MAC projection | D^MAC(σ G^MAC p) = D^MAC U^{n+½}, σ = 1/ρ; U^ADV = U^{n+½} − σ G^MAC p | S99 (27)–(28) | `FlowSolver::mac_project` |
| 2. Velocity advection term | edge states of u, v upwinded with U^ADV, then [(U·∇)U]^{n+½} | S99 (29)–(32) | `godunov_edge_states`, `godunov_advective_derivative` (`src/godunov.cpp`) |
| 3. Crank–Nicolson viscous solve | (U* − Uⁿ)/Δt = −[(U·∇)U] − Gp^{n−½}/ρ + (L* + Lⁿ)/(2ρ) + F | S99 (16); L in §3.3 | `FlowSolver::viscous_operator`, `FlowSolver::viscous_solve` |
| 4. Approximate nodal projection | V = U*/Δt + Gp^{n−½}/ρ; L_σ p^{n+½} = D V (Q1 finite-element stencil); U^{n+1} = Δt (V − σ G p^{n+½}) | S99 (17), §3.4; A98 §3.3 | `FlowSolver::nodal_project` |
| Lagged pressure gradient | Gp at cell centres from the four surrounding nodes | S99 (33) | `FlowSolver::pressure_gradient_from_nodes` |
| Δt | CFL on \|U\| (and on the wall velocity), plus the gravity constraint | S99 §3.1.1 | `FlowSolver::estimate_dt` |
| Initialisation | project U⁰, iterate the first step to obtain p^{½} | S99 §3.6 | `FlowSolver::initialise` |

Viscous stress: `viscous_operator` forms ∇·(μ∇U) implicitly (Crank–Nicolson) and the transpose part
∇·(μ∇Uᵀ) explicitly. For constant μ the transpose part is μ∇(∇·U), which vanishes for the exact solution;
the discrete velocity is only *approximately* divergence-free, so it is kept as the O(h²) term it is. That
also removes the need for the explicit-viscous time-step limit of S99 §3.1.1.

Linear solves: PCG (`pcg`, `src/linsolve.cpp`) preconditioned by geometric multigrid (`MGCell` for the
cell-centred MAC/viscous systems, `MGNode` for the Q1 nodal system, `src/mg.cpp`). With constant ρ and μ all
three operators have constant coefficients, but the variable-coefficient interfaces are kept unchanged.

## B. Diagnostics — `src/flow.cpp`, `src/main_flow.cpp`
| Quantity | Form | ls2d |
|---|---|---|
| Kinetic energy | ∫ ½ρ\|U\|² dx | `FlowSolver::kinetic_energy` |
| MAC divergence | max \|D^MAC U^ADV\|, at solver tolerance after step 1c | `FlowSolver::max_divergence_mac` |
| Nodal divergence | max \|D U\| at interior nodes (transpose of eq. 33); O(h²) for a smooth flow — this is what makes the projection *approximate* (A98 §3.3) | `FlowSolver::max_divergence_nodal` |
| Steadiness | max \|U^{n+1} − Uⁿ\|/Δt; the run stops when it drops below `ns.steady_tol` | `run_flow` (`src/main_flow.cpp`) |
| Centreline profiles | u(y) at x = ½(x_lo+x_hi) and v(x) at y = ½(y_lo+y_hi), linearly interpolated from the cell centres — the two cuts Ghia et al. tabulate | `column_profile`, `row_profile` (`src/main_flow.cpp`) |
| Taylor–Green errors | L2 and L∞ of U − U_exact; L2 of p − p_exact after removing the mean of each | `run_flow` |

## C. Boundary conditions (`src/bc.h`, `FlowSolver::fill_velocity_ghosts`)
Codes in `ns.lo_bc / ns.hi_bc`: 0 periodic · 2 outflow · 3 symmetry · 4 slip wall · 5 no-slip wall.

Ghost cells reflect about the wall face: the normal component is odd, the tangential component is odd for a
no-slip wall and even for a slip/symmetry wall; outflow is zero-gradient. A no-slip wall may slide
tangentially with a velocity w (`ns.wall_vel_lo / ns.wall_vel_hi`), in which case the tangential ghost is
2w − interior, so that the interpolated wall value is w — this is the lid of the driven cavity. Because that
fill is affine rather than linear, the viscous solve evaluates it once on a zero field and moves the result
to its right-hand side, and the operator it hands to PCG uses the homogeneous fill
(`fill_velocity_ghosts(..., homogeneous = true)`).

Pressure (nodal): homogeneous Neumann at walls (only interior cells contribute to the node divergence, which
is the zero-normal-flux condition), Dirichlet 0 at outflow. Cell-centred scalars: zero-gradient at all
non-periodic boundaries (`fill_scalar_ghosts`).
