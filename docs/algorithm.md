# Algorithm ↔ code correspondence (2D, single grid)

Sources: **S99** = Sussman, Almgren, Bell, Colella, Howell, Welcome, JCP 148 (1999); **A98** = Almgren, Bell,
Colella, Howell, Welcome, JCP 142 (1998); **SF99** = Sussman & Fatemi, SIAM J. Sci. Comput. 20 (1999).
Equation numbers refer to the printed papers. Every routine named below carries the same reference in its header comment.

Staggering (S99 §3): `u, v, φ, ρ, μ` cell-centred; `u^ADV` face-centred (MAC); `p` node-centred.

## A. State and material properties (`src/levelset.*`, `src/flow.cpp`)
| Step | Form | Source | ls2d |
|---|---|---|---|
| Smoothed Heaviside | H_ε(φ) = 0 / ½[1 + φ/ε + sin(πφ/ε)/π] / 1, ε = α Δx (α = `ns.epsilon`, default 2) | S99 (50) | `heaviside_eps` |
| Density, viscosity | ρ = ρ₂ + (ρ₁ − ρ₂) H_ε(φ), same for μ | S99 (5)–(6) | `material_from_phi`, `FlowSolver::set_material_from_phi` |
| Mid-time φ, ρ, μ | φ^{n+½} = ½(φⁿ + φ^{n+1}) → ρ^{n+½}, μ^{n+½} | S99 (13)–(15) | `FlowSolver::advance` |

## B. One time step (S99 §3.1) — `FlowSolver::advance` in `src/flow.cpp`
| Step | Form | Source | ls2d |
|---|---|---|---|
| 1a. Predictor to faces | Taylor extrapolation of U (and φ) to faces at t^{n+½}; 4th-order limited normal slopes; upwinded transverse terms; forcing term | S99 (20)–(24); A98 §3.2 | `plm_slope`, `godunov_edge_states` (`src/godunov.cpp`) |
| 1b. Upwind face velocity | choose left/right state by the sign rule | S99 (25) | `godunov_edge_states` |
| 1c. MAC projection | D^MAC(1/ρⁿ G^MAC p) = D^MAC U^{n+½}; U^ADV = U^{n+½} − (1/ρⁿ) G^MAC p | S99 (27)–(28) | `FlowSolver::predict_mac_velocity`, `FlowSolver::mac_project` |
| 2. Level-set advection | φ^{n+1} = φⁿ − Δt [U^ADV·∇φ]^{n+½}, non-conservative | S99 (12), (29)–(30) | `godunov_advective_derivative` (`src/godunov.cpp`) |
| 3. Velocity advection term | [(U·∇)U]^{n+½} from face states and U^ADV | S99 (31)–(32) | `godunov_advective_derivative` |
| 4. Surface tension | M = σ κ ∇H_ε, κ from nodal normals | S99 §3.3 (34)–(39) | `FlowSolver::surface_tension` |
| 5. Crank–Nicolson viscous solve | (U* − Uⁿ)/Δt = −[(U·∇)U] − Gp^{n−½}/ρ^{n+½} + (L* + Lⁿ)/(2ρ^{n+½}) + F | S99 (16); L in §3.3 | `FlowSolver::viscous_operator`, `FlowSolver::viscous_solve` |
| 6. Approximate nodal projection | V = U*/Δt + Gp^{n−½}/ρ^{n+½}; L_ρ p^{n+½} = D V (Q1 finite-element stencil); U^{n+1} = Δt (V − Gp^{n+½}/ρ^{n+½}) | S99 (17), §3.4; A98 §3.3 | `FlowSolver::nodal_project` |
| 7. Redistance | see C; applied right after step 2 (switch `reinit_after_projection = 1` restores the S99 order) | S99 §3.5 | `redistance` |
| 8. Δt | min of CFL, gravity, per-phase viscous and capillary constraints | S99 §3.1.1 | `FlowSolver::estimate_dt` |
| Initialisation | project U⁰, iterate the first step to obtain p^{½} | S99 §3.6 | `FlowSolver::initialise` |

Linear solves: PCG (`pcg`, `src/linsolve.cpp`) preconditioned by geometric multigrid (`MGCell` for the cell-centred
MAC/viscous systems, `MGNode` for the Q1 nodal system, `src/mg.cpp`).

## C. Redistancing (S99 §3.5, SF99 §4) — `redistance` in `src/levelset.cpp`
| Step | Form | Source | ls2d |
|---|---|---|---|
| Init | d⁰ = φ^{n+1} | S99 step 1 | `redistance` |
| Sign function | S = 2(H_ε(d⁰) − ½), frozen during the iteration (`reinit.sign = peng` recomputes it from the current d) | S99 (57) | `sign_eps` |
| Pseudo-time step | Δτ = Δx/2; N = `ns.number_of_reinit` iterations (default 4 → τ = ε) | S99 after (55) | `ReinitParams` |
| One-sided ENO2 derivatives | D^L = D⁻d + (Δx/2) m(D⁺D⁻d_i, D⁺D⁻d_{i−1}); D^R = D⁺d − (Δx/2) m(D⁺D⁻d_i, D⁺D⁻d_{i+1}); m = minmod | S99 (58)–(59), (63)–(65) | `eno2_onesided` |
| Godunov upwind choice | w^L = S D^L, w^R = S D^R; D^L if w^L > 0 and w^L + w^R > 0; D^R if w^R < 0 and w^L + w^R < 0; else 0 | S99 (60)–(62) | `upwind_gradient` |
| Spatial operator | L(d) = S (1 − |∇d|) | S99 (56) | `redistance` |
| Heun / RK2 | d^{(1)} = d^k + Δτ L(d^k); d^{k+1} = d^k + ½Δτ [L(d^k) + L(d^{(1)})] | S99 (54)–(55) | `redistance` |
| Volume constraint | λ = −∫H'_ε(d⁰)(d̃^k − d⁰)/(τ^k − τ⁰) / ∫H'_ε(d⁰)²; d^k = d̃^k + λ (τ^k − τ⁰) H'_ε(d⁰), per cell | S99 (66)–(71); SF99 (4.7) | `redistance` (`volume_fix` block) |
| Cell integrals | ∫_Ω g ≈ h²/24 (16 g_ij + Σ 8 neighbours) | SF99 (4.8) | `vf_weight_centre = 16`, `vf_weight_neigh = 1` |
| H'_ε | ½[1 + cos(πd/ε)]/ε for |d| ≤ ε, else 0; width `reinit.alpha_delta` Δx (default = α) | S99 (68) | `delta_eps` |

## D. Diagnostics (S99 §5.1) — `src/levelset.cpp`, `src/main.cpp`, `src/main_flow.cpp`
| Quantity | Form | Source | ls2d |
|---|---|---|---|
| Interface error | E = ∫|H(φ_exact) − H(φ)| dx, midpoint rule on sub-cells with bilinear φ; reported as E/L | S99 (80); Enright 2002 (14) | `interface_error` (main.cpp) |
| Phase area | ∫H(φ) dx (sharp, sub-cell sampled) and ∫H_ε(φ) dx | S99 (81) | `volume_sharp`, `volume_heaviside` |
| Distance-function quality | max / mean of ||∇φ| − 1| in the band |φ| < ε | — | `grad_norm_error` |
| Bubble quantities | centroid y_c, rise velocity v_c, circularity | Hysing 2009 (eq. 1–3) | `run_flow` (main_flow.cpp) |

## E. Boundary conditions (`src/bc.h`, `FlowSolver::fill_velocity_ghosts`)
Codes in `ns.lo_bc / ns.hi_bc`: 0 periodic · 2 outflow · 3 symmetry · 4 slip wall · 5 no-slip wall.
φ: zero-gradient at all non-periodic walls. p (nodal): homogeneous Neumann at walls, Dirichlet 0 at outflow.
