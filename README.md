# ls2d-baseline (`single-phase`) — a readable single-core 2D incompressible Navier–Stokes solver

Single-phase variant of `ls2d-baseline`: the two-phase level-set machinery of the `main` branch
(level set, redistancing, variable density and viscosity, surface tension, the prescribed-velocity
interface tests) has been removed, leaving the constant-density projection solver on its own — a
compact reference implementation that can be read end-to-end and is verified against published
results, small enough to serve as a starting point for further method development.

A from-scratch, dependency-free C++17 implementation of the approximate projection method of Almgren
et al. (JCP 142, 1998), discretised as in §3 of Sussman et al. (JCP 148, 1999) — that paper's §3
reduces to exactly this scheme when the two fluids are identical, and its equation numbers are the
ones the routine comments cite. `docs/algorithm.md` maps paper equations to functions,
`docs/results.md` compares the code with the literature.

Solved on a uniform, single-level Cartesian grid, ρ and μ constant:

    U_t + (U·∇)U = -(1/ρ) ∇p + (1/ρ) ∇·( μ (∇U + ∇Uᵀ) ) + g e_y ,     ∇·U = 0

## Build and run
```
make                               # -> ./ls2d   (g++ -O2 -std=c++17 -Wall -Wextra, no libraries)
./ls2d tests/inputs.taylorgreen    # decaying vortex, analytic solution
./ls2d tests/inputs.cavity         # lid-driven cavity at Re = 100
python3 scripts/compare_ghia.py 100 out/cavity_re100_96      # compare with Ghia et al. 1982
```
Outputs go to `out_dir`: `diag.csv` time series, `flow_*.vtk` snapshots (legacy VTK, open in
ParaView), `centreline_u.txt` / `centreline_v.txt` (u along the vertical centreline, v along the
horizontal one) and `summary.txt` with the final numbers and the wall time.

## Input file
`key = value` lines, `#` starts a comment, unknown keys are ignored. All of the keys the code reads:

| key | meaning | default |
|---|---|---|
| `prob.type` | `taylorgreen` or `cavity` | `taylorgreen` |
| `prob.lid_velocity` | lid speed of the `cavity` problem | 1.0 |
| `geometry.prob_lo`, `geometry.prob_hi` | domain corners, two numbers each | `0 0`, `1 1` |
| `geometry.is_periodic` | periodicity in x, y | `0 0` |
| `amr.n_cell` | cells in x, y | `64 64` |
| `amr.plot_int` | write a VTK snapshot every N steps (0 = only first/last) | 0 |
| `ns.rho`, `ns.mu` | constant density and dynamic viscosity (`ns.mu = 0` is inviscid) | 1.0, 0.0 |
| `ns.gravity` | acceleration in y (negative = downward) | 0.0 |
| `ns.lo_bc`, `ns.hi_bc` | velocity BC per side: 0 periodic, 2 outflow, 3 symmetry, 4 slip wall, 5 no-slip wall | `0 0` |
| `ns.wall_vel_lo`, `ns.wall_vel_hi` | tangential velocity of a no-slip wall (x-side, y-side) — the cavity lid | `0 0` |
| `ns.cfl` | CFL factor on the time step (0.5 is the validated value; 0.9 is unstable) | 0.5 |
| `ns.fixed_dt` | fixed time step, overrides the CFL estimate | 0 (off) |
| `ns.init_iter` | initial pressure iterations, S99 §3.6 | 3 |
| `ns.steady_tol` | stop when max \|dU/dt\| falls below this | 0 (off) |
| `ns.v` | per-step progress printing | 1 |
| `stop_time`, `max_step` | end of the run | 0, 10⁶ |
| `diag.interval` | write a `diag.csv` row every N steps | 10 |
| `out_dir` | output directory | `out` |
| `lin.rtol`, `lin.maxiter` | PCG tolerance and iteration cap | 1e-10, 5000 |
| `lin.precond` | `mg` (multigrid V-cycle) or `jacobi` | `mg` |
| `lin.mg_nu`, `lin.mg_omega` | multigrid smoother sweeps and Jacobi damping | 2, 0.8 |
| `godunov.fourth_order_slopes`, `godunov.use_transverse` | PLM predictor options | 1, 1 |

## Layout
```
src/grid.h          Field2D with ghost cells, Geometry, index conventions
src/params.h        key = value input reader
src/bc.h            boundary-condition codes and scalar ghost-cell fill (periodic / zero-gradient)
src/godunov.*       unsplit Godunov-PLM predictor-corrector: edge states, eq. (31)-(32) advective derivatives
src/linsolve.*      preconditioned conjugate gradients
src/mg.*            geometric multigrid V-cycles (cell-centred 5-point, node-centred Q1)
src/flow.*          the projection-method time step: MAC predictor + projection, Crank–Nicolson
                    viscous solve, nodal approximate projection, time-step estimate, wall BCs
src/problems_flow.* initial conditions and input parsing: taylorgreen, cavity
src/io.*            legacy VTK output
src/main_flow.*     driver: time loop, diagnostics, centreline profiles, main()
tests/inputs.*      one input file per case in docs/results.md
scripts/            run_table.sh (grid sweeps), compare_ghia.py (cavity vs. Ghia et al. 1982)
docs/               algorithm.md (equations <-> code), results.md (comparison with the literature)
out/                outputs of the runs reported in docs/results.md (VTK snapshots not tracked; rerun the inputs to regenerate)
references/         validation sources (see references/README.md)
```

## Validation (details and sources in `docs/results.md`)
| case | reference | result |
|---|---|---|
| Taylor–Green vortex | analytic | 2nd order in velocity; L2(u) = 4.81e-4 / 9.08e-5 / 2.37e-5 on 32² / 64² / 128² — identical to the two-phase code with equal densities |
| lid-driven cavity, Re = 100 | Ghia, Ghia & Shin 1982, Tables I, II | max deviation 0.0088 on 96² (210 s single-core) |
| lid-driven cavity, Re = 400 | Ghia, Ghia & Shin 1982, Tables I, II | max deviation 0.0032 on 96² (478 s single-core) |
| lid-driven cavity, Re = 1000 | Ghia, Ghia & Shin 1982, Tables I, II | max deviation 0.0106 on 96² (778 s single-core) |
