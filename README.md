# ls2d-baseline — a readable single-core 2D level-set two-phase Navier–Stokes solver

Baseline code for two FrontierPhysics tasks (multiphase flow, fluid–structure interaction).

A from-scratch, dependency-free C++17 implementation of the variable-density approximate projection
method of Almgren et al. (JCP 142, 1998) coupled to the level-set method of Sussman et al. (JCP 148,
1999), with the redistancing / volume constraint of Sussman & Fatemi (SIAM J. Sci. Comput. 20, 1999).
Every routine cites the equation it discretises; `docs/algorithm.md` maps paper equations to functions,
`docs/results.md` compares the code with the literature.

Purpose: a baseline that a beginner can read end-to-end and that is verified against published
results, from which the FrontierPhysics tasks (multiphase flow, fluid–structure interaction) are built.

## Build and run
```
make                               # -> ./ls2d   (g++ -O2 -std=c++17, no libraries)
./ls2d tests/inputs.zalesak        # interface transport with a prescribed velocity
./ls2d tests/inputs.hysing1        # two-phase Navier–Stokes (prob.flow = 1)
```
Outputs go to `out_dir` from the input file: `diag.csv` time series, `*.vtk` snapshots (legacy VTK, open
in ParaView), `summary.txt` with the final numbers. Input files are `key = value` lines; the main keys are
`ns.epsilon` (Heaviside width in cells), `ns.number_of_reinit`, `ns.rho_w / ns.rho_a / ns.mu_w / ns.mu_a`
(phi > 0 / phi < 0 fluids), `ns.gravity`, `ns.sigma`, `ns.lo_bc / ns.hi_bc`, `reinit.*` (redistancing
options, see `src/levelset.h`), `lin.*` (linear solver).

## Layout
```
src/grid.h          Field2D with ghost cells, Geometry, index conventions
src/params.h        key = value input reader
src/bc.h            scalar ghost-cell fill (periodic / zero-gradient)
src/levelset.*      H_eps, delta_eps, sign, rho/mu from phi, redistancing (ENO2 + Heun + volume constraint), diagnostics
src/godunov.*       unsplit Godunov-PLM predictor-corrector: edge states, eq. (31)-(32) advective derivatives
src/advect.*        simpler ENO2/RK2 advection, kept for a side-by-side with the Godunov scheme
src/linsolve.*      preconditioned conjugate gradients
src/mg.*            geometric multigrid V-cycles (cell-centred 5-point, node-centred Q1)
src/flow.*          the projection-method time step: MAC predictor + projection, surface tension,
                    Crank–Nicolson viscous solve, nodal approximate projection, time-step estimate
src/problems*.cpp   initial conditions: disk, zalesak, vortex (prescribed velocity); taylorgreen, hydrostatic, rt, bubble (flow)
src/io.*            VTK / ASCII output
src/main*.cpp       drivers (interface tests / flow tests)
tests/inputs.*      one input file per case in docs/results.md
scripts/            run_table.sh (grid sweeps), plot_contours.py, digitize_gs2009.py (Rayleigh–Taylor reference)
docs/               algorithm.md (equations <-> code), results.md (comparison with the literature)
out/                outputs of the runs reported in docs/results.md (VTK snapshots not tracked; rerun the inputs to regenerate)
references/         validation sources (see references/README.md)
```

## Validation (details and sources in `docs/results.md`)
| case | reference | result |
|---|---|---|
| redistancing of a circle | analytic | zero contour fixed to 2e-5, \|∇φ\| = 1 to 1e-3 |
| Zalesak disk, 1 revolution | Enright 2002 Table 1 | err/L 0.40 (100²), 0.036 (200²) vs 0.61, 0.08 |
| single vortex, T = 8 | Enright 2002 Table 3 | err/L 0.026 (128²), 0.0079 (256²) vs 0.031, 0.008 |
| Taylor–Green | analytic | 2nd order in velocity |
| hydrostatic interface, ρ 1:1000 | analytic | spurious velocity 2e-12 |
| Rayleigh–Taylor, At 0.5, Re 1000 | Guermond–Salgado 2009 Fig. 1 | spike tip within 0.015, bubble tip within 0.04 (128×512) |
| static drop | Laplace law | Δp within 0.8 % |
| rising bubble, Hysing 2009 case 1 | FeatFlow reference | circularity, rise velocity, centroid within 0.3 % (80×160) |
