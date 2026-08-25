# ls2d-baseline — a readable single-core 2D level-set two-phase Navier–Stokes solver

A from-scratch, dependency-free C++17 implementation of the variable-density approximate projection
method of Almgren et al. (JCP 142, 1998) coupled to the level-set method of Sussman et al. (JCP 148,
1999), with the redistancing / volume constraint of Sussman & Fatemi (SIAM J. Sci. Comput. 20, 1999).
Every routine cites the equation it discretises; `docs/algorithm.md` maps paper equations to functions,
`docs/results.md` compares the code with the literature.

Purpose: a compact reference implementation that can be read end-to-end and that is verified against
published results, small enough to serve as a starting point for further method development.

## What is in this repository, and on which branch

The repository holds one solver in two configurations. They share the same projection-method core —
the Godunov predictor, the MAC and nodal projections, the Crank–Nicolson viscous solve and the
geometric multigrid — and differ only in whether an interface between two fluids exists at all.
Pick the branch that matches the problem you are starting from; each is validated on its own and
neither is a work-in-progress version of the other.

| branch | what it is | validated against |
|---|---|---|
| `main` | the **two-phase** solver: a level set with Sussman–Fatemi redistancing, density and viscosity varying across the interface, surface tension, and the prescribed-velocity interface tests | Zalesak disk and single vortex (Enright et al. 2002), Taylor–Green, hydrostatic balance at a density ratio of 1000, Rayleigh–Taylor (Guermond & Salgado 2009), Laplace's law, and the rising-bubble benchmark of Hysing et al. (2009) |
| `single-phase` | the **constant-density** solver: the level set, redistancing, variable material properties, surface tension and the interface tests are removed, leaving the projection method on its own (`src/` goes from 2032 to 1427 lines). Adds a lid-driven cavity problem | Taylor–Green, reproducing `main`'s errors bit for bit, which is the regression showing the strip changed nothing; and the lid-driven cavity at Re = 100, 400 and 1000 against Ghia, Ghia & Shin (1982) |

Use `main` if the problem has an interface: a rising bubble, a dam break, a droplet, a wave. Use
`single-phase` if it does not, or if you want the smallest readable starting point — a
constant-property incompressible solver with no interface machinery to work around.

The two branches are kept in step: a fix to the shared core is applied to both, and the Taylor–Green
regression above is what checks that they have not drifted apart.

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
