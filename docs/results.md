# ls2d (single-phase) — validation results against the literature

All numbers here are produced by the code in `src/` (unsplit Godunov-PLM advection, MAC projection +
Crank–Nicolson viscous step + Q1 nodal approximate projection, geometric multigrid, constant ρ and μ).
Every row can be reproduced from the listed input file; `scripts/run_table.sh` runs the grid sequences and
`scripts/compare_ghia.py` produces the cavity comparison tables. Wall times are single-core, `g++ -O2`,
on a machine also running other jobs.

## 1. Taylor–Green vortex (`tests/inputs.taylorgreen`) — analytic solution

u = sin2πx cos2πy e^{−8π²νt}, v = −cos2πx sin2πy e^{−8π²νt}, periodic unit box, ν = 0.01, t = 0.2, CFL 0.5.
Grids 32² and 128² via `scripts/run_table.sh tests/inputs.taylorgreen tg "" 32 64 128`.

| grid | L2(u) error | order | max nodal divergence |
|---|---|---|---|
| 32²  | 4.81e-4 | –    | 7.5e-3 |
| 64²  | 9.08e-5 | 2.40 | 8.9e-4 |
| 128² | 2.37e-5 | 1.94 | 1.7e-4 |

Second-order velocity; the O(h²) residual nodal divergence is the expected signature of the approximate
projection (Almgren et al. 1998).

**Regression against the two-phase code.** These are the numbers of the `main` branch (two-phase solver run
with equal densities and viscosities, `ns.do_phi = 0`). Removing the level-set machinery reproduces them
*bit for bit*: L2(u) = 4.8097373e-4 / 9.0829289e-5 / 2.3725758e-5, L2(v) = 4.7106217e-4 / 8.5033621e-5 /
2.2327919e-5, L2(p) = 9.501992e-4 / 5.7014875e-4 / 1.7821359e-4, max nodal divergence 7.5368108e-3 /
8.9474169e-4 / 1.6681381e-4 and the same 12 / 24 / 48 time steps on 32² / 64² / 128² in both branches. The
strip is therefore exact on this path, not merely accurate to the three digits printed above.

## 2. Lid-driven cavity — Ghia, Ghia & Shin, JCP 48 (1982) 387, Tables I and II

`tests/inputs.cavity` (Re = 100), `tests/inputs.cavity_re400`, `tests/inputs.cavity_re1000`.
Unit square, fluid initially at rest, no-slip on all four walls, the top wall sliding at u = 1;
ρ = 1, μ = 1/Re, so Re = ρUL/μ. Grid 96², CFL 0.5. The run stops when the flow is steady,
max |dU/dt| < `ns.steady_tol` = 1e-4, or at `stop_time`.

Reference: the 129×129 solution of Ghia et al., Table I (u at the 17 tabulated y on the vertical centreline
x = 0.5) and Table II (v at the 17 tabulated x on the horizontal centreline y = 0.5), transcribed in
`references/ghia1982_cavity.txt`. ls2d writes those two cuts to `centreline_u.txt` / `centreline_v.txt`;
`scripts/compare_ghia.py <Re> <run dir>` interpolates them linearly to Ghia's locations and prints the
deviations. The two wall points of each table (u = 0/1 at y = 0/1, v = 0 at x = 0/1) are imposed exactly by
the boundary condition and are excluded, as is the Re = 400, x = 0.9063 entry of Table II, which is a
misprint in the original (see the header of `references/ghia1982_cavity.txt`; ls2d gives −0.3838 there,
between the neighbouring −0.4503 and −0.2301, where the table prints −0.23827).

| Re | grid | steps | t_final | steady | wall time | max \|ls2d − Ghia\| | where |
|---|---|---|---|---|---|---|---|
| 100  | 96² | 2607 | 13.58 | yes, max\|dU/dt\| = 1.0e-4 | 210 s | 0.0088 | v(x = 0.8594) |
| 400  | 96² | 6205 | 32.32 | yes, max\|dU/dt\| = 1.0e-4 | 478 s | 0.0032 | u(y = 0.1016) |
| 1000 | 96² | 8640 | 45.00 | see below | 778 s | 0.0106 | u(y = 0.0703) |

Grid-refinement check at Re = 100 (`tests/inputs.cavity_128`, 128², CFL 0.25 — see the remark on the time
step below): 6950 steps, t = 13.57, steady, 815 s, max deviation again 0.0088 at v(x = 0.8594). The two grids
agree with each other to 4.7e-4 at every one of Ghia's locations (e.g. v(0.8594) = −0.23324 on 96² and −0.23325
on 128²), so the 96² results are grid-converged and the remaining difference from Table I/II is not a
resolution effect.

Away from the listed worst points the agreement is much closer: at Re = 400 every compared value is
within 0.0033 and the whole v-profile within 0.0017; at Re = 1000 every value is within 0.0024
except the four points inside the bottom boundary layer (y ≤ 0.1016), where 96² has about three cells across
the layer and ls2d under-predicts |u| by ~0.010. Full tables:
`python3 scripts/compare_ghia.py 100 out/cavity_re100_96 400 out/cavity_re400_96 1000 out/cavity_re1000_96`.

### Steadiness at Re = 1000

At Re = 1000 the two cells at the upper corners settle into a small limit cycle of period 60 steps, so the
*maximum* of |dU/dt| never falls below the tolerance (it cycles between 0.35 and 1.03, with the maximum
always in cell (94,94) — `max_dudt_cell` in `summary.txt`). The bulk flow is steady: over the last time unit
of the run the kinetic energy changes by 2.5e-4 of itself, and the domain rms of |dU/dt| is 0.015, comparable to
the 0.010 that a single cell of the 9216 oscillating with amplitude 1 would produce on its own. That run therefore ends at
`stop_time` = 45 rather than on the steadiness test. At Re = 100 and
Re = 400 the maximum sits at an interior cell and the test triggers normally.

### Remarks

The lid velocity is discontinuous at the two upper corners, so the flow is not smooth there. Two consequences
are visible in the diagnostics:

* the residual nodal divergence of the approximate projection, which is O(h²) for the Taylor–Green vortex, is
  O(1) in the two corner cells (`diag.csv` column `max_div_nodal`: 1.41 / 1.82 / 2.98 at Re = 100 / 400 /
  1000) and decays away from them;
* the stable time step is smaller than the CFL estimate suggests. CFL 0.5 is stable on 96² at all three
  Reynolds numbers, but on 128² the corner cells at Re = 100 drive a saturated oscillation of the whole
  cavity; CFL 0.25 runs 128² cleanly (max nodal divergence 1.74, monotone decay of max |dU/dt|). CFL 0.9 is
  unstable already on 64². `ns.cfl = 0.5` is the validated setting for the grids reported above.

The moving-wall condition itself is exact: with periodic sides (`geometry.is_periodic = 1 0`,
`ns.lo_bc = 0 5`, `ns.hi_bc = 0 5`, `ns.wall_vel_hi = 0 1`) the same problem is plane Couette flow, and ls2d
reproduces u(y) = y with max |u − y| = 1.0e-8 on 32² — exactly the level at which that run was stopped
(`ns.steady_tol = 1e-8`), so the error is the residual transient and not the discretisation.

The Re = 100 deviation is dominated by the v-profile on the right half of the horizontal centreline, where
ls2d is consistently 3–4 % larger in magnitude than Table II (v(0.8594) = −0.2332 against Ghia's −0.22445).
As the 96²/128² comparison above shows, that value is grid-converged, so the difference is a difference
between two converged solutions, not a resolution effect on this side. Two candidates remain and nothing
measured here settles which dominates: the O(1) discretisation error at the singular upper corners, which at
Re = 100 is spread through the whole cavity by viscous diffusion (the diffusive length over the run,
sqrt(νt) ≈ 0.37, is comparable to the cavity side, whereas at Re = 400 and Re = 1000 it stays local and the
profiles agree to ~0.002 outside the bottom boundary layer); and the accuracy of the 1982 129×129 reference
itself, which is not quantified in the paper.

## 3. Summary

| case | reference | status |
|---|---|---|
| Taylor–Green | analytic | 2nd order in velocity; bit-identical to the two-phase branch with equal fluids |
| lid-driven cavity, Re = 100 | Ghia et al. 1982, Tables I, II | max deviation 0.0088 (96², 210 s) |
| lid-driven cavity, Re = 400 | Ghia et al. 1982, Tables I, II | max deviation 0.0032 (96², 478 s) |
| lid-driven cavity, Re = 1000 | Ghia et al. 1982, Tables I, II | max deviation 0.0106 (96², 778 s) |
| plane Couette flow | analytic | u(y) = y to 1.0e-8, the level the run was stopped at |
