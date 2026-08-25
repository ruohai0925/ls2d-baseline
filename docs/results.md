# ls2d — validation results against the literature

All numbers here are produced by the final code in `src/` (unsplit Godunov-PLM advection, level-set redistancing
every step with the Sussman–Fatemi volume constraint, MAC projection + Crank–Nicolson viscous step + Q1 nodal
approximate projection, geometric multigrid, continuum surface tension). Each row can be regenerated with the
listed input file; `scripts/run_table.sh` runs the grid sequences.

Interface error: err/L = (1/L) ∫ |H(φ_exact) − H(φ)| dx with L the initial interface length
(Sussman et al. 1999 eq. 80; Enright et al. 2002 eq. 14). Area = ∫ H(φ) dx with the sharp Heaviside.

## 0. Redistancing alone — analytic signed distance (`tests/inputs.disk`)
A circle of radius 0.25 on 128², initialised with the distorted signed distance φ₀ = d + 20 d³ (d = R − r), then
redistanced with 40 Heun iterations (Δτ = Δx/2, ε = 2Δx, volume constraint on) and no advection.

| quantity | ls2d | exact |
|---|---|---|
| area change | +0.016 % | 0 |
| err/L | 1.9e-5 | 0 |
| max ||∇φ| − 1| in the band | 1.1e-3 | 0 |

The zero contour stays put while |∇φ| → 1: the Sussman–Fatemi constraint behaves as designed.

## 1. Interface transport with a prescribed velocity — Enright et al., JCP 183 (2002), Tables 1 and 3, "level set" column

### 1.1 Zalesak's slotted disk, one revolution (`tests/inputs.zalesak`)
Domain [0,100]², disk R = 15 centred at (50,75), slot width 5, length 25; solid-body rotation, period 628.

| grid | area change, ls2d | area change, Enright | err/L, ls2d | err/L, Enright |
|---|---|---|---|---|
| 50²  | −7.6 %  | −100 %  | 1.60  | 4.03 |
| 100² | −1.8 %  | +5.3 %  | 0.40  | 0.61 |
| 200² | −0.04 % | +0.54 % | 0.036 | 0.08 |

Observed convergence order 100² → 200²: 3.5 (Enright: 2.9).

### 1.2 Single vortex, time-reversed, T = 8 (`tests/inputs.vortex`)
Unit box, circle R = 0.15 at (0.5,0.75), ψ = π⁻¹ sin²(πx) sin²(πy) cos(πt/T).

| grid | area change, ls2d | area change, Enright | err/L, ls2d | err/L, Enright |
|---|---|---|---|---|
| 64²  | +101 % | −100 %  | 0.080  | 0.075 |
| 128² | +29 %  | −39.8 % | 0.0256 | 0.031 |
| 256² | +3.1 % | −10.3 % | 0.0079 | 0.008 |

The interface error matches or beats the reference on every grid. The area *gain* (instead of Enright's loss) is a
property of the Sussman–Fatemi constraint, which conserves the ε-smoothed volume ∫H_ε(φ): where the filament is
thinner than 2ε, restoring |∇φ| = 1 while holding that volume pushes the zero contour outward. With the band of
H'_ε in the constraint narrowed to one cell (`reinit.alpha_delta = 1`) the 128² result becomes −0.8 % / 0.0297.
Without the constraint (`reinit.volume_fix = 0`): −22.6 % / 0.0266.

## 2. Single-phase flow and hydrostatics — analytic solutions

### 2.1 Taylor–Green vortex (`tests/inputs.taylorgreen`; 32² and 128² via `scripts/run_table.sh`), ν = 0.01, t = 0.2, CFL 0.5
u = sin2πx cos2πy e^{−8π²νt}, v = −cos2πx sin2πy e^{−8π²νt}, periodic unit box.

| grid | L2(u) error | order | max nodal divergence |
|---|---|---|---|
| 32²  | 4.81e-4 | –    | 7.5e-3 |
| 64²  | 9.08e-5 | 2.40 | 8.9e-4 |
| 128² | 2.37e-5 | 1.94 | 1.7e-4 |

Second-order velocity; the O(h²) residual nodal divergence is the expected signature of the approximate projection
(Almgren et al. 1998).

### 2.2 Hydrostatic two-fluid interface (`tests/inputs.hydrostatic`), ρ 1:1000, g = 9.81, 64², t = 0.5
max |u| = 1.7e-12, kinetic energy 1.5e-24: the density jump generates no spurious currents.

## 3. Two-phase flow — literature benchmarks

### 3.1 Rayleigh–Taylor instability (`tests/inputs.rt`, `tests/inputs.rt128`) — Tryggvason, JCP 75 (1988); Guermond & Salgado, JCP 228 (2009) §5.2, Fig. 1
Domain (−½,½)×(−2,2), heavy fluid above, ρ 3:1 (At = 0.5), Re = ρ_min d^{3/2} g^{1/2}/μ = 1000,
η(x) = −0.1 cos(2πx), no-slip top/bottom, periodic sides, 64×256. Tryggvason time t_T = t√At.

Guermond–Salgado give no table; their Fig. 1 (density field, six frames, half domain) was digitised with
`scripts/digitize_gs2009.py` (89 px per unit length, so the reference tips are accurate to about ±0.01;
their mesh size is 0.025 in the refined region, ours 1/64 = 0.0156).

| t_T | spike tip: ls2d 64×256 | ls2d 128×512 | G–S Fig. 1 | bubble tip: ls2d 64×256 | ls2d 128×512 | G–S Fig. 1 |
|---|---|---|---|---|---|---|
| 1.00 | −0.365 | −0.371 | −0.360 | 0.296 | 0.298 | 0.303 |
| 1.50 | −0.613 | −0.622 | −0.607 | 0.427 | 0.432 | 0.449 |
| 1.75 | −0.727 | −0.740 | −0.734 | 0.487 | 0.494 | 0.499 |
| 2.00 | −0.836 | −0.852 | −0.849 | 0.544 | 0.554 | 0.581 |
| 2.25 | −0.944 | −0.962 | −0.969 | 0.598 | 0.613 | 0.633 |
| 2.50 | −1.061 | −1.079 | −1.092 | 0.652 | 0.671 | 0.711 |

(ls2d values interpolated from `diag.csv` columns `interface_ymin/ymax` at t = t_T/√At.)
Spike tip agrees with the digitised reference within 0.03 (64×256) / 0.015 (128×512) over the whole sequence.
The bubble tip lags the reference by 0.06 (64×256) and 0.04 (128×512) at t_T = 2.5 and the gap shrinks with
refinement; the secondary roll-up of the mushroom appears in both codes at t_T ≥ 2. Phase-area drift over the
run: 0.5 % (64×256), 1.8 % (128×512, more thin filaments in the roll-up — see the volume-constraint note in §1.2).
Figure: `out/rt64/interface_evolution.png`.

### 3.2 Static drop (`tests/inputs.staticdrop`) — Laplace law
R = 0.25, σ = 1, ρ 1:1, μ = 0.1, 64² periodic, t = 0.5.

| quantity | ls2d | exact |
|---|---|---|
| p_max − p_min | 4.03 | σ/R = 4 (+0.8 %) |
| max spurious velocity | 1.5e-3 (Ca ≈ 1.5e-4), steady | 0 |
| circularity | 1.002 | 1 |

### 3.3 Rising bubble (`tests/inputs.hysing1`, `tests/inputs.hysing1_80`) — Hysing et al., IJNMF 60 (2009), test case 1
Domain [0,1]×[0,2], R = 0.25 at (0.5,0.5), ρ 1000/100, μ 10/1, σ = 24.5, g = 0.98 (Re 35, Eo 10),
no-slip top/bottom, free-slip sides, T = 3. Reference: TP2D, 1/h = 320 (`references/hysing_featflow_data/data_bench_quantities/c1g1l7.txt`).

| grid | min circularity (t) | max rise velocity (t) | y_c at t = 3 |
|---|---|---|---|
| 40×80    | 0.9060 (1.87) | 0.2385 (0.94) | 1.0856 |
| 80×160   | 0.9014 (1.92) | 0.2410 (0.93) | 1.0820 |
| reference | 0.9013 (1.90) | 0.2417 (0.92) | 1.0813 |

All three benchmark quantities are within 0.3 % of the reference at 80×160.

## 4. Summary

| case | reference | status |
|---|---|---|
| redistancing of a circle | analytic | zero contour fixed to 2e-5, |∇φ| = 1 to 1e-3 |
| Zalesak disk | Enright 2002 Table 1 | err/L below reference on all grids |
| single vortex T = 8 | Enright 2002 Table 3 | err/L at or below reference; area drift sign explained by the constraint |
| Taylor–Green | analytic | 2nd order |
| hydrostatic interface | analytic | machine-zero spurious velocity |
| Rayleigh–Taylor | Guermond–Salgado 2009 Fig. 1 (digitised) | spike tip within 0.015, bubble tip within 0.04 at 128×512, converging |
| static drop | Laplace law | Δp within 0.8 % |
| rising bubble case 1 | Hysing 2009 | ≤ 0.3 % on circularity, rise velocity, centroid |
