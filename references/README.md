# references/ — validation sources used by ls2d-baseline, `single-phase` branch

Two Elsevier PDFs are kept locally but not redistributed in the public repository;
DOIs: 10.1006/jcph.1998.6106 (S99), 10.1006/jcph.1998.5890 (A98).
The lid-driven-cavity benchmark data are transcribed into a plain-text file that *is* tracked.

| File | Paper | Used for (see `docs/results.md`) |
|---|---|---|
| `An Adaptive Level Set Approach ...pdf` (not tracked) | Sussman, Almgren, Bell, Colella, Howell, Welcome, JCP 148 (1999) — "S99" | the equation numbering the code comments cite: predictor eq. (20)-(25), MAC projection eq. (27)-(28), advective derivatives eq. (31)-(32), lagged pressure gradient eq. (33), viscous step eq. (16), nodal projection eq. (17), time step §3.1.1. Its §3 reduces to the present single-phase scheme when the two fluids are identical |
| `A Conservative Adaptive Projection Method ...pdf` (not tracked) | Almgren, Bell, Colella, Howell, Welcome, JCP 142 (1998) — "A98" | the approximate-projection method itself: Godunov-PLM predictor §3.2, Q1 nodal approximate projection §3.3 |
| `ghia1982_cavity.txt` | Ghia, Ghia & Shin, J. Comput. Phys. 48 (1982) 387-411, DOI 10.1016/0021-9991(82)90058-4 | §2: lid-driven cavity at Re = 100, 400, 1000 — Table I (u on the vertical centreline) and Table II (v on the horizontal centreline), 129×129 reference solution. The file header lists the sources the values were taken from and cross-checked against, and documents one known misprint in the original |

The Taylor–Green vortex (`docs/results.md` §1) is compared with its analytic solution and needs no reference file.
