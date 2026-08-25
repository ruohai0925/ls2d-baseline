# references/ — validation sources used by ls2d-baseline (2D)

The Elsevier/Wiley PDFs (S99, A98, Guermond–Salgado 2009, Hysing 2009) are kept locally but not redistributed in
the public repository; DOIs: 10.1006/jcph.1998.6106 (S99), 10.1006/jcph.1998.5890 (A98), 10.1016/j.jcp.2008.12.036,
10.1002/fld.1934. The open-access items (Sussman–Fatemi HAL copy, Enright 2002 author copy, FeatFlow data, one
digitised figure) are included.

| File | Paper | Used for (see `ls2d-baseline/docs/results.md`) |
|---|---|---|
| `An Adaptive Level Set Approach ...pdf` | Sussman, Almgren, Bell, Colella, Howell, Welcome, JCP 148 (1999) | the algorithm itself; interface-error metric eq. (80) |
| `A Conservative Adaptive Projection Method ...pdf` | Almgren et al., JCP 142 (1998) | Godunov / approximate projection details |
| `redistance.pdf` | Sussman & Fatemi, SIAM J. Sci. Comput. 20 (1999) | redistancing with volume constraint; nine-point weights eq. (4.8) (deviation D5, confirmed) |
| `Enright2002_hybrid_particle_level_set.pdf` | Enright, Fedkiw, Ferziger, Mitchell, JCP 183 (2002) | §1: Zalesak disk (Table 1) and single vortex T = 8 (Table 3), "level set" columns; error measure eq. (14) |
| `GuermondSalgado2009_splitting_variable_density.pdf`, `gs2009_p11-11.png` | Guermond & Salgado, JCP 228 (2009) §5.2 | §3.1: Rayleigh–Taylor (Tryggvason setup, Re 1000); Fig. 1 digitised with `ls2d-baseline/scripts/digitize_gs2009.py` |
| `Hysing2009_bubble_benchmark.pdf`, `hysing_featflow_data/` | Hysing et al., IJNMF 60 (2009); reference data from featflow.de | §3.3: rising bubble test case 1 (only TP2D level 7, `data_bench_quantities/c1g1l7.txt`, is kept) |

Not obtainable (paywalled) and therefore not used: Rider & Kothe 1998, LeVeque 1996, Tryggvason 1988, Zalesak 1979 — their setups are reproduced in Enright 2002 / Guermond–Salgado 2009.
