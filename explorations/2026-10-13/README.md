# Exploration, 2026-10-13: review checks, screened rings, shell BECs, Rossby strips

Exploratory, NUMERICAL work done with Claude after an adversarial review of the repository. Nothing here is part of
the papers or the `havelock` package. Every script runs with numpy/scipy (and mpmath for `review/`). Run each
script from its own folder (`cd rings && python ring_bands.py`).

## review/ — independent checks of the papers' numbers

| Script | What it checks | Result |
|---|---|---|
| `my_c4.py N` | Level-set quartic c4 for the MT family, from scratch (Cartesian KKT, 60 digits) | Matches the closed forms for n = 4, 5, 7, 8, 9, 12, 13, 15, 17 (relative difference 3e-12 to 6e-8) |
| `slice_c4.py N` | Same quartic on the LINEAR slice ker dJ (no level-set curvature) | Differs from c4: 0.2307 vs 0.25 (n=4), 2.2779 vs 1.6964 (n=8), 2.5781 vs 1.5833 (n=9), 10.9142 vs 6.4904 (n=13). Both positive. The paper should say which of these Montaldi–Tokieda computed. |
| `my_tanh.py N xi...` | tanh-model ring: full Cartesian linearisation and Hessian, no block formulas | Brackets xi*(8) = 0.06735622300…; N=9 Krein window width ~3.2e-8 (paper: 3.22e-8); N=7 definite at all tested radii |

## rings/ — point-vortex rings with other kernels and geometries

| Script | Question | Result |
|---|---|---|
| `ring_bands.py` | (1) Rings with the screened kernel K0(r/lambda) in the plane. (2) Ring of N around a central vortex of strength q. | (1) K0 rings: N <= 5 stable at all radii tested (to 40 lambda); N = 6 stable only for R < 1.288 lambda; N >= 7 unstable at every radius. Screening destabilises; it does not give the tanh model's critical spacing. (2) Log kernel: stable iff q_lo <= q <= (N-1)^2/4 with q_lo = (N^2-8N+8)/16 (even N), (N-1)(N-7)/16 (odd N) — matches these closed forms to 4 digits (believed classical; original reference not checked). At q = 1 the largest stable ring is N = 8 (N = 9 marginal). With K0 screening, N = 8 at q = 1 needs R < 1.597 lambda. |
| `sphere.py`, `scan_sphere.py STAGGER q...` | Shell BEC: ring + mirrored antiring (+ pole charges +-q), zero net charge | Validated on MT sphere thresholds (n = 4, 5, 6). Near the pole the bands reduce to the planar ones. White (PRA 109, 013301) reports a charge-22 pole vortex with 19 ring vortices; the band for q = 22 is N in [11, 22], and N = 19 is linearly stable for ring colatitudes up to ~40 deg. On an ideal sphere, trap rotation only adds a multiple of the conserved momentum, so it cannot change point-vortex linear stability. |

## strips/ — Rossby/edge waves on vorticity strips, roll-up, and Saturn

| Script | Question | Result |
|---|---|---|
| `kuo.py` | Rayleigh–Kuo eigenproblem on plane / rotating sphere / H^2 | Validated: Kelvin annulus growth 0.1690 vs 0.1690 (delta = 0.85, m = 6); Rossby–Haurwitz frequencies exact; H^2 grid-converged to 5 digits |
| `jet.py` | Negative result: point-vortex "jets" (ring + opposite ring) | Fastest mode is always the grid scale m = N/2 (vortex-sheet Kelvin–Helmholtz), so point vortices cannot select a polygon number. (theta0 = 15, d = 6 has no co-rotating strength in the bracket; ignore that row.) |
| `h2strip.py h2` / `h2strip.py plane`, `h2check.py` | A strip rolls up into its fastest m* vortices; does that ring survive? | Plane: thin strips give m* > 7 -> ring unstable -> merging. H^2: rings survive when the roll-up spacing (~6-7 D) exceeds l* = 0.844 curvature radii, i.e. D >~ 0.13, at any radius: e.g. R = 1, D = 0.1 -> 10-gon; R = 2, D = 0.2 -> 17-gon; R = 3, D = 0.2 -> 46-gon (plane at R = 2, D = 0.2 gives m* = 9 -> merges). |
| `saturn.py LAT` | Barotropic Gaussian jet with Saturn's rotation | U = 100 m/s: m* = 6 at 76 deg with width 1.5 deg; m* = 10 at 63 deg with width 2.0 deg. m* scales like circumference / jet width (the 63-deg circle is 1.88x the 76-deg circle). beta stabilises wide jets, more strongly at lower latitude. Pattern speeds come out ~+32 m/s (observed ~0 and ~2.5 m/s): the barotropic model gets drift wrong. |
| `saturn_qg.py` | Same jets with a finite deformation radius Ld (equivalent-barotropic QG) | m stays 6-7 (hexagon) and 10-11 (decagon) for Ld from infinity down to 1500 km. Shrinking Ld slows the pattern: hexagon +34 -> +9 m/s at Ld = 1000 km; decagon +33 -> +12 m/s at 1500 km, and the decagon jet becomes stable at 1000 km. Right trend toward the observed near-stationary patterns, not a match. |

## rollup/ — nonlinear roll-up of a strip, H^2 vs plane

| Script | What it does | Result |
|---|---|---|
| `rollup.py`, `clumps.py`, `plot_rollup.py` (figure: `rollup_R2.png`) | Vortex-blob simulation of a uniform strip (960-1440 blobs), clump count by single-linkage clustering | R = 1, D = 0.1: plane rolls into 9 clumps and merges to 6 by t = 60; H^2 rolls into 11 and holds 9-10. R = 2, D = 0.2: plane rolls into 11 and merges to 5 by t = 47 (then steady); H^2 rolls into 20, partially merges, and holds 10-13 clumps in a ring to t = 100, with about half the circulation stripped into filaments. So H^2 keeps polygons with N > 7 that the plane cannot, but the linear prediction (survival of the full m*-gon, 17 here) is too optimistic: unequal clumps still merge. Single seed, coarse discretisation. |
| `patchring.py r a T N...` | Ring of N equal patches (or point vortices, a = 0) on H^2 at radius r, nonlinear | **Point vortices, r = 1.22** (certified edge between N = 13, r* = 1.178, and N = 14, r* = 1.269): N = 10, 12, 13 hold (radial deviation ~0.004 to t = 300); N = 14, 15, 16 break up by t ~ 5. The tanh thresholds hold for the full nonlinear dynamics. **Patches:** a = 0.06: N = 12, 13 hold, N = 14 breaks (same as points). a = 0.15: with 19 blobs/patch N = 13 broke by t = 40, but with 37 blobs/patch it held to t = 40 (max deviation 0.036) and broke between t = 60 and 120; N = 14 breaks. a = 0.30: even N = 6 breaks, by the patches themselves shearing apart (rms patch size 0.24 -> 1.1), not by the ring mode. A lone patch does not drift (checked), so this is not a regularisation artefact, but a 37-blob cluster is a poor model of a smooth patch. Finite-size effects are therefore NOT settled; this needs contour dynamics on H^2. |

## acoustic/ — the acoustic spacetime of a vortex ring (analogue gravity)

Uniform 2D condensate; lengths in l_e = kappa/(2 pi c) (~1.4 healing lengths), speeds in c.

| Script | Question | Result |
|---|---|---|
| `ergo.py` | When do the ring vortices' individual ergoregions (\|v\| > c) merge into one annulus? | Below R_annulus(N): 1.59 (N=3), 3.82 (6), 5.49 (8), 9.01 (12), 15.46 (19), 20.18 (24); R_annulus/N rises from 0.53 to 0.84 toward the giant-vortex limit 1. At merger the vortex spacing is 2.7-5.3 l_e (4-7 healing lengths), so the point-vortex/uniform-density picture is still usable. Havelock stability (N <= 7 alone) then says which merged geometries last. |
| `plot_ergo.py` (figure: `ergo.png`) | Lab frame vs co-rotating frame | The rotating ring is stationary only along chi = d_t + Omega d_phi. In that frame the lab-frame annulus disappears; what remains is an outer acoustic light cylinder r_LC = [1 + sqrt(1 + 4 Omega Q)]/(2 Omega), Q = q + N, and, with a central charge q, an inner ergoregion r_in = [-1 + sqrt(1 + 4 Omega q)]/(2 Omega). Both closed forms match the grid (6.77 vs 6.78 for N=6, R=3; 21.26 and 9.51 for q=22, N=19, R=15). The ring sits in the stationary band r_in < R < r_LC. |

## Caveats

- All of this is NUMERICAL and exploratory; no specialist has reviewed it.
- The linear survival rule assumes roll-up into exactly m* equal vortices; rollup/ shows partial merging on H^2 anyway.
- The Saturn jet widths are a parameter scan, not measured profiles. The decagon latitude (63 deg) is from press
  coverage; the Science Advances paper (doi:10.1126/sciadv.aee4251) was not read.
