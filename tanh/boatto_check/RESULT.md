# Boatto (2008) vs our tanh thresholds: result

## Verdict: (b)
- Boatto's **stated equations** are exactly our tanh model with the same criterion:
  - eq. (14): Ψ_K = −(1/2π) ln tanh(√|K| r/2);
  - eq. (17): H = c Σ ln((d−1)/(d+1)), with d = cosh(√|K| D) from eq. (19)/(22);
  - she monitors the minimum non-rotational eigenvalue of the full Hessian of H̃ = H + νM in canonical (q, p).
- Implemented faithfully from those equations, her computation gives **our** thresholds to all digits.
- Her **Fig. 5a/5b on the hyperbolic side (K < 0)** were not produced with eq. (17). They match, point by point, the pair law
  c Σ ln((d−1)(d+1)) = −(1/2π) Σ ln sinh(√|K| D),
  which is eq. (17) with the sign of the H2 = −c Σ ln(1+d) term flipped. The pair interaction is −ln sinh of the **full** distance.
- This is a slip in her numerics for K < 0. Her equations do not imply it.
- Her K > 0 (punctured-sphere, tan) curves agree with her equations.
- The discrepancy is **not** (a), a mapping error: the map K = −r² is correct.
- It is **not** (c): our blocks are confirmed independently once more.
- It is **not** the Montaldi–Tokieda sinh(D/2) model either. The "4ξ*_sinh" closeness was a near-coincidence that misses by 2%, about 10× the digitizing error.

## Evidence

### 1. Map (NOTES.md)
- She fixes the geodesic radius r_o = 1 (p. 362, Remark 2) and varies K. So the curvature −1 radius is r = √(−K), and
  ξ = tanh²(√(−K)/2), K = −(2 artanh √ξ)².
- She adds no Robin or background term on H². M_z is linear in p, so Hess H̃ = Hess H in (q, p) at the critical ring.

### 2. Faithful reimplementation from her equations (`boatto_hessian.py`, `map_and_verify.py`)
- Method: pair function f(φ, p_a, p_b) from eqs. (17) and (22), differentiated symbolically, assembled into the 2N×2N Hessian, at 40 digits. It shares no code with `../blocks.py`.

| N | K* from her eqs. (17)/(22) | ξ = tanh²(√−K*/2) | ours (`blocks.py`) |
|---|---|---|---|
| 8 | −0.28218454857650183188 | 0.06735622300408126227 | 0.06735622300408126227 (diff 1e-41) |
| 9 | −0.48869837218039686475 | 0.11287109543513 | ξ_EM(9) = ξ* + 3.22e-8 (the note's certified Krein-window gap) |
| 10 | −0.74890006954582178601 | 0.16611424152623911808 | 0.16611424152623911808 (diff 1e-42) |

- For odd N her criterion (definiteness) gives ξ_EM, not the spectral ξ*. Our note already distinguishes the two. The gap is 3.22e-8 for N=9, reproduced exactly here.
- **Mode-by-mode congruence** (K = −0.25, N = 8, 9, all m), comparing her Fourier blocks Q_m, P_m, G_m with our b_m, a_m, g_m:
  - Q_m/b_m = 0.159155 = 1/2π for every m;
  - P_m/a_m = 0.134906 for every m, i.e. (1/2π)(ds/dp)²;
  - |G_m|/g_m = √(product of the two ratios) = 0.14653.
  So her (correctly implemented) Hessian is our block Hessian in other coordinates, entry for entry.
- At K = −0.25, N = 8: a_4 = −0.115177, and her own matrix gives 4πP_4 = −0.1953 < 0. The ring is unstable there in her model too.

### 3. What her figure actually shows (`scan_fig5a.py`, `beta_family.py`, `compare_detail.py`, `fig5b_digitize.py`)
- Test family for K < 0: H_β = c Σ [ln(d−1) + β ln(d+1)]. Here β = −1 is her eq. (17) (tanh), β = 0 is the sinh(D/2) model, and β = +1 is sinh(D).

| β (K<0 law) | N=8 crossing | N=9 crossing | median abs. residual vs Fig. 5a (N=6..9, K<0) |
|---|---|---|---|
| −1 (tanh, her eq. 17) | −0.28218 | −0.48870 | 0.405 |
| 0 (sinh(D/2), M–T 1.12) | −0.26202 | −0.43360 | 0.199 |
| 0.75 | −0.25010 | −0.40498 | 0.057 |
| **+1 (sinh D)** | **−0.24656** | **−0.39692** | **0.0077** |
| 1.25 | −0.24321 | −0.38943 | 0.041 |
| digitized Fig. 5a | −0.2465 | −0.3964 | — |

- With β = +1, every clean digitized point of Fig. 5a (N = 5..9, K < 0) agrees to ≤ 0.015 in eigenvalue units. Fig. 5b (N = 6) agrees too:
  - K = −0.75: plotted min 4.27 vs model 4.24 (tanh would give 2.86); plotted max 19.33 vs 19.33 (tanh 17.67);
  - K = −0.5: plotted min 3.20 vs 3.18 (tanh 2.51); plotted max 19.24 vs 19.24 (tanh 18.06).
- For K > 0 her curves equal her own punctured-sphere (tan) law:
  - N=6 crossing: 0.2024 digitized vs 0.202445 computed;
  - Fig. 5b at K = 0.1..0.7: agreement ≤ 0.03.
- Exact values of the competing laws:
  - sinh(D/2): ξ* = 0.0627461 = 8 − 3√7 (N=8) and 0.1010205 (N=9), so 4ξ* = 0.2510 and 0.4041;
  - the figure's crossings are 0.2465 and 0.3964;
  - sinh(D) gives 0.24656 and 0.39692, a match.
- **Plotting artefact (does not move any zero crossing).** Plotted values follow min(λ, 2λ) for the minimum and 2λ for the maximum, in her 4π-normalized units. This produces the kinks exactly at zero in Fig. 5a/5b, and it holds for N = 6..9 and in Fig. 5b to about 1%. The likely mechanism is a list containing both λ_i and 2λ_i; the cause is unknown.
- **Unexplained residual.** Her N=5 curve departs from her own tan law near K ≈ 0.4–0.6, crossing at ≈ 0.449 instead of 0.4711. This is on the sphere side, irrelevant to H². Over-weighting the coupling G does not explain it.

### 4. Could our blocks be wrong (c)? No.
- Our blocks were already checked against a Cartesian brute-force Hessian (`../check_jacobian.py`, `../review/indep_lin.py`).
- This adds a third, structurally different check: Boatto's canonical (q, p) variables, her d-formula, and a full 2N×2N Hessian with the rotation zero removed. Every mode is congruent to (a_m, b_m, g_m) with mode-independent factors, and the crossings agree to 1e-41.

## Do any of our tanh claims change?
- **No.** ξ*(N), ξ_EM(N), the block formulas, the algebraic claims and the K = −r_*² conversion all stand.
- The note's sentence "Her criterion is Hessian definiteness, i.e. our energy–momentum threshold" is correct. For odd N, the K values from her equations are the EM ones: K_EM(9) = −0.4886984, versus −0.4886982 at ξ*.

## What the note should say about Boatto (replace the "unreconciled" text in §motivation and §Limitations)
> Boatto [Boatto] formulated exactly this problem and model on the hyperbolic surface: Green's function −(1/2π) ln tanh(√|K| r/2), Hamiltonian (17) c Σ ln((d−1)/(d+1)), ring of fixed geodesic radius r_o = 1. She also gave the Q/P/G circulant block structure with a coupling block and the Hessian-definiteness criterion, and she observed numerically that N ≥ 7 rings can be stabilized only by sufficiently negative curvature. Credit for the question, the block structure and the qualitative picture is hers. Her hyperbolic blocks are not written out, and her Fig. 5a puts the N = 8, 9 crossings at K ≈ −0.247, −0.396. Implementing her eqs. (17), (19) and (38) as stated gives K = −0.28218 and −0.48870, i.e. our ξ_*(8) and ξ_EM(9) under ξ = tanh²(√−K/2). Her plotted K < 0 curves are reproduced to plotting accuracy (≤ 0.015) by the same computation with the sign of the ln(d+1) term reversed, i.e. by the pair law −ln sinh(D). We therefore attribute the difference to a slip in the numerics behind the figure, not to the model. Her K > 0 (punctured-sphere) curves agree with her equations.

Keep the tone factual. Credit stays as already written.

## Files
- `NOTES.md`: quotes with pages and equations, the map, digitization notes.
- `boatto_hessian.py`: reimplementation from eqs. (11), (17), (22), (29), (38).
- `scan_fig5a.py` → `scan_fig5a.out`: Fig. 5a curves and crossings for the tanh and sinh(D/2) laws.
- `modes.py` → `modes_boatto.out`: per-mode Q, P, G.
- `map_and_verify.py` → `map_and_verify.out`: precise crossings, the map, and congruence with `../blocks.py`.
- `digitize_fig5a.py` → `fig5a_digitized.txt`/`.npy`.
- `n5_check.py` → `n5_check.out`; `n5_coupling.py` → `n5_coupling.out`.
- `beta_family.py` → `beta_family.out`.
- `compare_detail.py` → `compare_beta1_gamma-1.out`.
- `fig5b_check.py` → `fig5b_check.out`; `fig5b_digitize.py` → `fig5b_digitize.out`.
- `dl/`: both PDFs. `boatto2008.txt`, `mparc06-311.txt`: text extracts. `fig/`: renders.
