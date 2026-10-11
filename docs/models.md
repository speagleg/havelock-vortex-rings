# Models

| name | Hamiltonian | parameter | stable side |
|---|---|---|---|
| `mt` | Montaldi–Tokieda family [MT (1.12)], chart radius 1 | `lam_r2` = λr₀² ∈ (−1, 1] | below b_n |
| `sphere`, `plane`, `sinh` | the same family restricted to λr₀² ∈ (0, 1], {0}, (−1, 0) | `lam_r2` | below b_n |
| `tanh` | Kimura's Hamiltonian −Σ ln tanh(d/2) on H² | `xi` = tanh²(r/2) ∈ (0, 1) | above ξ* |

Every threshold also reports the geodesic ring radius `r` (curvature normalised to ±1; the chart radius 1 for the
plane), and on H² also `xi` = tanh²(r/2) and `cosh_r`.

## Montaldi–Tokieda family (`mt`, `sphere`, `plane`, `sinh`)

- Hamiltonian H = −Σ log(|z_i − z_j|²/(σ_iσ_j)), σ = 1 + λ|z|², symplectic form Σ 2/σ² dx∧dy
  (conj13/paper/main.tex, eq. (H)). Fourier blocks a_m = C1 − m(n−m), b_m = m(n−m), g_m = 0 with
  C1 = 2(n−1)(1+b²)/(1+b)² (paper Proposition prop:blocks).
- **Threshold, n ≥ 4: PROVED.** b_n from paper eq. (bn); critical mode ⌊n/2⌋. Lean: `Havelock.Conj13.bn_spec`.
  The `sinh` view cites cosh² r* = q(N) (sinh/note.tex; Lean `Havelock.Sinh`). A view that does not contain b_n
  returns `NoThreshold` ("stable at every radius" or "unstable at every radius"), PROVED.
- **c4 at b_n, n ≥ 4: PROVED.** Paper Theorems thm:even / thm:odd; Lean `Havelock.Conj13.c4_at_bn` checks the algebra
  downstream of the series coefficients. Verdict "stable" (Lyapunov stability modulo rotations, paper Theorem
  thm:main).
- **n = 2, 3 on the sphere / full family: OutOfValidityRange.** Eq. (bn) has its only root at the equator b = 1, and
  c4 = 0 at n = 3; the case rests on Montaldi–Tokieda's own argument.
- **stability(): PROVED** from exact signs of the blocks at the (exact rational) value.
- **simulate(): NUMERICAL** always.

## tanh model (`tanh`)

- Hamiltonian H = −Σ ln tanh(d/2) on H² (Kimura's Hamiltonian; tanh/note_tanh.tex eq. (H)), Poincaré disk, curvature
  −1; parameter `xi` = tanh²(r/2) ∈ (0, 1); stable for xi above the threshold.
- **Blocks: PROVED** (note Thm 2), with a radial–tangential coupling g_m. The Boatto (2008) discrepancy was resolved
  in tanh/boatto_check/RESULT.md (her equations give our thresholds; her figure used a different pair law).
- **Threshold, any N <= 64: CERTIFIED** at run time: Arb root isolation of the irreducible factor of A_⌊N/2⌋, exact
  Sturm counts, exact signs of every block on both sides, and no root of any a_m, b_m on the stable side. The basis
  adds "Lean-checked" (`Havelock.Tanh.cert_N{N}`) for 8 <= N <= 16 and `stable_N{N}` for N <= 7. N <= 7: stable at
  every radius (CERTIFIED). Above N = 64 (`hv.limits()["tanh_exact_max_N"]`) `strict=True` raises
  `CertificationFailed`; `strict=False` returns a NUMERICAL value with `fallback_from = CERTIFIED`.
- **Spectral versus energy–momentum threshold.** The certified threshold ξ* is the *spectral* one: the value above
  which every block has a_m, b_m > 0 and below which the critical mode has a_m b_m < 0. For even N it is also the
  energy–momentum (EM) threshold. For odd N the EM threshold lies strictly above ξ*: on the Krein window
  (ξ*, ξ_EM) the ring is linearly stable (a_m, b_m > 0) but the critical block has determinant D_m = a_m b_m − g_m² < 0,
  so the Hessian is indefinite there; nonlinear stability on that window is undecided by this package. At ξ* itself
  a_m = 0, so the critical block is indefinite (D_m = −g_m² < 0).
- **stability(): CERTIFIED** (exact signs). For odd N the window (ξ*, ξ_EM) gives the verdict
  "spectrally-stable-only".
- **c4: NUMERICAL** (level-set engine), even N only. Odd N raises OutOfValidityRange: at ξ* the critical block is
  indefinite (Krein window, above), so the quartic test does not apply. No closed form is known.
- **simulate(): NUMERICAL** always.
- **Lean export:** `Certificate.to_lean(path)` writes a module importing `Havelock.Tanh.General` and proving
  `Havelock.Tanh.cert_N{N}` (8 <= N <= 16 emitted and tested; the N = 8 export is byte-identical to the committed
  `lean/Havelock/Tanh/CertN8.lean` body and is compiled in CI, `.github/workflows/lean-export.yml`). MT-family
  certificates raise `InvalidInput`: those thresholds are proved in the Lean library itself.
