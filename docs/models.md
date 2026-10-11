# Models

| name | Hamiltonian | parameter | stable side |
|---|---|---|---|
| `mt` | Montaldi–Tokieda family [MT (1.12)], chart radius 1 | `lam_r2` = λr₀² ∈ (−1, 1] | below b_n |
| `sphere`, `plane`, `sinh` | the same family restricted to λr₀² ∈ (0, 1], {0}, (−1, 0) | `lam_r2` | below b_n |

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
