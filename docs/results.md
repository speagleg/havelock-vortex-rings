# Results and their status

Every claim in this repository carries one of these labels:

| Label | Meaning |
|---|---|
| **PROVED** | A pen-and-paper proof in one of the papers. Where marked, its algebra is also checked in Lean (CI-built). |
| **CERTIFIED** | Established by a rigorous computation: exact rational root isolation (sympy), Arb ball arithmetic, PARI/GP factorisation, or a Lean certificate. Of the research computations, only the Lean and PARI results leave a certificate that can be re-checked without re-running the computation; the others are re-run by `scripts/reproduce.py` (see the Certificates section of [REPRODUCE.md](../REPRODUCE.md)). The `havelock` package's threshold certificates (`havelock cert`) use exact Sturm counts and are re-checked independently by `havelock verify` ([certificates.md](certificates.md)). |
| **NUMERICAL** | High-precision computation with no rigorous error bound. |

No human specialist has reviewed any of this yet. Independent AI readers and the computer checks listed below are the
only verification so far. One of those checks, in [`independent_check/`](../independent_check/README.md), was written
from scratch and shares no code, series coefficients or Fourier blocks with the papers. It works from the Hamiltonian
in Cartesian coordinates at 40–50 digits. Its results are NUMERICAL and corroborate the claims marked below; they do
not replace the proofs or certificates.

---

## 1. Montaldi–Tokieda Conjecture 13

Paper: [`conj13/paper/main.pdf`](../conj13/paper/main.pdf).

**Setting.** Montaldi and Tokieda (arXiv:1210.5662) put point vortices on the sphere, the plane and the hyperbolic
plane into one family. The regular ring of n vortices loses stability at an explicit curvature value b_n. Their
Conjecture 13 says that the ring is still stable *at* b_n.

**What they proved.**

- Even n: they showed the cubic term vanishes and checked the sign of the quartic coefficient with Maple.
- n = 3, 5, 7, 9, 11: they evaluated the odd-n coefficients numerically.

**What is new here.**

| Claim | Status |
|---|---|
| The first proof for odd n ≥ 13. Our argument is uniform in n, so it also re-proves 4 ≤ n ≤ 11 by a second route. | PROVED (paper, Theorem 1) |
| Stability is Lyapunov stability modulo rotations, via an energy–momentum argument on the momentum level set. It holds for perturbations of nearby momentum too. | PROVED (paper, §5) |
| Closed form for the quartic coefficient on the level set, even n: c₄ = (n−2)(m⁴+40m²+48m+16)/(192n(n−1)), m = n−2 | PROVED; the algebra from the series coefficients onward is Lean-checked. Also reproduced numerically for even n = 8–16 by the [independent check](../independent_check/README.md) (level-set quartic from the Cartesian Hamiltonian, relative difference ≤ 5·10⁻¹⁰) |
| Closed form, odd n ≥ 5: c₄ = (n−1)(n−3)(n²−9n+38)/(128n). Positivity is immediate (discriminant −71). The formula vanishes at n = 3, which is why n = 3 rests on Montaldi–Tokieda's own argument. | PROVED; the algebra is Lean-checked. Also reproduced numerically for n = 9–21 (odd) by the [independent check](../independent_check/README.md) (relative difference ≤ 7·10⁻⁹), and 9/14 in the plane at n = 7 |
| Exact table of c₄ for n = 4–25 | PROVED (Lean, exact over ℚ), plus an independent check modulo three 61-bit primes |
| The series coefficients h₂, h₄, f₂, E_ij that the closed forms start from | Computer algebra (sympy). Derived twice independently and printed in the paper; **not** Lean-checked |
| A comparison with intermediate steps of the published argument (paper §7, notes P1–P3), on our reading of it: their computation evaluates the quartic on a linear slice rather than the level set (6/7 against 9/14 at n = 7), and with the weights of their Lemma 8 the key combination is δ − 3γ²/β. None of this affects their conclusions, which our results confirm. | Our direct recomputation reproduces their printed β, γ, δ to every printed digit at n = 7, 9, 11 (`conj13/review/p2_mt_coeffs.py`). The normalisation-free ratio test is `conj13/review/p2_ratio_test.py` |

## 2. The tanh (Kimura) model on the hyperbolic plane

Note: [`tanh/note_tanh.pdf`](../tanh/note_tanh.pdf).

**Prior work.** Boatto (*Discrete Contin. Dyn. Syst. B* 10, 2008) posed the ring-stability question for this
Hamiltonian, set up the same block structure, and computed stability numerically for N = 5–9. The question and the
qualitative picture are hers. Hwang–Kim (2013) is closed access and unread; it may overlap.

**What is added here.**

| Claim | Status |
|---|---|
| In log-polar coordinates the energy–momentum Hessian splits into 2×2 Fourier blocks whose entries are explicit rational functions of ξ = tanh²(r/2) | PROVED (note, Thm 2); the block identities are Lean-checked |
| Reciprocity: under ξ ↦ 1/ξ, a_m and b_m are invariant and g_m changes sign | PROVED; Lean-checked |
| N ≤ 7: the ring is stable at every radius | CERTIFIED (exact rational root isolation), and Lean-checked (Bernstein certificates) |
| N = 8–16: the threshold ξ*(N) is the unique root of an explicit integer polynomial f_N in (0, 1), and the critical mode is ⌊N/2⌋ | CERTIFIED (exact root isolation plus Arb enclosure, 25 digits), and Lean-checked (uniqueness, sign change and an interval of width 10⁻²⁰). For N = 8–12 and 16, also reproduced to at least 20 digits by bisection on the full 2N×2N linearisation in the [independent check](../independent_check/README.md) |
| f_N is irreducible and palindromic; degrees 10, 16, 12, 20, 16, 24, 18, 28, 22 for N = 8–16 | CERTIFIED (PARI/GP `polisirreducible`, plus a second check by factorisation patterns mod p; `tanh/irreducibility/`). Not in Lean. Irreducibility over ℚ and palindromy also confirmed with sympy by the [independent check](../independent_check/README.md) |
| ξ*(8) is an algebraic unit of degree 10; ξ*(9)–ξ*(16) are not algebraic integers | Unit property: Lean. Degree and non-integrality: PARI |
| For odd N = 9, 11, 13, 15, a "Krein window" lies between the spectral threshold and the energy–momentum threshold. In it the ring is linearly (spectrally) stable; nonlinear stability is undecided. Width 3.22×10⁻⁸ at N = 9. | Non-emptiness for N = 9, 11, 13, 15 is Lean-checked; the widths are CERTIFIED. At N = 9 the [independent check](../independent_check/README.md) finds the reduced Hessian indefinite at ξ* + 3.0·10⁻⁸ and definite at ξ* + 3.4·10⁻⁸, consistent with the width. Non-emptiness for every odd N ≥ 9 is now PROVED (row below) |
| Exact symbol form: with τ = −(N/2) ln ξ and x = 1 − 2m/N, b_m/N² and g_m/N² are explicit functions of (x, τ), and a_m/N² is one plus a term that depends on N and τ only (of order N⁻²) | PROVED (note, Thm 11); checked against the block formulas to 10⁻⁵⁸ (`tanh/large_n.out`) |
| b_m > 0 and the critical mode is ⌊N/2⌋, for **every** N and every radius (so the ring is spectrally stable iff a_⌊N/2⌋ ≥ 0) | PROVED (note, Lemmas 12–13, Cor. 14). Previously certified only for N ≤ 16 |
| Even N: exactly one spectral threshold for every even N ≥ 8, and stability at every radius for even N ≤ 6. Closed form 8a_{N/2}/N² = 2μ tanh η − (coth η − μ)² − 4/N², η = τ/2, μ = (2/N) coth(2η/N) | PROVED (note, Prop. 15). For odd N ≥ 17 uniqueness is NUMERICAL (odd N ≤ 15 is certified) |
| Large N: cosh r* = κ∞N + O(1/N) for both parities, κ∞ = 1/(2η∞) = 0.1343122142726…, where η∞ is the root of (coth η − 1/η)² = 2 tanh(η)/η. Equivalently, the arc length along the ring between neighbours at threshold tends to ℓ* = π/η∞ = 0.8439 curvature radii | PROVED (note, Thm 16). The limit function is (1/8)[2 tanh(η)/η − L(η)²], L the Langevin function, which is strictly decreasing |
| Second order: cosh r* = κ∞N + C₁/N + O(N⁻³), C₁ = 0.5412188 (even N), 0.4010382 (odd N) | Derived, without written remainder bounds; matches the exact thresholds to 2·10⁻⁶ at N = 800, 801. The even value explains the earlier fitted 0.541 |
| Energy–momentum threshold: equal to the spectral threshold for **every** even N (no Krein window); for **every** odd N ≥ 9 a non-empty Krein window, of width ξ_EM − ξ* = O(N⁻³); for both parities cosh r_EM = κ∞N + O(1/N). The key inequality: in the symbol variables, the EM determinant is N⁴β[Φ_N(τ) + G(x, τ)] with G ≥ 0 | PROVED (note, Lemma 17 and Thm 18). Previously certified only for N ≤ 16. The limiting window constant N³(ξ_EM − ξ*) → 0.0037 is derived, not proved |
| Physical picture: the large-N limit is an infinite vortex row on a horocycle; the critical mode is the row's pairing mode; the row is spectrally stable iff its spacing is at least ℓ*, with definite Hessian iff the spacing exceeds ℓ* (no Krein window) | PROVED (note, Prop. 19 and Cor. 20, by four classical lattice sums incl. one Poisson summation); confirmed against direct lattice sums to 4·10⁻¹² (`tanh/large_n.out`) |
| Boatto's equations, implemented as written, reproduce our ξ*(8) to 10⁻⁴¹ and ξ_EM(9) exactly. Her plotted K < 0 curves (Fig. 5a/5b) appear to match the pair law −ln sinh D instead, i.e. her eq. (17) with one sign flipped; her qualitative conclusions are unaffected. | NUMERICAL comparison (digitised figure, ±0.002 in K); `tanh/boatto_check/RESULT.md` |

Thresholds for N = 8–16 (ξ = tanh²(r/2), r the hyperbolic ring radius, curvature −1):

| N | ξ* | cosh r* | degree | leading coefficient |
|---|---|---|---|---|
| 8 | 0.06735622300408126226… | 1.1444414784 | 10 | 1 |
| 9 | 0.11287106323660274043… | 1.2544637167 | 16 | 4 |
| 10 | 0.16611424152623911808… | 1.3984100696 | 12 | 7 |
| 11 | 0.20460117688266929063… | 1.5144618547 | 20 | 10 |
| 12 | 0.24741777680406770917… | 1.6575169309 | 16 | 7 |
| 13 | 0.27986138864390668702… | 1.7772431147 | 24 | 18 |
| 14 | 0.31493896311550662339… | 1.9194478920 | 18 | 23 |
| 15 | 0.34245858351217025660… | 2.0416335000 | 28 | 28 |
| 16 | 0.37168181390604545663… | 2.1831006077 | 22 | 17 |

## 3. The sinh model (expository)

Note: [`sinh/note.pdf`](../sinh/note.pdf). The thresholds are known (Montaldi–Tokieda 2013). The note writes them in
the form cosh² r* = q(N), with q(N) = (N−3)/4 for odd N and (N−2)²/(4(N−1)) for even N. It also classifies when the
threshold is an algebraic unit: c(N) = ξ* + 1/ξ* is an integer exactly for N ∈ {8, 9, 11, 15, 23}. That arithmetic is
Lean-checked. The [independent check](../independent_check/README.md) reproduces cosh² r* = q(N) for N = 8, 9 and 13 by bisection on the full linearisation.

## Open items

- **A specialist reader.** Montaldi and Tokieda are the natural readers for both results.
- **Hwang–Kim (2013).** It needs journal access.
- **Odd-n, n ≥ 5 limit.** The odd closed form is derived for n ≥ 5. The case n = 3 relies on Montaldi–Tokieda.
- **Not formalised.** The derivation of the series coefficients (Conjecture 13) and of the Fourier blocks (tanh) are
  computer algebra. The dynamics are not formalised in either case.
