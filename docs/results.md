# Results and their status

The status of each claim is one of these three labels. Where a cell instead names a method (for example
"Computer algebra (sympy)", or "Derived" for a calculation without written remainder bounds), that method is the
support for the claim, and the claim is not labelled PROVED.

| Label | Meaning |
|---|---|
| **PROVED** | A pen-and-paper proof in one of the papers. Where marked, its algebra is also checked in Lean (CI-built). |
| **CERTIFIED** | Established by a rigorous computation: exact rational root isolation (sympy), Arb ball arithmetic, PARI/GP factorisation, or a Lean certificate. Of the research computations, only the Lean and PARI results leave a certificate that can be re-checked without re-running the computation; the others are re-run by `scripts/reproduce.py` (see the Certificates section of [REPRODUCE.md](../REPRODUCE.md)). The `havelock` package's threshold certificates (`havelock cert`) use exact Sturm counts and are re-checked independently by `havelock verify` ([certificates.md](certificates.md)). |
| **NUMERICAL** | High-precision computation with no rigorous error bound. |

No human specialist has reviewed any of this yet; the computer checks listed below are the only verification so far.
One of them, in [`independent_check/`](../independent_check/README.md), is an independent from-scratch implementation
that shares no code, series coefficients or Fourier blocks with the papers. It works from the Hamiltonian
in Cartesian coordinates at 40–50 digits. Its results are NUMERICAL and corroborate the claims marked below; they do
not replace the proofs or certificates.

---

## 1. Montaldi–Tokieda Conjecture 13

Paper: [`conj13/paper/main.pdf`](../conj13/paper/main.pdf).

**Setting.** Montaldi and Tokieda (arXiv:1210.5662) put point vortices on the sphere, the plane and the hyperbolic
plane into one family. The regular ring of n vortices loses stability at an explicit curvature value b_n. Their
Conjecture 13 says that the ring is still stable *at* b_n.

**Their treatment.** Montaldi and Tokieda treat even n by a closed-form symbolic computation (Maple), odd
n = 5, 9, 11 numerically, n = 7 with Maple, and n = 3 by a separate remark.

**What is added here.** The decisive quartic coefficient is computed on the momentum level set for every n ≥ 4.

| Claim | Status |
|---|---|
| A proof of Conjecture 13 for every n ≥ 4, uniform in n; in particular it covers odd n ≥ 13. | PROVED (paper, Theorem 1) |
| Stability is Lyapunov stability modulo rotations, via an energy–momentum argument on the momentum level set. It holds for perturbations of nearby momentum too. | PROVED (paper, §5) |
| Closed form for the quartic coefficient on the level set, even n: c₄ = (n−2)(m⁴+40m²+48m+16)/(192n(n−1)), m = n−2 | PROVED (paper, Theorem 6); Lean-checked from the series coefficients onward, and the series coefficients themselves are derived in Lean (row below). Also reproduced numerically for even n = 8–16 by the [independent check](../independent_check/README.md) (level-set quartic from the Cartesian Hamiltonian, relative difference ≤ 5·10⁻¹⁰) |
| Closed form, odd n ≥ 5: c₄ = (n−1)(n−3)(n²−9n+38)/(128n). Positivity is immediate (discriminant −71). The formula vanishes at n = 3, which is why n = 3 rests on Montaldi–Tokieda's own argument. | PROVED (paper, Theorem 7); the algebra is Lean-checked. Also reproduced numerically for odd n = 9–21 by the [independent check](../independent_check/README.md) (relative difference ≤ 7·10⁻⁹), and at n = 7 in the plane, where it gives c₄ = 9/14 |
| Exact table of c₄ for n = 4–25 | PROVED (Lean, exact over ℚ), plus an independent check modulo three 61-bit primes |
| The series coefficients h₂, h₄, f₂, E_ij that the closed forms start from | Computer algebra (sympy), derived twice independently and stated in the paper. All of them (h₂, h₄ for even n; f₂ and the fourth-order E_ij along the odd ansatz), the J₁ constraint and both normalisations are also **derived from the Hamiltonian in Lean** (`Conj13Series.lean`). The structural steps (symmetric criticality, completeness of the odd ansatz) are not formalised |
| The Fourier blocks and spectrum at the ring: a_m = C1 − m(n−m), b_m = m(n−m), kernel exactly the critical mode, positive on the rest | PROVED (paper, Propositions 3 and 4). The Hessian computation is on paper; Lean checks the arithmetic that follows from the block formulas (the Fejér identity, the kernel mode and the spectral gap) |
| §5 nonlinear-stability argument | PROVED (paper). Its analytic core (Lemmas B1–B4, Theorem 17) is Lean-checked in an abstract setting (`Conj13Stability.lean`). There, Lemma B1 takes as hypotheses the existence of the slaving function β and its Taylor (slaving) bound, which the paper derives, and the hand-off from B3 (chart coordinates) to B4 (metric space) goes through a chart that is also a hypothesis. The instantiation to the vortex ring is not formalised |

## 2. The tanh (Kimura) model on the hyperbolic plane

Note: [`tanh/note_tanh.pdf`](../tanh/note_tanh.pdf).

**Prior work.** Boatto (*Discrete Contin. Dyn. Syst. B* 10, 2008) posed the ring-stability question for this
Hamiltonian, wrote the Hessian in the same circulant block form, and computed the stability ranges numerically for
N = 5–9. The question and the qualitative picture are hers.

**What is added here.**

| Claim | Status |
|---|---|
| In log-polar coordinates the energy–momentum Hessian splits into 2×2 Fourier blocks whose entries are explicit rational functions of ξ = tanh²(r/2) | PROVED (note, Thm 2); **derived from the Hamiltonian in Lean** (`TanhBlocks.lean`: criticality, the Hessian as circulant matrices, its Fourier blocks), and the full 2N-dimensional linearisation's spectrum is the union of the block spectra (`LinearFlow.lean`) |
| Reciprocity: under ξ ↦ 1/ξ, a_m and b_m are invariant and g_m changes sign | PROVED (note, Prop. 5); Lean-checked for a_m, b_m and g_m (`aBlock_recip`, `bBlock_recip`, `gBlock_recip`) |
| N ≤ 7: the ring is stable at every radius | CERTIFIED (exact rational root isolation), and Lean-checked (Bernstein certificates) |
| N = 8–16: the threshold ξ*(N) is the unique root of an explicit integer polynomial f_N in (0, 1), and the critical mode is ⌊N/2⌋ | CERTIFIED (exact root isolation plus Arb enclosure), and Lean-checked (uniqueness, sign change and an interval of width 10⁻²⁰). The table below prints 25 digits; `tanh/certify_v2.out` prints 55, from Arb balls of radius < 3·10⁻¹⁵³ that lie inside exact rational isolating intervals of width < 10⁻⁴⁵; the Lean intervals have width 10⁻²⁰. For N = 8–12 and 16, also reproduced to at least 20 digits by bisection on the full 2N×2N linearisation in the [independent check](../independent_check/README.md) |
| f_N is irreducible and palindromic; degrees 10, 16, 12, 20, 16, 24, 18, 28, 22 for N = 8–16 | CERTIFIED (PARI/GP `polisirreducible`, plus a second check by factorisation patterns mod p; `tanh/irreducibility/`). **Lean-checked** too (`TanhIrred.lean`: mod-p certificates for p = 7, 13, 107, 17, 11, 101, 3, 29, 163, checked by kernel reduction). Irreducibility over ℚ and palindromy also confirmed with sympy by the [independent check](../independent_check/README.md) |
| ξ*(8) is an algebraic unit of degree 10; ξ*(9)–ξ*(16) are not algebraic integers | CERTIFIED. Lean checks irreducibility of f_N, its coefficients (hence its degree and leading coefficient), its unique root in (0, 1), and that ξ*(8) is a unit. That f_N is then the minimal polynomial of ξ*(N) (so the degrees are those of ξ*), and that a leading coefficient ≠ ±1 makes ξ*(N) a non-integer, follow on paper (Gauss's lemma) |
| For odd N, a "Krein window" lies between the spectral threshold and the energy–momentum threshold. In it the ring is linearly (spectrally) stable; nonlinear stability is undecided. Width 3.22×10⁻⁸ at N = 9. | Non-empty for every odd N ≥ 9: PROVED (note, Thm 18(b); row below). For N = 9, 11, 13, 15 the widths are CERTIFIED and non-emptiness is also Lean-checked from the certificates. At N = 9 the [independent check](../independent_check/README.md) finds the reduced Hessian indefinite at ξ* + 3.0·10⁻⁸ and definite at ξ* + 3.4·10⁻⁸, consistent with the width |
| Exact symbol form: with τ = −(N/2) ln ξ and x = 1 − 2m/N, b_m/N² and g_m/N² are explicit functions of (x, τ), and a_m/N² is one plus a term that depends on N and τ only (of order N⁻²) | PROVED (note, Thm 11); checked against the block formulas to 10⁻⁵⁸ (`tanh/large_n.out`). **Lean-checked**, together with the rows below that are marked PROVED, except where a row says otherwise (`lean/Havelock/TanhLargeN.lean`; see [lean.md](lean.md)) |
| b_m > 0 and the critical mode is ⌊N/2⌋, for **every** N and every radius (so the ring is spectrally stable iff a_⌊N/2⌋ ≥ 0) | PROVED (note, Lemmas 12–13, Cor. 14) |
| Exactly one spectral threshold for every N ≥ 8, and stability at every radius for N ≤ 7. Even N: closed form 8a_{N/2}/N² = 2μ tanh η − (coth η − μ)² − 4/N², η = τ/2, μ = (2/N) coth(2η/N) | Even N ≥ 8 and even N ≤ 6: PROVED (note, Prop. 15). Odd N ≥ 17: PROVED (note, Prop. 21). Odd 9 ≤ N ≤ 15 and N = 3, 5, 7: CERTIFIED (note, Thm 7). Lean: `unique_threshold_of_ge_eight` (`AllN.lean`) combines them |
| Large N: cosh r* = κ∞N + O(1/N) for both parities, κ∞ = 1/(2η∞) = 0.1343122142726…, where η∞ is the root of (coth η − 1/η)² = 2 tanh(η)/η. Equivalently, the arc length along the ring between neighbours at threshold tends to ℓ* = π/η∞ = 0.8439 curvature radii | PROVED (note, Thm 16). The limit function is (1/8)[2 tanh(η)/η − L(η)²], L the Langevin function, which is strictly decreasing |
| Second order: cosh r* = κ∞N + C₁/N + O(N⁻³), C₁ = 0.5412188 (even N), 0.4010382 (odd N) | Derived, without written remainder bounds; matches the exact thresholds to 2·10⁻⁶ at N = 800, 801 |
| Energy–momentum threshold: equal to the spectral threshold for **every** even N (no Krein window); for **every** odd N ≥ 9 a non-empty Krein window, of width ξ_EM − ξ* = O(N⁻³); for both parities cosh r_EM = κ∞N + O(1/N). The key inequality: in the symbol variables, the EM determinant is N⁴β[Φ_N(τ) + G(x, τ)] with G ≥ 0 | PROVED (note, Lemma 17 and Thm 18). The limiting window constant N³(ξ_EM − ξ*) → 0.0037 is derived, not proved |
| Physical picture: the large-N limit is an infinite vortex row on a horocycle; the critical mode is the row's pairing mode; the row is spectrally stable iff its spacing is at least ℓ*, with definite Hessian iff the spacing exceeds ℓ* (no Krein window) | PROVED (note, Prop. 19 and Cor. 20, by four classical lattice sums incl. one Poisson summation); confirmed against direct lattice sums to 4·10⁻¹² (`tanh/large_n.out`). Lean checks the lattice sums, the Bloch symbol and the stability criterion stated on that symbol; the reduction of the row's flow to Bloch modes is on paper only |
| An independent implementation of Boatto's (2008) published equations agrees with the thresholds here: ξ*(8) to 10⁻⁴¹, and her definiteness threshold at N = 9 equals ξ_EM(9) | NUMERICAL; `tanh/boatto_check/` |

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
Lean-checked; the step from "ξ* is an algebraic integer" to "c(N) ∈ ℤ" for irrational ξ* is on paper only. The ring is
spectrally stable iff r ≥ r* and its Hessian is definite iff r > r*; at r = r* (N ≥ 8) it is not linearly stable, and
that bifurcation point is the subject of the Conjecture 13 paper. The [independent check](../independent_check/README.md) reproduces cosh² r* = q(N) for N = 8, 9 and 13 by bisection on the full linearisation.

## Open items

- **A specialist reader.** Montaldi and Tokieda are the natural readers for both results.
- **Scope: n = 3.** The odd closed form is derived for n ≥ 5. The case n = 3 is Montaldi–Tokieda's.
- **Not formalised.** Conjecture 13: the Hessian computation behind Propositions 3 and 4; the structural steps of the
  expansion (symmetric criticality, completeness of the odd ansatz); the hypotheses of the abstract §5 argument named
  in the table (the slaving function and its bound, the chart); and the instantiation of §5 to the vortex ring. Tanh:
  that κ₀⁻¹JH is the derivative of the vortex vector field; the Gauss-lemma step from irreducibility to minimal
  polynomial and non-integrality; the reduction of the horocyclic row's flow to Bloch modes. Sinh: the step from
  algebraic integer to c(N) ∈ ℤ. Every other entry above marked Lean-checked is in `lean/` ([lean.md](lean.md)).
