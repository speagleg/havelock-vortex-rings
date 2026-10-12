# Sinh note: second-reader review, core claims C1–C3

Reviewer: independent AI reader (fresh context), 2026-10-07. I re-derived everything from the note
(`sinh/note.tex`) before reading the author's log, the literature notes (`literature.md`) or the author's scripts.
File paths below are given in the current repository layout; the Lean file reviewed is now
`lean/Havelock/SinhArithmetic.lean`.
Scope: note §1–§4 (C1, C2, C3), the Lean file, the literature/attribution claims. C4 and C5 are not in scope.

## Verdict: **ACCEPT** (core C1–C3), with minor non-blocking recommendations

No mathematical error found. All formulas re-derived by hand and confirmed by an independent
computation that uses no finite differences and no Fourier reduction. The literature claims (C1/C2 KNOWN, Montaldi–Tokieda 2013)
check out against the primary text. The attribution warning (Kimura's H² Hamiltonian is the tanh one) is
supported by three more secondary sources than the author found, including Boatto's own preprint.

## 1. Re-derivation (by hand; all agree with the note)

- **Hamiltonian.** sinh(d/2) = |z−w|/√((1−|z|²)(1−|w|²)). Each k is in N−1 pairs, so
  −Σ ln sinh = −Σ ln|z_j−z_k| + ((N−1)/2) Σ ln(1−|z_k|²). ✓
- **Log-polar form.** z = √ξ e^{s+iθ}: dx∧dy = |z|² ds∧dθ, so κ(s) = 4ξe^{2s}/(1−ξe^{2s})². ✓
- **Pair expansion.** ln(2coshΔs − 2cos(φ+Δθ)). The Δθ² coefficient is cosφ/D − 2sin²φ/D². Using sin²φ = D(1+cosφ)/2
  this is −1/D, and the Δs² coefficient is 1/D. There is no ΔsΔθ term, and the linear θ-term cancels under l ↦ N−l.
  So the quadratic part is (Δs²−Δθ²)/D. ✓
- **Fourier lemma.** The proof in the note (|1−w^m|²/|1−w|² = |Σw^a|², orthogonality, |a−b| < m ≤ N) is correct, including m = N.
  It is also equivalent to M–T (2.4) minus its ℓ = 0 case. On mode m, Σ_{j<k}|x_j−x_k|²/D_{jk} = (N/2)·m(N−m), and
  Σ|x_k|² = N, so the Hessian entries are ∓m(N−m)/2. (The factor 2 between the quadratic form and the Hessian is consistent.) ✓
- **One-body terms.** W'(0) = −(N−1)ξ/(1−ξ). W''(0) = −2(N−1)ξ/(1−ξ)². J₁' = 4u/(1−u)². J₁'' = 8u(1+u)/(1−u)³.
  The equilibrium condition −(N−1)/2 + W'(0) = ΩJ₁'(0) gives Ω = −(N−1)(1+ξ)(1−ξ)/(8ξ).
  −ΩJ₁''(0) = (N−1)(1+ξ)²/(1−ξ)², so C₁ = (N−1)(1+ξ²)/(1−ξ)². ✓
- **Linearisation.** At a critical point of F the variation of κ drops out. κ₀ṡ = bθ and κ₀θ̇ = −as give σ² = −ab/κ₀². ✓
  b_m > 0 for 1 ≤ m ≤ N−1. a₁ = C₁ − (N−1)/2 > 0 always, and a₀ = C₁ > 0. So the Hessian is positive definite off the
  rotation direction iff C₁ > T. ✓
- **cosh² form.** ξ = tanh²(r/2) (the disk metric 4|dz|²/(1−|z|²)² gives r = 2 artanh r_E). cosh r = (1+ξ)/(1−ξ), and
  (1+cosh²r)/2 = (1+ξ²)/(1−ξ)². For odd N, T = (N²−1)/8 and q = (N−3)/4. For even N, T = N²/8 and q = (N−2)²/(4(N−1)). ✓
  q ≤ 1 for N ≤ 7 (q = 1/3, 1/2, 4/5, 1 for N = 4..7), and q > 1 for N ≥ 8. C₁ is increasing in ξ, so the ring is stable iff r > r*. ✓
- **ξ*+1/ξ*.** This equals 2(q+1)/(q−1) = 2 + 4/(q−1). For odd N it is 16/(N−7). For even N it is 16(N−1)/(N²−8N+8). ✓
- **Units.** For odd N, (N−7) | 16 with N−7 even gives N ∈ {9, 11, 15, 23} and c = 10, 6, 4, 3.
  For even N, N²−24N+24 has roots 12 ± √120 ≈ 1.05 and 22.95. So 0 < 4/(q−1) < 1 for N ≥ 24, while N = 22 gives −20 < 0,
  so the case check up to 22 is needed. I checked N = 10, …, 22 by hand: 144/28, 176/56, 208/92, 240/136, 272/188, 304/248 and 336/316
  are all non-integers. N = 8 gives 112/8 = 14, so c = 16. ✓
  Converse: if ξ* ∈ ℚ, an algebraic integer in (0,1) is impossible. Otherwise the minimal polynomial is x² − cx + 1, and
  integrality gives a minimal polynomial in ℤ[x], so c ∈ ℤ. Correct. ξ*(10) = (4/3−1)/(4/3+1) = 1/7 ✓. ξ*(8) = 8 − 3√7 ✓.
  Independent exact (Fraction) scan, 8 ≤ N < 5000: c ∈ ℤ only for N = 8, 9, 11, 15, 23 (c = 16, 10, 6, 4, 3).
- **Table.** I recomputed q, ξ*, √ξ* and r* = arccosh √q = 2 artanh √ξ* (both formulas agree) for N = 8..15. All entries match. The digits of ξ*
  in the note are truncated, not rounded (for example 0.20241086403185…), which is acceptable.

## 2. Independent numerics (NUMERICAL) — `sinh/review_1_stab.py` → `sinh/review_1_stab.out`

This is a different method from the author's `lin_stab.py`. The exact second derivatives of the pair term and of the one-body term
(H − ΩJ) are produced symbolically by sympy and assembled into the full 2N×2N Cartesian Hessian, evaluated at 50 digits. There is no finite differencing
and no Fourier or log-polar reduction. Ω is solved from ∂F/∂x₀ = 0, which is linear in Ω; the formula is not used.
The linearised field is ẋ = λ⁻¹∂_yF, ẏ = −λ⁻¹∂_xF, with λ = 4/(1−ξ)².
- (A) N = 8..15, ξ ∈ {0.03, 0.8ξ*, 1.2ξ*, 0.6}. The squared eigenvalues of the full Jacobian match {0,0} ∪ {−a_m b_m/κ₀²}
  to ≤ 1.5·10⁻⁴⁵. Ω matches −(N−1)(1+ξ)(1−ξ)/(8ξ) and the ring gradient vanishes to ≤ 5·10⁻⁴⁸.
- (B) The threshold is the sign change of the largest nonzero σ², found by bisection on [0.7ξ*, 1.3ξ*] (60 steps).
  It agrees with ξ*(N) to (1.6–7.0)·10⁻²⁰ for every N = 8..15. This is the bisection resolution.
- I also re-ran the author's `verify_core.py`. Its output is byte-identical to `verify_core.out` (ALL PASS).

## 3. Lean

`export PATH=$HOME/.elan/bin:$PATH; cd lean && lake env lean Havelock/SinhArithmetic.lean` exits 0 with no
output. There is no `sorry`, `admit`, `axiom` or `native_decide`.
`#print axioms` (in a temporary uncommitted copy of the file) gives [propext, Classical.choice, Quot.sound] for all 15
declarations.
Statement audit:
- `q` and `c` are exactly the note's q(N) and c(N). Lean's x/0 = 0 does not matter because q − 1 ≠ 0 for N > 7.
- `c_int_iff` holds for all N > 7, the hypothesis is not vacuous, and both directions are proved.
- `even_dvd_iff`/`odd_dvd_iff` are the note's case analysis, including the N ≥ 23 bound and `interval_cases` up to 22.
- `xi_add_inv` holds for real q > 1.
- `isIntegral_of_c_int` is general: ξ ≠ 0 and ξ + ξ⁻¹ = k ∈ ℤ give integrality and ξ(k−ξ) = 1.
- `xi8_eq` and `xi10_eq` are stated with literal √(9/7) and √(16/9); `q8` and `q10` link these to q separately.
Cosmetic gap: no single theorem composes `xi_add_inv` (ℝ) with `c` (ℚ). That composition is a trivial cast. The note
correctly lists what is not formalised: the dynamics, the lemma, and the converse.

## 4. Literature (checked myself from primary PDFs)

- **Montaldi–Tokieda 2013, arXiv:1210.5662: claims CONFIRMED.**
  - Green's function (1.12) is log|z−w|²/((1+λ|z|²)(1+λ|w|²)). At λ < 0 this is the sinh model (curvature normalisation 4λ, same z). Symplectic form (1.8) is 2dx∧dy/σ².
  - (2.2): the Hessian has no r–θ coupling.
  - (2.6): ε_r^(ℓ) = κ²/(4πr₀²)[2(n−1)(1+λ²r₀⁴)/σ² − ℓ(n−ℓ)] and ε_θ^(ℓ) = κ²ℓ(n−ℓ)/4π. With λ = −1 and r₀² = ξ this is 2a_ℓ/r₀², a scaled copy of the note's (a_ℓ, b_ℓ).
  - (2.11): (1+λ²r₀⁴)/(1+λr₀²)² > ⌊n²/4⌋/(2(n−1)), which is exactly C₁ > T.
  - Theorem 2: Lyapunov stability (G_µ sense) if (2.11) holds. Remark 4: linear instability if it is reversed.
  - Table: λr₀² = −0.0627, −0.101, −0.143, −0.172, −0.202, −0.225 for n = 8..13, which is −ξ*(N).
  - "Hyperbolic planes" paragraph: n ≤ 7 always stable, and "n = 15 stable for r₀² > (2−√3)/|λ|".
  - §4.2: the image Green's function log|z−w|²/|1+λzw̄|² "On the hyperbolic plane λ < 0 … reduces to the Hamiltonian
    adopted by Kimura [10]", and its stability problem "seems no longer tractable analytically".
- **Montaldi–Nava-Gaxiola 2014, arXiv:1403.2138: CONFIRMED.**
  - Eq. (12): H = −(1/4π)ΣΓᵢΓⱼ ln[(⟨Xᵢ,Xⱼ⟩+1)/(⟨Xᵢ,Xⱼ⟩−1)] "as constructed by Kimura [10]".
  - Eq. (1): ⟨X,Y⟩_H = x₁x₂ + y₁y₂ − z₁z₂. Eq. (2): ⟨X,Y⟩_H = −cosh d. So the ratio is tanh²(d/2).
  - Their remark that the pair energy tends to 0 as d → ∞ holds only for tanh.
  - Only N = 2, 3 are treated.
- **Additional sources found by this reviewer, all quoting Kimura's H² Green's function as tanh:**
  1. Boatto, "Curvature perturbations and stability of a ring of vortices", mp_arc 06-311 (2006 preprint of the DCDS-B
     2008 paper). The text has Type-3 fonts and is unextractable, so I OCR'd it (tesseract) and read p. 7 as an image.
     - Eq. (14), hyperbolic: 2πΨ = −ln tanh(√|K| r/2).
     - Eq. (17): H_{K<0} = c Σ ΓΓ ln((d_αβ−1)/(d_αβ+1)), with c = −1/4π and d = cosh(dist). That is the tanh model.
     - Eq. (18): H_{S²} uses ln(1−d), the sin/chord model.
     - Boatto says her formulation "was greatly inspired" by Kimura.
     - The H² ring stability result there is numerical only (Fig. 5, minimum Hessian eigenvalue versus K). "N > 7 can only be assured to be
       nonlinearly stable for sufficiently negative values of K."
     - The only closed-form criterion is for S² (Thm 5.3, attributed to Boatto–Simó).
     - **So Boatto (preprint version) does not contain C1/C2 for the sinh model.** The published 2008 version was not
       obtained and could differ.
  2. Garra, "Localization of point vortices under curvature", Math. Mech. Complex Syst. 1(1) (2013), msp.org. Eq. (3-3)
     attributes 2πG_{K<0} = −ln tanh(√|K| r/2) to Kimura 1999. Eq. (3-1) is the sphere −ln sin.
  3. Drivas–Khanikati–Khanikati, arXiv:2607.16490 (2026), eq. (11): "The Green's functions on these surfaces are given
     by [13 = Kimura 1999]": G_{H²} = −(1/2π) log tanh(d/2), with G_{S²} = −(1/2π) log sin(d/2).
  4. Nava-Gaxiola thesis (Manchester 2013), eq. (4.18): tanh form attributed to Kimura. It says Hwang–Kim 2009 uses
     different equations of motion (an extra factor in the denominator), and that Hwang–Kim 2013's Hamiltonian coincides with
     Kimura's (tanh).
- **Kimura 1999 itself:** not obtained. royalsocietypublishing.org returns 403 to both curl and WebFetch. Crossref confirms
  the DOI 10.1098/rspa.1999.0311 and pp. 245–259, with no abstract. Unpaywall shows no open-access copy.
- **Hwang–Kim 2009** (J. Geom. Phys. 59, doi 10.1016/j.geomphys.2009.01.003): Unpaywall lists a publisher open-access PDF, but
  ScienceDirect returns 403. Not read.
  **Hwang–Kim 2013** (J. Math. Phys. 54, 063101, doi 10.1063/1.4811454): only the Crossref abstract was available ("relative equilibria … the
  configuration of rings of vortices is investigated in detail"). It does not mention stability thresholds, and
  zbMATH's review is blocked. Per the thesis, its Hamiltonian is the tanh one. So even if it contains ring thresholds,
  they would bear on C5 (tanh), not on C1–C3.
  **C3 (unit classification): not found anywhere.** The novelty label POSSIBLY-NEW (immaterial) stands.

**Attribution verdict.** The note's caution is correct and, if anything, too weak. Six secondary sources (M–T 2013 §4.2,
M–NG 2014 eq. 12, Nava-Gaxiola thesis eq. 4.18, Boatto preprint eqs. 14/17, Garra 2013 eq. 3-3, plus Drivas et al. 2026
eq. 11) identify Kimura's H² Hamiltonian/Green's function as the **tanh** (image) one. None identifies it with ln sinh.
Kimura's original text remains unread, so the wording "per the secondary literature" is right. Do not call the sinh model "Kimura's".

## 5. Issues (none blocking)

1. **The project README (at the time of review) still said "Kimura's H² model"/"Kimura sinh model"** for the sinh
   Hamiltonian and had no entry recording this misattribution; the name came from earlier project material. The author
   had proposed the fix in the log but not yet applied it. It should be applied before the note is final. The README
   also called the Euclidean radii √ξ* "r*". The note uses r* for the hyperbolic radius, so the README should rename them.
2. Note, abstract: "linearly stable iff C₁ > T". This is correct if understood modulo the rotation direction. The m = 0 block
   (a₀ > 0, b₀ = 0) is a nilpotent drift, and at C₁ = T the critical block is nilpotent. Add one sentence saying that
   equality is degenerate (M–T §3 treat it by bifurcation analysis).
3. Note: Boatto is described as "may be in Boatto earlier (not verified)" and "text unreadable". Update it with §4 above: in the
   preprint, the H² model is tanh and the H² results are numerical; the closed form is only for S².
4. Note §2: the closing sentence on M–T Theorem 2 should mention that M–T's symplectic form and curvature are scaled
   (2dx∧dy/σ², curvature 4λ). This scaling does not change the thresholds in ξ = −λr₀². The note says "up to the curvature normalisation", which is adequate.
5. Lean (cosmetic): there is no single theorem composing `xi_add_inv` with `c`, and `xi8_eq` uses the literal √(9/7) rather than `q 8`.

## 6. Commands

```
cd sinh && python3 review_1_stab.py            # 22.6 min; output review_1_stab.out
cd sinh && python3 verify_core.py | diff - verify_core.out   # identical
export PATH=$HOME/.elan/bin:$PATH; cd lean && lake env lean Havelock/SinhArithmetic.lean   # exit 0
# plus a temporary copy of the file with #print axioms for all 15 decls
```
Primary PDFs: arxiv.org/pdf/1210.5662, arxiv.org/pdf/1403.2138, web.ma.utexas.edu/mp_arc/c/06/06-311.pdf (OCR),
msp.org/memocs/2013/1-1/memocs-v1-n1-p02-s.pdf, arxiv.org/pdf/2607.16490,
eprints.maths.manchester.ac.uk/2022/1/thesis_final.pdf. Local copies of the papers are not committed (copyright).
