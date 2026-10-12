> Independent AI reader's report (2026-10-12) on Proposition 19 and Corollary 20 of `tanh/note_tanh.tex`. All GAP and WORDING
> items below were applied before commit; line references are to the draft it reviewed. The reader's scripts are kept locally.

# Referee report: horocyclic vortex row (tanh/note_tanh.tex, Sec. "Large N", subsection "Physical interpretation")

Scope: the paragraph "The limit object is a vortex row on a horocycle", Prop. `prop:row` with its proof, and Cor. `cor:row` with its proof.
I re-derived everything by hand and checked it numerically with my own code (scripts and outputs are in `review_row/`).

**Verdict: ACCEPT WITH MINOR REVISIONS.** I found no mathematical errors. Every formula checks. The issues are two small gaps in the corollary's proof and some wording about edge modes and the meaning of "stability" for an infinite row.

## 1. Hand re-derivation

**Model.** On the upper half-plane, tanh(d/2) = |z-w|/|z-w̄|, so the pair term is -ln|z-w| + ln|z-w̄|. This is the same intrinsic H as the disk model (Cayley is an isometry). There is no self-term, consistent with "no Robin term".

**Hamilton's equations and momentum.** With area form dx dy/y² and the note's convention (xdot = λ⁻¹ ∂_y H, ydot = -λ⁻¹ ∂_x H, λ = y⁻²), the equations are xdot = y² ∂_y H and ydot = -y² ∂_x H. The x-translation field (1,0) needs y² ∂_y J = 1, so J = -1/y and J = -Σ 1/y_j. Correct.

**Relative equilibrium.** If dF = 0 with F = H - vJ, then xdot = v y² ∂_y J = v. So v is the drift velocity in +x.

**Force balance.**
- x-forces vanish by k ↔ -k symmetry.
- The pair part gives no y-force, since z0 - zk is real.
- The image part gives ∂_{y0} ln|z0 - z̄k| = 2/(k²ℓ² + 4).
- So v = Σ_{k≠0} 2/(k²ℓ² + 4) = (2/ℓ²)[(π/a) coth(πa) - 1/a²] with a = 2/ℓ, πa = τ. This equals (τ coth τ - 1)/2. Correct.
- Sign check: positive vortices above the absorbing boundary drift in +x, as with a wall image. Consistent.

**Hessian.** For h = Re f, Hess h = [[Re f'', -Im f''], [-Im f'', -Re f'']].
- Pair term, f = -log(z0 - zk): f'' = (kℓ)⁻². On-site this gives diag(1,-1)(kℓ)⁻². The coupling is minus that, i.e. diag(-1,1)(kℓ)⁻².
- Image term, g = log(z0 - z̄k): g'' = -p_k with p_k = (kℓ - 2i)⁻². On-site this gives [[-Re p, Im p], [Im p, Re p]]. The off-diagonal Im p_k is odd in k, so it cancels in pairs.
- Coupling for the image term: ∂_{xk} = -∂_X and ∂_{yk} = +∂_Y. This gives [[Re p, Im p], [-Im p, Re p]], exactly as stated.
- The -vJ term adds +2v/y³ = 2v to the on-site yy entry.
- Hermiticity: C_{-k} = C_kᵀ because p_{-k} = conj(p_k).
- Decay of the coupling blocks: the xx entry is O(k⁻⁴) (the pair and image parts cancel at order k⁻²), xy is O(k⁻³), and yy is O(k⁻²). All are absolutely summable, as claimed.

**Lattice sums.**
- S1 = τ²/12.
- Σ_{k≠0} cos(qk)/(kℓ)² = τ²/12 - τ²θ(1-θ)/2, valid on 0 ≤ q ≤ 2π.
- Poisson summation, with convention f̂(ω) = ∫ f(t) e^{-iωt} dt and s = q - 2πn:
  - For s > 0, close the contour upward. The double pole at ib gives 2πi · (is e^{-bs}) = -2π s e^{-bs}.
  - For s < 0 the transform is 0. For s = 0, ∫(t - ib)⁻² dt = 0.
- Endpoints: at q = 0 the n = 0 term vanishes. At q = 2π the extra n = 1 term has s = 0 and vanishes. So the formula holds on the closed interval [0, 2π] with n ≥ 0 after re-indexing.
- Validity: f = O(t⁻²) and f̂ is continuous with exponential decay, so Poisson summation applies.
- Result: Σ_{k∈Z} e^{iqk} p_k = -τ² e^{-2θτ}[θ/(1-X) + X/(1-X)²]. Removing p_0 = -1/4 gives P(θ). Correct.
- P is real. The cos and sin projections are (P(θ) + P(1-θ))/2 and (P(1-θ) - P(θ))/2. Correct.
- The bracket in P equals W_m/N with ξ^m = e^{-2θτ}. I checked this using 1 - X = 2e^{-τ} sinh τ.

**Reduction to β, α, γ.**
- M_xx = τ²[(1-x²)/8 + (1-Q)/(4 sinh²τ)] = τ²β.
- M_yy = τ coth τ - 1/2 - τ²(1-x²)/8 - τ²(1+Q)/(4 sinh²τ) = τ²α, using τ²E_∞ = τ coth τ - 1/2.
- M_xy = i Σ sin(qk) Im p_k = iτ²γ.
- All three match the note.

**Linearised flow, corollary item (i).** Write V = P(y) dH with P = y² J0. Then DV = P Hess H + (∂P)·dH. At the row dH = v dJ, and P dJ = (1,0) is constant. Hence (∂P)·dJ = -P Hess J, and DV = P Hess F.
- At y = 1 this is J0 M(θ) on Bloch modes. The eigenvalues are -iγ' ± sqrt(-M_xx M_yy) with γ' = τ²γ.
- The claim is right. The reason, though, is dF = 0, not only "the area form is standard at y = 1".

**Spectral stability for the infinite row, corollary item (ii).** The linear operator is a bounded block convolution on ℓ²(Z; C²). Its symbol is continuous, so its spectrum is the union over θ of spec J0 M(θ). The definition "every Bloch exponent purely imaginary" is the standard spectral-stability definition and is acceptable.
- For 0 < θ < 1, β > 0 (Lemma bpos). So stability holds iff α(x) ≥ 0 for all x.
- min α = α(0) = Φ (Lemma crit). Φ is strictly decreasing, so stability holds iff τ ≤ 2η_∞, i.e. ℓ ≥ ℓ* = π/η_∞.
- At θ = 0: β = γ = 0, M(0) = diag(0, τ²α(1,τ)), and the exponents are 0. The product condition M_xx M_yy ≥ 0 holds trivially.
- The boundary ℓ = ℓ* is included correctly: at θ = 1/2 the exponent is 0.

**Definiteness and Krein, corollary item (iii).**
- αβ - γ² = β(α - γ²/β) = β(Φ + G), with G as in Lemma em. That lemma applies verbatim because it is a pointwise identity in (x, τ).
- So M(θ) > 0 for all 0 < θ < 1 iff Φ > 0, i.e. iff ℓ > ℓ*.
- At ℓ = ℓ*, γ(0, τ) = 0 (R(0, τ) = 0) and M(1/2) = diag(β, 0). The only degeneracy is at θ = 1/2, since G > 0 elsewhere.
- The spectral and definiteness thresholds coincide, so there is no Krein window. Correct.

## 2. Numerical verification

All scripts are my own and run with the repo's .venv python.

- **`check_row.py` / `check_row.out`**
  - v: mpmath brute sum vs closed form agree to all 15 printed digits at ℓ = 0.3, 0.8439, 1.5, 4.
  - P(θ): brute sum with K = 2·10⁶ vs closed form at θ ∈ {0, 0.1, 0.25, 0.5, 0.77, 1}, including both endpoints. Agreement is better than 1·10⁻¹².
  - Full symbol M(θ): Hessian blocks from sympy second derivatives of the pair function (independent of the note's Re/Im bookkeeping), summed to 2·10⁶. Compared with τ²[[β, iγ], [-iγ, α]] at ℓ ∈ {0.4, ℓ*, 1, 2.5} and θ ∈ {0, 0.05, 0.2, 0.5, 0.63, 0.9, 1}.
    - For 0 < θ < 1 the difference is at most 3·10⁻¹².
    - At θ = 0 and 1 the difference is 1.2·10⁻⁵. This is exactly the truncation tail 4/(ℓ²K) of the non-oscillating yy coupling, not a formula error.
  - At ℓ = ℓ*, M(1/2) = diag(6.913, 0). Confirmed.
- **`dyn_row.py`, `dyn_row2.py` / `.out`** (independent of the Hessian and of claim (i))
  - Setup: a periodic row of 12 vortices with all images summed (|j| ≤ 4000, plus an analytic tail). The Jacobian is taken by finite differences of the true velocity field (y² ∂_y H, -y² ∂_x H).
  - Equilibrium: the mean x-velocity equals v to 1.5·10⁻⁸, and max |ydot| is about 10⁻¹².
  - Raw-Jacobian eigenvalues match the symbol prediction at θ = j/12 to 2·10⁻⁷, away from the θ = 0 Jordan pair.

| ℓ/ℓ* | 0.9 | 0.98 | 0.995 | 1.005 | 1.02 | 1.2 |
|---|---|---|---|---|---|---|
| max Re σ (raw) | 2.599 | 0.975 | 0.472 | 1e-10 | 2e-10 | 3e-10 |
| max Re σ (symbol) | 2.599 | 0.975 | 0.472 | 0 | 0 | 0 |

  - At ℓ = ℓ* exactly, the θ = 1/2 pair from the raw Jacobian is ±4·10⁻⁵ (real). This is the square root of the finite-difference error, the signature of a nilpotent 2×2 Jordan block (J0 diag(β, 0)).
  - At θ = 0 the raw Jacobian shows ±1.7·10⁻⁵ i. This is also a Jordan pair, for the same reason.
- **`edge_checks.py` / `.out`**
  - τ²α(1, τ) > 0 for every τ sampled, and β(1) = γ(1) = 0.
  - Near θ = 0: β/θ² → about 6.8. The smallest eigenvalue of M(θ)/τ² is about 4.9·10⁻⁶ at θ = 10⁻³, so it tends to 0 (acoustic branch). The exponents are about ±0.92θ i.

## 3. Issues

**Issue 1 (GAP, minor). Corollary proof: the reason for the formula σ = -iγ' ± sqrt(-M_xx M_yy).** The proof says only that "the area form is standard at y = 1". But the Poisson tensor y² J0 is not constant, and the linearisation drops the term (∂P)·dH only because dF = 0 at the row. That is the argument of Prop. lin.
- Fix: "Since dF = 0 at the row, the linearisation of xdot = y² ∂_y F, ydot = -y² ∂_x F at y = 1 is J0 Hess F."
- Also define γ' := τ²γ = -i M_xy. Right now γ' is used without a definition.

**Issue 2 (GAP, minor). The θ = 0 mode in the corollary proof.** "M_xx ∝ β > 0 (Lemma bpos)" holds only for 0 < θ < 1, because Lemma bpos is for |x| < 1. At θ = 0 (and θ = 1) we have β = γ = 0 and M(0) = diag(0, τ²α(1, τ)).
- Fix: add a sentence. At θ = 0 the exponents are 0, so that mode does not affect the criterion. The null direction is the x-translation symmetry.
- Also note τ²α(1, τ) = (τ coth τ - 1) + ½(1 - τ²/sinh²τ) > 0, so M(0) is positive semidefinite with kernel exactly the translation direction.

**Issue 3 (WORDING). "Spectrally stable" vs. linear stability, including the boundary.** The equivalence "stable iff ℓ ≥ ℓ*", boundary included, is correct for the stated definition. However:
- At ℓ = ℓ* the θ = 1/2 block J0 M = [[0, 0], [-M_xx, 0]] is a nontrivial Jordan block, so that mode grows linearly.
- At θ = 0 there is always a Jordan block: a uniform height change makes the row drift, because the velocity changes.
- So the row is spectrally stable but never linearly stable in the bounded sense.
- Fix: one sentence saying that "spectrally stable" means the spectrum is on iR, so there is no exponential growth; algebraic growth occurs at θ = 0 for every ℓ and at θ = 1/2 when ℓ = ℓ*.

**Issue 4 (WORDING). Meaning of "positive definite on every Bloch mode with 0 < θ < 1" and of "no Krein window".** The claim is correct pointwise in θ. But two things should be said:
1. θ = 0 is excluded for the same reason mode 0 is excluded for the ring: it holds the translation symmetry and the height (momentum) direction.
2. The definiteness is not uniform. Because β = O(θ²), the smallest eigenvalue of M(θ) tends to 0 as θ → 0 (acoustic branch). Also, H and J are infinite for the row.

So pointwise definiteness gives no energy–momentum or Lyapunov conclusion for the infinite row. "No Krein window" should be read as "the spectral and definiteness thresholds coincide, so all Krein signatures are definite for ℓ > ℓ*".
- Fix: add a sentence saying so, and that no nonlinear statement is made for the row.

**Issue 5 (WORDING). "Its Bloch symbol is the large-N limit of the ring blocks".** This is true up to reordering the basis and an orientation convention.
- The ring block in (s, θ) is [[a, ig], [-ig, b]]. In the order (tangential, normal) it becomes [[b, -ig], [ig, a]].
- In the limit that is the complex conjugate of the row's matrix, i.e. the row's M(1-θ). This comes from s ↔ -δy: the outward normal of the ring corresponds to -y for the horocycle bounding the horoball y > 1.
- The scalings also differ: N⁻² for the ring vs (ℓ/2π)² for the row.
- Fix: add "after reordering the basis to (tangential, normal) and up to the orientation θ ↦ 1 - θ (equivalently γ ↦ -γ)".

**Issue 6 (WORDING). The Poisson-summation sentence.** State the Fourier convention f̂(ω) = ∫ f e^{-iωt} dt and why Poisson summation applies (f = O(t⁻²), f̂ continuous and exponentially decaying). Also note the endpoint handling: the terms with s = 0 vanish, so the formula holds on the closed interval [0, 2π]. Everything stated is already correct; this only makes it self-contained.

**Issue 7 (WORDING, optional). The infinite row's equations of motion.** Separately, the pair force and the image force from site k decay only like 1/k and sum only conditionally. The tanh interaction taken as a unit gives forces of order k⁻³ (x) and k⁻² (y), which are absolutely convergent. The note's "absolutely summable" remark covers the Hessian. A half-sentence saying the vector field is also absolutely convergent would make "relative equilibrium" fully rigorous for the infinite system.

## 4. Correct use of earlier results

- **Lemma bpos:** used for 0 < θ < 1. Correct, with the θ = 0 caveat of Issue 2.
- **Lemma crit** (ψ even and strictly decreasing on [0, 1], for every τ): gives min α = Φ, with the minimum only at x = 0. Correct.
- **Lemma em:** a pointwise identity in (x, τ). It applies with Φ in place of Φ_N because αβ - γ² = β(Φ + G) follows from the definitions. Correct.
- **Monotonicity of Φ** (from eq. phi): gives the unique crossing at τ_∞. Correct.
- **Prop. lin:** used as a template. It applies, given the dF = 0 argument of Issue 1.
