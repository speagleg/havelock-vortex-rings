import Havelock.Tanh.BlocksCalc

/-! # Tanh model: derivation of the Fourier blocks (note Lemma `lem:re`, Thm `thm:blocks`)

`General.lean` takes the blocks `aBlock`, `bBlock`, `gBlock` as definitions. Here they are derived as the
Fourier symbol of the Hessian of `F = H − ΩJ` (`ringF`, `BlocksCalc.lean`) at the ring `s = θ = 0`, in the
log-polar coordinates `z_j = √ξ w^j e^{s_j + iθ_j}`.

**Notions used.** `hessQ s θ = d²/dt² F(ts, tθ)|_{t=0}` (iterated one-variable `deriv`) is the Hessian
quadratic form; the Hessian entries are obtained by polarisation,
`hessB x y = (hessQ(x+y) − hessQ(x−y))/4`, and
`hessSS a b = hessB (e_a,0) (e_b,0)` (`= ∂²F/∂s_a∂s_b`), `hessTT a b = hessB (0,e_a) (0,e_b)`
(`= ∂²F/∂θ_a∂θ_b`), `hessST a b = hessB (e_a,0) (0,e_b)` (`= ∂²F/∂s_a∂θ_b`),
`hessTS a b = hessB (0,e_a) (e_b,0)` (`= ∂²F/∂θ_a∂s_b`). (For a `C²` function these are the usual
second partial derivatives.) Criticality is stated for the directional derivative along every line through
the ring.

Main results (all for `0 < ξ < 1`):
* `ring_critical` (note Lemma `lem:re`): with `Ω = ringOmega N ξ`, i.e. `Ωκ₀ = omegaKappa N ξ`,
  `t ↦ F(ts, tθ)` has derivative `0` at `t = 0` for every direction `(s, θ)`;
  `grad_formula`: in general that derivative is `(omegaKappa N ξ − Ωκ₀) Σ_j s_j`.
* (`BlocksHess.lean`) `hessQ_eq`, `hessSS_eq`, `hessTT_eq`, `hessST_eq`, `hessTS_eq`: closed forms of the entries (circulant in `a − b`).
* (`BlocksFourier.lean`) `blocks_eigen` (note Thm `thm:blocks`): for `1 ≤ m ≤ N−1` and `v_b = w^{mb}`,
  `Σ_b ∂²F/∂s_a∂s_b v_b = a_m v_a`, `Σ_b ∂²F/∂θ_a∂θ_b v_b = b_m v_a`, `Σ_b ∂²F/∂s_a∂θ_b v_b = i g_m v_a`,
  `Σ_b ∂²F/∂θ_a∂s_b v_b = −i g_m v_a`: the Hessian is block-circulant and acts on the mode
  `(S,T) w^{mj}` as `M_m = [[a_m, i g_m], [−i g_m, b_m]]`. `blocks_eigen_zero`: on mode `0`, `b_0 = g_0 = 0`.
* `blocks_hermitian`: `Σ_{a,b} conj(w^{m'a}) ∂²F/∂s_a∂s_b w^{mb} = N a_m [m = m']` for `0 ≤ m, m' ≤ N−1`
  (`1 ≤ m`), and likewise `N b_m`, `N (i g_m)`, `N (−i g_m)` for the other three blocks
  (cross terms between different Fourier modes vanish).
-/

noncomputable section

namespace Havelock.Tanh

open scoped ComplexConjugate

/-! ## Symmetric double sums -/

theorem sum_lt_add_swap {N : ℕ} (f : Fin N → Fin N → ℝ) :
    ∑ j : Fin N, ∑ k : Fin N, (if j < k then f j k + f k j else 0)
      = ∑ j : Fin N, ∑ k : Fin N, (if j = k then 0 else f j k) := by
  have e1 : ∀ j k : Fin N, (if j < k then f j k + f k j else 0)
      = (if j < k then f j k else 0) + (if j < k then f k j else 0) := by
    intro j k; split_ifs <;> ring
  have e2 : ∀ j k : Fin N, (if j = k then 0 else f j k)
      = (if j < k then f j k else 0) + (if k < j then f j k else 0) := by
    intro j k
    rcases lt_trichotomy j k with h | h | h
    · simp [h, h.ne, not_lt.mpr h.le]
    · simp [h]
    · simp [h, h.ne', not_lt.mpr h.le]
  simp_rw [e1, e2, Finset.sum_add_distrib]
  congr 1
  rw [Finset.sum_comm]

/-! ## Criticality (note Lemma `lem:re`) -/

/-- `Ω` of note Lemma `lem:re`: `Ω κ₀ = −(N−1)/2 + ξ/(1−ξ) − N X/(1−X)` (`omegaKappa`). -/
def ringOmega (N : ℕ) (ξ : ℝ) : ℝ := omegaKappa N ξ / kappa0 ξ

theorem zeta_swap (N : ℕ) (j k : Fin N) : zeta N k j = conj (zeta N j k) := by
  unfold zeta; simp only [map_mul, Complex.conj_conj]; ring

/-- One half of the first-derivative pair term. -/
def uPair (N : ℕ) (ξ : ℝ) (s θ : Fin N → ℝ) (j k : Fin N) : ℝ :=
  -(epsC s θ j * (1 / (1 - conj (zeta N j k)) + rho1 ξ (zeta N j k))).re

theorem lPair_eq {N : ℕ} (hN : 0 < N) (ξ : ℝ) (s θ : Fin N → ℝ) {j k : Fin N} (hjk : j ≠ k) :
    lPair N ξ s θ j k = uPair N ξ s θ j k + uPair N ξ s θ k j := by
  have hw := ringW_pow_ne hN hjk
  have ha : ringW N ^ (j : ℕ) ≠ 0 := pow_ne_zero _ (ringW_ne_zero N)
  have hb : ringW N ^ (k : ℕ) ≠ 0 := pow_ne_zero _ (ringW_ne_zero N)
  have hca := conj_pow_ringW N j
  have hcb := conj_pow_ringW N k
  unfold lPair uPair
  rw [zeta_swap N j k, rho1_conj, Complex.conj_conj]
  generalize rho1 ξ (zeta N j k) = r
  have hz : zeta N j k = ringW N ^ (j : ℕ) * (ringW N ^ (k : ℕ))⁻¹ := by unfold zeta; rw [hcb]
  have hzc : conj (zeta N j k) = (ringW N ^ (j : ℕ))⁻¹ * ringW N ^ (k : ℕ) := by
    rw [hz, map_mul, hca, map_inv₀, hcb, inv_inv]
  rw [hzc, hz]
  set a := ringW N ^ (j : ℕ)
  set b := ringW N ^ (k : ℕ)
  have hab : a - b ≠ 0 := sub_ne_zero.mpr hw
  have hba : b - a ≠ 0 := sub_ne_zero.mpr (Ne.symm hw)
  have h1 : (1 : ℂ) - a⁻¹ * b ≠ 0 := by
    rw [show (1 : ℂ) - a⁻¹ * b = (a - b) / a by field_simp]; exact div_ne_zero hab ha
  have h2 : (1 : ℂ) - a * b⁻¹ ≠ 0 := by
    rw [show (1 : ℂ) - a * b⁻¹ = (b - a) / b by field_simp]; exact div_ne_zero hba hb
  have e1 : (a * epsC s θ j - b * epsC s θ k) / (a - b)
      = epsC s θ j * (1 / (1 - a⁻¹ * b)) + epsC s θ k * (1 / (1 - a * b⁻¹)) := by
    field_simp; ring
  rw [e1]
  generalize 1 / (1 - a⁻¹ * b) = p
  generalize 1 / (1 - a * b⁻¹) = q
  simp only [Complex.add_re, Complex.mul_re, Complex.conj_re, Complex.conj_im, Complex.add_im,
    ]
  ring

theorem sum_inv_one_sub_conj {N : ℕ} [NeZero N] :
    ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0 else 1 / (1 - conj (ringW N ^ (l : ℕ))))
      = ((N : ℂ) - 1) / 2 := by
  rw [← sum_inv_one_sub (N := N), ← sum_conj (fun ζ => if ζ = 1 then 0 else 1 / (1 - ζ))]
  refine Finset.sum_congr rfl fun l _ => ?_
  have e1 : (conj (ringW N ^ (l : ℕ)) = 1) ↔ (ringW N ^ (l : ℕ) = 1) := by
    constructor
    · intro h; have := congrArg conj h; simpa using this
    · intro h; rw [h, map_one]
  simp only [e1]

theorem rho1_one (ξ : ℝ) : rho1 ξ 1 = ((ξ / (1 - ξ) : ℝ) : ℂ) := by
  unfold rho1; push_cast; ring

theorem rho_one (ξ : ℝ) : rho ξ 1 = ((ξ / (1 - ξ) ^ 2 : ℝ) : ℂ) := by
  unfold rho; push_cast; ring

/-- `Σ_{j<k} lPair_{jk} = −((N−1)/2 + N X/(1−X) − ξ/(1−ξ)) Σ_j s_j` (proof of note Lemma `lem:re`). -/
theorem sum_lPair {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) (s θ : Fin N → ℝ) :
    ∑ j : Fin N, ∑ k : Fin N, (if j < k then lPair N ξ s θ j k else 0)
      = -(((N : ℝ) - 1) / 2 + N * ξ ^ N / (1 - ξ ^ N) - ξ / (1 - ξ)) * ∑ j : Fin N, s j := by
  have : NeZero N := ⟨hN.ne'⟩
  have e0 : ∀ j k : Fin N, (if j < k then lPair N ξ s θ j k else 0)
      = (if j < k then uPair N ξ s θ j k + uPair N ξ s θ k j else 0) := by
    intro j k; split_ifs with h
    · exact lPair_eq hN ξ s θ (ne_of_lt h)
    · rfl
  simp_rw [e0]
  rw [sum_lt_add_swap]
  set G : ℂ → ℂ := fun ζ => if ζ = 1 then 0 else 1 / (1 - conj ζ) + rho1 ξ ζ with hG
  have hK : ∑ l : Fin N, G (ringW N ^ (l : ℕ))
      = ((((N : ℝ) - 1) / 2 + N * ξ ^ N / (1 - ξ ^ N) - ξ / (1 - ξ) : ℝ) : ℂ) := by
    have e : ∀ l : Fin N, G (ringW N ^ (l : ℕ))
        = (if ringW N ^ (l : ℕ) = 1 then 0 else 1 / (1 - conj (ringW N ^ (l : ℕ))))
          + (if ringW N ^ (l : ℕ) = 1 then 0 else rho1 ξ (ringW N ^ (l : ℕ))) := by
      intro l; simp only [hG]; split_ifs <;> simp
    rw [Finset.sum_congr rfl (fun l _ => e l), Finset.sum_add_distrib, sum_inv_one_sub_conj,
      sum_ite_ne_one (fun ζ => rho1 ξ ζ), sum_rho1 hN h0 h1, rho1_one]
    push_cast; ring
  have hrow : ∀ j : Fin N, ∑ k : Fin N, (if j = k then 0 else uPair N ξ s θ j k)
      = -(epsC s θ j * ∑ k : Fin N, G (zeta N j k)).re := by
    intro j
    rw [Finset.mul_sum, Complex.re_sum, ← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    have hz := zeta_eq_one_iff j k
    unfold uPair
    simp only [hG]
    by_cases hjk : j = k
    · rw [if_pos hjk, if_pos (by unfold zeta; exact hz.mpr hjk)]; simp
    · rw [if_neg hjk, if_neg (by unfold zeta; exact fun h => hjk (hz.mp h))]
  simp_rw [hrow]
  have hz : ∀ j : Fin N, ∑ k : Fin N, G (zeta N j k) = ∑ l : Fin N, G (ringW N ^ (l : ℕ)) :=
    fun j => sum_zeta j G
  simp_rw [hz, hK]
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  simp only [epsC, Complex.mul_re, Complex.add_re, Complex.ofReal_re, Complex.mul_im,
    Complex.ofReal_im, Complex.I_re, Complex.I_im, Complex.add_im]
  ring

/-- The directional derivative of `F` at the ring along `(s, θ)`. -/
theorem grad_formula {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ)
    (s θ : Fin N → ℝ) :
    deriv (fun t : ℝ => ringF N ξ Ω (t • s) (t • θ)) 0
      = (omegaKappa N ξ - Ω * kappa0 ξ) * ∑ j : Fin N, s j := by
  rw [grad_eq hN h0 h1, sum_lPair hN h0.le h1, ← Finset.mul_sum]
  unfold omegaKappa
  ring

/-- **Note Lemma `lem:re`.** With `Ω = ringOmega N ξ` (`Ωκ₀ = −(N−1)/2 + ξ/(1−ξ) − NX/(1−X)`), the ring
`s = θ = 0` is a critical point of `F = H − ΩJ`: along every line `t ↦ (ts, tθ)` the derivative at `t = 0`
is `0`. -/
theorem ring_critical {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (s θ : Fin N → ℝ) :
    HasDerivAt (fun t : ℝ => ringF N ξ (ringOmega N ξ) (t • s) (t • θ)) 0 0 := by
  obtain ⟨f1, H, _⟩ := hasD2_ringF hN h0 h1 (ringOmega N ξ) s θ
  have hd := H.1.self_of_nhds
  have hk : kappa0 ξ ≠ 0 := by
    unfold kappa0
    have : (1 : ℝ) - ξ ≠ 0 := by linarith
    positivity
  have := grad_formula hN h0 h1 (ringOmega N ξ) s θ
  rw [hd.deriv] at this
  rw [this] at hd
  unfold ringOmega at hd
  rw [div_mul_cancel₀ _ hk, sub_self, zero_mul] at hd
  exact hd

/-- Uniqueness of `Ω` (for `N ≥ 1`): if the ring is critical then `Ωκ₀ = omegaKappa N ξ`. -/
theorem ring_critical_iff {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ) :
    (∀ s θ : Fin N → ℝ, deriv (fun t : ℝ => ringF N ξ Ω (t • s) (t • θ)) 0 = 0)
      ↔ Ω * kappa0 ξ = omegaKappa N ξ := by
  constructor
  · intro h
    have := h (fun _ => 1) 0
    rw [grad_formula hN h0 h1] at this
    simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one] at this
    have hN' : (N : ℝ) ≠ 0 := by exact_mod_cast hN.ne'
    rcases mul_eq_zero.mp this with h' | h'
    · linarith
    · exact absurd h' hN'
  · intro h s θ
    rw [grad_formula hN h0 h1, h, sub_self, zero_mul]

end Havelock.Tanh
