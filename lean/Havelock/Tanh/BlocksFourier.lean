import Havelock.Tanh.BlocksHess

/-! # Tanh model: the Fourier blocks as the symbol of the Hessian (note Thm `thm:blocks`)

With `Ω = ringOmega N ξ` (note Lemma `lem:re`) and the Hessian entries `hessSS`, `hessTT`, `hessST`,
`hessTS` of `F = H − ΩJ` at the ring (defined by polarisation of `d²/dt² F(ts,tθ)|₀`, see `Blocks.lean`),
for `0 < ξ < 1` and the Fourier vector `v_b = (w^b)^m = w^{mb}`:

* `blocks_eigen` (`1 ≤ m ≤ N−1`): `Σ_b ∂²F/∂s_a∂s_b v_b = a_m v_a`, `Σ_b ∂²F/∂θ_a∂θ_b v_b = b_m v_a`,
  `Σ_b ∂²F/∂s_a∂θ_b v_b = i g_m v_a`, `Σ_b ∂²F/∂θ_a∂s_b v_b = −i g_m v_a`, with `a_m = aBlock N m ξ`,
  `b_m = bBlock N m ξ`, `g_m = gBlock N m ξ` of `General.lean`. So the Hessian is block-circulant and acts on
  the mode `(S,T) w^{mj}` as `M_m = [[a_m, i g_m], [−i g_m, b_m]]`.
* `blocks_eigen_zero`: on mode `0` the `θθ`, `sθ`, `θs` blocks vanish (`b_0 = g_0 = 0`).
* `blocks_hermitian`: for `0 ≤ m, m' ≤ N−1`, `Σ_{a,b} conj(v'_a) ∂²F/∂s_a∂s_b v_b = N a_m [m = m']` (`m ≥ 1`),
  and likewise `N b_m`, `N (i g_m)`, `N (−i g_m)`; cross terms between different modes vanish.
-/

noncomputable section

namespace Havelock.Tanh

open scoped ComplexConjugate

theorem sum_zeta' {N : ℕ} [NeZero N] (j : Fin N) (G : ℂ → ℂ) :
    ∑ k : Fin N, G (zeta N j k) = ∑ l : Fin N, G (ringW N ^ (l : ℕ)) := sum_zeta j G

theorem zeta_eq_one_iff' {N : ℕ} [NeZero N] (j k : Fin N) : zeta N j k = 1 ↔ j = k :=
  zeta_eq_one_iff j k

/-- Circulant eigenvector identity. -/
theorem eig_circ {N : ℕ} [NeZero N] (a : Fin N) (m : ℕ) (D : ℂ) (h : ℂ → ℂ) :
    ∑ b : Fin N, (if a = b then D else h (zeta N a b)) * (ringW N ^ (b : ℕ)) ^ m
      = (D + ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0
          else h (ringW N ^ (l : ℕ)) * conj (ringW N ^ (l : ℕ)) ^ m)) * (ringW N ^ (a : ℕ)) ^ m := by
  have hv : ∀ b : Fin N,
      (ringW N ^ (b : ℕ)) ^ m = conj (zeta N a b) ^ m * (ringW N ^ (a : ℕ)) ^ m := by
    intro b
    unfold zeta
    rw [map_mul, Complex.conj_conj, ← mul_pow]
    congr 1
    have := ringW_pow_mul_conj N a
    linear_combination (-(ringW N ^ (b : ℕ))) * this
  set G : ℂ → ℂ := fun ζ => if ζ = 1 then D else h ζ * conj ζ ^ m with hG
  have step : ∀ b : Fin N, (if a = b then D else h (zeta N a b)) * (ringW N ^ (b : ℕ)) ^ m
      = G (zeta N a b) * (ringW N ^ (a : ℕ)) ^ m := by
    intro b
    by_cases hab : a = b
    · subst hab
      have h1 : zeta N a a = 1 := (zeta_eq_one_iff' a a).mpr rfl
      simp [hG, h1]
    · have h1 : zeta N a b ≠ 1 := fun h => hab ((zeta_eq_one_iff' a b).mp h)
      simp only [hG, if_neg hab, if_neg h1]
      rw [hv b]; ring
  rw [Finset.sum_congr rfl (fun b _ => step b), ← Finset.sum_mul, sum_zeta' a]
  congr 1
  rw [← sum_ite_one (N := N) (fun _ => D), ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun l _ => ?_
  simp only [hG]
  split_ifs <;> simp

theorem rowsum_zeta {N : ℕ} [NeZero N] (a : Fin N) (f : ℂ → ℝ) :
    (((∑ k : Fin N, (if a = k then 0 else f (zeta N a k)) : ℝ)) : ℂ)
      = ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0 else ((f (ringW N ^ (l : ℕ)) : ℝ) : ℂ)) := by
  rw [Complex.ofReal_sum, ← sum_zeta' a (fun ζ => if ζ = 1 then 0 else ((f ζ : ℝ) : ℂ))]
  refine Finset.sum_congr rfl fun k _ => ?_
  by_cases hak : a = k
  · have h1 : zeta N a k = 1 := (zeta_eq_one_iff' a k).mpr hak
    rw [if_pos hak, if_pos h1, Complex.ofReal_zero]
  · have h1 : zeta N a k ≠ 1 := fun h => hak ((zeta_eq_one_iff' a k).mp h)
    rw [if_neg hak, if_neg h1]

/-! ## The three lattice sums -/

/-- `Σ_{l≠0} α(w^l)(1 − conj(w^l)^m)`. -/
def sAlpha (N m : ℕ) : ℂ :=
  ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0
    else ((alphaF (ringW N ^ (l : ℕ)) : ℝ) : ℂ) * (1 - conj (ringW N ^ (l : ℕ)) ^ m))

/-- `Σ_{l≠0} Re ρ(w^l) conj(w^l)^m`. -/
def sBeta (N : ℕ) (ξ : ℝ) (m : ℕ) : ℂ :=
  ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0
    else (((rho ξ (ringW N ^ (l : ℕ))).re : ℝ) : ℂ) * conj (ringW N ^ (l : ℕ)) ^ m)

/-- `Σ_{l≠0} Im ρ(w^l) conj(w^l)^m`. -/
def sGamma (N : ℕ) (ξ : ℝ) (m : ℕ) : ℂ :=
  ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0
    else (((rho ξ (ringW N ^ (l : ℕ))).im : ℝ) : ℂ) * conj (ringW N ^ (l : ℕ)) ^ m)

theorem sAlpha_eq {N : ℕ} [NeZero N] (m : ℕ) (hm : m ≤ N) :
    sAlpha N m = (m : ℂ) * ((N : ℂ) - m) / 2 := by
  rw [← fejerC m hm]
  unfold sAlpha
  refine Finset.sum_congr rfl fun l _ => ?_
  split_ifs
  · rfl
  · unfold alphaF; push_cast; ring

/-- `Σ_l ρ(w^l) (w^l)^m = N W_{N−m}` and its conjugate. -/
theorem sum_conj_rho {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) (m : ℕ) (hm : m ≤ N) :
    ∑ l : Fin N, conj (rho ξ (ringW N ^ (l : ℕ))) * conj (ringW N ^ (l : ℕ)) ^ m
      = ((N * W N (N - m) ξ : ℝ) : ℂ) := by
  have e : ∀ l : Fin N, conj (rho ξ (ringW N ^ (l : ℕ))) * conj (ringW N ^ (l : ℕ)) ^ m
      = conj (rho ξ (ringW N ^ (l : ℕ)) * conj (ringW N ^ (l : ℕ)) ^ (N - m)) := by
    intro l
    rw [conj_pow_eq hN l (N - m) (Nat.sub_le N m), Nat.sub_sub_self hm]
    simp only [map_mul, map_pow]
  rw [Finset.sum_congr rfl (fun l _ => e l), ← map_sum, sum_rho hN h0 h1 (N - m) (Nat.sub_le N m),
    Complex.conj_ofReal]

theorem sBeta_eq {N : ℕ} [NeZero N] {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) (m : ℕ) (hm : m ≤ N) :
    sBeta N ξ m = (((N * W N m ξ + N * W N (N - m) ξ) / 2 - ξ / (1 - ξ) ^ 2 : ℝ) : ℂ) := by
  have hN : 0 < N := Nat.pos_of_ne_zero (NeZero.ne N)
  unfold sBeta
  rw [sum_ite_ne_one (fun ζ => (((rho ξ ζ).re : ℝ) : ℂ) * conj ζ ^ m)]
  simp only [map_one, one_pow, mul_one, rho_one, Complex.ofReal_re]
  have e : ∀ l : Fin N, (((rho ξ (ringW N ^ (l : ℕ))).re : ℝ) : ℂ) * conj (ringW N ^ (l : ℕ)) ^ m
      = (rho ξ (ringW N ^ (l : ℕ)) * conj (ringW N ^ (l : ℕ)) ^ m
        + conj (rho ξ (ringW N ^ (l : ℕ))) * conj (ringW N ^ (l : ℕ)) ^ m) / 2 := by
    intro l; rw [Complex.re_eq_add_conj]; ring
  rw [Finset.sum_congr rfl (fun l _ => e l), ← Finset.sum_div, Finset.sum_add_distrib,
    sum_rho hN h0 h1 m hm, sum_conj_rho hN h0 h1 m hm]
  push_cast; ring

theorem sGamma_eq {N : ℕ} [NeZero N] {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) (m : ℕ) (hm : m ≤ N) :
    sGamma N ξ m = ((N * W N m ξ - N * W N (N - m) ξ : ℝ) : ℂ) / (2 * Complex.I) := by
  have hN : 0 < N := Nat.pos_of_ne_zero (NeZero.ne N)
  unfold sGamma
  rw [sum_ite_ne_one (fun ζ => (((rho ξ ζ).im : ℝ) : ℂ) * conj ζ ^ m)]
  simp only [map_one, one_pow, mul_one, rho_one, Complex.ofReal_im, Complex.ofReal_zero, sub_zero]
  have e : ∀ l : Fin N, (((rho ξ (ringW N ^ (l : ℕ))).im : ℝ) : ℂ) * conj (ringW N ^ (l : ℕ)) ^ m
      = (rho ξ (ringW N ^ (l : ℕ)) * conj (ringW N ^ (l : ℕ)) ^ m
        - conj (rho ξ (ringW N ^ (l : ℕ))) * conj (ringW N ^ (l : ℕ)) ^ m) / (2 * Complex.I) := by
    intro l; rw [Complex.im_eq_sub_conj]; ring
  rw [Finset.sum_congr rfl (fun l _ => e l), ← Finset.sum_div, Finset.sum_sub_distrib,
    sum_rho hN h0 h1 m hm, sum_conj_rho hN h0 h1 m hm]
  push_cast; ring

theorem W_zero_eq (N : ℕ) (ξ : ℝ) (hX : ξ ^ N ≠ 1) : W N 0 ξ = W N N ξ := by
  have hne : (1 : ℝ) - ξ ^ N ≠ 0 := sub_ne_zero.mpr (Ne.symm hX)
  unfold W
  simp only [CharP.cast_eq_zero, pow_zero, zero_mul, zero_div, zero_add, pow_add]
  field_simp
  ring

/-! ## Eigenvector identities for the explicit matrices -/

section eig
variable {N : ℕ} [NeZero N] (ξ Ω : ℝ) (a : Fin N) (m : ℕ)

theorem eig_hSS :
    ∑ b : Fin N, ((hSS N ξ Ω a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
      = ((Vj ξ Ω : ℂ) - sAlpha N m - sBeta N ξ 0 - sBeta N ξ m) * (ringW N ^ (a : ℕ)) ^ m := by
  have e : ∀ b : Fin N, ((hSS N ξ Ω a b : ℝ) : ℂ)
      = if a = b then (((Vj ξ Ω - ∑ k : Fin N,
          (if a = k then 0 else alphaF (zeta N a k) + (rho ξ (zeta N a k)).re)) : ℝ) : ℂ)
        else (fun ζ => ((alphaF ζ - (rho ξ ζ).re : ℝ) : ℂ)) (zeta N a b) := by
    intro b; unfold hSS alphaM betaM; split_ifs <;> rfl
  rw [Finset.sum_congr rfl (fun b _ => by rw [e b]), eig_circ a m _ (fun ζ => ((alphaF ζ - (rho ξ ζ).re : ℝ) : ℂ))]
  congr 1
  rw [Complex.ofReal_sub, rowsum_zeta a (fun ζ => alphaF ζ + (rho ξ ζ).re)]
  unfold sAlpha sBeta
  simp only [pow_zero, mul_one]
  rw [sub_sub, sub_sub, sub_add, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib,
    ← Finset.sum_sub_distrib]
  congr 1
  refine Finset.sum_congr rfl fun l _ => ?_
  split_ifs
  · ring
  · push_cast; ring

theorem eig_hTT :
    ∑ b : Fin N, ((hTT N ξ a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
      = (sAlpha N m + sBeta N ξ 0 - sBeta N ξ m) * (ringW N ^ (a : ℕ)) ^ m := by
  have e : ∀ b : Fin N, ((hTT N ξ a b : ℝ) : ℂ)
      = if a = b then (((∑ k : Fin N,
          (if a = k then 0 else alphaF (zeta N a k) + (rho ξ (zeta N a k)).re)) : ℝ) : ℂ)
        else (fun ζ => ((-(alphaF ζ + (rho ξ ζ).re) : ℝ) : ℂ)) (zeta N a b) := by
    intro b; unfold hTT alphaM betaM; split_ifs <;> rfl
  rw [Finset.sum_congr rfl (fun b _ => by rw [e b]), eig_circ a m _ (fun ζ => ((-(alphaF ζ + (rho ξ ζ).re) : ℝ) : ℂ))]
  congr 1
  rw [rowsum_zeta a (fun ζ => alphaF ζ + (rho ξ ζ).re)]
  unfold sAlpha sBeta
  simp only [pow_zero, mul_one]
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun l _ => ?_
  split_ifs
  · ring
  · push_cast; ring

theorem eig_hST :
    ∑ b : Fin N, ((hST N ξ a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
      = (sGamma N ξ 0 - sGamma N ξ m) * (ringW N ^ (a : ℕ)) ^ m := by
  have e : ∀ b : Fin N, ((hST N ξ a b : ℝ) : ℂ)
      = if a = b then (((∑ k : Fin N, (if a = k then 0 else (rho ξ (zeta N a k)).im)) : ℝ) : ℂ)
        else (fun ζ => ((-(rho ξ ζ).im : ℝ) : ℂ)) (zeta N a b) := by
    intro b; unfold hST gammaM; split_ifs <;> rfl
  rw [Finset.sum_congr rfl (fun b _ => by rw [e b]), eig_circ a m _ (fun ζ => ((-(rho ξ ζ).im : ℝ) : ℂ))]
  congr 1
  rw [rowsum_zeta a (fun ζ => (rho ξ ζ).im)]
  unfold sGamma
  simp only [pow_zero, mul_one]
  rw [← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun l _ => ?_
  split_ifs
  · ring
  · push_cast; ring

theorem eig_hTS :
    ∑ b : Fin N, ((hST N ξ b a : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
      = (sGamma N ξ 0 + sGamma N ξ m) * (ringW N ^ (a : ℕ)) ^ m := by
  have e : ∀ b : Fin N, ((hST N ξ b a : ℝ) : ℂ)
      = if a = b then (((∑ k : Fin N, (if a = k then 0 else (rho ξ (zeta N a k)).im)) : ℝ) : ℂ)
        else (fun ζ => (((rho ξ ζ).im : ℝ) : ℂ)) (zeta N a b) := by
    intro b; unfold hST
    by_cases hab : a = b
    · subst hab; simp only [if_true]; unfold gammaM; rfl
    · rw [if_neg (Ne.symm hab), if_neg hab, gammaM_swap]; unfold gammaM; push_cast; ring
  rw [Finset.sum_congr rfl (fun b _ => by rw [e b]), eig_circ a m _ (fun ζ => (((rho ξ ζ).im : ℝ) : ℂ))]
  congr 1
  rw [rowsum_zeta a (fun ζ => (rho ξ ζ).im)]
  unfold sGamma
  simp only [pow_zero, mul_one]

end eig

/-! ## Evaluation -/

theorem Vj_ringOmega {N : ℕ} {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) :
    Vj ξ (ringOmega N ξ) = Vpp N ξ - 2 * ξ / (1 - ξ) ^ 2 := by
  have hne : (1 : ℝ) - ξ ≠ 0 := by linarith
  unfold Vj ringOmega kappa0 Vpp cR
  field_simp
  ring

theorem div_two_I (x : ℂ) : x / (2 * Complex.I) = -(Complex.I * x) / 2 := by
  have hI : Complex.I ≠ 0 := Complex.I_ne_zero
  rw [div_eq_iff (mul_ne_zero two_ne_zero hI)]
  ring_nf
  rw [Complex.I_sq]
  ring

/-- **Note Thm `thm:blocks` (Fourier blocks).** For `0 < ξ < 1`, `Ω = ringOmega N ξ` (note Lemma `lem:re`),
`1 ≤ m ≤ N−1`, and `v_b = (w^b)^m = w^{mb}`: the Hessian of `F = H − ΩJ` at the ring satisfies
`Σ_b ∂²F/∂s_a∂s_b v_b = a_m v_a`, `Σ_b ∂²F/∂θ_a∂θ_b v_b = b_m v_a`, `Σ_b ∂²F/∂s_a∂θ_b v_b = i g_m v_a`
and `Σ_b ∂²F/∂θ_a∂s_b v_b = −i g_m v_a`. -/
theorem blocks_eigen {N : ℕ} {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (m : ℕ) (hm1 : 1 ≤ m) (hmN : m < N)
    (a : Fin N) :
    ∑ b : Fin N, ((hessSS N ξ (ringOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
        = (aBlock N m ξ : ℂ) * (ringW N ^ (a : ℕ)) ^ m
      ∧ ∑ b : Fin N, ((hessTT N ξ (ringOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
        = (bBlock N m ξ : ℂ) * (ringW N ^ (a : ℕ)) ^ m
      ∧ ∑ b : Fin N, ((hessST N ξ (ringOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
        = Complex.I * (gBlock N m ξ : ℂ) * (ringW N ^ (a : ℕ)) ^ m
      ∧ ∑ b : Fin N, ((hessTS N ξ (ringOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
        = -Complex.I * (gBlock N m ξ : ℂ) * (ringW N ^ (a : ℕ)) ^ m := by
  have hN : 0 < N := by omega
  have : NeZero N := ⟨hN.ne'⟩
  have hX : ξ ^ N ≠ 1 := (pow_lt_one₀ h0.le h1 hN.ne').ne
  have hW0 := W_zero_eq N ξ hX
  have hsB0 := sBeta_eq (N := N) h0.le h1 0 (Nat.zero_le N)
  have hsG0 := sGamma_eq (N := N) h0.le h1 0 (Nat.zero_le N)
  rw [Nat.sub_zero, hW0] at hsB0 hsG0
  have hsBm := sBeta_eq (N := N) h0.le h1 m hmN.le
  have hsGm := sGamma_eq (N := N) h0.le h1 m hmN.le
  have hsA := sAlpha_eq (N := N) m hmN.le
  refine ⟨?_, ?_, ?_, ?_⟩
  · simp_rw [hessSS_eq hN h0 h1]
    rw [eig_hSS, hsA, hsB0, hsBm, Vj_ringOmega h0 h1,
      ← aSum_eq_aBlock N m hmN.le ξ h1.ne hX]
    unfold aSum
    push_cast
    ring
  · simp_rw [hessTT_eq hN h0 h1]
    rw [eig_hTT, hsA, hsB0, hsBm]
    unfold bBlock
    push_cast
    ring
  · simp_rw [hessST_eq hN h0 h1]
    rw [eig_hST, hsG0, hsGm, div_two_I, div_two_I]
    unfold gBlock
    push_cast
    ring
  · simp_rw [hessTS_eq hN h0 h1]
    rw [eig_hTS, hsG0, hsGm, div_two_I, div_two_I]
    unfold gBlock
    push_cast
    ring

/-- Mode `0` (note Thm `thm:blocks`: `b_0 = g_0 = 0`): the `θθ`, `sθ` and `θs` blocks annihilate the
constant vector. -/
theorem blocks_eigen_zero {N : ℕ} {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ) (hN : 0 < N)
    (a : Fin N) :
    ∑ b : Fin N, hessTT N ξ Ω a b = 0 ∧ ∑ b : Fin N, hessST N ξ Ω a b = 0
      ∧ ∑ b : Fin N, hessTS N ξ Ω a b = 0 := by
  have : NeZero N := ⟨hN.ne'⟩
  have hX : ξ ^ N ≠ 1 := (pow_lt_one₀ h0.le h1 hN.ne').ne
  have hsG0 := sGamma_eq (N := N) h0.le h1 0 (Nat.zero_le N)
  rw [Nat.sub_zero, W_zero_eq N ξ hX, sub_self, Complex.ofReal_zero, zero_div] at hsG0
  have hsA := sAlpha_eq (N := N) 0 (Nat.zero_le N)
  simp only [CharP.cast_eq_zero, zero_mul, zero_div] at hsA
  refine ⟨?_, ?_, ?_⟩
  · have := eig_hTT (N := N) ξ a 0
    simp only [pow_zero, mul_one, hsA, zero_add, sub_self] at this
    simp_rw [hessTT_eq hN h0 h1]
    exact_mod_cast this
  · have := eig_hST (N := N) ξ a 0
    simp only [pow_zero, mul_one, sub_self] at this
    simp_rw [hessST_eq hN h0 h1]
    exact_mod_cast this
  · have := eig_hTS (N := N) ξ a 0
    simp only [pow_zero, mul_one, hsG0, add_zero] at this
    simp_rw [hessTS_eq hN h0 h1]
    exact_mod_cast this

/-- Hermitian form from an eigenvector identity (orthogonality of characters). -/
theorem herm_of_eig {N : ℕ} (hN : 0 < N) {H : Fin N → Fin N → ℝ} {μ : ℂ} {m : ℕ}
    (hH : ∀ a : Fin N, ∑ b : Fin N, ((H a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
      = μ * (ringW N ^ (a : ℕ)) ^ m) (m' : ℕ) (hm : m < N) (hm' : m' < N) :
    ∑ a : Fin N, ∑ b : Fin N, conj ((ringW N ^ (a : ℕ)) ^ m') * ((H a b : ℝ) : ℂ)
        * (ringW N ^ (b : ℕ)) ^ m
      = if m = m' then (N : ℂ) * μ else 0 := by
  simp_rw [mul_assoc, ← Finset.mul_sum, hH]
  have e : ∀ a : Fin N, conj ((ringW N ^ (a : ℕ)) ^ m') * (μ * (ringW N ^ (a : ℕ)) ^ m)
      = μ * ((ringW N ^ (a : ℕ)) ^ m * conj (ringW N ^ (a : ℕ)) ^ m') := by
    intro a; rw [map_pow]; ring
  rw [Finset.sum_congr rfl (fun a _ => e a), ← Finset.mul_sum, rootSum hN, Nat.mod_eq_of_lt hm,
    Nat.mod_eq_of_lt hm']
  split_ifs <;> ring

/-- **Note Thm `thm:blocks`, Hermitian form.** For `1 ≤ m ≤ N−1`, `0 ≤ m' ≤ N−1`, `v_b = w^{mb}`,
`v'_a = w^{m'a}`: `Σ_{a,b} conj(v'_a) ∂²F/∂s_a∂s_b v_b = N a_m [m = m']`, and likewise for the
`θθ` (`N b_m`), `sθ` (`N i g_m`) and `θs` (`−N i g_m`) blocks. -/
theorem blocks_hermitian {N : ℕ} {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (m m' : ℕ) (hm1 : 1 ≤ m)
    (hmN : m < N) (hm' : m' < N) :
    (∑ a : Fin N, ∑ b : Fin N, conj ((ringW N ^ (a : ℕ)) ^ m')
        * ((hessSS N ξ (ringOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
      = if m = m' then (N : ℂ) * (aBlock N m ξ : ℂ) else 0)
    ∧ (∑ a : Fin N, ∑ b : Fin N, conj ((ringW N ^ (a : ℕ)) ^ m')
        * ((hessTT N ξ (ringOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
      = if m = m' then (N : ℂ) * (bBlock N m ξ : ℂ) else 0)
    ∧ (∑ a : Fin N, ∑ b : Fin N, conj ((ringW N ^ (a : ℕ)) ^ m')
        * ((hessST N ξ (ringOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
      = if m = m' then (N : ℂ) * (Complex.I * (gBlock N m ξ : ℂ)) else 0)
    ∧ (∑ a : Fin N, ∑ b : Fin N, conj ((ringW N ^ (a : ℕ)) ^ m')
        * ((hessTS N ξ (ringOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
      = if m = m' then (N : ℂ) * (-Complex.I * (gBlock N m ξ : ℂ)) else 0) := by
  have hN : 0 < N := by omega
  refine ⟨herm_of_eig hN (fun a => (blocks_eigen h0 h1 m hm1 hmN a).1) m' hmN hm',
    herm_of_eig hN (fun a => (blocks_eigen h0 h1 m hm1 hmN a).2.1) m' hmN hm',
    herm_of_eig hN (fun a => ?_) m' hmN hm', herm_of_eig hN (fun a => ?_) m' hmN hm'⟩
  · rw [(blocks_eigen h0 h1 m hm1 hmN a).2.2.1, mul_assoc]
  · rw [(blocks_eigen h0 h1 m hm1 hmN a).2.2.2, mul_assoc]

end Havelock.Tanh
