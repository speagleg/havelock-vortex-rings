import Havelock.Tanh.LargeN.Symbol

/-! # Tanh model, large `N`: tangential stiffness (note Lemma 12)

`β(x,τ) > 0` for `|x| < 1`, `τ > 0`, hence `b_m > 0` for every `N ≥ 2`, `1 ≤ m ≤ N−1`, `0 < ξ < 1`.

With `u = (1+x)τ/2`, `v = (1−x)τ/2` (the note's `a/2`, `b/2`) one has
`β = [u v sinh²(u+v) − (u+v)(v sinh²u + u sinh²v)] / (2 (u+v)² sinh²(u+v))`, and positivity of the
numerator is the superadditivity of `t ↦ sinh²t / t`, which follows from `sinh t / t` being strictly
increasing (the note's route). Here that monotonicity is proved as
`(u+v) sinh u < u sinh(u+v)`, from `sinh u < u cosh u` and `v < sinh v`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real

/-- `sinh u < u cosh u` for `u > 0` (i.e. `tanh u < u`). -/
theorem sinh_lt_mul_cosh {u : ℝ} (hu : 0 < u) : sinh u < u * cosh u := by
  set h : ℝ → ℝ := fun t => t * cosh t - sinh t with h_def
  have hd : ∀ t, HasDerivAt h (t * sinh t) t := by
    intro t
    have := ((hasDerivAt_id t).mul (hasDerivAt_cosh t)).sub (hasDerivAt_sinh t)
    convert this using 1
    · rw [h_def]; rfl
    · rw [id_eq]; ring
  have hmono : StrictMonoOn h (Set.Ici 0) := by
    refine strictMonoOn_of_deriv_pos (convex_Ici 0) ?_ ?_
    · exact (continuous_id.mul continuous_cosh).sub continuous_sinh |>.continuousOn
    · intro t ht
      rw [interior_Ici] at ht
      have ht' : 0 < t := ht
      rw [(hd t).deriv]
      exact mul_pos ht' (sinh_pos_iff.mpr ht')
  have := hmono (Set.mem_Ici.mpr le_rfl) (Set.mem_Ici.mpr hu.le) hu
  simp only [h, cosh_zero, sinh_zero, zero_mul, sub_zero] at this
  linarith

/-- `sinh t / t` is strictly increasing, in the form `(u+v) sinh u < u sinh(u+v)` (proof of note Lemma 12). -/
theorem sinh_div_strict {u v : ℝ} (hu : 0 < u) (hv : 0 < v) :
    (u + v) * sinh u < u * sinh (u + v) := by
  rw [sinh_add]
  have h1 := sinh_lt_mul_cosh hu
  have h2 : v < sinh v := self_lt_sinh_iff.mpr hv
  have h3 : 0 < sinh u := sinh_pos_iff.mpr hu
  have h4 : 1 ≤ cosh v := one_le_cosh v
  have h5 : 0 < sinh v := sinh_pos_iff.mpr hv
  nlinarith [mul_pos h5 (sub_pos.mpr h1), mul_pos h3 (sub_pos.mpr h2),
    mul_nonneg (mul_pos hu h3).le (sub_nonneg.mpr h4)]

/-- Superadditivity of `t ↦ sinh²t/t`, cleared of denominators (proof of note Lemma 12). -/
theorem sinh_sq_superadd {u v : ℝ} (hu : 0 < u) (hv : 0 < v) :
    (u + v) * (v * sinh u ^ 2 + u * sinh v ^ 2) < u * v * sinh (u + v) ^ 2 := by
  have hA := sinh_div_strict hu hv
  have hB := sinh_div_strict hv hu
  rw [add_comm v u] at hB
  have h3 : 0 < sinh u := sinh_pos_iff.mpr hu
  have h5 : 0 < sinh v := sinh_pos_iff.mpr hv
  have hs : 0 < u + v := by linarith
  have hA2 : ((u + v) * sinh u) ^ 2 < (u * sinh (u + v)) ^ 2 :=
    pow_lt_pow_left₀ hA (by positivity) two_ne_zero
  have hB2 : ((u + v) * sinh v) ^ 2 < (v * sinh (u + v)) ^ 2 :=
    pow_lt_pow_left₀ hB (by positivity) two_ne_zero
  have key : (u + v) * ((u + v) * (v * sinh u ^ 2 + u * sinh v ^ 2))
      < (u + v) * (u * v * sinh (u + v) ^ 2) := by
    nlinarith [mul_lt_mul_of_pos_left hA2 hv, mul_lt_mul_of_pos_left hB2 hu]
  exact lt_of_mul_lt_mul_left key hs.le

/-- `β` in the variables `u = (1+x)τ/2`, `v = (1−x)τ/2` (note Lemma 12, with `a = 2u`, `b = 2v`). -/
theorem beta_uv {u v : ℝ} (hu : 0 < u) (hv : 0 < v) :
    beta ((u - v) / (u + v)) (u + v) =
      (u * v * sinh (u + v) ^ 2 - (u + v) * (v * sinh u ^ 2 + u * sinh v ^ 2))
        / (2 * (u + v) ^ 2 * sinh (u + v) ^ 2) := by
  have hs : u + v ≠ 0 := by positivity
  unfold beta Q
  rw [show (u - v) / (u + v) * (u + v) = u - v by field_simp]
  simp only [cosh_eq, sinh_eq, exp_neg, exp_add, exp_sub]
  have hA1 : 1 < exp u := one_lt_exp_iff.mpr hu
  have hB1 : 1 < exp v := one_lt_exp_iff.mpr hv
  generalize exp u = A at *
  generalize exp v = B at *
  have hA0 : A ≠ 0 := by positivity
  have hB0 : B ≠ 0 := by positivity
  have hAB : 1 < A * B := by nlinarith
  have hAB2 : A ^ 2 * B ^ 2 - 1 ≠ 0 := by nlinarith
  have hAB2' : A * B * (A * B) - 1 ≠ 0 := by nlinarith
  have hAB3 : (A * B) ^ 2 - 1 ≠ 0 := by nlinarith
  field_simp
  ring

/-- Tangential stiffness: `β(x,τ) > 0` for `|x| < 1`, `τ > 0` (note Lemma 12). -/
theorem beta_pos {x τ : ℝ} (hx : |x| < 1) (hτ : 0 < τ) : 0 < beta x τ := by
  obtain ⟨hx1, hx2⟩ := abs_lt.mp hx
  set u := (1 + x) * τ / 2 with hu_def
  set v := (1 - x) * τ / 2 with hv_def
  have hu : 0 < u := by rw [hu_def]; nlinarith
  have hv : 0 < v := by rw [hv_def]; nlinarith
  have hτuv : u + v = τ := by rw [hu_def, hv_def]; ring
  have hxuv : (u - v) / (u + v) = x := by
    rw [hτuv, hu_def, hv_def]; field_simp; ring
  rw [← hxuv, ← hτuv, beta_uv hu hv]
  have hS : 0 < sinh (u + v) := sinh_pos_iff.mpr (by linarith)
  exact div_pos (by linarith [sinh_sq_superadd hu hv]) (by positivity)

/-- `b_m > 0` for every `N ≥ 2`, `1 ≤ m ≤ N−1` and `0 < ξ < 1` (note Lemma 12, "Hence"). -/
theorem bBlock_pos {N m : ℕ} {ξ : ℝ} (hN : 2 ≤ N) (hm1 : 1 ≤ m) (hmN : m < N) (hξ0 : 0 < ξ)
    (hξ1 : ξ < 1) : 0 < bBlock N m ξ := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  set τ := -((N : ℝ) / 2) * log ξ with hτ_def
  have hlog : log ξ < 0 := log_neg hξ0 hξ1
  have hτ : 0 < τ := by rw [hτ_def]; nlinarith
  have hξτ : xiOf N τ = ξ := by
    unfold xiOf
    rw [hτ_def, show -2 * (-((N : ℝ) / 2) * log ξ) / N = log ξ by field_simp, exp_log hξ0]
  have hx : |xOf N m| < 1 := by
    unfold xOf
    have hm1' : (1 : ℝ) ≤ m := by exact_mod_cast hm1
    have hmN' : (m : ℝ) < N := by exact_mod_cast hmN
    have h1 : 0 < 2 * (m : ℝ) / N := by positivity
    have h2 : 2 * (m : ℝ) / N < 2 := by rw [div_lt_iff₀ hNpos]; linarith
    rw [abs_lt]; constructor <;> linarith
  rw [← hξτ, bBlock_xiOf hm1 hmN hτ]
  exact mul_pos (by positivity) (beta_pos hx hτ)

end Havelock.Tanh.LargeN
