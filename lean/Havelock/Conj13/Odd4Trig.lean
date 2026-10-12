/-
Havelock — Conjecture 13, odd `n = 2M+1`, `M ≥ 2`: the one-index Fourier sums used for the `A⁴`
coefficient (paper `conj13/paper/main.tex`, proof of Thm odd, "Sums over roots of unity").

Notation: `η_i = 2πi/n + h` (any real phase `h`). For `M ≥ 2`:

* `B1`  `Σ sin²(Mη_i) = n/2`;            `B2`  `Σ sin²(Mη_i)cos²(Mη_i) = n/8`;
* `B3`  `Σ sin⁴(Mη_i) = 3n/8`;           `B4`  `Σ sin²(Mη_i)cos η_i = −(n/4)cos(nh)`;
* `B5`  `Σ sin(Mη_i)cos(Mη_i)sin η_i = −(n/4)cos(nh)`;
* `B6`  `Σ sin²η_i = n/2`;  `B7`  `Σ cos²η_i = n/2`;  `B8`  `Σ cos η_i = 0`;
* `B9`  `Σ sin²(Mη_i)cos(Mη_i) = 0`;      `B10` `Σ sin(Mη_i) sin η_i = 0`.

Also the one-body sums `Σ c_k³ = 0`, `Σ c_k⁴ = 3n/8` (`c_k = cos(2πMk/n)`).
-/
import Havelock.Conj13.Odd4Fejer
import Havelock.Conj13.OddConstraint

open Real Finset

namespace Havelock.Conj13.Odd4

open Havelock.Conj13.OddSums Havelock.Conj13.Odd

/-- `η_i = 2πi/n + h`, `n = 2M+1`. -/
noncomputable def eta (M i : ℕ) (h : ℝ) : ℝ := 2 * π * i / (2 * M + 1 : ℕ) + h

theorem nd_lt (M k : ℕ) (h0 : 0 < k) (h1 : k < 2 * M + 1) : ¬ (2 * M + 1) ∣ k :=
  fun h => by have := Nat.le_of_dvd h0 h; omega

theorem nd_lt2 (M k : ℕ) (h0 : 2 * M + 1 < k) (h1 : k < 2 * (2 * M + 1)) : ¬ (2 * M + 1) ∣ k := by
  rintro ⟨c, hc⟩
  have hc2 : c < 2 := by
    by_contra hcon
    push Not at hcon
    have : (2 * M + 1) * 2 ≤ (2 * M + 1) * c := Nat.mul_le_mul_left _ hcon
    omega
  interval_cases c <;> omega

theorem sum_cos_eta (M k : ℕ) (hk : ¬ (2 * M + 1) ∣ k) (h : ℝ) :
    ∑ i ∈ range (2 * M + 1), Real.cos ((k : ℝ) * eta M i h) = 0 := by
  rw [← sum_cos_add_eq_zero (2 * M + 1) k (by omega) hk ((k : ℝ) * h)]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  unfold eta; congr 1; ring

theorem sum_cos_eta_n (M : ℕ) (h : ℝ) :
    ∑ i ∈ range (2 * M + 1), Real.cos (((2 * M + 1 : ℕ) : ℝ) * eta M i h)
      = ((2 * M + 1 : ℕ) : ℝ) * Real.cos (((2 * M + 1 : ℕ) : ℝ) * h) := by
  have hn : ((2 * M + 1 : ℕ) : ℝ) ≠ 0 := by positivity
  have e : ∀ i ∈ range (2 * M + 1), Real.cos (((2 * M + 1 : ℕ) : ℝ) * eta M i h)
      = Real.cos (((2 * M + 1 : ℕ) : ℝ) * h) := by
    intro i _
    rw [show ((2 * M + 1 : ℕ) : ℝ) * eta M i h = ((2 * M + 1 : ℕ) : ℝ) * h + i * (2 * π) by
      unfold eta; field_simp; ring, Real.cos_add_nat_mul_two_pi]
  rw [Finset.sum_congr rfl e, Finset.sum_const, card_range, nsmul_eq_mul]

/-- Casting helpers. -/
theorem cast2M (M : ℕ) (y : ℝ) : ((2 * M : ℕ) : ℝ) * y = 2 * ((M : ℝ) * y) := by push_cast; ring
theorem cast4M (M : ℕ) (y : ℝ) : ((4 * M : ℕ) : ℝ) * y = 2 * (2 * ((M : ℝ) * y)) := by push_cast; ring
theorem cast3M (M : ℕ) (y : ℝ) : ((3 * M : ℕ) : ℝ) * y = 3 * ((M : ℝ) * y) := by push_cast; ring
theorem castn (M : ℕ) (y : ℝ) : ((2 * M + 1 : ℕ) : ℝ) * y = 2 * ((M : ℝ) * y) + y := by push_cast; ring
theorem cast2Mm1 (M : ℕ) (hM : 1 ≤ M) (y : ℝ) : ((2 * M - 1 : ℕ) : ℝ) * y = 2 * ((M : ℝ) * y) - y := by
  rw [Nat.cast_sub (by omega)]; push_cast; ring
theorem castMm1 (M : ℕ) (hM : 1 ≤ M) (y : ℝ) : ((M - 1 : ℕ) : ℝ) * y = (M : ℝ) * y - y := by
  rw [Nat.cast_sub hM]; push_cast; ring
theorem castMp1 (M : ℕ) (y : ℝ) : ((M + 1 : ℕ) : ℝ) * y = (M : ℝ) * y + y := by push_cast; ring

section
variable (M : ℕ) (hM : 2 ≤ M) (h : ℝ)
include hM

theorem B1 : ∑ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h) ^ 2
    = ((2 * M + 1 : ℕ) : ℝ) / 2 := by
  have e : ∀ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h) ^ 2
      = 1 / 2 - 1 / 2 * Real.cos (((2 * M : ℕ) : ℝ) * eta M i h) := by
    intro i _
    rw [cast2M, Real.cos_two_mul]
    linear_combination Real.sin_sq_add_cos_sq ((M : ℝ) * eta M i h)
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]
  rw [sum_cos_eta M (2 * M) (nd_lt M _ (by omega) (by omega))]
  simp; ring

theorem B2 : ∑ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h) ^ 2
    * Real.cos ((M : ℝ) * eta M i h) ^ 2 = ((2 * M + 1 : ℕ) : ℝ) / 8 := by
  have e : ∀ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h) ^ 2
      * Real.cos ((M : ℝ) * eta M i h) ^ 2
      = 1 / 8 - 1 / 8 * Real.cos (((4 * M : ℕ) : ℝ) * eta M i h) := by
    intro i _
    rw [cast4M, Real.cos_two_mul, Real.cos_two_mul]
    linear_combination (Real.cos ((M : ℝ) * eta M i h) ^ 2) * Real.sin_sq_add_cos_sq ((M : ℝ) * eta M i h)
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]
  rw [sum_cos_eta M (4 * M) (nd_lt2 M _ (by omega) (by omega))]
  simp; ring

theorem B3 : ∑ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h) ^ 4
    = 3 * ((2 * M + 1 : ℕ) : ℝ) / 8 := by
  have e : ∀ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h) ^ 4
      = 3 / 8 - 1 / 2 * Real.cos (((2 * M : ℕ) : ℝ) * eta M i h)
        + 1 / 8 * Real.cos (((4 * M : ℕ) : ℝ) * eta M i h) := by
    intro i _
    rw [cast2M, cast4M, Real.cos_two_mul, Real.cos_two_mul (2 * _), Real.cos_two_mul]
    linear_combination (Real.sin ((M : ℝ) * eta M i h) ^ 2 + 1 - Real.cos ((M : ℝ) * eta M i h) ^ 2)
      * Real.sin_sq_add_cos_sq ((M : ℝ) * eta M i h)
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]
  rw [sum_cos_eta M (2 * M) (nd_lt M _ (by omega) (by omega)),
    sum_cos_eta M (4 * M) (nd_lt2 M _ (by omega) (by omega))]
  simp; ring

theorem B4 : ∑ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h) ^ 2 * Real.cos (eta M i h)
    = -(((2 * M + 1 : ℕ) : ℝ) / 4) * Real.cos (((2 * M + 1 : ℕ) : ℝ) * h) := by
  have e : ∀ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h) ^ 2 * Real.cos (eta M i h)
      = 1 / 2 * Real.cos (((1 : ℕ) : ℝ) * eta M i h)
        - 1 / 4 * Real.cos (((2 * M + 1 : ℕ) : ℝ) * eta M i h)
        - 1 / 4 * Real.cos (((2 * M - 1 : ℕ) : ℝ) * eta M i h) := by
    intro i _
    rw [castn, cast2Mm1 M (by omega), Nat.cast_one, one_mul, Real.cos_add, Real.cos_sub,
      Real.cos_two_mul, Real.sin_two_mul]
    linear_combination (Real.cos (eta M i h)) * Real.sin_sq_add_cos_sq ((M : ℝ) * eta M i h)
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]
  rw [sum_cos_eta M 1 (nd_lt M _ (by omega) (by omega)),
    sum_cos_eta M (2 * M - 1) (nd_lt M _ (by omega) (by omega)), sum_cos_eta_n]
  ring

theorem B5 : ∑ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h)
    * Real.cos ((M : ℝ) * eta M i h) * Real.sin (eta M i h)
    = -(((2 * M + 1 : ℕ) : ℝ) / 4) * Real.cos (((2 * M + 1 : ℕ) : ℝ) * h) := by
  have e : ∀ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h)
      * Real.cos ((M : ℝ) * eta M i h) * Real.sin (eta M i h)
      = 1 / 4 * Real.cos (((2 * M - 1 : ℕ) : ℝ) * eta M i h)
        - 1 / 4 * Real.cos (((2 * M + 1 : ℕ) : ℝ) * eta M i h) := by
    intro i _
    rw [castn, cast2Mm1 M (by omega), Real.cos_add, Real.cos_sub, Real.sin_two_mul]
    ring
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]
  rw [sum_cos_eta M (2 * M - 1) (nd_lt M _ (by omega) (by omega)), sum_cos_eta_n]
  ring

theorem B6 : ∑ i ∈ range (2 * M + 1), Real.sin (eta M i h) ^ 2 = ((2 * M + 1 : ℕ) : ℝ) / 2 := by
  have e : ∀ i ∈ range (2 * M + 1), Real.sin (eta M i h) ^ 2
      = 1 / 2 - 1 / 2 * Real.cos (((2 : ℕ) : ℝ) * eta M i h) := by
    intro i _
    rw [Nat.cast_ofNat, Real.cos_two_mul]
    linear_combination Real.sin_sq_add_cos_sq (eta M i h)
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]
  rw [sum_cos_eta M 2 (nd_lt M _ (by omega) (by omega))]
  simp; ring

theorem B7 : ∑ i ∈ range (2 * M + 1), Real.cos (eta M i h) ^ 2 = ((2 * M + 1 : ℕ) : ℝ) / 2 := by
  have e : ∀ i ∈ range (2 * M + 1), Real.cos (eta M i h) ^ 2
      = 1 / 2 + 1 / 2 * Real.cos (((2 : ℕ) : ℝ) * eta M i h) := by
    intro i _
    rw [Nat.cast_ofNat, Real.cos_two_mul]
    ring
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]
  rw [sum_cos_eta M 2 (nd_lt M _ (by omega) (by omega))]
  simp; ring

theorem B8 : ∑ i ∈ range (2 * M + 1), Real.cos (eta M i h) = 0 := by
  have := sum_cos_eta M 1 (nd_lt M _ (by omega) (by omega)) h
  simpa using this

theorem B9 : ∑ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h) ^ 2
    * Real.cos ((M : ℝ) * eta M i h) = 0 := by
  have e : ∀ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h) ^ 2
      * Real.cos ((M : ℝ) * eta M i h)
      = 1 / 4 * Real.cos (((M : ℕ) : ℝ) * eta M i h) - 1 / 4 * Real.cos (((3 * M : ℕ) : ℝ) * eta M i h) := by
    intro i _
    rw [cast3M, Real.cos_three_mul]
    linear_combination (Real.cos ((M : ℝ) * eta M i h)) * Real.sin_sq_add_cos_sq ((M : ℝ) * eta M i h)
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]
  rw [sum_cos_eta M M (nd_lt M _ (by omega) (by omega)),
    sum_cos_eta M (3 * M) (nd_lt2 M _ (by omega) (by omega))]
  simp

theorem B10 : ∑ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h) * Real.sin (eta M i h) = 0 := by
  have e : ∀ i ∈ range (2 * M + 1), Real.sin ((M : ℝ) * eta M i h) * Real.sin (eta M i h)
      = 1 / 2 * Real.cos (((M - 1 : ℕ) : ℝ) * eta M i h) - 1 / 2 * Real.cos (((M + 1 : ℕ) : ℝ) * eta M i h) := by
    intro i _
    rw [castMm1 M (by omega), castMp1, Real.cos_add, Real.cos_sub]
    ring
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]
  rw [sum_cos_eta M (M - 1) (nd_lt M _ (by omega) (by omega)),
    sum_cos_eta M (M + 1) (nd_lt M _ (by omega) (by omega))]
  simp

/-- `Σ_k c_k³ = 0`. -/
theorem sum_cc_cube : ∑ k ∈ range (2 * M + 1), cc M k ^ 3 = 0 := by
  have := B9 M hM 0
  have e : ∀ k ∈ range (2 * M + 1), cc M k ^ 3 = Real.cos ((M : ℝ) * eta M k 0)
      - Real.sin ((M : ℝ) * eta M k 0) ^ 2 * Real.cos ((M : ℝ) * eta M k 0) := by
    intro k _
    have hc : cc M k = Real.cos ((M : ℝ) * eta M k 0) := by unfold cc eta; congr 1; ring
    rw [hc]
    linear_combination (Real.cos ((M : ℝ) * eta M k 0)) * Real.sin_sq_add_cos_sq ((M : ℝ) * eta M k 0)
  rw [Finset.sum_congr rfl e, Finset.sum_sub_distrib, this, sub_zero]
  have := sum_cos_eta M M (nd_lt M _ (by omega) (by omega)) 0
  simpa using this

/-- `Σ_k c_k⁴ = 3n/8`. -/
theorem sum_cc_four : ∑ k ∈ range (2 * M + 1), cc M k ^ 4 = 3 * ((2 * M + 1 : ℕ) : ℝ) / 8 := by
  have h3 := B3 M hM 0
  have h1 := B1 M hM 0
  have e : ∀ k ∈ range (2 * M + 1), cc M k ^ 4 = 1 - 2 * Real.sin ((M : ℝ) * eta M k 0) ^ 2
      + Real.sin ((M : ℝ) * eta M k 0) ^ 4 := by
    intro k _
    have hc : cc M k = Real.cos ((M : ℝ) * eta M k 0) := by unfold cc eta; congr 1; ring
    rw [hc]
    linear_combination (Real.cos ((M : ℝ) * eta M k 0) ^ 2 - Real.sin ((M : ℝ) * eta M k 0) ^ 2 + 1)
      * Real.sin_sq_add_cos_sq ((M : ℝ) * eta M k 0)
  rw [Finset.sum_congr rfl e, Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum, h1, h3]
  simp; ring

end

end Havelock.Conj13.Odd4
