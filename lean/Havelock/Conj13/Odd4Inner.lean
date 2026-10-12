/-
Havelock — Conjecture 13, odd `n = 2M+1`, `M ≥ 2`: the inner (fixed-offset) sums of the pair
coefficients of the `A⁴` expansion (paper `conj13/paper/main.tex`, proof of Thm odd).

For a real offset `x` put `h = πx/n` (half the angle `θ = 2πx/n`), `c_i = cos(2πMi/n)`,
`c_j = cos(2πM(i+x)/n)`, `C_i, C_j, S_i, S_j` the same with frequency `1` (cos, sin). Half-angle
factorisation: `c_i − c_j = 2 sin(Mη_i) sin(Mh)`, `c_i + c_j = 2 cos(Mη_i) cos(Mh)`,
`c_ic_j = cos²(Mη_i) − sin²(Mh)`, `C_i − C_j = 2 sin η_i sin h`, `C_i + C_j = 2 cos η_i cos h`,
`S_j − S_i = 2 cos η_i sin h`, with `η_i = 2πi/n + h`. With `S = sin(Mh)`, `Cm = cos(Mh)`,
`s₁ = sin h`, `c₁ = cos h`, `Z = cos(nh)`, the sums over `i ∈ [0, n)` are:

* `inner_a`  `Σ (c_i−c_j)²(C_i+C_j) = −2n S² c₁ Z`
* `inner_b`  `Σ (c_i−c_j)²(c_i²+c_j²) = 2nS²Cm² − nS² + 4nS⁴`
* `inner_c`  `Σ (c_i−c_j)² c_ic_j = nS²/2 − 2nS⁴`
* `inner_d`  `Σ (c_i−c_j)(C_i−C_j)(c_i+c_j) = −2n S Cm s₁ Z`
* `inner_e`  `Σ (C_i−C_j)² = 2n s₁²`,  `inner_s`  `Σ (S_j−S_i)² = 2n s₁²`
* `inner_f`  `Σ (c_i−c_j)⁴ = 6nS⁴`
* `inner_g`  `Σ (c_i−c_j)²(S_j−S_i) = −2n S² s₁ Z`
* `inner_lin`, `inner_3a`, `inner_3b`: `Σ (S_j−S_i) = Σ (c_i−c_j)²(c_i+c_j) = Σ (c_i−c_j)(C_i−C_j) = 0`.
-/
import Havelock.Conj13.Odd4Trig

open Real Finset

namespace Havelock.Conj13.Odd4

open Havelock.Conj13.OddSums Havelock.Conj13.Odd

/-! ### Half-angle factorisations -/

theorem hdiff (u v : ℝ) : Real.cos u - Real.cos (u + v) = 2 * Real.sin (u + v / 2) * Real.sin (v / 2) := by
  have h := Real.cos_sub (u + v / 2) (v / 2)
  have h2 := Real.cos_add (u + v / 2) (v / 2)
  rw [show u + v / 2 - v / 2 = u by ring] at h
  rw [show u + v / 2 + v / 2 = u + v by ring] at h2
  rw [h, h2]; ring

theorem hsum (u v : ℝ) : Real.cos u + Real.cos (u + v) = 2 * Real.cos (u + v / 2) * Real.cos (v / 2) := by
  have h := Real.cos_sub (u + v / 2) (v / 2)
  have h2 := Real.cos_add (u + v / 2) (v / 2)
  rw [show u + v / 2 - v / 2 = u by ring] at h
  rw [show u + v / 2 + v / 2 = u + v by ring] at h2
  rw [h, h2]; ring

theorem hprod (u v : ℝ) :
    Real.cos u * Real.cos (u + v) = Real.cos (u + v / 2) ^ 2 - Real.sin (v / 2) ^ 2 := by
  have h := Real.cos_sub (u + v / 2) (v / 2)
  have h2 := Real.cos_add (u + v / 2) (v / 2)
  rw [show u + v / 2 - v / 2 = u by ring] at h
  rw [show u + v / 2 + v / 2 = u + v by ring] at h2
  rw [h, h2]
  linear_combination (Real.cos (u + v / 2) ^ 2) * Real.sin_sq_add_cos_sq (v / 2)
    - Real.sin (v / 2) ^ 2 * Real.sin_sq_add_cos_sq (u + v / 2)

theorem hsdiff (u v : ℝ) : Real.sin (u + v) - Real.sin u = 2 * Real.cos (u + v / 2) * Real.sin (v / 2) := by
  have h := Real.sin_sub (u + v / 2) (v / 2)
  have h2 := Real.sin_add (u + v / 2) (v / 2)
  rw [show u + v / 2 - v / 2 = u by ring] at h
  rw [show u + v / 2 + v / 2 = u + v by ring] at h2
  rw [h, h2]; ring

/-- Half the offset angle, `πx/n`. -/
noncomputable def hh (M : ℕ) (x : ℝ) : ℝ := π * x / (2 * M + 1 : ℕ)

section
variable (M i : ℕ) (x : ℝ)

theorem argM : 2 * π * (M : ℝ) * (i : ℝ) / ((2 * M + 1 : ℕ) : ℝ)
      + 2 * π * (M : ℝ) * x / ((2 * M + 1 : ℕ) : ℝ) / 2 = (M : ℝ) * eta M i (hh M x) := by
  unfold eta hh; ring

theorem arg1 : 2 * π * (i : ℝ) / ((2 * M + 1 : ℕ) : ℝ)
      + 2 * π * x / ((2 * M + 1 : ℕ) : ℝ) / 2 = eta M i (hh M x) := by
  unfold eta hh; ring

theorem cc_diff : cc M i - cc M (i + x)
    = 2 * Real.sin ((M : ℝ) * eta M i (hh M x)) * Real.sin ((M : ℝ) * hh M x) := by
  have := hdiff (2 * π * M * i / (2 * M + 1 : ℕ)) (2 * π * M * x / (2 * M + 1 : ℕ))
  rw [argM, show 2 * π * (M : ℝ) * x / ((2 * M + 1 : ℕ) : ℝ) / 2 = (M : ℝ) * hh M x by unfold hh; ring]
    at this
  unfold cc
  rw [← this]; congr 2; ring

theorem cc_sum : cc M i + cc M (i + x)
    = 2 * Real.cos ((M : ℝ) * eta M i (hh M x)) * Real.cos ((M : ℝ) * hh M x) := by
  have := hsum (2 * π * M * i / (2 * M + 1 : ℕ)) (2 * π * M * x / (2 * M + 1 : ℕ))
  rw [argM, show 2 * π * (M : ℝ) * x / ((2 * M + 1 : ℕ) : ℝ) / 2 = (M : ℝ) * hh M x by unfold hh; ring]
    at this
  unfold cc
  rw [← this]; congr 2; ring

theorem cc_prod : cc M i * cc M (i + x)
    = Real.cos ((M : ℝ) * eta M i (hh M x)) ^ 2 - Real.sin ((M : ℝ) * hh M x) ^ 2 := by
  have := hprod (2 * π * M * i / (2 * M + 1 : ℕ)) (2 * π * M * x / (2 * M + 1 : ℕ))
  rw [argM, show 2 * π * (M : ℝ) * x / ((2 * M + 1 : ℕ) : ℝ) / 2 = (M : ℝ) * hh M x by unfold hh; ring]
    at this
  unfold cc
  rw [← this]; congr 2; ring

theorem CC_diff : CC M i - CC M (i + x) = 2 * Real.sin (eta M i (hh M x)) * Real.sin (hh M x) := by
  have := hdiff (2 * π * i / (2 * M + 1 : ℕ)) (2 * π * x / (2 * M + 1 : ℕ))
  rw [arg1, show 2 * π * x / ((2 * M + 1 : ℕ) : ℝ) / 2 = hh M x by unfold hh; ring] at this
  unfold CC
  rw [← this]; congr 2; ring

theorem CC_sum : CC M i + CC M (i + x) = 2 * Real.cos (eta M i (hh M x)) * Real.cos (hh M x) := by
  have := hsum (2 * π * i / (2 * M + 1 : ℕ)) (2 * π * x / (2 * M + 1 : ℕ))
  rw [arg1, show 2 * π * x / ((2 * M + 1 : ℕ) : ℝ) / 2 = hh M x by unfold hh; ring] at this
  unfold CC
  rw [← this]; congr 2; ring

theorem SS_diff : SS M (i + x) - SS M i = 2 * Real.cos (eta M i (hh M x)) * Real.sin (hh M x) := by
  have := hsdiff (2 * π * i / (2 * M + 1 : ℕ)) (2 * π * x / (2 * M + 1 : ℕ))
  rw [arg1, show 2 * π * x / ((2 * M + 1 : ℕ) : ℝ) / 2 = hh M x by unfold hh; ring] at this
  unfold SS
  rw [← this]; congr 2; ring

end

/-! ### Inner sums -/

section
variable (M : ℕ) (hM : 2 ≤ M) (x : ℝ)
include hM

theorem inner_a : ∑ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) ^ 2 * (CC M i + CC M (i + x))
    = -2 * ((2 * M + 1 : ℕ) : ℝ) * Real.sin ((M : ℝ) * hh M x) ^ 2 * Real.cos (hh M x)
      * Real.cos (((2 * M + 1 : ℕ) : ℝ) * hh M x) := by
  have e : ∀ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) ^ 2 * (CC M i + CC M (i + x))
      = (8 * Real.sin ((M : ℝ) * hh M x) ^ 2 * Real.cos (hh M x))
        * (Real.sin ((M : ℝ) * eta M i (hh M x)) ^ 2 * Real.cos (eta M i (hh M x))) := by
    intro i _; rw [cc_diff, CC_sum]; ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum, B4 M hM]; ring

theorem inner_b : ∑ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) ^ 2 * (cc M i ^ 2 + cc M (i + x) ^ 2)
    = 2 * ((2 * M + 1 : ℕ) : ℝ) * Real.sin ((M : ℝ) * hh M x) ^ 2 * Real.cos ((M : ℝ) * hh M x) ^ 2
      - ((2 * M + 1 : ℕ) : ℝ) * Real.sin ((M : ℝ) * hh M x) ^ 2
      + 4 * ((2 * M + 1 : ℕ) : ℝ) * Real.sin ((M : ℝ) * hh M x) ^ 4 := by
  have e : ∀ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) ^ 2 * (cc M i ^ 2 + cc M (i + x) ^ 2)
      = (16 * Real.sin ((M : ℝ) * hh M x) ^ 2 * Real.cos ((M : ℝ) * hh M x) ^ 2
          - 8 * Real.sin ((M : ℝ) * hh M x) ^ 2)
        * (Real.sin ((M : ℝ) * eta M i (hh M x)) ^ 2 * Real.cos ((M : ℝ) * eta M i (hh M x)) ^ 2)
        + (8 * Real.sin ((M : ℝ) * hh M x) ^ 4) * Real.sin ((M : ℝ) * eta M i (hh M x)) ^ 2 := by
    intro i _
    rw [show cc M i ^ 2 + cc M (i + x) ^ 2 = (cc M i + cc M (i + x)) ^ 2 - 2 * (cc M i * cc M (i + x))
      by ring, cc_diff, cc_sum, cc_prod]
    ring
  rw [Finset.sum_congr rfl e, Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum,
    B2 M hM, B1 M hM]
  ring

theorem inner_c : ∑ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) ^ 2 * (cc M i * cc M (i + x))
    = ((2 * M + 1 : ℕ) : ℝ) * Real.sin ((M : ℝ) * hh M x) ^ 2 / 2
      - 2 * ((2 * M + 1 : ℕ) : ℝ) * Real.sin ((M : ℝ) * hh M x) ^ 4 := by
  have e : ∀ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) ^ 2 * (cc M i * cc M (i + x))
      = (4 * Real.sin ((M : ℝ) * hh M x) ^ 2)
        * (Real.sin ((M : ℝ) * eta M i (hh M x)) ^ 2 * Real.cos ((M : ℝ) * eta M i (hh M x)) ^ 2)
        - (4 * Real.sin ((M : ℝ) * hh M x) ^ 4) * Real.sin ((M : ℝ) * eta M i (hh M x)) ^ 2 := by
    intro i _
    rw [cc_diff, cc_prod]
    ring
  rw [Finset.sum_congr rfl e, Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.mul_sum,
    B2 M hM, B1 M hM]
  ring

theorem inner_d : ∑ i ∈ range (2 * M + 1),
    (cc M i - cc M (i + x)) * (CC M i - CC M (i + x)) * (cc M i + cc M (i + x))
    = -2 * ((2 * M + 1 : ℕ) : ℝ) * Real.sin ((M : ℝ) * hh M x) * Real.cos ((M : ℝ) * hh M x)
      * Real.sin (hh M x) * Real.cos (((2 * M + 1 : ℕ) : ℝ) * hh M x) := by
  have e : ∀ i ∈ range (2 * M + 1),
      (cc M i - cc M (i + x)) * (CC M i - CC M (i + x)) * (cc M i + cc M (i + x))
      = (8 * Real.sin ((M : ℝ) * hh M x) * Real.cos ((M : ℝ) * hh M x) * Real.sin (hh M x))
        * (Real.sin ((M : ℝ) * eta M i (hh M x)) * Real.cos ((M : ℝ) * eta M i (hh M x))
          * Real.sin (eta M i (hh M x))) := by
    intro i _
    rw [cc_diff, CC_diff, cc_sum]
    ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum, B5 M hM]
  ring

theorem inner_e : ∑ i ∈ range (2 * M + 1), (CC M i - CC M (i + x)) ^ 2
    = 2 * ((2 * M + 1 : ℕ) : ℝ) * Real.sin (hh M x) ^ 2 := by
  have e : ∀ i ∈ range (2 * M + 1), (CC M i - CC M (i + x)) ^ 2
      = (4 * Real.sin (hh M x) ^ 2) * Real.sin (eta M i (hh M x)) ^ 2 := by
    intro i _; rw [CC_diff]; ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum, B6 M hM]
  ring

theorem inner_s : ∑ i ∈ range (2 * M + 1), (SS M (i + x) - SS M i) ^ 2
    = 2 * ((2 * M + 1 : ℕ) : ℝ) * Real.sin (hh M x) ^ 2 := by
  have e : ∀ i ∈ range (2 * M + 1), (SS M (i + x) - SS M i) ^ 2
      = (4 * Real.sin (hh M x) ^ 2) * Real.cos (eta M i (hh M x)) ^ 2 := by
    intro i _; rw [SS_diff]; ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum, B7 M hM]
  ring

theorem inner_lin : ∑ i ∈ range (2 * M + 1), (SS M (i + x) - SS M i) = 0 := by
  have e : ∀ i ∈ range (2 * M + 1), (SS M (i + x) - SS M i)
      = (2 * Real.sin (hh M x)) * Real.cos (eta M i (hh M x)) := by
    intro i _; rw [SS_diff]; ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum, B8 M hM, mul_zero]

theorem inner_3a : ∑ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) ^ 2 * (cc M i + cc M (i + x)) = 0 := by
  have e : ∀ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) ^ 2 * (cc M i + cc M (i + x))
      = (8 * Real.sin ((M : ℝ) * hh M x) ^ 2 * Real.cos ((M : ℝ) * hh M x))
        * (Real.sin ((M : ℝ) * eta M i (hh M x)) ^ 2 * Real.cos ((M : ℝ) * eta M i (hh M x))) := by
    intro i _; rw [cc_diff, cc_sum]; ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum, B9 M hM, mul_zero]

theorem inner_3b : ∑ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) * (CC M i - CC M (i + x)) = 0 := by
  have e : ∀ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) * (CC M i - CC M (i + x))
      = (4 * Real.sin ((M : ℝ) * hh M x) * Real.sin (hh M x))
        * (Real.sin ((M : ℝ) * eta M i (hh M x)) * Real.sin (eta M i (hh M x))) := by
    intro i _; rw [cc_diff, CC_diff]; ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum, B10 M hM, mul_zero]

theorem inner_f : ∑ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) ^ 4
    = 6 * ((2 * M + 1 : ℕ) : ℝ) * Real.sin ((M : ℝ) * hh M x) ^ 4 := by
  have e : ∀ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) ^ 4
      = (16 * Real.sin ((M : ℝ) * hh M x) ^ 4) * Real.sin ((M : ℝ) * eta M i (hh M x)) ^ 4 := by
    intro i _; rw [cc_diff]; ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum, B3 M hM]
  ring

theorem inner_g : ∑ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) ^ 2 * (SS M (i + x) - SS M i)
    = -2 * ((2 * M + 1 : ℕ) : ℝ) * Real.sin ((M : ℝ) * hh M x) ^ 2 * Real.sin (hh M x)
      * Real.cos (((2 * M + 1 : ℕ) : ℝ) * hh M x) := by
  have e : ∀ i ∈ range (2 * M + 1), (cc M i - cc M (i + x)) ^ 2 * (SS M (i + x) - SS M i)
      = (8 * Real.sin ((M : ℝ) * hh M x) ^ 2 * Real.sin (hh M x))
        * (Real.sin ((M : ℝ) * eta M i (hh M x)) ^ 2 * Real.cos (eta M i (hh M x))) := by
    intro i _; rw [cc_diff, SS_diff]; ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum, B4 M hM]
  ring

end

end Havelock.Conj13.Odd4
