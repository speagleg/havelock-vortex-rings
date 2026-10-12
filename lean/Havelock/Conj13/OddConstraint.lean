/-
Havelock — Conjecture 13, odd `n = 2M+1`: the `J₁` constraint along the paper's ansatz
(paper `conj13/paper/main.tex`, proof of Thm odd, "Constraint").

The ansatz in canonical coordinates: `w_k = w₀ + A cos(lφ_k) + A² p cos φ_k`,
`ϑ_k = φ_k + A² s sin φ_k` (`oddAnsatz`). Results:

* `J3_oddAnsatz`   `J₃ = Σ w_k` is exactly constant (`J₃` is linear in `w`, `Σ cos(lφ_k) = Σ cos φ_k = 0`).
* `J1_oddAnsatz_jet`  **`J₁ = (n/2) A² J1con(λ, p, s) + O(A³)`**, where
  `J1con λ p s = Q'(w₀)p − Q(w₀)s + Q''(w₀)/4` is the definition used in `Conj13Algebra.lean`
  (with `Q = √P`, `Q(w₀) = 1/(1+λ)`, `Q'(w₀) = (1−λ)/2`, `Q''(w₀) = −(1+λ)³/4`, see `Q_values`), for
  every odd `n ≥ 5`, `λ ∈ (−1, 1)` and all `p, s`.
-/
import Havelock.Conj13.OddQuadratic

open Real Asymptotics Filter Topology Finset

namespace Havelock.Conj13.Odd

open Havelock.Conj13.Series Havelock.Conj13.Config Havelock.Conj13.OddSums

/-! ### Jets of `cos`, `sin` and products -/

theorem cos_sub_one_bigO : (fun t : ℝ => Real.cos t - 1) =O[𝓝 0] (fun t => t ^ 2) := by
  refine IsBigO.of_bound 1 ?_
  filter_upwards [Metric.ball_mem_nhds (0 : ℝ) one_pos] with t ht
  rw [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs] at ht
  have h := Real.cos_bound ht.le
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow, one_mul]
  have h2 : |t| ^ 4 ≤ |t| ^ 2 := pow_le_pow_of_le_one (abs_nonneg t) ht.le (by norm_num)
  have h3 : |t ^ 2 / 2| = |t| ^ 2 / 2 := by rw [abs_div, abs_pow]; norm_num
  calc |Real.cos t - 1| = |(Real.cos t - (1 - t ^ 2 / 2)) - t ^ 2 / 2| := by ring_nf
    _ ≤ |Real.cos t - (1 - t ^ 2 / 2)| + |t ^ 2 / 2| := abs_sub _ _
    _ ≤ |t| ^ 2 := by rw [h3]; nlinarith [pow_nonneg (abs_nonneg t) 4]

theorem sin_sub_self_bigO : (fun t : ℝ => Real.sin t - t) =O[𝓝 0] (fun t => t ^ 3) := by
  refine IsBigO.of_bound 1 ?_
  filter_upwards [Metric.ball_mem_nhds (0 : ℝ) one_pos] with t ht
  rw [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs] at ht
  have h := Real.sin_bound ht.le
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow, one_mul]
  have h2 : |t| ^ 5 ≤ |t| ^ 3 := pow_le_pow_of_le_one (abs_nonneg t) ht.le (by norm_num)
  have h3 : |t ^ 3 / 6| = |t| ^ 3 / 6 := by rw [abs_div, abs_pow]; norm_num
  calc |Real.sin t - t| = |(Real.sin t - (t - t ^ 3 / 6)) - t ^ 3 / 6| := by ring_nf
    _ ≤ |Real.sin t - (t - t ^ 3 / 6)| + |t ^ 3 / 6| := abs_sub _ _
    _ ≤ |t| ^ 3 := by rw [h3]; nlinarith [pow_nonneg (abs_nonneg t) 5]

/-- Product of two second-order jets. -/
theorem mul_jet2 {f g : ℝ → ℝ} {a0 a1 a2 c0 c1 c2 : ℝ}
    (hf : (fun A => f A - (a0 + a1 * A + a2 * A ^ 2)) =O[𝓝 0] (fun A => A ^ 3))
    (hg : (fun A => g A - (c0 + c1 * A + c2 * A ^ 2)) =O[𝓝 0] (fun A => A ^ 3)) :
    (fun A => f A * g A - (a0 * c0 + (a0 * c1 + a1 * c0) * A + (a0 * c2 + a1 * c1 + a2 * c0) * A ^ 2))
      =O[𝓝 0] (fun A => A ^ 3) := by
  have h1 : (fun A : ℝ => A ^ 3 * ((a1 * c2 + a2 * c1) + a2 * c2 * A)) =O[𝓝 0] (fun A => A ^ 3) :=
    mul_cont (by fun_prop)
  have h2 := bigO_mul_cont hg (q := fun A => a0 + a1 * A + a2 * A ^ 2) (by fun_prop)
  have h3 := bigO_mul_cont hf (q := fun A => c0 + c1 * A + c2 * A ^ 2) (by fun_prop)
  have h4 := bigO_weaken (by norm_num : 3 ≤ 3 + 3) (mul_pow_bigO hf hg)
  refine (((h1.add h2).add h3).add h4).congr_left (fun A => ?_)
  ring

/-! ### The ansatz -/

/-- `sin(2πx/n)`. -/
noncomputable def SS (M : ℕ) (x : ℝ) : ℝ := Real.sin (2 * π * x / (2 * M + 1 : ℕ))

/-- The paper's odd ansatz in canonical coordinates:
`w_k = w₀ + A cos(lφ_k) + A² p cos φ_k`, `ϑ_k = φ_k + A² s sin φ_k`. -/
noncomputable def oddAnsatz (M : ℕ) (b p s A : ℝ) (k : Fin (2 * M + 1)) : ℂ :=
  (√(uOf b (w0 b + A * cc M ((k : ℕ) : ℝ) + A ^ 2 * p * CC M ((k : ℕ) : ℝ))) : ℂ)
    * Complex.exp (((2 * π * ((k : ℕ) : ℝ) / (2 * M + 1 : ℕ) + A ^ 2 * s * SS M ((k : ℕ) : ℝ) : ℝ) : ℂ)
      * Complex.I)

/-- In the chart (`|h| < w₀`): `0 < w₀ + h` and `0 < 1 − λ(w₀ + h)`. -/
theorem w_pos' (b h : ℝ) (hb : -1 < b) (hb1 : b < 1) (hh : |h| < w0 b) :
    0 < w0 b + h ∧ 0 < 1 - b * (w0 b + h) := by
  have h1 : 0 < 1 + b := by linarith
  have hb' : |b| ≤ 1 := abs_le.mpr ⟨by linarith, by linarith⟩
  have hbh : |b * h| ≤ |h| := by rw [abs_mul]; nlinarith [abs_nonneg b, abs_nonneg h]
  have e : 1 - b * (w0 b + h) = w0 b - b * h := by
    unfold w0; rw [mul_add, show b * (1 / (1 + b)) = 1 - 1 / (1 + b) by field_simp; ring]; ring
  rw [e]
  constructor
  · have := neg_abs_le h; linarith
  · have := le_abs_self (b * h); linarith

/-- `J₃ = Σ_k w_k` is exactly `J₃(ring)` along the ansatz (`|A| + A²|p| < w₀`). -/
theorem J3_oddAnsatz (M : ℕ) (hM : 1 ≤ M) (b p s A : ℝ) (hb : -1 < b) (hb1 : b < 1)
    (hA : |A| + A ^ 2 * |p| < w0 b) :
    EvenQuartic.J3 b (oddAnsatz M b p s A) = (2 * M + 1 : ℕ) * w0 b := by
  have hterm : ∀ k : Fin (2 * M + 1),
      Complex.normSq (oddAnsatz M b p s A k) / sig b (oddAnsatz M b p s A k)
        = w0 b + A * cc M ((k : ℕ) : ℝ) + A ^ 2 * p * CC M ((k : ℕ) : ℝ) := by
    intro k
    have hh : |A * cc M ((k : ℕ) : ℝ) + A ^ 2 * p * CC M ((k : ℕ) : ℝ)| < w0 b := by
      have e1 : |A * cc M ((k : ℕ) : ℝ)| ≤ |A| := by
        rw [abs_mul]; exact mul_le_of_le_one_right (abs_nonneg _) (cc_abs_le M _)
      have e2 : |A ^ 2 * p * CC M ((k : ℕ) : ℝ)| ≤ A ^ 2 * |p| := by
        rw [abs_mul, abs_mul, abs_of_nonneg (sq_nonneg A)]
        exact mul_le_of_le_one_right (by positivity) (Real.abs_cos_le_one _)
      calc |A * cc M ((k : ℕ) : ℝ) + A ^ 2 * p * CC M ((k : ℕ) : ℝ)|
          ≤ |A * cc M ((k : ℕ) : ℝ)| + |A ^ 2 * p * CC M ((k : ℕ) : ℝ)| := abs_add_le _ _
        _ ≤ |A| + A ^ 2 * |p| := add_le_add e1 e2
        _ < w0 b := hA
    obtain ⟨hw, h1⟩ := w_pos' b _ hb hb1 hh
    rw [← add_assoc] at hw h1
    have hu : 0 ≤ uOf b (w0 b + A * cc M ((k : ℕ) : ℝ) + A ^ 2 * p * CC M ((k : ℕ) : ℝ)) :=
      (div_pos hw h1).le
    obtain ⟨E1, E2⟩ := uOf_rel b _ h1
    have hs : 1 + b * uOf b (w0 b + A * cc M ((k : ℕ) : ℝ) + A ^ 2 * p * CC M ((k : ℕ) : ℝ)) ≠ 0 :=
      fun h => by rw [h, zero_mul] at E1; exact zero_ne_one E1
    unfold oddAnsatz sig
    rw [normSq_polar, Real.sq_sqrt hu, div_eq_iff hs]
    linear_combination (-(uOf b (w0 b + A * cc M ((k : ℕ) : ℝ) + A ^ 2 * p * CC M ((k : ℕ) : ℝ)))) * E1
      + (1 + b * uOf b (w0 b + A * cc M ((k : ℕ) : ℝ) + A ^ 2 * p * CC M ((k : ℕ) : ℝ))) * E2
  have hC : ∑ k ∈ range (2 * M + 1), CC M (k : ℝ) = 0 := by
    have := sum_cos_eq_zero (2 * M + 1) 1 (by omega)
      (fun h => by have := Nat.le_of_dvd one_pos h; omega)
    unfold CC; simpa using this
  unfold EvenQuartic.J3
  rw [Finset.sum_congr rfl (fun k _ => hterm k), Finset.sum_add_distrib, Finset.sum_add_distrib,
    ← Finset.mul_sum, ← Finset.mul_sum,
    Fin.sum_univ_eq_sum_range (fun k => cc M (k : ℝ)) (2 * M + 1),
    Fin.sum_univ_eq_sum_range (fun k => CC M (k : ℝ)) (2 * M + 1), sum_cc M hM, hC]
  simp

/-! ### Fourier sums for the constraint -/

theorem sum_CC (M : ℕ) (hM : 1 ≤ M) : ∑ k ∈ range (2 * M + 1), CC M (k : ℝ) = 0 := by
  have := sum_cos_eq_zero (2 * M + 1) 1 (by omega) (fun h => by have := Nat.le_of_dvd one_pos h; omega)
  unfold CC; simpa using this

theorem sum_CC_sq (M : ℕ) (hM : 1 ≤ M) :
    ∑ k ∈ range (2 * M + 1), CC M (k : ℝ) ^ 2 = ((2 * M + 1 : ℕ) : ℝ) / 2 := by
  have h := sum_cos_eq_zero (2 * M + 1) 2 (by omega) (fun h => by have := Nat.le_of_dvd (by norm_num) h; omega)
  have e : ∀ k ∈ range (2 * M + 1), CC M (k : ℝ) ^ 2
      = 1 / 2 + Real.cos (2 * π * ((2 : ℕ) : ℝ) * k / ((2 * M + 1 : ℕ) : ℝ)) / 2 := by
    intro k _; unfold CC; rw [Real.cos_sq]; congr 2; push_cast; ring_nf
  rw [Finset.sum_congr rfl e, Finset.sum_add_distrib, ← Finset.sum_div, ← Finset.sum_div, h]
  simp

theorem sum_SS_sq (M : ℕ) (hM : 1 ≤ M) :
    ∑ k ∈ range (2 * M + 1), SS M (k : ℝ) ^ 2 = ((2 * M + 1 : ℕ) : ℝ) / 2 := by
  have e : ∀ k ∈ range (2 * M + 1), SS M (k : ℝ) ^ 2 = 1 - CC M (k : ℝ) ^ 2 := by
    intro k _; unfold SS CC; rw [Real.sin_sq]
  rw [Finset.sum_congr rfl e, Finset.sum_sub_distrib, sum_CC_sq M hM]
  simp; ring

theorem sum_cc_CC (M : ℕ) (hM : 2 ≤ M) :
    ∑ k ∈ range (2 * M + 1), cc M (k : ℝ) * CC M (k : ℝ) = 0 := by
  have h1 := sum_cos_eq_zero (2 * M + 1) (M + 1) (by omega)
    (fun h => by have := Nat.le_of_dvd (by omega) h; omega)
  have h2 := sum_cos_eq_zero (2 * M + 1) (M - 1) (by omega)
    (fun h => by have := Nat.le_of_dvd (by omega) h; omega)
  have e : ∀ k ∈ range (2 * M + 1), cc M (k : ℝ) * CC M (k : ℝ)
      = Real.cos (2 * π * ((M + 1 : ℕ) : ℝ) * k / ((2 * M + 1 : ℕ) : ℝ)) / 2
        + Real.cos (2 * π * ((M - 1 : ℕ) : ℝ) * k / ((2 * M + 1 : ℕ) : ℝ)) / 2 := by
    intro k _
    unfold cc CC
    rw [Nat.cast_sub (by omega : 1 ≤ M)]
    set x := 2 * π * (M : ℝ) * k / ((2 * M + 1 : ℕ) : ℝ)
    set y := 2 * π * (k : ℝ) / ((2 * M + 1 : ℕ) : ℝ)
    have ex : 2 * π * ((M + 1 : ℕ) : ℝ) * k / ((2 * M + 1 : ℕ) : ℝ) = x + y := by
      simp only [x, y]; push_cast; ring
    have ey : 2 * π * ((M : ℝ) - ((1 : ℕ) : ℝ)) * k / ((2 * M + 1 : ℕ) : ℝ) = x - y := by
      simp only [x, y]; push_cast; ring
    rw [ex, ey, Real.cos_add, Real.cos_sub]; ring
  rw [Finset.sum_congr rfl e, Finset.sum_add_distrib, ← Finset.sum_div, ← Finset.sum_div, h1, h2]
  simp

theorem sum_cc_sq_CC (M : ℕ) (hM : 1 ≤ M) :
    ∑ k ∈ range (2 * M + 1), cc M (k : ℝ) ^ 2 * CC M (k : ℝ) = ((2 * M + 1 : ℕ) : ℝ) / 4 := by
  have h := sum_cos_sq_mul_cos (2 * M + 1) (by omega) (by omega)
  rw [show (2 * M + 1 - 1) / 2 = M by omega] at h
  unfold cc CC
  exact h

/-! ### `J₁` along the ansatz -/

/-- `Re(z/σ) = √P(w) cos ϑ` in canonical coordinates. -/
theorem re_div_sig (b w θ : ℝ) (hw : 0 < w) (h1 : 0 < 1 - b * w) :
    (((√(uOf b w) : ℂ) * Complex.exp ((θ : ℂ) * Complex.I))
        / (sig b ((√(uOf b w) : ℂ) * Complex.exp ((θ : ℂ) * Complex.I)) : ℂ)).re
      = √(Pw b w) * Real.cos θ := by
  have hu : 0 ≤ uOf b w := (div_pos hw h1).le
  have hr : √(uOf b w) * (1 - b * w) = √(Pw b w) := by
    rw [← Real.sqrt_sq h1.le, ← Real.sqrt_mul hu]
    congr 1; unfold uOf Pw; field_simp
  obtain ⟨E1, -⟩ := uOf_rel b w h1
  have hsig : sig b ((√(uOf b w) : ℂ) * Complex.exp ((θ : ℂ) * Complex.I)) = 1 / (1 - b * w) := by
    unfold sig; rw [normSq_polar, Real.sq_sqrt hu, eq_div_iff h1.ne', E1]
  rw [hsig, show (((1 / (1 - b * w) : ℝ)) : ℂ) = ((1 / (1 - b * w) : ℝ) : ℂ) from rfl,
    div_eq_mul_inv, ← Complex.ofReal_inv, one_div, inv_inv, mul_comm, ← mul_assoc,
    ← Complex.ofReal_mul, mul_comm (1 - b * w), hr, Complex.re_ofReal_mul,
    Complex.exp_ofReal_mul_I_re]

/-- One vortex: `√P(w_k) cos ϑ_k` to second order in `A`. -/
theorem J1_term_jet (b p s c C S φ : ℝ) (hb : -1 < b) (hC : Real.cos φ = C) (hS : Real.sin φ = S) :
    (fun A => √(Pw b (w0 b + A * c + A ^ 2 * p * C)) * Real.cos (φ + A ^ 2 * s * S)
      - (w0 b * C + ((1 - b) / (1 + b) * c / (2 * w0 b)) * C * A
        + (w0 b * (-(s * S ^ 2)) + (((1 - b) / (1 + b) * p * C - b * c ^ 2)
            - ((1 - b) / (1 + b) * c / (2 * w0 b)) ^ 2) / (2 * w0 b) * C) * A ^ 2))
      =O[𝓝 0] (fun A => A ^ 3) := by
  have h1 : 0 < 1 + b := by linarith
  set P0 := (1 / (1 + b)) ^ 2 with hP0
  set P1 := (1 - b) / (1 + b) with hP1
  have hPw : ∀ h : ℝ, Pw b (w0 b + h) = P0 + P1 * h - b * h ^ 2 := by
    intro h; unfold Pw w0; rw [hP0, hP1]; field_simp; ring
  have hg : (fun A => Pw b (w0 b + A * c + A ^ 2 * p * C)
      - (P0 + P1 * c * A + (P1 * p * C - b * c ^ 2) * A ^ 2)) =O[𝓝 0] (fun A => A ^ 3) := by
    have : (fun A : ℝ => A ^ 3 * (-(b * (2 * c * p * C + A * p ^ 2 * C ^ 2))))
        =O[𝓝 0] (fun A => A ^ 3) := mul_cont (by fun_prop)
    refine this.congr_left (fun A => ?_)
    rw [add_assoc, hPw]; ring
  have hsq := sqrt_jet2 (by positivity : 0 < P0) hg
  have hs0 : √P0 = w0 b := by rw [hP0, Real.sqrt_sq (by positivity)]; unfold w0; ring
  rw [hs0] at hsq
  -- `cos(φ + A²sS) = C − sS² A² + O(A³)`
  have hcos : (fun A => Real.cos (φ + A ^ 2 * s * S) - (C + 0 * A + (-(s * S ^ 2)) * A ^ 2))
      =O[𝓝 0] (fun A => A ^ 3) := by
    have ht : (fun A : ℝ => A ^ 2 * s * S) =O[𝓝 0] (fun A => A ^ 2) := by
      have : (fun A : ℝ => A ^ 2 * (s * S)) =O[𝓝 0] (fun A => A ^ 2) := mul_cont continuousAt_const
      exact this.congr_left (fun A => by ring)
    have hc1 := bigO_weaken (by norm_num : 3 ≤ 2 * 2) (comp_bigO (by norm_num) cos_sub_one_bigO ht)
    have hs1 := bigO_weaken (by norm_num : 3 ≤ 2 * 3) (comp_bigO (by norm_num) sin_sub_self_bigO ht)
    refine ((hc1.const_mul_left C).sub (hs1.const_mul_left S)).congr_left (fun A => ?_)
    rw [Real.cos_add, hC, hS]; ring
  have := mul_jet2 hsq hcos
  refine this.congr_left (fun A => ?_)
  ring

/-- **The `J₁` constraint** (paper, proof of Thm odd): along the ansatz,
`J₁ = Re Σ z_k/σ_k = (n/2) A² (Q'(w₀)p − Q(w₀)s + Q''(w₀)/4) + O(A³)`, i.e. `(n/2) A² · J1con λ p s`,
for every odd `n = 2M+1 ≥ 5`, `λ ∈ (−1,1)`, `p`, `s`. -/
theorem J1_oddAnsatz_jet (M : ℕ) (hM : 2 ≤ M) (b p s : ℝ) (hb : -1 < b) (hb1 : b < 1) :
    (fun A => (EvenQuartic.Jp b (oddAnsatz M b p s A)).re
        - ((2 * M + 1 : ℕ) : ℝ) / 2 * J1con b p s * A ^ 2) =O[𝓝 0] (fun A => A ^ 3) := by
  have h1 : 0 < 1 + b := by linarith
  have hw0 : 0 < w0 b := by unfold w0; exact div_pos one_pos h1
  set n := 2 * M + 1 with hn
  -- the per-vortex jets
  have hk : ∀ k ∈ (Finset.univ : Finset (Fin n)),
      (fun A => √(Pw b (w0 b + A * cc M ((k : ℕ) : ℝ) + A ^ 2 * p * CC M ((k : ℕ) : ℝ)))
          * Real.cos (2 * π * ((k : ℕ) : ℝ) / (2 * M + 1 : ℕ) + A ^ 2 * s * SS M ((k : ℕ) : ℝ))
        - (w0 b * CC M ((k : ℕ) : ℝ)
          + ((1 - b) / (1 + b) * cc M ((k : ℕ) : ℝ) / (2 * w0 b)) * CC M ((k : ℕ) : ℝ) * A
          + (w0 b * (-(s * SS M ((k : ℕ) : ℝ) ^ 2)) + (((1 - b) / (1 + b) * p * CC M ((k : ℕ) : ℝ)
              - b * cc M ((k : ℕ) : ℝ) ^ 2) - ((1 - b) / (1 + b) * cc M ((k : ℕ) : ℝ) / (2 * w0 b)) ^ 2)
              / (2 * w0 b) * CC M ((k : ℕ) : ℝ)) * A ^ 2))
        =O[𝓝 0] (fun A => A ^ 3) :=
    fun k _ => J1_term_jet b p s _ _ _ _ hb rfl rfl
  have hsum := IsBigO.fun_sum hk
  -- the sum of the polynomial parts
  have hpoly : ∀ A : ℝ, (∑ k : Fin n, (w0 b * CC M ((k : ℕ) : ℝ)
          + ((1 - b) / (1 + b) * cc M ((k : ℕ) : ℝ) / (2 * w0 b)) * CC M ((k : ℕ) : ℝ) * A
          + (w0 b * (-(s * SS M ((k : ℕ) : ℝ) ^ 2)) + (((1 - b) / (1 + b) * p * CC M ((k : ℕ) : ℝ)
              - b * cc M ((k : ℕ) : ℝ) ^ 2) - ((1 - b) / (1 + b) * cc M ((k : ℕ) : ℝ) / (2 * w0 b)) ^ 2)
              / (2 * w0 b) * CC M ((k : ℕ) : ℝ)) * A ^ 2))
        = ((2 * M + 1 : ℕ) : ℝ) / 2 * J1con b p s * A ^ 2 := by
    intro A
    have e : ∀ k : Fin n, (w0 b * CC M ((k : ℕ) : ℝ)
          + ((1 - b) / (1 + b) * cc M ((k : ℕ) : ℝ) / (2 * w0 b)) * CC M ((k : ℕ) : ℝ) * A
          + (w0 b * (-(s * SS M ((k : ℕ) : ℝ) ^ 2)) + (((1 - b) / (1 + b) * p * CC M ((k : ℕ) : ℝ)
              - b * cc M ((k : ℕ) : ℝ) ^ 2) - ((1 - b) / (1 + b) * cc M ((k : ℕ) : ℝ) / (2 * w0 b)) ^ 2)
              / (2 * w0 b) * CC M ((k : ℕ) : ℝ)) * A ^ 2)
        = w0 b * CC M ((k : ℕ) : ℝ)
          + ((1 - b) / (1 + b) / (2 * w0 b) * A) * (cc M ((k : ℕ) : ℝ) * CC M ((k : ℕ) : ℝ))
          + (-(w0 b * s) * A ^ 2) * SS M ((k : ℕ) : ℝ) ^ 2
          + ((1 - b) / (1 + b) * p / (2 * w0 b) * A ^ 2) * CC M ((k : ℕ) : ℝ) ^ 2
          + (-(b + ((1 - b) / (1 + b) / (2 * w0 b)) ^ 2) / (2 * w0 b) * A ^ 2)
              * (cc M ((k : ℕ) : ℝ) ^ 2 * CC M ((k : ℕ) : ℝ)) := by
      intro k; ring
    rw [Finset.sum_congr rfl (fun k _ => e k)]
    simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
    rw [Fin.sum_univ_eq_sum_range (fun k => CC M (k : ℝ)) n,
      Fin.sum_univ_eq_sum_range (fun k => cc M (k : ℝ) * CC M (k : ℝ)) n,
      Fin.sum_univ_eq_sum_range (fun k => SS M (k : ℝ) ^ 2) n,
      Fin.sum_univ_eq_sum_range (fun k => CC M (k : ℝ) ^ 2) n,
      Fin.sum_univ_eq_sum_range (fun k => cc M (k : ℝ) ^ 2 * CC M (k : ℝ)) n,
      sum_CC M (by omega), sum_cc_CC M hM, sum_SS_sq M (by omega), sum_CC_sq M (by omega),
      sum_cc_sq_CC M (by omega)]
    unfold J1con w0
    field_simp
    ring
  -- `Re(Jp) = Σ √P_k cos ϑ_k` near `A = 0`
  have hev : ∀ᶠ A in 𝓝 (0 : ℝ), |A| + A ^ 2 * |p| < w0 b := by
    have hc : Continuous (fun A : ℝ => |A| + A ^ 2 * |p|) := by fun_prop
    exact hc.continuousAt.eventually_lt continuousAt_const (by simpa using hw0)
  refine hsum.congr' ?_ (Eventually.of_forall fun _ => rfl)
  filter_upwards [hev] with A hA
  rw [← hpoly A, Finset.sum_sub_distrib]
  congr 1
  unfold EvenQuartic.Jp
  rw [Complex.re_sum]
  refine Finset.sum_congr rfl (fun k _ => ?_)
  have hh : |A * cc M ((k : ℕ) : ℝ) + A ^ 2 * p * CC M ((k : ℕ) : ℝ)| < w0 b := by
    have e1 : |A * cc M ((k : ℕ) : ℝ)| ≤ |A| := by
      rw [abs_mul]; exact mul_le_of_le_one_right (abs_nonneg _) (cc_abs_le M _)
    have e2 : |A ^ 2 * p * CC M ((k : ℕ) : ℝ)| ≤ A ^ 2 * |p| := by
      rw [abs_mul, abs_mul, abs_of_nonneg (sq_nonneg A)]
      exact mul_le_of_le_one_right (by positivity) (Real.abs_cos_le_one _)
    calc |A * cc M ((k : ℕ) : ℝ) + A ^ 2 * p * CC M ((k : ℕ) : ℝ)|
        ≤ |A * cc M ((k : ℕ) : ℝ)| + |A ^ 2 * p * CC M ((k : ℕ) : ℝ)| := abs_add_le _ _
      _ ≤ |A| + A ^ 2 * |p| := add_le_add e1 e2
      _ < w0 b := hA
  obtain ⟨hw, hc1⟩ := w_pos' b _ hb hb1 hh
  rw [← add_assoc] at hw hc1
  unfold oddAnsatz
  rw [re_div_sig b _ _ hw hc1]

end Havelock.Conj13.Odd
