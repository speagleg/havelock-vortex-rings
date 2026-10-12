/-
Havelock — Conjecture 13: elementary Taylor expansions at `0` with `O(xᵏ)` remainders.

These are the one-variable series used to expand the reduced Hamiltonian along the even-`n`
symmetric curve (paper `conj13/paper/main.tex`, proof of Thm even). All remainders are stated as
`Asymptotics.IsBigO` at `𝓝 0`; the expansions are derived from Mathlib's bounds
`Real.abs_log_sub_add_sum_range_le` (logarithm) and `Real.exp_bound` (exponential).

* `log_one_sub_jet`   `log(1−t) = −t − t²/2 + O(t³)`.
* `log_one_sub_sq_jet` `log(1−(cx)²) = −c²x² − c⁴x⁴/2 + O(x⁶)`.
* `artanh_jet`        `(log(1+x) − log(1−x))/2 = x + x³/3 + O(x⁵)`.
* `log_cosh_jet`      `log(cosh y) = y²/2 − y⁴/12 + O(y⁶)`.
* `log_cosh_comp_jet` the same composed with `y(x) = c₁x + c₃x³ + O(x⁵)`.
-/
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.Asymptotics.Lemmas

open Asymptotics Filter Topology Real

namespace Havelock.Conj13.Series

/-- `xᵏ·q(x) = O(xᵏ)` at `0` when `q` is continuous at `0`. -/
theorem mul_cont {k : ℕ} {q : ℝ → ℝ} (hq : ContinuousAt q 0) :
    (fun x => x ^ k * q x) =O[𝓝 0] (fun x => x ^ k) := by
  have := (isBigO_refl (fun x : ℝ => x ^ k) (𝓝 0)).mul (hq.tendsto.isBigO_one ℝ)
  simpa using this

/-- `xᵏ = O(xʲ)` at `0` for `j ≤ k`. -/
theorem pow_bigO_pow {j k : ℕ} (h : j ≤ k) : (fun x : ℝ => x ^ k) =O[𝓝 0] (fun x => x ^ j) := by
  have e : (fun x : ℝ => x ^ k) = fun x => x ^ j * x ^ (k - j) := by
    funext x; rw [← pow_add]; congr 1; omega
  rw [e]; exact mul_cont (by fun_prop)

/-- `(c·x)ᵏ = O(xᵏ)`. -/
theorem cpow_bigO (c : ℝ) (k : ℕ) : (fun x : ℝ => (c * x) ^ k) =O[𝓝 0] (fun x => x ^ k) := by
  have e : (fun x : ℝ => (c * x) ^ k) = fun x => x ^ k * c ^ k := by funext x; ring
  rw [e]; exact mul_cont continuousAt_const

/-- If `f = O(xʲ)` and `g = O(xᵏ)` then `f·g = O(xʲ⁺ᵏ)`. -/
theorem mul_pow_bigO {f g : ℝ → ℝ} {j k : ℕ} (hf : f =O[𝓝 0] (fun x => x ^ j))
    (hg : g =O[𝓝 0] (fun x => x ^ k)) : (fun x => f x * g x) =O[𝓝 0] (fun x => x ^ (j + k)) := by
  have := hf.mul hg
  simpa [← pow_add] using this

/-- If `f = O(xʲ)` then `fᵐ = O(x^{jm})`. -/
theorem pow_of_bigO {f : ℝ → ℝ} {j : ℕ} (hf : f =O[𝓝 0] (fun x => x ^ j)) (m : ℕ) :
    (fun x => f x ^ m) =O[𝓝 0] (fun x => x ^ (j * m)) := by
  have := hf.pow m
  simpa [← pow_mul] using this

/-- `f = O(xʲ)` and `j ≥ 1` give `f → 0`. -/
theorem tendsto_of_bigO {f : ℝ → ℝ} {j : ℕ} (hj : 1 ≤ j) (hf : f =O[𝓝 0] (fun x => x ^ j)) :
    Tendsto f (𝓝 0) (𝓝 0) := by
  have h0 : Tendsto (fun x : ℝ => x ^ j) (𝓝 0) (𝓝 0) := by
    have : Continuous (fun x : ℝ => x ^ j) := by fun_prop
    simpa [zero_pow (by omega : j ≠ 0)] using this.tendsto 0
  exact hf.trans_tendsto h0

/-- Composition: if `F(y) = O(yᵐ)` and `y(x) = O(xʲ)` with `y → 0`, then `F(y(x)) = O(x^{jm})`. -/
theorem comp_bigO {F y : ℝ → ℝ} {j m : ℕ} (hj : 1 ≤ j) (hF : F =O[𝓝 0] (fun t => t ^ m))
    (hy : y =O[𝓝 0] (fun x => x ^ j)) : (fun x => F (y x)) =O[𝓝 0] (fun x => x ^ (j * m)) :=
  (hF.comp_tendsto (tendsto_of_bigO hj hy)).trans (pow_of_bigO hy m)

/-! ### Logarithm -/

/-- `log(1−t) + Σ_{i<N} t^{i+1}/(i+1) = O(t^{N+1})` (from `Real.abs_log_sub_add_sum_range_le`). -/
theorem log_one_sub_series (N : ℕ) :
    (fun t : ℝ => (∑ i ∈ Finset.range N, t ^ (i + 1) / (i + 1)) + log (1 - t))
      =O[𝓝 0] (fun t => t ^ (N + 1)) := by
  refine IsBigO.of_bound 2 ?_
  filter_upwards [Metric.ball_mem_nhds (0 : ℝ) (by norm_num : (0 : ℝ) < 1 / 2)] with t ht
  rw [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs] at ht
  have h := Real.abs_log_sub_add_sum_range_le (by linarith : |t| < 1) N
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
  refine h.trans ?_
  rw [div_le_iff₀ (by linarith)]
  have := pow_nonneg (abs_nonneg t) (N + 1)
  nlinarith

/-- `log(1−t) = −t − t²/2 + O(t³)`. -/
theorem log_one_sub_jet : (fun t : ℝ => log (1 - t) + t + t ^ 2 / 2) =O[𝓝 0] (fun t => t ^ 3) := by
  refine (log_one_sub_series 2).congr_left (fun t => ?_)
  simp [Finset.sum_range_succ]; ring

/-- `log(1−t) = −t − t²/2 − t³/3 − t⁴/4 + O(t⁵)`. -/
theorem log_one_sub_jet4 :
    (fun t : ℝ => log (1 - t) + t + t ^ 2 / 2 + t ^ 3 / 3 + t ^ 4 / 4) =O[𝓝 0] (fun t => t ^ 5) := by
  refine (log_one_sub_series 4).congr_left (fun t => ?_)
  simp [Finset.sum_range_succ]; ring

/-- `log(1−(cx)²) = −c²x² − c⁴x⁴/2 + O(x⁶)`. -/
theorem log_one_sub_sq_jet (c : ℝ) :
    (fun x : ℝ => log (1 - (c * x) ^ 2) + c ^ 2 * x ^ 2 + c ^ 4 * x ^ 4 / 2)
      =O[𝓝 0] (fun x => x ^ 6) := by
  have hy : (fun x : ℝ => (c * x) ^ 2) =O[𝓝 0] (fun x => x ^ 2) := cpow_bigO c 2
  have := comp_bigO (by norm_num) log_one_sub_jet hy
  refine this.congr_left (fun x => ?_)
  ring

/-- `artanh x = (log(1+x) − log(1−x))/2 = x + x³/3 + O(x⁵)`. -/
theorem artanh_jet :
    (fun x : ℝ => (log (1 + x) - log (1 - x)) / 2 - x - x ^ 3 / 3) =O[𝓝 0] (fun x => x ^ 5) := by
  have h1 := log_one_sub_jet4
  have hy : (fun x : ℝ => -x) =O[𝓝 0] (fun x => x ^ 1) := by
    simpa using (isBigO_refl (fun x : ℝ => x) (𝓝 0)).neg_left
  have h2 := comp_bigO (by norm_num) log_one_sub_jet4 hy
  have := (h2.sub h1).const_mul_left (1 / 2)
  refine this.congr_left (fun x => ?_)
  simp only [sub_neg_eq_add, neg_sq]
  ring_nf

/-- `f = O(xᵏ)` and `q` continuous at `0` give `f·q = O(xᵏ)`. -/
theorem bigO_mul_cont {f q : ℝ → ℝ} {k : ℕ} (hf : f =O[𝓝 0] (fun x => x ^ k))
    (hq : ContinuousAt q 0) : (fun x => f x * q x) =O[𝓝 0] (fun x => x ^ k) := by
  have := hf.mul (hq.tendsto.isBigO_one ℝ)
  simpa using this

/-- `f = O(xʲ)`, `j ≤ k`... the weakening `O(xᵏ) ⊆ O(xʲ)`. -/
theorem bigO_weaken {f : ℝ → ℝ} {j k : ℕ} (h : j ≤ k) (hf : f =O[𝓝 0] (fun x => x ^ k)) :
    f =O[𝓝 0] (fun x => x ^ j) :=
  hf.trans (pow_bigO_pow h)

/-! ### Exponential and hyperbolic functions -/

/-- `exp x = Σ_{m<6} xᵐ/m! + O(x⁶)` (from `Real.exp_bound`). -/
theorem exp_jet :
    (fun x : ℝ => exp x - (1 + x + x ^ 2 / 2 + x ^ 3 / 6 + x ^ 4 / 24 + x ^ 5 / 120))
      =O[𝓝 0] (fun x => x ^ 6) := by
  refine IsBigO.of_bound 1 ?_
  filter_upwards [Metric.ball_mem_nhds (0 : ℝ) (by norm_num : (0 : ℝ) < 1)] with x hx
  rw [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs] at hx
  have h := Real.exp_bound hx.le (n := 6) (by norm_num)
  simp [Finset.sum_range_succ, Nat.factorial] at h
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow]
  refine le_trans (le_of_eq_of_le ?_ h) ?_
  · congr 1
  · have := pow_nonneg (abs_nonneg x) 6
    norm_num
    nlinarith

/-- `cosh y = 1 + y²/2 + y⁴/24 + O(y⁶)`. -/
theorem cosh_jet :
    (fun y : ℝ => cosh y - 1 - y ^ 2 / 2 - y ^ 4 / 24) =O[𝓝 0] (fun y => y ^ 6) := by
  have hy : (fun x : ℝ => -x) =O[𝓝 0] (fun x => x ^ 1) := by
    simpa using (isBigO_refl (fun x : ℝ => x) (𝓝 0)).neg_left
  have h2 := comp_bigO (by norm_num) exp_jet hy
  have := (exp_jet.add h2).const_mul_left (1 / 2)
  refine this.congr' ?_ (Eventually.of_forall fun x => by ring)
  exact Eventually.of_forall fun y => by dsimp only; rw [Real.cosh_eq]; ring

/-- `log(cosh y) = y²/2 − y⁴/12 + O(y⁶)`. -/
theorem log_cosh_jet :
    (fun y : ℝ => log (cosh y) - y ^ 2 / 2 + y ^ 4 / 12) =O[𝓝 0] (fun y => y ^ 6) := by
  set e : ℝ → ℝ := fun y => cosh y - 1 - y ^ 2 / 2 - y ^ 4 / 24 with he_def
  have he : e =O[𝓝 0] (fun y => y ^ 6) := cosh_jet
  set t : ℝ → ℝ := fun y => 1 - cosh y with ht_def
  have ht : t =O[𝓝 0] (fun y => y ^ 2) := by
    have h1 : (fun y : ℝ => y ^ 2 * (-(1 / 2 + y ^ 2 / 24))) =O[𝓝 0] (fun y => y ^ 2) :=
      mul_cont (by fun_prop)
    refine (h1.sub (bigO_weaken (by norm_num) he)).congr_left (fun y => ?_)
    simp only [he_def, ht_def]; ring
  have hL := comp_bigO (by norm_num) log_one_sub_jet ht
  have hq : ContinuousAt (fun y : ℝ => 1 - (y ^ 2 / 2 + y ^ 4 / 24)) 0 := by fun_prop
  have hee : (fun y => e y * e y) =O[𝓝 0] (fun y => y ^ 6) :=
    bigO_weaken (by norm_num) (mul_pow_bigO he he)
  have hy6 : (fun y : ℝ => y ^ 6 * (1 / 48 + y ^ 2 / 1152)) =O[𝓝 0] (fun y => y ^ 6) :=
    mul_cont (by fun_prop)
  have := ((hL.add (bigO_mul_cont he hq)).sub (hee.const_mul_left (1 / 2))).sub hy6
  refine this.congr_left (fun y => ?_)
  have hc : log (1 - t y) = log (cosh y) := by simp [ht_def]
  simp only [hc]
  simp only [he_def, ht_def]
  ring

/-- Composition: if `Y(x) = c₁x + c₃x³ + O(x⁵)` then
`log cosh Y(x) = (c₁²/2) x² + (c₁c₃ − c₁⁴/12) x⁴ + O(x⁶)`. -/
theorem log_cosh_comp_jet {Y : ℝ → ℝ} {c₁ c₃ : ℝ}
    (hY : (fun x => Y x - c₁ * x - c₃ * x ^ 3) =O[𝓝 0] (fun x => x ^ 5)) :
    (fun x => log (cosh (Y x)) - c₁ ^ 2 / 2 * x ^ 2 - (c₁ * c₃ - c₁ ^ 4 / 12) * x ^ 4)
      =O[𝓝 0] (fun x => x ^ 6) := by
  set r : ℝ → ℝ := fun x => Y x - c₁ * x - c₃ * x ^ 3 with hr
  have hY1 : Y =O[𝓝 0] (fun x => x ^ 1) := by
    have h1 : (fun x : ℝ => x ^ 1 * (c₁ + c₃ * x ^ 2)) =O[𝓝 0] (fun x => x ^ 1) :=
      mul_cont (by fun_prop)
    refine (h1.add (bigO_weaken (by norm_num) hY)).congr_left (fun x => ?_)
    ring
  have hA := comp_bigO (by norm_num) log_cosh_jet hY1
  -- `Y² − c₁²x² − 2c₁c₃x⁴ = r·(Y + c₁x + c₃x³) + c₃²x⁶`
  have hp : (fun x : ℝ => Y x + c₁ * x + c₃ * x ^ 3) =O[𝓝 0] (fun x => x ^ 1) := by
    have h1 : (fun x : ℝ => x ^ 1 * (c₁ + c₃ * x ^ 2)) =O[𝓝 0] (fun x => x ^ 1) :=
      mul_cont (by fun_prop)
    refine (hY1.add h1).congr_left (fun x => ?_)
    ring
  have hB : (fun x => r x * (Y x + c₁ * x + c₃ * x ^ 3)) =O[𝓝 0] (fun x => x ^ 6) :=
    mul_pow_bigO hY hp
  have hB6 : (fun x : ℝ => x ^ 6 * c₃ ^ 2) =O[𝓝 0] (fun x => x ^ 6) := mul_cont continuousAt_const
  -- `Y⁴ − c₁⁴x⁴ = (Y − c₁x)(Y + c₁x)(Y² + c₁²x²)`
  have hm : (fun x => Y x - c₁ * x) =O[𝓝 0] (fun x => x ^ 3) := by
    have h1 : (fun x : ℝ => x ^ 3 * c₃) =O[𝓝 0] (fun x => x ^ 3) := mul_cont continuousAt_const
    refine (h1.add (bigO_weaken (by norm_num) hY)).congr_left (fun x => ?_)
    ring
  have hpl : (fun x => Y x + c₁ * x) =O[𝓝 0] (fun x => x ^ 1) := by
    have h1 : (fun x : ℝ => x ^ 1 * c₁) =O[𝓝 0] (fun x => x ^ 1) := mul_cont continuousAt_const
    refine (hY1.add h1).congr_left (fun x => ?_)
    ring
  have hsq : (fun x => Y x ^ 2 + c₁ ^ 2 * x ^ 2) =O[𝓝 0] (fun x => x ^ 2) := by
    have h1 : (fun x : ℝ => x ^ 2 * c₁ ^ 2) =O[𝓝 0] (fun x => x ^ 2) := mul_cont continuousAt_const
    refine ((pow_of_bigO hY1 2).add h1).congr_left (fun x => ?_)
    ring
  have hC : (fun x => (Y x - c₁ * x) * (Y x + c₁ * x) * (Y x ^ 2 + c₁ ^ 2 * x ^ 2))
      =O[𝓝 0] (fun x => x ^ 6) := mul_pow_bigO (mul_pow_bigO hm hpl) hsq
  have := (hA.add (((hB.add hB6)).const_mul_left (1 / 2))).sub (hC.const_mul_left (1 / 12))
  refine this.congr_left (fun x => ?_)
  simp only [hr]
  ring

end Havelock.Conj13.Series
