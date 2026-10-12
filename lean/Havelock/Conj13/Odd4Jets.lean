/-
Havelock — Conjecture 13: a small calculus of fourth-order jets at `0`.

`Jet4 f a₀ a₁ a₂ a₃ a₄` means `f(A) = a₀ + a₁A + a₂A² + a₃A³ + a₄A⁴ + O(A⁵)` as `A → 0`. Closure
under `+`, `−`, scalar multiples, products, finite sums, `log` (`a₀ > 0`), `√` (`a₀ > 0`) and
inverses (`a₀ ≠ 0`), and the angular jet `1 − cos(θ + A²d)`. Used for the `A⁴` coefficient of the
Hamiltonian along the odd ansatz (paper `conj13/paper/main.tex`, proof of Thm odd).
-/
import Havelock.Conj13.Series

open Real Asymptotics Filter Topology Finset

namespace Havelock.Conj13.Odd4

open Havelock.Conj13.Series

/-- `f(A) = a₀ + a₁A + a₂A² + a₃A³ + a₄A⁴ + O(A⁵)` at `A → 0`. -/
def Jet4 (f : ℝ → ℝ) (a0 a1 a2 a3 a4 : ℝ) : Prop :=
  (fun A => f A - (a0 + a1 * A + a2 * A ^ 2 + a3 * A ^ 3 + a4 * A ^ 4)) =O[𝓝 0] (fun A => A ^ 5)

namespace Jet4

variable {f g : ℝ → ℝ} {a0 a1 a2 a3 a4 c0 c1 c2 c3 c4 : ℝ}

theorem of_coeff (h : Jet4 f a0 a1 a2 a3 a4) {b0 b1 b2 b3 b4 : ℝ} (e0 : a0 = b0) (e1 : a1 = b1)
    (e2 : a2 = b2) (e3 : a3 = b3) (e4 : a4 = b4) : Jet4 f b0 b1 b2 b3 b4 := by
  subst e0 e1 e2 e3 e4; exact h

theorem congr' (h : Jet4 f a0 a1 a2 a3 a4) (hfg : ∀ᶠ A in 𝓝 0, f A = g A) : Jet4 g a0 a1 a2 a3 a4 := by
  unfold Jet4 at h ⊢
  refine h.congr' ?_ (Eventually.of_forall fun _ => rfl)
  filter_upwards [hfg] with A hA
  rw [hA]

theorem congr (h : Jet4 f a0 a1 a2 a3 a4) (hfg : ∀ A, f A = g A) : Jet4 g a0 a1 a2 a3 a4 :=
  h.congr' (Eventually.of_forall hfg)

/-- A polynomial of degree `≤ 4` is its own jet. -/
theorem poly (a0 a1 a2 a3 a4 : ℝ) :
    Jet4 (fun A => a0 + a1 * A + a2 * A ^ 2 + a3 * A ^ 3 + a4 * A ^ 4) a0 a1 a2 a3 a4 := by
  unfold Jet4
  simp only [sub_self]
  exact isBigO_zero _ _

theorem const (c : ℝ) : Jet4 (fun _ => c) c 0 0 0 0 :=
  (poly c 0 0 0 0).congr (fun A => by ring)

theorem add (hf : Jet4 f a0 a1 a2 a3 a4) (hg : Jet4 g c0 c1 c2 c3 c4) :
    Jet4 (fun A => f A + g A) (a0 + c0) (a1 + c1) (a2 + c2) (a3 + c3) (a4 + c4) := by
  unfold Jet4 at *
  exact (hf.add hg).congr_left (fun A => by ring)

theorem sub (hf : Jet4 f a0 a1 a2 a3 a4) (hg : Jet4 g c0 c1 c2 c3 c4) :
    Jet4 (fun A => f A - g A) (a0 - c0) (a1 - c1) (a2 - c2) (a3 - c3) (a4 - c4) := by
  unfold Jet4 at *
  exact (hf.sub hg).congr_left (fun A => by ring)

theorem const_mul (c : ℝ) (hf : Jet4 f a0 a1 a2 a3 a4) :
    Jet4 (fun A => c * f A) (c * a0) (c * a1) (c * a2) (c * a3) (c * a4) := by
  unfold Jet4 at *
  exact (hf.const_mul_left c).congr_left (fun A => by ring)

/-- Adding an `O(A⁵)` function does not change the jet. -/
theorem add_bigO (hf : Jet4 f a0 a1 a2 a3 a4) (hg : g =O[𝓝 0] (fun A => A ^ 5)) :
    Jet4 (fun A => f A + g A) a0 a1 a2 a3 a4 := by
  unfold Jet4 at *
  exact (hf.add hg).congr_left (fun A => by ring)

theorem of_bigO (hg : (fun A => g A - f A) =O[𝓝 0] (fun A => A ^ 5)) (hf : Jet4 f a0 a1 a2 a3 a4) :
    Jet4 g a0 a1 a2 a3 a4 :=
  (hf.add_bigO hg).congr (fun A => by ring)

theorem zero_bigO (hf : Jet4 f 0 0 0 0 0) : f =O[𝓝 0] (fun A => A ^ 5) := by
  unfold Jet4 at hf
  exact hf.congr_left (fun A => by ring)

/-- The remainder of a polynomial truncation. -/
theorem mul (hf : Jet4 f a0 a1 a2 a3 a4) (hg : Jet4 g c0 c1 c2 c3 c4) :
    Jet4 (fun A => f A * g A) (a0 * c0) (a0 * c1 + a1 * c0) (a0 * c2 + a1 * c1 + a2 * c0)
      (a0 * c3 + a1 * c2 + a2 * c1 + a3 * c0) (a0 * c4 + a1 * c3 + a2 * c2 + a3 * c1 + a4 * c0) := by
  unfold Jet4 at *
  have h1 : (fun A : ℝ => A ^ 5 * ((a1 * c4 + a2 * c3 + a3 * c2 + a4 * c1)
      + (a2 * c4 + a3 * c3 + a4 * c2) * A + (a3 * c4 + a4 * c3) * A ^ 2 + a4 * c4 * A ^ 3))
      =O[𝓝 0] (fun A => A ^ 5) := mul_cont (by fun_prop)
  have h2 := bigO_mul_cont hg (q := fun A => a0 + a1 * A + a2 * A ^ 2 + a3 * A ^ 3 + a4 * A ^ 4)
    (by fun_prop)
  have h3 := bigO_mul_cont hf (q := fun A => c0 + c1 * A + c2 * A ^ 2 + c3 * A ^ 3 + c4 * A ^ 4)
    (by fun_prop)
  have h4 := bigO_weaken (by norm_num : 5 ≤ 5 + 5) (mul_pow_bigO hf hg)
  refine (((h1.add h2).add h3).add h4).congr_left (fun A => ?_)
  ring

theorem sum {ι : Type*} (s : Finset ι) {F : ι → ℝ → ℝ} {b0 b1 b2 b3 b4 : ι → ℝ}
    (h : ∀ i ∈ s, Jet4 (F i) (b0 i) (b1 i) (b2 i) (b3 i) (b4 i)) :
    Jet4 (fun A => ∑ i ∈ s, F i A) (∑ i ∈ s, b0 i) (∑ i ∈ s, b1 i) (∑ i ∈ s, b2 i) (∑ i ∈ s, b3 i)
      (∑ i ∈ s, b4 i) := by
  unfold Jet4 at *
  refine (IsBigO.fun_sum h).congr_left (fun A => ?_)
  simp only [Finset.sum_sub_distrib, Finset.sum_add_distrib, Finset.sum_mul]

/-- A jet with vanishing constant term is `O(A)`. -/
theorem bigO_one (hf : Jet4 f 0 a1 a2 a3 a4) : f =O[𝓝 0] (fun A => A ^ 1) := by
  unfold Jet4 at hf
  have h1 : (fun A : ℝ => A ^ 1 * (a1 + a2 * A + a3 * A ^ 2 + a4 * A ^ 3)) =O[𝓝 0] (fun A => A ^ 1) :=
    mul_cont (by fun_prop)
  refine (h1.add (bigO_weaken (by norm_num) hf)).congr_left (fun A => ?_)
  ring

theorem tendsto (hf : Jet4 f a0 a1 a2 a3 a4) : Tendsto f (𝓝 0) (𝓝 a0) := by
  have h0 : Jet4 (fun A => f A - a0) 0 a1 a2 a3 a4 :=
    (hf.sub (const a0)).of_coeff (by ring) (by ring) (by ring) (by ring) (by ring)
  have := tendsto_of_bigO (le_refl 1) h0.bigO_one
  have := this.add_const a0
  simpa using this

theorem eventually_pos (hf : Jet4 f a0 a1 a2 a3 a4) (h0 : 0 < a0) : ∀ᶠ A in 𝓝 0, 0 < f A :=
  hf.tendsto.eventually (lt_mem_nhds h0)

/-- `log` of a jet with positive constant term. -/
theorem log (hf : Jet4 f a0 a1 a2 a3 a4) (h0 : 0 < a0) :
    Jet4 (fun A => Real.log (f A)) (Real.log a0) (a1 / a0) (a2 / a0 - a1 ^ 2 / (2 * a0 ^ 2))
      (a3 / a0 - a1 * a2 / a0 ^ 2 + a1 ^ 3 / (3 * a0 ^ 3))
      (a4 / a0 - (a2 ^ 2 + 2 * a1 * a3) / (2 * a0 ^ 2) + a1 ^ 2 * a2 / a0 ^ 3
        - a1 ^ 4 / (4 * a0 ^ 4)) := by
  set t : ℝ → ℝ := fun A => 1 - f A / a0 with ht_def
  have ht : Jet4 t 0 (-(a1 / a0)) (-(a2 / a0)) (-(a3 / a0)) (-(a4 / a0)) := by
    have := (const 1).sub (hf.const_mul (1 / a0))
    refine (this.congr (fun A => by simp only [ht_def]; ring)).of_coeff ?_ ?_ ?_ ?_ ?_ <;>
      field_simp <;> ring
  have hL := comp_bigO (le_refl 1) (log_one_sub_series 4) ht.bigO_one
  have h2 := ht.mul ht
  have h3 := h2.mul ht
  have h4 := h3.mul ht
  have hS := (((ht.add (h2.const_mul (1 / 2))).add (h3.const_mul (1 / 3))).add
    (h4.const_mul (1 / 4)))
  have hlog := ((const (Real.log a0)).sub hS).add_bigO hL
  refine (hlog.congr' ?_).of_coeff ?_ ?_ ?_ ?_ ?_
  · filter_upwards [hf.eventually_pos h0] with A hA
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, ht_def]
    have : Real.log (f A) = Real.log a0 + Real.log (1 - (1 - f A / a0)) := by
      rw [sub_sub_cancel, ← Real.log_mul h0.ne' (div_pos hA h0).ne', mul_div_cancel₀ _ h0.ne']
    rw [this]
    push_cast
    ring
  · ring
  · field_simp; ring
  · field_simp; ring
  · field_simp; ring
  · field_simp; ring

/-- `√` of a jet with positive constant term (coefficients given implicitly). -/
theorem sqrt (hf : Jet4 f a0 a1 a2 a3 a4) {s0 s1 s2 s3 s4 : ℝ} (hs0 : 0 < s0) (e0 : s0 ^ 2 = a0)
    (e1 : 2 * s0 * s1 = a1) (e2 : 2 * s0 * s2 + s1 ^ 2 = a2) (e3 : 2 * s0 * s3 + 2 * s1 * s2 = a3)
    (e4 : 2 * s0 * s4 + 2 * s1 * s3 + s2 ^ 2 = a4) :
    Jet4 (fun A => √(f A)) s0 s1 s2 s3 s4 := by
  set p : ℝ → ℝ := fun A => s0 + s1 * A + s2 * A ^ 2 + s3 * A ^ 3 + s4 * A ^ 4 with hp
  have hpp := (poly s0 s1 s2 s3 s4).mul (poly s0 s1 s2 s3 s4)
  have hz : (fun A => f A - p A * p A) =O[𝓝 0] (fun A => A ^ 5) :=
    zero_bigO ((hf.sub hpp).of_coeff (by rw [← e0]; ring) (by rw [← e1]; ring) (by rw [← e2]; ring)
      (by rw [← e3]; ring) (by rw [← e4]; ring))
  have hft := hf.tendsto
  have hpt : Tendsto p (𝓝 0) (𝓝 s0) := (poly s0 s1 s2 s3 s4).tendsto
  have hden : Tendsto (fun A => √(f A) + p A) (𝓝 0) (𝓝 (2 * s0)) := by
    have h1 := (Real.continuous_sqrt.tendsto a0).comp hft
    have := h1.add hpt
    rwa [← e0, Real.sqrt_sq hs0.le, ← two_mul] at this
  have hinv : (fun A => (√(f A) + p A)⁻¹) =O[𝓝 0] (fun _ => (1 : ℝ)) :=
    (hden.inv₀ (by positivity)).isBigO_one ℝ
  have := hz.mul hinv
  simp only [mul_one] at this
  unfold Jet4
  refine this.congr' ?_ (Eventually.of_forall fun _ => rfl)
  have hfpos : ∀ᶠ A in 𝓝 0, 0 < f A := hf.eventually_pos (by rw [← e0]; positivity)
  have hdpos : ∀ᶠ A in 𝓝 0, 0 < √(f A) + p A := hden.eventually (lt_mem_nhds (by positivity))
  filter_upwards [hfpos, hdpos] with A hA hd
  have hsq : √(f A) ^ 2 = f A := Real.sq_sqrt hA.le
  rw [mul_inv_eq_iff_eq_mul₀ hd.ne']
  simp only [hp]
  linear_combination (-1 : ℝ) * hsq

/-- Inverse of a jet with nonzero constant term (coefficients given implicitly). -/
theorem inv (hf : Jet4 f a0 a1 a2 a3 a4) (ha0 : a0 ≠ 0) {i0 i1 i2 i3 i4 : ℝ} (e0 : a0 * i0 = 1)
    (e1 : a0 * i1 + a1 * i0 = 0) (e2 : a0 * i2 + a1 * i1 + a2 * i0 = 0)
    (e3 : a0 * i3 + a1 * i2 + a2 * i1 + a3 * i0 = 0)
    (e4 : a0 * i4 + a1 * i3 + a2 * i2 + a3 * i1 + a4 * i0 = 0) :
    Jet4 (fun A => (f A)⁻¹) i0 i1 i2 i3 i4 := by
  set p : ℝ → ℝ := fun A => i0 + i1 * A + i2 * A ^ 2 + i3 * A ^ 3 + i4 * A ^ 4 with hp
  have hfp := hf.mul (poly i0 i1 i2 i3 i4)
  have hz : (fun A => f A * p A - 1) =O[𝓝 0] (fun A => A ^ 5) :=
    zero_bigO ((hfp.sub (const 1)).of_coeff (by rw [e0]; ring) (by linear_combination e1)
      (by linear_combination e2) (by linear_combination e3) (by linear_combination e4))
  have hinv : (fun A => (f A)⁻¹) =O[𝓝 0] (fun _ => (1 : ℝ)) :=
    (hf.tendsto.inv₀ ha0).isBigO_one ℝ
  have := (hz.mul hinv).neg_left
  simp only [mul_one] at this
  unfold Jet4
  refine this.congr' ?_ (Eventually.of_forall fun _ => rfl)
  have hne : ∀ᶠ A in 𝓝 0, f A ≠ 0 := hf.tendsto.eventually_ne ha0
  filter_upwards [hne] with A hA
  simp only [hp]
  field_simp
  ring

end Jet4

/-! ### The angular jet -/

theorem cos_sub_quad_bigO : (fun t : ℝ => Real.cos t - 1 + t ^ 2 / 2) =O[𝓝 0] (fun t => t ^ 4) := by
  refine IsBigO.of_bound 1 ?_
  filter_upwards [Metric.ball_mem_nhds (0 : ℝ) one_pos] with t ht
  rw [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs] at ht
  have h := Real.cos_bound ht.le
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow, one_mul]
  calc |Real.cos t - 1 + t ^ 2 / 2| = |Real.cos t - (1 - t ^ 2 / 2)| := by ring_nf
    _ ≤ |t| ^ 4 * (5 / 96) := h
    _ ≤ |t| ^ 4 := by nlinarith [pow_nonneg (abs_nonneg t) 4]

theorem sin_sub_self_bigO' : (fun t : ℝ => Real.sin t - t) =O[𝓝 0] (fun t => t ^ 3) := by
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

/-- `1 − cos(θ + A²d) = (1 − cos θ) + (d sin θ) A² + (d² cos θ/2) A⁴ + O(A⁵)`. -/
theorem jet_one_sub_cos (θ d : ℝ) :
    Jet4 (fun A => 1 - Real.cos (θ + A ^ 2 * d)) (1 - Real.cos θ) 0 (Real.sin θ * d) 0
      (Real.cos θ * d ^ 2 / 2) := by
  have hx : (fun A : ℝ => A ^ 2 * d) =O[𝓝 0] (fun A => A ^ 2) := by
    have : (fun A : ℝ => A ^ 2 * d) =O[𝓝 0] (fun A => A ^ 2) := mul_cont continuousAt_const
    exact this
  have hc := bigO_weaken (by norm_num : 5 ≤ 2 * 4) (comp_bigO (by norm_num) cos_sub_quad_bigO hx)
  have hs := bigO_weaken (by norm_num : 5 ≤ 2 * 3) (comp_bigO (by norm_num) sin_sub_self_bigO' hx)
  unfold Jet4
  refine ((hs.const_mul_left (Real.sin θ)).sub (hc.const_mul_left (Real.cos θ))).congr_left
    (fun A => ?_)
  simp only [Real.cos_add]
  ring

end Havelock.Conj13.Odd4
