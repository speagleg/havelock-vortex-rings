/-
Havelock — Conjecture 13: the finite-dimensional analysis behind the nonlinear stability argument
(paper `conj13/paper/main.tex`, Section "Nonlinear stability for a degenerate minimum",
`sec:nonlinear`, Lemmas B1–B3), stated in chart coordinates.

Coordinates (paper): `w = (a, b) ∈ A × B` with `A = K₀` (the Hessian kernel) and `B = E` (its
complement) on `M_μ`; `j = J₃ − μ₃` the momentum coordinate on `S/K`.

* Tools: `quadratic_bound` (a `C¹` map vanishing to first order at `0` with derivative
  differentiable at `0` is `O(‖w‖²)`, with explicit constants), `quadratic_bound_deriv` (its
  one-variable form), `taylor_lower` / `taylor_lower_abs` (second-order Taylor lower bound from a
  lower bound on the second derivative).
* Lemma B1: `quartic_form_coercive` (a continuous positive-definite quartic form is `≥ m‖a‖⁴`),
  `quartic_lower` (`g = q₄ + o(‖a‖⁴)` ⇒ `g ≥ (m/2)‖a‖⁴`), `slave_taylor_of_line` (Taylor in `b`
  about the slaving point), `B1_combine` (the final inequality), and `lemma_B1` (assembled).
* Lemma B2: `B2_absorb` (the AM–GM absorption step) and `lemma_B2` (assembled from the Taylor
  bound in `j`).
* Lemma B3: `lemma_B3` (the penalty function `V = H − ω₀J₃ + K_c(J₃−μ₃)²` has a strict local
  minimum at the origin of the `(j, w)` chart for `K_c` large).

What is *not* here (see the report): the existence of the chart (`S/K` a manifold, the slaving map
`β` from the implicit function theorem, the identification `{j = 0} = M_μ`), and the Taylor
expansion of `h` itself; these enter as hypotheses in the form stated in the paper.
-/
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Comp
import Mathlib.Analysis.Normed.Module.FiniteDimension

open Metric Set Filter Topology Asymptotics

namespace Havelock.Conj13.Lyapunov

/-! ### Tools -/

section Tools

variable {E F : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [NormedAddCommGroup F]
  [NormedSpace ℝ F]

/-- A map with `f 0 = 0`, derivative `f'` near `0`, `f' 0 = 0` and `f'` differentiable at `0`
satisfies `‖f w‖ ≤ L ‖w‖²` near `0`. (Used for the slaving map `β`, for `∂ⱼΦ(0, ·)` and for `ψ`.) -/
theorem quadratic_bound {f : E → F} {f' : E → E →L[ℝ] F} {ρ : ℝ} (hρ : 0 < ρ)
    (hf : ∀ w, ‖w‖ < ρ → HasFDerivAt f (f' w) w) (hf0 : f 0 = 0) (hf'0 : f' 0 = 0)
    (hf'd : DifferentiableAt ℝ f' 0) :
    ∃ L > 0, ∃ ρ' > 0, ∀ w, ‖w‖ < ρ' → ‖f w‖ ≤ L * ‖w‖ ^ 2 := by
  obtain ⟨L, hL, hLw⟩ := hf'd.isBigO_sub.exists_pos
  have hev := hLw.bound
  simp only [hf'0, sub_zero] at hev
  obtain ⟨ε, hε, hεb⟩ := Metric.eventually_nhds_iff.1 hev
  refine ⟨L, hL, min ρ ε, lt_min hρ hε, fun w hw => ?_⟩
  have hball : ∀ v ∈ closedBall (0 : E) ‖w‖, ‖v‖ < min ρ ε := fun v hv =>
    lt_of_le_of_lt (by simpa using hv) hw
  have key := (convex_closedBall (0 : E) ‖w‖).norm_image_sub_le_of_norm_hasFDerivWithin_le
    (f := f) (f' := f') (C := L * ‖w‖)
    (fun v hv => (hf v ((hball v hv).trans_le (min_le_left _ _))).hasFDerivWithinAt)
    (fun v hv => by
      have h1 : ‖f' v‖ ≤ L * ‖v‖ := hεb (by
        simpa [dist_zero_right] using (hball v hv).trans_le (min_le_right _ _))
      have h2 : ‖v‖ ≤ ‖w‖ := by simpa using hv
      calc ‖f' v‖ ≤ L * ‖v‖ := h1
        _ ≤ L * ‖w‖ := by gcongr)
    (mem_closedBall_self (norm_nonneg w)) (by simp : w ∈ closedBall (0 : E) ‖w‖)
  simp only [hf0, sub_zero] at key
  calc ‖f w‖ ≤ L * ‖w‖ * ‖w‖ := key
    _ = L * ‖w‖ ^ 2 := by ring

/-- One-variable form of `quadratic_bound`: if `ψ 0 = 0`, `ψ' = d` near `0`, `d 0 = 0` and `d` is
differentiable at `0`, then `|ψ j| ≤ L j²` near `0`. -/
theorem quadratic_bound_deriv {ψ d : ℝ → ℝ} {ρ : ℝ} (hρ : 0 < ρ)
    (hψ : ∀ j, |j| < ρ → HasDerivAt ψ (d j) j) (hψ0 : ψ 0 = 0) (hd0 : d 0 = 0)
    (hdd : DifferentiableAt ℝ d 0) :
    ∃ L > 0, ∃ ρ' > 0, ∀ j, |j| < ρ' → |ψ j| ≤ L * j ^ 2 := by
  obtain ⟨L, hL, hLw⟩ := hdd.isBigO_sub.exists_pos
  have hev := hLw.bound
  simp only [hd0, sub_zero] at hev
  obtain ⟨ε, hε, hεb⟩ := Metric.eventually_nhds_iff.1 hev
  refine ⟨L, hL, min ρ ε, lt_min hρ hε, fun j hj => ?_⟩
  have hball : ∀ v ∈ closedBall (0 : ℝ) |j|, |v| < min ρ ε := fun v hv =>
    lt_of_le_of_lt (by simpa using hv) hj
  have key := (convex_closedBall (0 : ℝ) |j|).norm_image_sub_le_of_norm_hasDerivWithin_le
    (f := ψ) (f' := d) (C := L * |j|)
    (fun v hv => (hψ v ((hball v hv).trans_le (min_le_left _ _))).hasDerivWithinAt)
    (fun v hv => by
      have h1 : ‖d v‖ ≤ L * ‖v‖ := hεb (by
        simpa [dist_zero_right] using (hball v hv).trans_le (min_le_right _ _))
      have h2 : |v| ≤ |j| := by simpa using hv
      rw [Real.norm_eq_abs] at h1
      calc ‖d v‖ ≤ L * |v| := h1
        _ ≤ L * |j| := by gcongr)
    (mem_closedBall_self (abs_nonneg j)) (by simp : j ∈ closedBall (0 : ℝ) |j|)
  simp only [hψ0, sub_zero, Real.norm_eq_abs] at key
  calc |ψ j| ≤ L * |j| * |j| := key
    _ = L * j ^ 2 := by rw [mul_assoc, ← sq, sq_abs]

/-- Second-order Taylor lower bound on `[0, T]`: if `φ'' ≥ c` there, then
`φ T ≥ φ 0 + φ'(0) T + (c/2) T²`. -/
theorem taylor_lower {φ φ' φ'' : ℝ → ℝ} {T c : ℝ} (hT : 0 ≤ T)
    (h1 : ∀ s ∈ Icc 0 T, HasDerivAt φ (φ' s) s) (h2 : ∀ s ∈ Icc 0 T, HasDerivAt φ' (φ'' s) s)
    (hc : ∀ s ∈ Icc 0 T, c ≤ φ'' s) : φ 0 + φ' 0 * T + c / 2 * T ^ 2 ≤ φ T := by
  -- `χ' = φ' - φ'(0) - c s` is monotone, hence `≥ 0`
  have hd1 : ∀ s ∈ Icc 0 T,
      HasDerivAt (fun s => φ' s - φ' 0 - c * s) (φ'' s - c) s := fun s hs => by
    exact (((h2 s hs).sub_const (φ' 0)).sub ((hasDerivAt_id' s).const_mul c)).congr_deriv
      (by ring)
  have hmono1 : MonotoneOn (fun s => φ' s - φ' 0 - c * s) (Icc 0 T) :=
    monotoneOn_of_hasDerivWithinAt_nonneg (convex_Icc 0 T)
      (fun s hs => (hd1 s hs).continuousAt.continuousWithinAt)
      (fun s hs => (hd1 s (interior_subset hs)).hasDerivWithinAt)
      (fun s hs => by have := hc s (interior_subset hs); linarith)
  have hpos : ∀ s ∈ Icc 0 T, 0 ≤ φ' s - φ' 0 - c * s := fun s hs => by
    have := hmono1 (left_mem_Icc.2 hT) hs hs.1
    simpa using this
  have hd0 : ∀ s ∈ Icc 0 T, HasDerivAt (fun s => φ s - φ 0 - φ' 0 * s - c / 2 * s ^ 2)
      (φ' s - φ' 0 - c * s) s := fun s hs => by
    have := (((h1 s hs).sub_const (φ 0)).sub ((hasDerivAt_id' s).const_mul (φ' 0))).sub
      ((hasDerivAt_pow 2 s).const_mul (c / 2))
    convert this using 1
    simp; ring
  have hmono0 : MonotoneOn (fun s => φ s - φ 0 - φ' 0 * s - c / 2 * s ^ 2) (Icc 0 T) :=
    monotoneOn_of_hasDerivWithinAt_nonneg (convex_Icc 0 T)
      (fun s hs => (hd0 s hs).continuousAt.continuousWithinAt)
      (fun s hs => (hd0 s (interior_subset hs)).hasDerivWithinAt)
      (fun s hs => hpos s (interior_subset hs))
  have := hmono0 (left_mem_Icc.2 hT) (right_mem_Icc.2 hT) hT
  simp only at this
  linarith

/-- Two-sided form of `taylor_lower`: for `|j| ≤ T`, `φ j ≥ φ 0 + φ'(0) j + (c/2) j²` when
`φ'' ≥ c` on `[-T, T]`. -/
theorem taylor_lower_abs {φ φ' φ'' : ℝ → ℝ} {T c : ℝ}
    (h1 : ∀ s, |s| ≤ T → HasDerivAt φ (φ' s) s) (h2 : ∀ s, |s| ≤ T → HasDerivAt φ' (φ'' s) s)
    (hc : ∀ s, |s| ≤ T → c ≤ φ'' s) (j : ℝ) (hj : |j| ≤ T) :
    φ 0 + φ' 0 * j + c / 2 * j ^ 2 ≤ φ j := by
  rcases le_or_gt 0 j with hj0 | hj0
  · have hsub : ∀ s ∈ Icc 0 j, |s| ≤ T := fun s hs => by
      rw [abs_of_nonneg hs.1]; exact hs.2.trans (le_abs_self j |>.trans hj)
    exact taylor_lower hj0 (fun s hs => h1 s (hsub s hs)) (fun s hs => h2 s (hsub s hs))
      (fun s hs => hc s (hsub s hs))
  · have hsub : ∀ s ∈ Icc 0 (-j), |-s| ≤ T := fun s hs => by
      rw [abs_neg, abs_of_nonneg hs.1]; exact hs.2.trans (neg_le_abs j |>.trans hj)
    have := taylor_lower (φ := fun s => φ (-s)) (φ' := fun s => -φ' (-s))
      (φ'' := fun s => φ'' (-s)) (T := -j) (c := c) (by linarith)
      (fun s hs => by
        exact ((h1 (-s) (hsub s hs)).comp s (hasDerivAt_neg' s)).congr_deriv (by ring))
      (fun s hs => by
        exact ((h2 (-s) (hsub s hs)).comp s (hasDerivAt_neg' s)).neg.congr_deriv (by ring))
      (fun s hs => hc (-s) (hsub s hs))
    simp only [neg_zero, neg_neg] at this
    have e : -φ' 0 * -j = φ' 0 * j := by ring
    rw [e, neg_sq] at this
    exact this

end Tools

/-! ### Lemma B1: quartic lower bound on `M_μ` -/

section B1

variable {A B : Type*} [NormedAddCommGroup A] [NormedSpace ℝ A] [NormedAddCommGroup B]
  [NormedSpace ℝ B]

/-- A continuous, `4`-homogeneous function on a finite-dimensional space, positive away from `0`
(the positive definite quartic Taylor term `q₄` of (A4)), satisfies `q₄ a ≥ m ‖a‖⁴` for some
`m > 0` (`m = min_{|a|=1} q₄`). -/
theorem quartic_form_coercive [FiniteDimensional ℝ A] {q : A → ℝ} (hqc : Continuous q)
    (hhom : ∀ (t : ℝ) (a : A), q (t • a) = t ^ 4 * q a) (hpos : ∀ a, a ≠ 0 → 0 < q a) :
    ∃ m > 0, ∀ a, m * ‖a‖ ^ 4 ≤ q a := by
  have hq0 : q 0 = 0 := by simpa using hhom 0 0
  have hunit : ∀ a : A, a ≠ 0 → ‖a‖⁻¹ • a ∈ sphere (0 : A) 1 := fun a ha => by
    simp [norm_smul, norm_ne_zero_iff.2 ha]
  have hscale : ∀ a : A, a ≠ 0 → q a = ‖a‖ ^ 4 * q (‖a‖⁻¹ • a) := fun a ha => by
    rw [← hhom, smul_smul, mul_inv_cancel₀ (norm_ne_zero_iff.2 ha), one_smul]
  rcases (sphere (0 : A) 1).eq_empty_or_nonempty with h | h
  · refine ⟨1, one_pos, fun a => ?_⟩
    by_cases ha : a = 0
    · simp [ha, hq0]
    · have := hunit a ha; rw [h] at this; exact this.elim
  · obtain ⟨u, hu, humin⟩ := (isCompact_sphere (0 : A) 1).exists_isMinOn h hqc.continuousOn
    have hu0 : u ≠ 0 := by
      intro h0; simp [h0] at hu
    refine ⟨q u, hpos u hu0, fun a => ?_⟩
    by_cases ha : a = 0
    · simp [ha, hq0]
    · rw [hscale a ha]
      have h1 : q u ≤ q (‖a‖⁻¹ • a) := humin (hunit a ha)
      have h2 : 0 ≤ ‖a‖ ^ 4 := by positivity
      nlinarith

omit [NormedSpace ℝ A] in
/-- (A4) ⇒ quartic lower bound for `g`: if `g = q₄ + o(‖a‖⁴)` and `q₄ ≥ m‖a‖⁴`, then
`g a ≥ (m/2) ‖a‖⁴` near `0` (paper Lemma B1, "`g(a) ≥ ½ min_{|a|=1} q₄ |a|⁴` for small `a`"). -/
theorem quartic_lower {g q : A → ℝ} {m : ℝ} (hm : 0 < m) (hq : ∀ a, m * ‖a‖ ^ 4 ≤ q a)
    (hgq : (fun a => g a - q a) =o[𝓝 0] (fun a => ‖a‖ ^ 4)) :
    ∃ ρ > 0, ∀ a, ‖a‖ < ρ → m / 2 * ‖a‖ ^ 4 ≤ g a := by
  obtain ⟨ε, hε, hεb⟩ := Metric.eventually_nhds_iff.1 (hgq.def (half_pos hm))
  refine ⟨ε, hε, fun a ha => ?_⟩
  have h1 := hεb (y := a) (by simpa [dist_zero_right] using ha)
  simp only [Real.norm_eq_abs, abs_pow, abs_norm] at h1
  have h2 := hq a
  have h3 := neg_abs_le (g a - q a)
  linarith

omit [NormedAddCommGroup A] [NormedSpace ℝ A] in
/-- Taylor's theorem in `b` about the slaving point (paper Lemma B1): along the segment
`t ↦ f a (β a + t • (b - β a))`, if the first derivative vanishes at `t = 0`
(`∂_b f(a, β a) = 0`) and the second derivative is `≥ 2 c₁ ‖b − β a‖²` (uniform positivity of
`∂²_b f` near `0`), then `f a b ≥ f a (β a) + c₁ ‖b − β a‖²`. -/
theorem slave_taylor_of_line (f : A → B → ℝ) (β : A → B) (a : A) (b : B) (c₁ : ℝ)
    (φ' φ'' : ℝ → ℝ)
    (h1 : ∀ t ∈ Icc (0 : ℝ) 1, HasDerivAt (fun t => f a (β a + t • (b - β a))) (φ' t) t)
    (h2 : ∀ t ∈ Icc (0 : ℝ) 1, HasDerivAt φ' (φ'' t) t) (hcrit : φ' 0 = 0)
    (hconv : ∀ t ∈ Icc (0 : ℝ) 1, 2 * c₁ * ‖b - β a‖ ^ 2 ≤ φ'' t) :
    f a (β a) + c₁ * ‖b - β a‖ ^ 2 ≤ f a b := by
  have := taylor_lower zero_le_one h1 h2 hconv
  have e : β a + (b - β a) = b := by abel
  simp only [zero_smul, add_zero, one_smul, hcrit, one_pow, mul_one, e] at this
  linarith

omit [NormedSpace ℝ A] [NormedSpace ℝ B] in
/-- The final step of paper Lemma B1: from `f(a,b) − g(a) ≥ c₁‖b − β(a)‖²`, `‖β(a)‖ ≤ L‖a‖²` and
`g(a) ≥ c₂‖a‖⁴` we get `f(a,b) ≥ κ(‖a‖⁴ + ‖b‖²)` with `κ = min(c₁/2, c₂/(1+2L²))`. -/
theorem B1_combine {f : A → B → ℝ} {g : A → ℝ} {β : A → B} {c₁ c₂ L ρ : ℝ} (hc₁ : 0 < c₁)
    (hc₂ : 0 < c₂)
    (h1 : ∀ a b, ‖a‖ < ρ → ‖b‖ < ρ → g a + c₁ * ‖b - β a‖ ^ 2 ≤ f a b)
    (h2 : ∀ a, ‖a‖ < ρ → ‖β a‖ ≤ L * ‖a‖ ^ 2)
    (h3 : ∀ a, ‖a‖ < ρ → c₂ * ‖a‖ ^ 4 ≤ g a) :
    ∀ a b, ‖a‖ < ρ → ‖b‖ < ρ →
      min (c₁ / 2) (c₂ / (1 + 2 * L ^ 2)) * (‖a‖ ^ 4 + ‖b‖ ^ 2) ≤ f a b := by
  intro a b ha hb
  set κ := min (c₁ / 2) (c₂ / (1 + 2 * L ^ 2))
  set s := ‖b - β a‖
  have hD : 0 < 1 + 2 * L ^ 2 := by positivity
  have hκ1 : κ ≤ c₁ / 2 := min_le_left _ _
  have hκ2 : κ * (1 + 2 * L ^ 2) ≤ c₂ := (le_div_iff₀ hD).1 (min_le_right _ _)
  have hκ0 : 0 ≤ κ := le_min (by positivity) (by positivity)
  have htri : ‖b‖ ≤ s + L * ‖a‖ ^ 2 := by
    calc ‖b‖ = ‖(b - β a) + β a‖ := by rw [sub_add_cancel]
      _ ≤ s + ‖β a‖ := norm_add_le _ _
      _ ≤ s + L * ‖a‖ ^ 2 := by gcongr; exact h2 a ha
  have hb2 : ‖b‖ ^ 2 ≤ 2 * s ^ 2 + 2 * L ^ 2 * ‖a‖ ^ 4 := by
    have hs : 0 ≤ s := norm_nonneg _
    have := pow_le_pow_left₀ (norm_nonneg b) htri 2
    nlinarith [sq_nonneg (s - L * ‖a‖ ^ 2)]
  have hf := h1 a b ha hb
  have hg := h3 a ha
  have e1 : κ * ‖b‖ ^ 2 ≤ κ * (2 * s ^ 2 + 2 * L ^ 2 * ‖a‖ ^ 4) :=
    mul_le_mul_of_nonneg_left hb2 hκ0
  have e2 : κ * (2 * s ^ 2) ≤ c₁ * s ^ 2 := by nlinarith [sq_nonneg s]
  have e3 : κ * (1 + 2 * L ^ 2) * ‖a‖ ^ 4 ≤ c₂ * ‖a‖ ^ 4 :=
    mul_le_mul_of_nonneg_right hκ2 (by positivity)
  nlinarith

/-- **Paper Lemma B1 (quartic lower bound on `M_μ`)**, in chart coordinates `w = (a, b)`,
`a ∈ A = K₀` (finite-dimensional), `b ∈ B = E`, with `f = h − h(0)`. Hypotheses:
* the slaving map `β` (from the implicit function theorem) is differentiable near `0` with
  `β 0 = 0`, `dβ(0) = 0` (the cross block of `d²f(0)` vanishes) and `dβ` differentiable at `0`;
* Taylor in `b` about `β(a)`: `f(a,b) ≥ f(a,β(a)) + c₁‖b − β(a)‖²` (see `slave_taylor_of_line`);
* (A4): `g(a) = f(a, β(a)) = q₄(a) + o(‖a‖⁴)` with `q₄` a continuous `4`-homogeneous form,
  positive away from `0`.
Conclusion: `f(a,b) ≥ κ(‖a‖⁴ + ‖b‖²)` for `‖a‖, ‖b‖ < ρ'`. -/
theorem lemma_B1 [FiniteDimensional ℝ A] (f : A → B → ℝ) (β : A → B) (β' : A → A →L[ℝ] B)
    (q : A → ℝ) {c₁ ρ : ℝ} (hc₁ : 0 < c₁) (hρ : 0 < ρ)
    (hβ : ∀ a, ‖a‖ < ρ → HasFDerivAt β (β' a) a) (hβ0 : β 0 = 0) (hβ'0 : β' 0 = 0)
    (hβ'd : DifferentiableAt ℝ β' 0)
    (hslave : ∀ a b, ‖a‖ < ρ → ‖b‖ < ρ → f a (β a) + c₁ * ‖b - β a‖ ^ 2 ≤ f a b)
    (hqc : Continuous q) (hhom : ∀ (t : ℝ) (a : A), q (t • a) = t ^ 4 * q a)
    (hpos : ∀ a, a ≠ 0 → 0 < q a)
    (hg : (fun a => f a (β a) - q a) =o[𝓝 0] (fun a => ‖a‖ ^ 4)) :
    ∃ κ > 0, ∃ ρ' > 0, ∀ a b, ‖a‖ < ρ' → ‖b‖ < ρ' → κ * (‖a‖ ^ 4 + ‖b‖ ^ 2) ≤ f a b := by
  obtain ⟨m, hm, hqm⟩ := quartic_form_coercive hqc hhom hpos
  obtain ⟨ρ₁, hρ₁, hg4⟩ := quartic_lower hm hqm hg
  obtain ⟨L, hL, ρ₂, hρ₂, hβL⟩ := quadratic_bound hρ hβ hβ0 hβ'0 hβ'd
  set ρ' := min ρ (min ρ₁ ρ₂)
  have hρ'1 : ρ' ≤ ρ := min_le_left _ _
  have hρ'2 : ρ' ≤ ρ₁ := (min_le_right _ _).trans (min_le_left _ _)
  have hρ'3 : ρ' ≤ ρ₂ := (min_le_right _ _).trans (min_le_right _ _)
  refine ⟨min (c₁ / 2) (m / 2 / (1 + 2 * L ^ 2)), lt_min (by positivity) (by positivity), ρ',
    lt_min hρ (lt_min hρ₁ hρ₂), ?_⟩
  exact B1_combine (g := fun a => f a (β a)) hc₁ (half_pos hm)
    (fun a b ha hb => hslave a b (ha.trans_le hρ'1) (hb.trans_le hρ'1))
    (fun a ha => hβL a (ha.trans_le hρ'3)) (fun a ha => hg4 a (ha.trans_le hρ'2))

end B1

/-! ### Lemma B2: uniformity in nearby momentum -/

/-- The absorption step of paper Lemma B2 (AM–GM on the `j‖a‖²` term, direct absorption of the
`j‖b‖²` term once `L'|j| ≤ κ/2`). Here `A = ‖a‖`, `Bn = ‖b‖`, `Φ₀ = Φ(0, w)`, `Φ = Φ(j, w)`. -/
theorem B2_absorb {κ L' C₂ A Bn j Φ Φ₀ : ℝ} (hκ : 0 < κ)
    (hj : L' * |j| ≤ κ / 2) (h0 : κ * (A ^ 4 + Bn ^ 2) ≤ Φ₀)
    (hT : Φ₀ - L' * |j| * (A ^ 2 + Bn ^ 2) - C₂ * j ^ 2 ≤ Φ) :
    κ / 2 * (A ^ 4 + Bn ^ 2) - (L' ^ 2 / (2 * κ) + C₂) * j ^ 2 ≤ Φ := by
  -- AM–GM: `L'|j|A² ≤ (κ/2)A⁴ + (L'²/(2κ)) j²`
  have amgm : L' * |j| * A ^ 2 ≤ κ / 2 * A ^ 4 + L' ^ 2 / (2 * κ) * j ^ 2 := by
    have e : κ / 2 * A ^ 4 + L' ^ 2 / (2 * κ) * j ^ 2
        = (κ ^ 2 * A ^ 4 + L' ^ 2 * j ^ 2) / (2 * κ) := by field_simp
    rw [e, le_div_iff₀ (by positivity)]
    nlinarith [sq_nonneg (κ * A ^ 2 - L' * |j|), sq_abs j]
  have absorb : L' * |j| * Bn ^ 2 ≤ κ / 2 * Bn ^ 2 :=
    mul_le_mul_of_nonneg_right hj (sq_nonneg Bn)
  nlinarith

/-- **Paper Lemma B2 (uniformity in nearby momentum)**, in chart coordinates `(j, w)`,
`w = (a, b) ∈ A × B` (product norm), for `Φ(j, w) = H(j, w) − H(j, 0)`. Hypotheses:
* (Lemma B1) `Φ(0, w) ≥ κ(‖a‖⁴ + ‖b‖²)` near `0`;
* `j ↦ Φ(j, w)` has derivatives `D₁(j, w)`, `D₂(j, w)` near `0` with `D₂ ≥ −2C₂`;
* `w ↦ D₁(0, w) = ∂ⱼΦ(0, w)` vanishes at `0` with vanishing derivative there
  (`∂ⱼΦ(0,0) = 0`, `∂_w∂ⱼΦ(0,0) = 0`), its derivative being differentiable at `0`.
Conclusion: `Φ(j, w) ≥ (κ/2)(‖a‖⁴ + ‖b‖²) − C j²` near `0`. -/
theorem lemma_B2 {A B : Type*} [NormedAddCommGroup A] [NormedSpace ℝ A] [NormedAddCommGroup B]
    [NormedSpace ℝ B] (Φ D₁ D₂ : ℝ → A × B → ℝ) (D₁' : A × B → A × B →L[ℝ] ℝ)
    {κ C₂ ρ : ℝ} (hκ : 0 < κ) (hρ : 0 < ρ)
    (hB1 : ∀ a b, ‖a‖ < ρ → ‖b‖ < ρ → κ * (‖a‖ ^ 4 + ‖b‖ ^ 2) ≤ Φ 0 (a, b))
    (hD₁ : ∀ w : A × B, ‖w‖ < ρ → ∀ j, |j| ≤ ρ → HasDerivAt (fun j => Φ j w) (D₁ j w) j)
    (hD₂ : ∀ w : A × B, ‖w‖ < ρ → ∀ j, |j| ≤ ρ → HasDerivAt (fun j => D₁ j w) (D₂ j w) j)
    (hC₂ : ∀ w : A × B, ‖w‖ < ρ → ∀ j, |j| ≤ ρ → -(2 * C₂) ≤ D₂ j w)
    (hD₁w : ∀ w : A × B, ‖w‖ < ρ → HasFDerivAt (D₁ 0) (D₁' w) w) (hD₁0 : D₁ 0 0 = 0)
    (hD₁'0 : D₁' 0 = 0) (hD₁'d : DifferentiableAt ℝ D₁' 0) :
    ∃ C, ∃ ρ' > 0, ∀ j (a : A) (b : B), |j| < ρ' → ‖a‖ < ρ' → ‖b‖ < ρ' →
      κ / 2 * (‖a‖ ^ 4 + ‖b‖ ^ 2) - C * j ^ 2 ≤ Φ j (a, b) := by
  obtain ⟨L, hL, ρ₁, hρ₁, hDL⟩ := quadratic_bound hρ hD₁w hD₁0 hD₁'0 hD₁'d
  set ρ' := min ρ (min ρ₁ (κ / (2 * L)))
  have hρ'1 : ρ' ≤ ρ := min_le_left _ _
  have hρ'2 : ρ' ≤ ρ₁ := (min_le_right _ _).trans (min_le_left _ _)
  have hρ'3 : ρ' ≤ κ / (2 * L) := (min_le_right _ _).trans (min_le_right _ _)
  refine ⟨L ^ 2 / (2 * κ) + C₂, ρ', lt_min hρ (lt_min hρ₁ (by positivity)),
    fun j a b hj ha hb => ?_⟩
  set w : A × B := (a, b)
  have hw : ‖w‖ < ρ' := by
    rw [Prod.norm_def]; exact max_lt ha hb
  have hw2 : ‖w‖ ^ 2 ≤ ‖a‖ ^ 2 + ‖b‖ ^ 2 := by
    rw [Prod.norm_def]
    rcases le_total ‖a‖ ‖b‖ with h | h
    · rw [max_eq_right h]; nlinarith [sq_nonneg ‖a‖]
    · rw [max_eq_left h]; nlinarith [sq_nonneg ‖b‖]
  have hwρ : ‖w‖ < ρ := hw.trans_le hρ'1
  -- Taylor in `j`
  have hT := taylor_lower_abs (φ := fun j => Φ j w) (φ' := fun j => D₁ j w)
    (φ'' := fun j => D₂ j w) (T := ρ') (c := -(2 * C₂))
    (fun s hs => hD₁ w hwρ s (hs.trans hρ'1)) (fun s hs => hD₂ w hwρ s (hs.trans hρ'1))
    (fun s hs => hC₂ w hwρ s (hs.trans hρ'1)) j hj.le
  -- `|∂ⱼΦ(0, w)| ≤ L ‖w‖² ≤ L (‖a‖² + ‖b‖²)`
  have hD := hDL w (hw.trans_le hρ'2)
  rw [Real.norm_eq_abs] at hD
  have hlin : -(L * |j| * (‖a‖ ^ 2 + ‖b‖ ^ 2)) ≤ D₁ 0 w * j := by
    have h1 : |D₁ 0 w * j| ≤ L * |j| * (‖a‖ ^ 2 + ‖b‖ ^ 2) := by
      rw [abs_mul]
      calc |D₁ 0 w| * |j| ≤ L * ‖w‖ ^ 2 * |j| := by gcongr
        _ ≤ L * (‖a‖ ^ 2 + ‖b‖ ^ 2) * |j| := by gcongr
        _ = L * |j| * (‖a‖ ^ 2 + ‖b‖ ^ 2) := by ring
    exact (abs_le.1 h1).1
  have hj' : L * |j| ≤ κ / 2 := by
    have := hj.le.trans hρ'3
    rw [le_div_iff₀ (by positivity)] at this
    linarith
  exact B2_absorb hκ hj'
    (hB1 a b (ha.trans_le hρ'1) (hb.trans_le hρ'1)) (by linarith)

/-! ### Lemma B3: the penalty function -/

/-- **Paper Lemma B3 (conserved function with a strict minimum modulo `K`)**, in chart
coordinates `(j, w)` on `S/K` (`j = J₃ − μ₃`). `Hc j w` is `H` in these coordinates and
`V j w = Hc j w − ω₀ j + K_c j²` is `H − ω₀J₃ + K_c(J₃ − μ₃)²` up to the additive constant
`−ω₀μ₃`. Hypotheses:
* (Lemma B2) `Φ(j, w) = Hc j w − Hc j 0 ≥ (κ/2)(‖a‖⁴ + ‖b‖²) − C j²` near `0`;
* (A3) along the family of relative equilibria `w = 0`, `d/dj H(x_j) = ω(j)`, with `ω`
  differentiable at `0` and `ω(0) = ω₀`.
Conclusion: there is `K₀` such that for `K_c > K₀`, `V j w − V 0 0 ≥ (κ/2)(‖a‖⁴ + ‖b‖²) +
(K_c − K₀) j²` near `0`, and hence `V j w > V 0 0` for `(j, w) ≠ 0` near `0`. -/
theorem lemma_B3 {A B : Type*} [NormedAddCommGroup A] [NormedAddCommGroup B]
    (Hc : ℝ → A × B → ℝ) (ω : ℝ → ℝ) {ω₀ κ C ρ : ℝ} (hκ : 0 < κ) (hρ : 0 < ρ)
    (hB2 : ∀ j (a : A) (b : B), |j| < ρ → ‖a‖ < ρ → ‖b‖ < ρ →
      κ / 2 * (‖a‖ ^ 4 + ‖b‖ ^ 2) - C * j ^ 2 ≤ Hc j (a, b) - Hc j 0)
    (hfam : ∀ j, |j| < ρ → HasDerivAt (fun j => Hc j 0) (ω j) j) (hω0 : ω 0 = ω₀)
    (hωd : DifferentiableAt ℝ ω 0) :
    ∃ K₀, ∀ K_c, K₀ < K_c → ∃ ρ' > 0, ∀ j (a : A) (b : B), |j| < ρ' → ‖a‖ < ρ' → ‖b‖ < ρ' →
      κ / 2 * (‖a‖ ^ 4 + ‖b‖ ^ 2) + (K_c - K₀) * j ^ 2 ≤
        (Hc j (a, b) - ω₀ * j + K_c * j ^ 2) - (Hc 0 0 - ω₀ * 0 + K_c * 0 ^ 2) ∧
      ((j, a, b) ≠ 0 →
        Hc 0 0 - ω₀ * 0 + K_c * 0 ^ 2 < Hc j (a, b) - ω₀ * j + K_c * j ^ 2) := by
  -- `ψ(j) = H(x_j) − H(x₀) − ω₀ j` is `O(j²)`
  have hψ := quadratic_bound_deriv (ψ := fun j => Hc j 0 - Hc 0 0 - ω₀ * j)
    (d := fun j => ω j - ω₀) hρ
    (fun j hj => by
      exact (((hfam j hj).sub_const (Hc 0 0)).sub ((hasDerivAt_id' j).const_mul ω₀)).congr_deriv
        (by ring))
    (by simp) (by simp [hω0]) (hωd.sub_const ω₀)
  obtain ⟨C₃, -, ρ₃, hρ₃, hψb⟩ := hψ
  refine ⟨C + C₃, fun K_c hK => ⟨min ρ ρ₃, lt_min hρ hρ₃, fun j a b hj ha hb => ?_⟩⟩
  have hj1 : |j| < ρ := hj.trans_le (min_le_left _ _)
  have ha1 : ‖a‖ < ρ := ha.trans_le (min_le_left _ _)
  have hb1 : ‖b‖ < ρ := hb.trans_le (min_le_left _ _)
  have h2 := hB2 j a b hj1 ha1 hb1
  have h3 := hψb j (hj.trans_le (min_le_right _ _))
  have h3' := neg_abs_le (Hc j 0 - Hc 0 0 - ω₀ * j)
  have hmain : κ / 2 * (‖a‖ ^ 4 + ‖b‖ ^ 2) + (K_c - (C + C₃)) * j ^ 2 ≤
      (Hc j (a, b) - ω₀ * j + K_c * j ^ 2) - (Hc 0 0 - ω₀ * 0 + K_c * 0 ^ 2) := by
    nlinarith
  refine ⟨hmain, fun hne => ?_⟩
  have hpos : 0 < κ / 2 * (‖a‖ ^ 4 + ‖b‖ ^ 2) + (K_c - (C + C₃)) * j ^ 2 := by
    by_cases hj0 : j = 0
    · have hab : a ≠ 0 ∨ b ≠ 0 := by
        by_contra h
        push Not at h
        exact hne (by simp [hj0, h.1, h.2])
      have : 0 < ‖a‖ ^ 4 + ‖b‖ ^ 2 := by
        rcases hab with h | h
        · have := norm_pos_iff.2 h; positivity
        · have := norm_pos_iff.2 h; positivity
      subst hj0
      simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero, add_zero]
      positivity
    · have : 0 < j ^ 2 := by positivity
      have : 0 < K_c - (C + C₃) := by linarith
      have : 0 ≤ κ / 2 * (‖a‖ ^ 4 + ‖b‖ ^ 2) := by positivity
      nlinarith
  linarith

end Havelock.Conj13.Lyapunov
