import Havelock.Tanh.LargeN.Defs
import Mathlib.Analysis.Calculus.Deriv.MeanValue

/-! # Tanh model, large `N`: the energy–momentum core (note Lemma 17, analytic part)

For `|x| < 1`, `τ > 0` (and `β(x,τ) ≠ 0`):
`β G = (cosh τ + 1)² K / (32 sinh⁴ τ)`, and `K(x,τ) > 0` for `0 < |x| < 1`, `τ > 0`.

The sign follows the note: `K(x,0) = 0` and
`∂_τ K = sinh τ (1 − x²) [x²(cosh τ − 1) − (cosh xτ − 1)]`, where the bracket is positive because
its own `τ`-derivative is `x (x sinh τ − sinh xτ)`, which is positive because *its* derivative is
`x² (cosh τ − cosh xτ) > 0`. Each step uses `strictMonoOn_of_deriv_pos` on `[0, ∞)`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real Set

/-- A differentiable function vanishing at `0` with positive derivative on `(0,∞)` is positive
on `(0,∞)`. -/
theorem pos_of_hasDerivAt_pos {f f' : ℝ → ℝ} (hf : ∀ t, HasDerivAt f (f' t) t) (h0 : f 0 = 0)
    (hpos : ∀ t, 0 < t → 0 < f' t) {t : ℝ} (ht : 0 < t) : 0 < f t := by
  have hmono : StrictMonoOn f (Ici 0) :=
    strictMonoOn_of_deriv_pos (convex_Ici 0)
      (fun y _ => (hf y).continuousAt.continuousWithinAt)
      (by
        intro y hy
        rw [interior_Ici] at hy
        rw [(hf y).deriv]
        exact hpos y hy)
  have := hmono (mem_Ici.mpr le_rfl) (mem_Ici.mpr ht.le) ht
  rwa [h0] at this

/-- `x (x sinh τ − sinh xτ) > 0` for `0 < |x| < 1`, `τ > 0` (i.e. `sinh xτ < x sinh τ` for `x > 0`). -/
theorem x_mul_sinh_gap_pos {x τ : ℝ} (hx0 : 0 < |x|) (hx1 : |x| < 1) (hτ : 0 < τ) :
    0 < x * (x * sinh τ - sinh (x * τ)) := by
  have hx2 : 0 < x ^ 2 := by rw [← sq_abs]; positivity
  refine pos_of_hasDerivAt_pos (f := fun t => x * (x * sinh t - sinh (x * t)))
    (f' := fun t => x ^ 2 * (cosh t - cosh (x * t))) ?_ (by simp) ?_ hτ
  · intro t
    have h1 := ((hasDerivAt_id t).const_mul x).sinh
    have h2 := ((Real.hasDerivAt_sinh t).const_mul x).sub h1
    have h3 := h2.const_mul x
    convert h3 using 1
    · funext y; simp only [id, Pi.sub_apply]
    · simp only [id]; ring
  · intro t ht
    have : cosh (x * t) < cosh t := by
      rw [cosh_lt_cosh, abs_mul, abs_of_pos ht]
      nlinarith
    have : 0 < cosh t - cosh (x * t) := by linarith
    positivity

/-- `cosh xτ − 1 < x² (cosh τ − 1)` for `0 < |x| < 1`, `τ > 0`. -/
theorem cosh_gap_pos {x τ : ℝ} (hx0 : 0 < |x|) (hx1 : |x| < 1) (hτ : 0 < τ) :
    0 < x ^ 2 * (cosh τ - 1) - (cosh (x * τ) - 1) := by
  refine pos_of_hasDerivAt_pos (f := fun t => x ^ 2 * (cosh t - 1) - (cosh (x * t) - 1))
    (f' := fun t => x * (x * sinh t - sinh (x * t))) ?_ (by simp) ?_ hτ
  · intro t
    have h1 := ((hasDerivAt_id t).const_mul x).cosh
    have h2 := (((Real.hasDerivAt_cosh t).sub_const 1).const_mul (x ^ 2)).sub (h1.sub_const 1)
    convert h2 using 1
    · funext y; simp only [id, Pi.sub_apply]
    · simp only [id]; ring
  · intro t ht
    exact x_mul_sinh_gap_pos hx0 hx1 ht

/-- `∂_τ K(x,τ) = sinh τ (1 − x²) [x²(cosh τ − 1) − (cosh xτ − 1)]` (note Lemma 17, proof). -/
theorem hasDerivAt_K (x τ : ℝ) :
    HasDerivAt (fun t => K x t)
      (sinh τ * (1 - x ^ 2) * (x ^ 2 * (cosh τ - 1) - (cosh (x * τ) - 1))) τ := by
  have hC := ((hasDerivAt_id τ).const_mul x).cosh
  have hS := ((hasDerivAt_id τ).const_mul x).sinh
  have hc := Real.hasDerivAt_cosh τ
  have hs := Real.hasDerivAt_sinh τ
  have h := ((hc.mul ((hasDerivAt_const τ (1 : ℝ)).sub hC)).add
      ((hs.const_mul x).mul hS)).sub
      ((((hc.sub_const 1).const_mul (x ^ 2)).div_const 2).mul
        ((hasDerivAt_const τ (2 : ℝ)).sub ((hc.sub_const 1).const_mul (1 - x ^ 2))))
  convert h using 1
  · funext t
    simp only [K, id, Pi.add_apply, Pi.sub_apply, Pi.mul_apply]
  · simp only [id, Pi.sub_apply]
    ring

/-- note Lemma 17: `K(x,τ) > 0` for `0 < |x| < 1`, `τ > 0`. -/
theorem K_pos {x τ : ℝ} (hx0 : 0 < |x|) (hx1 : |x| < 1) (hτ : 0 < τ) : 0 < K x τ := by
  have hx2 : x ^ 2 < 1 := by
    rw [← sq_abs]; nlinarith [abs_nonneg x]
  refine pos_of_hasDerivAt_pos (f := fun t => K x t)
    (f' := fun t => sinh t * (1 - x ^ 2) * (x ^ 2 * (cosh t - 1) - (cosh (x * t) - 1)))
    (hasDerivAt_K x) (by simp [K]) ?_ hτ
  intro t ht
  have h1 : 0 < sinh t := sinh_pos_iff.mpr ht
  have h2 : 0 < 1 - x ^ 2 := by linarith
  have h3 := cosh_gap_pos hx0 hx1 ht
  positivity

/-- note Lemma 17, the identity: `β G = (cosh τ + 1)² K / (32 sinh⁴ τ)` for `τ > 0`
(given `β ≠ 0`, which holds for `|x| < 1` by note Lemma 12). -/
theorem beta_mul_G {x τ : ℝ} (hτ : 0 < τ) (hβ : beta x τ ≠ 0) :
    beta x τ * G x τ = (cosh τ + 1) ^ 2 * K x τ / (32 * sinh τ ^ 4) := by
  have hs : sinh τ ≠ 0 := (sinh_pos_iff.mpr hτ).ne'
  have hG : beta x τ * G x τ = beta x τ * (psi 0 τ - psi x τ) - gamma x τ ^ 2 := by
    rw [G]; field_simp
  set c := cosh τ
  set s := sinh τ
  set C := cosh (x * τ)
  set S := sinh (x * τ)
  have hc : c ^ 2 - s ^ 2 = 1 := cosh_sq_sub_sinh_sq τ
  have hC : C ^ 2 - S ^ 2 = 1 := cosh_sq_sub_sinh_sq (x * τ)
  have hB : beta x τ = ((1 - x ^ 2) * s ^ 2 + 2 - 2 * (c * C - x * s * S)) / (8 * s ^ 2) := by
    simp only [beta, Q]; field_simp; ring
  have hA : psi 0 τ - psi x τ = (x ^ 2 * s ^ 2 + 2 * c - 2 * (c * C - x * s * S)) / (8 * s ^ 2) := by
    simp only [psi, Q, zero_mul, cosh_zero, sinh_zero]; field_simp; ring
  have hR : gamma x τ = (c * S - x * s * C) / (4 * s ^ 2) := by
    rfl
  have hK : K x τ = c * (1 - C) + x * s * S - x ^ 2 * (c - 1) / 2 * (2 - (1 - x ^ 2) * (c - 1)) := by
    rfl
  have key : ((1 - x ^ 2) * s ^ 2 + 2 - 2 * (c * C - x * s * S))
        * (x ^ 2 * s ^ 2 + 2 * c - 2 * (c * C - x * s * S)) - 4 * (c * S - x * s * C) ^ 2
      = 2 * (c + 1) ^ 2
        * (c * (1 - C) + x * s * S - x ^ 2 * (c - 1) / 2 * (2 - (1 - x ^ 2) * (c - 1))) := by
    linear_combination
      (-(-4*C^2*x^2 - 2*C*c + 4*S^2*x^2 + 2*S*s*x - c^2*x^4 + c^2*x^2 - 2*c*x^2 + 2*c - s^2*x^4
        + s^2*x^2 + x^4 + x^2)) * hc + (-(4*c^2*x^2 - 4*c^2 - 4*x^2)) * hC
  rw [hG, hB, hA, hR, hK]
  rw [show ((1 - x ^ 2) * s ^ 2 + 2 - 2 * (c * C - x * s * S)) / (8 * s ^ 2) *
        ((x ^ 2 * s ^ 2 + 2 * c - 2 * (c * C - x * s * S)) / (8 * s ^ 2)) -
      ((c * S - x * s * C) / (4 * s ^ 2)) ^ 2
      = (((1 - x ^ 2) * s ^ 2 + 2 - 2 * (c * C - x * s * S))
        * (x ^ 2 * s ^ 2 + 2 * c - 2 * (c * C - x * s * S)) - 4 * (c * S - x * s * C) ^ 2)
        / (64 * s ^ 4) by field_simp; ring, key]
  field_simp
  ring

end Havelock.Tanh.LargeN
