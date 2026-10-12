/-
Havelock — Conjecture 13, odd `n`: fourth-order expansion of one pair term and of the one-body
term in canonical coordinates (paper `conj13/paper/main.tex`, proof of Thm odd).

For a pair of vortices with `w_i = w₀ + e_iA + f_iA²`, `w_j = w₀ + e_jA + f_jA²` (`w₀ = 1/(1+λ)`)
and angle difference `θ + A²d`, the pair quantity
`D = w_i + w_j − 2λw_iw_j − 2√(P(w_i)P(w_j)) cos(θ + A²d)`, `P(w) = w(1−λw)`, factors as
`D = 2√(P_iP_j)·(K + Y)` with `K = 1 − cos(θ + A²d)` and
`Y = (w_i + w_j − 2λw_iw_j)/(2√(P_iP_j)) − 1 = O(A²)`. Hence

* `jet_logPw`  `log P(w₀ + eA + fA²) = log w₀² + L₁e A + (L₁f + L₂e²)A² + (2L₂ef + L₃e³)A³
                 + (L₂f² + 3L₃e²f + L₄e⁴)A⁴ + O(A⁵)`, `Lₖ` the Taylor coefficients of
                 `log(1+(1+λ)x) + log(1−λ(1+λ)x)`;
* `jet_pair`   `log D − ½log P_i − ½log P_j = log 2 + log K₀ + ((K₂+Y₂)/K₀)A² + (Y₃/K₀)A³
                 + ((K₄+Y₄)/K₀ − (K₂+Y₂)²/(2K₀²))A⁴ + O(A⁵)`, with `K₀ = 1 − cos θ`,
                 `K₂ = d sin θ`, `K₄ = d² cos θ/2` and the explicit polynomials `Y₂, Y₃, Y₄`
                 (`Y2c`, `Y3c`, `Y4c`).
-/
import Havelock.Conj13.Odd4Jets
import Havelock.Conj13.OddQuadratic

open Real Asymptotics Filter Topology Finset

namespace Havelock.Conj13.Odd4

open Havelock.Conj13.Series Havelock.Conj13.Odd

namespace Jet4

variable {f : ℝ → ℝ} {a1 a2 a3 a4 : ℝ}

/-- `√(1 + x)` for a jet `x = O(A)`. -/
theorem sqrt_one (hf : Jet4 f 1 a1 a2 a3 a4) :
    Jet4 (fun A => √(f A)) 1 (a1 / 2) (a2 / 2 - a1 ^ 2 / 8) (a3 / 2 - a1 * a2 / 4 + a1 ^ 3 / 16)
      (a4 / 2 - a1 * a3 / 4 - a2 ^ 2 / 8 + 3 * a1 ^ 2 * a2 / 16 - 5 * a1 ^ 4 / 128) :=
  hf.sqrt one_pos (by ring) (by ring) (by ring) (by ring) (by ring)

/-- `1/(1 + x)` for a jet `x = O(A)`. -/
theorem inv_one (hf : Jet4 f 1 a1 a2 a3 a4) :
    Jet4 (fun A => (f A)⁻¹) 1 (-a1) (a1 ^ 2 - a2) (-a3 + 2 * a1 * a2 - a1 ^ 3)
      (-a4 + a2 ^ 2 + 2 * a1 * a3 - 3 * a1 ^ 2 * a2 + a1 ^ 4) :=
  hf.inv one_ne_zero (by ring) (by ring) (by ring) (by ring) (by ring)

end Jet4

/-! ### The one-body term -/

/-- `w₀(1 − λw₀) = w₀²`, i.e. `1 − λw₀ = w₀`. -/
theorem one_sub_b_w0 (b : ℝ) (hb : -1 < b) : 1 - b * w0 b = w0 b := by
  have h1 : (1 + b) ≠ 0 := by linarith
  unfold w0; field_simp; ring

theorem w0_pos (b : ℝ) (hb : -1 < b) : 0 < w0 b := by
  unfold w0; exact div_pos one_pos (by linarith)

/-- The Taylor coefficients of `log(1 + (1+λ)x) + log(1 − λ(1+λ)x)`. -/
noncomputable def L1 (b : ℝ) : ℝ := (1 + b) * (1 - b)
noncomputable def L2 (b : ℝ) : ℝ := -(1 + b) ^ 2 * (1 + b ^ 2) / 2
noncomputable def L3 (b : ℝ) : ℝ := (1 + b) ^ 3 * (1 - b ^ 3) / 3
noncomputable def L4 (b : ℝ) : ℝ := -(1 + b) ^ 4 * (1 + b ^ 4) / 4

/-- **One-body jet.** `log P(w₀ + eA + fA²)` to fourth order. -/
theorem jet_logPw (b e f : ℝ) (hb : -1 < b) :
    Jet4 (fun A => Real.log (Pw b (w0 b + A * e + A ^ 2 * f))) (Real.log (w0 b ^ 2)) (L1 b * e)
      (L1 b * f + L2 b * e ^ 2) (2 * L2 b * e * f + L3 b * e ^ 3)
      (L2 b * f ^ 2 + 3 * L3 b * e ^ 2 * f + L4 b * e ^ 4) := by
  have h1 : (1 + b) ≠ 0 := by linarith
  have hw := w0_pos b hb
  have hP : Jet4 (fun A => Pw b (w0 b + A * e + A ^ 2 * f)) (w0 b ^ 2) (w0 b * (1 - b) * e)
      (w0 b * (1 - b) * f - b * e ^ 2) (-2 * b * e * f) (-b * f ^ 2) := by
    refine (Jet4.poly _ _ _ _ _).congr (fun A => ?_)
    unfold Pw
    linear_combination (-(w0 b + A * e + A ^ 2 * f)) * one_sub_b_w0 b hb
  refine (hP.log (by positivity)).of_coeff rfl ?_ ?_ ?_ ?_ <;>
    · simp only [L1, L2, L3, L4, w0]; field_simp; try ring

/-! ### One pair -/

/-- The pair quantity `w_i + w_j − 2λw_iw_j − 2√(P(w_i)P(w_j)) cos(θ + A²d)` with
`w_i = w₀ + e_iA + f_iA²`, `w_j = w₀ + e_jA + f_jA²`. -/
noncomputable def Dg4 (b ei fi ej fj θ d A : ℝ) : ℝ :=
  (w0 b + A * ei + A ^ 2 * fi) + (w0 b + A * ej + A ^ 2 * fj)
    - 2 * b * (w0 b + A * ei + A ^ 2 * fi) * (w0 b + A * ej + A ^ 2 * fj)
    - 2 * √(Pw b (w0 b + A * ei + A ^ 2 * fi) * Pw b (w0 b + A * ej + A ^ 2 * fj))
      * Real.cos (θ + A ^ 2 * d)

/-- `Y₂` (coefficient of `A²` in `Y`). -/
noncomputable def Y2c (b ei ej : ℝ) : ℝ := (1 + b) ^ 4 * (ei - ej) ^ 2 / 8

/-- `Y₃`. -/
noncomputable def Y3c (b ei fi ej fj : ℝ) : ℝ :=
  (b - 1) * (1 + b) ^ 5 / 8 * ((ei - ej) ^ 2 * (ei + ej)) + (1 + b) ^ 4 / 4 * ((ei - ej) * (fi - fj))

/-- `Y₄`. -/
noncomputable def Y4c (b ei fi ej fj : ℝ) : ℝ :=
  3 * (1 + b) ^ 6 * (5 * b ^ 2 - 6 * b + 5) / 128 * ((ei - ej) ^ 2 * (ei ^ 2 + ej ^ 2))
    + (1 + b) ^ 6 * (9 * b ^ 2 - 14 * b + 9) / 64 * ((ei - ej) ^ 2 * (ei * ej))
    + (b - 1) * (1 + b) ^ 5 / 8 * ((ei - ej) ^ 2 * (fi + fj))
    + (b - 1) * (1 + b) ^ 5 / 4 * ((ei - ej) * (fi - fj) * (ei + ej))
    + (1 + b) ^ 4 / 8 * (fi - fj) ^ 2

/-- `g(A) = (1+u)(1−λu)`, `u = (1+λ)(eA + fA²)`, so that `P(w₀ + eA + fA²) = w₀² g(A)`. -/
noncomputable def gfac (b e f A : ℝ) : ℝ :=
  (1 + (1 + b) * (A * e + A ^ 2 * f)) * (1 - b * ((1 + b) * (A * e + A ^ 2 * f)))

theorem jet_gfac (b e f : ℝ) :
    Jet4 (gfac b e f) 1 ((1 - b) * (1 + b) * e) ((1 - b) * (1 + b) * f - b * (1 + b) ^ 2 * e ^ 2)
      (-2 * b * (1 + b) ^ 2 * e * f) (-b * (1 + b) ^ 2 * f ^ 2) :=
  (Jet4.poly _ _ _ _ _).congr (fun A => by unfold gfac; ring)

theorem Pw_eq_gfac (b e f A : ℝ) (hb : -1 < b) :
    Pw b (w0 b + A * e + A ^ 2 * f) = w0 b ^ 2 * gfac b e f A := by
  have h1 : (1 + b) ≠ 0 := by linarith
  unfold Pw gfac w0; field_simp; ring

/-- `N(A) = 2 + (1−λ)(u_i+u_j) − 2λu_iu_j`, so that `w_i + w_j − 2λw_iw_j = w₀² N(A)`. -/
noncomputable def Nfac (b ei fi ej fj A : ℝ) : ℝ :=
  2 + (1 - b) * ((1 + b) * (A * ei + A ^ 2 * fi) + (1 + b) * (A * ej + A ^ 2 * fj))
    - 2 * b * ((1 + b) * (A * ei + A ^ 2 * fi)) * ((1 + b) * (A * ej + A ^ 2 * fj))

theorem jet_Nfac (b ei fi ej fj : ℝ) :
    Jet4 (Nfac b ei fi ej fj) 2 ((1 - b) * (1 + b) * (ei + ej))
      ((1 - b) * (1 + b) * (fi + fj) - 2 * b * (1 + b) ^ 2 * ei * ej)
      (-2 * b * (1 + b) ^ 2 * (ei * fj + fi * ej)) (-2 * b * (1 + b) ^ 2 * fi * fj) :=
  (Jet4.poly _ _ _ _ _).congr (fun A => by unfold Nfac; ring)

theorem num_eq_Nfac (b ei fi ej fj A : ℝ) (hb : -1 < b) :
    (w0 b + A * ei + A ^ 2 * fi) + (w0 b + A * ej + A ^ 2 * fj)
      - 2 * b * (w0 b + A * ei + A ^ 2 * fi) * (w0 b + A * ej + A ^ 2 * fj)
      = w0 b ^ 2 * Nfac b ei fi ej fj A := by
  have h1 : (1 + b) ≠ 0 := by linarith
  unfold Nfac w0; field_simp; ring

/-- `Y = N/(2√(g_ig_j)) − 1`. -/
noncomputable def Yfun (b ei fi ej fj A : ℝ) : ℝ :=
  Nfac b ei fi ej fj A * (√(gfac b ei fi A * gfac b ej fj A))⁻¹ / 2 - 1

/-- **`Y = Y₂A² + Y₃A³ + Y₄A⁴ + O(A⁵)`.** -/
theorem jet_Yfun (b ei fi ej fj : ℝ) :
    Jet4 (Yfun b ei fi ej fj) 0 0 (Y2c b ei ej) (Y3c b ei fi ej fj) (Y4c b ei fi ej fj) := by
  have hG := ((jet_gfac b ei fi).mul (jet_gfac b ej fj)).of_coeff (b0 := 1) (by ring) rfl rfl rfl rfl
  have hR := (hG.sqrt_one).inv_one
  have hY := (((jet_Nfac b ei fi ej fj).mul hR).const_mul (1 / 2)).sub (Jet4.const 1)
  refine (hY.congr (fun A => by unfold Yfun; ring)).of_coeff ?_ ?_ ?_ ?_ ?_
  · ring
  · ring
  · unfold Y2c; ring
  · unfold Y3c; ring
  · unfold Y4c; ring

/-- The factorisation `D = 2w₀²√(g_ig_j)·(K + Y)`. -/
theorem Dg4_eq (b ei fi ej fj θ d A : ℝ) (hb : -1 < b) (hgi : 0 < gfac b ei fi A)
    (hgj : 0 < gfac b ej fj A) :
    Dg4 b ei fi ej fj θ d A = 2 * w0 b ^ 2 * √(gfac b ei fi A * gfac b ej fj A)
      * ((1 - Real.cos (θ + A ^ 2 * d)) + Yfun b ei fi ej fj A) := by
  have hG : 0 < √(gfac b ei fi A * gfac b ej fj A) := Real.sqrt_pos.mpr (mul_pos hgi hgj)
  have hs : √(Pw b (w0 b + A * ei + A ^ 2 * fi) * Pw b (w0 b + A * ej + A ^ 2 * fj))
      = w0 b ^ 2 * √(gfac b ei fi A * gfac b ej fj A) := by
    rw [Pw_eq_gfac b ei fi A hb, Pw_eq_gfac b ej fj A hb,
      show w0 b ^ 2 * gfac b ei fi A * (w0 b ^ 2 * gfac b ej fj A)
        = (w0 b ^ 2) ^ 2 * (gfac b ei fi A * gfac b ej fj A) by ring,
      Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq (sq_nonneg _)]
  unfold Dg4
  rw [num_eq_Nfac b ei fi ej fj A hb, hs]
  unfold Yfun
  field_simp
  ring

/-- **One pair, fourth order.** With `K₀ = 1 − cos θ > 0`:
`log D − ½log P(w_i) − ½log P(w_j) = log 2 + log K₀ + ((d sin θ + Y₂)/K₀)A² + (Y₃/K₀)A³
 + ((d² cos θ/2 + Y₄)/K₀ − (d sin θ + Y₂)²/(2K₀²))A⁴ + O(A⁵)`. -/
theorem jet_pair (b ei fi ej fj θ d : ℝ) (hb : -1 < b) (hθ : Real.cos θ < 1) :
    Jet4 (fun A => Real.log (Dg4 b ei fi ej fj θ d A)
        - Real.log (Pw b (w0 b + A * ei + A ^ 2 * fi)) / 2
        - Real.log (Pw b (w0 b + A * ej + A ^ 2 * fj)) / 2)
      (Real.log 2 + Real.log (1 - Real.cos θ)) 0
      ((Real.sin θ * d + Y2c b ei ej) / (1 - Real.cos θ))
      (Y3c b ei fi ej fj / (1 - Real.cos θ))
      ((Real.cos θ * d ^ 2 / 2 + Y4c b ei fi ej fj) / (1 - Real.cos θ)
        - (Real.sin θ * d + Y2c b ei ej) ^ 2 / (2 * (1 - Real.cos θ) ^ 2)) := by
  have hK0 : 0 < 1 - Real.cos θ := by linarith
  have hKY := (jet_one_sub_cos θ d).add (jet_Yfun b ei fi ej fj)
  have hK0' : 0 < 1 - Real.cos θ + 0 := by linarith
  have hL := (Jet4.const (Real.log 2)).add (hKY.log hK0')
  have hgi := (jet_gfac b ei fi).eventually_pos one_pos
  have hgj := (jet_gfac b ej fj).eventually_pos one_pos
  have hpos := hKY.eventually_pos hK0'
  have hw : 0 < w0 b := w0_pos b hb
  refine (hL.congr' ?_).of_coeff ?_ ?_ ?_ ?_ ?_
  · filter_upwards [hgi, hgj, hpos] with A h1 h2 h3
    have hG : 0 < √(gfac b ei fi A * gfac b ej fj A) := Real.sqrt_pos.mpr (mul_pos h1 h2)
    rw [Dg4_eq b ei fi ej fj θ d A hb h1 h2, Pw_eq_gfac b ei fi A hb, Pw_eq_gfac b ej fj A hb,
      Real.log_mul (by positivity) h3.ne', Real.log_mul (by positivity) hG.ne',
      Real.log_mul (by positivity) (by positivity), Real.log_sqrt (mul_pos h1 h2).le,
      Real.log_mul h1.ne' h2.ne', Real.log_mul (by positivity) h1.ne',
      Real.log_mul (by positivity) h2.ne']
    ring
  · simp
  · simp
  all_goals (have := hK0.ne'; field_simp; try ring)

end Havelock.Conj13.Odd4
