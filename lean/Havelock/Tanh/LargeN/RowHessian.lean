import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Havelock.Tanh.LargeN.RowSymbol

/-! # Tanh model, large `N`: derivation of the row Hessian (proof of note Prop. `prop:row`)

`RowSymbol.lean` takes the per-pair Hessian entries `onXX … cpYY` and the velocity `rowVel` of the horocyclic
row as definitions copied from the proof of note Prop. `prop:row`. Here they are derived by calculus.

* `pairG x₀ y₀ x₁ y₁ = −log|z₀ − z₁| + log|z₀ − z̄₁|`, `z = x + iy` (the tanh pair term in the upper half-plane).
* Second partial derivatives are taken as iterated one-variable `deriv`s of restrictions:
  `∂²g/∂s∂t` at a point `P` means `deriv (fun s ↦ deriv (fun t ↦ g(…s…t…)) t_P) s_P`, the inner
  derivative being taken first (in `t`) at every value of the outer variable `s`, then differentiated in `s`.
  The point is the row pair `(x₀, y₀, x₁, y₁) = (0, 1, kℓ, 1)`; the defs `hessOnXX … hessCpYY` (below) spell out
  exactly which restriction is used.
* `hessOnXX_eq … hessCpYY_eq`: these equal `onXX … cpYY` of `RowSymbol` (for `kℓ ≠ 0`).
* `pair_dy0`: `∂g/∂y₀ = 2/((kℓ)² + 4)`; `pair_dx0`: `∂g/∂x₀ = 1/(kℓ) − kℓ/((kℓ)² + 4)` (odd in `k`).
* Momentum `J = −∑ 1/y_j`, single-site term `rowJ y = −1/y`: `∂J/∂y = 1` and `∂²(−vJ)/∂y² = 2v` at `y = 1`.
* `row_force_balance_y`: `∑_{k≠0} ∂g/∂y₀ = rowVel ℓ · ∂J/∂y₀`, i.e. `∂(H − vJ)/∂y₀ = 0` at the row with
  `v = rowVel ℓ`; `row_force_balance_x`: `∑_{k≠0} ∂g/∂x₀ = 0` (and `J` does not depend on `x`).
* `row_bloch_from_hessian`: `rowMxx … rowMyy` are the `ℤ`-sums over `k ≠ 0` of
  (derived on-site second derivative) `+ e^{2πiθk}` (derived coupling second derivative), plus, for `M_yy`,
  the derived second derivative `∂²(−vJ)/∂y₀²` at `y₀ = 1` with `v = rowVel ℓ`.
-/

noncomputable section
namespace Havelock.Tanh.LargeN

open Real Filter Topology
open Complex (I)

/-- `z = x + i y`. -/
def rowPt (x y : ℝ) : ℂ := (x : ℂ) + (y : ℂ) * I

/-- The tanh pair term in the upper half-plane, `g = −log|z₀ − z₁| + log|z₀ − z̄₁|`, `z_j = x_j + i y_j`
(note, paragraph "The limit object is a vortex row on a horocycle"). -/
def pairG (x₀ y₀ x₁ y₁ : ℝ) : ℝ :=
  -Real.log ‖rowPt x₀ y₀ - rowPt x₁ y₁‖ + Real.log ‖rowPt x₀ y₀ - (starRingEnd ℂ) (rowPt x₁ y₁)‖

/-- `log(u² + v²)`. -/
def rowLq (u v : ℝ) : ℝ := Real.log (u ^ 2 + v ^ 2)

/-- `u/(u² + v²)`. -/
def rowPhi (u v : ℝ) : ℝ := u / (u ^ 2 + v ^ 2)

theorem log_norm_rowPt (u v : ℝ) : Real.log ‖rowPt u v‖ = rowLq u v / 2 := by
  rw [Complex.norm_def, Real.log_sqrt (Complex.normSq_nonneg _), rowLq]
  congr 2
  simp [rowPt, Complex.normSq_apply, sq]

/-- `g = −½ log((x₀−x₁)² + (y₀−y₁)²) + ½ log((x₀−x₁)² + (y₀+y₁)²)`. -/
theorem pairG_eq (x₀ y₀ x₁ y₁ : ℝ) :
    pairG x₀ y₀ x₁ y₁ = -(rowLq (x₀ - x₁) (y₀ - y₁)) / 2 + rowLq (x₀ - x₁) (y₀ + y₁) / 2 := by
  have h1 : rowPt x₀ y₀ - rowPt x₁ y₁ = rowPt (x₀ - x₁) (y₀ - y₁) := by
    simp only [rowPt]; push_cast; ring
  have h2 : rowPt x₀ y₀ - (starRingEnd ℂ) (rowPt x₁ y₁) = rowPt (x₀ - x₁) (y₀ + y₁) := by
    simp only [rowPt, map_add, map_mul, Complex.conj_ofReal, Complex.conj_I]; push_cast; ring
  rw [pairG, h1, h2, log_norm_rowPt, log_norm_rowPt]
  ring

/-! ## One-variable derivative helpers -/

theorem hasDerivAt_rowLq {U V : ℝ → ℝ} {U' V' t : ℝ} (hU : HasDerivAt U U' t)
    (hV : HasDerivAt V V' t) (h : U t ^ 2 + V t ^ 2 ≠ 0) :
    HasDerivAt (fun s => rowLq (U s) (V s)) (2 * (U t * U' + V t * V') / (U t ^ 2 + V t ^ 2)) t := by
  have := ((hU.fun_pow 2).add (hV.fun_pow 2)).log h
  convert this using 1
  · rfl
  · rw [show (2 : ℕ) - 1 = 1 from rfl]
    simp only [Pi.add_apply, pow_one, Nat.cast_ofNat]
    ring

theorem hasDerivAt_rowPhi {U V : ℝ → ℝ} {U' V' t : ℝ} (hU : HasDerivAt U U' t)
    (hV : HasDerivAt V V' t) (h : U t ^ 2 + V t ^ 2 ≠ 0) :
    HasDerivAt (fun s => rowPhi (U s) (V s))
      ((U' * (V t ^ 2 - U t ^ 2) - 2 * V' * U t * V t) / (U t ^ 2 + V t ^ 2) ^ 2) t := by
  have := hU.fun_div ((hU.fun_pow 2).add (hV.fun_pow 2)) h
  convert this using 1
  · rfl
  · rw [show (2 : ℕ) - 1 = 1 from rfl]
    simp only [Pi.add_apply, pow_one, Nat.cast_ofNat]
    ring

/-! ## First partial derivatives of `g` -/

/-- `∂g/∂x₀ = −u/(u²+s²) + u/(u²+t²)`, `u = x₀−x₁`, `s = y₀−y₁`, `t = y₀+y₁` (also `∂g/∂x₁ = −∂g/∂x₀`). -/
def pairGx (x₀ y₀ x₁ y₁ : ℝ) : ℝ := -rowPhi (x₀ - x₁) (y₀ - y₁) + rowPhi (x₀ - x₁) (y₀ + y₁)
/-- `∂g/∂y₀ = −s/(s²+u²) + t/(t²+u²)`. -/
def pairGy (x₀ y₀ x₁ y₁ : ℝ) : ℝ := -rowPhi (y₀ - y₁) (x₀ - x₁) + rowPhi (y₀ + y₁) (x₀ - x₁)
/-- `∂g/∂y₁ = s/(s²+u²) + t/(t²+u²)`. -/
def pairGy₁ (x₀ y₀ x₁ y₁ : ℝ) : ℝ := rowPhi (y₀ - y₁) (x₀ - x₁) + rowPhi (y₀ + y₁) (x₀ - x₁)

section first
variable {x₀ y₀ x₁ y₁ : ℝ}
  (h1 : (x₀ - x₁) ^ 2 + (y₀ - y₁) ^ 2 ≠ 0) (h2 : (x₀ - x₁) ^ 2 + (y₀ + y₁) ^ 2 ≠ 0)
include h1 h2

theorem hasDerivAt_pairG_x₀ : HasDerivAt (fun x => pairG x y₀ x₁ y₁) (pairGx x₀ y₀ x₁ y₁) x₀ := by
  have e : (fun x => pairG x y₀ x₁ y₁) =
      fun x => -(rowLq (x - x₁) (y₀ - y₁)) / 2 + rowLq (x - x₁) (y₀ + y₁) / 2 :=
    funext fun x => pairG_eq _ _ _ _
  rw [e]
  have hU : HasDerivAt (fun x : ℝ => x - x₁) 1 x₀ := (hasDerivAt_id' x₀).sub_const x₁
  have := ((hasDerivAt_rowLq hU (hasDerivAt_const x₀ (y₀ - y₁)) h1).neg.div_const 2).add
    ((hasDerivAt_rowLq hU (hasDerivAt_const x₀ (y₀ + y₁)) h2).div_const 2)
  convert this using 1
  simp only [pairGx, rowPhi]
  ring

theorem hasDerivAt_pairG_y₀ : HasDerivAt (fun y => pairG x₀ y x₁ y₁) (pairGy x₀ y₀ x₁ y₁) y₀ := by
  have e : (fun y => pairG x₀ y x₁ y₁) =
      fun y => -(rowLq (x₀ - x₁) (y - y₁)) / 2 + rowLq (x₀ - x₁) (y + y₁) / 2 :=
    funext fun y => pairG_eq _ _ _ _
  rw [e]
  have := ((hasDerivAt_rowLq (hasDerivAt_const y₀ (x₀ - x₁)) ((hasDerivAt_id' y₀).sub_const y₁)
    h1).neg.div_const 2).add
    ((hasDerivAt_rowLq (hasDerivAt_const y₀ (x₀ - x₁)) ((hasDerivAt_id' y₀).add_const y₁) h2).div_const 2)
  convert this using 1
  simp only [pairGy, rowPhi]
  ring

theorem hasDerivAt_pairG_x₁ : HasDerivAt (fun x => pairG x₀ y₀ x y₁) (-pairGx x₀ y₀ x₁ y₁) x₁ := by
  have e : (fun x => pairG x₀ y₀ x y₁) =
      fun x => -(rowLq (x₀ - x) (y₀ - y₁)) / 2 + rowLq (x₀ - x) (y₀ + y₁) / 2 :=
    funext fun x => pairG_eq _ _ _ _
  rw [e]
  have hU : HasDerivAt (fun x : ℝ => x₀ - x) (-1) x₁ := by
    simpa using (hasDerivAt_id' x₁).const_sub x₀
  have := ((hasDerivAt_rowLq hU (hasDerivAt_const x₁ (y₀ - y₁)) h1).neg.div_const 2).add
    ((hasDerivAt_rowLq hU (hasDerivAt_const x₁ (y₀ + y₁)) h2).div_const 2)
  convert this using 1
  simp only [pairGx, rowPhi]
  ring

theorem hasDerivAt_pairG_y₁ : HasDerivAt (fun y => pairG x₀ y₀ x₁ y) (pairGy₁ x₀ y₀ x₁ y₁) y₁ := by
  have e : (fun y => pairG x₀ y₀ x₁ y) =
      fun y => -(rowLq (x₀ - x₁) (y₀ - y)) / 2 + rowLq (x₀ - x₁) (y₀ + y) / 2 :=
    funext fun y => pairG_eq _ _ _ _
  rw [e]
  have hV : HasDerivAt (fun y : ℝ => y₀ - y) (-1) y₁ := by
    simpa using (hasDerivAt_id' y₁).const_sub y₀
  have := ((hasDerivAt_rowLq (hasDerivAt_const y₁ (x₀ - x₁)) hV h1).neg.div_const 2).add
    ((hasDerivAt_rowLq (hasDerivAt_const y₁ (x₀ - x₁)) ((hasDerivAt_id' y₁).const_add y₀)
      h2).div_const 2)
  convert this using 1
  simp only [pairGy₁, rowPhi]
  ring

end first

/-! ## The row pair `(0, 1, a, 1)`, `a = kℓ` -/

section row
variable {a : ℝ}

theorem rowD1_ne {x : ℝ} (hx : x ≠ a) : (x - a) ^ 2 + ((1 : ℝ) - 1) ^ 2 ≠ 0 := by
  have := sq_pos_of_ne_zero (sub_ne_zero.mpr hx); nlinarith
theorem rowD2_ne (x : ℝ) : (x - a) ^ 2 + ((1 : ℝ) + 1) ^ 2 ≠ 0 := by
  nlinarith [sq_nonneg (x - a)]
theorem rowD1_ne' (ha : a ≠ 0) (y : ℝ) : ((0 : ℝ) - a) ^ 2 + (y - 1) ^ 2 ≠ 0 := by
  have := sq_pos_of_ne_zero ha; nlinarith [sq_nonneg (y - 1)]
theorem rowD2_ne' (ha : a ≠ 0) (y : ℝ) : ((0 : ℝ) - a) ^ 2 + (y + 1) ^ 2 ≠ 0 := by
  have := sq_pos_of_ne_zero ha; nlinarith [sq_nonneg (y + 1)]

theorem rowD1_ne'' (ha : a ≠ 0) (y : ℝ) : ((0 : ℝ) - a) ^ 2 + ((1 : ℝ) - y) ^ 2 ≠ 0 := by
  have := sq_pos_of_ne_zero ha; nlinarith [sq_nonneg (1 - y)]
theorem rowD2_ne'' (ha : a ≠ 0) (y : ℝ) : ((0 : ℝ) - a) ^ 2 + ((1 : ℝ) + y) ^ 2 ≠ 0 := by
  have := sq_pos_of_ne_zero ha; nlinarith [sq_nonneg (1 + y)]

/-- Inner `x₀`-derivative near `x₀ = 0` (valid for `x₀ ≠ a`). -/
theorem ev_inner_x₀ (ha : a ≠ 0) :
    (fun x => deriv (fun x' => pairG x' 1 a 1) x) =ᶠ[𝓝 0] fun x => pairGx x 1 a 1 := by
  filter_upwards [eventually_ne_nhds (show (0 : ℝ) ≠ a from fun h => ha h.symm)] with x hx
  exact (hasDerivAt_pairG_x₀ (rowD1_ne hx) (rowD2_ne x)).deriv

end row

/-! ## Second partial derivatives at the row: definitions -/

/-- `∂²g/∂x₀²` at `(0, 1, a, 1)`. -/
def hessOnXX (a : ℝ) : ℝ := deriv (fun x => deriv (fun x' => pairG x' 1 a 1) x) 0
/-- `∂²g/∂x₀∂y₀` at `(0, 1, a, 1)` (inner `∂/∂y₀`, outer `∂/∂x₀`). -/
def hessOnXY (a : ℝ) : ℝ := deriv (fun x => deriv (fun y => pairG x y a 1) 1) 0
/-- `∂²g/∂y₀∂x₀` at `(0, 1, a, 1)` (inner `∂/∂x₀`, outer `∂/∂y₀`). -/
def hessOnYX (a : ℝ) : ℝ := deriv (fun y => deriv (fun x => pairG x y a 1) 0) 1
/-- `∂²g/∂y₀²` at `(0, 1, a, 1)`. -/
def hessOnYY (a : ℝ) : ℝ := deriv (fun y => deriv (fun y' => pairG 0 y' a 1) y) 1
/-- `∂²g/∂x₀∂x₁` at `(0, 1, a, 1)` (inner `∂/∂x₁`, outer `∂/∂x₀`). -/
def hessCpXX (a : ℝ) : ℝ := deriv (fun x => deriv (fun x' => pairG x 1 x' 1) a) 0
/-- `∂²g/∂x₀∂y₁` at `(0, 1, a, 1)` (inner `∂/∂y₁`, outer `∂/∂x₀`). -/
def hessCpXY (a : ℝ) : ℝ := deriv (fun x => deriv (fun y => pairG x 1 a y) 1) 0
/-- `∂²g/∂y₀∂x₁` at `(0, 1, a, 1)` (inner `∂/∂x₁`, outer `∂/∂y₀`). -/
def hessCpYX (a : ℝ) : ℝ := deriv (fun y => deriv (fun x' => pairG 0 y x' 1) a) 1
/-- `∂²g/∂y₀∂y₁` at `(0, 1, a, 1)` (inner `∂/∂y₁`, outer `∂/∂y₀`). -/
def hessCpYY (a : ℝ) : ℝ := deriv (fun y => deriv (fun y' => pairG 0 y a y') 1) 1

/-! ## Second partial derivatives at the row: values -/

section values
variable {a : ℝ} (ha : a ≠ 0)
include ha

theorem hessOnXX_val : hessOnXX a = 1 / a ^ 2 - (a ^ 2 - 4) / (a ^ 2 + 4) ^ 2 := by
  rw [hessOnXX, (ev_inner_x₀ ha).deriv_eq]
  have hU : HasDerivAt (fun x : ℝ => x - a) 1 0 := (hasDerivAt_id' 0).sub_const a
  have h := (hasDerivAt_rowPhi hU (hasDerivAt_const 0 ((1 : ℝ) - 1))
    (rowD1_ne (fun h => ha h.symm))).neg.add
    (hasDerivAt_rowPhi hU (hasDerivAt_const 0 ((1 : ℝ) + 1)) (rowD2_ne 0))
  rw [HasDerivAt.deriv (f := fun x => pairGx x 1 a 1) h]
  have : a ^ 2 + 4 ≠ 0 := by positivity
  field_simp
  ring

theorem hessOnXY_val : hessOnXY a = 4 * a / (a ^ 2 + 4) ^ 2 := by
  have hev : (fun x => deriv (fun y => pairG x y a 1) 1) =ᶠ[𝓝 0] fun x => pairGy x 1 a 1 := by
    filter_upwards [eventually_ne_nhds (show (0 : ℝ) ≠ a from fun h => ha h.symm)] with x hx
    exact (hasDerivAt_pairG_y₀ (rowD1_ne hx) (rowD2_ne x)).deriv
  rw [hessOnXY, hev.deriv_eq]
  have hU : HasDerivAt (fun x : ℝ => x - a) 1 0 := (hasDerivAt_id' 0).sub_const a
  have h := (hasDerivAt_rowPhi (hasDerivAt_const 0 ((1 : ℝ) - 1)) hU
    (by rw [add_comm]; exact rowD1_ne (fun h => ha h.symm))).neg.add
    (hasDerivAt_rowPhi (hasDerivAt_const 0 ((1 : ℝ) + 1)) hU (by rw [add_comm]; exact rowD2_ne 0))
  rw [HasDerivAt.deriv (f := fun x => pairGy x 1 a 1) h]
  have : a ^ 2 + 4 ≠ 0 := by positivity
  field_simp
  ring

theorem hessOnYX_val : hessOnYX a = 4 * a / (a ^ 2 + 4) ^ 2 := by
  have hev : (fun y => deriv (fun x => pairG x y a 1) 0) = fun y => pairGx 0 y a 1 :=
    funext fun y => (hasDerivAt_pairG_x₀ (rowD1_ne' ha y) (rowD2_ne' ha y)).deriv
  rw [hessOnYX, hev]
  have h := (hasDerivAt_rowPhi (hasDerivAt_const 1 ((0 : ℝ) - a)) ((hasDerivAt_id' 1).sub_const 1)
    (rowD1_ne' ha 1)).neg.add
    (hasDerivAt_rowPhi (hasDerivAt_const 1 ((0 : ℝ) - a)) ((hasDerivAt_id' 1).add_const 1)
      (rowD2_ne' ha 1))
  rw [HasDerivAt.deriv (f := fun y => pairGx 0 y a 1) h]
  have : a ^ 2 + 4 ≠ 0 := by positivity
  field_simp
  ring

theorem hessOnYY_val : hessOnYY a = -(1 / a ^ 2) + (a ^ 2 - 4) / (a ^ 2 + 4) ^ 2 := by
  have hev : (fun y => deriv (fun y' => pairG 0 y' a 1) y) = fun y => pairGy 0 y a 1 :=
    funext fun y => (hasDerivAt_pairG_y₀ (rowD1_ne' ha y) (rowD2_ne' ha y)).deriv
  rw [hessOnYY, hev]
  have h := (hasDerivAt_rowPhi ((hasDerivAt_id' 1).sub_const 1) (hasDerivAt_const 1 ((0 : ℝ) - a))
    (by rw [add_comm]; exact rowD1_ne' ha 1)).neg.add
    (hasDerivAt_rowPhi ((hasDerivAt_id' 1).add_const 1) (hasDerivAt_const 1 ((0 : ℝ) - a))
      (by rw [add_comm]; exact rowD2_ne' ha 1))
  rw [HasDerivAt.deriv (f := fun y => pairGy 0 y a 1) h]
  have : a ^ 2 + 4 ≠ 0 := by positivity
  field_simp
  ring

theorem hessCpXX_val : hessCpXX a = -(1 / a ^ 2) + (a ^ 2 - 4) / (a ^ 2 + 4) ^ 2 := by
  have hev : (fun x => deriv (fun x' => pairG x 1 x' 1) a) =ᶠ[𝓝 0] fun x => -pairGx x 1 a 1 := by
    filter_upwards [eventually_ne_nhds (show (0 : ℝ) ≠ a from fun h => ha h.symm)] with x hx
    exact (hasDerivAt_pairG_x₁ (rowD1_ne hx) (rowD2_ne x)).deriv
  rw [hessCpXX, hev.deriv_eq]
  have hU : HasDerivAt (fun x : ℝ => x - a) 1 0 := (hasDerivAt_id' 0).sub_const a
  have h := ((hasDerivAt_rowPhi hU (hasDerivAt_const 0 ((1 : ℝ) - 1))
    (rowD1_ne (fun h => ha h.symm))).neg.add
    (hasDerivAt_rowPhi hU (hasDerivAt_const 0 ((1 : ℝ) + 1)) (rowD2_ne 0))).neg
  rw [HasDerivAt.deriv (f := fun x => -pairGx x 1 a 1) h]
  have : a ^ 2 + 4 ≠ 0 := by positivity
  field_simp
  ring

theorem hessCpXY_val : hessCpXY a = 4 * a / (a ^ 2 + 4) ^ 2 := by
  have hev : (fun x => deriv (fun y => pairG x 1 a y) 1) =ᶠ[𝓝 0] fun x => pairGy₁ x 1 a 1 := by
    filter_upwards [eventually_ne_nhds (show (0 : ℝ) ≠ a from fun h => ha h.symm)] with x hx
    exact (hasDerivAt_pairG_y₁ (rowD1_ne hx) (rowD2_ne x)).deriv
  rw [hessCpXY, hev.deriv_eq]
  have hU : HasDerivAt (fun x : ℝ => x - a) 1 0 := (hasDerivAt_id' 0).sub_const a
  have h := (hasDerivAt_rowPhi (hasDerivAt_const 0 ((1 : ℝ) - 1)) hU
    (by rw [add_comm]; exact rowD1_ne (fun h => ha h.symm))).add
    (hasDerivAt_rowPhi (hasDerivAt_const 0 ((1 : ℝ) + 1)) hU (by rw [add_comm]; exact rowD2_ne 0))
  rw [HasDerivAt.deriv (f := fun x => pairGy₁ x 1 a 1) h]
  have : a ^ 2 + 4 ≠ 0 := by positivity
  field_simp
  ring

theorem hessCpYX_val : hessCpYX a = -(4 * a / (a ^ 2 + 4) ^ 2) := by
  have hev : (fun y => deriv (fun x' => pairG 0 y x' 1) a) = fun y => -pairGx 0 y a 1 :=
    funext fun y => (hasDerivAt_pairG_x₁ (rowD1_ne' ha y) (rowD2_ne' ha y)).deriv
  rw [hessCpYX, hev]
  have h := ((hasDerivAt_rowPhi (hasDerivAt_const 1 ((0 : ℝ) - a)) ((hasDerivAt_id' 1).sub_const 1)
    (rowD1_ne' ha 1)).neg.add
    (hasDerivAt_rowPhi (hasDerivAt_const 1 ((0 : ℝ) - a)) ((hasDerivAt_id' 1).add_const 1)
      (rowD2_ne' ha 1))).neg
  rw [HasDerivAt.deriv (f := fun y => -pairGx 0 y a 1) h]
  have : a ^ 2 + 4 ≠ 0 := by positivity
  field_simp
  ring

theorem hessCpYY_val : hessCpYY a = 1 / a ^ 2 + (a ^ 2 - 4) / (a ^ 2 + 4) ^ 2 := by
  have hev : (fun y => deriv (fun y' => pairG 0 y a y') 1) = fun y => pairGy₁ 0 y a 1 :=
    funext fun y => (hasDerivAt_pairG_y₁ (rowD1_ne' ha y) (rowD2_ne' ha y)).deriv
  rw [hessCpYY, hev]
  have h := (hasDerivAt_rowPhi ((hasDerivAt_id' 1).sub_const 1) (hasDerivAt_const 1 ((0 : ℝ) - a))
    (by rw [add_comm]; exact rowD1_ne' ha 1)).add
    (hasDerivAt_rowPhi ((hasDerivAt_id' 1).add_const 1) (hasDerivAt_const 1 ((0 : ℝ) - a))
      (by rw [add_comm]; exact rowD2_ne' ha 1))
  rw [HasDerivAt.deriv (f := fun y => pairGy₁ 0 y a 1) h]
  have : a ^ 2 + 4 ≠ 0 := by positivity
  field_simp
  ring

/-- First derivative: `∂g/∂y₀ = 2/(a² + 4)` at `(0, 1, a, 1)` (the `−log|z₀ − z₁|` part contributes `0`). -/
theorem pair_dy0 : deriv (fun y => pairG 0 y a 1) 1 = 2 / (a ^ 2 + 4) := by
  rw [(hasDerivAt_pairG_y₀ (rowD1_ne' ha 1) (rowD2_ne' ha 1)).deriv]
  simp only [pairGy, rowPhi]
  have : a ^ 2 + 4 ≠ 0 := by positivity
  field_simp
  ring

/-- First derivative: `∂g/∂x₀ = 1/a − a/(a² + 4)` at `(0, 1, a, 1)` (odd in `a`). -/
theorem pair_dx0 : deriv (fun x => pairG x 1 a 1) 0 = 1 / a - a / (a ^ 2 + 4) := by
  rw [(hasDerivAt_pairG_x₀ (rowD1_ne (fun h => ha h.symm)) (rowD2_ne 0)).deriv]
  simp only [pairGx, rowPhi]
  have : a ^ 2 + 4 ≠ 0 := by positivity
  field_simp
  ring

end values

/-! ## `Re p_k`, `Im p_k` -/

theorem rowP_eq (ℓ : ℝ) (k : ℤ) :
    rowP ℓ k = ((((k : ℝ) * ℓ) ^ 2 - 4 : ℝ) + ((4 * ((k : ℝ) * ℓ)) : ℝ) * I) /
      (((((k : ℝ) * ℓ) ^ 2 + 4) ^ 2 : ℝ) : ℂ) := by
  set a : ℝ := (k : ℝ) * ℓ with hadef
  have hD : ((((a ^ 2 + 4) ^ 2 : ℝ)) : ℂ) ≠ 0 := by
    have : (a ^ 2 + 4) ^ 2 ≠ 0 := by positivity
    exact_mod_cast this
  have hw : (k : ℂ) * ℓ - 2 * I = (a : ℂ) - 2 * I := by rw [hadef]; push_cast; ring
  have hw0 : (a : ℂ) - 2 * I ≠ 0 := by
    intro h; have := congrArg Complex.im h; simp at this
  rw [rowP, hw, div_eq_div_iff (pow_ne_zero 2 hw0) hD]
  push_cast
  linear_combination (-(4 * ((a : ℂ) ^ 2 - 4 + 4 * a * I) - 16 * (a : ℂ) ^ 2)) * Complex.I_sq

theorem rowP_re (ℓ : ℝ) (k : ℤ) :
    (rowP ℓ k).re = (((k : ℝ) * ℓ) ^ 2 - 4) / (((k : ℝ) * ℓ) ^ 2 + 4) ^ 2 := by
  rw [rowP_eq, Complex.div_ofReal_re, Complex.add_re, Complex.ofReal_re, Complex.re_ofReal_mul,
    Complex.I_re, mul_zero, add_zero]

theorem rowP_im (ℓ : ℝ) (k : ℤ) :
    (rowP ℓ k).im = 4 * ((k : ℝ) * ℓ) / (((k : ℝ) * ℓ) ^ 2 + 4) ^ 2 := by
  rw [rowP_eq, Complex.div_ofReal_im, Complex.add_im, Complex.ofReal_im, Complex.im_ofReal_mul,
    Complex.I_im, mul_one, zero_add]

/-! ## The derived entries are those of `RowSymbol` -/

section matchRow
variable {ℓ : ℝ} {k : ℤ} (hkℓ : (k : ℝ) * ℓ ≠ 0)
include hkℓ

theorem hessOnXX_eq : hessOnXX (k * ℓ) = onXX ℓ k := by
  rw [hessOnXX_val hkℓ, onXX, rowS, rowP_re]
theorem hessOnXY_eq : hessOnXY (k * ℓ) = onXY ℓ k := by
  rw [hessOnXY_val hkℓ, onXY, rowP_im]
theorem hessOnYX_eq : hessOnYX (k * ℓ) = onXY ℓ k := by
  rw [hessOnYX_val hkℓ, onXY, rowP_im]
theorem hessOnYY_eq : hessOnYY (k * ℓ) = onYY ℓ k := by
  rw [hessOnYY_val hkℓ, onYY, rowS, rowP_re]
theorem hessCpXX_eq : hessCpXX (k * ℓ) = cpXX ℓ k := by
  rw [hessCpXX_val hkℓ, cpXX, rowS, rowP_re]
theorem hessCpXY_eq : hessCpXY (k * ℓ) = cpXY ℓ k := by
  rw [hessCpXY_val hkℓ, cpXY, rowP_im]
theorem hessCpYX_eq : hessCpYX (k * ℓ) = cpYX ℓ k := by
  rw [hessCpYX_val hkℓ, cpYX, rowP_im]
theorem hessCpYY_eq : hessCpYY (k * ℓ) = cpYY ℓ k := by
  rw [hessCpYY_val hkℓ, cpYY, rowS, rowP_re]

end matchRow

theorem kℓ_ne {ℓ : ℝ} (hℓ : ℓ ≠ 0) {k : ℤ} (hk : k ≠ 0) : (k : ℝ) * ℓ ≠ 0 :=
  mul_ne_zero (Int.cast_ne_zero.mpr hk) hℓ

/-- note Prop. `prop:row`, proof: the per-pair on-site and coupling Hessian entries used in `RowSymbol` are
the second partial derivatives of the pair term `g` at the row pair `(z₀, z_k) = (i, kℓ + i)`, `kℓ ≠ 0`.
(`∂²g/∂x₀∂y₀ = ∂²g/∂y₀∂x₀ = onXY`.) -/
theorem pair_hessian_row {ℓ : ℝ} (hℓ : ℓ ≠ 0) {k : ℤ} (hk : k ≠ 0) :
    hessOnXX (k * ℓ) = onXX ℓ k ∧ hessOnXY (k * ℓ) = onXY ℓ k ∧ hessOnYX (k * ℓ) = onXY ℓ k ∧
    hessOnYY (k * ℓ) = onYY ℓ k ∧ hessCpXX (k * ℓ) = cpXX ℓ k ∧ hessCpXY (k * ℓ) = cpXY ℓ k ∧
    hessCpYX (k * ℓ) = cpYX ℓ k ∧ hessCpYY (k * ℓ) = cpYY ℓ k :=
  have h := kℓ_ne hℓ hk
  ⟨hessOnXX_eq h, hessOnXY_eq h, hessOnYX_eq h, hessOnYY_eq h, hessCpXX_eq h, hessCpXY_eq h,
    hessCpYX_eq h, hessCpYY_eq h⟩

/-! ## Momentum term and force balance -/

/-- Single-site momentum term `J(y) = −1/y` (`J = −∑_j 1/y_j`). -/
def rowJ (y : ℝ) : ℝ := -1 / y

theorem hasDerivAt_rowJ {y : ℝ} (hy : y ≠ 0) : HasDerivAt rowJ (1 / y ^ 2) y := by
  have := ((hasDerivAt_id' y).inv hy).neg
  convert this using 1
  · funext s; simp [rowJ, div_eq_mul_inv]
  · field_simp

/-- `∂J/∂y = 1/y² = 1` at `y = 1`. -/
theorem deriv_rowJ_one : deriv rowJ 1 = 1 := by
  rw [(hasDerivAt_rowJ one_ne_zero).deriv]; norm_num

/-- The momentum part of `F = H − vJ`: `∂²(−vJ)/∂y² = 2v` at `y = 1`. -/
theorem momentum_hess (v : ℝ) : deriv (fun y => deriv (fun y' => -v * rowJ y') y) 1 = 2 * v := by
  have hev : (fun y => deriv (fun y' => -v * rowJ y') y) =ᶠ[𝓝 1] fun y => -v * y⁻¹ ^ 2 := by
    filter_upwards [eventually_ne_nhds (one_ne_zero : (1 : ℝ) ≠ 0)] with y hy
    rw [((hasDerivAt_rowJ hy).const_mul (-v)).deriv]; field_simp
  rw [hev.deriv_eq]
  have h := (((hasDerivAt_id' (1 : ℝ)).inv one_ne_zero).fun_pow 2).const_mul (-v)
  rw [HasDerivAt.deriv (f := fun y => -v * y⁻¹ ^ 2) h]
  norm_num
  ring

/-- note Prop. `prop:row`, force balance in `y`: `∑_{k≠0} ∂g/∂y₀ (0,1,kℓ,1) = v · ∂J/∂y₀(1)` with
`v = rowVel ℓ`, i.e. `∂(H − vJ)/∂y₀ = 0` at the row. -/
theorem row_force_balance_y {ℓ : ℝ} (hℓ : 0 < ℓ) :
    HasSum (fun k : ℤ => if k = 0 then 0 else deriv (fun y => pairG 0 y (k * ℓ) 1) 1)
      (rowVel ℓ * deriv rowJ 1) := by
  rw [deriv_rowJ_one, mul_one, rowVel, (hasSum_rowVel hℓ).tsum_eq]
  convert hasSum_rowVel hℓ using 2 with k
  by_cases hk : k = 0
  · simp [hk]
  · simp only [hk, ite_false, pair_dy0 (kℓ_ne hℓ.ne' hk)]

/-- `rowVel` is the velocity fixed by the `y`-force balance: `rowVel ℓ = ∑_{k≠0} ∂g/∂y₀ (0,1,kℓ,1)`. -/
theorem rowVel_eq_force {ℓ : ℝ} (hℓ : 0 < ℓ) :
    rowVel ℓ = ∑' k : ℤ, if k = 0 then 0 else deriv (fun y => pairG 0 y (k * ℓ) 1) 1 := by
  rw [(row_force_balance_y hℓ).tsum_eq, deriv_rowJ_one, mul_one]

/-- Force balance in `x`: `∑_{k≠0} ∂g/∂x₀ (0,1,kℓ,1) = 0` (the summand is odd in `k`; `J` is
independent of `x`). -/
theorem row_force_balance_x {ℓ : ℝ} (hℓ : ℓ ≠ 0) :
    ∑' k : ℤ, (if k = 0 then 0 else deriv (fun x => pairG x 1 (k * ℓ) 1) 0) = 0 := by
  set f : ℤ → ℝ := fun k => if k = 0 then 0 else deriv (fun x => pairG x 1 (k * ℓ) 1) 0 with hf
  have hodd : ∀ k, f (-k) = -f k := by
    intro k
    by_cases hk : k = 0
    · simp [hf, hk]
    · simp only [hf, neg_eq_zero, hk, ite_false, pair_dx0 (kℓ_ne hℓ hk),
        pair_dx0 (kℓ_ne hℓ (neg_ne_zero.mpr hk))]
      push_cast; ring
  have h : ∑' k, f k = -∑' k, f k := by
    rw [← tsum_neg, ← (Equiv.neg ℤ).tsum_eq]
    exact tsum_congr fun k => by simp [hodd]
  linarith

/-! ## The Bloch symbol from the derived Hessian -/

/-- note Prop. `prop:row`: the Bloch entries of `RowSymbol` are the `ℤ`-sums over `k ≠ 0` of the derived
second derivatives of `F = H − vJ` at the row: on-site `∂²/∂(x₀,y₀)²` of the pair term `(0,k)` plus
`e^{2πiθk}` times the coupling `∂²/∂(x₀,y₀)∂(x_k,y_k)`; the momentum term `−vJ` (which depends only on the `y_j`)
adds its derived on-site `∂²/∂y₀² = 2v` to `M_yy`, with `v = rowVel ℓ` the force-balance velocity
(`rowVel_eq_force`). -/
theorem row_bloch_from_hessian {ℓ : ℝ} (hℓ : ℓ ≠ 0) (θ : ℝ) :
    rowMxx ℓ θ = ∑' k : ℤ, (if k = 0 then 0 else
      (hessOnXX (k * ℓ) : ℂ) + blochPhase θ k * (hessCpXX (k * ℓ) : ℂ)) ∧
    rowMxy ℓ θ = ∑' k : ℤ, (if k = 0 then 0 else
      (hessOnXY (k * ℓ) : ℂ) + blochPhase θ k * (hessCpXY (k * ℓ) : ℂ)) ∧
    rowMyx ℓ θ = ∑' k : ℤ, (if k = 0 then 0 else
      (hessOnYX (k * ℓ) : ℂ) + blochPhase θ k * (hessCpYX (k * ℓ) : ℂ)) ∧
    rowMyy ℓ θ = (∑' k : ℤ, (if k = 0 then 0 else
      (hessOnYY (k * ℓ) : ℂ) + blochPhase θ k * (hessCpYY (k * ℓ) : ℂ))) +
      ((deriv (fun y => deriv (fun y' => -rowVel ℓ * rowJ y') y) 1 : ℝ) : ℂ) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [rowMxx]; congr 1; funext k
    by_cases hk : k = 0
    · simp [termXX, hk]
    · obtain ⟨h1, -, -, -, h5, -, -, -⟩ := pair_hessian_row hℓ hk
      simp only [termXX, hk, ite_false, h1, h5]
  · rw [rowMxy]; congr 1; funext k
    by_cases hk : k = 0
    · simp [termXY, hk]
    · obtain ⟨-, h2, -, -, -, h6, -, -⟩ := pair_hessian_row hℓ hk
      simp only [termXY, hk, ite_false, h2, h6]
  · rw [rowMyx]; congr 1; funext k
    by_cases hk : k = 0
    · simp [termYX, hk]
    · obtain ⟨-, -, h3, -, -, -, h7, -⟩ := pair_hessian_row hℓ hk
      simp only [termYX, hk, ite_false, h3, h7]
  · rw [rowMyy, momentum_hess]; push_cast; congr 1; congr 1; funext k
    by_cases hk : k = 0
    · simp [termYY, hk]
    · obtain ⟨-, -, -, h4, -, -, -, h8⟩ := pair_hessian_row hℓ hk
      simp only [termYY, hk, ite_false, h4, h8]

end Havelock.Tanh.LargeN
