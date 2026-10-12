import Havelock.Tanh.LargeN.RowFourier

/-! # Tanh model, large `N`: the horocyclic row and its Bloch symbol (note Prop. `prop:row`)

Upper half-plane (curvature `−1`), tanh pair term `−ln|z − w| + ln|z − w̄|`, the row
`z_k = kℓ + i`, `k ∈ ℤ`, on the horocycle `Im z = 1`; `F = H − vJ` with `J = −∑ 1/y_j`.
With `p_k = 1/(kℓ − 2i)²` (`rowP`), the note's proof of Prop. `prop:row` gives the second derivatives of the
pair term between sites `0` and `k ≠ 0` at the row:

* on-site block (derivatives in `x₀, y₀`): `[[ (kℓ)⁻² − Re p_k, Im p_k ], [ Im p_k, −(kℓ)⁻² + Re p_k ]]`
  (`onXX`, `onXY`, `onYY`; the off-diagonal entries cancel in pairs `k ↔ −k`);
* coupling block (derivatives in `(x₀, y₀)` and `(x_k, y_k)`):
  `[[ −(kℓ)⁻² + Re p_k, Im p_k ], [ −Im p_k, (kℓ)⁻² + Re p_k ]]` (`cpXX`, `cpXY`, `cpYX`, `cpYY`);
* the momentum term contributes `diag(0, 2v)` on site, `v = ∑_{k≠0} 2/(k²ℓ² + 4)` (`rowVel`, from the
  `y`-force balance).

These entries are taken here as the DEFINITION of the per-pair Hessian (the derivative computation itself,
`Hess Re f = [[Re f'', −Im f''], [−Im f'', −Re f'']]` applied to `f = −log(z − z_k)` and `log(z − z̄_k)`,
is not formalised). The Bloch symbol on `(δx_k, δy_k) = e^{2πiθk}(X, Y)` is
`M(θ) = ∑_{k≠0} [on(k) + e^{2πiθk} cp(k)] + diag(0, 2v)` (`rowMxx`, `rowMxy`, `rowMyx`, `rowMyy`,
`ℤ`-indexed tsums whose `k = 0` term is set to `0`).

Main results (`row_bloch_symbol`): for `ℓ > 0`, `0 ≤ θ ≤ 1`, `τ = 2π/ℓ`, `x = 1 − 2θ`,
`v = (τ coth τ − 1)/2`, `(ℓ/2π)² M_xx = β(x,τ)`, `(ℓ/2π)² M_xy = i γ(x,τ)`, `(ℓ/2π)² M_yx = −i γ(x,τ)`,
`(ℓ/2π)² M_yy = E_∞(τ) − 1/(4 sinh²τ) − ψ(x,τ)`. Each is proved as a `HasSum` (so the series converge).
Also `hasSum_note_sum4`: sum (4) of the note in its literal form.
-/

noncomputable section
namespace Havelock.Tanh.LargeN

open Real
open Complex (I)

/-- `p_k = 1/(kℓ − 2i)²` (proof of note Prop. `prop:row`). -/
def rowP (ℓ : ℝ) (k : ℤ) : ℂ := 1 / ((k : ℂ) * ℓ - 2 * I) ^ 2

/-- `(kℓ)⁻²` (zero at `k = 0` by Lean's convention `1/0 = 0`). -/
def rowS (ℓ : ℝ) (k : ℤ) : ℝ := 1 / ((k : ℝ) * ℓ) ^ 2

/-- The Bloch phase `e^{iqk}`, `q = 2πθ`. -/
def blochPhase (θ : ℝ) (k : ℤ) : ℂ := Complex.exp (2 * π * I * k * θ)

/-- The closed form `P(θ) = 1/4 − τ² w(θ)`, `τ = 2π/ℓ`, of `∑_{k≠0} e^{iqk} p_k`. -/
def rowPsum (ℓ θ : ℝ) : ℝ := 1 / 4 - (2 * π / ℓ) ^ 2 * rowW (2 * π / ℓ) θ

/-- `p_{−k} = conj p_k`. -/
theorem rowP_neg (ℓ : ℝ) (k : ℤ) : rowP ℓ (-k) = (starRingEnd ℂ) (rowP ℓ k) := by
  simp only [rowP, map_div₀, map_one, map_pow, map_sub, map_mul, Complex.conj_ofReal,
    map_intCast, Complex.conj_I, map_ofNat]
  push_cast
  ring

/-- `e^{−iqk} = e^{i(2π−q)k}`. -/
theorem blochPhase_neg (θ : ℝ) (k : ℤ) : blochPhase θ (-k) = blochPhase (1 - θ) k := by
  simp only [blochPhase]
  have : 2 * π * I * (k : ℂ) * ((1 - θ : ℝ) : ℂ) =
      (k : ℂ) * (2 * π * I) + 2 * π * I * ((-k : ℤ) : ℂ) * θ := by push_cast; ring
  rw [this, Complex.exp_add, Complex.exp_int_mul_two_pi_mul_I, one_mul]

theorem blochPhase_zero (k : ℤ) : blochPhase 0 k = 1 := by simp [blochPhase]

/-- Sum (4) with the `k = 0` term removed: `∑_{k≠0} e^{iqk} p_k = P(θ)` for `0 ≤ θ ≤ 1`. -/
theorem hasSum_bloch_rowP {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (fun k : ℤ => if k = 0 then 0 else blochPhase θ k * rowP ℓ k) (rowPsum ℓ θ : ℂ) := by
  have h := hasSum_ite_zero (hasSum_rowP hℓ hθ)
  convert h using 1
  · funext k
    by_cases hk : k = 0
    · simp [hk]
    · simp only [hk, ite_false, blochPhase, rowP]; ring
  · simp only [rowPsum, Int.cast_zero, zero_mul, zero_sub, mul_zero,
      Complex.exp_zero, mul_one]
    push_cast
    rw [neg_sq, mul_pow, Complex.I_sq]
    ring

/-- The reflected sum `∑_{k≠0} e^{iqk} p_{−k} = P(1−θ)`. -/
theorem hasSum_bloch_rowP_neg {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (fun k : ℤ => if k = 0 then 0 else blochPhase θ k * rowP ℓ (-k))
      (rowPsum ℓ (1 - θ) : ℂ) := by
  have hθ' : 1 - θ ∈ Set.Icc (0 : ℝ) 1 := ⟨by linarith [hθ.2], by linarith [hθ.1]⟩
  rw [← (Equiv.neg ℤ).hasSum_iff]
  convert hasSum_bloch_rowP hℓ hθ' using 1
  funext k
  simp only [Function.comp_apply, Equiv.neg_apply, neg_eq_zero, neg_neg]
  by_cases hk : k = 0
  · simp [hk]
  · simp only [hk, ite_false, blochPhase_neg]

/-- `∑_{k≠0} e^{iqk} Re p_k = (P(θ) + P(1−θ))/2`. -/
theorem hasSum_bloch_re {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (fun k : ℤ => if k = 0 then 0 else blochPhase θ k * ((rowP ℓ k).re : ℂ))
      (((rowPsum ℓ θ + rowPsum ℓ (1 - θ)) / 2 : ℝ) : ℂ) := by
  have h := ((hasSum_bloch_rowP hℓ hθ).add (hasSum_bloch_rowP_neg hℓ hθ)).div_const 2
  convert h using 1
  · funext k
    by_cases hk : k = 0
    · simp [hk]
    · simp only [hk, ite_false, Complex.re_eq_add_conj, rowP_neg]; ring
  · push_cast; ring

/-- `∑_{k≠0} e^{iqk} Im p_k = (P(θ) − P(1−θ))/(2i)`. -/
theorem hasSum_bloch_im {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (fun k : ℤ => if k = 0 then 0 else blochPhase θ k * ((rowP ℓ k).im : ℂ))
      (((rowPsum ℓ θ - rowPsum ℓ (1 - θ)) : ℝ) / (2 * I)) := by
  have h := ((hasSum_bloch_rowP hℓ hθ).sub (hasSum_bloch_rowP_neg hℓ hθ)).div_const (2 * I)
  convert h using 1
  · funext k
    by_cases hk : k = 0
    · simp [hk]
    · simp only [hk, ite_false, Complex.im_eq_sub_conj, rowP_neg]; ring
  · push_cast; ring

/-- Sums (1)–(2) scaled: `∑_{k≠0} e^{iqk} (kℓ)⁻² = (τ²/2)(θ² − θ + 1/6)`, `τ = 2π/ℓ`. -/
theorem hasSum_bloch_rowS {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (fun k : ℤ => if k = 0 then 0 else blochPhase θ k * (rowS ℓ k : ℂ))
      (((2 * π / ℓ) ^ 2 / 2 * (θ ^ 2 - θ + 1 / 6) : ℝ) : ℂ) := by
  have h := (hasSum_inv_sq_fourier hθ).div_const ((ℓ : ℂ) ^ 2)
  have hℓ' : (ℓ : ℂ) ≠ 0 := by exact_mod_cast hℓ.ne'
  convert h using 1
  · funext k
    by_cases hk : k = 0
    · simp [hk]
    · simp only [hk, ite_false, blochPhase, rowS]; push_cast; field_simp
  · push_cast; field_simp

/-- Force balance / sum (3): `v = ∑_{k≠0} 2/(k²ℓ² + 4) = (τ coth τ − 1)/2`, `τ = 2π/ℓ`. -/
theorem hasSum_rowVel {ℓ : ℝ} (hℓ : 0 < ℓ) :
    HasSum (fun k : ℤ => if k = 0 then 0 else 2 / (((k : ℝ) * ℓ) ^ 2 + 4))
      (((2 * π / ℓ) * coth (2 * π / ℓ) - 1) / 2) := by
  set τ := 2 * π / ℓ with hτdef
  have hτ : 0 < τ := by positivity
  have h := hasSum_ite_zero ((hasSum_inv_sq_add_sq hτ).mul_left (τ ^ 2 / 2))
  have hℓτ : ℓ = 2 * π / τ := by rw [hτdef]; field_simp
  convert h using 1
  · funext k
    by_cases hk : k = 0
    · simp [hk]
    · simp only [hk, ite_false]
      rw [hℓτ]
      field_simp
      ring
  · simp only [Int.cast_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero,
      zero_add]
    field_simp

/-! ## Hyperbolic form of `w` -/

/-- `w(θ) = e^{xτ}(2θ sinh τ + e^{−τ})/(4 sinh²τ)`, `x = 1 − 2θ` (using `1 − e^{−2τ} = 2e^{−τ} sinh τ`). -/
theorem rowW_eq {τ : ℝ} (hτ : 0 < τ) (θ : ℝ) :
    rowW τ θ = Real.exp ((1 - 2 * θ) * τ) * (2 * θ * sinh τ + Real.exp (-τ)) / (4 * sinh τ ^ 2) := by
  have hs : sinh τ ≠ 0 := (sinh_pos_iff.mpr hτ).ne'
  have he : Real.exp (-τ) ≠ 0 := (Real.exp_pos _).ne'
  have h1 : 1 - Real.exp (-2 * τ) = 2 * Real.exp (-τ) * sinh τ := by
    rw [Real.sinh_eq, show -2 * τ = -τ + -τ by ring, Real.exp_add, Real.exp_neg]
    field_simp
  have h2 : Real.exp (-2 * τ) = Real.exp (-τ) ^ 2 := by
    rw [sq, ← Real.exp_add]; congr 1; ring
  have h3 : Real.exp (-2 * τ * θ) = Real.exp (-τ) * Real.exp ((1 - 2 * θ) * τ) := by
    rw [← Real.exp_add]; congr 1; ring
  rw [rowW, h1, h2, h3]
  field_simp
  ring

/-- `w(θ) + w(1−θ) = Q(x,τ)/(2 sinh²τ)`, `x = 1 − 2θ` (cf. note Thm 11, `W_m`). -/
theorem rowW_add {τ : ℝ} (hτ : 0 < τ) (θ : ℝ) :
    rowW τ θ + rowW τ (1 - θ) = Q (1 - 2 * θ) τ / (2 * sinh τ ^ 2) := by
  have hs : sinh τ ≠ 0 := (sinh_pos_iff.mpr hτ).ne'
  have hc : cosh τ = sinh τ + Real.exp (-τ) := by rw [Real.cosh_eq, Real.sinh_eq]; ring
  rw [rowW_eq hτ, rowW_eq hτ, show (1 - 2 * (1 - θ)) * τ = -((1 - 2 * θ) * τ) by ring, Q, hc,
    Real.cosh_eq ((1 - 2 * θ) * τ), Real.sinh_eq ((1 - 2 * θ) * τ)]
  have hy : Real.exp ((1 - 2 * θ) * τ) ≠ 0 := (Real.exp_pos _).ne'
  rw [Real.exp_neg ((1 - 2 * θ) * τ)]
  field_simp
  ring

/-- `w(θ) − w(1−θ) = R(x,τ)/(2 sinh²τ)`, `x = 1 − 2θ`. -/
theorem rowW_sub {τ : ℝ} (hτ : 0 < τ) (θ : ℝ) :
    rowW τ θ - rowW τ (1 - θ) = R (1 - 2 * θ) τ / (2 * sinh τ ^ 2) := by
  have hs : sinh τ ≠ 0 := (sinh_pos_iff.mpr hτ).ne'
  have hc : cosh τ = sinh τ + Real.exp (-τ) := by rw [Real.cosh_eq, Real.sinh_eq]; ring
  rw [rowW_eq hτ, rowW_eq hτ, show (1 - 2 * (1 - θ)) * τ = -((1 - 2 * θ) * τ) by ring, R, hc,
    Real.cosh_eq ((1 - 2 * θ) * τ), Real.sinh_eq ((1 - 2 * θ) * τ)]
  have hy : Real.exp ((1 - 2 * θ) * τ) ≠ 0 := (Real.exp_pos _).ne'
  rw [Real.exp_neg ((1 - 2 * θ) * τ)]
  field_simp
  ring

/-- `Q(1,τ) = 1`. -/
theorem Q_one (τ : ℝ) : Q 1 τ = 1 := by
  simp only [Q, one_mul]; nlinarith [cosh_sq_sub_sinh_sq τ]

/-! ## The Bloch symbol of the row -/

/-- On-site `xx` entry of the Hessian of the pair term `−ln|z₀ − z_k| + ln|z₀ − z̄_k|` with respect to
`z₀ = x₀ + i y₀`, at the row `z_j = jℓ + i`: `(kℓ)⁻² − Re p_k`. -/
def onXX (ℓ : ℝ) (k : ℤ) : ℝ := rowS ℓ k - (rowP ℓ k).re
/-- On-site `xy` (= `yx`) entry of the same Hessian: `Im p_k` (these cancel in pairs `k ↔ −k`). -/
def onXY (ℓ : ℝ) (k : ℤ) : ℝ := (rowP ℓ k).im
/-- On-site `yy` entry of the same Hessian: `−(kℓ)⁻² + Re p_k`. -/
def onYY (ℓ : ℝ) (k : ℤ) : ℝ := -rowS ℓ k + (rowP ℓ k).re
/-- Coupling `∂²/∂x₀∂x_k` of the pair term at the row: `−(kℓ)⁻² + Re p_k`. -/
def cpXX (ℓ : ℝ) (k : ℤ) : ℝ := -rowS ℓ k + (rowP ℓ k).re
/-- Coupling `∂²/∂x₀∂y_k` of the pair term at the row: `Im p_k`. -/
def cpXY (ℓ : ℝ) (k : ℤ) : ℝ := (rowP ℓ k).im
/-- Coupling `∂²/∂y₀∂x_k` of the pair term at the row: `−Im p_k`. -/
def cpYX (ℓ : ℝ) (k : ℤ) : ℝ := -(rowP ℓ k).im
/-- Coupling `∂²/∂y₀∂y_k` of the pair term at the row: `(kℓ)⁻² + Re p_k`. -/
def cpYY (ℓ : ℝ) (k : ℤ) : ℝ := rowS ℓ k + (rowP ℓ k).re

/-- The translation velocity of the row from the `y`-force balance: `v = ∑_{k≠0} 2/(k²ℓ² + 4)`
(the `y`-derivative of `ln|z₀ − z̄_k|` at the row is `2/(k²ℓ² + 4)`; that of `−ln|z₀ − z_k|` is `0`). -/
def rowVel (ℓ : ℝ) : ℝ := ∑' k : ℤ, if k = 0 then 0 else 2 / (((k : ℝ) * ℓ) ^ 2 + 4)

/-- Summand of `M_xx(θ)`: pair `k ≠ 0` contributes `on_xx(k) + e^{iqk} cp_xx(k)`, `q = 2πθ`. -/
def termXX (ℓ θ : ℝ) (k : ℤ) : ℂ :=
  if k = 0 then 0 else (onXX ℓ k : ℂ) + blochPhase θ k * (cpXX ℓ k : ℂ)
/-- Summand of `M_xy(θ)`. -/
def termXY (ℓ θ : ℝ) (k : ℤ) : ℂ :=
  if k = 0 then 0 else (onXY ℓ k : ℂ) + blochPhase θ k * (cpXY ℓ k : ℂ)
/-- Summand of `M_yx(θ)`. -/
def termYX (ℓ θ : ℝ) (k : ℤ) : ℂ :=
  if k = 0 then 0 else (onXY ℓ k : ℂ) + blochPhase θ k * (cpYX ℓ k : ℂ)
/-- Summand of the pair part of `M_yy(θ)`. -/
def termYY (ℓ θ : ℝ) (k : ℤ) : ℂ :=
  if k = 0 then 0 else (onYY ℓ k : ℂ) + blochPhase θ k * (cpYY ℓ k : ℂ)

/-- Bloch symbol `M_xx(θ) = ∑_{k≠0} [on_xx(k) + e^{2πiθk} cp_xx(k)]` of the Hessian of `F = H − vJ`
at the row `z_k = kℓ + i` on the Bloch wave `(δx_k, δy_k) = e^{2πiθk}(X, Y)`. -/
def rowMxx (ℓ θ : ℝ) : ℂ := ∑' k : ℤ, termXX ℓ θ k
/-- Bloch symbol `M_xy(θ) = ∑_{k≠0} [on_xy(k) + e^{2πiθk} cp_xy(k)]`. -/
def rowMxy (ℓ θ : ℝ) : ℂ := ∑' k : ℤ, termXY ℓ θ k
/-- Bloch symbol `M_yx(θ) = ∑_{k≠0} [on_xy(k) + e^{2πiθk} cp_yx(k)]`. -/
def rowMyx (ℓ θ : ℝ) : ℂ := ∑' k : ℤ, termYX ℓ θ k
/-- Bloch symbol `M_yy(θ) = ∑_{k≠0} [on_yy(k) + e^{2πiθk} cp_yy(k)] + 2v`; the last term is the
Hessian `∂²_y (v/y) = 2v` at `y = 1` of the momentum term `−vJ = v ∑_j 1/y_j`. -/
def rowMyy (ℓ θ : ℝ) : ℂ := (∑' k : ℤ, termYY ℓ θ k) + 2 * (rowVel ℓ : ℂ)

/-- `∑_{k≠0} (kℓ)⁻² = τ²/12` (sum (1)). -/
theorem hasSum_rowS0 {ℓ : ℝ} (hℓ : 0 < ℓ) :
    HasSum (fun k : ℤ => if k = 0 then 0 else (rowS ℓ k : ℂ))
      (((2 * π / ℓ) ^ 2 / 2 * ((0 : ℝ) ^ 2 - 0 + 1 / 6) : ℝ) : ℂ) := by
  convert hasSum_bloch_rowS hℓ (θ := 0) (by simp) using 1
  funext k; simp [blochPhase_zero]

/-- `∑_{k≠0} Re p_k = (P(0) + P(1))/2`. -/
theorem hasSum_re0 {ℓ : ℝ} (hℓ : 0 < ℓ) :
    HasSum (fun k : ℤ => if k = 0 then 0 else ((rowP ℓ k).re : ℂ))
      (((rowPsum ℓ 0 + rowPsum ℓ (1 - 0)) / 2 : ℝ) : ℂ) := by
  convert hasSum_bloch_re hℓ (θ := 0) (by simp) using 1
  funext k; simp [blochPhase_zero]

/-- `∑_{k≠0} Im p_k = (P(0) − P(1))/(2i)`. -/
theorem hasSum_im0 {ℓ : ℝ} (hℓ : 0 < ℓ) :
    HasSum (fun k : ℤ => if k = 0 then 0 else ((rowP ℓ k).im : ℂ))
      (((rowPsum ℓ 0 - rowPsum ℓ (1 - 0)) : ℝ) / (2 * I)) := by
  convert hasSum_bloch_im hℓ (θ := 0) (by simp) using 1
  funext k; simp [blochPhase_zero]

/-- note Prop. `prop:row`, `xx` entry: `∑_{k≠0}[on_xx + e^{iqk} cp_xx] = τ² β(x,τ)`. -/
theorem hasSum_termXX {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (termXX ℓ θ) (((2 * π / ℓ) ^ 2 * beta (1 - 2 * θ) (2 * π / ℓ) : ℝ) : ℂ) := by
  have h := (((hasSum_rowS0 hℓ).sub (hasSum_re0 hℓ)).sub (hasSum_bloch_rowS hℓ hθ)).add
    (hasSum_bloch_re hℓ hθ)
  convert h using 1
  · funext k
    by_cases hk : k = 0
    · simp [termXX, hk]
    · simp only [termXX, hk, ite_false, onXX, cpXX]; push_cast; ring
  · set τ := 2 * π / ℓ with hτdef
    have hτ : 0 < τ := by positivity
    have hA := rowW_add hτ θ
    have hA0 := rowW_add hτ 0
    have hQ1 : Q (1 - 2 * 0) τ = 1 := by rw [mul_zero, sub_zero]; exact Q_one τ
    have key : τ ^ 2 / 2 * ((0 : ℝ) ^ 2 - 0 + 1 / 6) - (rowPsum ℓ 0 + rowPsum ℓ (1 - 0)) / 2
        - τ ^ 2 / 2 * (θ ^ 2 - θ + 1 / 6) + (rowPsum ℓ θ + rowPsum ℓ (1 - θ)) / 2
        = τ ^ 2 * beta (1 - 2 * θ) τ := by
      simp only [rowPsum, ← hτdef, beta]
      linear_combination (-τ ^ 2 / 2) * hA + (τ ^ 2 / 2) * hA0 + (τ ^ 2 / (4 * sinh τ ^ 2)) * hQ1
    rw [← key]; push_cast; ring

/-- note Prop. `prop:row`, `yy` entry: `∑_{k≠0}[on_yy + e^{iqk} cp_yy] + 2v = τ² α(x,τ)`,
`α = E_∞ − 1/(4 sinh²τ) − ψ`. -/
theorem hasSum_termYY {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (fun k => termYY ℓ θ k + if k = 0 then 2 * (rowVel ℓ : ℂ) else 0)
      (((2 * π / ℓ) ^ 2 * (Einf (2 * π / ℓ) - 1 / (4 * sinh (2 * π / ℓ) ^ 2)
        - psi (1 - 2 * θ) (2 * π / ℓ)) : ℝ) : ℂ) := by
  have h := ((((hasSum_re0 hℓ).sub (hasSum_rowS0 hℓ)).add (hasSum_bloch_rowS hℓ hθ)).add
    (hasSum_bloch_re hℓ hθ)).add (hasSum_ite_eq (0 : ℤ) (2 * (rowVel ℓ : ℂ)))
  convert h using 1
  · funext k
    by_cases hk : k = 0
    · simp [termYY, hk]
    · simp only [termYY, hk, ite_false, onYY, cpYY]; push_cast; ring
  · set τ := 2 * π / ℓ with hτdef
    have hτ : 0 < τ := by positivity
    have hA := rowW_add hτ θ
    have hA0 := rowW_add hτ 0
    have hQ1 : Q (1 - 2 * 0) τ = 1 := by rw [mul_zero, sub_zero]; exact Q_one τ
    have hv : rowVel ℓ = (τ * coth τ - 1) / 2 := (hasSum_rowVel hℓ).tsum_eq
    have hE : τ ^ 2 * Einf τ = τ * coth τ - 1 / 2 := by rw [Einf]; field_simp
    have key : (rowPsum ℓ 0 + rowPsum ℓ (1 - 0)) / 2 - τ ^ 2 / 2 * ((0 : ℝ) ^ 2 - 0 + 1 / 6)
        + τ ^ 2 / 2 * (θ ^ 2 - θ + 1 / 6) + (rowPsum ℓ θ + rowPsum ℓ (1 - θ)) / 2
        + 2 * ((τ * coth τ - 1) / 2)
        = τ ^ 2 * (Einf τ - 1 / (4 * sinh τ ^ 2) - psi (1 - 2 * θ) τ) := by
      simp only [rowPsum, ← hτdef, psi]
      linear_combination (-τ ^ 2 / 2) * hA + (-τ ^ 2 / 2) * hA0
        + (-τ ^ 2 / (4 * sinh τ ^ 2)) * hQ1 + (-1) * hE
    rw [hv, ← key]; push_cast; ring

/-- `P(0) = P(1)`, so `∑_{k≠0} Im p_k = 0`. -/
theorem rowPsum_sub_zero {ℓ : ℝ} (hℓ : 0 < ℓ) : rowPsum ℓ 0 - rowPsum ℓ (1 - 0) = 0 := by
  have hτ : 0 < 2 * π / ℓ := by positivity
  simp only [rowPsum, sub_zero, rowW_zero_eq_one hτ]; ring

/-- `P(θ) − P(1−θ) = −2τ²γ(x,τ)`. -/
theorem rowPsum_sub {ℓ : ℝ} (hℓ : 0 < ℓ) (θ : ℝ) :
    rowPsum ℓ θ - rowPsum ℓ (1 - θ) = -2 * ((2 * π / ℓ) ^ 2 * gamma (1 - 2 * θ) (2 * π / ℓ)) := by
  set τ := 2 * π / ℓ with hτdef
  have hτ : 0 < τ := by positivity
  have hB := rowW_sub hτ θ
  simp only [rowPsum, ← hτdef, gamma]
  linear_combination (-τ ^ 2) * hB

/-- note Prop. `prop:row`, `xy` entry: `∑_{k≠0}[on_xy + e^{iqk} cp_xy] = i τ² γ(x,τ)`. -/
theorem hasSum_termXY {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (termXY ℓ θ) (I * (((2 * π / ℓ) ^ 2 * gamma (1 - 2 * θ) (2 * π / ℓ) : ℝ) : ℂ)) := by
  have h := (hasSum_im0 hℓ).add (hasSum_bloch_im hℓ hθ)
  convert h using 1
  · funext k
    by_cases hk : k = 0
    · simp [termXY, hk]
    · simp only [termXY, hk, ite_false, onXY, cpXY]
  · rw [rowPsum_sub_zero hℓ, rowPsum_sub hℓ, ← add_div, eq_div_iff (by simp)]
    push_cast
    linear_combination (2 * ((2 * (π : ℂ) / ℓ) ^ 2 * (gamma (1 - 2 * θ) (2 * π / ℓ) : ℂ)))
      * Complex.I_sq

/-- note Prop. `prop:row`, `yx` entry: `∑_{k≠0}[on_xy + e^{iqk} cp_yx] = −i τ² γ(x,τ)`. -/
theorem hasSum_termYX {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (termYX ℓ θ) (-I * (((2 * π / ℓ) ^ 2 * gamma (1 - 2 * θ) (2 * π / ℓ) : ℝ) : ℂ)) := by
  have h := (hasSum_im0 hℓ).sub (hasSum_bloch_im hℓ hθ)
  convert h using 1
  · funext k
    by_cases hk : k = 0
    · simp [termYX, hk]
    · simp only [termYX, hk, ite_false, onXY, cpYX]; push_cast; ring
  · rw [rowPsum_sub_zero hℓ, rowPsum_sub hℓ, ← sub_div, eq_div_iff (by simp)]
    push_cast
    linear_combination (-2 * ((2 * (π : ℂ) / ℓ) ^ 2 * (gamma (1 - 2 * θ) (2 * π / ℓ) : ℂ)))
      * Complex.I_sq

/-! ## Main statements (note Prop. `prop:row`) -/

/-- note Prop. `prop:row` (force balance): the row translates with
`v = ∑_{k≠0} 2/(k²ℓ² + 4) = (τ coth τ − 1)/2`, `τ = 2π/ℓ`. -/
theorem rowVel_eq {ℓ : ℝ} (hℓ : 0 < ℓ) :
    rowVel ℓ = ((2 * π / ℓ) * coth (2 * π / ℓ) - 1) / 2 :=
  (hasSum_rowVel hℓ).tsum_eq

/-- `(ℓ/2π)² τ² a = a`. -/
theorem scale_mul {ℓ : ℝ} (hℓ : 0 < ℓ) (a : ℝ) :
    ((ℓ / (2 * π)) ^ 2 : ℂ) * (((2 * π / ℓ) ^ 2 * a : ℝ) : ℂ) = a := by
  have hℓ' : (ℓ : ℂ) ≠ 0 := by exact_mod_cast hℓ.ne'
  have hπ : (π : ℂ) ≠ 0 := by exact_mod_cast pi_pos.ne'
  push_cast; field_simp

/-- note Prop. `prop:row`, `xx` entry: `(ℓ/2π)² M_xx(θ) = β(x,τ)`, `x = 1 − 2θ`, `τ = 2π/ℓ`. -/
theorem rowMxx_eq {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    ((ℓ / (2 * π)) ^ 2 : ℂ) * rowMxx ℓ θ = (beta (1 - 2 * θ) (2 * π / ℓ) : ℂ) := by
  rw [rowMxx, (hasSum_termXX hℓ hθ).tsum_eq, scale_mul hℓ]

/-- note Prop. `prop:row`, `xy` entry: `(ℓ/2π)² M_xy(θ) = i γ(x,τ)`. -/
theorem rowMxy_eq {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    ((ℓ / (2 * π)) ^ 2 : ℂ) * rowMxy ℓ θ = I * (gamma (1 - 2 * θ) (2 * π / ℓ) : ℂ) := by
  rw [rowMxy, (hasSum_termXY hℓ hθ).tsum_eq, mul_left_comm, scale_mul hℓ]

/-- note Prop. `prop:row`, `yx` entry: `(ℓ/2π)² M_yx(θ) = −i γ(x,τ)`. -/
theorem rowMyx_eq {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    ((ℓ / (2 * π)) ^ 2 : ℂ) * rowMyx ℓ θ = -I * (gamma (1 - 2 * θ) (2 * π / ℓ) : ℂ) := by
  rw [rowMyx, (hasSum_termYX hℓ hθ).tsum_eq, mul_left_comm, scale_mul hℓ]

/-- note Prop. `prop:row`, `yy` entry: `(ℓ/2π)² M_yy(θ) = α(x,τ) = E_∞(τ) − 1/(4 sinh²τ) − ψ(x,τ)`. -/
theorem rowMyy_eq {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    ((ℓ / (2 * π)) ^ 2 : ℂ) * rowMyy ℓ θ =
      ((Einf (2 * π / ℓ) - 1 / (4 * sinh (2 * π / ℓ) ^ 2) - psi (1 - 2 * θ) (2 * π / ℓ) : ℝ) : ℂ) := by
  have h := (hasSum_termYY hℓ hθ).sub (hasSum_ite_eq (0 : ℤ) (2 * (rowVel ℓ : ℂ)))
  simp only [add_sub_cancel_right] at h
  rw [rowMyy, h.tsum_eq, sub_add_cancel, scale_mul hℓ]

/-- note Prop. `prop:row`, combined: for `ℓ > 0`, `0 ≤ θ ≤ 1`, `τ = 2π/ℓ`, `x = 1 − 2θ`,
`(ℓ/2π)² M(θ) = [[β, iγ], [−iγ, α]](x,τ)` with `α = E_∞ − 1/(4 sinh²τ) − ψ`, and the row translates
with `v = (τ coth τ − 1)/2`. -/
theorem row_bloch_symbol {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    rowVel ℓ = ((2 * π / ℓ) * coth (2 * π / ℓ) - 1) / 2 ∧
    ((ℓ / (2 * π)) ^ 2 : ℂ) * rowMxx ℓ θ = (beta (1 - 2 * θ) (2 * π / ℓ) : ℂ) ∧
    ((ℓ / (2 * π)) ^ 2 : ℂ) * rowMxy ℓ θ = I * (gamma (1 - 2 * θ) (2 * π / ℓ) : ℂ) ∧
    ((ℓ / (2 * π)) ^ 2 : ℂ) * rowMyx ℓ θ = -I * (gamma (1 - 2 * θ) (2 * π / ℓ) : ℂ) ∧
    ((ℓ / (2 * π)) ^ 2 : ℂ) * rowMyy ℓ θ =
      ((Einf (2 * π / ℓ) - 1 / (4 * sinh (2 * π / ℓ) ^ 2) - psi (1 - 2 * θ) (2 * π / ℓ) : ℝ) : ℂ) :=
  ⟨rowVel_eq hℓ, rowMxx_eq hℓ hθ, rowMxy_eq hℓ hθ, rowMyx_eq hℓ hθ, rowMyy_eq hℓ hθ⟩

/-! ## Sum (4) in the note's form -/

/-- The note's right-hand side of sum (4): `∑_{n≥0}(q + 2πn) e^{−b(q+2πn)} = 2π w(θ)` with
`τ = πb`, `θ = q/(2π)` (geometric series). -/
theorem hasSum_note_rhs4 {b q : ℝ} (hb : 0 < b) :
    HasSum (fun n : ℕ => (q + 2 * π * n) * Real.exp (-b * (q + 2 * π * n)))
      (2 * π * rowW (π * b) (q / (2 * π))) := by
  set r := Real.exp (-2 * (π * b)) with hr
  have hr0 : 0 ≤ r := (Real.exp_pos _).le
  have hr1 : r < 1 := by
    have := Real.exp_lt_exp.mpr (show -2 * (π * b) < 0 by nlinarith [pi_pos])
    rwa [Real.exp_zero, ← hr] at this
  have hrn : ‖r‖ < 1 := by rw [Real.norm_eq_abs, abs_of_nonneg hr0]; exact hr1
  have h := (((hasSum_geometric_of_lt_one hr0 hr1).mul_left q).add
    ((hasSum_coe_mul_geometric_of_norm_lt_one hrn).mul_left (2 * π))).mul_left (Real.exp (-b * q))
  have h1r : 1 - r ≠ 0 := by linarith
  convert h using 1
  · funext n
    have : Real.exp (-b * (q + 2 * π * n)) = Real.exp (-b * q) * r ^ n := by
      rw [hr, ← Real.exp_nat_mul, ← Real.exp_add]; congr 1; ring
    rw [this]; ring
  · have hπ : π ≠ 0 := pi_pos.ne'
    have he : Real.exp (-2 * (π * b) * (q / (2 * π))) = Real.exp (-b * q) := by
      congr 1; field_simp
    rw [rowW, he, ← hr]
    field_simp

/-- Lattice sum (4) of the note, literally: for `b > 0` and `0 ≤ q ≤ 2π`,
`∑_{k∈ℤ} e^{iqk}/(k − ib)² = −2π ∑_{n≥0}(q + 2πn) e^{−b(q+2πn)}`. -/
theorem hasSum_note_sum4 {b q : ℝ} (hb : 0 < b) (hq0 : 0 ≤ q) (hq : q ≤ 2 * π) :
    HasSum (fun k : ℤ => Complex.exp (I * q * k) / ((k : ℂ) - I * b) ^ 2)
      ((-2 * π * ∑' n : ℕ, (q + 2 * π * n) * Real.exp (-b * (q + 2 * π * n)) : ℝ) : ℂ) := by
  set ℓ : ℝ := 2 / b with hℓdef
  have hℓ : 0 < ℓ := by positivity
  have hθ : q / (2 * π) ∈ Set.Icc (0 : ℝ) 1 :=
    ⟨by positivity, (div_le_one (by positivity)).mpr hq⟩
  have h := (hasSum_rowP hℓ hθ).mul_left ((ℓ : ℂ) ^ 2)
  have hτ : 2 * π / ℓ = π * b := by rw [hℓdef]; field_simp
  rw [(hasSum_note_rhs4 hb).tsum_eq]
  have hb' : (b : ℂ) ≠ 0 := by exact_mod_cast hb.ne'
  have hπ : (π : ℂ) ≠ 0 := by exact_mod_cast pi_pos.ne'
  convert h using 1
  · funext k
    have hk : (k : ℂ) - I * b ≠ 0 := by
      intro h0; have := congrArg Complex.im h0; simp at this; linarith
    have e1 : (k : ℂ) * ℓ - 2 * I = ℓ * ((k : ℂ) - I * b) := by
      rw [hℓdef]; push_cast; field_simp
    have e2 : 2 * π * I * k * ((q / (2 * π) : ℝ) : ℂ) = I * q * k := by
      push_cast; field_simp
    rw [e1, e2]
    have hℓ' : (ℓ : ℂ) ≠ 0 := by exact_mod_cast hℓ.ne'
    field_simp
  · rw [hτ, hℓdef]
    push_cast
    field_simp

end Havelock.Tanh.LargeN
