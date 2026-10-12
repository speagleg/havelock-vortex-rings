import Havelock.Tanh.General
import Havelock.Tanh.LargeN.Defs

/-! # Tanh model, large `N`: the exact symbol (note Thm 11)

With `ξ = exp(−2τ/N)` and `x = 1 − 2m/N`, the blocks `b_m`, `a_m`, `g_m` of note Thm 2 are `N²` times
`β(x,τ)`, `E_N(τ) − 1/(4 sinh²τ) − ψ(x,τ)` and `γ(x,τ)`. Proof: `ξ^m = e^{xτ} e^{−τ}`,
`ξ^{N−m} = e^{−xτ} e^{−τ}`, `ξ^N = e^{−2τ}`; every identity is then rational in `e^τ`, `e^{xτ}`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real

/-- `0 < xiOf N τ` (note Thm 11). -/
theorem xiOf_pos (N : ℕ) (τ : ℝ) : 0 < xiOf N τ := exp_pos _

/-- `xiOf N τ < 1` for `N ≥ 1`, `τ > 0` (note Thm 11). -/
theorem xiOf_lt_one {N : ℕ} {τ : ℝ} (hN : 1 ≤ N) (hτ : 0 < τ) : xiOf N τ < 1 := by
  unfold xiOf
  rw [exp_lt_one_iff]
  have : (0 : ℝ) < N := by exact_mod_cast hN
  have : 0 < 2 * τ / N := by positivity
  linarith [show -2 * τ / N = -(2 * τ / N) by ring]

/-- `c = cosh r = cR ξ = coth(τ/N)` (note §`sec:asym`). -/
theorem cR_xiOf {N : ℕ} {τ : ℝ} (hN : 1 ≤ N) (hτ : 0 < τ) : cR (xiOf N τ) = coth (τ / N) := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  set z := exp (τ / N) with hz
  have hz1 : 1 < z := one_lt_exp_iff.mpr (by positivity)
  have hξ : xiOf N τ = (z ^ 2)⁻¹ := by
    unfold xiOf
    rw [hz, ← exp_nat_mul, ← exp_neg]
    congr 1
    push_cast
    ring
  have hz0 : z ≠ 0 := by positivity
  have hz2 : z ^ 2 - 1 ≠ 0 := by nlinarith
  have hz2' : z ^ 2 + 1 ≠ 0 := by positivity
  unfold coth cR
  rw [hξ, cosh_eq, sinh_eq, exp_neg, ← hz]
  field_simp

/-- Powers of `ξ = xiOf N τ` in terms of `E = e^τ`, `y = e^{xτ}` (proof of note Thm 11). -/
theorem xiOf_pow_m {N m : ℕ} (hN : 1 ≤ N) (τ : ℝ) :
    xiOf N τ ^ m = exp (xOf N m * τ) * (exp τ)⁻¹ := by
  have hNpos : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  unfold xiOf xOf
  rw [← exp_nat_mul, ← exp_neg, ← exp_add]
  congr 1
  field_simp
  ring

theorem xiOf_pow_sub {N m : ℕ} (hN : 1 ≤ N) (hm : m ≤ N) (τ : ℝ) :
    xiOf N τ ^ (N - m) = (exp (xOf N m * τ))⁻¹ * (exp τ)⁻¹ := by
  have hNpos : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  unfold xiOf xOf
  rw [← exp_nat_mul, ← exp_neg, ← exp_neg, ← exp_add, Nat.cast_sub hm]
  congr 1
  field_simp
  ring

theorem xiOf_pow_N {N : ℕ} (hN : 1 ≤ N) (τ : ℝ) :
    xiOf N τ ^ N = ((exp τ)⁻¹) ^ 2 := by
  have hNpos : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  unfold xiOf
  rw [← exp_nat_mul, ← exp_neg, ← exp_nat_mul]
  congr 1
  field_simp
  push_cast
  ring

theorem cast_m_eq {N m : ℕ} (hN : 1 ≤ N) : (m : ℝ) = N * (1 - xOf N m) / 2 := by
  have hNpos : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  unfold xOf
  field_simp
  ring

/-- `b_m = N² β(x,τ)` (note Thm 11). -/
theorem bBlock_xiOf {N m : ℕ} (hm1 : 1 ≤ m) (hmN : m < N) {τ : ℝ} (hτ : 0 < τ) :
    bBlock N m (xiOf N τ) = (N : ℝ) ^ 2 * beta (xOf N m) τ := by
  have hN : 1 ≤ N := by omega
  have hNpos : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  unfold bBlock W beta Q
  simp only [pow_add]
  rw [xiOf_pow_N hN, xiOf_pow_m (m := m) hN, xiOf_pow_sub hN hmN.le, Nat.cast_sub hmN.le, cast_m_eq (m := m) hN]
  simp only [cosh_eq, sinh_eq, exp_neg]
  have hE1 : 1 < exp τ := one_lt_exp_iff.mpr hτ
  have hy0 : exp (xOf N m * τ) ≠ 0 := (exp_pos _).ne'
  generalize exp (xOf N m * τ) = y at *
  generalize exp τ = E at *
  generalize xOf N m = x at *
  have hE0 : E ≠ 0 := by positivity
  have hE2 : E ^ 2 - 1 ≠ 0 := by nlinarith
  have hE2' : E * E - 1 ≠ 0 := by nlinarith
  field_simp
  ring

/-- `a_m = N² (E_N(τ) − 1/(4 sinh²τ) − ψ(x,τ))` (note Thm 11). -/
theorem aBlock_xiOf {N m : ℕ} (hm1 : 1 ≤ m) (hmN : m < N) {τ : ℝ} (hτ : 0 < τ) :
    aBlock N m (xiOf N τ) =
      (N : ℝ) ^ 2 * (EN N τ - 1 / (4 * sinh τ ^ 2) - psi (xOf N m) τ) := by
  have hN : 1 ≤ N := by omega
  have hNpos : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have h2 : 2 * N - m = N + (N - m) := by omega
  unfold aBlock EN kN
  rw [← cR_xiOf hN hτ, h2]
  unfold coth psi Q
  simp only [pow_add]
  rw [xiOf_pow_N hN, xiOf_pow_m (m := m) hN, xiOf_pow_sub hN hmN.le, cast_m_eq (m := m) hN]
  simp only [cosh_eq, sinh_eq, exp_neg]
  have hE1 : 1 < exp τ := one_lt_exp_iff.mpr hτ
  have hy0 : exp (xOf N m * τ) ≠ 0 := (exp_pos _).ne'
  generalize cR (xiOf N τ) = c
  generalize exp (xOf N m * τ) = y at *
  generalize exp τ = E at *
  generalize xOf N m = x at *
  have hE0 : E ≠ 0 := by positivity
  have hE2 : E ^ 2 - 1 ≠ 0 := by nlinarith
  have hE2' : E * E - 1 ≠ 0 := by nlinarith
  field_simp
  ring

/-- `g_m = N² γ(x,τ)` with `γ = R/(4 sinh²τ)` (note Thm 11). -/
theorem gBlock_xiOf {N m : ℕ} (hm1 : 1 ≤ m) (hmN : m < N) {τ : ℝ} (hτ : 0 < τ) :
    gBlock N m (xiOf N τ) = (N : ℝ) ^ 2 * gamma (xOf N m) τ := by
  have hN : 1 ≤ N := by omega
  have hNpos : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  unfold gBlock W gamma R
  simp only [pow_add]
  rw [xiOf_pow_N hN, xiOf_pow_m (m := m) hN, xiOf_pow_sub hN hmN.le, Nat.cast_sub hmN.le,
    cast_m_eq (m := m) hN]
  simp only [cosh_eq, sinh_eq, exp_neg]
  have hE1 : 1 < exp τ := one_lt_exp_iff.mpr hτ
  have hy0 : exp (xOf N m * τ) ≠ 0 := (exp_pos _).ne'
  generalize exp (xOf N m * τ) = y at *
  generalize exp τ = E at *
  generalize xOf N m = x at *
  have hE0 : E ≠ 0 := by positivity
  have hE2 : E ^ 2 - 1 ≠ 0 := by nlinarith
  have hE2' : E * E - 1 ≠ 0 := by nlinarith
  field_simp
  ring

end Havelock.Tanh.LargeN
