import Havelock.Tanh.LargeN.RowSymbol
import Havelock.Tanh.LargeN.PhiLimit
import Havelock.Tanh.LargeN.CritMode

/-! # Tanh model, large `N`: critical spacing of the horocyclic row (note Cor. `cor:row`)

For `0 < θ < 1` the Bloch matrix of the row (note Prop. `prop:row`, `RowSymbol.lean`) has real diagonal
`M_xx = τ²β`, `M_yy = τ²α` and purely imaginary `M_xy = iτ²γ` (`τ = 2π/ℓ`). With `ℓ_* = 2π/τ_∞`:

* `row_stable_iff`: `M_xx M_yy ≥ 0` for every `θ ∈ (0,1)` iff `ℓ ≥ ℓ_*` (the note's spectral-stability
  criterion; the reduction of the linearised dynamics to this sign condition is not formalised);
* `row_definite_iff`: `M_xx M_yy − |M_xy|² > 0` for every `θ ∈ (0,1)` iff `ℓ > ℓ_*` (with `M_xx > 0`
  this is positive definiteness of the Hermitian matrix `M(θ)`);
* `rowMyy_half_crit`: at `ℓ = ℓ_*` the pairing mode `θ = 1/2` is degenerate (`M_yy = 0`).

Ingredients: `β > 0` (note Lemma 12), `ψ(·,τ)` even and strictly decreasing on `[0,1]` (note Lemma 13),
`βG = (cosh τ + 1)² K/(32 sinh⁴τ)` with `K ≥ 0` (note Lemma 17), and `Φ` strictly decreasing with
`Φ(τ_∞) = 0`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real
open Complex (I)

theorem rowMxx_eq_ofReal {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    rowMxx ℓ θ = (((2 * π / ℓ) ^ 2 * beta (1 - 2 * θ) (2 * π / ℓ) : ℝ) : ℂ) :=
  (hasSum_termXX hℓ hθ).tsum_eq

theorem rowMxy_eq_ofReal {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    rowMxy ℓ θ = I * (((2 * π / ℓ) ^ 2 * gamma (1 - 2 * θ) (2 * π / ℓ) : ℝ) : ℂ) :=
  (hasSum_termXY hℓ hθ).tsum_eq

theorem rowMyy_eq_ofReal {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    rowMyy ℓ θ = (((2 * π / ℓ) ^ 2 * (Einf (2 * π / ℓ) - 1 / (4 * sinh (2 * π / ℓ) ^ 2)
        - psi (1 - 2 * θ) (2 * π / ℓ)) : ℝ) : ℂ) := by
  have h := (hasSum_termYY hℓ hθ).sub (hasSum_ite_eq (0 : ℤ) (2 * (rowVel ℓ : ℂ)))
  simp only [add_sub_cancel_right] at h
  rw [rowMyy, h.tsum_eq, sub_add_cancel]

/-- `M_xx(θ)` and `M_yy(θ)` are real and `M_xy(θ)` is purely imaginary. -/
theorem rowM_real_parts {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    (rowMxx ℓ θ).im = 0 ∧ (rowMyy ℓ θ).im = 0 ∧ (rowMxy ℓ θ).re = 0 := by
  rw [rowMxx_eq_ofReal hℓ hθ, rowMyy_eq_ofReal hℓ hθ, rowMxy_eq_ofReal hℓ hθ]
  refine ⟨Complex.ofReal_im _, Complex.ofReal_im _, ?_⟩
  simp only [Complex.mul_re, Complex.I_re, Complex.I_im, Complex.ofReal_re, Complex.ofReal_im]
  ring

/-- `ψ(x,τ) ≤ ψ(0,τ)` for `|x| ≤ 1`, with equality only at `x = 0` (note Lemma 13). -/
theorem psi_le_psi_zero_row {x τ : ℝ} (hτ : 0 < τ) (hx : |x| ≤ 1) : psi x τ ≤ psi 0 τ := by
  have h : psi x τ = psi |x| τ := by
    rcases abs_cases x with ⟨h1, _⟩ | ⟨h1, _⟩
    · rw [h1]
    · rw [h1, psi_even]
  rw [h]
  exact (psi_strictAntiOn hτ).antitoneOn ⟨le_rfl, zero_le_one⟩ ⟨abs_nonneg x, hx⟩ (abs_nonneg x)

/-- `Φ(τ) ≥ 0 ↔ τ ≤ τ_∞` on `τ > 0`. -/
theorem Phi_nonneg_iff {τ : ℝ} (hτ : 0 < τ) : 0 ≤ Phi τ ↔ τ ≤ tauInf := by
  constructor
  · intro h
    by_contra hlt
    replace hlt := not_le.mp hlt
    have := Phi_strictAntiOn (show tauInf ∈ Set.Ioi 0 from tauInf_pos) (show τ ∈ Set.Ioi 0 from hτ) hlt
    rw [Phi_tauInf] at this
    linarith
  · intro h
    rw [← Phi_tauInf]
    exact Phi_strictAntiOn.antitoneOn hτ (show tauInf ∈ Set.Ioi 0 from tauInf_pos) h

/-- `Φ(τ) > 0 ↔ τ < τ_∞` on `τ > 0`. -/
theorem Phi_pos_iff {τ : ℝ} (hτ : 0 < τ) : 0 < Phi τ ↔ τ < tauInf := by
  constructor
  · intro h
    by_contra hle
    replace hle := not_lt.mp hle
    have := Phi_strictAntiOn.antitoneOn (show tauInf ∈ Set.Ioi 0 from tauInf_pos) hτ hle
    rw [Phi_tauInf] at this
    linarith
  · intro h
    rw [← Phi_tauInf]
    exact Phi_strictAntiOn hτ (show tauInf ∈ Set.Ioi 0 from tauInf_pos) h

theorem tau_le_iff {ℓ : ℝ} (hℓ : 0 < ℓ) : 2 * π / ℓ ≤ tauInf ↔ 2 * π / tauInf ≤ ℓ := by
  have ht := tauInf_pos
  rw [div_le_iff₀ hℓ, div_le_iff₀ ht, mul_comm tauInf ℓ]

theorem tau_lt_iff {ℓ : ℝ} (hℓ : 0 < ℓ) : 2 * π / ℓ < tauInf ↔ 2 * π / tauInf < ℓ := by
  have ht := tauInf_pos
  rw [div_lt_iff₀ hℓ, div_lt_iff₀ ht, mul_comm tauInf ℓ]

theorem abs_one_sub_two_mul_lt {θ : ℝ} (hθ : θ ∈ Set.Ioo (0 : ℝ) 1) : |1 - 2 * θ| < 1 := by
  rw [abs_lt]; constructor <;> linarith [hθ.1, hθ.2]

theorem gamma_zero (τ : ℝ) : gamma 0 τ = 0 := by simp [gamma, R]

/-- `α(x,τ) = Φ(τ) + ψ(0,τ) − ψ(x,τ)`. -/
theorem alpha_eq (x τ : ℝ) :
    Einf τ - 1 / (4 * sinh τ ^ 2) - psi x τ = Phi τ + (psi 0 τ - psi x τ) := by
  rw [Phi]; ring

/-- `αβ − γ² = β(Φ + G)` and `βG ≥ 0` for `|x| < 1` (note Lemma 17), hence `αβ − γ² ≥ βΦ`. -/
theorem det_ge {x τ : ℝ} (hx : |x| < 1) (hτ : 0 < τ) :
    beta x τ * Phi τ ≤
      beta x τ * (Einf τ - 1 / (4 * sinh τ ^ 2) - psi x τ) - gamma x τ ^ 2 := by
  have hβ := beta_pos hx hτ
  have hBG := beta_mul_G hτ hβ.ne'
  have hK : 0 ≤ K x τ := by
    rcases eq_or_ne x 0 with h0 | h0
    · subst h0; simp [K]
    · exact (K_pos (abs_pos.mpr h0) hx hτ).le
  have hG : 0 ≤ beta x τ * G x τ := by rw [hBG]; positivity
  have : beta x τ * G x τ = beta x τ * (psi 0 τ - psi x τ) - gamma x τ ^ 2 := by
    rw [G]; field_simp
  rw [alpha_eq]
  nlinarith

/-- note Cor. `cor:row`, stability part: `M_xx(θ) M_yy(θ) ≥ 0` for every `0 < θ < 1` iff
`ℓ ≥ ℓ_* = 2π/τ_∞`. -/
theorem row_stable_iff {ℓ : ℝ} (hℓ : 0 < ℓ) :
    (∀ θ ∈ Set.Ioo (0 : ℝ) 1, 0 ≤ (rowMxx ℓ θ).re * (rowMyy ℓ θ).re) ↔ 2 * π / tauInf ≤ ℓ := by
  set τ := 2 * π / ℓ with hτdef
  have hτ : 0 < τ := by positivity
  rw [← tau_le_iff hℓ, ← hτdef, ← Phi_nonneg_iff hτ]
  constructor
  · intro h
    have h2 := h (1 / 2) ⟨by norm_num, by norm_num⟩
    have hθ : (1 / 2 : ℝ) ∈ Set.Icc (0 : ℝ) 1 := ⟨by norm_num, by norm_num⟩
    rw [rowMxx_eq_ofReal hℓ hθ, rowMyy_eq_ofReal hℓ hθ, Complex.ofReal_re, Complex.ofReal_re,
      ← hτdef, show (1 : ℝ) - 2 * (1 / 2) = 0 by norm_num, alpha_eq, sub_self, add_zero] at h2
    have hβ := beta_pos (x := 0) (by simp) hτ
    have hτ2 : 0 < τ ^ 2 := by positivity
    by_contra hneg
    replace hneg := not_le.mp hneg
    have : τ ^ 2 * beta 0 τ * (τ ^ 2 * Phi τ) < 0 :=
      mul_neg_of_pos_of_neg (by positivity) (mul_neg_of_pos_of_neg hτ2 hneg)
    linarith
  · intro hPhi θ hθ
    have hθ' : θ ∈ Set.Icc (0 : ℝ) 1 := Set.Ioo_subset_Icc_self hθ
    have hx := abs_one_sub_two_mul_lt hθ
    rw [rowMxx_eq_ofReal hℓ hθ', rowMyy_eq_ofReal hℓ hθ', Complex.ofReal_re, Complex.ofReal_re,
      ← hτdef, alpha_eq]
    have hβ := beta_pos hx hτ
    have hψ := psi_le_psi_zero_row (x := 1 - 2 * θ) hτ hx.le
    have : 0 ≤ Phi τ + (psi 0 τ - psi (1 - 2 * θ) τ) := by linarith
    positivity

/-- note Cor. `cor:row`, definiteness part: `M_xx(θ) M_yy(θ) − |M_xy(θ)|² > 0` for every `0 < θ < 1` iff
`ℓ > ℓ_* = 2π/τ_∞`. -/
theorem row_definite_iff {ℓ : ℝ} (hℓ : 0 < ℓ) :
    (∀ θ ∈ Set.Ioo (0 : ℝ) 1,
        0 < (rowMxx ℓ θ).re * (rowMyy ℓ θ).re - Complex.normSq (rowMxy ℓ θ)) ↔
      2 * π / tauInf < ℓ := by
  set τ := 2 * π / ℓ with hτdef
  have hτ : 0 < τ := by positivity
  rw [← tau_lt_iff hℓ, ← hτdef, ← Phi_pos_iff hτ]
  have key : ∀ θ ∈ Set.Icc (0 : ℝ) 1,
      (rowMxx ℓ θ).re * (rowMyy ℓ θ).re - Complex.normSq (rowMxy ℓ θ) =
        τ ^ 4 * (beta (1 - 2 * θ) τ * (Einf τ - 1 / (4 * sinh τ ^ 2) - psi (1 - 2 * θ) τ)
          - gamma (1 - 2 * θ) τ ^ 2) := by
    intro θ hθ
    rw [rowMxx_eq_ofReal hℓ hθ, rowMyy_eq_ofReal hℓ hθ, rowMxy_eq_ofReal hℓ hθ,
      Complex.ofReal_re, Complex.ofReal_re, Complex.normSq_mul, Complex.normSq_I,
      Complex.normSq_ofReal, ← hτdef]
    ring
  constructor
  · intro h
    have hθ : (1 / 2 : ℝ) ∈ Set.Icc (0 : ℝ) 1 := ⟨by norm_num, by norm_num⟩
    have h2 := h (1 / 2) ⟨by norm_num, by norm_num⟩
    rw [key _ hθ, show (1 : ℝ) - 2 * (1 / 2) = 0 by norm_num, alpha_eq, sub_self, add_zero,
      gamma_zero] at h2
    have hβ := beta_pos (x := 0) (by simp) hτ
    by_contra hneg
    replace hneg := not_lt.mp hneg
    have : τ ^ 4 * (beta 0 τ * Phi τ - 0 ^ 2) ≤ 0 := by
      have : beta 0 τ * Phi τ ≤ 0 := mul_nonpos_of_nonneg_of_nonpos hβ.le hneg
      have : (0 : ℝ) < τ ^ 4 := by positivity
      nlinarith
    linarith
  · intro hPhi θ hθ
    have hθ' : θ ∈ Set.Icc (0 : ℝ) 1 := Set.Ioo_subset_Icc_self hθ
    have hx := abs_one_sub_two_mul_lt hθ
    rw [key _ hθ']
    have hβ := beta_pos hx hτ
    have hd := det_ge hx hτ
    have : 0 < beta (1 - 2 * θ) τ * Phi τ := mul_pos hβ hPhi
    have : (0 : ℝ) < τ ^ 4 := by positivity
    nlinarith

/-- note Cor. `cor:row`: at the critical spacing `ℓ = ℓ_*` the pairing mode `θ = 1/2` is degenerate,
`M_yy(1/2) = 0` (so `M(1/2)` is singular, `M_xy(1/2) = 0`). -/
theorem rowMyy_half_crit :
    rowMyy (2 * π / tauInf) (1 / 2) = 0 ∧ rowMxy (2 * π / tauInf) (1 / 2) = 0 := by
  have ht := tauInf_pos
  have hℓ : 0 < 2 * π / tauInf := by positivity
  have hθ : (1 / 2 : ℝ) ∈ Set.Icc (0 : ℝ) 1 := ⟨by norm_num, by norm_num⟩
  have hτ : 2 * π / (2 * π / tauInf) = tauInf := by field_simp
  rw [rowMyy_eq_ofReal hℓ hθ, rowMxy_eq_ofReal hℓ hθ, hτ,
    show (1 : ℝ) - 2 * (1 / 2) = 0 by norm_num, alpha_eq, sub_self, add_zero, Phi_tauInf,
    gamma_zero]
  simp

end Havelock.Tanh.LargeN
