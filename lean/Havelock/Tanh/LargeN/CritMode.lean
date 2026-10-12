import Havelock.Tanh.LargeN.Stiffness
import Mathlib.Analysis.Convex.SpecificFunctions.Basic

/-! # Tanh model, large `N`: the critical mode (note Lemma 13)

`ψ(·,τ)` is even and strictly decreasing on `[0,1]` for every `τ > 0`. Hence, for every `N ≥ 2`,
`1 ≤ m ≤ N−1` and `0 < ξ < 1`, `a_m ≥ a_{⌊N/2⌋}`, with equality only for `m ∈ {⌊N/2⌋, N − ⌊N/2⌋}`.

The proof of the monotonicity differs from the note's (which uses `λ = ln(sinh t/t)` and a case
split at `τ = √3 e^{3/2}`). With `h = xτ ∈ (0,τ)` and `d = τ − h > 0`, the note's quantity
`τ · 4 sinh²τ · ∂ₓψ = τ sinh τ · v(h)` is negative by two elementary inequalities:
* (A) `τ cosh τ sinh h − h cosh h sinh τ ≤ d sinh τ sinh h`, which after `cosh t − sinh t = e^{−t}`
  is `τ (e^{2h} − 1) ≤ h (e^{2τ} − 1)`, i.e. convexity of `exp`;
* (B) `h sinh τ > (d − 1) τ sinh h`, from `sinh τ = sinh h cosh d + cosh h sinh d`,
  `sinh h < h cosh h`, `cosh d ≥ 1 + d²/2` and `sinh d ≥ d (1 + d²/8)`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real Set

/-- `1 + y²/2 ≤ cosh y` (from `cosh y = 1 + 2 sinh²(y/2)` and `sinh u ≥ u` for `u ≥ 0`). -/
theorem one_add_sq_half_le_cosh {y : ℝ} (hy : 0 ≤ y) : 1 + y ^ 2 / 2 ≤ cosh y := by
  have h2 : cosh y = cosh (y / 2) ^ 2 + sinh (y / 2) ^ 2 := by
    rw [← cosh_two_mul]; ring_nf
  have hs : y / 2 ≤ sinh (y / 2) := self_le_sinh_iff.mpr (by linarith)
  have hc := cosh_sq (y / 2)
  nlinarith

/-- `d (1 + d²/8) ≤ sinh d` for `d ≥ 0`. -/
theorem cubic_le_sinh {d : ℝ} (hd : 0 ≤ d) : d * (1 + d ^ 2 / 8) ≤ sinh d := by
  have h2 : sinh d = 2 * sinh (d / 2) * cosh (d / 2) := by
    rw [← sinh_two_mul]; ring_nf
  have hs : d / 2 ≤ sinh (d / 2) := self_le_sinh_iff.mpr (by linarith)
  have hc : 1 + (d / 2) ^ 2 / 2 ≤ cosh (d / 2) := one_add_sq_half_le_cosh (by linarith)
  have hc0 : 0 ≤ 1 + (d / 2) ^ 2 / 2 := by positivity
  rw [h2]
  have := mul_le_mul hs hc hc0 (by linarith)
  nlinarith

/-- Inequality (A): `τ cosh τ sinh h − h cosh h sinh τ ≤ (τ − h) sinh τ sinh h` for `0 < h ≤ τ`
(convexity of `exp`). -/
theorem ineqA {h τ : ℝ} (hh : 0 < h) (hhτ : h ≤ τ) :
    τ * cosh τ * sinh h - h * cosh h * sinh τ ≤ (τ - h) * sinh τ * sinh h := by
  have hτ : 0 < τ := lt_of_lt_of_le hh hhτ
  -- convexity: exp (2h) ≤ (h/τ) exp (2τ) + (1 − h/τ)
  have hconv := convexOn_exp.2 (mem_univ (2 * τ)) (mem_univ 0) (div_nonneg hh.le hτ.le)
    (by rw [sub_nonneg, div_le_one hτ]; exact hhτ) (add_sub_cancel (h / τ) 1)
  simp only [smul_eq_mul, mul_zero, add_zero, exp_zero, mul_one] at hconv
  rw [show h / τ * (2 * τ) = h + h by field_simp; ring] at hconv
  rw [show 2 * τ = τ + τ by ring, exp_add, exp_add] at hconv
  set a := exp h with ha_def
  set b := exp τ with hb_def
  have ha : 0 < a := exp_pos h
  have hb : 0 < b := exp_pos τ
  have key : τ * (a ^ 2 - 1) ≤ h * (b ^ 2 - 1) := by
    have h1 : τ * (a * a) ≤ τ * (h / τ * (b * b) + (1 - h / τ)) :=
      mul_le_mul_of_nonneg_left hconv hτ.le
    have h2 : τ * (h / τ * (b * b) + (1 - h / τ)) = h * (b * b) + (τ - h) := by
      field_simp
    nlinarith
  have hid : (τ - h) * sinh τ * sinh h - (τ * cosh τ * sinh h - h * cosh h * sinh τ)
      = (h * (b ^ 2 - 1) - τ * (a ^ 2 - 1)) / (2 * a * b) := by
    rw [sinh_eq, sinh_eq, cosh_eq, cosh_eq, exp_neg, exp_neg, ← ha_def, ← hb_def]
    field_simp; ring
  have : 0 ≤ (h * (b ^ 2 - 1) - τ * (a ^ 2 - 1)) / (2 * a * b) :=
    div_nonneg (by linarith) (by positivity)
  linarith

/-- Inequality (B): `(τ − h − 1) τ sinh h < h sinh τ` for `0 < h ≤ τ`. -/
theorem ineqB {h τ : ℝ} (hh : 0 < h) (hhτ : h ≤ τ) :
    (τ - h - 1) * τ * sinh h < h * sinh τ := by
  set d := τ - h with hd_def
  have hd : 0 ≤ d := by rw [hd_def]; linarith
  have hτ : τ = h + d := by rw [hd_def]; ring
  have hsh : 0 < sinh h := sinh_pos_iff.mpr hh
  have hsd : 0 ≤ sinh d := sinh_nonneg_iff.mpr hd
  have hhc : sinh h < h * cosh h := sinh_lt_mul_cosh hh
  have hcd : 1 + d ^ 2 / 2 ≤ cosh d := one_add_sq_half_le_cosh hd
  have hsd' : d * (1 + d ^ 2 / 8) ≤ sinh d := cubic_le_sinh hd
  -- `h sinh τ ≥ sinh h (h cosh d + sinh d)`
  have e1 : h * sinh τ ≥ sinh h * (h * cosh d + sinh d) := by
    rw [hτ, sinh_add]
    have := mul_le_mul_of_nonneg_right hhc.le hsd
    nlinarith
  -- `h cosh d + sinh d > (d − 1)(h + d)`
  have e2 : (d - 1) * (h + d) < h * cosh d + sinh d := by
    have : h * (d - 1) < h * cosh d := by
      apply mul_lt_mul_of_pos_left _ hh
      nlinarith [sq_nonneg (d - 1)]
    nlinarith [sq_nonneg (d - 4)]
  have e3 := mul_lt_mul_of_pos_left e2 hsh
  rw [← hτ] at e3
  nlinarith

/-- `x`-derivative of `ψ(x,τ)`. -/
theorem hasDerivAt_psi (x τ : ℝ) :
    HasDerivAt (fun y => psi y τ)
      (-(x / 4) + (τ * cosh τ * sinh (x * τ) - sinh τ * sinh (x * τ)
        - x * τ * sinh τ * cosh (x * τ)) / (4 * sinh τ ^ 2)) x := by
  have h1 : HasDerivAt (fun y : ℝ => y * τ) τ x := by
    simpa using (hasDerivAt_id x).mul_const τ
  have hc := (hasDerivAt_cosh (x * τ)).comp x h1
  have hs := (hasDerivAt_sinh (x * τ)).comp x h1
  have hQ := (hc.const_mul (cosh τ)).sub (((hasDerivAt_id x).mul_const (sinh τ)).mul hs)
  have hP := (((hasDerivAt_id x).pow 2).const_sub 1).div_const 8
  have := hP.add (hQ.div_const (4 * sinh τ ^ 2))
  convert this using 1
  · funext y
    simp only [psi, Q, Function.comp_apply, id, Pi.add_apply, Pi.sub_apply, Pi.mul_apply,
      Pi.pow_apply]
  · simp only [id, Function.comp_apply, Nat.cast_ofNat]
    ring

/-- The numerator of `∂ₓψ` is negative on `0 < x < 1` (`τ > 0`). -/
theorem psi_deriv_neg {x τ : ℝ} (hτ : 0 < τ) (hx0 : 0 < x) (hx1 : x < 1) :
    -(x / 4) + (τ * cosh τ * sinh (x * τ) - sinh τ * sinh (x * τ)
        - x * τ * sinh τ * cosh (x * τ)) / (4 * sinh τ ^ 2) < 0 := by
  have hst : 0 < sinh τ := sinh_pos_iff.mpr hτ
  set h := x * τ with h_def
  have hh : 0 < h := mul_pos hx0 hτ
  have hhτ : h ≤ τ := by rw [h_def]; nlinarith
  have hA := ineqA hh hhτ
  have hB := ineqB hh hhτ
  have hsh : 0 < sinh h := sinh_pos_iff.mpr hh
  -- `τ · numerator < 0`, where numerator `= −x sinh²τ + (…)`
  have hnum : τ * (τ * cosh τ * sinh h - sinh τ * sinh h - h * sinh τ * cosh h)
      < h * sinh τ ^ 2 := by
    have s1 : τ * cosh τ * sinh h - sinh τ * sinh h - h * sinh τ * cosh h
        ≤ (τ - h - 1) * sinh τ * sinh h := by nlinarith
    have s2 := mul_le_mul_of_nonneg_left s1 hτ.le
    have s3 := mul_lt_mul_of_pos_right hB hst
    nlinarith
  have hxτ : x * τ * sinh τ * cosh h = h * sinh τ * cosh h := by rw [h_def]
  rw [hxτ]
  have hden : 0 < 4 * sinh τ ^ 2 := by positivity
  rw [neg_add_lt_iff_lt_add, add_zero, div_lt_iff₀ hden]
  have : h * sinh τ ^ 2 = τ * (x * sinh τ ^ 2) := by rw [h_def]; ring
  rw [this] at hnum
  have := lt_of_mul_lt_mul_left hnum hτ.le
  linarith

/-- note Lemma 13: `ψ(·,τ)` is even. -/
theorem psi_even (x τ : ℝ) : psi (-x) τ = psi x τ := by
  simp only [psi, Q, neg_mul, cosh_neg, sinh_neg]
  ring

/-- note Lemma 13: `ψ(·,τ)` is strictly decreasing on `[0,1]` for every `τ > 0`. -/
theorem psi_strictAntiOn {τ : ℝ} (hτ : 0 < τ) :
    StrictAntiOn (fun x => psi x τ) (Icc 0 1) := by
  refine strictAntiOn_of_deriv_neg (convex_Icc 0 1)
    (fun y _ => (hasDerivAt_psi y τ).continuousAt.continuousWithinAt) ?_
  intro y hy
  rw [interior_Icc] at hy
  rw [(hasDerivAt_psi y τ).deriv]
  exact psi_deriv_neg hτ hy.1 hy.2

/-- `|x_{⌊N/2⌋}| ≤ |x_m|`, strictly unless `m ∈ {⌊N/2⌋, N − ⌊N/2⌋}`, as integers. -/
theorem abs_N_sub_two_mul {N m : ℕ} :
    ((N : ℤ) - 2 * ((N / 2 : ℕ) : ℤ)) ≤ |(N : ℤ) - 2 * (m : ℤ)| ∧
      (m ≠ N / 2 → m ≠ N - N / 2 → m ≤ N →
        ((N : ℤ) - 2 * ((N / 2 : ℕ) : ℤ)) < |(N : ℤ) - 2 * (m : ℤ)|) := by
  constructor
  · rcases abs_cases ((N : ℤ) - 2 * (m : ℤ)) with ⟨h1, h2⟩ | ⟨h1, h2⟩ <;> rw [h1] <;> omega
  · intro h1 h2 h3
    rcases abs_cases ((N : ℤ) - 2 * (m : ℤ)) with ⟨e1, e2⟩ | ⟨e1, e2⟩ <;> rw [e1] <;> omega

/-- `|xOf N m| = |N − 2m| / N`. -/
theorem abs_xOf {N m : ℕ} (hN : 1 ≤ N) :
    |xOf N m| = (|(N : ℤ) - 2 * (m : ℤ)| : ℤ) / (N : ℝ) := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  have : xOf N m = (((N : ℤ) - 2 * (m : ℤ) : ℤ) : ℝ) / N := by
    unfold xOf; push_cast; field_simp
  rw [this, abs_div, abs_of_pos hNpos, Int.cast_abs]

/-- `|xOf N m| < 1` for `1 ≤ m < N`. -/
theorem abs_xOf_lt_one {N m : ℕ} (hm1 : 1 ≤ m) (hmN : m < N) : |xOf N m| < 1 := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  unfold xOf
  have hm1' : (1 : ℝ) ≤ m := by exact_mod_cast hm1
  have hmN' : (m : ℝ) < N := by exact_mod_cast hmN
  have h1 : 0 < 2 * (m : ℝ) / N := by positivity
  have h2 : 2 * (m : ℝ) / N < 2 := by rw [div_lt_iff₀ hNpos]; linarith
  rw [abs_lt]; constructor <;> linarith

/-- `ψ(x_{⌊N/2⌋}) ≥ ψ(x_m)`, strictly unless `m ∈ {⌊N/2⌋, N − ⌊N/2⌋}`. -/
theorem psi_xOf_le {N m : ℕ} (hN : 2 ≤ N) (hm1 : 1 ≤ m) (hmN : m < N) {τ : ℝ} (hτ : 0 < τ) :
    psi (xOf N m) τ ≤ psi (xOf N (N / 2)) τ ∧
      (m ≠ N / 2 → m ≠ N - N / 2 → psi (xOf N m) τ < psi (xOf N (N / 2)) τ) := by
  have hN1 : 1 ≤ N := by omega
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  have hk1 : 1 ≤ N / 2 := by omega
  have hkN : N / 2 < N := by omega
  have ev : ∀ y, psi y τ = psi |y| τ := by
    intro y
    rcases abs_cases y with ⟨h, _⟩ | ⟨h, _⟩
    · rw [h]
    · rw [h, psi_even]
  have memk : |xOf N (N / 2)| ∈ Icc (0 : ℝ) 1 :=
    ⟨abs_nonneg _, (abs_xOf_lt_one hk1 hkN).le⟩
  have memm : |xOf N m| ∈ Icc (0 : ℝ) 1 := ⟨abs_nonneg _, (abs_xOf_lt_one hm1 hmN).le⟩
  have hk_eq : (|(N : ℤ) - 2 * ((N / 2 : ℕ) : ℤ)| : ℤ) = (N : ℤ) - 2 * ((N / 2 : ℕ) : ℤ) :=
    abs_of_nonneg (by omega)
  obtain ⟨hle, hlt⟩ := abs_N_sub_two_mul (N := N) (m := m)
  rw [ev (xOf N m), ev (xOf N (N / 2))]
  constructor
  · apply (psi_strictAntiOn hτ).antitoneOn memk memm
    rw [abs_xOf hN1, abs_xOf hN1, hk_eq]
    apply div_le_div_of_nonneg_right _ hNpos.le
    exact_mod_cast hle
  · intro h1 h2
    apply psi_strictAntiOn hτ memk memm
    rw [abs_xOf hN1, abs_xOf hN1, hk_eq]
    apply div_lt_div_of_pos_right _ hNpos
    exact_mod_cast hlt h1 h2 hmN.le

/-- note Lemma 13 (consequence for the blocks): for `N ≥ 2`, `1 ≤ m ≤ N−1`, `0 < ξ < 1`,
`a_{⌊N/2⌋} ≤ a_m`. -/
theorem aBlock_half_le {N m : ℕ} {ξ : ℝ} (hN : 2 ≤ N) (hm1 : 1 ≤ m) (hmN : m < N)
    (hξ0 : 0 < ξ) (hξ1 : ξ < 1) : aBlock N (N / 2) ξ ≤ aBlock N m ξ := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  set τ := -((N : ℝ) / 2) * log ξ with hτ_def
  have hlog : log ξ < 0 := log_neg hξ0 hξ1
  have hτ : 0 < τ := by rw [hτ_def]; nlinarith
  have hξτ : xiOf N τ = ξ := by
    unfold xiOf
    rw [hτ_def, show -2 * (-((N : ℝ) / 2) * log ξ) / N = log ξ by field_simp, exp_log hξ0]
  rw [← hξτ, aBlock_xiOf hm1 hmN hτ, aBlock_xiOf (by omega) (by omega) hτ]
  have := (psi_xOf_le hN hm1 hmN hτ).1
  have hN2 : (0 : ℝ) < (N : ℝ) ^ 2 := by positivity
  nlinarith

/-- note Lemma 13 (consequence for the blocks): for `N ≥ 2`, `1 ≤ m ≤ N−1`, `0 < ξ < 1`,
`a_{⌊N/2⌋} < a_m` unless `m = ⌊N/2⌋` or `m = N − ⌊N/2⌋`. -/
theorem aBlock_half_lt {N m : ℕ} {ξ : ℝ} (hN : 2 ≤ N) (hm1 : 1 ≤ m) (hmN : m < N)
    (hξ0 : 0 < ξ) (hξ1 : ξ < 1) (hne1 : m ≠ N / 2) (hne2 : m ≠ N - N / 2) :
    aBlock N (N / 2) ξ < aBlock N m ξ := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  set τ := -((N : ℝ) / 2) * log ξ with hτ_def
  have hlog : log ξ < 0 := log_neg hξ0 hξ1
  have hτ : 0 < τ := by rw [hτ_def]; nlinarith
  have hξτ : xiOf N τ = ξ := by
    unfold xiOf
    rw [hτ_def, show -2 * (-((N : ℝ) / 2) * log ξ) / N = log ξ by field_simp, exp_log hξ0]
  rw [← hξτ, aBlock_xiOf hm1 hmN hτ, aBlock_xiOf (by omega) (by omega) hτ]
  have := (psi_xOf_le hN hm1 hmN hτ).2 hne1 hne2
  have hN2 : (0 : ℝ) < (N : ℝ) ^ 2 := by positivity
  nlinarith

end Havelock.Tanh.LargeN
