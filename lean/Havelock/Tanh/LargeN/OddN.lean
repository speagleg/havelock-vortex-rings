import Havelock.Tanh.LargeN.EvenN
import Havelock.Tanh.LargeN.ThresholdLimit

/-! # Tanh model, large `N`: odd `N ≥ 17` has a unique threshold (note Prop oddN)

Note, end of section "Large N", subsection "Odd N: the threshold is unique", Proposition
`prop:oddN`. With `τ = −(N/2) ln ξ`, `x = x_{⌊N/2⌋} ∈ [0, 1/N]` (`= 1/N` for odd `N`) and
`a_crit = a_{⌊N/2⌋}`, the note writes `a_crit/N² = Φ + Δ_N + ε_N` (`Fcrit`), with
`Φ_N = Φ + Δ_N` strictly decreasing (note Prop 15 holds for all `N ≥ 3`: `PhiN_strictAntiOn`) and
`ε_N = ψ(0,τ) − ψ(x,τ) ∈ [0, x²/4]` (`psi_zero_sub_le_quarter`).

Main results, for every `N ≥ 17` (the proof never uses the parity, since only `0 ≤ x ≤ 1/N`
enters; for even `N` it re-proves part of note Prop 15):
* `aCrit_pos_of_le_five`: `a_crit(e^{−2τ/N}) > 0` for `0 < τ ≤ 5`;
* `aCrit_neg_of_ten_le`: `a_crit(e^{−2τ/N}) < 0` for `τ ≥ 10`;
* `aCrit_strictAntiOn`: `τ ↦ a_crit(e^{−2τ/N})` is strictly decreasing on `[5,10]`;
* `aBlock_half_threshold_of_ge_17`: a unique zero `ξ₀ ∈ (0,1)` of `a_crit`, negative on `(0,ξ₀)`,
  positive on `(ξ₀,1)` (the shape of `aBlock_half_threshold` for even `N`);
* `xiStar_unique_of_ge_17`: `ξ_*(N) = e^{−2τ₀/N}` with `τ₀ ∈ (5,10)`, `a_crit` vanishes on `(0,1)`
  only at `ξ_*(N)`, and `specStable N ξ ↔ ξ_*(N) ≤ ξ` for `0 < ξ < 1`.

Constants used (cruder than or equal to the note's, but enough for `N ≥ 17`):
`Φ(5) ≥ 1/20`, `Φ(10) < −3/100`, `Δ_N(5) ≥ −13/N²`, `Δ_N(10) ≤ 10/(3N²)`, and on `[5,10]`:
`Φ' ≤ −1/250` (first term of `Φ'` only), `Δ_N' ≤ (1/2 + 1/15)/N²`, `ε_N' ≤ (3/10) x²`,
so `(a_crit/N²)' ≤ −1/250 + (13/15)/N² < 0` for `N² ≥ 289`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real Set Filter Topology

/-- `𝓛(u) ≤ u/3` for `u > 0`. -/
theorem langevin_le_third {u : ℝ} (hu : 0 < u) : langevin u ≤ u / 3 := by
  have hs : 0 < sinh u := sinh_pos_iff.mpr hu
  have hg : 0 < u ^ 2 * sinh u - 3 * (u * cosh u - sinh u) := by
    refine pos_of_hasDerivAt_pos (f := fun t => t ^ 2 * sinh t - 3 * (t * cosh t - sinh t))
      (f' := fun t => t * (t * cosh t - sinh t)) ?_ (by simp) ?_ hu
    · intro t
      have h := (((hasDerivAt_id t).pow 2).mul (Real.hasDerivAt_sinh t)).sub
        ((((hasDerivAt_id t).mul (Real.hasDerivAt_cosh t)).sub (Real.hasDerivAt_sinh t)).const_mul 3)
      convert h using 1
      · funext y; simp only [id, Pi.sub_apply, Pi.pow_apply, Pi.mul_apply]
      · simp only [id, Pi.pow_apply, Nat.cast_ofNat]; ring
    · intro t ht
      exact mul_pos ht (mul_cosh_sub_sinh_pos ht)
  have : langevin u = (u * cosh u - sinh u) / (u * sinh u) := by
    simp only [langevin, coth]; field_simp
  rw [this, div_le_iff₀ (by positivity)]
  nlinarith

/-- `sinh u ≤ u + (2/9) u³` for `0 ≤ u ≤ 1` (from `Real.exp_bound`). -/
theorem sinh_le_cubic {u : ℝ} (h0 : 0 ≤ u) (h1 : u ≤ 1) : sinh u ≤ u + 2 / 9 * u ^ 3 := by
  have hp := Real.exp_bound (x := u) (by rw [abs_of_nonneg h0]; exact h1) (n := 3) (by norm_num)
  have hm := Real.exp_bound (x := -u) (by rw [abs_neg, abs_of_nonneg h0]; exact h1) (n := 3)
    (by norm_num)
  rw [abs_neg, abs_of_nonneg h0] at hm
  rw [abs_of_nonneg h0] at hp
  norm_num [Finset.sum_range_succ, Nat.factorial] at hp hm
  have h1 := (abs_le.mp hp).2
  have h2 := (abs_le.mp hm).1
  rw [sinh_eq]
  nlinarith

/-- note Prop 15 for every `N ≥ 1`: `8 Φ_N(τ) = evenRHS N τ` (no parity needed). -/
theorem eight_mul_PhiN {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    8 * PhiN N τ = evenRHS N τ := by
  have hN0 : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  unfold PhiN evenRHS muEven EN kN
  generalize coth (τ / N) = c
  obtain ⟨η, rfl⟩ : ∃ η, τ = 2 * η := ⟨τ / 2, by ring⟩
  have hη : 0 < η := by linarith
  have hs : sinh η ≠ 0 := (sinh_pos_iff.mpr hη).ne'
  have hc : cosh η ≠ 0 := (cosh_pos _).ne'
  have hpy := cosh_sq_sub_sinh_sq η
  rw [show 2 * η / 2 = η by ring]
  simp only [psi, Q, coth, tanh_eq_sinh_div_cosh, zero_mul, cosh_zero, sinh_zero, mul_zero, sub_zero,
    mul_one, cosh_two_mul, sinh_two_mul]
  field_simp
  linear_combination (16 * cosh η ^ 2 * (N : ℝ) ^ 2 + 8 * (N : ℝ) ^ 2) * hpy

/-- `e⁵ ≥ 148`. -/
theorem exp_five_ge : (148 : ℝ) ≤ exp 5 := by
  refine le_trans ?_ (Real.sum_le_exp_of_nonneg (by norm_num) 15)
  norm_num [Finset.sum_range_succ, Nat.factorial]

/-- `e^{5/2} ≥ 12`. -/
theorem exp_five_half_ge : (12 : ℝ) ≤ exp (5 / 2) := by
  refine le_trans ?_ (Real.sum_le_exp_of_nonneg (by norm_num) 10)
  norm_num [Finset.sum_range_succ, Nat.factorial]


/-! ### Certified values of `Φ` -/

/-- `Φ(5) ≥ 1/20` (numerically `Φ(5) = 0.0516…`). -/
theorem Phi_five_ge : (1 / 20 : ℝ) ≤ Phi 5 := by
  have hE := exp_five_ge
  have h8 := eight_mul_Phi_t (τ := 5) (by norm_num)
  set E := exp (5 : ℝ)
  have ht0 : (147 / 149 : ℝ) ≤ (E - 1) / (E + 1) := by
    rw [div_le_div_iff₀ (by norm_num) (by linarith)]; linarith
  have ht1 : (E - 1) / (E + 1) < 1 := by
    rw [div_lt_one (by linarith)]; linarith
  have hg := g_mono (η := 5 / 2) (by norm_num) (by norm_num) ht0 ht1
  have hnum : (2 / 5 : ℝ) ≤ 2 * (147 / 149) / (5 / 2) - (1 / (147 / 149) - 1 / (5 / 2)) ^ 2 := by
    norm_num
  norm_num at h8 hg
  linarith

/-- `Φ(10) < −3/100` (numerically `Φ(10) = −0.0300…`). -/
theorem Phi_ten_lt : Phi 10 < -3 / 100 := by
  have h8 := eight_mul_Phi_t (τ := 10) (by norm_num)
  have hE1 : 1 < exp (10 : ℝ) := one_lt_exp_iff.mpr (by norm_num)
  set E := exp (10 : ℝ)
  set t := (E - 1) / (E + 1) with ht
  have ht0 : 0 < t := div_pos (by linarith) (by linarith)
  have ht1 : t < 1 := by rw [ht, div_lt_one (by linarith)]; linarith
  have h1t : 1 < 1 / t := by rw [lt_div_iff₀ ht0]; linarith
  have hsq : (4 / 5 : ℝ) ^ 2 < (1 / t - 1 / (10 / 2)) ^ 2 := by
    apply pow_lt_pow_left₀ (by linarith) (by norm_num) (by norm_num)
  have : 2 * t / (10 / 2) < 2 / 5 := by norm_num; linarith
  linarith

/-! ### `Φ'` on `[5,10]` -/

/-- note Prop oddN (iii): `Φ'(τ) ≤ −1/250` on `[5,10]` (the note has `−Φ' > 0.0052`). -/
theorem PhiDeriv_le {τ : ℝ} (h5 : 5 ≤ τ) (h10 : τ ≤ 10) : PhiDeriv τ ≤ -1 / 250 := by
  set η := τ / 2 with hηdef
  have hη : 0 < η := by rw [hηdef]; linarith
  have hη1 : 5 / 2 ≤ η := by rw [hηdef]; linarith
  have hη2 : η ≤ 5 := by rw [hηdef]; linarith
  have hs : 0 < sinh η := sinh_pos_iff.mpr hη
  have hc : 0 < cosh η := cosh_pos η
  have hsη : η < sinh η := self_lt_sinh_iff.mpr hη
  have h3 : 0 < 1 / η ^ 2 - 1 / sinh η ^ 2 := by
    have : η ^ 2 < sinh η ^ 2 := by nlinarith
    have := one_div_lt_one_div_of_lt (by positivity) this
    linarith
  have h4 := mul_pos (langevin_pos hη) h3
  have hE : 12 ≤ exp η := le_trans exp_five_half_ge (exp_le_exp.mpr hη1)
  have hF : 0 < exp (-η) := exp_pos _
  have hEF : exp η * exp (-η) = 1 := by rw [← exp_add]; simp
  have hS : sinh η = (exp η - exp (-η)) / 2 := sinh_eq η
  have hC : cosh η = (exp η + exp (-η)) / 2 := cosh_eq η
  have hsc : 49 / 50 * cosh η ≤ sinh η := by
    rw [hS, hC]; nlinarith
  have hc6 : 6 ≤ cosh η := by rw [hC]; linarith
  have key : η - sinh η * cosh η ≤ -(83 / 100) * cosh η ^ 2 := by nlinarith
  have hfirst : (η - sinh η * cosh η) / (η ^ 2 * cosh η ^ 2) ≤ -(83 / 100) / 25 := by
    rw [div_le_iff₀ (by positivity)]
    have : η ^ 2 ≤ 25 := by nlinarith
    nlinarith [sq_nonneg (cosh η)]
  simp only [PhiDeriv]
  rw [← hηdef]
  linarith

/-! ### `Δ_N = E_N − E_∞` and its derivative -/

/-- `k = 1/τ + 𝓛(τ/N)/N`. -/
theorem kN_eq {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    kN N τ = 1 / τ + langevin (τ / N) / N := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  unfold kN langevin
  field_simp
  ring

/-- `𝓛' = 1/t² − 1/sinh²t`. -/
theorem hasDerivAt_langevin {t : ℝ} (ht : t ≠ 0) :
    HasDerivAt langevin (-1 / sinh t ^ 2 + 1 / t ^ 2) t := by
  have h := (hasDerivAt_cosh_div_sinh ht).sub ((hasDerivAt_const t (1 : ℝ)).div (hasDerivAt_id t) ht)
  have hf : langevin = fun y => cosh y / sinh y - 1 / id y := by funext y; rfl
  rw [hf]
  convert h using 1
  try simp only [id]
  field_simp
  ring

/-- `δ = k − 1/τ` has derivative `1/τ² − 1/(N² sinh²(τ/N))`. -/
theorem hasDerivAt_delta {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    HasDerivAt (fun t => kN N t - 1 / t) (-1 / sinh (τ / N) ^ 2 / N / N + 1 / τ ^ 2) τ := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hk : HasDerivAt (kN N) (-1 / sinh (τ / N) ^ 2 / N / N) τ :=
    (hasDerivAt_coth_div hNpos hτ).div_const (N : ℝ)
  have h := hk.sub (hasDerivAt_inv hτ.ne')
  have e : (fun t => kN N t - 1 / t) = kN N - fun y => y⁻¹ := by funext t; simp [one_div]
  rw [e]
  refine h.congr_deriv ?_
  ring

/-- The derivative of `Δ_N = E_N − E_∞`: `Δ_N' = δ' (𝓛(τ) − δ) + δ 𝓛'(τ)` with `δ = k − 1/τ`
(proof of note Prop oddN, (iii)). -/
theorem hasDerivAt_Delta {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    HasDerivAt (fun t => EN N t - Einf t)
      ((-1 / sinh (τ / N) ^ 2 / N / N + 1 / τ ^ 2) * (langevin τ - (kN N τ - 1 / τ))
        + (kN N τ - 1 / τ) * (-1 / sinh τ ^ 2 + 1 / τ ^ 2)) τ := by
  have hd := hasDerivAt_delta hN hτ
  have hL := hasDerivAt_langevin hτ.ne'
  have h := (hd.mul (hL.sub (hd.div_const 2))).sub_const (1 / (2 * (N : ℝ) ^ 2))
  have heq : (fun t => EN N t - Einf t) =ᶠ[𝓝 τ]
      (fun t => (kN N t - 1 / t) * (langevin t - (kN N t - 1 / t) / 2) - 1 / (2 * (N : ℝ) ^ 2)) := by
    filter_upwards [Ioi_mem_nhds hτ] with t ht
    exact EN_sub_Einf N ht
  convert h.congr_of_eventuallyEq heq using 1
  simp only [Pi.sub_apply]
  ring

/-! ### `ε = ψ(0,τ) − ψ(x,τ)` and its derivative -/

/-- `∂_τ Q(x,τ) = (1 − x²) sinh τ cosh xτ`. -/
theorem hasDerivAt_Q_tau (x τ : ℝ) :
    HasDerivAt (fun t => Q x t) ((1 - x ^ 2) * sinh τ * cosh (x * τ)) τ := by
  have hx : HasDerivAt (fun t : ℝ => x * t) x τ := by
    simpa using (hasDerivAt_id τ).const_mul x
  have h := ((Real.hasDerivAt_cosh τ).mul ((Real.hasDerivAt_cosh (x * τ)).comp τ hx)).sub
    (((Real.hasDerivAt_sinh τ).const_mul x).mul ((Real.hasDerivAt_sinh (x * τ)).comp τ hx))
  convert h using 1
  · funext t; simp [Q, Function.comp]
  · simp only [Function.comp_apply]; ring

/-- `ε = ψ(0,τ) − ψ(x,τ) = x²/8 + (cosh τ − Q)/(4 sinh²τ)` (note Prop oddN). -/
theorem psi_zero_sub_eq (x t : ℝ) :
    psi 0 t - psi x t = x ^ 2 / 8 + (cosh t - Q x t) / (4 * sinh t ^ 2) := by
  have : Q 0 t = cosh t := by simp [Q]
  unfold psi
  rw [this]
  ring

/-- The `τ`-derivative of `ε(τ) = ψ(0,τ) − ψ(x,τ)`. -/
def epsDeriv (x τ : ℝ) : ℝ :=
  (sinh τ * (1 - (1 - x ^ 2) * cosh (x * τ)) * (4 * sinh τ ^ 2)
    - (cosh τ - Q x τ) * (8 * sinh τ * cosh τ)) / (4 * sinh τ ^ 2) ^ 2

theorem hasDerivAt_eps (x : ℝ) {τ : ℝ} (hτ : 0 < τ) :
    HasDerivAt (fun t => psi 0 t - psi x t) (epsDeriv x τ) τ := by
  have hs : sinh τ ≠ 0 := (sinh_pos_iff.mpr hτ).ne'
  have hD := (Real.hasDerivAt_cosh τ).sub (hasDerivAt_Q_tau x τ)
  have hden := ((Real.hasDerivAt_sinh τ).pow 2).const_mul (4 : ℝ)
  have h := (hD.div hden (mul_ne_zero (by norm_num) (pow_ne_zero 2 hs))).const_add (x ^ 2 / 8)
  have hf : (fun t => psi 0 t - psi x t) = fun t => x ^ 2 / 8 + (cosh t - Q x t) / (4 * sinh t ^ 2) := by
    funext t; exact psi_zero_sub_eq x t
  rw [hf]
  convert h using 1
  simp only [epsDeriv, Pi.pow_apply, Pi.sub_apply, Nat.cast_ofNat]
  ring


/-! ### Elementary bounds on `[5,10]` -/

/-- For `τ ≥ 5`: `cosh τ ≤ (11/10) sinh τ` and `1 + cosh τ ≤ (6/5) sinh τ`. -/
theorem cosh_sinh_bounds {τ : ℝ} (h5 : 5 ≤ τ) :
    cosh τ ≤ 11 / 10 * sinh τ ∧ 1 + cosh τ ≤ 6 / 5 * sinh τ := by
  have hE : 148 ≤ exp τ := le_trans exp_five_ge (exp_le_exp.mpr h5)
  have hF1 : exp (-τ) ≤ 1 := exp_le_one_iff.mpr (by linarith)
  have hF : 0 < exp (-τ) := exp_pos _
  rw [sinh_eq, cosh_eq]
  constructor <;> nlinarith

/-- `𝓛(τ) ≤ 1` on `[5,10]`. -/
theorem langevin_le_one {τ : ℝ} (h5 : 5 ≤ τ) (h10 : τ ≤ 10) : langevin τ ≤ 1 := by
  have hτ : 0 < τ := by linarith
  have hs : 0 < sinh τ := sinh_pos_iff.mpr hτ
  have hc := (cosh_sinh_bounds h5).1
  have h1 : cosh τ / sinh τ ≤ 11 / 10 := by rw [div_le_iff₀ hs]; linarith
  have h2 : 1 / 10 ≤ 1 / τ := one_div_le_one_div_of_le hτ h10
  simp only [langevin, coth]
  linarith

/-- `δ = k − 1/τ = 𝓛(τ/N)/N` satisfies `0 < δ ≤ τ/(3N²)` (proof of note Prop oddN). -/
theorem delta_bounds {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    0 < kN N τ - 1 / τ ∧ kN N τ - 1 / τ ≤ τ / (3 * (N : ℝ) ^ 2) := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hu : 0 < τ / N := by positivity
  rw [kN_eq hN hτ, add_sub_cancel_left]
  refine ⟨div_pos (langevin_pos hu) hNpos, ?_⟩
  have h := langevin_le_third hu
  calc langevin (τ / N) / N ≤ τ / N / 3 / N := div_le_div_of_nonneg_right h hNpos.le
    _ = τ / (3 * (N : ℝ) ^ 2) := by field_simp

/-- `Δ_N ≤ τ/(3N²)` where `𝓛(τ) ≤ 1` (proof of note Prop oddN, (ii)). -/
theorem Delta_le {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 < τ) (hL : langevin τ ≤ 1) :
    EN N τ - Einf τ ≤ τ / (3 * (N : ℝ) ^ 2) := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  obtain ⟨hd0, hd1⟩ := delta_bounds hN hτ
  have hL0 := langevin_pos hτ
  rw [EN_sub_Einf N hτ]
  have : 0 < 1 / (2 * (N : ℝ) ^ 2) := by positivity
  nlinarith

/-- `0 ≤ δ' ≤ 1/(2N²)` for `0 < τ ≤ N` (proof of note Prop oddN, (iii)). -/
theorem ddelta_bounds {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 < τ) (hτN : τ ≤ N) :
    0 ≤ -1 / sinh (τ / N) ^ 2 / N / N + 1 / τ ^ 2 ∧
      -1 / sinh (τ / N) ^ 2 / N / N + 1 / τ ^ 2 ≤ 1 / (2 * (N : ℝ) ^ 2) := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  set u := τ / N with hu_def
  have hu : 0 < u := by positivity
  have hu1 : u ≤ 1 := by rw [hu_def, div_le_one hNpos]; exact hτN
  have hτu : τ = N * u := by rw [hu_def]; field_simp
  have hs : 0 < sinh u := sinh_pos_iff.mpr hu
  have hsu : u ≤ sinh u := self_le_sinh_iff.mpr hu.le
  have hsu2 := sinh_le_cubic hu.le hu1
  have heq : -1 / sinh u ^ 2 / N / N + 1 / τ ^ 2 = (1 / u ^ 2 - 1 / sinh u ^ 2) / (N : ℝ) ^ 2 := by
    rw [hτu]; field_simp; ring
  have hlo : 0 ≤ 1 / u ^ 2 - 1 / sinh u ^ 2 := by
    have := one_div_le_one_div_of_le (by positivity) (pow_le_pow_left₀ hu.le hsu 2)
    linarith
  have hhi : 1 / u ^ 2 - 1 / sinh u ^ 2 ≤ 1 / 2 := by
    rw [div_sub_div _ _ (by positivity) (by positivity), div_le_iff₀ (by positivity)]
    have h1 : sinh u - u ≤ 2 / 9 * u ^ 3 := by linarith
    have hu3 : u ^ 3 ≤ u := by nlinarith
    have h2 : sinh u + u ≤ 20 / 9 * u := by nlinarith
    have h3 : (sinh u - u) * (sinh u + u) ≤ 2 / 9 * u ^ 3 * (20 / 9 * u) :=
      mul_le_mul h1 h2 (by linarith) (by positivity)
    have h4 : u ^ 4 ≤ u ^ 2 * sinh u ^ 2 := by
      have := pow_le_pow_left₀ hu.le hsu 2
      nlinarith
    nlinarith
  rw [heq]
  constructor
  · positivity
  · rw [div_le_div_iff₀ (by positivity) (by positivity)]; nlinarith

/-- `Δ_N' ≤ (17/30)/N²` on `[5,10]` for `N ≥ 10` (proof of note Prop oddN, (iii)). -/
theorem dDelta_le {N : ℕ} (hN : 10 ≤ N) {τ : ℝ} (h5 : 5 ≤ τ) (h10 : τ ≤ 10) :
    (-1 / sinh (τ / N) ^ 2 / N / N + 1 / τ ^ 2) * (langevin τ - (kN N τ - 1 / τ))
        + (kN N τ - 1 / τ) * (-1 / sinh τ ^ 2 + 1 / τ ^ 2) ≤ 17 / 30 / (N : ℝ) ^ 2 := by
  have hτ : 0 < τ := by linarith
  have hNr : (10 : ℝ) ≤ N := by exact_mod_cast hN
  have hNpos : (0 : ℝ) < N := by linarith
  obtain ⟨hdd0, hdd1⟩ := ddelta_bounds (N := N) (by omega) hτ (by linarith)
  obtain ⟨hd0, hd1⟩ := delta_bounds (N := N) (by omega) hτ
  have hL := langevin_le_one h5 h10
  set dd := -1 / sinh (τ / N) ^ 2 / N / N + 1 / τ ^ 2
  set d := kN N τ - 1 / τ
  have e1 : dd * (langevin τ - d) ≤ dd := by nlinarith
  have hsn : 0 ≤ 1 / sinh τ ^ 2 := by positivity
  have e2 : d * (-1 / sinh τ ^ 2 + 1 / τ ^ 2) ≤ d * (1 / τ ^ 2) := by
    apply mul_le_mul_of_nonneg_left _ hd0.le
    rw [neg_div]; linarith
  have e3 : d * (1 / τ ^ 2) ≤ τ / (3 * (N : ℝ) ^ 2) * (1 / τ ^ 2) :=
    mul_le_mul_of_nonneg_right hd1 (by positivity)
  have e4 : τ / (3 * (N : ℝ) ^ 2) * (1 / τ ^ 2) ≤ 1 / (15 * (N : ℝ) ^ 2) := by
    rw [div_mul_div_comm, div_le_div_iff₀ (by positivity) (by positivity)]
    nlinarith [mul_nonneg (mul_nonneg (sq_nonneg (N : ℝ)) hτ.le) (sub_nonneg.mpr h5)]
  have e5 : 1 / (2 * (N : ℝ) ^ 2) + 1 / (15 * (N : ℝ) ^ 2) = 17 / 30 / (N : ℝ) ^ 2 := by
    field_simp; ring
  linarith

/-- `ε' ≤ (3/10) x²` on `τ ≥ 5` for `0 ≤ x ≤ 1` (proof of note Prop oddN, (iii)). -/
theorem epsDeriv_le {x τ : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) (h5 : 5 ≤ τ) :
    epsDeriv x τ ≤ 3 / 10 * x ^ 2 := by
  have hτ : 0 < τ := by linarith
  have hs : 0 < sinh τ := sinh_pos_iff.mpr hτ
  have hc : 0 < cosh τ := cosh_pos τ
  have hcw : 1 ≤ cosh (x * τ) := one_le_cosh _
  have hD : -(x ^ 2 * sinh τ ^ 2 / 2) ≤ cosh τ - Q x τ := by
    have h := psi_le_psi_zero (x := x) (by rw [abs_of_nonneg hx0]; exact hx1) hτ
    have h2 := psi_zero_sub_eq x τ
    have h3 : -(x ^ 2 / 8) ≤ (cosh τ - Q x τ) / (4 * sinh τ ^ 2) := by linarith
    rw [le_div_iff₀ (by positivity)] at h3
    linarith
  have e3 := (cosh_sinh_bounds h5).2
  have hx2 : 0 ≤ 1 - x ^ 2 := by nlinarith
  have e1 : sinh τ * (1 - (1 - x ^ 2) * cosh (x * τ)) ≤ sinh τ * x ^ 2 :=
    mul_le_mul_of_nonneg_left (by nlinarith) hs.le
  have e1' := mul_le_mul_of_nonneg_right e1 (by positivity : (0 : ℝ) ≤ 4 * sinh τ ^ 2)
  have e2 : -(cosh τ - Q x τ) * (8 * sinh τ * cosh τ)
      ≤ x ^ 2 * sinh τ ^ 2 / 2 * (8 * sinh τ * cosh τ) :=
    mul_le_mul_of_nonneg_right (by linarith) (by positivity)
  have e4 := mul_le_mul_of_nonneg_left e3 (by positivity : (0 : ℝ) ≤ 4 * x ^ 2 * sinh τ ^ 3)
  unfold epsDeriv
  rw [div_le_iff₀ (by positivity)]
  nlinarith

/-- note Prop oddN (i): `ψ(0,τ) − ψ(x,τ) ≤ x²/4` for `0 ≤ x ≤ 1`, `τ > 0`. -/
theorem psi_zero_sub_le_quarter {x τ : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) (hτ : 0 < τ) :
    psi 0 τ - psi x τ ≤ x ^ 2 / 4 := by
  rcases hx0.eq_or_lt with rfl | hx
  · simp
  have hs : 0 < sinh τ := sinh_pos_iff.mpr hτ
  have hD : 0 < x ^ 2 / 2 * sinh τ ^ 2 - (cosh τ - Q x τ) := by
    refine pos_of_hasDerivAt_pos (f := fun t => x ^ 2 / 2 * sinh t ^ 2 - (cosh t - Q x t))
      (f' := fun t => x ^ 2 / 2 * (2 * sinh t * cosh t)
        - (sinh t - (1 - x ^ 2) * sinh t * cosh (x * t))) ?_ ?_ ?_ hτ
    · intro t
      have h := (((Real.hasDerivAt_sinh t).pow 2).const_mul (x ^ 2 / 2)).sub
        ((Real.hasDerivAt_cosh t).sub (hasDerivAt_Q_tau x t))
      refine h.congr_deriv ?_
      simp only [Nat.cast_ofNat]
      ring
    · simp [Q]
    · intro t ht
      have hst : 0 < sinh t := sinh_pos_iff.mpr ht
      have hxt : 0 < x * t := mul_pos hx ht
      have hcw1 : 1 < cosh (x * t) := one_lt_cosh.mpr hxt.ne'
      have hcw2 : cosh (x * t) ≤ cosh t := by
        rw [cosh_le_cosh, abs_of_pos hxt, abs_of_pos ht]; nlinarith
      have e : x ^ 2 / 2 * (2 * sinh t * cosh t) - (sinh t - (1 - x ^ 2) * sinh t * cosh (x * t))
          = sinh t * (x ^ 2 * (cosh t - cosh (x * t)) + (cosh (x * t) - 1)) := by ring
      rw [e]
      apply mul_pos hst
      nlinarith
  rw [psi_zero_sub_eq]
  have : (cosh τ - Q x τ) / (4 * sinh τ ^ 2) ≤ x ^ 2 / 8 := by
    rw [div_le_iff₀ (by positivity)]; nlinarith
  linarith

/-! ### The critical block for `N ≥ 17` -/

/-- `a_crit/N² = Φ + Δ_N + ε_N` as a function of `τ` (note Prop oddN). -/
def Fcrit (N : ℕ) (τ : ℝ) : ℝ :=
  Phi τ + (EN N τ - Einf τ) + (psi 0 τ - psi (xOf N (N / 2)) τ)

/-- `a_crit(e^{−2τ/N}) = N² · Fcrit N τ` (note Thm 16 proof / Prop oddN). -/
theorem aBlock_half_eq_Fcrit {N : ℕ} (hN : 2 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    aBlock N (N / 2) (xiOf N τ) = (N : ℝ) ^ 2 * Fcrit N τ :=
  aBlock_half_decomp hN hτ

/-- The derivative of `a_crit/N²` in `τ`: `Φ' + Δ_N' + ε_N'`. -/
theorem hasDerivAt_Fcrit {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    HasDerivAt (Fcrit N)
      (PhiDeriv τ + ((-1 / sinh (τ / N) ^ 2 / N / N + 1 / τ ^ 2) * (langevin τ - (kN N τ - 1 / τ))
        + (kN N τ - 1 / τ) * (-1 / sinh τ ^ 2 + 1 / τ ^ 2)) + epsDeriv (xOf N (N / 2)) τ) τ :=
  ((hasDerivAt_Phi hτ).add (hasDerivAt_Delta hN hτ)).add (hasDerivAt_eps _ hτ)

/-- note Prop oddN (iii): `d/dτ (a_crit/N²) < 0` on `[5,10]` for `N ≥ 17`. -/
theorem Fcrit_deriv_neg {N : ℕ} (hN : 17 ≤ N) {τ : ℝ} (h5 : 5 ≤ τ) (h10 : τ ≤ 10) :
    PhiDeriv τ + ((-1 / sinh (τ / N) ^ 2 / N / N + 1 / τ ^ 2) * (langevin τ - (kN N τ - 1 / τ))
        + (kN N τ - 1 / τ) * (-1 / sinh τ ^ 2 + 1 / τ ^ 2)) + epsDeriv (xOf N (N / 2)) τ < 0 := by
  have hNr : (17 : ℝ) ≤ N := by exact_mod_cast hN
  have hNpos : (0 : ℝ) < N := by linarith
  have h1 := PhiDeriv_le h5 h10
  have h2 := dDelta_le (N := N) (by omega) h5 h10
  obtain ⟨hx0, hx1⟩ := xOf_half_bounds (N := N) (by omega)
  have hx1' : xOf N (N / 2) ≤ 1 := le_trans hx1 (by rw [div_le_one hNpos]; linarith)
  have h3 := epsDeriv_le hx0 hx1' h5
  have hx2 : xOf N (N / 2) ^ 2 ≤ 1 / (N : ℝ) ^ 2 := by
    rw [one_div, ← inv_pow, ← one_div]; exact pow_le_pow_left₀ hx0 hx1 2
  have hN2 : (289 : ℝ) ≤ (N : ℝ) ^ 2 := by nlinarith
  have e : 17 / 30 / (N : ℝ) ^ 2 + 3 / 10 * (1 / (N : ℝ) ^ 2) ≤ 26 / 30 / 289 := by
    rw [show 17 / 30 / (N : ℝ) ^ 2 + 3 / 10 * (1 / (N : ℝ) ^ 2) = 26 / 30 / (N : ℝ) ^ 2 by ring]
    exact div_le_div_of_nonneg_left (by norm_num) (by norm_num) hN2
  nlinarith


/-- note Prop oddN: `a_crit/N²` is strictly decreasing on `[5,10]` for `N ≥ 17`. -/
theorem Fcrit_strictAntiOn {N : ℕ} (hN : 17 ≤ N) : StrictAntiOn (Fcrit N) (Icc 5 10) := by
  refine strictAntiOn_of_deriv_neg (convex_Icc 5 10)
    (fun t ht => (hasDerivAt_Fcrit (N := N) (by omega)
      (by linarith [ht.1] : (0 : ℝ) < t)).continuousAt.continuousWithinAt) ?_
  intro t ht
  rw [interior_Icc] at ht
  rw [(hasDerivAt_Fcrit (N := N) (by omega) (by linarith [ht.1] : (0 : ℝ) < t)).deriv]
  exact Fcrit_deriv_neg hN ht.1.le ht.2.le

/-- `a_crit/N² = Φ_N + ε_N`. -/
theorem Fcrit_eq_PhiN (N : ℕ) (τ : ℝ) :
    Fcrit N τ = PhiN N τ + (psi 0 τ - psi (xOf N (N / 2)) τ) := by
  rw [Fcrit, PhiN_eq]

/-- note Prop 15 (which holds for every `N ≥ 3`): `Φ_N` is strictly decreasing on `(0,∞)`. -/
theorem PhiN_strictAntiOn {N : ℕ} (hN : 3 ≤ N) : StrictAntiOn (PhiN N) (Ioi 0) := by
  intro a ha b hb hab
  have h := evenRHS_strictAntiOn hN ha hb hab
  rw [← eight_mul_PhiN (by omega) ha, ← eight_mul_PhiN (by omega) hb] at h
  linarith

/-- `|x_{⌊N/2⌋}| ≤ 1`. -/
theorem xOf_half_abs_le {N : ℕ} (hN : 1 ≤ N) : |xOf N (N / 2)| ≤ 1 := by
  have hNr : (1 : ℝ) ≤ N := by exact_mod_cast hN
  obtain ⟨hx0, hx1⟩ := xOf_half_bounds (N := N) hN
  rw [abs_of_nonneg hx0]
  exact le_trans hx1 (by rw [div_le_one (by linarith)]; exact hNr)

/-- note Prop oddN (ii): `a_crit/N² > 0` for `0 < τ ≤ 5`, `N ≥ 17`. -/
theorem Fcrit_pos {N : ℕ} (hN : 17 ≤ N) {τ : ℝ} (h0 : 0 < τ) (h5 : τ ≤ 5) : 0 < Fcrit N τ := by
  have hNr : (17 : ℝ) ≤ N := by exact_mod_cast hN
  have heps := psi_le_psi_zero (xOf_half_abs_le (N := N) (by omega)) h0
  have hmono : PhiN N 5 ≤ PhiN N τ := (PhiN_strictAntiOn (N := N) (by omega)).antitoneOn
    (mem_Ioi.mpr h0) (mem_Ioi.mpr (by norm_num)) h5
  have hΔ := (EN_sub_Einf_bounds (N := N) (by omega) (τ := 5) (by norm_num)).1
  have hP := Phi_five_ge
  have hN2 : (289 : ℝ) ≤ (N : ℝ) ^ 2 := by nlinarith
  have hb : -(26 / 578 : ℝ) ≤ -((5 : ℝ) ^ 2 + 1) / (2 * (N : ℝ) ^ 2) := by
    rw [neg_div, neg_le_neg_iff, div_le_div_iff₀ (by positivity) (by norm_num)]; nlinarith
  have h5v := PhiN_eq N 5
  rw [Fcrit_eq_PhiN]
  linarith

/-- note Prop oddN (ii): `a_crit/N² < 0` for `τ ≥ 10`, `N ≥ 17`. -/
theorem Fcrit_neg {N : ℕ} (hN : 17 ≤ N) {τ : ℝ} (h10 : 10 ≤ τ) : Fcrit N τ < 0 := by
  have hNr : (17 : ℝ) ≤ N := by exact_mod_cast hN
  have hτ : 0 < τ := by linarith
  obtain ⟨hx0, hx1⟩ := xOf_half_bounds (N := N) (by omega)
  have hx1' : xOf N (N / 2) ≤ 1 := by
    have := xOf_half_abs_le (N := N) (by omega); rwa [abs_of_nonneg hx0] at this
  have heps := psi_zero_sub_le_quarter hx0 hx1' hτ
  have hmono : PhiN N τ ≤ PhiN N 10 := (PhiN_strictAntiOn (N := N) (by omega)).antitoneOn
    (mem_Ioi.mpr (by norm_num)) (mem_Ioi.mpr hτ) h10
  have hΔ := Delta_le (N := N) (by omega) (τ := 10) (by norm_num) (langevin_le_one (by norm_num) le_rfl)
  have hP := Phi_ten_lt
  have hN2 : (289 : ℝ) ≤ (N : ℝ) ^ 2 := by nlinarith
  have hx2 : xOf N (N / 2) ^ 2 ≤ 1 / (N : ℝ) ^ 2 := by
    rw [one_div, ← inv_pow, ← one_div]; exact pow_le_pow_left₀ hx0 hx1 2
  have e1 : 10 / (3 * (N : ℝ) ^ 2) ≤ 10 / (3 * 289) :=
    div_le_div_of_nonneg_left (by norm_num) (by norm_num) (by linarith)
  have e2 : 1 / (N : ℝ) ^ 2 ≤ 1 / 289 := div_le_div_of_nonneg_left (by norm_num) (by norm_num) hN2
  have h10v := PhiN_eq N 10
  rw [Fcrit_eq_PhiN]
  linarith

/-! ### Main statements (note Prop oddN) -/

/-- note Prop oddN: for `N ≥ 17` (in particular every odd `N ≥ 17`) and `0 < τ ≤ 5`,
`a_crit = a_{⌊N/2⌋}(e^{−2τ/N}) > 0`. -/
theorem aCrit_pos_of_le_five {N : ℕ} (hN : 17 ≤ N) {τ : ℝ} (h0 : 0 < τ) (h5 : τ ≤ 5) :
    0 < aBlock N (N / 2) (xiOf N τ) := by
  rw [aBlock_half_eq_Fcrit (by omega) h0]
  exact mul_pos (pow_pos (Nat.cast_pos.mpr (by omega)) 2) (Fcrit_pos hN h0 h5)

/-- note Prop oddN: for `N ≥ 17` and `τ ≥ 10`, `a_crit = a_{⌊N/2⌋}(e^{−2τ/N}) < 0`. -/
theorem aCrit_neg_of_ten_le {N : ℕ} (hN : 17 ≤ N) {τ : ℝ} (h10 : 10 ≤ τ) :
    aBlock N (N / 2) (xiOf N τ) < 0 := by
  rw [aBlock_half_eq_Fcrit (by omega) (by linarith)]
  exact mul_neg_of_pos_of_neg (pow_pos (Nat.cast_pos.mpr (by omega)) 2) (Fcrit_neg hN h10)

/-- note Prop oddN: for `N ≥ 17`, `τ ↦ a_{⌊N/2⌋}(e^{−2τ/N})` is strictly decreasing on `[5,10]`. -/
theorem aCrit_strictAntiOn {N : ℕ} (hN : 17 ≤ N) :
    StrictAntiOn (fun τ => aBlock N (N / 2) (xiOf N τ)) (Icc 5 10) := by
  intro a ha b hb hab
  simp only
  have ha0 : (0 : ℝ) < a := by linarith [(mem_Icc.mp ha).1]
  have hb0 : (0 : ℝ) < b := by linarith [(mem_Icc.mp hb).1]
  rw [aBlock_half_eq_Fcrit (by omega) ha0, aBlock_half_eq_Fcrit (by omega) hb0]
  exact mul_lt_mul_of_pos_left (Fcrit_strictAntiOn hN ha hb hab)
    (pow_pos (Nat.cast_pos.mpr (by omega)) 2)

/-- For `N ≥ 17`, `a_crit/N²` has exactly one sign change in `τ > 0`, at some `τ₀ ∈ (5,10)`. -/
theorem Fcrit_sign {N : ℕ} (hN : 17 ≤ N) :
    ∃ τ₀ ∈ Ioo (5 : ℝ) 10, Fcrit N τ₀ = 0 ∧ (∀ τ, 0 < τ → τ < τ₀ → 0 < Fcrit N τ) ∧
      (∀ τ, τ₀ < τ → Fcrit N τ < 0) := by
  have hcont : ContinuousOn (Fcrit N) (Icc 5 10) := fun t ht =>
    (hasDerivAt_Fcrit (N := N) (by omega)
      (by linarith [ht.1] : (0 : ℝ) < t)).continuousAt.continuousWithinAt
  have h5 := Fcrit_pos hN (by norm_num : (0 : ℝ) < 5) le_rfl
  have h10 := Fcrit_neg hN (le_refl (10 : ℝ))
  obtain ⟨τ₀, hτ₀, hF0⟩ := intermediate_value_Icc' (by norm_num : (5 : ℝ) ≤ 10) hcont
    ⟨h10.le, h5.le⟩
  have h50 : 5 < τ₀ := lt_of_le_of_ne hτ₀.1 (by rintro rfl; linarith)
  have h010 : τ₀ < 10 := lt_of_le_of_ne hτ₀.2 (by rintro rfl; linarith)
  have hA := Fcrit_strictAntiOn hN
  refine ⟨τ₀, ⟨h50, h010⟩, hF0, ?_, ?_⟩
  · intro τ h0 hlt
    rcases le_or_gt τ 5 with h | h
    · exact Fcrit_pos hN h0 h
    · have := hA ⟨h.le, by linarith⟩ hτ₀ hlt
      rwa [hF0] at this
  · intro τ hlt
    rcases le_or_gt 10 τ with h | h
    · exact Fcrit_neg hN h
    · have := hA hτ₀ ⟨by linarith, h.le⟩ hlt
      rwa [hF0] at this

/-- note Prop oddN, core form: for `N ≥ 17` there is `τ₀ ∈ (5,10)` such that, with
`ξ₀ = e^{−2τ₀/N}`, `a_{⌊N/2⌋}(ξ)` vanishes on `(0,1)` exactly at `ξ₀`, is negative on `(0,ξ₀)` and
positive on `(ξ₀,1)`. -/
theorem aCrit_threshold_tau {N : ℕ} (hN : 17 ≤ N) :
    ∃ τ₀ ∈ Ioo (5 : ℝ) 10,
      (∀ ξ ∈ Ioo (0 : ℝ) 1, aBlock N (N / 2) ξ = 0 ↔ ξ = xiOf N τ₀) ∧
      (∀ ξ ∈ Ioo 0 (xiOf N τ₀), aBlock N (N / 2) ξ < 0) ∧
      (∀ ξ ∈ Ioo (xiOf N τ₀) 1, 0 < aBlock N (N / 2) ξ) := by
  have h1N : 1 ≤ N := by omega
  have hNsq : (0 : ℝ) < (N : ℝ) ^ 2 := pow_pos (Nat.cast_pos.mpr (by omega)) 2
  obtain ⟨τ₀, hτ₀, hF0, hpos, hneg⟩ := Fcrit_sign hN
  have hτ0 : 0 < τ₀ := by linarith [hτ₀.1]
  have hval : ∀ ξ ∈ Ioo (0 : ℝ) 1,
      aBlock N (N / 2) ξ = (N : ℝ) ^ 2 * Fcrit N (tauOf N ξ) := fun ξ hξ => by
    rw [← aBlock_half_eq_Fcrit (by omega) (tauOf_pos h1N hξ.1 hξ.2), xiOf_tauOf h1N hξ.1]
  have hξ0 : xiOf N τ₀ ∈ Ioo (0 : ℝ) 1 := ⟨xiOf_pos N τ₀, xiOf_lt_one h1N hτ0⟩
  have hN' : ∀ ξ ∈ Ioo 0 (xiOf N τ₀), aBlock N (N / 2) ξ < 0 := by
    intro ξ ⟨h0, h⟩
    have h1 : ξ < 1 := lt_trans h hξ0.2
    rw [hval ξ ⟨h0, h1⟩]
    have ht : τ₀ < tauOf N ξ := by
      rw [← xiOf_lt_iff h1N, xiOf_tauOf h1N h0]; exact h
    exact mul_neg_of_pos_of_neg hNsq (hneg _ ht)
  have hP' : ∀ ξ ∈ Ioo (xiOf N τ₀) 1, 0 < aBlock N (N / 2) ξ := by
    intro ξ ⟨h, h1⟩
    have h0 : 0 < ξ := lt_trans hξ0.1 h
    rw [hval ξ ⟨h0, h1⟩]
    have ht : tauOf N ξ < τ₀ := by
      rw [← xiOf_lt_iff h1N, xiOf_tauOf h1N h0]; exact h
    exact mul_pos hNsq (hpos _ (tauOf_pos h1N h0 h1) ht)
  refine ⟨τ₀, hτ₀, ?_, hN', hP'⟩
  intro ξ ⟨h0, h1⟩
  constructor
  · intro hz
    by_contra hne
    rcases lt_or_gt_of_ne hne with hlt | hgt
    · exact (hN' ξ ⟨h0, hlt⟩).ne hz
    · exact (hP' ξ ⟨hgt, h1⟩).ne' hz
  · rintro rfl
    rw [aBlock_half_eq_Fcrit (by omega) hτ0, hF0, mul_zero]

/-- note Prop oddN: for every `N ≥ 17` (in particular every odd `N ≥ 17`), `a_{⌊N/2⌋}` has exactly
one zero `ξ₀ ∈ (0,1)`, where it changes sign from negative to positive (same shape as the even case
`aBlock_half_threshold`). -/
theorem aBlock_half_threshold_of_ge_17 {N : ℕ} (hN : 17 ≤ N) :
    ∃ ξ₀ ∈ Ioo (0 : ℝ) 1, (∀ ξ ∈ Ioo (0 : ℝ) 1, aBlock N (N / 2) ξ = 0 ↔ ξ = ξ₀) ∧
      (∀ ξ ∈ Ioo 0 ξ₀, aBlock N (N / 2) ξ < 0) ∧ (∀ ξ ∈ Ioo ξ₀ 1, 0 < aBlock N (N / 2) ξ) := by
  obtain ⟨τ₀, hτ₀, hz, hn, hp⟩ := aCrit_threshold_tau hN
  exact ⟨xiOf N τ₀, ⟨xiOf_pos N τ₀, xiOf_lt_one (by omega) (by linarith [hτ₀.1])⟩, hz, hn, hp⟩

/-- note Prop oddN (spectral threshold): for every `N ≥ 17` (in particular every odd `N ≥ 17`),
`ξ_*(N) = e^{−2τ₀/N}` for some `τ₀ ∈ (5,10)`, `a_{⌊N/2⌋}` vanishes on `(0,1)` exactly at `ξ_*(N)`,
and for `0 < ξ < 1` the ring is spectrally stable iff `ξ_*(N) ≤ ξ`: the threshold is unique. -/
theorem xiStar_unique_of_ge_17 {N : ℕ} (hN : 17 ≤ N) :
    (∃ τ₀ ∈ Ioo (5 : ℝ) 10, xiStar N = xiOf N τ₀) ∧ 0 < xiStar N ∧ xiStar N < 1 ∧
      (∀ ξ ∈ Ioo (0 : ℝ) 1, aBlock N (N / 2) ξ = 0 ↔ ξ = xiStar N) ∧
      (∀ ξ ∈ Ioo (0 : ℝ) 1, specStable N ξ ↔ xiStar N ≤ ξ) := by
  have h2N : 2 ≤ N := by omega
  obtain ⟨τ₀, hτ₀, hz, hn, hp⟩ := aCrit_threshold_tau hN
  set ξ₀ := xiOf N τ₀ with hξ₀def
  have hξ0 : ξ₀ ∈ Ioo (0 : ℝ) 1 := ⟨xiOf_pos N τ₀, xiOf_lt_one (by omega) (by linarith [hτ₀.1])⟩
  have hle : xiStar N ≤ ξ₀ :=
    thr_le (specStable N) hξ0.1.le hξ0.2.le fun ξ hξ =>
      (specStable_iff h2N (hξ0.1.trans hξ.1) hξ.2).mpr (hp ξ hξ).le
  have hge : ξ₀ ≤ xiStar N := by
    by_contra hlt
    replace hlt := not_le.mp hlt
    have hs0 : 0 ≤ xiStar N := thr_nonneg (specStable N)
    obtain ⟨ξ, hξdef⟩ : ∃ ξ : ℝ, ξ = (xiStar N + ξ₀) / 2 := ⟨_, rfl⟩
    have hξa : xiStar N < ξ := by linarith
    have hξb : ξ < ξ₀ := by linarith
    have hξ1 : ξ < 1 := lt_trans hξb hξ0.2
    have hξpos : 0 < ξ := by linarith [hξ0.1]
    have hst : specStable N ξ := holds_of_thr_lt (specStable N) hξa hξ1
    have := (specStable_iff h2N hξpos hξ1).mp hst
    exact absurd (hn ξ ⟨hξpos, hξb⟩) (not_lt.mpr this)
  have heq : xiStar N = ξ₀ := le_antisymm hle hge
  refine ⟨⟨τ₀, hτ₀, heq⟩, heq ▸ hξ0.1, heq ▸ hξ0.2, heq ▸ hz, ?_⟩
  intro ξ ⟨h0, h1⟩
  rw [specStable_iff h2N h0 h1, heq]
  constructor
  · intro ha
    by_contra hlt
    replace hlt := not_le.mp hlt
    exact absurd (hn ξ ⟨h0, hlt⟩) (not_lt.mpr ha)
  · intro hle'
    rcases hle'.eq_or_lt with heq' | hlt
    · exact ((hz ξ ⟨h0, h1⟩).mpr heq'.symm).ge
    · exact (hp ξ ⟨hlt, h1⟩).le

end Havelock.Tanh.LargeN
