import Havelock.Tanh.LargeN.EMCore
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.Calculus.Deriv.Slope
import Mathlib.Topology.Order.IntermediateValue

/-! # Tanh model, large `N`: the limit function `Φ` and the constant `τ∞` (note eq. (phi))

With `η = τ/2`, `8 Φ(τ) = 2 tanh η / η − 𝓛(η)²` (note eq. (phi)). We prove: `Φ` has a negative
derivative at every `τ > 0` (so it is strictly decreasing on `(0,∞)`), `Φ → 1/4` as `τ → 0⁺`,
`Φ → −1/8` as `τ → ∞`, and `Φ` has a unique zero `τ∞ > 0`, with `7.4453 < τ∞ < 7.4454`
(numerically `τ∞ = 7.44533924494901868722…`). We define `tauInf` and `kappaInf = 1/τ∞`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real Set Filter Topology

/-- The right-hand side of note eq. (phi), divided by 8, as a function of `η = τ/2`
(with `tanh η = sinh η / cosh η`, `𝓛(η) = cosh η / sinh η − 1/η`). -/
def PhiEta (η : ℝ) : ℝ := (2 * (sinh η / cosh η) / η - (cosh η / sinh η - 1 / η) ^ 2) / 8

/-- The derivative of `Φ` (note eq. (phi)): with `η = τ/2`,
`Φ'(τ) = [(η − sinh η cosh η)/(η² cosh² η) − 𝓛(η) (1/η² − 1/sinh² η)] / 8`. -/
def PhiDeriv (τ : ℝ) : ℝ :=
  ((τ / 2 - sinh (τ / 2) * cosh (τ / 2)) / ((τ / 2) ^ 2 * cosh (τ / 2) ^ 2)
    - langevin (τ / 2) * (1 / (τ / 2) ^ 2 - 1 / sinh (τ / 2) ^ 2)) / 8

/-- note eq. (phi): `8 Φ(τ) = 2 tanh(τ/2)/(τ/2) − 𝓛(τ/2)²` for `τ > 0`. -/
theorem eight_mul_Phi {τ : ℝ} (hτ : 0 < τ) :
    8 * Phi τ = 2 * tanh (τ / 2) / (τ / 2) - langevin (τ / 2) ^ 2 := by
  have hη : 0 < τ / 2 := by positivity
  have hs : sinh (τ / 2) ≠ 0 := (sinh_pos_iff.mpr hη).ne'
  have hc : cosh (τ / 2) ≠ 0 := (cosh_pos _).ne'
  have hS : sinh τ ≠ 0 := (sinh_pos_iff.mpr hτ).ne'
  have h2 : τ = 2 * (τ / 2) := by ring
  have hpy := cosh_sq_sub_sinh_sq (τ / 2)
  simp only [Phi, Einf, psi, Q, coth, langevin, tanh_eq_sinh_div_cosh, zero_mul, cosh_zero,
    sinh_zero, mul_zero, sub_zero, mul_one]
  rw [h2] at hS ⊢
  rw [cosh_two_mul, sinh_two_mul] at *
  field_simp
  linear_combination (16 * cosh (τ / 2) ^ 2 * τ ^ 2 + 8 * τ ^ 2) * hpy

theorem Phi_eq_PhiEta {τ : ℝ} (hτ : 0 < τ) : Phi τ = PhiEta (τ / 2) := by
  have h := eight_mul_Phi hτ
  simp only [PhiEta, langevin, coth, tanh_eq_sinh_div_cosh] at h ⊢
  linarith


/-! ## The derivative -/

theorem hasDerivAt_sinh_div_cosh (η : ℝ) :
    HasDerivAt (fun y => sinh y / cosh y) (1 / cosh η ^ 2) η := by
  have h := (Real.hasDerivAt_sinh η).div (Real.hasDerivAt_cosh η) (cosh_pos η).ne'
  convert h using 1
  rw [show cosh η * cosh η - sinh η * sinh η = 1 by linear_combination cosh_sq_sub_sinh_sq η]

theorem hasDerivAt_cosh_div_sinh {η : ℝ} (hη : η ≠ 0) :
    HasDerivAt (fun y => cosh y / sinh y) (-1 / sinh η ^ 2) η := by
  have hs : sinh η ≠ 0 := fun h => hη (sinh_eq_zero.mp h)
  have h := (Real.hasDerivAt_cosh η).div (Real.hasDerivAt_sinh η) hs
  convert h using 1
  rw [show sinh η * sinh η - cosh η * cosh η = -1 by linear_combination -cosh_sq_sub_sinh_sq η]

theorem hasDerivAt_PhiEta {η : ℝ} (hη : 0 < η) :
    HasDerivAt PhiEta
      (((η - sinh η * cosh η) / (η ^ 2 * cosh η ^ 2)
        - langevin η * (1 / η ^ 2 - 1 / sinh η ^ 2)) / 4) η := by
  have hη0 : η ≠ 0 := hη.ne'
  have hs : sinh η ≠ 0 := (sinh_pos_iff.mpr hη).ne'
  have hc : cosh η ≠ 0 := (cosh_pos η).ne'
  have h := ((((hasDerivAt_sinh_div_cosh η).const_mul 2).div (hasDerivAt_id η) hη0).sub
    (((hasDerivAt_cosh_div_sinh hη0).sub
      ((hasDerivAt_const η (1 : ℝ)).div (hasDerivAt_id η) hη0)).pow 2)).div_const 8
  convert h using 1
  · funext y; simp only [PhiEta, id, Pi.sub_apply, Pi.div_apply, Pi.pow_apply]
  · simp only [langevin, coth, id, Pi.sub_apply, Pi.div_apply, Nat.cast_ofNat, Nat.add_one_sub_one, pow_one]
    field_simp
    ring

/-- `Φ` is differentiable at every `τ > 0` with derivative `PhiDeriv τ` (note eq. (phi)). -/
theorem hasDerivAt_Phi {τ : ℝ} (hτ : 0 < τ) : HasDerivAt Phi (PhiDeriv τ) τ := by
  have hη : 0 < τ / 2 := by positivity
  have h1 : HasDerivAt (fun t : ℝ => t / 2) (1 / 2) τ := by
    simpa using (hasDerivAt_id τ).div_const 2
  have h := (hasDerivAt_PhiEta hη).comp τ h1
  have heq : Phi =ᶠ[𝓝 τ] (PhiEta ∘ fun t : ℝ => t / 2) := by
    filter_upwards [Ioi_mem_nhds hτ] with t ht
    exact Phi_eq_PhiEta ht
  convert h.congr_of_eventuallyEq heq using 1
  simp only [PhiDeriv]
  ring

/-- For `η > 0`: `η cosh η − sinh η > 0` (so `𝓛(η) > 0`). -/
theorem mul_cosh_sub_sinh_pos {η : ℝ} (hη : 0 < η) : 0 < η * cosh η - sinh η := by
  refine pos_of_hasDerivAt_pos (f := fun t => t * cosh t - sinh t) (f' := fun t => t * sinh t)
    ?_ (by simp) ?_ hη
  · intro t
    have h := ((hasDerivAt_id t).mul (Real.hasDerivAt_cosh t)).sub (Real.hasDerivAt_sinh t)
    convert h using 1
    · funext y; simp only [id, Pi.sub_apply, Pi.mul_apply]
    · simp only [id]; ring
  · intro t ht
    have := sinh_pos_iff.mpr ht
    positivity

/-- `𝓛(η) > 0` for `η > 0`. -/
theorem langevin_pos {η : ℝ} (hη : 0 < η) : 0 < langevin η := by
  have hs : 0 < sinh η := sinh_pos_iff.mpr hη
  have h := mul_cosh_sub_sinh_pos hη
  have : langevin η = (η * cosh η - sinh η) / (η * sinh η) := by
    simp only [langevin, coth]; field_simp
  rw [this]; positivity

/-- `𝓛(η) ≤ η` for `η > 0`. -/
theorem langevin_le_self {η : ℝ} (hη : 0 < η) : langevin η ≤ η := by
  have hs : 0 < sinh η := sinh_pos_iff.mpr hη
  have hg : 0 < η ^ 2 * sinh η - η * cosh η + sinh η := by
    refine pos_of_hasDerivAt_pos (f := fun t => t ^ 2 * sinh t - t * cosh t + sinh t)
      (f' := fun t => t * sinh t + t ^ 2 * cosh t) ?_ (by simp) ?_ hη
    · intro t
      have h := ((((hasDerivAt_id t).pow 2).mul (Real.hasDerivAt_sinh t)).sub
        ((hasDerivAt_id t).mul (Real.hasDerivAt_cosh t))).add (Real.hasDerivAt_sinh t)
      convert h using 1
      · funext y; simp only [id, Pi.sub_apply, Pi.pow_apply, Pi.mul_apply, Pi.add_apply]
      · simp only [id, Pi.pow_apply, Nat.cast_ofNat]; ring
    · intro t ht
      have := sinh_pos_iff.mpr ht
      positivity
  have : langevin η = (η * cosh η - sinh η) / (η * sinh η) := by
    simp only [langevin, coth]; field_simp
  rw [this, div_le_iff₀ (by positivity)]
  nlinarith

/-- note eq. (phi): `Φ'(τ) < 0` for every `τ > 0`. -/
theorem PhiDeriv_neg {τ : ℝ} (hτ : 0 < τ) : PhiDeriv τ < 0 := by
  set η := τ / 2 with hηdef
  have hη : 0 < η := by positivity
  have hs : 0 < sinh η := sinh_pos_iff.mpr hη
  have hsη : η < sinh η := self_lt_sinh_iff.mpr hη
  have hc : 1 ≤ cosh η := one_le_cosh η
  have h1 : η - sinh η * cosh η < 0 := by nlinarith
  have h2 : (η - sinh η * cosh η) / (η ^ 2 * cosh η ^ 2) < 0 :=
    div_neg_of_neg_of_pos h1 (by positivity)
  have h3 : 0 < 1 / η ^ 2 - 1 / sinh η ^ 2 := by
    have : η ^ 2 < sinh η ^ 2 := by nlinarith
    have := one_div_lt_one_div_of_lt (by positivity) this
    linarith
  have h4 := mul_pos (langevin_pos hη) h3
  simp only [PhiDeriv]
  rw [← hηdef]
  linarith

/-- note eq. (phi): `Φ` has a negative derivative at every `τ > 0`. -/
theorem hasDerivAt_Phi_neg {τ : ℝ} (hτ : 0 < τ) :
    ∃ d, HasDerivAt Phi d τ ∧ d < 0 :=
  ⟨PhiDeriv τ, hasDerivAt_Phi hτ, PhiDeriv_neg hτ⟩

theorem continuousOn_Phi : ContinuousOn Phi (Ioi 0) :=
  fun _ ht => (hasDerivAt_Phi ht).continuousAt.continuousWithinAt

/-- note "The limit": `Φ` is strictly decreasing on `(0,∞)`. -/
theorem Phi_strictAntiOn : StrictAntiOn Phi (Ioi 0) := by
  refine strictAntiOn_of_deriv_neg (convex_Ioi 0) continuousOn_Phi ?_
  intro t ht
  rw [interior_Ioi] at ht
  rw [(hasDerivAt_Phi ht).deriv]
  exact PhiDeriv_neg ht

/-! ## The limits -/

/-- `Φ(τ) → 1/4` as `τ → 0⁺` (note "The limit"). -/
theorem tendsto_Phi_zero : Tendsto Phi (𝓝[>] 0) (𝓝 (1 / 4)) := by
  -- `sinh η / η → 1`
  have hsl : Tendsto (fun η : ℝ => sinh η / η) (𝓝[>] 0) (𝓝 1) := by
    have h := hasDerivAt_iff_tendsto_slope.mp (Real.hasDerivAt_sinh 0)
    rw [cosh_zero] at h
    refine (h.mono_left (nhdsWithin_mono _ ?_)).congr' ?_
    · intro y hy; exact ne_of_gt hy
    · filter_upwards [self_mem_nhdsWithin] with y hy
      simp [slope_def_field, div_eq_inv_mul]
  -- `𝓛(η) → 0`
  have hL : Tendsto langevin (𝓝[>] 0) (𝓝 0) := by
    have hid : Tendsto (fun η : ℝ => η) (𝓝[>] 0) (𝓝 0) := nhdsWithin_le_nhds
    refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hid ?_ ?_
    · filter_upwards [self_mem_nhdsWithin] with y hy
      exact (langevin_pos hy).le
    · filter_upwards [self_mem_nhdsWithin] with y hy
      exact langevin_le_self hy
  have hcosh : Tendsto cosh (𝓝[>] 0) (𝓝 1) := by
    have := (Real.continuous_cosh.tendsto 0).mono_left (nhdsWithin_le_nhds (s := Ioi 0))
    simpa using this
  have hE : Tendsto PhiEta (𝓝[>] 0) (𝓝 (1 / 4)) := by
    have h := (((hsl.div hcosh one_ne_zero).const_mul 2).sub (hL.pow 2)).div_const 8
    have h' : (2 * (1 / 1) - (0 : ℝ) ^ 2) / 8 = 1 / 4 := by norm_num
    rw [h'] at h
    refine h.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with y hy
    simp only [PhiEta, langevin, coth, Pi.div_apply]
    have : y ≠ 0 := ne_of_gt hy
    have : cosh y ≠ 0 := (cosh_pos y).ne'
    field_simp
  have hhalf : Tendsto (fun t : ℝ => t / 2) (𝓝[>] 0) (𝓝[>] 0) := by
    refine tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within _ ?_ ?_
    · have := ((continuous_id.div_const (2 : ℝ)).tendsto 0).mono_left
        (nhdsWithin_le_nhds (s := Ioi 0))
      simpa using this
    · filter_upwards [self_mem_nhdsWithin] with y hy
      show (0 : ℝ) < y / 2; exact half_pos hy
  refine (hE.comp hhalf).congr' ?_
  filter_upwards [self_mem_nhdsWithin] with y hy
  exact (Phi_eq_PhiEta hy).symm

/-- `Φ(τ) → −1/8` as `τ → ∞` (note "The limit"). -/
theorem tendsto_Phi_atTop : Tendsto Phi atTop (𝓝 (-1 / 8)) := by
  have hu : Tendsto (fun η : ℝ => exp (-η)) atTop (𝓝 0) := Real.tendsto_exp_neg_atTop_nhds_zero
  have hu2 : Tendsto (fun η : ℝ => exp (-η) ^ 2) atTop (𝓝 0) := by simpa using hu.pow 2
  have htanh : Tendsto (fun η : ℝ => (1 - exp (-η) ^ 2) / (1 + exp (-η) ^ 2)) atTop (𝓝 1) := by
    have := (tendsto_const_nhds (x := (1 : ℝ))).sub hu2 |>.div
      ((tendsto_const_nhds (x := (1 : ℝ))).add hu2) (by norm_num)
    simp only [sub_zero, add_zero, div_one] at this
    exact this
  have hcoth : Tendsto (fun η : ℝ => (1 + exp (-η) ^ 2) / (1 - exp (-η) ^ 2)) atTop (𝓝 1) := by
    have := (tendsto_const_nhds (x := (1 : ℝ))).add hu2 |>.div
      ((tendsto_const_nhds (x := (1 : ℝ))).sub hu2) (by norm_num)
    simp only [sub_zero, add_zero, div_one] at this
    exact this
  have hinv : Tendsto (fun η : ℝ => η⁻¹) atTop (𝓝 0) := tendsto_inv_atTop_zero
  have hE : Tendsto PhiEta atTop (𝓝 (-1 / 8)) := by
    have h := (((htanh.mul hinv).const_mul 2).sub ((hcoth.sub hinv).pow 2)).div_const 8
    have h' : (2 * (1 * 0) - (1 - (0 : ℝ)) ^ 2) / 8 = -1 / 8 := by norm_num
    rw [h'] at h
    refine h.congr' ?_
    filter_upwards [eventually_gt_atTop 0] with y hy
    have hw : 0 < exp y := exp_pos y
    have hw1 : 1 < exp y := one_lt_exp_iff.mpr hy
    have hneg : exp (-y) = (exp y)⁻¹ := exp_neg y
    have hs : sinh y ≠ 0 := (sinh_pos_iff.mpr hy).ne'
    have hd : 1 - ((exp y)⁻¹) ^ 2 ≠ 0 := by
      have : (exp y)⁻¹ < 1 := inv_lt_one_of_one_lt₀ hw1
      have : 0 < (exp y)⁻¹ := inv_pos.mpr hw
      nlinarith
    have hd' : exp y ^ 2 - 1 ≠ 0 := by nlinarith
    simp only [PhiEta, sinh_eq, cosh_eq, hneg]
    field_simp
  have hhalf : Tendsto (fun t : ℝ => t / 2) atTop atTop :=
    tendsto_id.atTop_div_const (by norm_num)
  refine (hE.comp hhalf).congr' ?_
  filter_upwards [eventually_gt_atTop 0] with y hy
  exact (Phi_eq_PhiEta hy).symm

/-! ## The zero `τ∞` and its rational bracket -/

theorem tanh_half_eq (τ : ℝ) : tanh (τ / 2) = (exp τ - 1) / (exp τ + 1) := by
  have hw : 0 < exp (τ / 2) := exp_pos _
  have h1 : exp τ = exp (τ / 2) ^ 2 := by rw [← exp_nat_mul]; norm_num; ring_nf
  have h2 : exp (-(τ / 2)) = (exp (τ / 2))⁻¹ := exp_neg _
  rw [tanh_eq_sinh_div_cosh, sinh_eq, cosh_eq, h2, h1]
  field_simp

/-- `8 Φ(τ) = 2 t/(τ/2) − (1/t − 1/(τ/2))²` with `t = tanh(τ/2) = (e^τ − 1)/(e^τ + 1)`. -/
theorem eight_mul_Phi_t {τ : ℝ} (hτ : 0 < τ) :
    8 * Phi τ = 2 * ((exp τ - 1) / (exp τ + 1)) / (τ / 2)
      - (1 / ((exp τ - 1) / (exp τ + 1)) - 1 / (τ / 2)) ^ 2 := by
  rw [eight_mul_Phi hτ, ← tanh_half_eq, langevin, coth, tanh_eq_sinh_div_cosh]
  simp only [one_div_div]

/-- The function `t ↦ 2t/η − (1/t − 1/η)²` is monotone on `(0,1)` for `η ≥ 1`. -/
theorem g_mono {η t₀ t : ℝ} (hη : 1 ≤ η) (h0 : 0 < t₀) (h01 : t₀ ≤ t) (ht1 : t < 1) :
    2 * t₀ / η - (1 / t₀ - 1 / η) ^ 2 ≤ 2 * t / η - (1 / t - 1 / η) ^ 2 := by
  have hη0 : 0 < η := by linarith
  have ht : 0 < t := lt_of_lt_of_le h0 h01
  have a : 2 * t₀ / η ≤ 2 * t / η := by gcongr
  have b1 : 1 / η ≤ 1 := by rw [div_le_one hη0]; exact hη
  have b2 : 1 < 1 / t := by rw [lt_div_iff₀ ht]; linarith
  have b3 : 1 / t ≤ 1 / t₀ := one_div_le_one_div_of_le h0 h01
  have b : (1 / t - 1 / η) ^ 2 ≤ (1 / t₀ - 1 / η) ^ 2 := by
    apply pow_le_pow_left₀ (by linarith); linarith
  linarith

theorem exp_lower : (1711 : ℝ) ≤ exp (74453 / 10000) := by
  have h1 : (25361 / 10000 : ℝ) ≤ exp (74453 / 80000) := by
    refine le_trans ?_ (Real.sum_le_exp_of_nonneg (by norm_num) 9)
    norm_num [Finset.sum_range_succ, Nat.factorial]
  have h2 : exp (74453 / 10000 : ℝ) = exp (74453 / 80000) ^ 8 := by
    rw [← exp_nat_mul]; norm_num
  rw [h2]
  calc (1711 : ℝ) ≤ (25361 / 10000) ^ 8 := by norm_num
    _ ≤ _ := pow_le_pow_left₀ (by norm_num) h1 8

theorem exp_upper : exp (74454 / 10000) ≤ (1713 : ℝ) := by
  have h1 : exp (74454 / 80000) ≤ (25363 / 10000 : ℝ) := by
    refine le_trans (Real.exp_bound' (by norm_num) (by norm_num) (n := 9) (by norm_num)) ?_
    norm_num [Finset.sum_range_succ, Nat.factorial]
  have h2 : exp (74454 / 10000 : ℝ) = exp (74454 / 80000) ^ 8 := by
    rw [← exp_nat_mul]; norm_num
  rw [h2]
  calc exp (74454 / 80000 : ℝ) ^ 8 ≤ (25363 / 10000) ^ 8 :=
        pow_le_pow_left₀ (exp_pos _).le h1 8
    _ ≤ 1713 := by norm_num

/-- `Φ(7.4453) > 0`. -/
theorem Phi_lower_pos : 0 < Phi (74453 / 10000) := by
  have hE := exp_lower
  have h8 := eight_mul_Phi_t (τ := 74453 / 10000) (by norm_num)
  set E := exp (74453 / 10000 : ℝ)
  have ht0 : (1710 / 1712 : ℝ) ≤ (E - 1) / (E + 1) := by
    rw [div_le_div_iff₀ (by norm_num) (by linarith)]; linarith
  have ht1 : (E - 1) / (E + 1) < 1 := by
    rw [div_lt_one (by linarith)]; linarith
  have hg := g_mono (η := 74453 / 10000 / 2) (by norm_num) (by norm_num) ht0 ht1
  have hnum : (0 : ℝ) < 2 * (1710 / 1712) / (74453 / 10000 / 2)
      - (1 / (1710 / 1712) - 1 / (74453 / 10000 / 2)) ^ 2 := by norm_num
  linarith

/-- `Φ(7.4454) < 0`. -/
theorem Phi_upper_neg : Phi (74454 / 10000) < 0 := by
  have hE := exp_upper
  have hE1 : 1 < exp (74454 / 10000 : ℝ) := one_lt_exp_iff.mpr (by norm_num)
  have h8 := eight_mul_Phi_t (τ := 74454 / 10000) (by norm_num)
  set E := exp (74454 / 10000 : ℝ)
  have ht0 : 0 < (E - 1) / (E + 1) := div_pos (by linarith) (by linarith)
  have ht1 : (E - 1) / (E + 1) ≤ 1712 / 1714 := by
    rw [div_le_div_iff₀ (by linarith) (by norm_num)]; linarith
  have hg := g_mono (η := 74454 / 10000 / 2) (by norm_num) ht0 ht1 (by norm_num)
  have hnum : 2 * (1712 / 1714) / (74454 / 10000 / 2)
      - (1 / (1712 / 1714) - 1 / (74454 / 10000 / 2)) ^ 2 < (0 : ℝ) := by norm_num
  linarith

/-- note "The limit": `Φ` has a zero in `(7.4453, 7.4454)`. -/
theorem exists_Phi_zero_bracket :
    ∃ τ : ℝ, 74453 / 10000 < τ ∧ τ < 74454 / 10000 ∧ Phi τ = 0 := by
  have hcont : ContinuousOn Phi (Icc (74453 / 10000) (74454 / 10000)) :=
    continuousOn_Phi.mono (fun t ht => lt_of_lt_of_le (by norm_num) ht.1)
  have h := intermediate_value_Icc' (by norm_num) hcont
  obtain ⟨τ, ⟨h1, h2⟩, hτ⟩ :=
    h ⟨Phi_upper_neg.le, Phi_lower_pos.le⟩
  refine ⟨τ, lt_of_le_of_ne h1 ?_, lt_of_le_of_ne h2 ?_, hτ⟩
  · rintro rfl; exact Phi_lower_pos.ne' hτ
  · rintro rfl; exact Phi_upper_neg.ne hτ

/-- `τ∞`, the unique zero of `Φ` on `(0,∞)` (note eq. (kappa): `τ∞ = 2 η∞`). -/
def tauInf : ℝ := Classical.choose exists_Phi_zero_bracket

/-- `κ∞ = 1/τ∞` (note eq. (kappa)). -/
def kappaInf : ℝ := 1 / tauInf

/-- `7.4453 < τ∞ < 7.4454`. -/
theorem tauInf_bracket : 74453 / 10000 < tauInf ∧ tauInf < 74454 / 10000 :=
  ⟨(Classical.choose_spec exists_Phi_zero_bracket).1,
    (Classical.choose_spec exists_Phi_zero_bracket).2.1⟩

theorem tauInf_pos : 0 < tauInf := lt_trans (by norm_num) tauInf_bracket.1

theorem Phi_tauInf : Phi tauInf = 0 := (Classical.choose_spec exists_Phi_zero_bracket).2.2

/-- note "The limit": `τ∞` is the unique zero of `Φ` on `(0,∞)`. -/
theorem Phi_eq_zero_iff {τ : ℝ} (hτ : 0 < τ) : Phi τ = 0 ↔ τ = tauInf := by
  constructor
  · intro h
    by_contra hne
    rcases lt_or_gt_of_ne hne with hlt | hlt
    · have := Phi_strictAntiOn (mem_Ioi.mpr hτ) (mem_Ioi.mpr tauInf_pos) hlt
      rw [h, Phi_tauInf] at this; exact lt_irrefl _ this
    · have := Phi_strictAntiOn (mem_Ioi.mpr tauInf_pos) (mem_Ioi.mpr hτ) hlt
      rw [h, Phi_tauInf] at this; exact lt_irrefl _ this
  · rintro rfl; exact Phi_tauInf

/-- note "The limit": existence and uniqueness of `τ∞ > 0` with `Φ(τ∞) = 0`. -/
theorem existsUnique_Phi_zero : ∃! τ : ℝ, 0 < τ ∧ Phi τ = 0 :=
  ⟨tauInf, ⟨tauInf_pos, Phi_tauInf⟩, fun _ hτ => (Phi_eq_zero_iff hτ.1).mp hτ.2⟩

/-- `1/7.4454 < κ∞ < 1/7.4453`. -/
theorem kappaInf_bracket : 10000 / 74454 < kappaInf ∧ kappaInf < 10000 / 74453 := by
  obtain ⟨h1, h2⟩ := tauInf_bracket
  have hp := tauInf_pos
  constructor
  · rw [kappaInf, lt_div_iff₀ hp]
    rw [div_mul_eq_mul_div, div_lt_one (by norm_num)]; linarith
  · rw [kappaInf, div_lt_div_iff₀ hp (by norm_num)]; linarith

end Havelock.Tanh.LargeN
