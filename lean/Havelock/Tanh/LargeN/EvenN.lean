import Havelock.Tanh.LargeN.Stiffness
import Havelock.Tanh.LargeN.PhiLimit

/-! # Tanh model, large `N`: even `N` has one threshold (note Prop 15)

For even `N`, with `η = τ/2` and `μ = (2/N) coth(2η/N) = (2/N) coth(τ/N)`,
`8 a_{N/2}/N² = 2 μ tanh η − (coth η − μ)² − 4/N²` (`eight_mul_aBlock_half`). For `N = 2` this is
identically `1` (`evenRHS_two`). For `N ≥ 3` (in particular even `N ≥ 4`) it is strictly decreasing in
`τ > 0` (`evenRHS_strictAntiOn`): following the note, `μ tanh η` is strictly decreasing because
`t/sinh t` is strictly decreasing, and `coth η − μ` is non-negative (`t coth t` is increasing) and
non-decreasing (`t/sinh t` decreasing again); here this is done on the derivative `evenRHSDeriv`.
Limits: `2 − 4/N²` as `τ → 0⁺` (`r → ∞`) and `−(N² − 8N + 8)/N²` as `τ → ∞` (`r → 0`).
Consequences for `a_{N/2}` as a function of `ξ ∈ (0,1)` (via `τ = −(N/2) ln ξ`, a decreasing bijection
`(0,1) → (0,∞)`): `a_{N/2} > 0` for even `N ≤ 6` (`aBlock_half_pos_of_le_six`), and a unique zero
`ξ₀ ∈ (0,1)` with `a_{N/2} < 0` on `(0,ξ₀)` and `> 0` on `(ξ₀,1)` for even `N ≥ 8`
(`aBlock_half_threshold`). -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real Set Filter Topology

/-- `μ = (2/N) coth(τ/N)` (note Prop 15, with `τ/N = 2η/N`). -/
def muEven (N : ℕ) (τ : ℝ) : ℝ := 2 / N * coth (τ / N)

/-- The right-hand side of note Prop 15, `2 μ tanh η − (coth η − μ)² − 4/N²` with `η = τ/2`. -/
def evenRHS (N : ℕ) (τ : ℝ) : ℝ :=
  2 * muEven N τ * tanh (τ / 2) - (coth (τ / 2) - muEven N τ) ^ 2 - 4 / (N : ℝ) ^ 2

theorem cast_half_of_even {N : ℕ} (hE : Even N) : ((N / 2 : ℕ) : ℝ) = (N : ℝ) / 2 := by
  obtain ⟨n, rfl⟩ := hE
  rw [show (n + n) / 2 = n by omega]
  push_cast; ring

theorem xOf_half {N : ℕ} (hE : Even N) (hN : 1 ≤ N) : xOf N (N / 2) = 0 := by
  have hN0 : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  unfold xOf
  rw [cast_half_of_even hE]
  field_simp
  ring

/-- note Prop 15: for even `N ≥ 2` and `τ > 0`,
`8 a_{N/2}(ξ)/N² = 2 μ tanh η − (coth η − μ)² − 4/N²` with `ξ = e^{−2τ/N}`, `η = τ/2`,
`μ = (2/N) coth(τ/N)`. -/
theorem eight_mul_aBlock_half {N : ℕ} (hE : Even N) (hN : 2 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    8 * aBlock N (N / 2) (xiOf N τ) / (N : ℝ) ^ 2 = evenRHS N τ := by
  have hm1 : 1 ≤ N / 2 := by omega
  have hmN : N / 2 < N := by omega
  have hN0 : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  rw [aBlock_xiOf hm1 hmN hτ, xOf_half hE (by omega)]
  unfold evenRHS muEven EN kN
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
  linear_combination 8 * (N : ℝ) ^ 2 * (2 * cosh η ^ 2 + 1) * hpy


/-- note Prop 15: for `N = 2` the right-hand side is identically `1`. -/
theorem evenRHS_two {τ : ℝ} (hτ : 0 < τ) : evenRHS 2 τ = 1 := by
  have hη : 0 < τ / 2 := by positivity
  have hs : sinh (τ / 2) ≠ 0 := (sinh_pos_iff.mpr hη).ne'
  have hc : cosh (τ / 2) ≠ 0 := (cosh_pos _).ne'
  simp only [evenRHS, muEven, coth, tanh_eq_sinh_div_cosh, Nat.cast_ofNat]
  field_simp
  ring

/-! ## Derivative -/

theorem hasDerivAt_coth_div {a τ : ℝ} (ha : 0 < a) (hτ : 0 < τ) :
    HasDerivAt (fun t => coth (t / a)) (-1 / sinh (τ / a) ^ 2 / a) τ := by
  have h1 : HasDerivAt (fun t : ℝ => t / a) (1 / a) τ := by
    simpa using (hasDerivAt_id τ).div_const a
  have h := (hasDerivAt_cosh_div_sinh (η := τ / a) (by positivity)).comp τ h1
  convert h using 1
  · rfl
  · ring

theorem hasDerivAt_tanh_div (a τ : ℝ) :
    HasDerivAt (fun t => tanh (t / a)) (1 / cosh (τ / a) ^ 2 / a) τ := by
  have h1 : HasDerivAt (fun t : ℝ => t / a) (1 / a) τ := by
    simpa using (hasDerivAt_id τ).div_const a
  have h := (hasDerivAt_sinh_div_cosh (τ / a)).comp τ h1
  have hf : (fun t => tanh (t / a)) = (fun y => sinh y / cosh y) ∘ (fun t : ℝ => t / a) := by
    funext t; simp [tanh_eq_sinh_div_cosh]
  rw [hf]
  convert h using 1
  ring

def evenRHSDeriv (N : ℕ) (τ : ℝ) : ℝ :=
  2 * ((2 / N * (-1 / sinh (τ / N) ^ 2 / N)) * tanh (τ / 2)
      + muEven N τ * (1 / cosh (τ / 2) ^ 2 / 2))
    - 2 * (coth (τ / 2) - muEven N τ)
      * (-1 / sinh (τ / 2) ^ 2 / 2 - 2 / N * (-1 / sinh (τ / N) ^ 2 / N))

theorem hasDerivAt_evenRHS {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    HasDerivAt (evenRHS N) (evenRHSDeriv N τ) τ := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hm : HasDerivAt (muEven N) (2 / N * (-1 / sinh (τ / N) ^ 2 / N)) τ :=
    (hasDerivAt_coth_div hNpos hτ).const_mul (2 / (N : ℝ))
  have ht := hasDerivAt_tanh_div 2 τ
  have hk := hasDerivAt_coth_div (a := 2) two_pos hτ
  have h := (((hm.mul ht).const_mul 2).sub ((hk.sub hm).pow 2)).sub_const (4 / (N : ℝ) ^ 2)
  convert h using 1
  · funext x; simp only [evenRHS, Pi.sub_apply, Pi.mul_apply, Pi.pow_apply]; ring
  · simp only [evenRHSDeriv, Nat.cast_ofNat, Pi.sub_apply]; ring


/-! ## Sign of the derivative -/

/-- `t ↦ t coth t` is monotone on `(0,∞)`, cleared of denominators. -/
theorem mul_coth_le_evenN {a b : ℝ} (ha : 0 < a) (hab : a ≤ b) :
    a * cosh a * sinh b ≤ b * cosh b * sinh a := by
  have hb : 0 < b := lt_of_lt_of_le ha hab
  have hd : ∀ t : ℝ, 0 < t → HasDerivAt (fun y => y * (cosh y / sinh y))
      (1 * (cosh t / sinh t) + t * (-1 / sinh t ^ 2)) t := fun t ht =>
    (hasDerivAt_id t).mul (hasDerivAt_cosh_div_sinh ht.ne')
  have hmono : StrictMonoOn (fun y => y * (cosh y / sinh y)) (Ioi 0) := by
    refine strictMonoOn_of_deriv_pos (convex_Ioi 0)
      (fun y hy => (hd y hy).continuousAt.continuousWithinAt) ?_
    intro t ht
    rw [interior_Ioi] at ht
    have ht' : 0 < t := ht
    rw [(hd t ht').deriv]
    have hs : 0 < sinh t := sinh_pos_iff.mpr ht'
    have hst : t < sinh t := self_lt_sinh_iff.mpr ht'
    have hc : 1 ≤ cosh t := one_le_cosh t
    have hnum : 0 < cosh t * sinh t - t := by nlinarith
    have : 1 * (cosh t / sinh t) + t * (-1 / sinh t ^ 2) = (cosh t * sinh t - t) / sinh t ^ 2 := by
      field_simp
      ring
    rw [this]; positivity
  have h := hmono.monotoneOn (mem_Ioi.mpr ha) (mem_Ioi.mpr hb) hab
  have hsa : 0 < sinh a := sinh_pos_iff.mpr ha
  have hsb : 0 < sinh b := sinh_pos_iff.mpr hb
  rw [mul_div_assoc', mul_div_assoc', div_le_div_iff₀ hsa hsb] at h
  linarith

theorem evenRHSDeriv_neg {N : ℕ} (hN : 3 ≤ N) {τ : ℝ} (hτ : 0 < τ) : evenRHSDeriv N τ < 0 := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  have hN3 : (3 : ℝ) ≤ N := by exact_mod_cast hN
  set a := τ / N with ha_def
  set b := τ / 2 with hb_def
  have ha : 0 < a := by positivity
  have hb : 0 < b := by positivity
  have hτa : τ = N * a := by rw [ha_def]; field_simp
  have hτb : τ = 2 * b := by rw [hb_def]; ring
  have hab : a < b := by nlinarith
  have hsa : 0 < sinh a := sinh_pos_iff.mpr ha
  have hsb : 0 < sinh b := sinh_pos_iff.mpr hb
  have hca : 0 < cosh a := cosh_pos a
  have hcb : 0 < cosh b := cosh_pos b
  -- key 1: `N sinh a cosh a < 2 sinh b cosh b` (`t/sinh t` strictly decreasing)
  have k1 : N * sinh a * cosh a < 2 * sinh b * cosh b := by
    have h := sinh_div_strict (u := 2 * a) (v := τ - 2 * a) (by positivity) (by nlinarith)
    rw [show 2 * a + (τ - 2 * a) = τ by ring, sinh_two_mul] at h
    have hsτ : sinh τ = 2 * sinh b * cosh b := by rw [hτb, sinh_two_mul]
    rw [hsτ] at h
    have h2 : 2 * a * (N * sinh a * cosh a) < 2 * a * (2 * sinh b * cosh b) := by
      rw [hτa] at h; linarith
    exact lt_of_mul_lt_mul_left h2 (by positivity)
  -- key 2: `N sinh a ≤ 2 sinh b`
  have k2 : N * sinh a ≤ 2 * sinh b := by
    have h := sinh_div_strict (u := a) (v := b - a) ha (by linarith)
    rw [show a + (b - a) = b by ring] at h
    have h2 : a * (N * sinh a) < a * (2 * sinh b) := by
      have : 2 * b = N * a := by rw [← hτa, hτb]
      nlinarith
    exact (lt_of_mul_lt_mul_left h2 ha.le).le
  -- key 3: `2 sinh b cosh a ≤ N sinh a cosh b` (`t coth t` monotone)
  have k3 : 2 * sinh b * cosh a ≤ N * sinh a * cosh b := by
    have h := mul_coth_le_evenN ha hab.le
    have h2 : a * (2 * sinh b * cosh a) ≤ a * (N * sinh a * cosh b) := by
      have : 2 * b = N * a := by rw [← hτa, hτb]
      have e : a * (N * sinh a * cosh b) = 2 * (b * cosh b * sinh a) := by
        rw [show a * (N * sinh a * cosh b) = (N * a) * sinh a * cosh b by ring, ← this]; ring
      rw [e]; linarith
    exact le_of_mul_le_mul_left h2 ha
  have hT1 : 2 / N * (-1 / sinh a ^ 2 / N) * tanh b + muEven N τ * (1 / cosh b ^ 2 / 2)
      = (N * sinh a * cosh a - 2 * sinh b * cosh b) / (N ^ 2 * sinh a ^ 2 * cosh b ^ 2) := by
    simp only [muEven, coth, tanh_eq_sinh_div_cosh, ← ha_def]
    field_simp
    ring
  have hD : coth b - muEven N τ = (N * sinh a * cosh b - 2 * sinh b * cosh a) / (N * sinh a * sinh b) := by
    simp only [muEven, coth, ← ha_def]
    field_simp
  have hD' : -1 / sinh b ^ 2 / 2 - 2 / N * (-1 / sinh a ^ 2 / N)
      = (4 * sinh b ^ 2 - N ^ 2 * sinh a ^ 2) / (2 * N ^ 2 * sinh a ^ 2 * sinh b ^ 2) := by
    field_simp
    ring
  have e1 : (N * sinh a * cosh a - 2 * sinh b * cosh b) / (N ^ 2 * sinh a ^ 2 * cosh b ^ 2) < 0 :=
    div_neg_of_neg_of_pos (by linarith) (by positivity)
  have e2 : 0 ≤ (N * sinh a * cosh b - 2 * sinh b * cosh a) / (N * sinh a * sinh b) :=
    div_nonneg (by linarith) (by positivity)
  have e3 : 0 ≤ (4 * sinh b ^ 2 - N ^ 2 * sinh a ^ 2) / (2 * N ^ 2 * sinh a ^ 2 * sinh b ^ 2) := by
    apply div_nonneg _ (by positivity)
    have := mul_le_mul k2 k2 (by positivity) (by positivity)
    nlinarith
  simp only [evenRHSDeriv, ← ha_def, ← hb_def]
  rw [hT1, hD, hD']
  nlinarith [mul_nonneg e2 e3]

theorem continuousOn_evenRHS {N : ℕ} (hN : 1 ≤ N) : ContinuousOn (evenRHS N) (Ioi 0) :=
  fun _ ht => (hasDerivAt_evenRHS hN ht).continuousAt.continuousWithinAt

/-- note Prop 15: for `N ≥ 3` the right-hand side is strictly decreasing in `τ` on `(0,∞)`. -/
theorem evenRHS_strictAntiOn {N : ℕ} (hN : 3 ≤ N) : StrictAntiOn (evenRHS N) (Ioi 0) := by
  refine strictAntiOn_of_deriv_neg (convex_Ioi 0) (continuousOn_evenRHS (by omega)) ?_
  intro t ht
  rw [interior_Ioi] at ht
  rw [(hasDerivAt_evenRHS (by omega) ht).deriv]
  exact evenRHSDeriv_neg hN ht


/-- note Prop 15: for even `N ≥ 4` the right-hand side is strictly decreasing in `τ` on `(0,∞)`. -/
theorem evenRHS_strictAntiOn_of_even {N : ℕ} (_hE : Even N) (hN : 4 ≤ N) :
    StrictAntiOn (evenRHS N) (Ioi 0) :=
  evenRHS_strictAntiOn (by omega)

/-! ## Limits -/

theorem tendsto_coth_atTop : Tendsto coth atTop (𝓝 1) := by
  have hu : Tendsto (fun η : ℝ => exp (-η)) atTop (𝓝 0) := Real.tendsto_exp_neg_atTop_nhds_zero
  have hu2 : Tendsto (fun η : ℝ => exp (-η) ^ 2) atTop (𝓝 0) := by simpa using hu.pow 2
  have hcoth : Tendsto (fun η : ℝ => (1 + exp (-η) ^ 2) / (1 - exp (-η) ^ 2)) atTop (𝓝 1) := by
    have := (tendsto_const_nhds (x := (1 : ℝ))).add hu2 |>.div
      ((tendsto_const_nhds (x := (1 : ℝ))).sub hu2) (by norm_num)
    simp only [sub_zero, add_zero, div_one] at this
    exact this
  refine hcoth.congr' ?_
  filter_upwards [eventually_gt_atTop 0] with y hy
  have hw : 0 < exp y := exp_pos y
  have hw1 : 1 < exp y := one_lt_exp_iff.mpr hy
  have hneg : exp (-y) = (exp y)⁻¹ := exp_neg y
  have hd : 1 - ((exp y)⁻¹) ^ 2 ≠ 0 := by
    have : (exp y)⁻¹ < 1 := inv_lt_one_of_one_lt₀ hw1
    have : 0 < (exp y)⁻¹ := inv_pos.mpr hw
    nlinarith
  have hd' : exp y ^ 2 - 1 ≠ 0 := by nlinarith
  have hd'' : exp y - (exp y)⁻¹ ≠ 0 := by
    have : (exp y)⁻¹ < 1 := inv_lt_one_of_one_lt₀ hw1
    linarith
  simp only [coth, sinh_eq, cosh_eq, hneg]
  field_simp

theorem tendsto_tanh_atTop : Tendsto tanh atTop (𝓝 1) := by
  have h := tendsto_coth_atTop.inv₀ one_ne_zero
  rw [inv_one] at h
  refine h.congr' (Eventually.of_forall fun y => ?_)
  simp only [coth, tanh_eq_sinh_div_cosh, inv_div]

/-- note Prop 15: as `τ → ∞` (`r → 0`) the right-hand side tends to `−(N² − 8N + 8)/N²`. -/
theorem tendsto_evenRHS_atTop {N : ℕ} (hN : 1 ≤ N) :
    Tendsto (evenRHS N) atTop (𝓝 (-((N : ℝ) ^ 2 - 8 * N + 8) / (N : ℝ) ^ 2)) := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hcN : Tendsto (fun τ : ℝ => coth (τ / N)) atTop (𝓝 1) :=
    tendsto_coth_atTop.comp (tendsto_id.atTop_div_const hNpos)
  have hc2 : Tendsto (fun τ : ℝ => coth (τ / 2)) atTop (𝓝 1) :=
    tendsto_coth_atTop.comp (tendsto_id.atTop_div_const two_pos)
  have ht2 : Tendsto (fun τ : ℝ => tanh (τ / 2)) atTop (𝓝 1) :=
    tendsto_tanh_atTop.comp (tendsto_id.atTop_div_const two_pos)
  have hmu : Tendsto (muEven N) atTop (𝓝 (2 / N * 1)) := hcN.const_mul _
  have h := (((hmu.const_mul 2).mul ht2).sub ((hc2.sub hmu).pow 2)).sub_const (4 / (N : ℝ) ^ 2)
  have hv : 2 * (2 / (N : ℝ) * 1) * 1 - (1 - 2 / N * 1) ^ 2 - 4 / (N : ℝ) ^ 2
      = -((N : ℝ) ^ 2 - 8 * N + 8) / (N : ℝ) ^ 2 := by
    field_simp; ring
  rw [hv] at h
  exact h

theorem tendsto_langevin_zero : Tendsto langevin (𝓝[>] 0) (𝓝 0) := by
  have hid : Tendsto (fun η : ℝ => η) (𝓝[>] 0) (𝓝 0) := nhdsWithin_le_nhds
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hid ?_ ?_
  · filter_upwards [self_mem_nhdsWithin] with y hy
    exact (langevin_pos hy).le
  · filter_upwards [self_mem_nhdsWithin] with y hy
    exact langevin_le_self hy

theorem tendsto_tanh_div_self_zero : Tendsto (fun y : ℝ => tanh y / y) (𝓝[>] 0) (𝓝 1) := by
  have h := hasDerivAt_iff_tendsto_slope.mp (hasDerivAt_sinh_div_cosh 0)
  rw [show (1 : ℝ) / cosh 0 ^ 2 = 1 by simp] at h
  refine (h.mono_left (nhdsWithin_mono _ ?_)).congr' ?_
  · intro y hy; exact ne_of_gt hy
  · filter_upwards [self_mem_nhdsWithin] with y hy
    simp [slope_def_field, tanh_eq_sinh_div_cosh]

theorem tendsto_tanh_zero : Tendsto tanh (𝓝[>] 0) (𝓝 0) := by
  have h := tendsto_tanh_div_self_zero.mul (nhdsWithin_le_nhds (a := (0 : ℝ)) (s := Ioi 0)
    |> tendsto_id.mono_left)
  rw [mul_zero] at h
  refine h.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with y hy
  have : y ≠ 0 := ne_of_gt hy
  simp only [id]
  field_simp

theorem tendsto_div_nhdsGT {a : ℝ} (ha : 0 < a) :
    Tendsto (fun t : ℝ => t / a) (𝓝[>] 0) (𝓝[>] 0) := by
  refine tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within _ ?_ ?_
  · have := ((continuous_id.div_const a).tendsto 0).mono_left (nhdsWithin_le_nhds (s := Ioi 0))
    simpa using this
  · filter_upwards [self_mem_nhdsWithin] with y hy
    show (0 : ℝ) < y / a
    exact div_pos hy ha

/-- For `τ > 0`: `evenRHS` in terms of `𝓛`: `μ tanh η = (2/N) 𝓛(τ/N) tanh η + tanh η / η` and
`coth η − μ = 𝓛(η) − (2/N) 𝓛(τ/N)`. -/
theorem evenRHS_eq_langevin {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    evenRHS N τ = 2 * (2 / N * langevin (τ / N) * tanh (τ / 2) + tanh (τ / 2) / (τ / 2))
      - (langevin (τ / 2) - 2 / N * langevin (τ / N)) ^ 2 - 4 / (N : ℝ) ^ 2 := by
  have hN0 : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have hτ0 : τ ≠ 0 := hτ.ne'
  simp only [evenRHS, muEven, langevin]
  field_simp
  ring

/-- note Prop 15: as `τ → 0⁺` (`r → ∞`) the right-hand side tends to `2 − 4/N²`. -/
theorem tendsto_evenRHS_zero {N : ℕ} (hN : 1 ≤ N) :
    Tendsto (evenRHS N) (𝓝[>] 0) (𝓝 (2 - 4 / (N : ℝ) ^ 2)) := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hLN : Tendsto (fun τ : ℝ => langevin (τ / N)) (𝓝[>] 0) (𝓝 0) :=
    tendsto_langevin_zero.comp (tendsto_div_nhdsGT hNpos)
  have hL2 : Tendsto (fun τ : ℝ => langevin (τ / 2)) (𝓝[>] 0) (𝓝 0) :=
    tendsto_langevin_zero.comp (tendsto_div_nhdsGT two_pos)
  have hT2 : Tendsto (fun τ : ℝ => tanh (τ / 2)) (𝓝[>] 0) (𝓝 0) :=
    tendsto_tanh_zero.comp (tendsto_div_nhdsGT two_pos)
  have hQ2 : Tendsto (fun τ : ℝ => tanh (τ / 2) / (τ / 2)) (𝓝[>] 0) (𝓝 1) :=
    tendsto_tanh_div_self_zero.comp (tendsto_div_nhdsGT two_pos)
  have h := ((((hLN.const_mul (2 / (N : ℝ))).mul hT2).add hQ2).const_mul 2).sub
    ((hL2.sub (hLN.const_mul (2 / (N : ℝ)))).pow 2) |>.sub_const (4 / (N : ℝ) ^ 2)
  have hv : (2 : ℝ) * (2 / N * 0 * 0 + 1) - (0 - 2 / N * 0) ^ 2 - 4 / (N : ℝ) ^ 2
      = 2 - 4 / (N : ℝ) ^ 2 := by ring
  rw [hv] at h
  refine h.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with y hy
  exact (evenRHS_eq_langevin hN hy).symm


/-! ## Consequences for the block `a_{N/2}` as a function of `ξ` -/

/-- `τ = −(N/2) ln ξ`, the inverse of `xiOf N`. -/
def tauOf (N : ℕ) (ξ : ℝ) : ℝ := -((N : ℝ) / 2) * log ξ

theorem tauOf_pos {N : ℕ} (hN : 1 ≤ N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) : 0 < tauOf N ξ := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hl : log ξ < 0 := log_neg h0 h1
  unfold tauOf
  nlinarith

theorem xiOf_tauOf {N : ℕ} (hN : 1 ≤ N) {ξ : ℝ} (h0 : 0 < ξ) : xiOf N (tauOf N ξ) = ξ := by
  have hN0 : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  unfold xiOf tauOf
  rw [show -2 * (-((N : ℝ) / 2) * log ξ) / N = log ξ by field_simp, exp_log h0]

theorem tauOf_xiOf {N : ℕ} (hN : 1 ≤ N) (τ : ℝ) : tauOf N (xiOf N τ) = τ := by
  have hN0 : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  unfold xiOf tauOf
  rw [log_exp]
  field_simp

theorem tauOf_lt_tauOf {N : ℕ} (hN : 1 ≤ N) {ξ₁ ξ₂ : ℝ} (h0 : 0 < ξ₁) (h : ξ₁ < ξ₂) :
    tauOf N ξ₂ < tauOf N ξ₁ := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hl := log_lt_log h0 h
  unfold tauOf
  nlinarith

/-- `a_{N/2}(ξ) = (N²/8) · evenRHS N τ` with `τ = −(N/2) ln ξ` (note Prop 15). -/
theorem aBlock_half_eq {N : ℕ} (hE : Even N) (hN : 2 ≤ N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) :
    aBlock N (N / 2) ξ = (N : ℝ) ^ 2 / 8 * evenRHS N (tauOf N ξ) := by
  have hN0 : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have h := eight_mul_aBlock_half hE hN (tauOf_pos (N := N) (by omega) h0 h1)
  rw [xiOf_tauOf (by omega) h0] at h
  rw [← h]
  field_simp

/-- For `2 ≤ N ≤ 6` the right-hand side is positive for every `τ > 0`. -/
theorem evenRHS_pos_of_le_six {N : ℕ} (hN2 : 2 ≤ N) (hN6 : N ≤ 6) {τ : ℝ} (hτ : 0 < τ) :
    0 < evenRHS N τ := by
  rcases (show N = 2 ∨ 3 ≤ N by omega) with rfl | h3
  · rw [evenRHS_two hτ]; norm_num
  · have hL := tendsto_evenRHS_atTop (N := N) (by omega)
    have h3' : (3 : ℝ) ≤ N := by exact_mod_cast h3
    have h6' : (N : ℝ) ≤ 6 := by exact_mod_cast hN6
    have hLpos : 0 < -((N : ℝ) ^ 2 - 8 * N + 8) / (N : ℝ) ^ 2 := by
      apply div_pos _ (by positivity)
      nlinarith
    have hle : -((N : ℝ) ^ 2 - 8 * N + 8) / (N : ℝ) ^ 2 ≤ evenRHS N (τ + 1) := by
      refine le_of_tendsto hL ?_
      filter_upwards [eventually_gt_atTop (τ + 1)] with t ht
      exact (evenRHS_strictAntiOn h3 (mem_Ioi.mpr (by linarith)) (mem_Ioi.mpr (by linarith)) ht).le
    have := evenRHS_strictAntiOn h3 (mem_Ioi.mpr hτ) (mem_Ioi.mpr (by linarith)) (lt_add_one τ)
    linarith

/-- note Prop 15: every even `N ≤ 6` is stable at every radius: `a_{N/2}(ξ) > 0` on `(0,1)`. -/
theorem aBlock_half_pos_of_le_six {N : ℕ} (hE : Even N) (hN2 : 2 ≤ N) (hN6 : N ≤ 6) {ξ : ℝ}
    (h0 : 0 < ξ) (h1 : ξ < 1) : 0 < aBlock N (N / 2) ξ := by
  have hN0 : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  rw [aBlock_half_eq hE hN2 h0 h1]
  exact mul_pos (by positivity) (evenRHS_pos_of_le_six hN2 hN6 (tauOf_pos (by omega) h0 h1))

/-- For `N ≥ 8` the right-hand side has a zero `τ₀ > 0` (intermediate value theorem). -/
theorem evenRHS_exists_zero {N : ℕ} (hN : 8 ≤ N) :
    ∃ τ₀ : ℝ, 0 < τ₀ ∧ evenRHS N τ₀ = 0 := by
  have h1N : 1 ≤ N := by omega
  have hN' : (8 : ℝ) ≤ N := by exact_mod_cast hN
  have hc0 : (0 : ℝ) < 2 - 4 / (N : ℝ) ^ 2 := by
    have : 4 / (N : ℝ) ^ 2 ≤ 4 / 64 := by
      apply div_le_div_of_nonneg_left (by norm_num) (by norm_num); nlinarith
    linarith
  have hL : -((N : ℝ) ^ 2 - 8 * N + 8) / (N : ℝ) ^ 2 < 0 :=
    div_neg_of_neg_of_pos (by nlinarith) (by positivity)
  obtain ⟨τ₁, hF1, hτ1⟩ :=
    (((tendsto_evenRHS_zero h1N).eventually (lt_mem_nhds hc0)).and self_mem_nhdsWithin).exists
  obtain ⟨τ₂, hF2, hτ12⟩ :=
    (((tendsto_evenRHS_atTop h1N).eventually (gt_mem_nhds hL)).and (eventually_gt_atTop τ₁)).exists
  have hτ1' : 0 < τ₁ := hτ1
  have hcont : ContinuousOn (evenRHS N) (Icc τ₁ τ₂) :=
    (continuousOn_evenRHS h1N).mono (fun t ht => mem_Ioi.mpr (lt_of_lt_of_le hτ1' ht.1))
  obtain ⟨τ₀, hτ0, hF0⟩ := intermediate_value_Icc' hτ12.le hcont ⟨hF2.le, hF1.le⟩
  exact ⟨τ₀, lt_of_lt_of_le hτ1' hτ0.1, hF0⟩

/-- note Prop 15: every even `N ≥ 8` has exactly one spectral threshold `ξ₀ ∈ (0,1)`, at a sign
change of `a_{N/2}`: `a_{N/2} < 0` on `(0, ξ₀)` and `a_{N/2} > 0` on `(ξ₀, 1)`. -/
theorem aBlock_half_threshold {N : ℕ} (hE : Even N) (hN : 8 ≤ N) :
    ∃ ξ₀ ∈ Ioo (0 : ℝ) 1, (∀ ξ ∈ Ioo (0 : ℝ) 1, aBlock N (N / 2) ξ = 0 ↔ ξ = ξ₀) ∧
      (∀ ξ ∈ Ioo 0 ξ₀, aBlock N (N / 2) ξ < 0) ∧ (∀ ξ ∈ Ioo ξ₀ 1, 0 < aBlock N (N / 2) ξ) := by
  have h1N : 1 ≤ N := by omega
  have h2N : 2 ≤ N := by omega
  have h3N : 3 ≤ N := by omega
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  obtain ⟨τ₀, hτ0, hF0⟩ := evenRHS_exists_zero hN
  have hξ0 : xiOf N τ₀ ∈ Ioo (0 : ℝ) 1 := ⟨xiOf_pos N τ₀, xiOf_lt_one h1N hτ0⟩
  have hneg : ∀ ξ ∈ Ioo 0 (xiOf N τ₀), aBlock N (N / 2) ξ < 0 := by
    intro ξ ⟨h0, h⟩
    have h1 : ξ < 1 := lt_trans h hξ0.2
    rw [aBlock_half_eq hE h2N h0 h1]
    have ht := tauOf_lt_tauOf h1N h0 h
    rw [tauOf_xiOf h1N] at ht
    have := evenRHS_strictAntiOn h3N (mem_Ioi.mpr hτ0) (mem_Ioi.mpr (tauOf_pos h1N h0 h1)) ht
    rw [hF0] at this
    exact mul_neg_of_pos_of_neg (by positivity) this
  have hpos : ∀ ξ ∈ Ioo (xiOf N τ₀) 1, 0 < aBlock N (N / 2) ξ := by
    intro ξ ⟨h, h1⟩
    have h0 : 0 < ξ := lt_trans hξ0.1 h
    rw [aBlock_half_eq hE h2N h0 h1]
    have ht := tauOf_lt_tauOf h1N hξ0.1 h
    rw [tauOf_xiOf h1N] at ht
    have := evenRHS_strictAntiOn h3N (mem_Ioi.mpr (tauOf_pos h1N h0 h1)) (mem_Ioi.mpr hτ0) ht
    rw [hF0] at this
    exact mul_pos (by positivity) this
  refine ⟨xiOf N τ₀, hξ0, ?_, hneg, hpos⟩
  intro ξ ⟨h0, h1⟩
  constructor
  · intro hz
    by_contra hne
    rcases lt_or_gt_of_ne hne with hlt | hgt
    · exact (hneg ξ ⟨h0, hlt⟩).ne hz
    · exact (hpos ξ ⟨hgt, h1⟩).ne' hz
  · rintro rfl
    rw [aBlock_half_eq hE h2N h0 h1, tauOf_xiOf h1N, hF0, mul_zero]

end Havelock.Tanh.LargeN
