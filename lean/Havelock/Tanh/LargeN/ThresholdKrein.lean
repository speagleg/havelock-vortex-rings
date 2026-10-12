import Havelock.Tanh.LargeN.ThresholdLimit

/-! # Tanh model, large `N`: the Krein window for odd `N` (note Thm 18(b))

For odd `N ≥ 9`: the spectral threshold `ξ_*(N)` lies in `(0,1)` (the planar value
`a_crit(0) = −(N−1)(N−7)/8 < 0`, and `Φ_N > 0` for small `τ`), `a_crit(ξ_*) = 0` by continuity,
and there `a_crit b_crit − g_crit² = −g_crit² < 0` because `R(1/N,τ) > 0`
(`∂_τ R = (1−x²) sinh τ sinh xτ > 0`, the derivative form of `tanh xτ > x tanh τ`). Hence the
critical block is indefinite on a right-neighbourhood of `ξ_*`, and `ξ_*(N) < ξ_EM(N)`.
The width bound `ξ_EM − ξ_* = O(N⁻³)` is `xiEM_sub_xiStar_le` in `ThresholdLimit`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real Set Filter Topology

/-! ### Continuity of the blocks on `[0,1)` -/

theorem continuousAt_aBlock {N m : ℕ} (hN : 1 ≤ N) {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) :
    ContinuousAt (fun ξ => aBlock N m ξ) ξ := by
  have hd1 : 1 - ξ ≠ 0 := by linarith
  have hdN : 1 - ξ ^ N ≠ 0 := by
    have := pow_lt_one₀ h0 h1 (by omega : N ≠ 0); linarith
  have hdN2 : (1 - ξ ^ N) ^ 2 ≠ 0 := pow_ne_zero 2 hdN
  unfold aBlock cR
  fun_prop (disch := assumption)

theorem continuousAt_bBlock {N m : ℕ} (hN : 1 ≤ N) {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) :
    ContinuousAt (fun ξ => bBlock N m ξ) ξ := by
  have hdN : 1 - ξ ^ N ≠ 0 := by
    have := pow_lt_one₀ h0 h1 (by omega : N ≠ 0); linarith
  have hdN2 : (1 - ξ ^ N) ^ 2 ≠ 0 := pow_ne_zero 2 hdN
  unfold bBlock W
  fun_prop (disch := assumption)

theorem continuousAt_gBlock {N m : ℕ} (hN : 1 ≤ N) {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) :
    ContinuousAt (fun ξ => gBlock N m ξ) ξ := by
  have hdN : 1 - ξ ^ N ≠ 0 := by
    have := pow_lt_one₀ h0 h1 (by omega : N ≠ 0); linarith
  have hdN2 : (1 - ξ ^ N) ^ 2 ≠ 0 := pow_ne_zero 2 hdN
  unfold gBlock W
  fun_prop (disch := assumption)

theorem continuousAt_det {N m : ℕ} (hN : 1 ≤ N) {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) :
    ContinuousAt (fun ξ => aBlock N m ξ * bBlock N m ξ - gBlock N m ξ ^ 2) ξ :=
  ((continuousAt_aBlock hN h0 h1).mul (continuousAt_bBlock hN h0 h1)).sub
    ((continuousAt_gBlock hN h0 h1).pow 2)

/-! ### The planar value and `R > 0` -/

/-- The planar (`ξ = 0`, i.e. `r → 0`) value `a_m(0) = (N−1) − m(N−m)/2`. -/
theorem aBlock_zero {N m : ℕ} (hm1 : 1 ≤ m) (hmN : m < N) :
    aBlock N m 0 = ((N : ℝ) - 1) - m * ((N : ℝ) - m) / 2 := by
  have e1 : (0 : ℝ) ^ N = 0 := zero_pow (by omega)
  have e2 : (0 : ℝ) ^ m = 0 := zero_pow (by omega)
  have e3 : (0 : ℝ) ^ (N - m) = 0 := zero_pow (by omega)
  have e4 : (0 : ℝ) ^ (N + m) = 0 := zero_pow (by omega)
  have e5 : (0 : ℝ) ^ (2 * N - m) = 0 := zero_pow (by omega)
  unfold aBlock cR
  rw [e1, e2, e3, e4, e5]
  ring

/-- note Thm 18(b): the planar limit of `a_crit` for odd `N`, `8 a_{(N−1)/2}(0) = −(N−1)(N−7)`. -/
theorem aBlock_half_zero_odd {N : ℕ} (hN : 2 ≤ N) (hodd : Odd N) :
    8 * aBlock N (N / 2) 0 = -(((N : ℝ) - 1) * ((N : ℝ) - 7)) := by
  rw [aBlock_zero (by omega) (by omega)]
  obtain ⟨k, rfl⟩ := hodd
  have hk : (2 * k + 1) / 2 = k := by omega
  rw [hk]
  push_cast
  ring

/-- For odd `N`, the critical mode `⌊N/2⌋` has `x = 1/N`. -/
theorem xOf_half_odd {N : ℕ} (hodd : Odd N) : xOf N (N / 2) = 1 / N := by
  obtain ⟨k, rfl⟩ := hodd
  have hk : (2 * k + 1) / 2 = k := by omega
  rw [hk]
  unfold xOf
  have : (0 : ℝ) < 2 * k + 1 := by positivity
  push_cast
  field_simp
  ring

theorem hasDerivAt_R (x τ : ℝ) :
    HasDerivAt (fun t => R x t) ((1 - x ^ 2) * sinh τ * sinh (x * τ)) τ := by
  have hC := ((hasDerivAt_id τ).const_mul x).cosh
  have hS := ((hasDerivAt_id τ).const_mul x).sinh
  have hc := Real.hasDerivAt_cosh τ
  have hs := Real.hasDerivAt_sinh τ
  have h := (hc.mul hS).sub ((hs.const_mul x).mul hC)
  convert h using 1
  · funext t
    simp only [R, id, Pi.sub_apply, Pi.mul_apply]
  · simp only [id]
    ring

/-- note Thm 18(b): `R(x,τ) > 0` for `0 < x < 1`, `τ > 0` (equivalently `tanh xτ > x tanh τ`). -/
theorem R_pos {x τ : ℝ} (hx0 : 0 < x) (hx1 : x < 1) (hτ : 0 < τ) : 0 < R x τ := by
  refine pos_of_hasDerivAt_pos (f := fun t => R x t)
    (f' := fun t => (1 - x ^ 2) * sinh t * sinh (x * t)) (fun t => hasDerivAt_R x t)
    (by simp [R]) (fun t ht => ?_) hτ
  have h1 : 0 < 1 - x ^ 2 := by nlinarith
  have h2 : 0 < sinh t := sinh_pos_iff.mpr ht
  have h3 : 0 < sinh (x * t) := sinh_pos_iff.mpr (mul_pos hx0 ht)
  positivity

/-- note Thm 18(b): for odd `N ≥ 2`, `g_crit ≠ 0` at every radius. -/
theorem gBlock_half_odd_pos {N : ℕ} (hN : 2 ≤ N) (hodd : Odd N) {ξ : ℝ} (h0 : 0 < ξ)
    (h1 : ξ < 1) : 0 < gBlock N (N / 2) ξ := by
  obtain ⟨τ, hτ, rfl⟩ := exists_tau (N := N) (by omega) h0 h1
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  have hN1 : (1 : ℝ) < N := by exact_mod_cast (show 1 < N by omega)
  rw [gBlock_xiOf (by omega) (by omega) hτ, xOf_half_odd hodd]
  unfold gamma
  have hx1 : 1 / (N : ℝ) < 1 := by rw [div_lt_one hNpos]; exact hN1
  have := R_pos (one_div_pos.mpr hNpos) hx1 hτ
  have hs : 0 < sinh τ := sinh_pos_iff.mpr hτ
  positivity

/-! ### note Thm 18(b) -/

/-- For `N ≥ 3`, the ring is positive definite (hence stable) near `ξ = 1`: some `ξ₁ < 1` has
`emDefinite` on `(ξ₁, 1)`. -/
theorem exists_emDefinite_near_one {N : ℕ} (hN : 3 ≤ N) :
    ∃ ξ₁ : ℝ, 0 < ξ₁ ∧ ξ₁ < 1 ∧ ∀ ξ ∈ Ioo ξ₁ 1, emDefinite N ξ := by
  have hev : ∀ᶠ τ in 𝓝[>] (0 : ℝ), 1 / 8 < Phi τ :=
    Filter.Tendsto.eventually_const_lt (by norm_num) tendsto_Phi_zero
  obtain ⟨τ₀, hΦ0, hτ0⟩ := (hev.and (Ioo_mem_nhdsGT (zero_lt_one' ℝ))).exists
  have hN1 : 1 ≤ N := by omega
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN1
  have hN3 : (3 : ℝ) ≤ N := by exact_mod_cast hN
  refine ⟨xiOf N τ₀, xiOf_pos N τ₀, xiOf_lt_one hN1 hτ0.1, fun ξ hξ => ?_⟩
  obtain ⟨τ, hτ, rfl⟩ := exists_tau hN1 ((xiOf_pos N _).trans hξ.1) hξ.2
  have hττ0 : τ < τ₀ := (xiOf_lt_iff hN1).mp hξ.1
  refine emDefinite_of_PhiN_pos (by omega) hτ ?_
  rw [PhiN_eq]
  have hmono : Phi τ₀ ≤ Phi τ :=
    Phi_strictAntiOn.antitoneOn (show (0 : ℝ) < τ from hτ) (show (0 : ℝ) < τ₀ from hτ0.1) hττ0.le
  have hΔ := (EN_sub_Einf_bounds hN1 hτ).1
  have hb : -(τ ^ 2 + 1) / (2 * (N : ℝ) ^ 2) ≥ -(1 / 9) := by
    rw [ge_iff_le, neg_div, neg_le_neg_iff, div_le_iff₀ (by positivity)]
    have : τ ^ 2 < 1 := by nlinarith [hτ0.2]
    nlinarith
  linarith

/-- note Thm 18(b): for every odd `N ≥ 9`,
* `0 < ξ_*(N) < 1` and `a_crit(ξ_*) = 0` (`crit = ⌊N/2⌋`);
* the ring is spectrally stable on `(ξ_*, 1)`;
* the critical block is indefinite (`a b − g² < 0`) on a right-neighbourhood of `ξ_*`;
* the Krein window is non-empty: `ξ_*(N) < ξ_EM(N)`. -/
theorem krein_window_odd {N : ℕ} (hN : 9 ≤ N) (hodd : Odd N) :
    0 < xiStar N ∧ xiStar N < 1 ∧ aBlock N (N / 2) (xiStar N) = 0 ∧
      (∀ ξ ∈ Ioo (xiStar N) 1, specStable N ξ) ∧
      (∃ η > 0, ∀ ξ ∈ Ioo (xiStar N) (xiStar N + η),
        aBlock N (N / 2) ξ * bBlock N (N / 2) ξ - gBlock N (N / 2) ξ ^ 2 < 0) ∧
      xiStar N < xiEM N := by
  have hN1 : 1 ≤ N := by omega
  have hN2 : 2 ≤ N := by omega
  have hNr : (9 : ℝ) ≤ N := by exact_mod_cast hN
  -- (1) `0 < ξ_*`: the planar value is negative.
  have ha0 : aBlock N (N / 2) 0 < 0 := by
    have := aBlock_half_zero_odd hN2 hodd
    nlinarith
  have hpos : 0 < xiStar N := by
    have hev : ∀ᶠ ξ in 𝓝 (0 : ℝ), aBlock N (N / 2) ξ < 0 :=
      Filter.Tendsto.eventually_lt_const ha0 (continuousAt_aBlock hN1 le_rfl zero_lt_one)
    obtain ⟨ε, hε, hball⟩ := Metric.eventually_nhds_iff.mp hev
    set ξ₁ := min (ε / 2) (1 / 2)
    have hξ₁0 : 0 < ξ₁ := lt_min (by positivity) (by norm_num)
    have hξ₁1 : ξ₁ < 1 := lt_of_le_of_lt (min_le_right _ _) (by norm_num)
    have hξ₁ε : ξ₁ < ε := lt_of_le_of_lt (min_le_left _ _) (by linarith)
    have hneg := hball (y := ξ₁) (by rw [Real.dist_eq, sub_zero, abs_of_pos hξ₁0]; exact hξ₁ε)
    have hns : ¬ specStable N ξ₁ := fun hs =>
      absurd ((specStable_iff hN2 hξ₁0 hξ₁1).mp hs) (not_le.mpr hneg)
    exact lt_of_lt_of_le hξ₁0 (le_thr (specStable N) hξ₁1 hns)
  -- (2) `ξ_* < 1`.
  have hlt1 : xiStar N < 1 := by
    obtain ⟨ξ₁, h0, h1, hP⟩ := exists_emDefinite_near_one (N := N) (by omega)
    exact lt_of_le_of_lt (thr_le (specStable N) h0.le h1.le fun ξ hξ =>
      specStable_of_emDefinite hN2 (h0.trans hξ.1) hξ.2 (hP ξ hξ)) h1
  set s := xiStar N with hs_def
  have hcont := continuousAt_aBlock (m := N / 2) hN1 hpos.le hlt1
  have habove : ∀ ξ ∈ Ioo s 1, specStable N ξ := fun ξ hξ =>
    holds_of_thr_lt (specStable N) hξ.1 hξ.2
  -- (3) `a_crit(ξ_*) = 0`.
  have hzero : aBlock N (N / 2) s = 0 := by
    rcases lt_trichotomy (aBlock N (N / 2) s) 0 with h | h | h
    · exfalso
      obtain ⟨ε, hε, hball⟩ := Metric.eventually_nhds_iff.mp
        (Filter.Tendsto.eventually_lt_const h hcont)
      set ξ := min (s + ε / 2) ((s + 1) / 2)
      have hξ1 : s < ξ := lt_min (by linarith) (by linarith)
      have hξ2 : ξ < 1 := lt_of_le_of_lt (min_le_right _ _) (by linarith)
      have hξ3 : ξ < s + ε := lt_of_le_of_lt (min_le_left _ _) (by linarith)
      have hneg := hball (y := ξ) (by rw [Real.dist_eq, abs_of_pos (by linarith)]; linarith)
      have := (specStable_iff hN2 (hpos.trans hξ1) hξ2).mp (habove ξ ⟨hξ1, hξ2⟩)
      linarith
    · exact h
    · exfalso
      obtain ⟨ε, hε, hball⟩ := Metric.eventually_nhds_iff.mp
        (Filter.Tendsto.eventually_const_lt h hcont)
      set ξ₀ := max (s - ε / 2) (s / 2)
      have hξ₀1 : ξ₀ < s := max_lt (by linarith) (by linarith)
      have hξ₀0 : 0 ≤ ξ₀ := le_trans (by linarith) (le_max_right _ _)
      have hξ₀2 : s - ε < ξ₀ := lt_of_lt_of_le (by linarith) (le_max_left _ _)
      have hle : s ≤ ξ₀ := by
        refine thr_le (specStable N) hξ₀0 (by linarith) fun ξ hξ => ?_
        rcases le_or_gt ξ s with hξs | hξs
        · have hξpos : 0 < ξ := lt_of_le_of_lt hξ₀0 hξ.1
          have hp := hball (y := ξ) (by
            rw [Real.dist_eq, abs_of_nonpos (by linarith)]; linarith [hξ.1])
          exact (specStable_iff hN2 hξpos hξ.2).mpr hp.le
        · exact habove ξ ⟨hξs, hξ.2⟩
      linarith
  -- (4) the determinant of the critical block is `−g² < 0` at `ξ_*`.
  have hg := gBlock_half_odd_pos hN2 hodd hpos hlt1
  have hdet : aBlock N (N / 2) s * bBlock N (N / 2) s - gBlock N (N / 2) s ^ 2 < 0 := by
    rw [hzero]; nlinarith
  obtain ⟨ε, hε, hball⟩ := Metric.eventually_nhds_iff.mp
    (Filter.Tendsto.eventually_lt_const hdet (continuousAt_det (m := N / 2) hN1 hpos.le hlt1))
  have hind : ∀ ξ ∈ Ioo s (s + ε),
      aBlock N (N / 2) ξ * bBlock N (N / 2) ξ - gBlock N (N / 2) ξ ^ 2 < 0 := fun ξ hξ =>
    hball (by rw [Real.dist_eq, abs_of_pos (by linarith [hξ.1])]; linarith [hξ.2])
  refine ⟨hpos, hlt1, hzero, habove, ⟨ε, hε, hind⟩, ?_⟩
  -- (5) `ξ_* < ξ_EM`.
  set ξ := min (s + ε / 2) ((s + 1) / 2)
  have hξ1 : s < ξ := lt_min (by linarith) (by linarith)
  have hξ2 : ξ < 1 := lt_of_le_of_lt (min_le_right _ _) (by linarith)
  have hξ3 : ξ < s + ε := lt_of_le_of_lt (min_le_left _ _) (by linarith)
  have hnot : ¬ emDefinite N ξ := fun h =>
    absurd (h (N / 2) (by omega) (by omega)).2 (not_lt.mpr (hind ξ ⟨hξ1, hξ3⟩).le)
  exact lt_of_lt_of_le hξ1 (le_thr (emDefinite N) hξ2 hnot)

end Havelock.Tanh.LargeN
