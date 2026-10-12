import Havelock.Tanh.LargeN.Threshold
import Havelock.Tanh.LargeN.PhiLimit
import Mathlib.Analysis.SpecificLimits.Basic

/-! # Tanh model, large `N`: the limit of the thresholds (note Thm 16 and Thm 18(c))

Following the proof of note Thm 16: at `ξ = xiOf N τ`,
`a_{⌊N/2⌋}/N² = Φ(τ) + Δ_N(τ) + ε_N(τ)` with `Δ_N = E_N − E_∞`, `ε_N = ψ(0,τ) − ψ(x_N,τ) ≥ 0`.
Here `|Δ_N| = O(1/N²)` uniformly on `(0,8]` (via `δ = k − 1/τ = 𝓛(τ/N)/N ∈ (0, τ/N²]`) and
`ε_N ≤ (1/8 + (1+τ²)/4)/N²` (an explicit Taylor-type bound). With `Φ(τ∞ ∓ d) ≷ ±c₀ d` and
`d = K/N²`, the ring is positive definite (hence spectrally stable) for `τ ≤ τ∞ − d` and unstable
at `τ∞ + d`; so `|τ_*(N) − τ∞| ≤ K/N²`, and the same holds for `τ_EM(N)`. Finally
`cosh r = cR ξ = coth(τ/N)`. Main results: `cR_xiStar_rate`, `tendsto_cR_xiStar` (note Thm 16) and
`cR_xiEM_rate`, `tendsto_cR_xiEM` (note Thm 18(c)). -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real Set Filter Topology Asymptotics

/-! ### `Δ_N = E_N − E_∞` -/

/-- `k = coth(τ/N)/N = 1/τ + δ` with `0 < δ ≤ τ/N²` (proof of note Thm 16). -/
theorem kN_bounds {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    1 / τ < kN N τ ∧ kN N τ ≤ 1 / τ + τ / (N : ℝ) ^ 2 := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hu : 0 < τ / N := by positivity
  have h1 := langevin_pos hu
  have h2 := langevin_le_self hu
  have hk : kN N τ = 1 / τ + langevin (τ / N) / N := by
    unfold kN langevin
    field_simp
    ring
  rw [hk]
  constructor
  · have : 0 < langevin (τ / N) / N := by positivity
    linarith
  · have : langevin (τ / N) / N ≤ τ / (N : ℝ) ^ 2 := by
      rw [div_le_iff₀ hNpos]
      calc langevin (τ / N) ≤ τ / N := h2
        _ = τ / (N : ℝ) ^ 2 * N := by field_simp
    linarith

/-- `Δ_N = δ (coth τ − 1/τ − δ/2) − 1/(2N²)` with `δ = k − 1/τ` (proof of note Thm 16). -/
theorem EN_sub_Einf (N : ℕ) {τ : ℝ} (hτ : 0 < τ) :
    EN N τ - Einf τ = (kN N τ - 1 / τ) * (langevin τ - (kN N τ - 1 / τ) / 2)
      - 1 / (2 * (N : ℝ) ^ 2) := by
  unfold EN Einf langevin
  field_simp
  ring

/-- note Thm 16 (proof): `−(τ²+1)/(2N²) ≤ Δ_N(τ) ≤ τ²/N²` for `N ≥ 1`, `τ > 0`. -/
theorem EN_sub_Einf_bounds {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    -(τ ^ 2 + 1) / (2 * (N : ℝ) ^ 2) ≤ EN N τ - Einf τ ∧
      EN N τ - Einf τ ≤ τ ^ 2 / (N : ℝ) ^ 2 := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  have hN2 : (1 : ℝ) ≤ (N : ℝ) ^ 2 := by nlinarith
  obtain ⟨hk1, hk2⟩ := kN_bounds hN hτ
  have hL := langevin_pos hτ
  have hL2 := langevin_le_self hτ
  rw [EN_sub_Einf N hτ]
  set δ := kN N τ - 1 / τ with hδ
  set n2 := (N : ℝ) ^ 2
  have hδ0 : 0 < δ := by rw [hδ]; linarith
  have hδ1 : δ ≤ τ / n2 := by rw [hδ]; linarith
  have hn2 : 0 < n2 := by positivity
  have hτn : τ / n2 ≤ τ := div_le_self hτ.le hN2
  constructor
  · -- `δ²/2 ≤ τ²/(2 n2)`
    have hsq : δ ^ 2 ≤ τ ^ 2 / n2 := by
      have : δ ^ 2 ≤ (τ / n2) * τ := by nlinarith
      calc δ ^ 2 ≤ (τ / n2) * τ := this
        _ = τ ^ 2 / n2 := by ring
    have : -(τ ^ 2 + 1) / (2 * n2) = -(τ ^ 2 / n2) / 2 - 1 / (2 * n2) := by
      field_simp; ring
    rw [this]
    nlinarith
  · have h3 : δ * (langevin τ - δ / 2) ≤ δ * langevin τ := by nlinarith
    have h4 : δ * langevin τ ≤ τ / n2 * τ := by
      apply mul_le_mul hδ1 hL2 hL.le (by positivity)
    have h5 : 0 < 1 / (2 * n2) := by positivity
    have : τ / n2 * τ = τ ^ 2 / n2 := by ring
    linarith

/-! ### `ε_N = ψ(0,τ) − ψ(x_N,τ)` -/

/-- `τ coth τ ≤ 1 + τ²` for `τ > 0` (from `𝓛(τ) ≤ τ`). -/
theorem mul_coth_le {τ : ℝ} (hτ : 0 < τ) : τ * cosh τ ≤ (1 + τ ^ 2) * sinh τ := by
  have h := langevin_le_self hτ
  have hs : 0 < sinh τ := sinh_pos_iff.mpr hτ
  unfold langevin coth at h
  have : cosh τ / sinh τ ≤ τ + 1 / τ := by linarith
  rw [div_le_iff₀ hs] at this
  have h' : τ * cosh τ ≤ τ * ((τ + 1 / τ) * sinh τ) := mul_le_mul_of_nonneg_left this hτ.le
  calc τ * cosh τ ≤ τ * ((τ + 1 / τ) * sinh τ) := h'
    _ = (1 + τ ^ 2) * sinh τ := by field_simp; ring

/-- note Thm 16 (proof), explicit form of `ε_N ≤ C/N²`:
`ψ(0,τ) − ψ(x,τ) ≤ x² (1/8 + (1+τ²)/4)` for `0 ≤ x ≤ 1`, `τ > 0`. -/
theorem psi_zero_sub_le {x τ : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) (hτ : 0 < τ) :
    psi 0 τ - psi x τ ≤ x ^ 2 * (1 / 8 + (1 + τ ^ 2) / 4) := by
  have hs : 0 < sinh τ := sinh_pos_iff.mpr hτ
  have hxτ : 0 ≤ x * τ := mul_nonneg hx0 hτ.le
  have hsx : sinh (x * τ) ≤ x * τ * cosh (x * τ) := by
    rcases hxτ.eq_or_lt with h | h
    · rw [← h]; simp
    · exact (sinh_lt_mul_cosh h).le
  have hcx : cosh (x * τ) ≤ cosh τ := by
    rw [cosh_le_cosh, abs_of_nonneg hxτ, abs_of_pos hτ]
    nlinarith
  have hc1 := one_le_cosh (x * τ)
  have hcτ : 0 < cosh τ := cosh_pos τ
  have hsxn : 0 ≤ sinh (x * τ) := sinh_nonneg_iff.mpr hxτ
  have key : cosh τ - Q x τ ≤ x ^ 2 * sinh τ * (τ * cosh τ) := by
    unfold Q
    have e1 : cosh τ * (1 - cosh (x * τ)) ≤ 0 := by nlinarith
    have e2 : x * sinh τ * sinh (x * τ) ≤ x * sinh τ * (x * τ * cosh τ) := by
      apply mul_le_mul_of_nonneg_left _ (by positivity)
      calc sinh (x * τ) ≤ x * τ * cosh (x * τ) := hsx
        _ ≤ x * τ * cosh τ := mul_le_mul_of_nonneg_left hcx hxτ
    nlinarith
  have hmc := mul_coth_le hτ
  have key2 : cosh τ - Q x τ ≤ x ^ 2 * (1 + τ ^ 2) * sinh τ ^ 2 := by
    have : x ^ 2 * sinh τ * (τ * cosh τ) ≤ x ^ 2 * sinh τ * ((1 + τ ^ 2) * sinh τ) :=
      mul_le_mul_of_nonneg_left hmc (by positivity)
    nlinarith
  have hpsi : psi 0 τ - psi x τ = x ^ 2 / 8 + (cosh τ - Q x τ) / (4 * sinh τ ^ 2) := by
    unfold psi
    have : Q 0 τ = cosh τ := by simp [Q]
    rw [this]
    field_simp
    ring
  rw [hpsi]
  have : (cosh τ - Q x τ) / (4 * sinh τ ^ 2) ≤ x ^ 2 * (1 + τ ^ 2) / 4 := by
    rw [div_le_iff₀ (by positivity)]
    nlinarith
  nlinarith

/-- The critical mode `m = ⌊N/2⌋` has `0 ≤ x_N ≤ 1/N`. -/
theorem xOf_half_bounds {N : ℕ} (hN : 1 ≤ N) :
    0 ≤ xOf N (N / 2) ∧ xOf N (N / 2) ≤ 1 / N := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have h1 : (2 * ((N / 2 : ℕ) : ℝ)) ≤ N := by exact_mod_cast (show 2 * (N / 2) ≤ N by omega)
  have h2 : (N : ℝ) ≤ 2 * ((N / 2 : ℕ) : ℝ) + 1 := by
    exact_mod_cast (show N ≤ 2 * (N / 2) + 1 by omega)
  have hx : xOf N (N / 2) = ((N : ℝ) - 2 * ((N / 2 : ℕ) : ℝ)) / N := by
    unfold xOf; field_simp
  rw [hx]
  exact ⟨div_nonneg (by linarith) hNpos.le, div_le_div_of_nonneg_right (by linarith) hNpos.le⟩

/-- note Thm 16 (proof): `a_{⌊N/2⌋}/N² = Φ + Δ_N + ε_N`. -/
theorem aBlock_half_decomp {N : ℕ} (hN : 2 ≤ N) {τ : ℝ} (hτ : 0 < τ) :
    aBlock N (N / 2) (xiOf N τ) = (N : ℝ) ^ 2 *
      (Phi τ + (EN N τ - Einf τ) + (psi 0 τ - psi (xOf N (N / 2)) τ)) := by
  rw [aBlock_xiOf (by omega) (by omega) hτ]
  unfold Phi
  ring

/-- `Φ_N = Φ + Δ_N`. -/
theorem PhiN_eq (N : ℕ) (τ : ℝ) : PhiN N τ = Phi τ + (EN N τ - Einf τ) := by
  unfold PhiN Phi; ring

/-! ### Linear behaviour of `Φ` at `τ∞` -/

/-- From `Φ'(τ∞) < 0`: `Φ(τ∞ − d) ≥ c₀ d` and `Φ(τ∞ + d) ≤ −c₀ d` for `0 < d ≤ d₀`
(proof of note Thm 16, `c₀ = |Φ'(τ∞)|/2`). -/
theorem Phi_linear_at_tauInf :
    ∃ c₀ > 0, ∃ d₀ > 0, ∀ d : ℝ, 0 < d → d ≤ d₀ →
      c₀ * d ≤ Phi (tauInf - d) ∧ Phi (tauInf + d) ≤ -(c₀ * d) := by
  have hD := hasDerivAt_Phi tauInf_pos
  have hneg := PhiDeriv_neg tauInf_pos
  set f' := PhiDeriv tauInf
  set c₀ := -f' / 2 with hc₀
  have hc0 : 0 < c₀ := by rw [hc₀]; linarith
  have hlo := hD.isLittleO.def hc0
  obtain ⟨ε, hε, hball⟩ := Metric.eventually_nhds_iff.mp hlo
  refine ⟨c₀, hc0, ε / 2, by positivity, fun d hd hdε => ?_⟩
  have hdε' : d < ε := by linarith
  have e1 := hball (y := tauInf - d) (by
    rw [Real.dist_eq, show tauInf - d - tauInf = -d by ring, abs_neg, abs_of_pos hd]; exact hdε')
  have e2 := hball (y := tauInf + d) (by
    rw [Real.dist_eq, show tauInf + d - tauInf = d by ring, abs_of_pos hd]; exact hdε')
  simp only [Phi_tauInf, smul_eq_mul, Real.norm_eq_abs, sub_zero] at e1 e2
  rw [show tauInf - d - tauInf = -d by ring, abs_neg, abs_of_pos hd] at e1
  rw [show tauInf + d - tauInf = d by ring, abs_of_pos hd] at e2
  have e1' := (abs_le.mp e1).1
  have e2' := (abs_le.mp e2).2
  constructor
  · nlinarith
  · nlinarith

/-! ### The window `[τ∞ − K/N², τ∞ + K/N²]` -/

/-- note Thm 16 (proof): there are `K > 0` and `N₀` such that for `N ≥ N₀`, with `d = K/N²`:
`Φ_N > 0` on `(0, τ∞ − d]` and `a_{⌊N/2⌋} < 0` at `τ∞ + d`. -/
theorem crit_window :
    ∃ Kc > 0, ∃ N₀ : ℕ, ∀ N : ℕ, N₀ ≤ N →
      2 ≤ N ∧ Kc / (N : ℝ) ^ 2 ≤ 1 / 2 ∧
      (∀ τ : ℝ, 0 < τ → τ ≤ tauInf - Kc / (N : ℝ) ^ 2 → 0 < PhiN N τ) ∧
      aBlock N (N / 2) (xiOf N (tauInf + Kc / (N : ℝ) ^ 2)) < 0 := by
  obtain ⟨c₀, hc0, d₀, hd0, hlin⟩ := Phi_linear_at_tauInf
  set Kc := 100 / c₀ with hKc
  have hKpos : 0 < Kc := by positivity
  set D := min d₀ (1 / 2) with hD
  have hDpos : 0 < D := lt_min hd0 (by norm_num)
  refine ⟨Kc, hKpos, max 2 ⌈Kc / D⌉₊, fun N hN => ?_⟩
  have hN2 : 2 ≤ N := le_trans (le_max_left _ _) hN
  have hNc : Kc / D ≤ N := le_trans (Nat.le_ceil _) (by exact_mod_cast le_trans (le_max_right _ _) hN)
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast (show 1 ≤ N by omega)
  set n2 := (N : ℝ) ^ 2 with hn2
  have hn2pos : 0 < n2 := by positivity
  have hn2N : (N : ℝ) ≤ n2 := by rw [hn2]; nlinarith
  set d := Kc / n2 with hd
  have hdpos : 0 < d := by positivity
  have hdD : d ≤ D := by
    rw [hd, div_le_iff₀ hn2pos]
    rw [div_le_iff₀ hDpos] at hNc
    nlinarith
  have hdd0 : d ≤ d₀ := le_trans hdD (min_le_left _ _)
  have hdh : d ≤ 1 / 2 := le_trans hdD (min_le_right _ _)
  have hc0d : c₀ * d = 100 / n2 := by rw [hd, hKc]; field_simp
  obtain ⟨hlo, hhi⟩ := hlin d hdpos hdd0
  have hT := tauInf_bracket
  refine ⟨hN2, hdh, ?_, ?_⟩
  · intro τ hτ hτd
    have hmono : Phi (tauInf - d) ≤ Phi τ :=
      Phi_strictAntiOn.antitoneOn (show 0 < τ from hτ) (show 0 < tauInf - d by linarith) hτd
    have hΔ := (EN_sub_Einf_bounds (show 1 ≤ N by omega) hτ).1
    have hτ8 : τ ^ 2 + 1 ≤ 66 := by nlinarith
    have hΔ' : -(33 / n2) ≤ EN N τ - Einf τ := by
      have : -(τ ^ 2 + 1) / (2 * n2) ≥ -(33 / n2) := by
        rw [ge_iff_le, neg_div, neg_le_neg_iff, div_le_div_iff₀ (by positivity) hn2pos]
        nlinarith
      linarith
    rw [PhiN_eq]
    have : 0 < 100 / n2 - 33 / n2 := by
      rw [← sub_div]; positivity
    linarith
  · have hτ1 : 0 < tauInf + d := by linarith
    rw [aBlock_half_decomp hN2 hτ1]
    have hΔ := (EN_sub_Einf_bounds (show 1 ≤ N by omega) hτ1).2
    obtain ⟨hx0, hx1⟩ := xOf_half_bounds (show 1 ≤ N by omega)
    have hε := psi_zero_sub_le hx0 (le_trans hx1 (by
      rw [div_le_one hNpos]; exact hN1)) hτ1
    have hx2 : xOf N (N / 2) ^ 2 ≤ 1 / n2 := by
      have : xOf N (N / 2) ^ 2 ≤ (1 / (N : ℝ)) ^ 2 := pow_le_pow_left₀ hx0 hx1 2
      rw [hn2]; simpa [div_pow] using this
    have hτsq : (tauInf + d) ^ 2 ≤ 64 := by nlinarith
    have hε' : psi 0 (tauInf + d) - psi (xOf N (N / 2)) (tauInf + d) ≤ 17 / n2 := by
      have hc : 1 / 8 + (1 + (tauInf + d) ^ 2) / 4 ≤ 17 := by nlinarith
      calc _ ≤ xOf N (N / 2) ^ 2 * (1 / 8 + (1 + (tauInf + d) ^ 2) / 4) := hε
        _ ≤ (1 / n2) * 17 := mul_le_mul hx2 hc (by positivity) (by positivity)
        _ = 17 / n2 := by ring
    have hΔ' : EN N (tauInf + d) - Einf (tauInf + d) ≤ 64 / n2 := by
      calc _ ≤ (tauInf + d) ^ 2 / n2 := hΔ
        _ ≤ 64 / n2 := div_le_div_of_nonneg_right hτsq hn2pos.le
    have hsum : Phi (tauInf + d) + (EN N (tauInf + d) - Einf (tauInf + d))
        + (psi 0 (tauInf + d) - psi (xOf N (N / 2)) (tauInf + d)) < 0 := by
      have : 64 / n2 + 17 / n2 - 100 / n2 < 0 := by
        rw [← add_div, ← sub_div, div_neg_iff]; right; exact ⟨by norm_num, hn2pos⟩
      linarith
    exact mul_neg_of_pos_of_neg hn2pos hsum

/-! ### From the window to the thresholds -/

theorem xiOf_le_iff {N : ℕ} (hN : 1 ≤ N) {a b : ℝ} : xiOf N a ≤ xiOf N b ↔ b ≤ a := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  unfold xiOf
  rw [exp_le_exp, div_le_div_iff_of_pos_right hNpos]
  constructor <;> intro h <;> linarith

theorem xiOf_lt_iff {N : ℕ} (hN : 1 ≤ N) {a b : ℝ} : xiOf N a < xiOf N b ↔ b < a := by
  rw [← not_le, xiOf_le_iff hN, not_le]

/-- Window bounds for a threshold predicate `P` that holds wherever `Φ_N > 0` and fails wherever
`a_{⌊N/2⌋} < 0` (proof of note Thm 16): `xiOf N (τ∞ + d) ≤ ξ_P(N) ≤ xiOf N (τ∞ − d)`. -/
theorem thr_bounds (P : ℕ → ℝ → Prop)
    (hP1 : ∀ (N : ℕ) (τ : ℝ), 2 ≤ N → 0 < τ → 0 < PhiN N τ → P N (xiOf N τ))
    (hP2 : ∀ (N : ℕ) (ξ : ℝ), 2 ≤ N → 0 < ξ → ξ < 1 → aBlock N (N / 2) ξ < 0 → ¬ P N ξ)
    {Kc : ℝ} (hKpos : 0 < Kc) {N : ℕ}
    (hwin : 2 ≤ N ∧ Kc / (N : ℝ) ^ 2 ≤ 1 / 2 ∧
      (∀ τ : ℝ, 0 < τ → τ ≤ tauInf - Kc / (N : ℝ) ^ 2 → 0 < PhiN N τ) ∧
      aBlock N (N / 2) (xiOf N (tauInf + Kc / (N : ℝ) ^ 2)) < 0) :
    xiOf N (tauInf + Kc / (N : ℝ) ^ 2)
        ≤ sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P N ξ} ∧
      sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P N ξ}
        ≤ xiOf N (tauInf - Kc / (N : ℝ) ^ 2) := by
  obtain ⟨hN2, hdh, hpos, hneg⟩ := hwin
  have hN1 : 1 ≤ N := by omega
  have hT := tauInf_bracket
  have hlo_pos : 0 < tauInf - Kc / (N : ℝ) ^ 2 := by linarith
  have hhi_pos : 0 < tauInf + Kc / (N : ℝ) ^ 2 := by
    have : 0 < Kc / (N : ℝ) ^ 2 := by
      have : (0 : ℝ) < N := by exact_mod_cast hN1
      positivity
    linarith
  constructor
  · exact le_thr (P N) (xiOf_lt_one hN1 hhi_pos)
      (hP2 N _ hN2 (xiOf_pos N _) (xiOf_lt_one hN1 hhi_pos) hneg)
  · refine thr_le (P N) (xiOf_pos N _).le (xiOf_lt_one hN1 hlo_pos).le fun ξ hξ => ?_
    obtain ⟨τ, hτ, rfl⟩ := exists_tau hN1 ((xiOf_pos N _).trans hξ.1) hξ.2
    have : τ < tauInf - Kc / (N : ℝ) ^ 2 := (xiOf_lt_iff hN1).mp hξ.1
    exact hP1 N τ hN2 hτ (hpos τ hτ this.le)

/-- Generic rate for a threshold predicate `P` that holds wherever `Φ_N > 0` and fails wherever
`a_{⌊N/2⌋} < 0` (proof of note Thm 16): `|cR(ξ_P(N)) − κ∞ N| ≤ C/N` for `N ≥ N₀`. -/
theorem cR_thr_rate (P : ℕ → ℝ → Prop)
    (hP1 : ∀ (N : ℕ) (τ : ℝ), 2 ≤ N → 0 < τ → 0 < PhiN N τ → P N (xiOf N τ))
    (hP2 : ∀ (N : ℕ) (ξ : ℝ), 2 ≤ N → 0 < ξ → ξ < 1 → aBlock N (N / 2) ξ < 0 → ¬ P N ξ) :
    ∃ C : ℝ, 0 < C ∧ ∃ N₀ : ℕ, ∀ N : ℕ, N₀ ≤ N →
      |cR (sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P N ξ}) - kappaInf * N|
        ≤ C / N := by
  obtain ⟨Kc, hKpos, N₀, hwin⟩ := crit_window
  refine ⟨8 + Kc, by positivity, N₀, fun N hN => ?_⟩
  obtain ⟨ht_ge, ht_le⟩ := thr_bounds P hP1 hP2 hKpos (hwin N hN)
  obtain ⟨hN2, hdh, -, -⟩ := hwin N hN
  have hN1 : 1 ≤ N := by omega
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN1
  set d := Kc / (N : ℝ) ^ 2 with hd
  have hdpos : 0 < d := by positivity
  have hT := tauInf_bracket
  set t := sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P N ξ} with ht
  have hlo_pos : 0 < tauInf - d := by linarith
  have ht0 : 0 < t := lt_of_lt_of_le (xiOf_pos N _) ht_ge
  have ht1 : t < 1 := lt_of_le_of_lt ht_le (xiOf_lt_one hN1 hlo_pos)
  obtain ⟨τs, hτs, hτs_eq⟩ := exists_tau hN1 ht0 ht1
  rw [← hτs_eq] at ht_le ht_ge ⊢
  have hτs1 : tauInf - d ≤ τs := (xiOf_le_iff hN1).mp ht_le
  have hτs2 : τs ≤ tauInf + d := (xiOf_le_iff hN1).mp ht_ge
  rw [cR_xiOf hN1 hτs]
  -- `coth(τs/N) = N/τs + 𝓛(τs/N)`, `0 < 𝓛(τs/N) ≤ τs/N`
  have hu : 0 < τs / N := by positivity
  have hL1 := langevin_pos hu
  have hL2 := langevin_le_self hu
  have hcoth : coth (τs / N) = N / τs + langevin (τs / N) := by
    unfold langevin; field_simp; ring
  have hκ : kappaInf * N = N / tauInf := by unfold kappaInf; ring
  rw [hcoth, hκ]
  have hTpos : 0 < tauInf := tauInf_pos
  have hdiff : N / τs - N / tauInf = N * (tauInf - τs) / (τs * tauInf) := by
    field_simp
  have hprod : 1 ≤ τs * tauInf := by nlinarith
  have habs : |tauInf - τs| ≤ d := abs_le.mpr ⟨by linarith, by linarith⟩
  have hB : |N / τs - N / tauInf| ≤ Kc / N := by
    rw [hdiff, abs_div, abs_mul, abs_of_pos hNpos, abs_of_pos (by positivity : 0 < τs * tauInf)]
    calc (N : ℝ) * |tauInf - τs| / (τs * tauInf) ≤ (N : ℝ) * |tauInf - τs| :=
          div_le_self (by positivity) hprod
      _ ≤ N * d := mul_le_mul_of_nonneg_left habs hNpos.le
      _ = Kc / N := by rw [hd]; field_simp
  have hA : |langevin (τs / N)| ≤ 8 / N := by
    rw [abs_of_pos hL1]
    calc langevin (τs / N) ≤ τs / N := hL2
      _ ≤ 8 / N := div_le_div_of_nonneg_right (by linarith) hNpos.le
  calc |N / τs + langevin (τs / N) - N / tauInf|
      = |(N / τs - N / tauInf) + langevin (τs / N)| := by ring_nf
    _ ≤ |N / τs - N / tauInf| + |langevin (τs / N)| := abs_add_le _ _
    _ ≤ Kc / N + 8 / N := add_le_add hB hA
    _ = (8 + Kc) / N := by ring

/-- The rate implies the limit `cR(ξ_P(N))/N → κ∞`. -/
theorem tendsto_of_rate (f : ℕ → ℝ) (C : ℝ) (N₀ : ℕ)
    (h : ∀ N : ℕ, N₀ ≤ N → |f N - kappaInf * N| ≤ C / N) :
    Tendsto (fun N : ℕ => f N / N) atTop (𝓝 kappaInf) := by
  rw [tendsto_iff_dist_tendsto_zero]
  refine squeeze_zero' (Eventually.of_forall fun _ => dist_nonneg) ?_
    ((tendsto_const_div_atTop_nhds_zero_nat C).div_atTop tendsto_natCast_atTop_atTop)
  filter_upwards [eventually_ge_atTop (max N₀ 1)] with N hN
  have hN0 : N₀ ≤ N := le_trans (le_max_left _ _) hN
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (lt_of_lt_of_le zero_lt_one (le_trans (le_max_right _ _) hN))
  rw [Real.dist_eq, show f N / N - kappaInf = (f N - kappaInf * N) / N by field_simp,
    abs_div, abs_of_pos hNpos]
  exact div_le_div_of_nonneg_right (h N hN0) hNpos.le

/-! ### note Thm 16 and Thm 18(c) -/

theorem specStable_hP1 (N : ℕ) (τ : ℝ) (hN : 2 ≤ N) (hτ : 0 < τ) (hΦ : 0 < PhiN N τ) :
    specStable N (xiOf N τ) :=
  specStable_of_emDefinite hN (xiOf_pos N τ) (xiOf_lt_one (by omega) hτ)
    (emDefinite_of_PhiN_pos hN hτ hΦ)

theorem specStable_hP2 (N : ℕ) (ξ : ℝ) (hN : 2 ≤ N) (h0 : 0 < ξ) (h1 : ξ < 1)
    (ha : aBlock N (N / 2) ξ < 0) : ¬ specStable N ξ := fun hs =>
  absurd ((specStable_iff hN h0 h1).mp hs) (not_le.mpr ha)

theorem emDefinite_hP2 (N : ℕ) (ξ : ℝ) (hN : 2 ≤ N) (_ : 0 < ξ) (_ : ξ < 1)
    (ha : aBlock N (N / 2) ξ < 0) : ¬ emDefinite N ξ := fun hs =>
  absurd (hs (N / 2) (by omega) (by omega)).1 (not_lt.mpr ha.le)


/-- note Thm 16 (rate): `cosh r_*(N) = κ∞ N + O(1/N)`, both parities. -/
theorem cR_xiStar_rate :
    ∃ C : ℝ, 0 < C ∧ ∃ N₀ : ℕ, ∀ N : ℕ, N₀ ≤ N → |cR (xiStar N) - kappaInf * N| ≤ C / N :=
  cR_thr_rate specStable specStable_hP1 specStable_hP2

/-- note Thm 16: `cosh r_*(N)/N → κ∞`, both parities. -/
theorem tendsto_cR_xiStar :
    Tendsto (fun N : ℕ => cR (xiStar N) / N) atTop (𝓝 kappaInf) := by
  obtain ⟨C, -, N₀, h⟩ := cR_xiStar_rate
  exact tendsto_of_rate (fun N => cR (xiStar N)) C N₀ h

/-- note Thm 18(c) (rate): `cosh r_EM(N) = κ∞ N + O(1/N)`, both parities. -/
theorem cR_xiEM_rate :
    ∃ C : ℝ, 0 < C ∧ ∃ N₀ : ℕ, ∀ N : ℕ, N₀ ≤ N → |cR (xiEM N) - kappaInf * N| ≤ C / N :=
  cR_thr_rate emDefinite (fun _ _ hN hτ hΦ => emDefinite_of_PhiN_pos hN hτ hΦ) emDefinite_hP2

/-- note Thm 18(c): `cosh r_EM(N)/N → κ∞`, both parities. -/
theorem tendsto_cR_xiEM :
    Tendsto (fun N : ℕ => cR (xiEM N) / N) atTop (𝓝 kappaInf) := by
  obtain ⟨C, -, N₀, h⟩ := cR_xiEM_rate
  exact tendsto_of_rate (fun N => cR (xiEM N)) C N₀ h

/-- `e^{−2a/N} − e^{−2b/N} ≤ 2(b − a)/N` for `0 ≤ a ≤ b`. -/
theorem xiOf_sub_le {N : ℕ} (hN : 1 ≤ N) {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) :
    xiOf N a - xiOf N b ≤ 2 * (b - a) / N := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  unfold xiOf
  set u := -2 * a / N
  set v := -2 * b / N
  have hu : u ≤ 0 := by
    have : 0 ≤ 2 * a / N := by positivity
    simp only [u]; rw [show -2 * a / N = -(2 * a / N) by ring]; linarith
  have h1 := add_one_le_exp (v - u)
  rw [exp_sub] at h1
  have heu : 0 < exp u := exp_pos u
  have heu1 : exp u ≤ 1 := exp_le_one_iff.mpr hu
  have h2 : exp u * (v - u + 1) ≤ exp v := by
    rw [le_div_iff₀ heu] at h1; linarith
  have huv : u - v = 2 * (b - a) / N := by simp only [u, v]; field_simp; ring
  have hba : 0 ≤ b - a := by linarith
  have hnn : 0 ≤ u - v := by rw [huv]; positivity
  nlinarith

/-- note Thm 18(b), width (here for both parities): `0 ≤ ξ_EM(N) − ξ_*(N) = O(N⁻³)`. -/
theorem xiEM_sub_xiStar_le :
    ∃ C : ℝ, 0 < C ∧ ∃ N₀ : ℕ, ∀ N : ℕ, N₀ ≤ N →
      0 ≤ xiEM N - xiStar N ∧ xiEM N - xiStar N ≤ C / (N : ℝ) ^ 3 := by
  obtain ⟨Kc, hKpos, N₀, hwin⟩ := crit_window
  refine ⟨4 * Kc, by positivity, N₀, fun N hN => ?_⟩
  have hS := thr_bounds specStable specStable_hP1 specStable_hP2 hKpos (hwin N hN)
  have hE := thr_bounds emDefinite (fun _ _ hN hτ hΦ => emDefinite_of_PhiN_pos hN hτ hΦ)
    emDefinite_hP2 hKpos (hwin N hN)
  obtain ⟨hN2, hdh, -, -⟩ := hwin N hN
  have hN1 : 1 ≤ N := by omega
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN1
  have hT := tauInf_bracket
  have hdpos : 0 < Kc / (N : ℝ) ^ 2 := by positivity
  refine ⟨by linarith [xiStar_le_xiEM hN2], ?_⟩
  have hw := xiOf_sub_le hN1 (a := tauInf - Kc / (N : ℝ) ^ 2) (b := tauInf + Kc / (N : ℝ) ^ 2)
    (by linarith) (by linarith)
  have : 2 * (tauInf + Kc / (N : ℝ) ^ 2 - (tauInf - Kc / (N : ℝ) ^ 2)) / N = 4 * Kc / (N : ℝ) ^ 3 := by
    field_simp; ring
  have h1 : xiEM N ≤ xiOf N (tauInf - Kc / (N : ℝ) ^ 2) := hE.2
  have h2 : xiOf N (tauInf + Kc / (N : ℝ) ^ 2) ≤ xiStar N := hS.1
  linarith

end Havelock.Tanh.LargeN
