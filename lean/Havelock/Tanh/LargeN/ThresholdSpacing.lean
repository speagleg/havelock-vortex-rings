import Havelock.Tanh.LargeN.ThresholdLimit

/-! # Tanh model, large `N`: the critical spacing (note Thm 16)

With `cosh r_*(N) = cR(ξ_*(N))` and `sinh r_*(N) = √(cosh² r_*(N) − 1)`, the arc-length spacing
`ℓ_N = 2π sinh r_*(N)/N` satisfies `ℓ_N = 2π κ∞ + O(N⁻²)` (`spacing_rate`) and hence
`ℓ_N → ℓ_* = 2π κ∞` (`tendsto_spacing`). Both follow from `cR_xiStar_rate`
(`cosh r_* = κ∞ N + O(1/N)`) and `|√(c² − 1) − c| ≤ 1/c` for `c ≥ 1`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real Set Filter Topology

/-! ### note Thm 16: the critical spacing -/

/-- For `c ≥ 1`: `|√(c² − 1) − c| ≤ 1/c`. -/
theorem abs_sqrt_sq_sub_one_sub_le {c : ℝ} (hc : 1 ≤ c) : |Real.sqrt (c ^ 2 - 1) - c| ≤ 1 / c := by
  have hc0 : 0 < c := by linarith
  have hup : Real.sqrt (c ^ 2 - 1) ≤ c := by
    calc Real.sqrt (c ^ 2 - 1) ≤ Real.sqrt (c ^ 2) := Real.sqrt_le_sqrt (by linarith)
      _ = c := Real.sqrt_sq hc0.le
  have hinv : 1 / c ≤ 1 := by rw [div_le_one hc0]; exact hc
  have hinv0 : 0 < 1 / c := by positivity
  have hlo0 : 0 ≤ c - 1 / c := by linarith
  have hsq : (c - 1 / c) ^ 2 ≤ c ^ 2 - 1 := by
    have : (c - 1 / c) ^ 2 = c ^ 2 - 2 + (1 / c) ^ 2 := by field_simp; ring
    rw [this]
    nlinarith
  have hlo : c - 1 / c ≤ Real.sqrt (c ^ 2 - 1) := by
    calc c - 1 / c = Real.sqrt ((c - 1 / c) ^ 2) := (Real.sqrt_sq hlo0).symm
      _ ≤ Real.sqrt (c ^ 2 - 1) := Real.sqrt_le_sqrt hsq
  rw [abs_le]
  constructor <;> linarith

/-- note Thm 16 (spacing, rate): with `cosh r_*(N) = cR(ξ_*(N))` and
`sinh r_*(N) = √(cosh² r_*(N) − 1)`, the spacing `ℓ_N = 2π sinh r_*(N)/N` satisfies
`|ℓ_N − 2π κ∞| ≤ C/N²` for `N ≥ N₀`. -/
theorem spacing_rate :
    ∃ C : ℝ, 0 < C ∧ ∃ N₀ : ℕ, ∀ N : ℕ, N₀ ≤ N →
      |2 * π * Real.sqrt (cR (xiStar N) ^ 2 - 1) / N - 2 * π * kappaInf| ≤ C / (N : ℝ) ^ 2 := by
  obtain ⟨C, hC, N₀, h⟩ := cR_xiStar_rate
  have hκ : 0 < kappaInf := lt_trans (by norm_num) kappaInf_bracket.1
  refine ⟨2 * π * (2 / kappaInf + C), by positivity,
    max (max N₀ 1) ⌈2 * (C + 1) / kappaInf⌉₊, fun N hN => ?_⟩
  have hN0 : N₀ ≤ N := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hN
  have hN1 : 1 ≤ N := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hN
  have hNc : 2 * (C + 1) / kappaInf ≤ N :=
    le_trans (Nat.le_ceil _) (by exact_mod_cast le_trans (le_max_right _ _) hN)
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN1
  have hN1' : (1 : ℝ) ≤ N := by exact_mod_cast hN1
  have hκN : 2 * (C + 1) ≤ kappaInf * N := by
    rw [div_le_iff₀ hκ] at hNc; linarith
  set c := cR (xiStar N) with hc_def
  have hcN := h N hN0
  rw [← hc_def] at hcN
  have hCN : C / N ≤ C := div_le_self hC.le hN1'
  have hcl : kappaInf * N / 2 ≤ c := by
    have := (abs_le.mp hcN).1
    linarith
  have hc1 : 1 ≤ c := by linarith
  have hc0 : 0 < c := by linarith
  have hs := abs_sqrt_sq_sub_one_sub_le hc1
  have h1c : 1 / c ≤ 2 / kappaInf / N := by
    rw [div_div, div_le_div_iff₀ hc0 (by positivity)]
    linarith
  have htot : |Real.sqrt (c ^ 2 - 1) - kappaInf * N| ≤ (2 / kappaInf + C) / N := by
    calc |Real.sqrt (c ^ 2 - 1) - kappaInf * N|
        = |(Real.sqrt (c ^ 2 - 1) - c) + (c - kappaInf * N)| := by ring_nf
      _ ≤ |Real.sqrt (c ^ 2 - 1) - c| + |c - kappaInf * N| := abs_add_le _ _
      _ ≤ 2 / kappaInf / N + C / N := add_le_add (hs.trans h1c) hcN
      _ = (2 / kappaInf + C) / N := by ring
  have heq : 2 * π * Real.sqrt (c ^ 2 - 1) / N - 2 * π * kappaInf
      = 2 * π * (Real.sqrt (c ^ 2 - 1) - kappaInf * N) / N := by
    field_simp
  rw [heq, abs_div, abs_mul, abs_of_pos hNpos, abs_of_pos (by positivity : (0 : ℝ) < 2 * π)]
  rw [div_le_div_iff₀ hNpos (by positivity)]
  have := mul_le_mul_of_nonneg_left htot (by positivity : (0 : ℝ) ≤ 2 * π)
  calc 2 * π * |Real.sqrt (c ^ 2 - 1) - kappaInf * N| * (N : ℝ) ^ 2
      ≤ 2 * π * ((2 / kappaInf + C) / N) * (N : ℝ) ^ 2 :=
        mul_le_mul_of_nonneg_right this (by positivity)
    _ = 2 * π * (2 / kappaInf + C) * N := by field_simp

/-- note Thm 16 (spacing): `ℓ_N = 2π sinh r_*(N)/N → ℓ_* = 2π κ∞`. -/
theorem tendsto_spacing :
    Tendsto (fun N : ℕ => 2 * π * Real.sqrt (cR (xiStar N) ^ 2 - 1) / N) atTop
      (𝓝 (2 * π * kappaInf)) := by
  obtain ⟨C, hC, N₀, h⟩ := spacing_rate
  rw [tendsto_iff_dist_tendsto_zero]
  refine squeeze_zero' (Eventually.of_forall fun _ => dist_nonneg) ?_
    (tendsto_const_div_atTop_nhds_zero_nat C)
  filter_upwards [eventually_ge_atTop (max N₀ 1)] with N hN
  have hN1' : (1 : ℝ) ≤ N := by exact_mod_cast le_trans (le_max_right _ _) hN
  have hNpos : (0 : ℝ) < N := by linarith
  rw [Real.dist_eq]
  refine (h N (le_trans (le_max_left _ _) hN)).trans ?_
  rw [div_le_div_iff₀ (by positivity) hNpos]
  have hNN : (N : ℝ) ≤ (N : ℝ) ^ 2 := by nlinarith
  exact mul_le_mul_of_nonneg_left hNN hC.le

end Havelock.Tanh.LargeN
