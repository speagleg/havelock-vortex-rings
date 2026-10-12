import Havelock.Tanh.LargeN.EvenN
import Havelock.Tanh.LargeN.Threshold

/-! # Tanh model, large `N`: even `N` (note Prop 15, Thm 18(a))

* `aBlock_half_strictMonoOn`: for even `N ≥ 4`, `ξ ↦ a_{N/2}(ξ)` is strictly increasing on `(0,1)`
  (note Prop 15: `8 a_{N/2}/N²` is strictly decreasing in `τ = −(N/2) ln ξ`).
* `xiEM_eq_xiStar_of_even`: for every even `N ≥ 2`, `ξ_EM(N) = ξ_*(N)` (note Thm 18(a)).
* `specStable_all_of_even_le_six`: for `N ∈ {2,4,6}` every radius is spectrally stable with a
  positive definite Hessian (note Prop 15).

The spacing statement of note Thm 16 is in `ThresholdSpacing.lean`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real Set Filter Topology

/-! ### note Prop 15: monotonicity of `a_{N/2}` in `ξ` -/

/-- note Prop 15: for even `N ≥ 4`, `ξ ↦ a_{N/2}(ξ)` is strictly increasing on `(0,1)`. -/
theorem aBlock_half_strictMonoOn {N : ℕ} (hE : Even N) (hN : 4 ≤ N) :
    StrictMonoOn (fun ξ => aBlock N (N / 2) ξ) (Ioo 0 1) := by
  intro ξ₁ h₁ ξ₂ h₂ h12
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  simp only
  rw [aBlock_half_eq hE (by omega) h₁.1 h₁.2, aBlock_half_eq hE (by omega) h₂.1 h₂.2]
  have ht := tauOf_lt_tauOf (N := N) (by omega) h₁.1 h12
  have hr := evenRHS_strictAntiOn_of_even hE hN
    (mem_Ioi.mpr (tauOf_pos (N := N) (by omega) h₂.1 h₂.2))
    (mem_Ioi.mpr (tauOf_pos (N := N) (by omega) h₁.1 h₁.2)) ht
  exact mul_lt_mul_of_pos_left hr (by positivity)

/-! ### note Thm 18(a) -/

/-- note Thm 18(a): for every even `N ≥ 2`, `ξ_EM(N) = ξ_*(N)` (no Krein window for even `N`). -/
theorem xiEM_eq_xiStar_of_even {N : ℕ} (hE : Even N) (hN : 2 ≤ N) : xiEM N = xiStar N := by
  rcases (show N = 2 ∨ 4 ≤ N by obtain ⟨n, rfl⟩ := hE; omega) with rfl | h4
  · exact xiEM_eq_xiStar_even le_rfl hE fun _ ξ₂ _ h12 h2 _ =>
      aBlock_half_pos_of_le_six hE le_rfl (by norm_num) (by linarith) h2
  · exact xiEM_eq_xiStar_of_strictMonoOn hN hE (aBlock_half_strictMonoOn hE h4)

/-! ### note Prop 15: `N ∈ {2,4,6}` -/

/-- note Prop 15: for even `N ∈ {2,4,6}` and every radius (`0 < ξ < 1`) the ring is spectrally
stable and the Hessian is positive definite on the non-rotational modes. -/
theorem specStable_all_of_even_le_six {N : ℕ} (hE : Even N) (hN2 : 2 ≤ N) (hN6 : N ≤ 6)
    {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) : specStable N ξ ∧ emDefinite N ξ := by
  have hem : emDefinite N ξ :=
    (emDefinite_iff_even hN2 hE h0 h1).mpr (aBlock_half_pos_of_le_six hE hN2 hN6 h0 h1)
  exact ⟨specStable_of_emDefinite hN2 h0 h1 hem, hem⟩

end Havelock.Tanh.LargeN
