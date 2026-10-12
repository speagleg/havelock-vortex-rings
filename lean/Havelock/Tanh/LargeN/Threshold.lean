import Havelock.Tanh.LargeN.CritMode
import Havelock.Tanh.LargeN.EMCore

/-! # Tanh model, large `N`: spectral and energy–momentum thresholds (note Cor 14, Lemma 17, Thm 18(a))

* `specStable N ξ`: every block `m = 1, …, N−1` has `a_m b_m ≥ 0` (note Prop 4: spectral stability).
* `emDefinite N ξ`: every block is positive definite, `a_m > 0` and `a_m b_m − g_m² > 0`
  (note Prop 4 / Thm 18: the Hessian is positive definite on the non-rotational modes).
* `xiStar N`, `xiEM N`: the thresholds `inf {ξ₀ ∈ [0,1] : P holds on (ξ₀, 1)}`.

Results: note Cor 14 (`specStable_iff`), note Lemma 17 (`det_xiOf`, `G_nonneg`,
`emDefinite_of_PhiN_pos`), note Thm 18(a) (`emDefinite_iff_even`; the equality of the two
thresholds for even `N` is `xiEM_eq_xiStar_of_strictMonoOn`, under the monotonicity of
`a_{N/2}` in `ξ` that note Prop 15 provides). -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real Set

/-- Spectral stability (note Prop 4 with Thm 2): `a_m b_m ≥ 0` for every `1 ≤ m ≤ N−1`. -/
def specStable (N : ℕ) (ξ : ℝ) : Prop :=
  ∀ m : ℕ, 1 ≤ m → m < N → 0 ≤ aBlock N m ξ * bBlock N m ξ

/-- Positive definiteness of the Hessian on the non-rotational modes (note Prop 4, Thm 18):
`a_m > 0` and `a_m b_m − g_m² > 0` for every `1 ≤ m ≤ N−1`. -/
def emDefinite (N : ℕ) (ξ : ℝ) : Prop :=
  ∀ m : ℕ, 1 ≤ m → m < N →
    0 < aBlock N m ξ ∧ 0 < aBlock N m ξ * bBlock N m ξ - gBlock N m ξ ^ 2

/-- The spectral threshold `ξ_*(N)` (note Thm 16). -/
def xiStar (N : ℕ) : ℝ :=
  sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, specStable N ξ}

/-- The energy–momentum threshold `ξ_EM(N)` (note Thm 18). -/
def xiEM (N : ℕ) : ℝ :=
  sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, emDefinite N ξ}

/-! ### Generic facts about thresholds `inf {ξ₀ ∈ [0,1] : P on (ξ₀,1)}` -/

section Generic

variable (P : ℝ → Prop)

theorem thr_one_mem : (1 : ℝ) ∈ {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P ξ} :=
  ⟨⟨zero_le_one, le_rfl⟩, fun _ h => absurd h.2 (not_lt.mpr h.1.le)⟩

theorem thr_bdd : BddBelow {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P ξ} :=
  ⟨0, fun _ h => h.1.1⟩

theorem thr_nonneg : 0 ≤ sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P ξ} :=
  le_csInf ⟨1, thr_one_mem P⟩ fun _ h => h.1.1

theorem thr_le_one : sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P ξ} ≤ 1 :=
  csInf_le (thr_bdd P) (thr_one_mem P)

/-- If `P` holds on `(ξ₁,1)`, the threshold is `≤ ξ₁`. -/
theorem thr_le {ξ₁ : ℝ} (h0 : 0 ≤ ξ₁) (h1 : ξ₁ ≤ 1) (hP : ∀ ξ ∈ Ioo ξ₁ 1, P ξ) :
    sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P ξ} ≤ ξ₁ :=
  csInf_le (thr_bdd P) ⟨⟨h0, h1⟩, hP⟩

/-- If `P` fails at `ξ₁ < 1`, the threshold is `≥ ξ₁`. -/
theorem le_thr {ξ₁ : ℝ} (h1 : ξ₁ < 1) (hn : ¬ P ξ₁) :
    ξ₁ ≤ sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P ξ} := by
  refine le_csInf ⟨1, thr_one_mem P⟩ fun ξ₀ h => ?_
  by_contra hlt
  exact hn (h.2 ξ₁ ⟨lt_of_not_ge hlt, h1⟩)

/-- `P` holds strictly above the threshold. -/
theorem holds_of_thr_lt {ξ : ℝ}
    (h : sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P ξ} < ξ) (h1 : ξ < 1) : P ξ := by
  obtain ⟨ξ₀, hξ₀, hlt⟩ := exists_lt_of_csInf_lt ⟨1, thr_one_mem P⟩ h
  exact hξ₀.2 ξ ⟨hlt, h1⟩

/-- Monotonicity of the threshold in the predicate. -/
theorem thr_mono (P' : ℝ → Prop) (hPP : ∀ ξ ∈ Ioo (0 : ℝ) 1, P' ξ → P ξ) :
    sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P ξ}
      ≤ sInf {ξ₀ | ξ₀ ∈ Icc (0 : ℝ) 1 ∧ ∀ ξ ∈ Ioo ξ₀ 1, P' ξ} := by
  refine csInf_le_csInf (thr_bdd P) ⟨1, thr_one_mem P'⟩ fun ξ₀ h => ⟨h.1, fun ξ hξ => ?_⟩
  exact hPP ξ ⟨lt_of_le_of_lt h.1.1 hξ.1, hξ.2⟩ (h.2 ξ hξ)

end Generic

/-! ### The radial variable -/

/-- Every `ξ ∈ (0,1)` is `xiOf N τ` for `τ = −(N/2) ln ξ > 0`. -/
theorem exists_tau {N : ℕ} {ξ : ℝ} (hN : 1 ≤ N) (h0 : 0 < ξ) (h1 : ξ < 1) :
    ∃ τ : ℝ, 0 < τ ∧ xiOf N τ = ξ := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hlog : log ξ < 0 := log_neg h0 h1
  refine ⟨-((N : ℝ) / 2) * log ξ, by nlinarith, ?_⟩
  unfold xiOf
  rw [show -2 * (-((N : ℝ) / 2) * log ξ) / N = log ξ by field_simp, exp_log h0]

/-- `ψ(x,τ) ≤ ψ(0,τ)` for `|x| ≤ 1` (note Lemma 13). -/
theorem psi_le_psi_zero {x τ : ℝ} (hx : |x| ≤ 1) (hτ : 0 < τ) : psi x τ ≤ psi 0 τ := by
  have hab : psi x τ = psi |x| τ := by
    rcases le_or_gt 0 x with h | h
    · rw [abs_of_nonneg h]
    · rw [abs_of_neg h, psi_even]
  rw [hab]
  exact (psi_strictAntiOn hτ).antitoneOn ⟨le_rfl, zero_le_one⟩ ⟨abs_nonneg x, hx⟩ (abs_nonneg x)

/-! ### note Cor 14 -/

/-- note Cor 14: for every `N ≥ 2` and every radius (`0 < ξ < 1`), the ring is spectrally stable
iff `a_{⌊N/2⌋} ≥ 0`. -/
theorem specStable_iff {N : ℕ} {ξ : ℝ} (hN : 2 ≤ N) (h0 : 0 < ξ) (h1 : ξ < 1) :
    specStable N ξ ↔ 0 ≤ aBlock N (N / 2) ξ := by
  constructor
  · intro h
    have hab := h (N / 2) (by omega) (by omega)
    have hb := bBlock_pos (m := N / 2) hN (by omega) (by omega) h0 h1
    by_contra hneg
    nlinarith [not_le.mp hneg]
  · intro h m hm1 hmN
    exact mul_nonneg (le_trans h (aBlock_half_le hN hm1 hmN h0 h1))
      (bBlock_pos hN hm1 hmN h0 h1).le

/-! ### note Lemma 17 -/

/-- note Lemma 17 (last display): `a_m b_m − g_m² = N⁴ β(x,τ) [Φ_N(τ) + G(x,τ)]` with
`ξ = xiOf N τ`, `x = xOf N m`, for `1 ≤ m ≤ N−1`, `τ > 0`. -/
theorem det_xiOf {N m : ℕ} (hm1 : 1 ≤ m) (hmN : m < N) {τ : ℝ} (hτ : 0 < τ) :
    aBlock N m (xiOf N τ) * bBlock N m (xiOf N τ) - gBlock N m (xiOf N τ) ^ 2
      = (N : ℝ) ^ 4 * beta (xOf N m) τ * (PhiN N τ + G (xOf N m) τ) := by
  rw [aBlock_xiOf hm1 hmN hτ, bBlock_xiOf hm1 hmN hτ, gBlock_xiOf hm1 hmN hτ]
  have hβ : beta (xOf N m) τ ≠ 0 := (beta_pos (abs_xOf_lt_one hm1 hmN) hτ).ne'
  unfold PhiN G
  field_simp
  ring

/-- note Lemma 17: `G(x,τ) ≥ 0` for `|x| < 1`, `τ > 0`. -/
theorem G_nonneg {x τ : ℝ} (hx : |x| < 1) (hτ : 0 < τ) : 0 ≤ G x τ := by
  rcases eq_or_ne x 0 with rfl | hx0
  · simp [G, gamma, R]
  · have hβ := beta_pos hx hτ
    have hK := K_pos (abs_pos.mpr hx0) hx hτ
    have hs : 0 < sinh τ := sinh_pos_iff.mpr hτ
    have hprod : 0 < beta x τ * G x τ := by
      rw [beta_mul_G hτ hβ.ne']
      positivity
    by_contra hneg
    nlinarith [not_le.mp hneg]

/-- note Lemma 17 / Thm 18: the Hessian is positive definite (every block) wherever `Φ_N(τ) > 0`. -/
theorem emDefinite_of_PhiN_pos {N : ℕ} {τ : ℝ} (hN : 2 ≤ N) (hτ : 0 < τ)
    (hΦ : 0 < PhiN N τ) : emDefinite N (xiOf N τ) := by
  intro m hm1 hmN
  have hx := abs_xOf_lt_one hm1 hmN
  have hN2 : (0 : ℝ) < (N : ℝ) ^ 2 := by
    have : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
    positivity
  refine ⟨?_, ?_⟩
  · rw [aBlock_xiOf hm1 hmN hτ]
    have := psi_le_psi_zero hx.le hτ
    unfold PhiN at hΦ
    have : 0 < EN N τ - 1 / (4 * sinh τ ^ 2) - psi (xOf N m) τ := by linarith
    positivity
  · rw [det_xiOf hm1 hmN hτ]
    have := G_nonneg hx hτ
    have hβ := beta_pos hx hτ
    have : 0 < PhiN N τ + G (xOf N m) τ := by linarith
    positivity

/-- Positive definiteness implies spectral stability (`b_m > 0`, note Lemma 12). -/
theorem specStable_of_emDefinite {N : ℕ} {ξ : ℝ} (hN : 2 ≤ N) (h0 : 0 < ξ) (h1 : ξ < 1)
    (h : emDefinite N ξ) : specStable N ξ := fun m hm1 hmN =>
  (mul_pos (h m hm1 hmN).1 (bBlock_pos hN hm1 hmN h0 h1)).le

/-- The threshold for positive definiteness is at least the spectral threshold. -/
theorem xiStar_le_xiEM {N : ℕ} (hN : 2 ≤ N) : xiStar N ≤ xiEM N :=
  thr_mono (specStable N) (emDefinite N) fun _ hξ h => specStable_of_emDefinite hN hξ.1 hξ.2 h

/-! ### note Thm 18(a) -/

/-- For even `N`, the critical mode `N/2` has `x = 0`. -/
theorem xOf_half_even {N : ℕ} (hN : 1 ≤ N) (hE : Even N) : xOf N (N / 2) = 0 := by
  obtain ⟨n, rfl⟩ := hE
  have hn : (n + n) / 2 = n := by omega
  rw [hn]
  unfold xOf
  have : (0 : ℝ) < n := by exact_mod_cast (show 0 < n by omega)
  push_cast
  field_simp
  ring

/-- For even `N ≥ 2`, `a_{N/2} = N² Φ_N(τ)` at `ξ = xiOf N τ` (note §"The energy–momentum
threshold"). -/
theorem aBlock_half_even {N : ℕ} (hN : 2 ≤ N) (hE : Even N) {τ : ℝ} (hτ : 0 < τ) :
    aBlock N (N / 2) (xiOf N τ) = (N : ℝ) ^ 2 * PhiN N τ := by
  rw [aBlock_xiOf (by omega) (by omega) hτ, xOf_half_even (by omega) hE]
  rfl

/-- note Thm 18(a): for every even `N ≥ 2` and every radius (`0 < ξ < 1`), the Hessian is
positive definite on the non-rotational modes iff `a_{N/2} > 0`. -/
theorem emDefinite_iff_even {N : ℕ} {ξ : ℝ} (hN : 2 ≤ N) (hE : Even N) (h0 : 0 < ξ)
    (h1 : ξ < 1) : emDefinite N ξ ↔ 0 < aBlock N (N / 2) ξ := by
  constructor
  · intro h
    exact (h (N / 2) (by omega) (by omega)).1
  · intro h
    obtain ⟨τ, hτ, rfl⟩ := exists_tau (N := N) (by omega) h0 h1
    rw [aBlock_half_even hN hE hτ] at h
    have hN2 : (0 : ℝ) < (N : ℝ) ^ 2 := by
      have : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
      positivity
    exact emDefinite_of_PhiN_pos hN hτ (pos_of_mul_pos_right h hN2.le)

/-- note Thm 18(a), threshold equality, under the hypothesis that `a_{N/2} ≥ 0` at some radius
forces `a_{N/2} > 0` at every larger `ξ` (this is supplied by the monotonicity of note Prop 15,
see `xiEM_eq_xiStar_of_strictMonoOn`). For even `N ≥ 2`: `ξ_EM(N) = ξ_*(N)`. -/
theorem xiEM_eq_xiStar_even {N : ℕ} (hN : 2 ≤ N) (hE : Even N)
    (hH : ∀ ξ₁ ξ₂ : ℝ, 0 < ξ₁ → ξ₁ < ξ₂ → ξ₂ < 1 →
      0 ≤ aBlock N (N / 2) ξ₁ → 0 < aBlock N (N / 2) ξ₂) :
    xiEM N = xiStar N := by
  unfold xiEM xiStar
  congr 1
  ext ξ₀
  constructor
  · rintro ⟨hI, hP⟩
    exact ⟨hI, fun ξ hξ =>
      specStable_of_emDefinite hN (lt_of_le_of_lt hI.1 hξ.1) hξ.2 (hP ξ hξ)⟩
  · rintro ⟨hI, hP⟩
    refine ⟨hI, fun ξ hξ => ?_⟩
    have hξ0 : 0 < ξ := lt_of_le_of_lt hI.1 hξ.1
    rw [emDefinite_iff_even hN hE hξ0 hξ.2]
    set ξ' := (ξ₀ + ξ) / 2
    have h1' : ξ₀ < ξ' := by simp only [ξ']; linarith [hξ.1]
    have h2' : ξ' < ξ := by simp only [ξ']; linarith [hξ.1]
    have h0' : 0 < ξ' := lt_of_le_of_lt hI.1 h1'
    have := (specStable_iff hN h0' (h2'.trans hξ.2)).mp (hP ξ' ⟨h1', h2'.trans hξ.2⟩)
    exact hH ξ' ξ h0' h2' hξ.2 this

/-- note Thm 18(a): for even `N ≥ 2`, if `ξ ↦ a_{N/2}(ξ)` is strictly increasing on `(0,1)`
(note Prop 15, `N ≥ 4`), then `ξ_EM(N) = ξ_*(N)`: no Krein window. -/
theorem xiEM_eq_xiStar_of_strictMonoOn {N : ℕ} (hN : 2 ≤ N) (hE : Even N)
    (hmono : StrictMonoOn (fun ξ => aBlock N (N / 2) ξ) (Ioo 0 1)) :
    xiEM N = xiStar N :=
  xiEM_eq_xiStar_even hN hE fun _ _ h0 h12 h2 h =>
    lt_of_le_of_lt h (hmono ⟨h0, h12.trans h2⟩ ⟨h0.trans h12, h2⟩ h12)

end Havelock.Tanh.LargeN
