import Havelock.Tanh.LargeN.OddN
import Havelock.Tanh.StableN2to7
import Havelock.Tanh.CertN9
import Havelock.Tanh.CertN11
import Havelock.Tanh.CertN13b
import Havelock.Tanh.CertN15b

/-! # Tanh model: every `N` (note Prop oddN, last sentence, with note Thm cert)

Combines the large-`N` results with the finite-`N` certificates:

* `unique_threshold_of_ge_eight`: for every `N ≥ 8` the critical block `a_{⌊N/2⌋}` has exactly one
  zero `ξ₀ ∈ (0,1)`, is negative on `(0,ξ₀)` and positive on `(ξ₀,1)`; moreover `ξ_*(N) = ξ₀` and,
  for `0 < ξ < 1`, the ring is spectrally stable iff `ξ_*(N) ≤ ξ` (one spectral threshold).
  Sources: even `N ≥ 8`, `aBlock_half_threshold` (note Prop 15); `N ≥ 17`,
  `aBlock_half_threshold_of_ge_17` (note Prop oddN); `N ∈ {9,11,13,15}`, the certificates
  `cert_N9`, `cert_N11`, `cert_N13`, `cert_N15` (note Thm cert(b)), where the zero of the critical
  block at the root `ξs` of `f_N` comes from the numerator identity `N*_linkA*` and the factorisation
  `N*_crit_factor`.
* `stable_all_radii_of_le_seven`: for `2 ≤ N ≤ 7` and every `ξ ∈ (0,1)` the ring is spectrally
  stable and the Hessian is positive definite on the non-rotational modes (`stable_N2..stable_N7`,
  note Thm cert(a)); hence `ξ_*(N) = ξ_EM(N) = 0`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real Set

/-! ### Generic: from the sign pattern of the critical block to the threshold -/

/-- If `a_{⌊N/2⌋}` vanishes at `ξ₀ ∈ (0,1)`, is negative on `(0,ξ₀)` and positive on `(ξ₀,1)`,
then `ξ₀` is its only zero in `(0,1)`. -/
theorem zero_iff_of_sign {N : ℕ} {ξ₀ : ℝ} (hz : aBlock N (N / 2) ξ₀ = 0)
    (hn : ∀ ξ ∈ Ioo 0 ξ₀, aBlock N (N / 2) ξ < 0) (hp : ∀ ξ ∈ Ioo ξ₀ 1, 0 < aBlock N (N / 2) ξ) :
    ∀ ξ ∈ Ioo (0 : ℝ) 1, aBlock N (N / 2) ξ = 0 ↔ ξ = ξ₀ := by
  intro ξ ⟨h0, h1⟩
  refine ⟨fun ha => ?_, fun h => h ▸ hz⟩
  rcases lt_trichotomy ξ ξ₀ with h | h | h
  · exact absurd ha (hn ξ ⟨h0, h⟩).ne
  · exact h
  · exact absurd ha (hp ξ ⟨h, h1⟩).ne'

/-- If `a_{⌊N/2⌋}` vanishes on `(0,1)` only at `ξ₀`, is negative on `(0,ξ₀)` and positive on
`(ξ₀,1)`, then `ξ_*(N) = ξ₀` and `specStable N ξ ↔ ξ₀ ≤ ξ` for `0 < ξ < 1`. -/
theorem xiStar_of_sign {N : ℕ} (hN : 2 ≤ N) {ξ₀ : ℝ} (hξ0 : ξ₀ ∈ Ioo (0 : ℝ) 1)
    (hz : ∀ ξ ∈ Ioo (0 : ℝ) 1, aBlock N (N / 2) ξ = 0 ↔ ξ = ξ₀)
    (hn : ∀ ξ ∈ Ioo 0 ξ₀, aBlock N (N / 2) ξ < 0) (hp : ∀ ξ ∈ Ioo ξ₀ 1, 0 < aBlock N (N / 2) ξ) :
    xiStar N = ξ₀ ∧ ∀ ξ ∈ Ioo (0 : ℝ) 1, specStable N ξ ↔ ξ₀ ≤ ξ := by
  have hle : xiStar N ≤ ξ₀ :=
    thr_le (specStable N) hξ0.1.le hξ0.2.le fun ξ hξ =>
      (specStable_iff hN (hξ0.1.trans hξ.1) hξ.2).mpr (hp ξ hξ).le
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
    have := (specStable_iff hN hξpos hξ1).mp hst
    exact absurd (hn ξ ⟨hξpos, hξb⟩) (not_lt.mpr this)
  refine ⟨le_antisymm hle hge, ?_⟩
  intro ξ ⟨h0, h1⟩
  rw [specStable_iff hN h0 h1]
  constructor
  · intro ha
    by_contra hlt
    exact absurd (hn ξ ⟨h0, not_le.mp hlt⟩) (not_lt.mpr ha)
  · intro hle'
    rcases hle'.eq_or_lt with heq' | hlt
    · exact ((hz ξ ⟨h0, h1⟩).mpr heq'.symm).ge
    · exact (hp ξ ⟨hlt, h1⟩).le

/-- The common shape of the threshold statement. -/
def HasUniqueThreshold (N : ℕ) : Prop :=
  ∃ ξ₀ ∈ Ioo (0 : ℝ) 1, (∀ ξ ∈ Ioo (0 : ℝ) 1, aBlock N (N / 2) ξ = 0 ↔ ξ = ξ₀) ∧
    (∀ ξ ∈ Ioo 0 ξ₀, aBlock N (N / 2) ξ < 0) ∧ (∀ ξ ∈ Ioo ξ₀ 1, 0 < aBlock N (N / 2) ξ)

/-- From a certificate-shaped statement (`a_crit < 0` below `ξs`, `> 0` above, `= 0` at `ξs`). -/
theorem hasUniqueThreshold_of_cert {N : ℕ} {ξs : ℝ} (h0 : 0 < ξs) (h1 : ξs < 1)
    (hz : aBlock N (N / 2) ξs = 0)
    (hn : ∀ x, 0 < x → x < 1 → x < ξs → aBlock N (N / 2) x < 0)
    (hp : ∀ x, 0 < x → x < 1 → ξs < x → 0 < aBlock N (N / 2) x) : HasUniqueThreshold N := by
  have hn' : ∀ ξ ∈ Ioo 0 ξs, aBlock N (N / 2) ξ < 0 := fun ξ h => hn ξ h.1 (h.2.trans h1) h.2
  have hp' : ∀ ξ ∈ Ioo ξs 1, 0 < aBlock N (N / 2) ξ := fun ξ h => hp ξ (h0.trans h.1) h.2 h.1
  exact ⟨ξs, ⟨h0, h1⟩, zero_iff_of_sign hz hn' hp', hn', hp'⟩

/-! ### The certificates for odd `N ∈ {9, 11, 13, 15}` -/

/-- Positivity of the denominator `2 (1−x)² (1−x^N)²` in the numerator identities. -/
theorem denom_pos {N : ℕ} (hN : N ≠ 0) {x : ℝ} (h0 : 0 < x) (h1 : x < 1) :
    0 < 2 * (1 - x) ^ 2 * (1 - x ^ N) ^ 2 := by
  have ha : 0 < 1 - x := sub_pos.mpr h1
  have hb : 0 < 1 - x ^ N := sub_pos.mpr (pow_lt_one₀ h0.le h1 hN)
  positivity

/-- note Thm cert(b), `N = 9`, in the threshold shape. -/
theorem hasUniqueThreshold_nine : HasUniqueThreshold 9 := by
  obtain ⟨ξs, hta, htb, hf0, -, hs⟩ := cert_N9
  have h0 : 0 < ξs := lt_trans (by norm_num) hta
  have h1 : ξs < 1 := htb.trans (by norm_num)
  have hz : aBlock 9 4 ξs = 0 := by
    have hl := N9_linkA4 ξs h0 h1
    rw [N9_crit_factor, hf0] at hl
    exact (mul_eq_zero.mp (hl.trans (by simp))).resolve_right (denom_pos (by norm_num) h0 h1).ne'
  exact hasUniqueThreshold_of_cert (N := 9) h0 h1 hz (fun x hx0 hx1 => (hs x hx0 hx1).2.1)
    (fun x hx0 hx1 hgt => (hs x hx0 hx1).2.2 hgt 4 (by norm_num) (by norm_num))

/-- note Thm cert(b), `N = 11`, in the threshold shape. -/
theorem hasUniqueThreshold_eleven : HasUniqueThreshold 11 := by
  obtain ⟨ξs, hta, htb, hf0, -, hs⟩ := cert_N11
  have h0 : 0 < ξs := lt_trans (by norm_num) hta
  have h1 : ξs < 1 := htb.trans (by norm_num)
  have hz : aBlock 11 5 ξs = 0 := by
    have hl := N11_linkA5 ξs h0 h1
    rw [N11_crit_factor, hf0] at hl
    exact (mul_eq_zero.mp (hl.trans (by simp))).resolve_right (denom_pos (by norm_num) h0 h1).ne'
  exact hasUniqueThreshold_of_cert (N := 11) h0 h1 hz (fun x hx0 hx1 => (hs x hx0 hx1).2.1)
    (fun x hx0 hx1 hgt => (hs x hx0 hx1).2.2 hgt 5 (by norm_num) (by norm_num))

/-- note Thm cert(b), `N = 13`, in the threshold shape. -/
theorem hasUniqueThreshold_thirteen : HasUniqueThreshold 13 := by
  obtain ⟨ξs, hta, htb, hf0, -, hs⟩ := cert_N13
  have h0 : 0 < ξs := lt_trans (by norm_num) hta
  have h1 : ξs < 1 := htb.trans (by norm_num)
  have hz : aBlock 13 6 ξs = 0 := by
    have hl := N13_linkA6 ξs h0 h1
    rw [N13_crit_factor, hf0] at hl
    exact (mul_eq_zero.mp (hl.trans (by simp))).resolve_right (denom_pos (by norm_num) h0 h1).ne'
  exact hasUniqueThreshold_of_cert (N := 13) h0 h1 hz (fun x hx0 hx1 => (hs x hx0 hx1).2.1)
    (fun x hx0 hx1 hgt => (hs x hx0 hx1).2.2 hgt 6 (by norm_num) (by norm_num))

/-- note Thm cert(b), `N = 15`, in the threshold shape. -/
theorem hasUniqueThreshold_fifteen : HasUniqueThreshold 15 := by
  obtain ⟨ξs, hta, htb, hf0, -, hs⟩ := cert_N15
  have h0 : 0 < ξs := lt_trans (by norm_num) hta
  have h1 : ξs < 1 := htb.trans (by norm_num)
  have hz : aBlock 15 7 ξs = 0 := by
    have hl := N15_linkA7 ξs h0 h1
    rw [N15_crit_factor, hf0] at hl
    exact (mul_eq_zero.mp (hl.trans (by simp))).resolve_right (denom_pos (by norm_num) h0 h1).ne'
  exact hasUniqueThreshold_of_cert (N := 15) h0 h1 hz (fun x hx0 hx1 => (hs x hx0 hx1).2.1)
    (fun x hx0 hx1 hgt => (hs x hx0 hx1).2.2 hgt 7 (by norm_num) (by norm_num))

/-! ### Main results -/

/-- note Prop oddN (last sentence) with note Thm cert(b) and Prop 15: every `N ≥ 8` has exactly one
spectral threshold. There is a unique `ξ₀ ∈ (0,1)` with `a_{⌊N/2⌋}(ξ₀) = 0`; `a_{⌊N/2⌋} < 0` on
`(0,ξ₀)` and `> 0` on `(ξ₀,1)`; `ξ_*(N) = ξ₀`; and for `0 < ξ < 1` the ring is spectrally stable iff
`ξ_*(N) ≤ ξ`. -/
theorem unique_threshold_of_ge_eight {N : ℕ} (hN : 8 ≤ N) :
    ∃ ξ₀ ∈ Ioo (0 : ℝ) 1, (∀ ξ ∈ Ioo (0 : ℝ) 1, aBlock N (N / 2) ξ = 0 ↔ ξ = ξ₀) ∧
      (∀ ξ ∈ Ioo 0 ξ₀, aBlock N (N / 2) ξ < 0) ∧ (∀ ξ ∈ Ioo ξ₀ 1, 0 < aBlock N (N / 2) ξ) ∧
      xiStar N = ξ₀ ∧ (∀ ξ ∈ Ioo (0 : ℝ) 1, specStable N ξ ↔ xiStar N ≤ ξ) := by
  have hT : HasUniqueThreshold N := by
    rcases Nat.even_or_odd N with hE | hO
    · exact aBlock_half_threshold hE hN
    · rcases le_or_gt 17 N with h17 | h17
      · exact aBlock_half_threshold_of_ge_17 h17
      · have : N = 9 ∨ N = 11 ∨ N = 13 ∨ N = 15 := by obtain ⟨k, hk⟩ := hO; omega
        rcases this with rfl | rfl | rfl | rfl
        · exact hasUniqueThreshold_nine
        · exact hasUniqueThreshold_eleven
        · exact hasUniqueThreshold_thirteen
        · exact hasUniqueThreshold_fifteen
  obtain ⟨ξ₀, hξ0, hz, hn, hp⟩ := hT
  obtain ⟨hstar, hiff⟩ := xiStar_of_sign (by omega) hξ0 hz hn hp
  exact ⟨ξ₀, hξ0, hz, hn, hp, hstar, hstar ▸ hiff⟩

/-- note Thm cert(a) and Prop oddN (last sentence): for `2 ≤ N ≤ 7` and every radius
(`0 < ξ < 1`) the ring is spectrally stable and the Hessian is positive definite on the
non-rotational modes (`a_m > 0` and `a_m b_m − g_m² > 0` for all `1 ≤ m ≤ N−1`). -/
theorem stable_all_radii_of_le_seven {N : ℕ} (hN2 : 2 ≤ N) (hN7 : N ≤ 7) {ξ : ℝ} (h0 : 0 < ξ)
    (h1 : ξ < 1) : specStable N ξ ∧ emDefinite N ξ := by
  have hem : emDefinite N ξ := by
    intro m hm1 hmN
    interval_cases N
    · exact ⟨(stable_N2 ξ h0 h1 m hm1 (by omega)).1, (stable_N2 ξ h0 h1 m hm1 (by omega)).2.2⟩
    · exact ⟨(stable_N3 ξ h0 h1 m hm1 (by omega)).1, (stable_N3 ξ h0 h1 m hm1 (by omega)).2.2⟩
    · exact ⟨(stable_N4 ξ h0 h1 m hm1 (by omega)).1, (stable_N4 ξ h0 h1 m hm1 (by omega)).2.2⟩
    · exact ⟨(stable_N5 ξ h0 h1 m hm1 (by omega)).1, (stable_N5 ξ h0 h1 m hm1 (by omega)).2.2⟩
    · exact ⟨(stable_N6 ξ h0 h1 m hm1 (by omega)).1, (stable_N6 ξ h0 h1 m hm1 (by omega)).2.2⟩
    · exact ⟨(stable_N7 ξ h0 h1 m hm1 (by omega)).1, (stable_N7 ξ h0 h1 m hm1 (by omega)).2.2⟩
  exact ⟨specStable_of_emDefinite hN2 h0 h1 hem, hem⟩

/-- For `2 ≤ N ≤ 7` both thresholds vanish: `ξ_*(N) = ξ_EM(N) = 0`. -/
theorem xiStar_eq_zero_of_le_seven {N : ℕ} (hN2 : 2 ≤ N) (hN7 : N ≤ 7) :
    xiStar N = 0 ∧ xiEM N = 0 := by
  constructor
  · exact le_antisymm (thr_le (specStable N) le_rfl zero_le_one fun ξ hξ =>
      (stable_all_radii_of_le_seven hN2 hN7 hξ.1 hξ.2).1) (thr_nonneg (specStable N))
  · exact le_antisymm (thr_le (emDefinite N) le_rfl zero_le_one fun ξ hξ =>
      (stable_all_radii_of_le_seven hN2 hN7 hξ.1 hξ.2).2) (thr_nonneg (emDefinite N))

end Havelock.Tanh.LargeN
