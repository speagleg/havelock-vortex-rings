/-
Havelock — Conjecture 13, odd `n`: the first- and second-order Fejér sums used for the `A⁴`
coefficient (paper `conj13/paper/main.tex`, proof of Thm odd, "Sums over roots of unity";
`conj13/derivation.txt` §8, item 5: "provable by the second-difference recursion in p").

With `θ_m = 2πm/n`, for every `n ≥ 1` and `0 ≤ q ≤ n`:

* `fejer1`  `Σ_{m=1}^{n−1} cos(qθ_m)/(1 − cos θ_m) = (n²−1)/6 − q(n−q)`;
* `fejer2`  `Σ_{m=1}^{n−1} cos(qθ_m)/(1 − cos θ_m)² = (n⁴ − 30n²q² + 10n² + 60nq³ − 60nq − 30q⁴ + 60q² − 11)/180`.

`fejer2` is proved by the second-difference recursion
`F₂(q+1) − 2F₂(q) + F₂(q−1) = −2F₁(q)`, `F₂(1) = F₂(0) − F₁(0)`, and the normalisation
`Σ_{q=0}^{n−1} F₂(q) = 0` (which gives `F₂(0) = (n²−1)(n²+11)/180`).
-/
import Havelock.Conj13.OddSums

open Real Finset

namespace Havelock.Conj13.Odd4

open Havelock.Conj13.OddSums

/-- `1 − cos(2πm/n) ≠ 0` for `0 < m < n`. -/
theorem one_sub_cos_ne (n m : ℕ) (hm : 0 < m) (hmn : m < n) :
    1 - Real.cos (2 * π * m / n) ≠ 0 := by
  intro h
  have h1 : Real.cos (2 * π * m / n) = 1 := by linarith
  obtain ⟨k, hk⟩ := (Real.cos_eq_one_iff _).mp h1
  have hn : (0 : ℝ) < n := by exact_mod_cast (lt_trans hm hmn)
  have hk' : (k : ℝ) * n = m := by
    field_simp at hk
    nlinarith [Real.pi_pos]
  have hz : (k * n : ℤ) = m := by exact_mod_cast hk'
  have hk0 : 0 < k := by
    by_contra hc
    push Not at hc
    have : (k : ℤ) * n ≤ 0 := Int.mul_nonpos_of_nonpos_of_nonneg hc (by positivity)
    omega
  have : (n : ℤ) ≤ k * n := le_mul_of_one_le_left (by positivity) hk0
  omega

theorem one_sub_cos_pos (n m : ℕ) (hm : 0 < m) (hmn : m < n) :
    0 < 1 - Real.cos (2 * π * m / n) :=
  lt_of_le_of_ne (by linarith [Real.cos_le_one (2 * π * m / n)]) (one_sub_cos_ne n m hm hmn).symm

/-- **First-order Fejér sum.** `Σ_{m=1}^{n−1} cos(2πqm/n)/(1 − cos(2πm/n)) = (n²−1)/6 − q(n−q)`. -/
theorem fejer1 (n q : ℕ) (hn : 0 < n) (hq : q ≤ n) :
    ∑ m ∈ Ico 1 n, Real.cos (2 * π * q * m / n) / (1 - Real.cos (2 * π * m / n))
      = ((n : ℝ) ^ 2 - 1) / 6 - (q : ℝ) * ((n : ℝ) - q) := by
  have hF := fejer n q hn hq
  have hT := sum_inv_one_sub_cos n hn
  have e : ∀ m ∈ Ico 1 n, Real.cos (2 * π * q * m / n) / (1 - Real.cos (2 * π * m / n))
      = 1 / (1 - Real.cos (2 * π * m / n))
        - (1 - Real.cos (2 * π * q * m / n)) / (1 - Real.cos (2 * π * m / n)) := by
    intro m _; ring
  rw [Finset.sum_congr rfl e, Finset.sum_sub_distrib, hT, hF]

/-- The second-order Fejér polynomial. -/
noncomputable def fejer2Poly (n q : ℝ) : ℝ :=
  (n ^ 4 - 30 * n ^ 2 * q ^ 2 + 10 * n ^ 2 + 60 * n * q ^ 3 - 60 * n * q - 30 * q ^ 4 + 60 * q ^ 2 - 11)
    / 180

/-- The second-order sum as a function of the frequency. -/
noncomputable def F2 (n q : ℕ) : ℝ :=
  ∑ m ∈ Ico 1 n, Real.cos (2 * π * q * m / n) / (1 - Real.cos (2 * π * m / n)) ^ 2

theorem F2_one (n : ℕ) (hn : 0 < n) :
    F2 n 1 = F2 n 0 - (((n : ℝ) ^ 2 - 1) / 6) := by
  have h1 := fejer1 n 0 hn (Nat.zero_le _)
  simp only [Nat.cast_zero, mul_zero, zero_mul, zero_div, Real.cos_zero, sub_zero] at h1
  unfold F2
  rw [← h1, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun m hm => ?_)
  have hm' := Finset.mem_Ico.mp hm
  have hne := one_sub_cos_ne n m (by omega) hm'.2
  simp only [Nat.cast_one, mul_one, Nat.cast_zero, mul_zero, zero_mul, zero_div, Real.cos_zero]
  field_simp
  ring

/-- The second-difference recursion. -/
theorem F2_rec (n q : ℕ) (hn : 0 < n) (hq : q + 1 ≤ n) :
    F2 n (q + 2) = 2 * F2 n (q + 1) - F2 n q
      - 2 * (((n : ℝ) ^ 2 - 1) / 6 - ((q + 1 : ℕ) : ℝ) * ((n : ℝ) - (q + 1 : ℕ))) := by
  have h1 := fejer1 n (q + 1) hn hq
  unfold F2
  rw [← h1, Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun m hm => ?_)
  have hm' := Finset.mem_Ico.mp hm
  have hne := one_sub_cos_ne n m (by omega) hm'.2
  set x := 2 * π * m / n
  have hx : ∀ k : ℕ, 2 * π * (k : ℝ) * m / n = k * x := by intro k; simp only [x]; ring
  simp only [hx]
  push_cast
  have e2 : ((q : ℝ) + 2) * x = ((q : ℝ) + 1) * x + x := by ring
  have e0 : (q : ℝ) * x = ((q : ℝ) + 1) * x - x := by ring
  rw [e2, e0, Real.cos_add, Real.cos_sub]
  field_simp
  ring

theorem fejer2Poly_rec (n q : ℝ) :
    fejer2Poly n (q + 2) = 2 * fejer2Poly n (q + 1) - fejer2Poly n q
      - 2 * ((n ^ 2 - 1) / 6 - (q + 1) * (n - (q + 1))) := by
  unfold fejer2Poly; ring

/-- `F₂(q) − F₂(0) = P(q) − P(0)` for `q ≤ n`. -/
theorem F2_sub (n : ℕ) (hn : 0 < n) :
    ∀ q, q + 1 ≤ n → F2 n q - F2 n 0 = fejer2Poly n q - fejer2Poly n 0
      ∧ F2 n (q + 1) - F2 n 0 = fejer2Poly n (q + 1) - fejer2Poly n 0 := by
  intro q
  induction q with
  | zero =>
    intro _
    refine ⟨by simp, ?_⟩
    rw [F2_one n hn]
    unfold fejer2Poly; push_cast; ring
  | succ k ih =>
    intro hk
    obtain ⟨-, h2⟩ := ih (by omega)
    refine ⟨by push_cast at h2 ⊢; exact h2, ?_⟩
    have hr := F2_rec n k hn (by omega)
    have hp := fejer2Poly_rec (n : ℝ) (k : ℝ)
    obtain ⟨h1, -⟩ := ih (by omega)
    rw [show k + 1 + 1 = k + 2 by ring, hr]
    push_cast at hp h1 h2 ⊢
    rw [show (k : ℝ) + 2 = (k : ℝ) + 1 + 1 by ring] at hp
    linarith

/-- `Σ_{q<n} cos(2πqm/n) = 0` for `0 < m < n`. -/
theorem sum_q_cos (n m : ℕ) (hm : 0 < m) (hmn : m < n) :
    ∑ q ∈ range n, Real.cos (2 * π * q * m / n) = 0 := by
  have h := sum_cos_eq_zero n m (by omega) (fun h => by have := Nat.le_of_dvd hm h; omega)
  rw [← h]
  refine Finset.sum_congr rfl (fun q _ => by ring_nf)

theorem sum_range_cube_real (n : ℕ) :
    ∑ p ∈ range n, (p : ℝ) ^ 3 = ((n : ℝ) * (n - 1) / 2) ^ 2 := by
  induction n with
  | zero => simp
  | succ k ih => rw [sum_range_succ, ih]; push_cast; ring

theorem sum_range_four_real (n : ℕ) :
    ∑ p ∈ range n, (p : ℝ) ^ 4 = (n : ℝ) * (n - 1) * (2 * n - 1) * (3 * n ^ 2 - 3 * n - 1) / 30 := by
  induction n with
  | zero => simp
  | succ k ih => rw [sum_range_succ, ih]; push_cast; ring

/-- `F₂(0) = Σ 1/(1−cos θ_m)² = (n²−1)(n²+11)/180`. -/
theorem F2_zero (n : ℕ) (hn : 0 < n) : F2 n 0 = fejer2Poly n 0 := by
  have hsum : ∑ q ∈ range n, F2 n q = 0 := by
    unfold F2
    rw [Finset.sum_comm]
    refine Finset.sum_eq_zero (fun m hm => ?_)
    have hm' := Finset.mem_Ico.mp hm
    rw [← Finset.sum_div, sum_q_cos n m (by omega) hm'.2, zero_div]
  have hq : ∀ q ∈ range n, F2 n q = F2 n 0 + (fejer2Poly n q - fejer2Poly n 0) := by
    intro q hq
    have := (F2_sub n hn q (by have := Finset.mem_range.mp hq; omega)).1
    linarith
  rw [Finset.sum_congr rfl hq, Finset.sum_add_distrib, Finset.sum_const, card_range,
    nsmul_eq_mul] at hsum
  have hP : ∑ q ∈ range n, (fejer2Poly n q - fejer2Poly n 0)
      = -((n : ℝ) * ((n : ℝ) ^ 4 + 10 * (n : ℝ) ^ 2 - 11)) / 180 := by
    have e : ∀ q ∈ range n, fejer2Poly n q - fejer2Poly n 0
        = (-30 * (n : ℝ) ^ 2 + 60) / 180 * (q : ℝ) ^ 2 + 60 * (n : ℝ) / 180 * (q : ℝ) ^ 3
          - 60 * (n : ℝ) / 180 * (q : ℝ) - 30 / 180 * (q : ℝ) ^ 4 := by
      intro q _; unfold fejer2Poly; ring
    rw [Finset.sum_congr rfl e, Finset.sum_sub_distrib, Finset.sum_sub_distrib,
      Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum,
      sum_range_sq_real, sum_range_cube_real, sum_range_id_real, sum_range_four_real]
    ring
  rw [hP] at hsum
  have hn' : (n : ℝ) ≠ 0 := by positivity
  have : F2 n 0 = ((n : ℝ) ^ 4 + 10 * (n : ℝ) ^ 2 - 11) / 180 := by
    have h2 : (n : ℝ) * (F2 n 0 * 180 - ((n : ℝ) ^ 4 + 10 * (n : ℝ) ^ 2 - 11)) = 0 := by linarith
    rcases mul_eq_zero.mp h2 with h | h
    · exact absurd h hn'
    · linarith
  rw [this]; unfold fejer2Poly; ring

/-- **Second-order Fejér sum** (`conj13/derivation.txt` §8): for `n ≥ 1`, `0 ≤ q ≤ n`,
`Σ_{m=1}^{n−1} cos(2πqm/n)/(1 − cos(2πm/n))² = (n⁴ − 30n²q² + 10n² + 60nq³ − 60nq − 30q⁴ + 60q² − 11)/180`. -/
theorem fejer2 (n q : ℕ) (hn : 0 < n) (hq : q ≤ n) :
    ∑ m ∈ Ico 1 n, Real.cos (2 * π * q * m / n) / (1 - Real.cos (2 * π * m / n)) ^ 2
      = fejer2Poly n q := by
  have h0 := F2_zero n hn
  rcases Nat.eq_zero_or_pos q with rfl | hq0
  · change F2 n 0 = _
    simpa using h0
  · have := (F2_sub n hn (q - 1) (by omega)).2
    rw [show q - 1 + 1 = q by omega] at this
    have e : ((q - 1 : ℕ) : ℝ) + 1 = q := by rw [Nat.cast_sub (by omega)]; push_cast; ring
    rw [e] at this
    change F2 n q = _
    linarith

end Havelock.Conj13.Odd4
