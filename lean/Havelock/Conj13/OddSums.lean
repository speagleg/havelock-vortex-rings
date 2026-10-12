/-
Havelock — Conjecture 13, odd `n = 2M+1`: the trigonometric pair sums behind the kernel coefficient
`f₂` (paper `conj13/paper/main.tex`, proof of Thm odd), for every odd `n ≥ 3`.

Notation: `φ_k = 2πk/n`, `c_k = cos(Mφ_k)` (`M = l = (n−1)/2`, the critical mode), `C_ij = cos(φ_j − φ_i)`.
All sums are over ordered pairs `i ≠ j` in `range n`.

* `sum_inv_one_sub_cos`  `Σ_{m=1}^{n−1} 1/(1 − cos(2πm/n)) = (n²−1)/6`   (from paper eq. (fejer)).
* `reindex_pairs`        `Σ_i Σ_{j≠i} G_i(j − i) = Σ_i Σ_{m=1}^{n−1} G_i(m)` for `n`-periodic `G_i`.
* `pairA1`  `Σ_{i≠j} c_i c_j/(1 − C_ij) = −n(n²−1)/24`.
* `pairA2`  `Σ_{i≠j} (c_i² + c_j²)/(1 − C_ij) = n(n²−1)/6`.
* `pairB2`  `Σ_{i≠j} (c_i² + c_j²) = n(n−1)`;   `pairB0`  `Σ_{i≠j} (c_i + c_j) = 0`.
-/
import Havelock.Conj13Algebra

open Real Finset

namespace Havelock.Conj13.OddSums

theorem sum_sin_eq_zero (n j : ℕ) (hn : 0 < n) (hj : ¬ n ∣ j) :
    ∑ k ∈ range n, Real.sin (2 * π * j * k / n) = 0 := by
  have hζ : IsPrimitiveRoot (Complex.exp (2 * π * Complex.I / n)) n :=
    Complex.isPrimitiveRoot_exp n (by omega)
  have h := sum_zpow_geom hζ (j : ℤ) (by exact_mod_cast hj)
  have h2 : ∀ k : ℕ, (Complex.exp (2 * π * Complex.I / n) ^ (j : ℤ)) ^ k
      = Complex.exp (((2 * π * j * k / n : ℝ)) * Complex.I) := by
    intro k
    rw [zpow_natCast, ← pow_mul, ← Complex.exp_nat_mul]
    congr 1; push_cast; ring
  simp_rw [h2] at h
  have := congrArg Complex.im h
  rw [Complex.im_sum] at this
  simp only [Complex.exp_ofReal_mul_I_im] at this
  exact this

/-- `Σ_k cos(2πtk/n + a) = 0` for `n ∤ t`. -/
theorem sum_cos_add_eq_zero (n t : ℕ) (hn : 0 < n) (ht : ¬ n ∣ t) (a : ℝ) :
    ∑ k ∈ range n, Real.cos (2 * π * t * k / n + a) = 0 := by
  simp_rw [Real.cos_add]
  rw [Finset.sum_sub_distrib, ← Finset.sum_mul, ← Finset.sum_mul, sum_cos_eq_zero n t hn ht,
    sum_sin_eq_zero n t hn ht]
  ring

theorem sum_range_id_real (n : ℕ) : ∑ p ∈ range n, (p : ℝ) = n * (n - 1) / 2 := by
  induction n with
  | zero => simp
  | succ k ih => rw [sum_range_succ, ih]; push_cast; ring

theorem sum_range_sq_real (n : ℕ) : ∑ p ∈ range n, (p : ℝ) ^ 2 = n * (n - 1) * (2 * n - 1) / 6 := by
  induction n with
  | zero => simp
  | succ k ih => rw [sum_range_succ, ih]; push_cast; ring

/-- `Σ_{m=1}^{n−1} 1/(1 − cos(2πm/n)) = (n²−1)/6`: sum paper eq. (fejer) over `0 ≤ p < n`. -/
theorem sum_inv_one_sub_cos (n : ℕ) (hn : 0 < n) :
    ∑ m ∈ Ico 1 n, 1 / (1 - Real.cos (2 * π * m / n)) = ((n : ℝ) ^ 2 - 1) / 6 := by
  have hF : ∑ p ∈ range n, ∑ m ∈ Ico 1 n,
      (1 - Real.cos (2 * π * p * m / n)) / (1 - Real.cos (2 * π * m / n))
        = ∑ p ∈ range n, (p : ℝ) * ((n : ℝ) - p) :=
    Finset.sum_congr rfl (fun p hp => fejer n p hn (Finset.mem_range.mp hp).le)
  rw [Finset.sum_comm] at hF
  have hinner : ∀ m ∈ Ico 1 n, ∑ p ∈ range n,
      (1 - Real.cos (2 * π * p * m / n)) / (1 - Real.cos (2 * π * m / n))
        = n * (1 / (1 - Real.cos (2 * π * m / n))) := by
    intro m hm
    have hm' := Finset.mem_Ico.mp hm
    have hz := sum_cos_eq_zero n m hn (fun h => by have := Nat.le_of_dvd (by omega) h; omega)
    simp_rw [sub_div, Finset.sum_sub_distrib, ← Finset.sum_div]
    have : ∑ p ∈ range n, Real.cos (2 * π * p * m / n) = 0 := by
      rw [← hz]; refine Finset.sum_congr rfl (fun p _ => by ring_nf)
    rw [this]; simp; ring
  rw [Finset.sum_congr rfl hinner, ← Finset.mul_sum] at hF
  have hR : ∑ p ∈ range n, (p : ℝ) * ((n : ℝ) - p) = n * ((n : ℝ) ^ 2 - 1) / 6 := by
    have : ∑ p ∈ range n, (p : ℝ) * ((n : ℝ) - p)
        = n * ∑ p ∈ range n, (p : ℝ) - ∑ p ∈ range n, (p : ℝ) ^ 2 := by
      rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
      exact Finset.sum_congr rfl (fun p _ => by ring)
    rw [this, sum_range_id_real, sum_range_sq_real]; ring
  rw [hR] at hF
  have hn' : (n : ℝ) ≠ 0 := by positivity
  field_simp at hF ⊢
  linarith

/-- Reindexing the pairs `(i, j)`, `j ≠ i`, by the offset `m = j − i mod n`. -/
theorem reindex_row (n : ℕ) (G : ℝ → ℝ) (hG : ∀ x, G (x - n) = G x) (i : ℕ) (hi : i < n) :
    ∑ j ∈ (range n).erase i, G ((j : ℝ) - i) = ∑ m ∈ Ico 1 n, G m := by
  refine Finset.sum_nbij' (fun j => if i < j then j - i else j + n - i)
    (fun m => if i + m < n then i + m else i + m - n) ?_ ?_ ?_ ?_ ?_
  · intro j hj
    simp only [Finset.mem_erase, Finset.mem_range] at hj
    simp only [Finset.mem_Ico]
    split_ifs <;> omega
  · intro m hm
    simp only [Finset.mem_Ico] at hm
    simp only [Finset.mem_erase, Finset.mem_range]
    split_ifs <;> omega
  · intro j hj
    simp only [Finset.mem_erase, Finset.mem_range] at hj
    split_ifs <;> omega
  · intro m hm
    simp only [Finset.mem_Ico] at hm
    split_ifs <;> omega
  · intro j hj
    simp only [Finset.mem_erase, Finset.mem_range] at hj
    split_ifs with h
    · rw [Nat.cast_sub h.le]
    · have h' : i ≤ j + n := by omega
      rw [Nat.cast_sub h', Nat.cast_add, ← hG ((j : ℝ) + n - i)]
      congr 1; ring

/-- The pair sum with the diagonal removed, reindexed by offsets. -/
theorem reindex_pairs (n : ℕ) (G : ℕ → ℝ → ℝ) (hG : ∀ i x, G i (x - n) = G i x) :
    ∑ i ∈ range n, ∑ j ∈ range n, (if i = j then 0 else G i ((j : ℝ) - i))
      = ∑ i ∈ range n, ∑ m ∈ Ico 1 n, G i m := by
  refine Finset.sum_congr rfl (fun i hi => ?_)
  have hi' := Finset.mem_range.mp hi
  rw [← Finset.sum_erase_add _ _ hi, if_pos rfl, add_zero, ← reindex_row n (G i) (hG i) i hi']
  refine Finset.sum_congr rfl (fun j hj => ?_)
  rw [if_neg (Finset.ne_of_mem_erase hj).symm]

/-! ### The odd-`n` data -/

variable (M : ℕ)

/-- `c_k = cos(Mφ_k)`, `φ_k = 2πk/n`, `n = 2M+1`. -/
noncomputable def cc (k : ℝ) : ℝ := Real.cos (2 * π * M * k / (2 * M + 1 : ℕ))

/-- `cos(2π x/n)`. -/
noncomputable def CC (x : ℝ) : ℝ := Real.cos (2 * π * x / (2 * M + 1 : ℕ))

theorem cc_per (k : ℝ) : cc M (k - (2 * M + 1 : ℕ)) = cc M k := by
  unfold cc
  have hn : ((2 * M + 1 : ℕ) : ℝ) ≠ 0 := by positivity
  rw [show 2 * π * M * (k - (2 * M + 1 : ℕ)) / (2 * M + 1 : ℕ) = 2 * π * M * k / (2 * M + 1 : ℕ)
      - (M : ℤ) * (2 * π) by push_cast; field_simp]
  exact Real.cos_sub_int_mul_two_pi _ _

theorem CC_per (x : ℝ) : CC M (x - (2 * M + 1 : ℕ)) = CC M x := by
  unfold CC
  have hn : ((2 * M + 1 : ℕ) : ℝ) ≠ 0 := by positivity
  rw [show 2 * π * (x - (2 * M + 1 : ℕ)) / (2 * M + 1 : ℕ) = 2 * π * x / (2 * M + 1 : ℕ)
      - ((1 : ℤ) : ℝ) * (2 * π) by push_cast; field_simp]
  exact Real.cos_sub_int_mul_two_pi _ _

theorem CC_neg (x : ℝ) : CC M (-x) = CC M x := by
  unfold CC; rw [show 2 * π * -x / (2 * M + 1 : ℕ) = -(2 * π * x / (2 * M + 1 : ℕ)) by ring,
    Real.cos_neg]

theorem not_dvd_M (hM : 1 ≤ M) : ¬ (2 * M + 1) ∣ M := fun h => by
  have := Nat.le_of_dvd (by omega) h; omega

theorem not_dvd_2M (hM : 1 ≤ M) : ¬ (2 * M + 1) ∣ 2 * M := fun h => by
  have := Nat.le_of_dvd (by omega) h; omega

theorem sum_cc (hM : 1 ≤ M) : ∑ k ∈ range (2 * M + 1), cc M k = 0 := by
  have := sum_cos_eq_zero (2 * M + 1) M (by omega) (not_dvd_M M hM)
  unfold cc; push_cast at this ⊢; exact this

theorem sum_cc_sq (hM : 1 ≤ M) : ∑ k ∈ range (2 * M + 1), cc M k ^ 2 = (2 * M + 1 : ℕ) / 2 := by
  have h := sum_cos_eq_zero (2 * M + 1) (2 * M) (by omega) (not_dvd_2M M hM)
  have e : ∀ k ∈ range (2 * M + 1), cc M k ^ 2
      = 1 / 2 + Real.cos (2 * π * ((2 * M : ℕ) : ℝ) * k / ((2 * M + 1 : ℕ) : ℝ)) / 2 := by
    intro k _
    unfold cc
    rw [Real.cos_sq]
    congr 2; push_cast; ring_nf
  rw [Finset.sum_congr rfl e, Finset.sum_add_distrib, ← Finset.sum_div, ← Finset.sum_div, h]
  simp

theorem sum_ite_zero (s : Finset ℕ) (i : ℕ) (hi : i ∈ s) (g : ℕ → ℝ) :
    ∑ j ∈ s, (if i = j then 0 else g j) = ∑ j ∈ s, g j - g i := by
  rw [← Finset.sum_erase_add _ _ hi, ← Finset.sum_erase_add s g hi, if_pos rfl, add_zero,
    add_sub_cancel_right]
  refine Finset.sum_congr rfl (fun j hj => ?_)
  rw [if_neg (Finset.ne_of_mem_erase hj).symm]

/-- `Σ_{i≠j} (f_i + f_j) = 2(n−1) Σ f`. -/
theorem offdiag_add (n : ℕ) (f : ℕ → ℝ) :
    ∑ i ∈ range n, ∑ j ∈ range n, (if i = j then 0 else f i + f j)
      = 2 * ((n : ℝ) - 1) * ∑ k ∈ range n, f k := by
  have h : ∀ i ∈ range n, ∑ j ∈ range n, (if i = j then 0 else f i + f j)
      = n * f i + ∑ k ∈ range n, f k - 2 * f i := by
    intro i hi
    rw [sum_ite_zero _ i hi, Finset.sum_add_distrib, Finset.sum_const, Finset.card_range,
      nsmul_eq_mul]
    ring
  rw [Finset.sum_congr rfl h, Finset.sum_sub_distrib, Finset.sum_add_distrib, ← Finset.mul_sum,
    ← Finset.mul_sum, Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  ring

/-- `Σ_{i≠j} (c_i + c_j) = 0`. -/
theorem pairB0 (hM : 1 ≤ M) :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else cc M i + cc M j) = 0 := by
  rw [offdiag_add (2 * M + 1) (fun k => cc M k), sum_cc M hM, mul_zero]

/-- `Σ_{i≠j} (c_i² + c_j²) = n(n−1)`. -/
theorem pairB2 (hM : 1 ≤ M) :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else cc M i ^ 2 + cc M j ^ 2)
        = (2 * M + 1 : ℕ) * ((2 * M + 1 : ℕ) - 1) := by
  rw [offdiag_add (2 * M + 1) (fun k => cc M k ^ 2), sum_cc_sq M hM]
  ring

theorem T1_CC (hM : 1 ≤ M) :
    ∑ m ∈ Ico 1 (2 * M + 1), 1 / (1 - CC M m) = (((2 * M + 1 : ℕ) : ℝ) ^ 2 - 1) / 6 :=
  sum_inv_one_sub_cos (2 * M + 1) (by omega)

/-- `Σ_{i≠j} c_i²/(1 − C_ij) = (n/2)(n²−1)/6`. -/
theorem pairA2_half (hM : 1 ≤ M) :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else cc M i ^ 2 / (1 - CC M ((j : ℝ) - i)))
        = ((2 * M + 1 : ℕ) : ℝ) / 2 * ((((2 * M + 1 : ℕ) : ℝ) ^ 2 - 1) / 6) := by
  rw [reindex_pairs (2 * M + 1) (fun i x => cc M i ^ 2 / (1 - CC M x))
    (fun i x => by simp only [CC_per])]
  simp_rw [div_eq_mul_inv (cc M _ ^ 2), ← Finset.mul_sum]
  have : ∑ m ∈ Ico 1 (2 * M + 1), (1 - CC M m)⁻¹ = (((2 * M + 1 : ℕ) : ℝ) ^ 2 - 1) / 6 := by
    rw [← T1_CC M hM]; simp
  rw [this, ← Finset.sum_mul, sum_cc_sq M hM]

/-- `Σ_{i≠j} (c_i² + c_j²)/(1 − C_ij) = n(n²−1)/6`. -/
theorem pairA2 (hM : 1 ≤ M) :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (cc M i ^ 2 + cc M j ^ 2) / (1 - CC M ((j : ℝ) - i)))
        = ((2 * M + 1 : ℕ) : ℝ) * ((((2 * M + 1 : ℕ) : ℝ) ^ 2 - 1) / 6) := by
  have hsplit : ∀ i j : ℕ, (if i = j then (0 : ℝ) else (cc M i ^ 2 + cc M j ^ 2) / (1 - CC M ((j : ℝ) - i)))
      = (if i = j then 0 else cc M i ^ 2 / (1 - CC M ((j : ℝ) - i)))
        + (if j = i then 0 else cc M j ^ 2 / (1 - CC M ((i : ℝ) - j))) := by
    intro i j
    by_cases h : i = j
    · simp [h]
    · rw [if_neg h, if_neg (Ne.symm h), if_neg h, add_div,
        show (i : ℝ) - j = -((j : ℝ) - i) by ring, CC_neg]
  simp_rw [hsplit, Finset.sum_add_distrib]
  have hswap : (∑ x ∈ range (2 * M + 1), ∑ y ∈ range (2 * M + 1),
      if y = x then (0 : ℝ) else cc M y ^ 2 / (1 - CC M ((x : ℝ) - y)))
      = ∑ x ∈ range (2 * M + 1), ∑ y ∈ range (2 * M + 1),
      if x = y then (0 : ℝ) else cc M x ^ 2 / (1 - CC M ((y : ℝ) - x)) := Finset.sum_comm
  rw [hswap, pairA2_half M hM]
  ring

/-- `Σ_i c_i c_{i+m} = (n/2) c_m`. -/
theorem sum_cc_mul_shift (hM : 1 ≤ M) (m : ℝ) :
    ∑ i ∈ range (2 * M + 1), cc M i * cc M (i + m) = ((2 * M + 1 : ℕ) : ℝ) / 2 * cc M m := by
  have h := sum_cos_add_eq_zero (2 * M + 1) (2 * M) (by omega) (not_dvd_2M M hM)
    (2 * π * M * m / (2 * M + 1 : ℕ))
  have e : ∀ i ∈ range (2 * M + 1), cc M i * cc M (i + m)
      = cc M m / 2 + Real.cos (2 * π * ((2 * M : ℕ) : ℝ) * i / ((2 * M + 1 : ℕ) : ℝ)
          + 2 * π * M * m / (2 * M + 1 : ℕ)) / 2 := by
    intro i _
    unfold cc
    set a := 2 * π * M * i / ((2 * M + 1 : ℕ) : ℝ)
    set b := 2 * π * M * m / ((2 * M + 1 : ℕ) : ℝ)
    have h1 : 2 * π * M * (i + m) / ((2 * M + 1 : ℕ) : ℝ) = a + b := by simp only [a, b]; ring
    have h2 : 2 * π * ((2 * M : ℕ) : ℝ) * i / ((2 * M + 1 : ℕ) : ℝ) + b = (a + b) + a := by
      simp only [a, b]; push_cast; ring
    rw [h1, h2]
    have e1 := Real.cos_add (a + b) a
    have e2 := Real.cos_sub (a + b) a
    rw [show a + b - a = b by ring] at e2
    rw [e1, e2]
    ring
  rw [Finset.sum_congr rfl e, Finset.sum_add_distrib, ← Finset.sum_div, ← Finset.sum_div, h]
  simp; ring

/-- `Σ_{i≠j} c_i c_j/(1 − C_ij) = −n(n²−1)/24`. -/
theorem pairA1 (hM : 1 ≤ M) :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else cc M i * cc M j / (1 - CC M ((j : ℝ) - i)))
        = -((2 * M + 1 : ℕ) : ℝ) * ((((2 * M + 1 : ℕ) : ℝ) ^ 2 - 1) / 24) := by
  have hrw : ∀ i j : ℕ, (if i = j then (0 : ℝ) else cc M i * cc M j / (1 - CC M ((j : ℝ) - i)))
      = (if i = j then 0 else (fun (i : ℕ) (x : ℝ) => cc M i * cc M (i + x) / (1 - CC M x)) i ((j : ℝ) - i)) := by
    intro i j; simp only [add_sub_cancel]
  simp_rw [hrw]
  rw [reindex_pairs (2 * M + 1) (fun (i : ℕ) (x : ℝ) => cc M i * cc M (i + x) / (1 - CC M x)) (fun i x => by
    simp only [CC_per, show (i : ℝ) + (x - (2 * M + 1 : ℕ)) = (i + x) - (2 * M + 1 : ℕ) by ring, cc_per])]
  rw [Finset.sum_comm]
  simp_rw [← Finset.sum_div, sum_cc_mul_shift M hM]
  have hF := fejer (2 * M + 1) M (by omega) (by omega)
  have hT := T1_CC M hM
  have e : ∀ m ∈ Ico 1 (2 * M + 1), ((2 * M + 1 : ℕ) : ℝ) / 2 * cc M m / (1 - CC M m)
      = ((2 * M + 1 : ℕ) : ℝ) / 2 * (1 / (1 - CC M m)
        - (1 - Real.cos (2 * π * M * m / (2 * M + 1 : ℕ))) / (1 - Real.cos (2 * π * m / (2 * M + 1 : ℕ)))) := by
    intro m _
    unfold cc CC
    ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum, Finset.sum_sub_distrib, hT, hF]
  push_cast
  ring

end Havelock.Conj13.OddSums
