/-
Havelock — Conjecture 13, odd `n = 2M+1`: the `A⁴` coefficient of the Hamiltonian along the paper's
odd ansatz (paper `conj13/paper/main.tex`, proof of Thm odd, "Result").

Along `oddAnsatz` (`w_k = w₀ + A cos(lφ_k) + A² p cos φ_k`, `ϑ_k = φ_k + A² s sin φ_k`, `l = M`):

**`Ham_oddAnsatz_jet`**: for every `M ≥ 2` (odd `n ≥ 5`), `λ ∈ (−1, 1)` and all `p, s`,
`H(oddAnsatz A) − H(ring) = f₂ A² + E(p, s) A⁴ + O(A⁵)`, where `f₂ = f2 n λ` and
`E(p,s) = Efun n λ p s = E₀₀ + E₁₀p + E₀₁s + E₂₀p² + E₀₂s²` are exactly the DEFINITIONS of
`Conj13Algebra.lean` (copied from the paper). In particular the `A¹` and `A³` coefficients vanish
(`Ham_oddAnsatz_jet4` states the full fourth-order jet `(0, 0, f₂, 0, E)`).

Proof: `H = −½ Σ_{i≠j} log D_ij`; `log D_ij = core_ij + ½log P_i + ½log P_j` with the pair jet
`Odd4Pair.jet_pair` and the one-body jet `Odd4Pair.jet_logPw`; the coefficient sums are the pair
sums of `Odd4PairSums` and the one-body Fourier sums.
-/
import Havelock.Conj13.Odd4Pair
import Havelock.Conj13.Odd4PairSums

open Real Asymptotics Filter Topology Finset

namespace Havelock.Conj13.Odd4

open Havelock.Conj13.Series Havelock.Conj13.Odd Havelock.Conj13.OddSums Havelock.Conj13.Config

theorem Jet4.eval_zero {f : ℝ → ℝ} {a0 a1 a2 a3 a4 : ℝ} (h : Jet4 f a0 a1 a2 a3 a4) : f 0 = a0 := by
  unfold Jet4 at h
  obtain ⟨C, hC⟩ := h.bound
  have := hC.self_of_nhds
  simp only [mul_zero, add_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, norm_zero,
    Real.norm_eq_abs] at this
  have h0 : |f 0 - a0| = 0 := le_antisymm this (abs_nonneg _)
  linarith [abs_eq_zero.mp h0]

theorem Jet4.ite {f : ℝ → ℝ} {a0 a1 a2 a3 a4 : ℝ} (c : Prop) [Decidable c]
    (h : ¬ c → Jet4 f a0 a1 a2 a3 a4) :
    Jet4 (fun A => if c then 0 else f A) (if c then 0 else a0) (if c then 0 else a1)
      (if c then 0 else a2) (if c then 0 else a3) (if c then 0 else a4) := by
  by_cases hc : c
  · simp only [hc, ite_true]; exact Jet4.const 0
  · simp only [hc, ite_false]; exact h hc

/-! ### The pair data along the ansatz -/

section
variable (M : ℕ) (b p s : ℝ)

/-- The angle `θ_ij = 2π(j−i)/n`. -/
noncomputable def thIJ (i j : ℕ) : ℝ := 2 * π * ((j : ℝ) - i) / (2 * M + 1 : ℕ)

/-- `d_ij = s(S_j − S_i)`. -/
noncomputable def dIJ (i j : ℕ) : ℝ := s * (SS M j - SS M i)

/-- The one-body term `log P(w_k)`. -/
noncomputable def gK (k : ℕ) (A : ℝ) : ℝ :=
  Real.log (Pw b (w0 b + A * cc M k + A ^ 2 * (p * CC M k)))

/-- The pair term `log D_ij` along the ansatz. -/
noncomputable def logD (i j : ℕ) (A : ℝ) : ℝ :=
  Real.log (Dg4 b (cc M i) (p * CC M i) (cc M j) (p * CC M j) (thIJ M i j) (dIJ M s i j) A)

/-- The pair core `log D_ij − ½log P_i − ½log P_j`. -/
noncomputable def coreIJ (i j : ℕ) (A : ℝ) : ℝ := logD M b p s i j A - gK M b p i A / 2 - gK M b p j A / 2

theorem cos_thIJ (i j : ℕ) : Real.cos (thIJ M i j) = CC M ((j : ℝ) - i) := rfl
theorem sin_thIJ (i j : ℕ) : Real.sin (thIJ M i j) = SS M ((j : ℝ) - i) := rfl

end

/-! ### The Hamiltonian along the ansatz as a pair sum -/

theorem ansatz_bound (M : ℕ) (b p A : ℝ) (k : ℕ) (hA : |A| + A ^ 2 * |p| < w0 b) :
    |A * cc M k + A ^ 2 * (p * CC M k)| < w0 b := by
  have e1 : |A * cc M k| ≤ |A| := by
    rw [abs_mul]; exact mul_le_of_le_one_right (abs_nonneg _) (cc_abs_le M _)
  have e2 : |A ^ 2 * (p * CC M k)| ≤ A ^ 2 * |p| := by
    rw [abs_mul, abs_mul, abs_of_nonneg (sq_nonneg A), ← mul_assoc]
    exact mul_le_of_le_one_right (by positivity) (Real.abs_cos_le_one _)
  calc |A * cc M k + A ^ 2 * (p * CC M k)| ≤ |A * cc M k| + |A ^ 2 * (p * CC M k)| := abs_add_le _ _
    _ ≤ |A| + A ^ 2 * |p| := add_le_add e1 e2
    _ < w0 b := hA

theorem pairTerm_oddAnsatz (M : ℕ) (b p s A : ℝ) (hb : -1 < b) (hb1 : b < 1)
    (hA : |A| + A ^ 2 * |p| < w0 b) (i j : Fin (2 * M + 1)) :
    pairTerm b (oddAnsatz M b p s A i) (oddAnsatz M b p s A j) = logD M b p s i j A := by
  obtain ⟨hwi, h1i⟩ := w_pos' b _ hb hb1 (ansatz_bound M b p A i hA)
  obtain ⟨hwj, h1j⟩ := w_pos' b _ hb hb1 (ansatz_bound M b p A j hA)
  rw [← add_assoc] at hwi h1i hwj h1j
  have ei : w0 b + A * cc M ((i : ℕ) : ℝ) + A ^ 2 * p * CC M ((i : ℕ) : ℝ)
      = w0 b + A * cc M ((i : ℕ) : ℝ) + A ^ 2 * (p * CC M ((i : ℕ) : ℝ)) := by ring
  have ej : w0 b + A * cc M ((j : ℕ) : ℝ) + A ^ 2 * p * CC M ((j : ℕ) : ℝ)
      = w0 b + A * cc M ((j : ℕ) : ℝ) + A ^ 2 * (p * CC M ((j : ℕ) : ℝ)) := by ring
  unfold oddAnsatz
  rw [ei, ej, pairTerm_canonical b _ _ _ _ hwi hwj h1i h1j]
  unfold logD Dg4
  congr 2
  rw [← Real.cos_neg]
  congr 1
  unfold thIJ dIJ
  ring

theorem Ham_oddAnsatz_eq (M : ℕ) (b p s A : ℝ) (hb : -1 < b) (hb1 : b < 1)
    (hA : |A| + A ^ 2 * |p| < w0 b) :
    Ham b (oddAnsatz M b p s A) = -(1 / 2) * ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else logD M b p s i j A) := by
  unfold Ham
  congr 1
  rw [← sum_fin_offdiag]
  refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => ?_))
  split_ifs with h
  · rfl
  · exact pairTerm_oddAnsatz M b p s A hb hb1 hA i j

/-- Splitting off the one-body part: `Σ_{i≠j} log D_ij = Σ_{i≠j} core_ij + (n−1) Σ_k log P_k`. -/
theorem sum_logD_split (M : ℕ) (b p s A : ℝ) :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1), (if i = j then 0 else logD M b p s i j A)
      = ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1), (if i = j then 0 else coreIJ M b p s i j A)
        + (((2 * M + 1 : ℕ) : ℝ) - 1) * ∑ k ∈ range (2 * M + 1), gK M b p k A := by
  have h := offdiag_add (2 * M + 1) (fun k => gK M b p k A / 2)
  have e : ∀ i j : ℕ, (if i = j then (0 : ℝ) else logD M b p s i j A)
      = (if i = j then 0 else coreIJ M b p s i j A)
        + (if i = j then 0 else gK M b p i A / 2 + gK M b p j A / 2) := by
    intro i j; split_ifs <;> simp [coreIJ]
  simp_rw [e, Finset.sum_add_distrib]
  rw [h, ← Finset.sum_div]
  ring

/-- `A = 0` is the unit ring. -/
theorem oddAnsatz_zero (M : ℕ) (b p s : ℝ) (hb : -1 < b) :
    oddAnsatz M b p s 0 = ringConf (2 * M + 1) := by
  have h1 : 0 < 1 + b := by linarith
  have hu : uOf b (w0 b) = 1 := by unfold uOf w0; field_simp; ring
  funext k
  unfold oddAnsatz ringConf
  simp only [zero_mul, add_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, hu,
    Real.sqrt_one, Complex.ofReal_one, one_mul]
  congr 1
  push_cast
  ring

/-! ### Coefficient sums -/

section Coeffs
variable (M : ℕ) (hM : 2 ≤ M) (b p s : ℝ)
include hM

local notation "nR" => ((2 * M + 1 : ℕ) : ℝ)

theorem coef2_pairs :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1), (if i = j then 0 else
      (Real.sin (thIJ M i j) * dIJ M s i j + Y2c b (cc M i) (cc M j)) / (1 - Real.cos (thIJ M i j)))
      = (1 + b) ^ 4 / 8 * (nR * (nR ^ 2 - 1) / 4) := by
  have e : ∀ i j : ℕ, (if i = j then (0 : ℝ) else
      (Real.sin (thIJ M i j) * dIJ M s i j + Y2c b (cc M i) (cc M j)) / (1 - Real.cos (thIJ M i j)))
      = s * (if i = j then 0 else (SS M j - SS M i) * SS M ((j : ℝ) - i) / (1 - CC M ((j : ℝ) - i)) ^ 1)
        + (1 + b) ^ 4 / 8 * (if i = j then 0 else (cc M i - cc M j) ^ 2 / (1 - CC M ((j : ℝ) - i)) ^ 1) := by
    intro i j
    split_ifs
    · simp
    · rw [sin_thIJ, cos_thIJ]; unfold dIJ Y2c; ring
  simp_rw [e, Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [pairLin M hM, pair2 M hM]
  ring

theorem coef3_pairs :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1), (if i = j then 0 else
      Y3c b (cc M i) (p * CC M i) (cc M j) (p * CC M j) / (1 - Real.cos (thIJ M i j))) = 0 := by
  have e : ∀ i j : ℕ, (if i = j then (0 : ℝ) else
      Y3c b (cc M i) (p * CC M i) (cc M j) (p * CC M j) / (1 - Real.cos (thIJ M i j)))
      = (b - 1) * (1 + b) ^ 5 / 8 * (if i = j then 0 else
          (cc M i - cc M j) ^ 2 * (cc M i + cc M j) / (1 - CC M ((j : ℝ) - i)) ^ 1)
        + (1 + b) ^ 4 / 4 * p * (if i = j then 0 else
          (cc M i - cc M j) * (CC M i - CC M j) / (1 - CC M ((j : ℝ) - i)) ^ 1) := by
    intro i j
    split_ifs
    · simp
    · rw [cos_thIJ]; unfold Y3c; ring
  simp_rw [e, Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [pair3a M hM, pair3b M hM]
  ring

/-- Pointwise split of the `A⁴` pair coefficient. -/
theorem c4_split (ci cj Ci Cj Si Sj co si : ℝ) (hK : 1 - co ≠ 0) (hpy : si ^ 2 + co ^ 2 = 1) :
    (co * (s * (Sj - Si)) ^ 2 / 2 + Y4c b ci (p * Ci) cj (p * Cj)) / (1 - co)
      - (si * (s * (Sj - Si)) + Y2c b ci cj) ^ 2 / (2 * (1 - co) ^ 2)
      = 3 * (1 + b) ^ 6 * (5 * b ^ 2 - 6 * b + 5) / 128 * ((ci - cj) ^ 2 * (ci ^ 2 + cj ^ 2) / (1 - co) ^ 1)
        + (1 + b) ^ 6 * (9 * b ^ 2 - 14 * b + 9) / 64 * ((ci - cj) ^ 2 * (ci * cj) / (1 - co) ^ 1)
        + (b - 1) * (1 + b) ^ 5 / 8 * p * ((ci - cj) ^ 2 * (Ci + Cj) / (1 - co) ^ 1)
        + (b - 1) * (1 + b) ^ 5 / 4 * p * ((ci - cj) * (Ci - Cj) * (ci + cj) / (1 - co) ^ 1)
        + (1 + b) ^ 4 / 8 * p ^ 2 * ((Ci - Cj) ^ 2 / (1 - co) ^ 1)
        - s ^ 2 / 2 * ((Sj - Si) ^ 2 / (1 - co) ^ 1)
        - s * (1 + b) ^ 4 / 8 * ((ci - cj) ^ 2 * (Sj - Si) * si / (1 - co) ^ 2)
        - (1 + b) ^ 8 / 128 * ((ci - cj) ^ 4 / (1 - co) ^ 2) := by
  have h1 : (si * (s * (Sj - Si)) + Y2c b ci cj) ^ 2
      = (1 - co ^ 2) * (s * (Sj - Si)) ^ 2 + 2 * si * (s * (Sj - Si)) * Y2c b ci cj
        + Y2c b ci cj ^ 2 := by
    linear_combination (s * (Sj - Si)) ^ 2 * hpy
  rw [h1]
  unfold Y4c Y2c
  field_simp
  ring

theorem coef4_pairs :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1), (if i = j then 0 else
      ((Real.cos (thIJ M i j) * dIJ M s i j ^ 2 / 2 + Y4c b (cc M i) (p * CC M i) (cc M j) (p * CC M j))
          / (1 - Real.cos (thIJ M i j))
        - (Real.sin (thIJ M i j) * dIJ M s i j + Y2c b (cc M i) (cc M j)) ^ 2
          / (2 * (1 - Real.cos (thIJ M i j)) ^ 2)))
      = 3 * (1 + b) ^ 6 * (5 * b ^ 2 - 6 * b + 5) / 128 * (nR * (nR - 1) * (3 * nR + 1) / 8)
        + (1 + b) ^ 6 * (9 * b ^ 2 - 14 * b + 9) / 64 * (-(nR * (nR - 1) * (3 * nR - 1) / 16))
        + (b - 1) * (1 + b) ^ 5 / 8 * p * (nR * (nR - 1) ^ 2 / 4)
        + (b - 1) * (1 + b) ^ 5 / 4 * p * (nR * (nR - 1) / 2)
        + (1 + b) ^ 4 / 8 * p ^ 2 * (nR * (nR - 1))
        - s ^ 2 / 2 * (nR * (nR - 1))
        - s * (1 + b) ^ 4 / 8 * (nR * (nR - 1) ^ 2 / 4)
        - (1 + b) ^ 8 / 128 * (nR * (nR - 1) * (nR + 1) * (nR ^ 2 + 3) / 32) := by
  have e : ∀ i ∈ range (2 * M + 1), ∀ j ∈ range (2 * M + 1), (if i = j then (0 : ℝ) else
      ((Real.cos (thIJ M i j) * dIJ M s i j ^ 2 / 2 + Y4c b (cc M i) (p * CC M i) (cc M j) (p * CC M j))
          / (1 - Real.cos (thIJ M i j))
        - (Real.sin (thIJ M i j) * dIJ M s i j + Y2c b (cc M i) (cc M j)) ^ 2
          / (2 * (1 - Real.cos (thIJ M i j)) ^ 2)))
      = 3 * (1 + b) ^ 6 * (5 * b ^ 2 - 6 * b + 5) / 128 * (if i = j then 0 else
          (cc M i - cc M j) ^ 2 * (cc M i ^ 2 + cc M j ^ 2) / (1 - CC M ((j : ℝ) - i)) ^ 1)
        + (1 + b) ^ 6 * (9 * b ^ 2 - 14 * b + 9) / 64 * (if i = j then 0 else
          (cc M i - cc M j) ^ 2 * (cc M i * cc M j) / (1 - CC M ((j : ℝ) - i)) ^ 1)
        + (b - 1) * (1 + b) ^ 5 / 8 * p * (if i = j then 0 else
          (cc M i - cc M j) ^ 2 * (CC M i + CC M j) / (1 - CC M ((j : ℝ) - i)) ^ 1)
        + (b - 1) * (1 + b) ^ 5 / 4 * p * (if i = j then 0 else
          (cc M i - cc M j) * (CC M i - CC M j) * (cc M i + cc M j) / (1 - CC M ((j : ℝ) - i)) ^ 1)
        + (1 + b) ^ 4 / 8 * p ^ 2 * (if i = j then 0 else
          (CC M i - CC M j) ^ 2 / (1 - CC M ((j : ℝ) - i)) ^ 1)
        - s ^ 2 / 2 * (if i = j then 0 else (SS M j - SS M i) ^ 2 / (1 - CC M ((j : ℝ) - i)) ^ 1)
        - s * (1 + b) ^ 4 / 8 * (if i = j then 0 else
          (cc M i - cc M j) ^ 2 * (SS M j - SS M i) * SS M ((j : ℝ) - i) / (1 - CC M ((j : ℝ) - i)) ^ 2)
        - (1 + b) ^ 8 / 128 * (if i = j then 0 else
          (cc M i - cc M j) ^ 4 / (1 - CC M ((j : ℝ) - i)) ^ 2) := by
    intro i hi j hj
    split_ifs with hij
    · simp
    · have hlt := CC_lt_one M i j (Finset.mem_range.mp hi) (Finset.mem_range.mp hj) hij
      rw [sin_thIJ, cos_thIJ]
      unfold dIJ
      exact c4_split M hM b p s _ _ _ _ _ _ _ _ (by linarith)
        (by rw [← sin_thIJ, ← cos_thIJ]; exact Real.sin_sq_add_cos_sq _)
  rw [Finset.sum_congr rfl (fun i hi => Finset.sum_congr rfl (fun j hj => e i hi j hj))]
  simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]
  rw [pairB M hM, pairC M hM, pairA M hM, pairD M hM, pairE M hM, pairS M hM, pairG M hM, pairF M hM]

end Coeffs

/-! ### The main theorem -/

/-- **Odd `n = 2M+1 ≥ 5`: the fourth-order jet of `H` along the paper's ansatz.**
`H(oddAnsatz A) − H(ring) = 0·A + f₂A² + 0·A³ + E(p,s)A⁴ + O(A⁵)`. -/
theorem Ham_oddAnsatz_jet4 (M : ℕ) (hM : 2 ≤ M) (b p s : ℝ) (hb : -1 < b) (hb1 : b < 1) :
    Jet4 (fun A => Ham b (oddAnsatz M b p s A) - Ham b (ringConf (2 * M + 1))) 0 0
      (f2 ((2 * M + 1 : ℕ) : ℝ) b) 0 (Efun ((2 * M + 1 : ℕ) : ℝ) b p s) := by
  set n := 2 * M + 1 with hn
  -- pair jets
  have hpair : ∀ i ∈ range n, ∀ j ∈ range n,
      Jet4 (fun A => if i = j then 0 else coreIJ M b p s i j A)
        (if i = j then 0 else Real.log 2 + Real.log (1 - Real.cos (thIJ M i j)))
        (if i = j then 0 else 0)
        (if i = j then 0 else (Real.sin (thIJ M i j) * dIJ M s i j + Y2c b (cc M i) (cc M j))
          / (1 - Real.cos (thIJ M i j)))
        (if i = j then 0 else Y3c b (cc M i) (p * CC M i) (cc M j) (p * CC M j)
          / (1 - Real.cos (thIJ M i j)))
        (if i = j then 0 else
          (Real.cos (thIJ M i j) * dIJ M s i j ^ 2 / 2 + Y4c b (cc M i) (p * CC M i) (cc M j) (p * CC M j))
            / (1 - Real.cos (thIJ M i j))
          - (Real.sin (thIJ M i j) * dIJ M s i j + Y2c b (cc M i) (cc M j)) ^ 2
            / (2 * (1 - Real.cos (thIJ M i j)) ^ 2)) := by
    intro i hi j hj
    refine Jet4.ite (i = j) (fun hij => ?_)
    have hlt := CC_lt_one M i j (Finset.mem_range.mp hi) (Finset.mem_range.mp hj) hij
    exact jet_pair b _ _ _ _ _ _ hb hlt
  have hrow : ∀ i ∈ range n, Jet4 (fun A => ∑ j ∈ range n, if i = j then 0 else coreIJ M b p s i j A)
      (∑ j ∈ range n, if i = j then 0 else Real.log 2 + Real.log (1 - Real.cos (thIJ M i j)))
      (∑ j ∈ range n, if i = j then (0 : ℝ) else 0)
      (∑ j ∈ range n, if i = j then 0 else (Real.sin (thIJ M i j) * dIJ M s i j + Y2c b (cc M i) (cc M j))
          / (1 - Real.cos (thIJ M i j)))
      (∑ j ∈ range n, if i = j then 0 else Y3c b (cc M i) (p * CC M i) (cc M j) (p * CC M j)
          / (1 - Real.cos (thIJ M i j)))
      (∑ j ∈ range n, if i = j then 0 else
          (Real.cos (thIJ M i j) * dIJ M s i j ^ 2 / 2 + Y4c b (cc M i) (p * CC M i) (cc M j) (p * CC M j))
            / (1 - Real.cos (thIJ M i j))
          - (Real.sin (thIJ M i j) * dIJ M s i j + Y2c b (cc M i) (cc M j)) ^ 2
            / (2 * (1 - Real.cos (thIJ M i j)) ^ 2)) :=
    fun i hi => Jet4.sum (range n) (fun j hj => hpair i hi j hj)
  have hall := Jet4.sum (range n) hrow
  -- one-body jets
  have hone := Jet4.sum (range n) (fun k (_ : k ∈ range n) => jet_logPw b (cc M k) (p * CC M k) hb)
  have hH := (hall.const_mul (-(1 / 2))).add (hone.const_mul (-(1 / 2) * ((n : ℝ) - 1)))
  -- identify with the Hamiltonian near `A = 0`
  have hw0 : 0 < w0 b := w0_pos b hb
  have hev : ∀ᶠ A in 𝓝 (0 : ℝ), |A| + A ^ 2 * |p| < w0 b := by
    have hc : Continuous (fun A : ℝ => |A| + A ^ 2 * |p|) := by fun_prop
    exact hc.continuousAt.eventually_lt continuousAt_const (by simpa using hw0)
  have hEq : ∀ᶠ A in 𝓝 (0 : ℝ), -(1 / 2) * (∑ i ∈ range n, ∑ j ∈ range n,
        if i = j then 0 else coreIJ M b p s i j A)
      + -(1 / 2) * ((n : ℝ) - 1) * (∑ k ∈ range n,
        Real.log (Pw b (w0 b + A * cc M k + A ^ 2 * (p * CC M k))))
      = Ham b (oddAnsatz M b p s A) := by
    filter_upwards [hev] with A hA
    rw [Ham_oddAnsatz_eq M b p s A hb hb1 hA, sum_logD_split]
    unfold gK
    ring
  have hHam := hH.congr' hEq
  have h0 : Ham b (ringConf n) = Ham b (oddAnsatz M b p s 0) := by rw [oddAnsatz_zero M b p s hb]
  have hval := hHam.eval_zero
  refine (hHam.sub (Jet4.const (Ham b (ringConf n)))).of_coeff ?_ ?_ ?_ ?_ ?_
  · rw [h0, hval]; ring
  · -- `A¹`
    simp only [Finset.sum_const_zero, ite_self, mul_zero, zero_add, sub_zero]
    rw [← Finset.mul_sum, sum_cc M (by omega)]
    ring
  · -- `A²`
    rw [coef2_pairs M hM b s]
    simp only [Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.sum_mul]
    rw [sum_CC M (by omega), sum_cc_sq M (by omega)]
    unfold f2 L1 L2
    simp only [hn]
    push_cast
    ring
  · -- `A³`
    rw [coef3_pairs M hM b p]
    have e : ∀ k ∈ range n, 2 * L2 b * cc M k * (p * CC M k) + L3 b * cc M k ^ 3
        = (2 * L2 b * p) * (cc M k * CC M k) + L3 b * cc M k ^ 3 := by intro k _; ring
    rw [Finset.sum_congr rfl e, Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum,
      sum_cc_CC M hM, sum_cc_cube M hM]
    ring
  · -- `A⁴`
    rw [coef4_pairs M hM b p s]
    have e : ∀ k ∈ range n, L2 b * (p * CC M k) ^ 2 + 3 * L3 b * cc M k ^ 2 * (p * CC M k)
        + L4 b * cc M k ^ 4
        = (L2 b * p ^ 2) * CC M k ^ 2 + (3 * L3 b * p) * (cc M k ^ 2 * CC M k) + L4 b * cc M k ^ 4 := by
      intro k _; ring
    rw [Finset.sum_congr rfl e, Finset.sum_add_distrib, Finset.sum_add_distrib, ← Finset.mul_sum,
      ← Finset.mul_sum, ← Finset.mul_sum, sum_CC_sq M (by omega), sum_cc_sq_CC M (by omega),
      sum_cc_four M hM]
    unfold Efun E00 E10 E01 E20 E02 L2 L3 L4
    simp only [hn]
    push_cast
    ring

/-- **Odd `n = 2M+1 ≥ 5` (paper, proof of Thm odd, "Result").** Along the ansatz
`w_k = w₀ + A cos(lφ_k) + A² p cos φ_k`, `ϑ_k = φ_k + A² s sin φ_k`, for every `λ ∈ (−1,1)` and all
`p, s`: `H − H(ring) = f₂ A² + E(p,s) A⁴ + O(A⁵)`, with `f₂ = f2 n λ` and
`E(p,s) = E₀₀ + E₁₀p + E₀₁s + E₂₀p² + E₀₂s²` the definitions of `Conj13Algebra.lean`
(no `A¹`, no `A³` term). -/
theorem Ham_oddAnsatz_jet (M : ℕ) (hM : 2 ≤ M) (b p s : ℝ) (hb : -1 < b) (hb1 : b < 1) :
    (fun A => Ham b (oddAnsatz M b p s A) - Ham b (ringConf (2 * M + 1))
        - f2 ((2 * M + 1 : ℕ) : ℝ) b * A ^ 2 - Efun ((2 * M + 1 : ℕ) : ℝ) b p s * A ^ 4)
      =O[𝓝 0] (fun A => A ^ 5) := by
  have h := Ham_oddAnsatz_jet4 M hM b p s hb hb1
  unfold Jet4 at h
  refine h.congr_left (fun A => ?_)
  ring

end Havelock.Conj13.Odd4
