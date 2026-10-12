/-
Havelock — Conjecture 13, odd `n = 2M+1`, `M ≥ 2`: the pair sums of the `A²`, `A³`, `A⁴`
coefficients (paper `conj13/paper/main.tex`, proof of Thm odd, "Sums over roots of unity").

Sums run over ordered pairs `i ≠ j` in `[0, n)`; `c_k = cos(2πMk/n)`, `C_k = cos(2πk/n)`,
`S_k = sin(2πk/n)`, `K_ij = 1 − cos(2π(j−i)/n)`, `σ_ij = sin(2π(j−i)/n)`.

* `pair2`  `Σ (c_i−c_j)²/K = n(n²−1)/4`
* `pairLin` `Σ (S_j−S_i)σ/K = 0`;  `pair3a` `Σ (c_i−c_j)²(c_i+c_j)/K = 0`;  `pair3b` `Σ (c_i−c_j)(C_i−C_j)/K = 0`
* `pairA`  `Σ (c_i−c_j)²(C_i+C_j)/K = n(n−1)²/4`
* `pairB`  `Σ (c_i−c_j)²(c_i²+c_j²)/K = n(n−1)(3n+1)/8`
* `pairC`  `Σ (c_i−c_j)²c_ic_j/K = −n(n−1)(3n−1)/16`
* `pairD`  `Σ (c_i−c_j)(C_i−C_j)(c_i+c_j)/K = n(n−1)/2`
* `pairE`  `Σ (C_i−C_j)²/K = n(n−1)`,   `pairS`  `Σ (S_j−S_i)²/K = n(n−1)`
* `pairF`  `Σ (c_i−c_j)⁴/K² = n(n−1)(n+1)(n²+3)/32`
* `pairG`  `Σ (c_i−c_j)²(S_j−S_i)σ/K² = n(n−1)²/4`

Each: reindex by the offset `m = j − i`, evaluate the inner sum (`Odd4Inner`), reduce it to
`Σ_q a_q cos(2πqm/n)` using `n·(πm/n) = πm`, and apply `fejer1` / `fejer2` (`Odd4Fejer`).
-/
import Havelock.Conj13.Odd4Inner

open Real Finset

namespace Havelock.Conj13.Odd4

open Havelock.Conj13.OddSums Havelock.Conj13.Odd

theorem SS_per (M : ℕ) (x : ℝ) : SS M (x - (2 * M + 1 : ℕ)) = SS M x := by
  unfold SS
  have hn : ((2 * M + 1 : ℕ) : ℝ) ≠ 0 := by positivity
  rw [show 2 * π * (x - (2 * M + 1 : ℕ)) / (2 * M + 1 : ℕ) = 2 * π * x / (2 * M + 1 : ℕ)
      - ((1 : ℤ) : ℝ) * (2 * π) by push_cast; field_simp]
  exact Real.sin_sub_int_mul_two_pi _ _

/-- Reindexing a pair sum by offsets. -/
theorem pair_offsets (M : ℕ) (F : ℝ → ℝ → ℝ)
    (hF : ∀ a b : ℝ, F a (b - (2 * M + 1 : ℕ)) = F a b) :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1), (if i = j then 0 else F i j)
      = ∑ m ∈ Ico 1 (2 * M + 1), ∑ i ∈ range (2 * M + 1), F i (i + m) := by
  have h := reindex_pairs (2 * M + 1) (fun (i : ℕ) (x : ℝ) => F i (i + x))
    (fun i x => by simp only [show (i : ℝ) + (x - (2 * M + 1 : ℕ)) = (i + x) - (2 * M + 1 : ℕ) by ring, hF])
  simp only [add_sub_cancel] at h
  rw [h, Finset.sum_comm]

/-- The Fejér atoms `cos(2πqm/n)` and weight `1 − cos(2πm/n)`. -/
noncomputable def fa (M q m : ℕ) : ℝ := Real.cos (2 * π * q * m / (2 * M + 1 : ℕ))
noncomputable def kw (M m : ℕ) : ℝ := 1 - Real.cos (2 * π * m / (2 * M + 1 : ℕ))

theorem F1 (M q : ℕ) (hq : q ≤ 2 * M + 1) :
    ∑ m ∈ Ico 1 (2 * M + 1), fa M q m / kw M m
      = ((((2 * M + 1 : ℕ) : ℝ)) ^ 2 - 1) / 6 - (q : ℝ) * (((2 * M + 1 : ℕ) : ℝ) - q) :=
  fejer1 (2 * M + 1) q (by omega) hq

theorem F2' (M q : ℕ) (hq : q ≤ 2 * M + 1) :
    ∑ m ∈ Ico 1 (2 * M + 1), fa M q m / kw M m ^ 2 = fejer2Poly ((2 * M + 1 : ℕ) : ℝ) q :=
  fejer2 (2 * M + 1) q (by omega) hq

theorem kw_pos (M m : ℕ) (hm : m ∈ Ico 1 (2 * M + 1)) : 0 < kw M m := by
  have := Finset.mem_Ico.mp hm
  exact one_sub_cos_pos (2 * M + 1) m (by omega) this.2

/-! ### Relations at an integer offset `m` (`Y = πm/n`, `X = MY`, `nY = πm`) -/

section Rel
variable (M m : ℕ)

theorem nY : ((2 * M + 1 : ℕ) : ℝ) * hh M m = π * m := by
  unfold hh; have : ((2 * M + 1 : ℕ) : ℝ) ≠ 0 := by positivity
  field_simp

theorem two_nY : 2 * (2 * ((M : ℝ) * hh M m)) + 2 * hh M m = (m : ℝ) * (2 * π) := by
  have h := nY M m; push_cast at h; linear_combination 2 * h

theorem rel_s2 : Real.sin ((M : ℝ) * hh M m) ^ 2 = (1 - Real.cos (2 * ((M : ℝ) * hh M m))) / 2 := by
  rw [Real.cos_two_mul]; linear_combination Real.sin_sq_add_cos_sq ((M : ℝ) * hh M m)

theorem rel_cm2 : Real.cos ((M : ℝ) * hh M m) ^ 2 = (1 + Real.cos (2 * ((M : ℝ) * hh M m))) / 2 := by
  rw [Real.cos_two_mul]; ring

theorem rel_t2 : Real.sin (hh M m) ^ 2 = (1 - Real.cos (2 * hh M m)) / 2 := by
  rw [Real.cos_two_mul]; linear_combination Real.sin_sq_add_cos_sq (hh M m)

theorem rel_c4 : Real.cos (2 * ((M : ℝ) * hh M m)) ^ 2 = (1 + Real.cos (2 * hh M m)) / 2 := by
  have h : Real.cos (2 * (2 * ((M : ℝ) * hh M m))) = Real.cos (2 * hh M m) := by
    rw [show 2 * (2 * ((M : ℝ) * hh M m)) = (m : ℝ) * (2 * π) - 2 * hh M m by
      linear_combination two_nY M m, Real.cos_nat_mul_two_pi_sub]
  rw [Real.cos_two_mul] at h
  linear_combination h / 2

theorem rel_r1 : Real.cos (hh M m) * Real.cos (((2 * M + 1 : ℕ) : ℝ) * hh M m)
    = Real.cos (2 * ((M : ℝ) * hh M m)) := by
  have e1 : ((2 * M + 1 : ℕ) : ℝ) * hh M m + hh M m = (m : ℝ) * (2 * π) - 2 * ((M : ℝ) * hh M m) := by
    have := two_nY M m; push_cast; linear_combination this
  have e2 : ((2 * M + 1 : ℕ) : ℝ) * hh M m - hh M m = 2 * ((M : ℝ) * hh M m) := by push_cast; ring
  have h1 := Real.cos_add (((2 * M + 1 : ℕ) : ℝ) * hh M m) (hh M m)
  have h2 := Real.cos_sub (((2 * M + 1 : ℕ) : ℝ) * hh M m) (hh M m)
  rw [e1, Real.cos_nat_mul_two_pi_sub] at h1
  rw [e2] at h2
  linear_combination (-1 / 2 : ℝ) * h1 - (1 / 2 : ℝ) * h2

theorem rel_r2 : Real.sin (hh M m) * Real.cos (((2 * M + 1 : ℕ) : ℝ) * hh M m)
    = -(2 * Real.sin ((M : ℝ) * hh M m) * Real.cos ((M : ℝ) * hh M m)) := by
  have e1 : ((2 * M + 1 : ℕ) : ℝ) * hh M m + hh M m = (m : ℝ) * (2 * π) - 2 * ((M : ℝ) * hh M m) := by
    have := two_nY M m; push_cast; linear_combination this
  have e2 : ((2 * M + 1 : ℕ) : ℝ) * hh M m - hh M m = 2 * ((M : ℝ) * hh M m) := by push_cast; ring
  have h1 := Real.sin_add (((2 * M + 1 : ℕ) : ℝ) * hh M m) (hh M m)
  have h2 := Real.sin_sub (((2 * M + 1 : ℕ) : ℝ) * hh M m) (hh M m)
  rw [e1, Real.sin_nat_mul_two_pi_sub, Real.sin_two_mul] at h1
  rw [e2, Real.sin_two_mul] at h2
  linear_combination (-1 / 2 : ℝ) * h1 + (1 / 2 : ℝ) * h2

theorem rel_r6 : Real.sin ((M : ℝ) * hh M m) * Real.cos ((M : ℝ) * hh M m) * Real.sin (2 * hh M m)
    = (Real.cos (2 * ((M : ℝ) * hh M m) - 2 * hh M m) - Real.cos (2 * ((M : ℝ) * hh M m))) / 4 := by
  have e : 2 * ((M : ℝ) * hh M m) + 2 * hh M m = (m : ℝ) * (2 * π) - 2 * ((M : ℝ) * hh M m) := by
    linear_combination two_nY M m
  have h1 := Real.cos_add (2 * ((M : ℝ) * hh M m)) (2 * hh M m)
  have h2 := Real.cos_sub (2 * ((M : ℝ) * hh M m)) (2 * hh M m)
  rw [e, Real.cos_nat_mul_two_pi_sub] at h1
  rw [Real.sin_two_mul ((M : ℝ) * hh M m)] at h1 h2
  linear_combination (1 / 4 : ℝ) * h1 - (1 / 4 : ℝ) * h2

theorem rel_r7 : Real.cos (2 * ((M : ℝ) * hh M m)) * Real.cos (2 * ((M : ℝ) * hh M m) - 2 * hh M m)
    = (Real.cos (2 * (2 * hh M m)) + Real.cos (2 * hh M m)) / 2 := by
  have e : 2 * ((M : ℝ) * hh M m) + (2 * ((M : ℝ) * hh M m) - 2 * hh M m)
      = (m : ℝ) * (2 * π) - 2 * (2 * hh M m) := by linear_combination two_nY M m
  have h1 := Real.cos_add (2 * ((M : ℝ) * hh M m)) (2 * ((M : ℝ) * hh M m) - 2 * hh M m)
  have h2 := Real.cos_sub (2 * ((M : ℝ) * hh M m)) (2 * ((M : ℝ) * hh M m) - 2 * hh M m)
  rw [e, Real.cos_nat_mul_two_pi_sub] at h1
  rw [show 2 * ((M : ℝ) * hh M m) - (2 * ((M : ℝ) * hh M m) - 2 * hh M m) = 2 * hh M m by ring] at h2
  linear_combination (-1 / 2 : ℝ) * h1 + (-1 / 2 : ℝ) * h2

theorem fa_M : fa M M m = Real.cos (2 * ((M : ℝ) * hh M m)) := by unfold fa hh; congr 1; ring
theorem fa_1 : fa M 1 m = Real.cos (2 * hh M m) := by unfold fa hh; congr 1; push_cast; ring
theorem fa_2 : fa M 2 m = Real.cos (2 * (2 * hh M m)) := by unfold fa hh; congr 1; push_cast; ring
theorem fa_0 : fa M 0 m = 1 := by unfold fa; simp
theorem fa_Mm1 (hM : 1 ≤ M) : fa M (M - 1) m = Real.cos (2 * ((M : ℝ) * hh M m) - 2 * hh M m) := by
  unfold fa hh; congr 1; rw [Nat.cast_sub hM]; push_cast; ring
theorem CC_m : CC M m = Real.cos (2 * hh M m) := by unfold CC hh; congr 1; ring
theorem SS_m : SS M m = Real.sin (2 * hh M m) := by unfold SS hh; congr 1; ring
theorem kw_eq : kw M m = 1 - CC M m := rfl

end Rel

/-! ### Pair sums -/

section Pairs
variable (M : ℕ) (hM : 2 ≤ M)
include hM

local notation "nR" => ((2 * M + 1 : ℕ) : ℝ)

omit hM in
theorem per1 (a b : ℝ) : CC M ((b - (2 * M + 1 : ℕ)) - a) = CC M (b - a) := by
  rw [show (b - (2 * M + 1 : ℕ)) - a = (b - a) - (2 * M + 1 : ℕ) by ring, CC_per]

omit hM in
theorem perS (a b : ℝ) : SS M ((b - (2 * M + 1 : ℕ)) - a) = SS M (b - a) := by
  rw [show (b - (2 * M + 1 : ℕ)) - a = (b - a) - (2 * M + 1 : ℕ) by ring, SS_per]

/-- Generic: a pair sum with weight `1/K` or `1/K²` reduces to the offset sum of the inner sums. -/
theorem pair_weight (N : ℝ → ℝ → ℝ) (hN : ∀ a b : ℝ, N a (b - (2 * M + 1 : ℕ)) = N a b) (r : ℕ) :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else N i j / (1 - CC M ((j : ℝ) - i)) ^ r)
      = ∑ m ∈ Ico 1 (2 * M + 1), (∑ i ∈ range (2 * M + 1), N i (i + m)) / kw M m ^ r := by
  rw [pair_offsets M (fun a b => N a b / (1 - CC M (b - a)) ^ r)
    (fun a b => by simp only [hN, per1 M a b])]
  refine Finset.sum_congr rfl (fun m _ => ?_)
  rw [Finset.sum_div]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [add_sub_cancel_left, kw_eq]

theorem pair2 :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (cc M i - cc M j) ^ 2 / (1 - CC M ((j : ℝ) - i)) ^ 1)
      = nR * (nR ^ 2 - 1) / 4 := by
  rw [pair_weight M hM (fun a b => (cc M a - cc M b) ^ 2) (fun a b => by simp only [cc_per]) 1]
  have e : ∀ m ∈ Ico 1 (2 * M + 1), (∑ i ∈ range (2 * M + 1), (cc M i - cc M (i + m)) ^ 2) / kw M m ^ 1
      = nR * (fa M 0 m / kw M m) - nR * (fa M M m / kw M m) := by
    intro m _
    have hs : ∑ i ∈ range (2 * M + 1), (cc M i - cc M (i + m)) ^ 2 = 2 * nR * Real.sin ((M : ℝ) * hh M m) ^ 2 := by
      have e2 : ∀ i ∈ range (2 * M + 1), (cc M i - cc M (i + m)) ^ 2
          = (4 * Real.sin ((M : ℝ) * hh M m) ^ 2) * Real.sin ((M : ℝ) * eta M i (hh M m)) ^ 2 := by
        intro i _; rw [cc_diff]; ring
      rw [Finset.sum_congr rfl e2, ← Finset.mul_sum, B1 M hM]; ring
    rw [hs, fa_0, fa_M, rel_s2]; ring
  rw [Finset.sum_congr rfl e, Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.mul_sum,
    F1 M 0 (by omega), F1 M M (by omega)]
  push_cast; ring


theorem pairA :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (cc M i - cc M j) ^ 2 * (CC M i + CC M j) / (1 - CC M ((j : ℝ) - i)) ^ 1)
      = nR * (nR - 1) ^ 2 / 4 := by
  rw [pair_weight M hM (fun a b => (cc M a - cc M b) ^ 2 * (CC M a + CC M b)) (fun a b => by simp only [cc_per, CC_per]) 1]
  have e : ∀ m ∈ Ico 1 (2 * M + 1), (∑ i ∈ range (2 * M + 1),
      (fun a b : ℝ => (cc M a - cc M b) ^ 2 * (CC M a + CC M b)) i (i + m)) / kw M m ^ 1 = (-nR) * (fa M M m / kw M m) + (nR / 2) * (fa M 0 m / kw M m) + (nR / 2) * (fa M 1 m / kw M m) := by
    intro m _
    have hc : ∑ i ∈ range (2 * M + 1), (fun a b : ℝ => (cc M a - cc M b) ^ 2 * (CC M a + CC M b)) i (i + m) = (-nR) * fa M M m + (nR / 2) * fa M 0 m + (nR / 2) * fa M 1 m := by
      rw [inner_a M hM m]
      have s2 := rel_s2 M m
      have cm2 := rel_cm2 M m
      have t2 := rel_t2 M m
      have c4 := rel_c4 M m
      have r1 := rel_r1 M m
      have r2 := rel_r2 M m
      rw [← fa_M, ← fa_1] at *
      rw [fa_0]
      linear_combination (-2*(nR)*(Real.sin ((M : ℝ) * hh M m))^2) * r1 + (-2*(nR)*(fa M M m)) * s2 + (nR) * c4
    rw [hc]; ring
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [F1 M M (by omega), F1 M 0 (by omega), F1 M 1 (by omega)]
  push_cast; ring

theorem pairB :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (cc M i - cc M j) ^ 2 * (cc M i ^ 2 + cc M j ^ 2) / (1 - CC M ((j : ℝ) - i)) ^ 1)
      = nR * (nR - 1) * (3 * nR + 1) / 8 := by
  rw [pair_weight M hM (fun a b => (cc M a - cc M b) ^ 2 * (cc M a ^ 2 + cc M b ^ 2)) (fun a b => by simp only [cc_per]) 1]
  have e : ∀ m ∈ Ico 1 (2 * M + 1), (∑ i ∈ range (2 * M + 1),
      (fun a b : ℝ => (cc M a - cc M b) ^ 2 * (cc M a ^ 2 + cc M b ^ 2)) i (i + m)) / kw M m ^ 1 = (5 * nR / 4) * (fa M 0 m / kw M m) + (-(3 * nR / 2)) * (fa M M m / kw M m) + (nR / 4) * (fa M 1 m / kw M m) := by
    intro m _
    have hc : ∑ i ∈ range (2 * M + 1), (fun a b : ℝ => (cc M a - cc M b) ^ 2 * (cc M a ^ 2 + cc M b ^ 2)) i (i + m) = (5 * nR / 4) * fa M 0 m + (-(3 * nR / 2)) * fa M M m + (nR / 4) * fa M 1 m := by
      rw [inner_b M hM m]
      have s2 := rel_s2 M m
      have cm2 := rel_cm2 M m
      have t2 := rel_t2 M m
      have c4 := rel_c4 M m
      have r1 := rel_r1 M m
      have r2 := rel_r2 M m
      rw [← fa_M, ← fa_1] at *
      rw [fa_0]
      linear_combination (2*(nR)*(Real.cos ((M : ℝ) * hh M m))^2 - (nR) + 4*(nR)*((Real.sin ((M : ℝ) * hh M m))^2 + (1 - (fa M M m))/2)) * s2 + ((nR)*(1 - (fa M M m))) * cm2 + ((nR)/2) * c4
    rw [hc]; ring
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [F1 M 0 (by omega), F1 M M (by omega), F1 M 1 (by omega)]
  push_cast; ring

theorem pairC :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (cc M i - cc M j) ^ 2 * (cc M i * cc M j) / (1 - CC M ((j : ℝ) - i)) ^ 1)
      = -(nR * (nR - 1) * (3 * nR - 1) / 16) := by
  rw [pair_weight M hM (fun a b => (cc M a - cc M b) ^ 2 * (cc M a * cc M b)) (fun a b => by simp only [cc_per]) 1]
  have e : ∀ m ∈ Ico 1 (2 * M + 1), (∑ i ∈ range (2 * M + 1),
      (fun a b : ℝ => (cc M a - cc M b) ^ 2 * (cc M a * cc M b)) i (i + m)) / kw M m ^ 1 = (-(nR / 2)) * (fa M 0 m / kw M m) + (3 * nR / 4) * (fa M M m / kw M m) + (-(nR / 4)) * (fa M 1 m / kw M m) := by
    intro m _
    have hc : ∑ i ∈ range (2 * M + 1), (fun a b : ℝ => (cc M a - cc M b) ^ 2 * (cc M a * cc M b)) i (i + m) = (-(nR / 2)) * fa M 0 m + (3 * nR / 4) * fa M M m + (-(nR / 4)) * fa M 1 m := by
      rw [inner_c M hM m]
      have s2 := rel_s2 M m
      have cm2 := rel_cm2 M m
      have t2 := rel_t2 M m
      have c4 := rel_c4 M m
      have r1 := rel_r1 M m
      have r2 := rel_r2 M m
      rw [← fa_M, ← fa_1] at *
      rw [fa_0]
      linear_combination ((nR)/2 - 2*(nR)*((Real.sin ((M : ℝ) * hh M m))^2 + (1 - (fa M M m))/2)) * s2 + (-(nR)/2) * c4
    rw [hc]; ring
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [F1 M 0 (by omega), F1 M M (by omega), F1 M 1 (by omega)]
  push_cast; ring

theorem pairD :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (cc M i - cc M j) * (CC M i - CC M j) * (cc M i + cc M j) / (1 - CC M ((j : ℝ) - i)) ^ 1)
      = nR * (nR - 1) / 2 := by
  rw [pair_weight M hM (fun a b => (cc M a - cc M b) * (CC M a - CC M b) * (cc M a + cc M b)) (fun a b => by simp only [cc_per, CC_per]) 1]
  have e : ∀ m ∈ Ico 1 (2 * M + 1), (∑ i ∈ range (2 * M + 1),
      (fun a b : ℝ => (cc M a - cc M b) * (CC M a - CC M b) * (cc M a + cc M b)) i (i + m)) / kw M m ^ 1 = (nR / 2) * (fa M 0 m / kw M m) + (-(nR / 2)) * (fa M 1 m / kw M m) := by
    intro m _
    have hc : ∑ i ∈ range (2 * M + 1), (fun a b : ℝ => (cc M a - cc M b) * (CC M a - CC M b) * (cc M a + cc M b)) i (i + m) = (nR / 2) * fa M 0 m + (-(nR / 2)) * fa M 1 m := by
      rw [inner_d M hM m]
      have s2 := rel_s2 M m
      have cm2 := rel_cm2 M m
      have t2 := rel_t2 M m
      have c4 := rel_c4 M m
      have r1 := rel_r1 M m
      have r2 := rel_r2 M m
      rw [← fa_M, ← fa_1] at *
      rw [fa_0]
      linear_combination (-2*(nR)*(Real.sin ((M : ℝ) * hh M m))*(Real.cos ((M : ℝ) * hh M m))) * r2 + (4*(nR)*(Real.cos ((M : ℝ) * hh M m))^2) * s2 + (4*(nR)*((1 - (fa M M m))/2)) * cm2 + (-(nR)) * c4
    rw [hc]; ring
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [F1 M 0 (by omega), F1 M 1 (by omega)]
  push_cast; ring

theorem pairE :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (CC M i - CC M j) ^ 2 / (1 - CC M ((j : ℝ) - i)) ^ 1)
      = nR * (nR - 1) := by
  rw [pair_weight M hM (fun a b => (CC M a - CC M b) ^ 2) (fun a b => by simp only [CC_per]) 1]
  have e : ∀ m ∈ Ico 1 (2 * M + 1), (∑ i ∈ range (2 * M + 1),
      (fun a b : ℝ => (CC M a - CC M b) ^ 2) i (i + m)) / kw M m ^ 1 = (nR) * (fa M 0 m / kw M m) + (-nR) * (fa M 1 m / kw M m) := by
    intro m _
    have hc : ∑ i ∈ range (2 * M + 1), (fun a b : ℝ => (CC M a - CC M b) ^ 2) i (i + m) = (nR) * fa M 0 m + (-nR) * fa M 1 m := by
      rw [inner_e M hM m]
      have s2 := rel_s2 M m
      have cm2 := rel_cm2 M m
      have t2 := rel_t2 M m
      have c4 := rel_c4 M m
      have r1 := rel_r1 M m
      have r2 := rel_r2 M m
      rw [← fa_M, ← fa_1] at *
      rw [fa_0]
      linear_combination (2*(nR)) * t2
    rw [hc]; ring
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [F1 M 0 (by omega), F1 M 1 (by omega)]
  push_cast; ring

theorem pairS :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (SS M j - SS M i) ^ 2 / (1 - CC M ((j : ℝ) - i)) ^ 1)
      = nR * (nR - 1) := by
  rw [pair_weight M hM (fun a b => (SS M b - SS M a) ^ 2) (fun a b => by simp only [SS_per]) 1]
  have e : ∀ m ∈ Ico 1 (2 * M + 1), (∑ i ∈ range (2 * M + 1),
      (fun a b : ℝ => (SS M b - SS M a) ^ 2) i (i + m)) / kw M m ^ 1 = (nR) * (fa M 0 m / kw M m) + (-nR) * (fa M 1 m / kw M m) := by
    intro m _
    have hc : ∑ i ∈ range (2 * M + 1), (fun a b : ℝ => (SS M b - SS M a) ^ 2) i (i + m) = (nR) * fa M 0 m + (-nR) * fa M 1 m := by
      rw [inner_s M hM m]
      have s2 := rel_s2 M m
      have cm2 := rel_cm2 M m
      have t2 := rel_t2 M m
      have c4 := rel_c4 M m
      have r1 := rel_r1 M m
      have r2 := rel_r2 M m
      rw [← fa_M, ← fa_1] at *
      rw [fa_0]
      linear_combination (2*(nR)) * t2
    rw [hc]; ring
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [F1 M 0 (by omega), F1 M 1 (by omega)]
  push_cast; ring

theorem pairF :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (cc M i - cc M j) ^ 4 / (1 - CC M ((j : ℝ) - i)) ^ 2)
      = nR * (nR - 1) * (nR + 1) * (nR ^ 2 + 3) / 32 := by
  rw [pair_weight M hM (fun a b => (cc M a - cc M b) ^ 4) (fun a b => by simp only [cc_per]) 2]
  have e : ∀ m ∈ Ico 1 (2 * M + 1), (∑ i ∈ range (2 * M + 1),
      (fun a b : ℝ => (cc M a - cc M b) ^ 4) i (i + m)) / kw M m ^ 2 = (9 * nR / 4) * (fa M 0 m / kw M m ^ 2) + (-(3 * nR)) * (fa M M m / kw M m ^ 2) + (3 * nR / 4) * (fa M 1 m / kw M m ^ 2) := by
    intro m _
    have hc : ∑ i ∈ range (2 * M + 1), (fun a b : ℝ => (cc M a - cc M b) ^ 4) i (i + m) = (9 * nR / 4) * fa M 0 m + (-(3 * nR)) * fa M M m + (3 * nR / 4) * fa M 1 m := by
      rw [inner_f M hM m]
      have s2 := rel_s2 M m
      have cm2 := rel_cm2 M m
      have t2 := rel_t2 M m
      have c4 := rel_c4 M m
      have r1 := rel_r1 M m
      have r2 := rel_r2 M m
      rw [← fa_M, ← fa_1] at *
      rw [fa_0]
      linear_combination (6*(nR)*((Real.sin ((M : ℝ) * hh M m))^2 + (1 - (fa M M m))/2)) * s2 + (3*(nR)/2) * c4
    rw [hc]; ring
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [F2' M 0 (by omega), F2' M M (by omega), F2' M 1 (by omega)]
  push_cast; unfold fejer2Poly; ring

theorem pairG :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (cc M i - cc M j) ^ 2 * (SS M j - SS M i) * SS M ((j : ℝ) - i)
        / (1 - CC M ((j : ℝ) - i)) ^ 2)
      = nR * (nR - 1) ^ 2 / 4 := by
  rw [pair_weight M hM (fun a b => (cc M a - cc M b) ^ 2 * (SS M b - SS M a) * SS M (b - a))
    (fun a b => by simp only [cc_per, SS_per, perS M a b]) 2]
  have e : ∀ m ∈ Ico 1 (2 * M + 1), (∑ i ∈ range (2 * M + 1),
      (fun a b : ℝ => (cc M a - cc M b) ^ 2 * (SS M b - SS M a) * SS M (b - a)) i (i + m)) / kw M m ^ 2
      = (nR / 2) * (fa M (M - 1) m / kw M m ^ 2) + (-(nR / 2)) * (fa M M m / kw M m ^ 2)
        + (nR / 4) * (fa M 0 m / kw M m ^ 2) + (-(nR / 4)) * (fa M 2 m / kw M m ^ 2) := by
    intro m _
    have hc : ∑ i ∈ range (2 * M + 1),
        (fun a b : ℝ => (cc M a - cc M b) ^ 2 * (SS M b - SS M a) * SS M (b - a)) i (i + m)
        = (nR / 2) * fa M (M - 1) m + (-(nR / 2)) * fa M M m + (nR / 4) * fa M 0 m
          + (-(nR / 4)) * fa M 2 m := by
      simp only [add_sub_cancel_left]
      rw [← Finset.sum_mul, inner_g M hM m]
      have s2 := rel_s2 M m
      have c4 := rel_c4 M m
      have r2 := rel_r2 M m
      have r6 := rel_r6 M m
      have r7 := rel_r7 M m
      rw [← SS_m, ← fa_M, ← fa_1, ← fa_Mm1 M m (by omega), ← fa_2] at *
      rw [fa_0]
      linear_combination (-2 * nR * Real.sin ((M : ℝ) * hh M m) ^ 2 * SS M m) * r2
        + (4 * nR * Real.sin ((M : ℝ) * hh M m) ^ 2) * r6 + (nR * (fa M (M - 1) m - fa M M m)) * s2
        + (-nR / 2) * r7 + (nR / 2) * c4
    rw [hc]; ring
  rw [Finset.sum_congr rfl e]
  simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [F2' M (M - 1) (by omega), F2' M M (by omega), F2' M 0 (by omega), F2' M 2 (by omega),
    Nat.cast_sub (by omega : 1 ≤ M)]
  push_cast; unfold fejer2Poly; ring

theorem pairLin :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (SS M j - SS M i) * SS M ((j : ℝ) - i) / (1 - CC M ((j : ℝ) - i)) ^ 1)
      = 0 := by
  rw [pair_weight M hM (fun a b => (SS M b - SS M a) * SS M (b - a))
    (fun a b => by simp only [SS_per, perS M a b]) 1]
  refine Finset.sum_eq_zero (fun m _ => ?_)
  simp only [add_sub_cancel_left]
  rw [← Finset.sum_mul, inner_lin M hM m]; simp

theorem pair3a :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (cc M i - cc M j) ^ 2 * (cc M i + cc M j) / (1 - CC M ((j : ℝ) - i)) ^ 1)
      = 0 := by
  rw [pair_weight M hM (fun a b => (cc M a - cc M b) ^ 2 * (cc M a + cc M b))
    (fun a b => by simp only [cc_per]) 1]
  refine Finset.sum_eq_zero (fun m _ => ?_)
  rw [inner_3a M hM m]; simp

theorem pair3b :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else (cc M i - cc M j) * (CC M i - CC M j) / (1 - CC M ((j : ℝ) - i)) ^ 1)
      = 0 := by
  rw [pair_weight M hM (fun a b => (cc M a - cc M b) * (CC M a - CC M b))
    (fun a b => by simp only [cc_per, CC_per]) 1]
  refine Finset.sum_eq_zero (fun m _ => ?_)
  rw [inner_3b M hM m]; simp

end Pairs

end Havelock.Conj13.Odd4
