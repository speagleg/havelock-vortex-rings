/-
Havelock — Conjecture 13, even `n = 2M`: the Hamiltonian (paper eq. (1)) on the configuration of two
staggered regular `M`-gons (paper, proof of Thm even, "Product formulas").

* `Ham`               paper eq. (1): `H = −Σ_{i<j} log(|z_i−z_j|²/(σ_iσ_j))`, `σ = 1 + λ|z|²`, written for
                      any finite index type as `−½ Σ_{i≠j}`; `Ham_fin` shows it is `−Σ_{i<j}` on `Fin n`,
                      `Ham_equiv` that it does not depend on the labelling.
* `prod_X_sub_odd`    `∏_{k<M} (X − ω^{2k+1}) = X^M + 1`, `ω = e^{iπ/M}` (the paper's
                      `∏_{ω^M=−1}(x − yω) = x^M + y^M`).
* `stag M ρ₊ ρ₋`      the configuration `{ρ₊η^j} ∪ {ρ₋ζη^j}` (`η = e^{2πi/M}`, `ζ = e^{iπ/M}`), indexed by
                      `Fin M × Fin 2`; `stag M 1 1` is the regular `2M`-gon (`stag_one_one`, relabelled by
                      `finProdFinEquiv`).
* `Ham_stag_sub`      `H(stag) − H(ring) = −(M(M−1)/2)(log u₊ + log u₋) − 2M log((ρ₊^M+ρ₋^M)/2)
                      + (2M−1)M(log σ₊ + log σ₋ − 2 log(1+λ))`, `u = ρ²`, `σ = 1 + λu`.
-/
import Havelock.Conj13.EvenCurve
import Mathlib.Algebra.Polynomial.Roots
import Mathlib.RingTheory.RootsOfUnity.Complex
import Mathlib.Analysis.SpecialFunctions.Complex.Log
import Mathlib.Logic.Equiv.Fin.Basic

open Real Complex Finset

namespace Havelock.Conj13.Config

/-- `σ(z) = 1 + λ|z|²`. -/
noncomputable def sig (b : ℝ) (z : ℂ) : ℝ := 1 + b * Complex.normSq z

/-- The pair term `log(|z−w|²/(σ(z)σ(w)))` of paper eq. (1). -/
noncomputable def pairTerm (b : ℝ) (z w : ℂ) : ℝ :=
  Real.log (Complex.normSq (z - w) / (sig b z * sig b w))

theorem pairTerm_comm (b : ℝ) (z w : ℂ) : pairTerm b z w = pairTerm b w z := by
  unfold pairTerm
  rw [← Complex.normSq_neg, neg_sub, mul_comm (sig b z)]

/-- Paper eq. (1), `H = −Σ_{i<j} log(|z_i−z_j|²/(σ_iσ_j))`, for an arbitrary finite index type,
written as `−½ Σ_{i ≠ j}` (see `Ham_fin` for the `i < j` form). -/
noncomputable def Ham {ι : Type*} [Fintype ι] [DecidableEq ι] (b : ℝ) (z : ι → ℂ) : ℝ :=
  -(1 / 2) * ∑ i, ∑ j, if i = j then 0 else pairTerm b (z i) (z j)

/-- On `Fin n`, `Ham` is literally `−Σ_{i<j}` (paper eq. (1)). -/
theorem Ham_fin {n : ℕ} (b : ℝ) (z : Fin n → ℂ) :
    Ham b z = -∑ i, ∑ j, if i < j then pairTerm b (z i) (z j) else 0 := by
  unfold Ham
  have hsplit : ∀ i j : Fin n, (if i = j then (0 : ℝ) else pairTerm b (z i) (z j))
      = (if i < j then pairTerm b (z i) (z j) else 0) + (if j < i then pairTerm b (z j) (z i) else 0) := by
    intro i j
    rcases lt_trichotomy i j with h | h | h
    · simp [h, h.ne, not_lt.mpr h.le]
    · simp [h]
    · simp [h, h.ne', not_lt.mpr h.le, pairTerm_comm]
  simp_rw [hsplit, Finset.sum_add_distrib]
  rw [Finset.sum_comm (f := fun i j => if j < i then pairTerm b (z j) (z i) else 0)]
  ring

/-- `Ham` is invariant under relabelling. -/
theorem Ham_equiv {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]
    (b : ℝ) (e : κ ≃ ι) (z : ι → ℂ) : Ham b (z ∘ e) = Ham b z := by
  unfold Ham
  congr 1
  rw [← e.sum_comp]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [← e.sum_comp]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  simp [e.injective.eq_iff]

/-! ### Roots of unity -/

/-- `ω = e^{iπ/M}`, a primitive `2M`-th root of unity. -/
noncomputable def om (M : ℕ) : ℂ := Complex.exp (↑(π / M) * I)

theorem normSq_om_pow (M k : ℕ) : Complex.normSq (om M ^ k) = 1 := by
  rw [map_pow, Complex.normSq_eq_norm_sq, om, Complex.norm_exp_ofReal_mul_I]; simp

theorem om_pow_M (M : ℕ) (hM : M ≠ 0) : om M ^ M = -1 := by
  rw [om, ← Complex.exp_nat_mul]
  have : (M : ℂ) * (↑(π / M) * I) = π * I := by
    push_cast; field_simp
  rw [this, Complex.exp_pi_mul_I]

theorem om_sq_prim (M : ℕ) (hM : M ≠ 0) : IsPrimitiveRoot (om M ^ 2) M := by
  have h := Complex.isPrimitiveRoot_exp M hM
  convert h using 1
  rw [om, ← Complex.exp_nat_mul]
  congr 1
  push_cast; field_simp

/-- `∏_{k<M} (X − ω^{2k+1}) = X^M + 1`: the `ω^{2k+1}` are the `M` roots of `t^M = −1`. -/
theorem prod_X_sub_odd (M : ℕ) (hM : M ≠ 0) (X : ℂ) :
    ∏ k : Fin M, (X - om M ^ (2 * (k : ℕ) + 1)) = X ^ M + 1 := by
  have hz := om_sq_prim M hM
  have hroots := hz.nthRoots_eq (α := om M) (a := -1) (om_pow_M M hM)
  set p : Polynomial ℂ := Polynomial.X ^ M - Polynomial.C (-1) with hp
  have hmon : p.Monic := Polynomial.monic_X_pow_sub_C _ hM
  have hr : p.roots = (Multiset.range M).map (fun k => (om M ^ 2) ^ k * om M) := by
    rw [← hroots]; rfl
  have hcard : Multiset.card p.roots = p.natDegree := by
    rw [hr, Multiset.card_map, Multiset.card_range, hp, Polynomial.natDegree_X_pow_sub_C]
  have hfac := Polynomial.prod_multiset_X_sub_C_of_monic_of_roots_card_eq hmon hcard
  have hev := congrArg (Polynomial.eval X) hfac
  rw [Polynomial.eval_multiset_prod, hr, Multiset.map_map, Multiset.map_map] at hev
  simp only [Function.comp_def, Polynomial.eval_sub, Polynomial.eval_X, Polynomial.eval_C] at hev
  rw [hp] at hev
  simp only [Polynomial.eval_sub, Polynomial.eval_pow, Polynomial.eval_X, Polynomial.eval_C,
    sub_neg_eq_add] at hev
  rw [← hev, Fin.prod_univ_eq_prod_range (fun k => X - om M ^ (2 * k + 1)) M]
  rw [Finset.prod_eq_multiset_prod]
  simp only [Finset.range_val]
  congr 1
  refine Multiset.map_congr rfl (fun k _ => ?_)
  rw [← pow_mul, ← pow_succ]

/-! ### The staggered configuration -/

/-- The radius of gon `s`: `ρ₊` for `s = 0`, `ρ₋` for `s = 1`. -/
noncomputable def rad (ρp ρm : ℝ) (s : Fin 2) : ℝ := if s = 0 then ρp else ρm

/-- Two staggered regular `M`-gons: vortex `(j, s)` (label `k = 2j + s`) sits at `ρ_s ω^{2j+s}`,
`ω = e^{iπ/M}`; i.e. `{ρ₊ η^j} ∪ {ρ₋ ζ η^j}` with `η = ω²`, `ζ = ω`. -/
noncomputable def stag (M : ℕ) (ρp ρm : ℝ) (p : Fin M × Fin 2) : ℂ :=
  (rad ρp ρm p.2 : ℂ) * om M ^ (2 * (p.1 : ℕ) + (p.2 : ℕ))

/-- The regular `n`-gon `z_k = e^{2πik/n}` (radius `r₀ = 1`). -/
noncomputable def ringConf (n : ℕ) (k : Fin n) : ℂ := Complex.exp (2 * π * I * (k : ℕ) / n)

/-- `stag M 1 1` is the regular `2M`-gon, relabelled by `(j, s) ↦ s + 2j`. -/
theorem stag_one_one (M : ℕ) (hM : M ≠ 0) : stag M 1 1 = ringConf (M * 2) ∘ finProdFinEquiv := by
  funext p
  simp only [stag, rad, Function.comp_apply, ringConf, finProdFinEquiv, Equiv.coe_fn_mk]
  rw [show (if p.2 = 0 then (1 : ℝ) else 1) = 1 by split_ifs <;> rfl, Complex.ofReal_one, one_mul,
    om, ← Complex.exp_nat_mul]
  congr 1
  push_cast
  field_simp
  ring

/-- Splitting the ordered-pair sum over `Fin M × Fin 2` into the two gons and the two cross terms. -/
theorem sum_pairs_split (M : ℕ) (F : Fin M × Fin 2 → Fin M × Fin 2 → ℝ) :
    (∑ a, ∑ c, if a = c then 0 else F a c)
      = (∑ j, ∑ j', if j = j' then 0 else F (j, 0) (j', 0))
        + (∑ j, ∑ j', if j = j' then 0 else F (j, 1) (j', 1))
        + (∑ j, ∑ j', F (j, 0) (j', 1)) + (∑ j, ∑ j', F (j, 1) (j', 0)) := by
  simp only [Fintype.sum_prod_type, Fin.sum_univ_two, Prod.mk.injEq, Finset.sum_add_distrib]
  simp
  ring

theorem sum_offdiag_const (M : ℕ) (c : ℝ) :
    (∑ j : Fin M, ∑ j' : Fin M, if j = j' then (0 : ℝ) else c) = M * (M - 1) * c := by
  have : ∀ j : Fin M, (∑ j' : Fin M, if j = j' then (0 : ℝ) else c) = (M - 1) * c := by
    intro j
    have h1 : ∀ j' : Fin M, (if j = j' then (0 : ℝ) else c) = c - (if j = j' then c else 0) := by
      intro j'; split_ifs <;> simp
    simp_rw [h1, Finset.sum_sub_distrib, Finset.sum_ite_eq]
    simp; ring
  simp_rw [this]; simp; ring

/-- The constant `Σ_{j≠j'} log|η^j − η^{j'}|²` of one regular `M`-gon. -/
noncomputable def Kgon (M : ℕ) : ℝ :=
  ∑ j : Fin M, ∑ j' : Fin M, if j = j' then 0 else
    Real.log (Complex.normSq (om M ^ (2 * (j : ℕ)) - om M ^ (2 * (j' : ℕ))))

theorem sig_rpow (b ρ : ℝ) (k : ℕ) (M : ℕ) : sig b ((ρ : ℂ) * om M ^ k) = 1 + b * ρ ^ 2 := by
  simp [sig, Complex.normSq_mul, normSq_om_pow, Complex.normSq_ofReal, sq]

/-- Pairs within one gon of radius `ρ`. -/
theorem same_gon (M : ℕ) (hM : M ≠ 0) (b ρ : ℝ) (hρ : 0 < ρ) (hσ : 0 < 1 + b * ρ ^ 2) (s : ℕ) :
    (∑ j : Fin M, ∑ j' : Fin M, if j = j' then 0 else
        pairTerm b ((ρ : ℂ) * om M ^ (2 * (j : ℕ) + s)) ((ρ : ℂ) * om M ^ (2 * (j' : ℕ) + s)))
      = M * (M - 1) * (Real.log (ρ ^ 2) - 2 * Real.log (1 + b * ρ ^ 2)) + Kgon M := by
  have hz := om_sq_prim M hM
  have hpt : ∀ j j' : Fin M, (if j = j' then 0 else
        pairTerm b ((ρ : ℂ) * om M ^ (2 * (j : ℕ) + s)) ((ρ : ℂ) * om M ^ (2 * (j' : ℕ) + s)))
      = (if j = j' then (0 : ℝ) else (Real.log (ρ ^ 2) - 2 * Real.log (1 + b * ρ ^ 2)))
        + (if j = j' then 0 else
            Real.log (Complex.normSq (om M ^ (2 * (j : ℕ)) - om M ^ (2 * (j' : ℕ))))) := by
    intro j j'
    split_ifs with h
    · simp
    have hne : om M ^ (2 * (j : ℕ)) - om M ^ (2 * (j' : ℕ)) ≠ 0 := by
      rw [sub_ne_zero, pow_mul, pow_mul]
      intro he
      exact h (Fin.ext (hz.pow_inj j.isLt j'.isLt he))
    have hd : Complex.normSq (om M ^ (2 * (j : ℕ)) - om M ^ (2 * (j' : ℕ))) ≠ 0 :=
      (Complex.normSq_pos.mpr hne).ne'
    unfold pairTerm
    rw [sig_rpow, sig_rpow]
    have he : (ρ : ℂ) * om M ^ (2 * (j : ℕ) + s) - (ρ : ℂ) * om M ^ (2 * (j' : ℕ) + s)
        = ((ρ : ℂ) * om M ^ s) * (om M ^ (2 * (j : ℕ)) - om M ^ (2 * (j' : ℕ))) := by
      rw [pow_add, pow_add]; ring
    rw [he, Complex.normSq_mul, Complex.normSq_mul, normSq_om_pow, Complex.normSq_ofReal, mul_one,
      Real.log_div (by positivity) (by positivity), Real.log_mul (by positivity) hd,
      Real.log_mul hσ.ne' hσ.ne', ← sq]
    ring
  simp_rw [hpt, Finset.sum_add_distrib, sum_offdiag_const]
  rfl

/-- The product over the opposite gon: `∏_{j'} (ρ₊ ω^{2j} − ρ₋ ω^{2j'+1}) = ρ₊^M + ρ₋^M`. -/
theorem cross_prod (M : ℕ) (hM : M ≠ 0) (ρp ρm : ℝ) (hρm : 0 < ρm) (j : ℕ) :
    ∏ j' : Fin M, ((ρp : ℂ) * om M ^ (2 * j) - (ρm : ℂ) * om M ^ (2 * (j' : ℕ) + 1))
      = ((ρp ^ M + ρm ^ M : ℝ) : ℂ) := by
  have hm : (ρm : ℂ) ≠ 0 := by exact_mod_cast hρm.ne'
  have h1 : ∀ j' : Fin M, ((ρp : ℂ) * om M ^ (2 * j) - (ρm : ℂ) * om M ^ (2 * (j' : ℕ) + 1))
      = (ρm : ℂ) * ((ρp : ℂ) * om M ^ (2 * j) / ρm - om M ^ (2 * (j' : ℕ) + 1)) := by
    intro j'; field_simp
  simp_rw [h1]
  rw [Finset.prod_mul_distrib, prod_X_sub_odd M hM, Finset.prod_const, Finset.card_univ,
    Fintype.card_fin]
  have hw : (om M ^ (2 * j)) ^ M = 1 := by
    rw [← pow_mul, show 2 * j * M = M * (2 * j) by ring, pow_mul, om_pow_M M hM, pow_mul]; simp
  rw [div_pow, mul_pow, hw]
  push_cast
  field_simp

/-- Pairs `(+, −)`: `Σ_{j,j'} log(|ρ₊ω^{2j} − ρ₋ω^{2j'+1}|²/(σ₊σ₋)) = 2M log(ρ₊^M+ρ₋^M) − M²(log σ₊ + log σ₋)`. -/
theorem cross_sum (M : ℕ) (hM : M ≠ 0) (b ρp ρm : ℝ) (hρp : 0 < ρp) (hρm : 0 < ρm)
    (hσp : 0 < 1 + b * ρp ^ 2) (hσm : 0 < 1 + b * ρm ^ 2) :
    (∑ j : Fin M, ∑ j' : Fin M,
        pairTerm b ((ρp : ℂ) * om M ^ (2 * (j : ℕ))) ((ρm : ℂ) * om M ^ (2 * (j' : ℕ) + 1)))
      = 2 * M * Real.log (ρp ^ M + ρm ^ M)
        - M ^ 2 * (Real.log (1 + b * ρp ^ 2) + Real.log (1 + b * ρm ^ 2)) := by
  have hS : 0 < ρp ^ M + ρm ^ M := by positivity
  have hrow : ∀ j : Fin M, (∑ j' : Fin M,
        pairTerm b ((ρp : ℂ) * om M ^ (2 * (j : ℕ))) ((ρm : ℂ) * om M ^ (2 * (j' : ℕ) + 1)))
      = 2 * Real.log (ρp ^ M + ρm ^ M) - M * (Real.log (1 + b * ρp ^ 2) + Real.log (1 + b * ρm ^ 2)) := by
    intro j
    have hprod := cross_prod M hM ρp ρm hρm j
    have hnorm : ∏ j' : Fin M,
        Complex.normSq ((ρp : ℂ) * om M ^ (2 * (j : ℕ)) - (ρm : ℂ) * om M ^ (2 * (j' : ℕ) + 1))
          = (ρp ^ M + ρm ^ M) ^ 2 := by
      rw [← map_prod, hprod, Complex.normSq_ofReal, sq]
    have hnz : ∀ j' ∈ (Finset.univ : Finset (Fin M)),
        Complex.normSq ((ρp : ℂ) * om M ^ (2 * (j : ℕ)) - (ρm : ℂ) * om M ^ (2 * (j' : ℕ) + 1)) ≠ 0 := by
      have : (∏ j' : Fin M, Complex.normSq
          ((ρp : ℂ) * om M ^ (2 * (j : ℕ)) - (ρm : ℂ) * om M ^ (2 * (j' : ℕ) + 1))) ≠ 0 := by
        rw [hnorm]; positivity
      exact Finset.prod_ne_zero_iff.mp this
    have hterm : ∀ j' : Fin M,
        pairTerm b ((ρp : ℂ) * om M ^ (2 * (j : ℕ))) ((ρm : ℂ) * om M ^ (2 * (j' : ℕ) + 1))
          = Real.log (Complex.normSq
              ((ρp : ℂ) * om M ^ (2 * (j : ℕ)) - (ρm : ℂ) * om M ^ (2 * (j' : ℕ) + 1)))
            - (Real.log (1 + b * ρp ^ 2) + Real.log (1 + b * ρm ^ 2)) := by
      intro j'
      unfold pairTerm
      rw [sig_rpow, sig_rpow, Real.log_div (hnz j' (Finset.mem_univ _)) (by positivity),
        Real.log_mul hσp.ne' hσm.ne']
    simp_rw [hterm]
    rw [Finset.sum_sub_distrib, ← Real.log_prod hnz, hnorm, Real.log_pow]
    simp; ring
  simp_rw [hrow]
  simp
  ring

/-- `H` on the staggered configuration (paper, proof of Thm even, "Product formulas"). -/
theorem Ham_stag (M : ℕ) (hM : M ≠ 0) (b ρp ρm : ℝ) (hρp : 0 < ρp) (hρm : 0 < ρm)
    (hσp : 0 < 1 + b * ρp ^ 2) (hσm : 0 < 1 + b * ρm ^ 2) :
    Ham b (stag M ρp ρm) = -(1 / 2) *
      (M * (M - 1) * (Real.log (ρp ^ 2) - 2 * Real.log (1 + b * ρp ^ 2))
        + M * (M - 1) * (Real.log (ρm ^ 2) - 2 * Real.log (1 + b * ρm ^ 2)) + 2 * Kgon M
        + 2 * (2 * M * Real.log (ρp ^ M + ρm ^ M)
          - M ^ 2 * (Real.log (1 + b * ρp ^ 2) + Real.log (1 + b * ρm ^ 2)))) := by
  unfold Ham
  rw [sum_pairs_split]
  have h0 := same_gon M hM b ρp hρp hσp 0
  have h1 := same_gon M hM b ρm hρm hσm 1
  have hc := cross_sum M hM b ρp ρm hρp hρm hσp hσm
  have hc' : (∑ j : Fin M, ∑ j' : Fin M,
      pairTerm b ((ρm : ℂ) * om M ^ (2 * (j : ℕ) + 1)) ((ρp : ℂ) * om M ^ (2 * (j' : ℕ))))
        = 2 * M * Real.log (ρp ^ M + ρm ^ M)
          - M ^ 2 * (Real.log (1 + b * ρp ^ 2) + Real.log (1 + b * ρm ^ 2)) := by
    rw [Finset.sum_comm, ← hc]
    simp_rw [pairTerm_comm b ((ρm : ℂ) * _)]
  simp only [add_zero] at h0
  simp only [stag, rad, Fin.isValue, Fin.val_zero, Fin.val_one, one_ne_zero,
    add_zero, ↓reduceIte]
  rw [h0, h1, hc, hc']
  ring

/-- `H(stag) − H(ring)`, the constant `Kgon` cancelling. -/
theorem Ham_stag_sub (M : ℕ) (hM : M ≠ 0) (b ρp ρm : ℝ) (hb : -1 < b) (hρp : 0 < ρp)
    (hρm : 0 < ρm) (hσp : 0 < 1 + b * ρp ^ 2) (hσm : 0 < 1 + b * ρm ^ 2) :
    Ham b (stag M ρp ρm) - Ham b (stag M 1 1)
      = -(M * (M - 1) / 2) * (Real.log (ρp ^ 2) + Real.log (ρm ^ 2))
        - 2 * M * Real.log ((ρp ^ M + ρm ^ M) / 2)
        + (2 * M - 1) * M * (Real.log (1 + b * ρp ^ 2) + Real.log (1 + b * ρm ^ 2)
          - 2 * Real.log (1 + b)) := by
  have hb1 : (0 : ℝ) < 1 + b * 1 ^ 2 := by linarith
  rw [Ham_stag M hM b ρp ρm hρp hρm hσp hσm, Ham_stag M hM b 1 1 one_pos one_pos hb1 hb1,
    Real.log_div (by positivity) (by norm_num)]
  simp only [one_pow, mul_one, Real.log_one]
  norm_num
  ring

end Havelock.Conj13.Config
