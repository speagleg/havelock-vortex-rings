import Mathlib.RingTheory.RootsOfUnity.Complex
import Mathlib.Analysis.SpecificLimits.Normed
import Havelock.Tanh.General

/-! # Tanh model, Fourier blocks: sums over the `N`-th roots of unity (proof of note Thm `thm:blocks`)

Finite sums over `l ∈ Fin N` of functions of `w^l`, `w = e^{2πi/N}`, used in the derivation of the
Fourier blocks (`Havelock/Tanh/Blocks.lean`):

* `rootSum`: `∑_l (w^l)^p conj(w^l)^k = N [p ≡ k (mod N)]` (roots-of-unity filter);
* `sum_zeta`: for fixed `j`, `k ↦ w^j conj(w^k)` runs over all `w^l` (circulant reindexing); `sum_conj`;
* `sum_rho1`: `∑_l ξw^l/(1−ξw^l) = N X/(1−X)`, `X = ξ^N` (proof of note Lemma `lem:re`);
* `sum_rho`: `∑_l ξw^l/(1−ξw^l)² · conj(w^l)^k = N W_k` for `0 ≤ k ≤ N` (`W_0 = W_N`) (the power-sum step of the
  proof of note Thm `thm:blocks`, via `General.W_hasSum`);
* `fejerC`: `∑_{l≠0} (1 − conj(w^l)^m)/|1−w^l|² = m(N−m)/2` (`0 ≤ m ≤ N`; Fejér-type identity, the pair
  term of note Thm `thm:blocks`); `sum_inv_one_sub`: `∑_{l≠0} 1/(1−w^l) = (N−1)/2` (pair term of
  note Lemma `lem:re`).
-/

noncomputable section

namespace Havelock.Tanh

open scoped ComplexConjugate

/-- `w = e^{2πi/N}` (note §"Coordinates and the relative equilibrium"). -/
def ringW (N : ℕ) : ℂ := Complex.exp (2 * Real.pi * Complex.I / N)

/-- `ρ(ζ) = ξζ/(1−ξζ)² = Σ_{p≥1} p (ξζ)^p` (second derivative of the image term). -/
def rho (ξ : ℝ) (ζ : ℂ) : ℂ := (ξ : ℂ) * ζ / (1 - (ξ : ℂ) * ζ) ^ 2

/-- `ρ₁(ζ) = ξζ/(1−ξζ) = Σ_{p≥1} (ξζ)^p` (first derivative of the image term). -/
def rho1 (ξ : ℝ) (ζ : ℂ) : ℂ := (ξ : ℂ) * ζ / (1 - (ξ : ℂ) * ζ)

theorem ringW_prim {N : ℕ} (hN : 0 < N) : IsPrimitiveRoot (ringW N) N :=
  Complex.isPrimitiveRoot_exp N hN.ne'

theorem ringW_ne_zero (N : ℕ) : ringW N ≠ 0 := Complex.exp_ne_zero _

theorem ringW_pow_N {N : ℕ} (hN : 0 < N) : ringW N ^ N = 1 := (ringW_prim hN).pow_eq_one

theorem conj_ringW (N : ℕ) : conj (ringW N) = (ringW N)⁻¹ := by
  unfold ringW
  rw [← Complex.exp_conj, ← Complex.exp_neg]
  congr 1
  simp only [map_div₀, map_mul, Complex.conj_ofReal, Complex.conj_I, map_natCast, map_ofNat]
  ring

theorem conj_pow_ringW (N l : ℕ) : conj (ringW N ^ l) = (ringW N ^ l)⁻¹ := by
  rw [map_pow, conj_ringW, inv_pow]

theorem ringW_pow_mul_conj (N l : ℕ) : ringW N ^ l * conj (ringW N ^ l) = 1 := by
  rw [conj_pow_ringW, mul_inv_cancel₀ (pow_ne_zero _ (ringW_ne_zero N))]

theorem ringW_pow_mod {N : ℕ} (hN : 0 < N) (p : ℕ) : ringW N ^ (p % N) = ringW N ^ p := by
  conv_rhs => rw [← Nat.mod_add_div p N, pow_add, pow_mul, ringW_pow_N hN, one_pow, mul_one]

theorem norm_ringW_pow {N : ℕ} (hN : 0 < N) (l : ℕ) : ‖ringW N ^ l‖ = 1 := by
  rw [norm_pow, (ringW_prim hN).norm'_eq_one hN.ne', one_pow]

theorem ringW_pow_eq_one_iff {N : ℕ} (hN : 0 < N) {l : ℕ} (hl : l < N) : ringW N ^ l = 1 ↔ l = 0 := by
  rw [(ringW_prim hN).pow_eq_one_iff_dvd]
  constructor
  · intro h; exact Nat.eq_zero_of_dvd_of_lt h hl
  · rintro rfl; exact dvd_zero N

/-- `conj(w^l)^k = (w^l)^{N−k}` for `k ≤ N`. -/
theorem conj_pow_eq {N : ℕ} (hN : 0 < N) (l k : ℕ) (hk : k ≤ N) :
    conj (ringW N ^ l) ^ k = (ringW N ^ l) ^ (N - k) := by
  have h1 : (ringW N ^ l) ^ (N - k) * (ringW N ^ l) ^ k = 1 := by
    rw [← pow_add, Nat.sub_add_cancel hk, ← pow_mul, mul_comm, pow_mul, ringW_pow_N hN, one_pow]
  rw [conj_pow_ringW, inv_pow, inv_eq_of_mul_eq_one_left h1]

theorem conj_ringW_pow_eq {N : ℕ} (hN : 0 < N) (k : ℕ) (hk : k ≤ N) :
    conj (ringW N ^ k) = ringW N ^ (N - k) := by
  have h := conj_pow_eq hN 1 k hk
  rw [pow_one, ← map_pow] at h
  exact h

/-- Roots-of-unity filter: `∑_{l<N} (w^l)^p conj(w^l)^k = N [p ≡ k mod N]`. -/
theorem rootSum {N : ℕ} (hN : 0 < N) (p k : ℕ) :
    ∑ l : Fin N, (ringW N ^ (l : ℕ)) ^ p * conj (ringW N ^ (l : ℕ)) ^ k
      = if p % N = k % N then (N : ℂ) else 0 := by
  set w := ringW N with hwdef
  have hw := ringW_prim hN
  have hw0 : w ≠ 0 := ringW_ne_zero N
  have key : ∀ l : ℕ, (w ^ l) ^ p * conj (w ^ l) ^ k = (w ^ p * (w ^ k)⁻¹) ^ l := by
    intro l
    rw [hwdef, conj_pow_ringW, mul_pow, ← pow_mul, ← pow_mul, inv_pow, inv_pow, ← pow_mul, ← pow_mul,
      mul_comm l p, mul_comm l k]
  rw [Finset.sum_congr rfl (fun (l : Fin N) _ => key l)]
  rw [Fin.sum_univ_eq_sum_range (fun l => (w ^ p * (w ^ k)⁻¹) ^ l) N]
  set u := w ^ p * (w ^ k)⁻¹ with hu
  have huN : u ^ N = 1 := by
    rw [hu, mul_pow, inv_pow, ← pow_mul, ← pow_mul, mul_comm p, mul_comm k, pow_mul, pow_mul,
      hwdef, ringW_pow_N hN]
    simp
  by_cases h : p % N = k % N
  · have hpk : w ^ p = w ^ k := by rw [hwdef, ← ringW_pow_mod hN, h, ringW_pow_mod hN]
    have : u = 1 := by rw [hu, hpk, mul_inv_cancel₀ (pow_ne_zero _ hw0)]
    simp [this, h]
  · have hu1 : u ≠ 1 := by
      intro h1
      have hpk : w ^ p = w ^ k := by
        have h2 := congrArg (· * w ^ k) h1
        simp only [hu, mul_assoc, inv_mul_cancel₀ (pow_ne_zero _ hw0), mul_one, one_mul] at h2
        exact h2
      rw [hwdef, ← ringW_pow_mod hN p, ← ringW_pow_mod hN k] at hpk
      exact h (hw.pow_inj (Nat.mod_lt _ hN) (Nat.mod_lt _ hN) hpk)
    rw [if_neg h, geom_sum_eq hu1, huN, sub_self, zero_div]

/-- `w^j conj(w^k) = w^{(j−k) mod N}`. -/
theorem zeta_eq {N : ℕ} [NeZero N] (j k : Fin N) :
    ringW N ^ (j : ℕ) * conj (ringW N ^ (k : ℕ)) = ringW N ^ ((j - k : Fin N) : ℕ) := by
  have hN : 0 < N := Nat.pos_of_ne_zero (NeZero.ne N)
  rw [Fin.val_sub, ringW_pow_mod hN, pow_add, conj_ringW_pow_eq hN k k.2.le, mul_comm]

theorem zeta_eq_one_iff {N : ℕ} [NeZero N] (j k : Fin N) :
    ringW N ^ (j : ℕ) * conj (ringW N ^ (k : ℕ)) = 1 ↔ j = k := by
  have hN : 0 < N := Nat.pos_of_ne_zero (NeZero.ne N)
  rw [zeta_eq, ringW_pow_eq_one_iff hN (j - k).2, ← sub_eq_zero (a := j), Fin.ext_iff]
  simp

/-- Circulant reindexing: for fixed `j`, `∑_k G(w^j conj(w^k)) = ∑_l G(w^l)`. -/
theorem sum_zeta {N : ℕ} [NeZero N] (j : Fin N) (G : ℂ → ℂ) :
    ∑ k : Fin N, G (ringW N ^ (j : ℕ) * conj (ringW N ^ (k : ℕ)))
      = ∑ l : Fin N, G (ringW N ^ (l : ℕ)) := by
  refine Fintype.sum_equiv (Equiv.subLeft j) _ _ (fun k => ?_)
  rw [zeta_eq]
  rfl

/-- `∑_l G(conj(w^l)) = ∑_l G(w^l)`. -/
theorem sum_conj {N : ℕ} [NeZero N] (G : ℂ → ℂ) :
    ∑ l : Fin N, G (conj (ringW N ^ (l : ℕ))) = ∑ l : Fin N, G (ringW N ^ (l : ℕ)) := by
  have := sum_zeta (0 : Fin N) G
  simpa using this

/-- Only `l = 0` has `w^l = 1`. -/
theorem sum_ite_one {N : ℕ} [NeZero N] (h : ℂ → ℂ) :
    ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then h (ringW N ^ (l : ℕ)) else 0) = h 1 := by
  have hN : 0 < N := Nat.pos_of_ne_zero (NeZero.ne N)
  rw [Finset.sum_eq_single (0 : Fin N)]
  · simp
  · intro l _ hl
    rw [if_neg]
    rw [ringW_pow_eq_one_iff hN l.2]
    intro h0; exact hl (Fin.ext (by simpa using h0))
  · intro h0; exact absurd (Finset.mem_univ _) h0

theorem sum_ite_ne_one {N : ℕ} [NeZero N] (h : ℂ → ℂ) :
    ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0 else h (ringW N ^ (l : ℕ)))
      = ∑ l : Fin N, h (ringW N ^ (l : ℕ)) - h 1 := by
  rw [← sum_ite_one (N := N) h, eq_sub_iff_add_eq, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun l _ => ?_
  split_ifs <;> simp

/-- Series version of the filter. -/
theorem hasSum_filter {N : ℕ} (hN : 0 < N) {c : ℕ → ℂ} {Φ : ℂ → ℂ}
    (h : ∀ l : Fin N, HasSum (fun p => c p * (ringW N ^ (l : ℕ)) ^ p) (Φ (ringW N ^ (l : ℕ)))) (k : ℕ) :
    HasSum (fun p => if p % N = k % N then (N : ℂ) * c p else 0)
      (∑ l : Fin N, Φ (ringW N ^ (l : ℕ)) * conj (ringW N ^ (l : ℕ)) ^ k) := by
  have := hasSum_sum (s := Finset.univ) (fun l _ => (h l).mul_right (conj (ringW N ^ (l : ℕ)) ^ k))
  convert this using 1
  funext p
  simp_rw [mul_assoc, ← Finset.mul_sum, rootSum hN]
  split_ifs <;> ring

theorem norm_xi_w {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) (l : ℕ) :
    ‖(ξ : ℂ) * ringW N ^ l‖ < 1 := by
  rw [norm_mul, norm_ringW_pow hN, mul_one, Complex.norm_real, Real.norm_of_nonneg h0]
  exact h1

/-- `∑_l ξw^l/(1−ξw^l)² · conj(w^l)^k = N W_k` for `0 ≤ k ≤ N` (`W_0 = W_N`) (note Thm `thm:blocks`, power sums). -/
theorem sum_rho {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) (k : ℕ) (hkN : k ≤ N) :
    ∑ l : Fin N, rho ξ (ringW N ^ (l : ℕ)) * conj (ringW N ^ (l : ℕ)) ^ k = ((N * W N k ξ : ℝ) : ℂ) := by
  have hser : ∀ l : Fin N, HasSum (fun p : ℕ => ((p : ℂ) * (ξ : ℂ) ^ p) * (ringW N ^ (l : ℕ)) ^ p)
      (rho ξ (ringW N ^ (l : ℕ))) := by
    intro l
    have := hasSum_coe_mul_geometric_of_norm_lt_one (norm_xi_w hN h0 h1 l)
    unfold rho
    convert this using 1
    funext p; rw [mul_pow]; ring
  have hA := hasSum_filter hN hser k
  -- the same series from `W_hasSum`
  have hW := (Complex.hasSum_ofReal.mpr (W_hasSum N k h0 h1 hN)).mul_left (N : ℂ)
  have hinj : Function.Injective (fun q : ℕ => k + q * N) := by
    intro a b hab
    have : a * N = b * N := by simpa using hab
    exact Nat.eq_of_mul_eq_mul_right hN this
  set f : ℕ → ℂ := fun p => if p % N = k % N then (N : ℂ) * ((p : ℂ) * (ξ : ℂ) ^ p) else 0 with hf
  have hcomp : (f ∘ fun q : ℕ => k + q * N) = fun q => (N : ℂ) * (((((k + q * N : ℕ) : ℝ) * ξ ^ (k + q * N)) : ℝ) : ℂ) := by
    funext q
    simp only [Function.comp, hf, Nat.add_mul_mod_self_right, if_true]
    push_cast; ring
  have hoff : ∀ p ∉ Set.range (fun q : ℕ => k + q * N), f p = 0 := by
    intro p hp
    simp only [hf]
    split_ifs with hpk
    · by_cases hle : k ≤ p
      · exfalso
        apply hp
        refine ⟨(p - k) / N, ?_⟩
        have hd : N ∣ p - k := (Nat.modEq_iff_dvd' hle).mp hpk.symm
        simp only
        rw [Nat.div_mul_cancel hd]; omega
      · -- p < k ≤ N
        rcases Nat.lt_or_ge k N with hk | hk
        · rw [Nat.mod_eq_of_lt (by omega), Nat.mod_eq_of_lt hk] at hpk; omega
        · have hkN' : k = N := le_antisymm hkN hk
          subst hkN'
          rw [Nat.mod_self, Nat.mod_eq_of_lt (by omega)] at hpk
          subst hpk
          simp
    · rfl
  have hB : HasSum f ((N : ℂ) * ((W N k ξ : ℝ) : ℂ)) := by
    rw [← hinj.hasSum_iff hoff, hcomp]
    exact hW
  rw [hA.unique hB]
  push_cast; ring

/-- `∑_l ξw^l/(1−ξw^l) = N X/(1−X)`, `X = ξ^N` (proof of note Lemma `lem:re`). -/
theorem sum_rho1 {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) :
    ∑ l : Fin N, rho1 ξ (ringW N ^ (l : ℕ)) = ((N * ξ ^ N / (1 - ξ ^ N) : ℝ) : ℂ) := by
  have hser : ∀ l : Fin N, HasSum (fun p : ℕ => ((ξ : ℂ) ^ p) * (ringW N ^ (l : ℕ)) ^ p)
      ((1 - (ξ : ℂ) * ringW N ^ (l : ℕ))⁻¹) := by
    intro l
    have := hasSum_geometric_of_norm_lt_one (norm_xi_w hN h0 h1 l)
    convert this using 1
    funext p; rw [mul_pow]
  have hA := hasSum_filter (Φ := fun ζ => (1 - (ξ : ℂ) * ζ)⁻¹) hN hser 0
  simp only [pow_zero, mul_one, Nat.zero_mod] at hA
  have hX0 : 0 ≤ ξ ^ N := pow_nonneg h0 N
  have hX1 : ξ ^ N < 1 := pow_lt_one₀ h0 h1 hN.ne'
  have hg := (Complex.hasSum_ofReal.mpr (hasSum_geometric_of_lt_one hX0 hX1)).mul_left (N : ℂ)
  have hinj : Function.Injective (fun q : ℕ => q * N) := fun a b hab =>
    Nat.eq_of_mul_eq_mul_right hN hab
  set f : ℕ → ℂ := fun p => if p % N = 0 then (N : ℂ) * (ξ : ℂ) ^ p else 0 with hf
  have hcomp : (f ∘ fun q : ℕ => q * N) = fun q => (N : ℂ) * ((((ξ ^ N) ^ q : ℝ)) : ℂ) := by
    funext q
    simp only [Function.comp, hf, Nat.mul_mod_left, if_true]
    push_cast; rw [← pow_mul, mul_comm q N]
  have hoff : ∀ p ∉ Set.range (fun q : ℕ => q * N), f p = 0 := by
    intro p hp
    simp only [hf]
    split_ifs with hpk
    · exfalso; apply hp; exact ⟨p / N, Nat.div_mul_cancel (Nat.dvd_of_mod_eq_zero hpk)⟩
    · rfl
  have hB : HasSum f ((N : ℂ) * (((1 - ξ ^ N)⁻¹ : ℝ) : ℂ)) := by
    rw [← hinj.hasSum_iff hoff, hcomp]; exact hg
  have hsum := hA.unique hB
  have hne : ∀ l : Fin N, (1 : ℂ) - (ξ : ℂ) * ringW N ^ (l : ℕ) ≠ 0 := by
    intro l h
    have := norm_xi_w hN h0 h1 l
    have h' : (ξ : ℂ) * ringW N ^ (l : ℕ) = 1 := (sub_eq_zero.mp h).symm
    rw [h', norm_one] at this; exact lt_irrefl _ this
  have hr : ∀ l : Fin N, rho1 ξ (ringW N ^ (l : ℕ)) = (1 - (ξ : ℂ) * ringW N ^ (l : ℕ))⁻¹ - 1 := by
    intro l; unfold rho1; field_simp [hne l]; ring
  rw [Finset.sum_congr rfl (fun l _ => hr l), Finset.sum_sub_distrib, hsum]
  have hX : (1 : ℂ) - (ξ : ℂ) ^ N ≠ 0 := by
    have : (1 : ℝ) - ξ ^ N ≠ 0 := by linarith
    exact_mod_cast this
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one]
  push_cast
  field_simp
  ring

/-- `∑_{l<N} |∑_{a<m} (w^l)^a|² = N m` for `m ≤ N`. -/
theorem sum_normSq_geomW {N : ℕ} (hN : 0 < N) (m : ℕ) (hm : m ≤ N) :
    ∑ l : Fin N, ((Complex.normSq (∑ a ∈ Finset.range m, (ringW N ^ (l : ℕ)) ^ a) : ℝ) : ℂ)
      = (N : ℂ) * m := by
  have e : ∀ l : Fin N, ((Complex.normSq (∑ a ∈ Finset.range m, (ringW N ^ (l : ℕ)) ^ a) : ℝ) : ℂ)
      = ∑ a ∈ Finset.range m, ∑ b ∈ Finset.range m,
          (ringW N ^ (l : ℕ)) ^ a * conj (ringW N ^ (l : ℕ)) ^ b := by
    intro l
    rw [← Complex.mul_conj, map_sum, Finset.sum_mul_sum]
    simp only [map_pow]
  rw [Finset.sum_congr rfl (fun l _ => e l), Finset.sum_comm]
  simp_rw [Finset.sum_comm (s := (Finset.univ : Finset (Fin N))) (t := Finset.range m)]
  simp_rw [rootSum hN]
  have inner : ∀ a ∈ Finset.range m,
      (∑ b ∈ Finset.range m, if a % N = b % N then (N : ℂ) else 0) = N := by
    intro a ha
    have ha' := Finset.mem_range.mp ha
    rw [Finset.sum_eq_single a]
    · simp
    · intro b hb hba
      have hb' := Finset.mem_range.mp hb
      rw [if_neg]
      rw [Nat.mod_eq_of_lt (by omega), Nat.mod_eq_of_lt (by omega)]
      exact fun h => hba h.symm
    · intro h; exact absurd ha h
  rw [Finset.sum_congr rfl inner]
  simp [mul_comm]

/-- Complex Fejér identity: `∑_{l≠0} (1 − conj(w^l)^m)/|1−w^l|² = m(N−m)/2` for `m ≤ N`
(pair term of note Thm `thm:blocks`; cf. `Conj13Algebra.fejer`). -/
theorem fejerC {N : ℕ} [NeZero N] (m : ℕ) (hm : m ≤ N) :
    ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0
        else (1 - conj (ringW N ^ (l : ℕ)) ^ m) / ((Complex.normSq (1 - ringW N ^ (l : ℕ)) : ℝ) : ℂ))
      = (m : ℂ) * ((N : ℂ) - m) / 2 := by
  have hN : 0 < N := Nat.pos_of_ne_zero (NeZero.ne N)
  set A := ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0
        else (1 - conj (ringW N ^ (l : ℕ)) ^ m) / ((Complex.normSq (1 - ringW N ^ (l : ℕ)) : ℝ) : ℂ))
    with hA
  -- conjugate reindexing
  have hA' : A = ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0
        else (1 - (ringW N ^ (l : ℕ)) ^ m) / ((Complex.normSq (1 - ringW N ^ (l : ℕ)) : ℝ) : ℂ)) := by
    rw [hA, ← sum_conj (fun ζ => if ζ = 1 then 0
        else (1 - conj ζ ^ m) / ((Complex.normSq (1 - ζ) : ℝ) : ℂ))]
    refine Finset.sum_congr rfl fun l _ => ?_
    have e1 : (conj (ringW N ^ (l : ℕ)) = 1) ↔ (ringW N ^ (l : ℕ) = 1) := by
      constructor
      · intro h; have := congrArg conj h; simpa using this
      · intro h; rw [h, map_one]
    simp only [e1, Complex.conj_conj]
    split_ifs
    · rfl
    · congr 2
      rw [← Complex.normSq_conj, map_sub, map_one, Complex.conj_conj]
  -- per-term identity
  have hterm : ∀ l : Fin N,
      (if ringW N ^ (l : ℕ) = 1 then 0
        else (1 - conj (ringW N ^ (l : ℕ)) ^ m) / ((Complex.normSq (1 - ringW N ^ (l : ℕ)) : ℝ) : ℂ))
      + (if ringW N ^ (l : ℕ) = 1 then 0
        else (1 - (ringW N ^ (l : ℕ)) ^ m) / ((Complex.normSq (1 - ringW N ^ (l : ℕ)) : ℝ) : ℂ))
      = ((Complex.normSq (∑ a ∈ Finset.range m, (ringW N ^ (l : ℕ)) ^ a) : ℝ) : ℂ)
        - (if ringW N ^ (l : ℕ) = 1 then (m : ℂ) ^ 2 else 0) := by
    intro l
    set ζ := ringW N ^ (l : ℕ) with hζ
    split_ifs with h1
    · rw [h1]; simp [Complex.normSq_natCast]; ring
    · have hne : (1 : ℂ) - ζ ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
      have hnsq : ((Complex.normSq (1 - ζ) : ℝ) : ℂ) ≠ 0 := by
        rw [Complex.ofReal_ne_zero]; exact (Complex.normSq_pos.mpr hne).ne'
      have hgeom : ∑ a ∈ Finset.range m, ζ ^ a = (ζ ^ m - 1) / (ζ - 1) := geom_sum_eq h1 m
      have hzc : ζ * conj ζ = 1 := ringW_pow_mul_conj N l
      rw [hgeom, Complex.normSq_div]
      push_cast
      rw [← Complex.mul_conj, ← Complex.mul_conj, ← Complex.mul_conj]
      have hc : conj (1 - ζ) ≠ 0 := by
        rw [map_sub, map_one]; intro h
        apply hne
        have := congrArg conj h; simpa using this
      have hc' : conj (ζ - 1) ≠ 0 := by
        rw [map_sub, map_one]; intro h
        apply hne
        have := congrArg conj h; simp at this; rw [← neg_sub]; simp [this]
      have hne' : ζ - 1 ≠ 0 := sub_ne_zero.mpr h1
      simp only [map_sub, map_one, map_pow] at hc hc' ⊢
      field_simp
      have hcm : conj ζ ^ m * ζ ^ m = 1 := by rw [← mul_pow, mul_comm, hzc, one_pow]
      linear_combination (1 - conj ζ) * (ζ - 1) * hcm
  have h2A : 2 * A = (N : ℂ) * m - (m : ℂ) ^ 2 := by
    have := Finset.sum_congr rfl (fun l (_ : l ∈ Finset.univ) => hterm l)
    rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, sum_normSq_geomW hN m hm, ← hA, ← hA'] at this
    rw [sum_ite_one (fun _ => (m : ℂ) ^ 2)] at this
    linear_combination this
  linear_combination h2A / 2

/-- `∑_{l≠0} 1/(1−w^l) = (N−1)/2` (pair term in the proof of note Lemma `lem:re`). -/
theorem sum_inv_one_sub {N : ℕ} [NeZero N] :
    ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0 else 1 / (1 - ringW N ^ (l : ℕ)))
      = ((N : ℂ) - 1) / 2 := by
  set T := ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0 else 1 / (1 - ringW N ^ (l : ℕ))) with hT
  have hT' : T = ∑ l : Fin N, (if ringW N ^ (l : ℕ) = 1 then 0 else 1 / (1 - conj (ringW N ^ (l : ℕ)))) := by
    rw [hT, ← sum_conj (fun ζ => if ζ = 1 then 0 else 1 / (1 - ζ))]
    refine Finset.sum_congr rfl fun l _ => ?_
    have e1 : (conj (ringW N ^ (l : ℕ)) = 1) ↔ (ringW N ^ (l : ℕ) = 1) := by
      constructor
      · intro h; have := congrArg conj h; simpa using this
      · intro h; rw [h, map_one]
    simp only [e1]
  have hterm : ∀ l : Fin N,
      (if ringW N ^ (l : ℕ) = 1 then 0 else 1 / (1 - ringW N ^ (l : ℕ)))
      + (if ringW N ^ (l : ℕ) = 1 then 0 else 1 / (1 - conj (ringW N ^ (l : ℕ))))
      = 1 - (if ringW N ^ (l : ℕ) = 1 then 1 else 0) := by
    intro l
    set ζ := ringW N ^ (l : ℕ) with hζ
    split_ifs with h1
    · ring
    · have hzc : ζ * conj ζ = 1 := ringW_pow_mul_conj N l
      have hne : (1 : ℂ) - ζ ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
      have hc : (1 : ℂ) - conj ζ ≠ 0 := by
        intro h; apply hne
        have := congrArg conj h; simpa using this
      field_simp
      linear_combination -hzc
  have h2 : 2 * T = (N : ℂ) - 1 := by
    have := Finset.sum_congr rfl (fun l (_ : l ∈ Finset.univ) => hterm l)
    rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← hT, ← hT'] at this
    rw [sum_ite_one (fun _ => (1 : ℂ))] at this
    simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one] at this
    linear_combination this
  linear_combination h2 / 2

end Havelock.Tanh
