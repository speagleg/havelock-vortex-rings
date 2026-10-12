/-
Havelock — Conjecture 13: the inputs `h₂, h₄` (even `n`) and `f₂` (odd `n`) as iterated derivatives of
the Hamiltonian (paper eq. (1)) along the stated curves at the ring:

* `even_iteratedDeriv`  `d^k/dτ^k [H(curve((1+λ)τ))]|₀ = 0, 2h₂, 0, 24h₄` for `k = 1, 2, 3, 4`;
* `odd_iteratedDeriv`   `d^k/dA^k [H(oddLine A)]|₀ = 0, 2f₂` for `k = 1, 2`;

with `h₂ = h2 M λ`, `h₄ = h4 M λ`, `f₂ = f2 n λ` the definitions of `Conj13Algebra.lean`. They follow
from the expansions `EvenQuartic.Ham_curve_jet`, `Odd.Ham_oddLine_jet` by uniqueness of Taylor
coefficients (`iteratedDeriv_eq_of_isLittleO`, from Mathlib's `taylor_isLittleO`).
-/
import Havelock.Conj13.OddQuadratic
import Mathlib.Analysis.Calculus.Taylor
import Mathlib.Analysis.SpecialFunctions.Trigonometric.DerivHyp

open Real Asymptotics Filter Topology Finset Set

namespace Havelock.Conj13.Taylor

open Havelock.Conj13.Series

/-- A polynomial of degree `≤ n` that is `o(xⁿ)` on the punctured neighbourhood of `0` vanishes. -/
theorem poly_isLittleO_zero (n : ℕ) :
    ∀ c : ℕ → ℝ, (fun x : ℝ => ∑ k ∈ range (n + 1), c k * x ^ k) =o[𝓝[≠] 0] (fun x => x ^ n) →
      ∀ k ≤ n, c k = 0 := by
  induction n with
  | zero =>
    intro c h k hk
    have hk0 : k = 0 := by omega
    subst hk0
    simp only [zero_add, range_one, sum_singleton, pow_zero, mul_one] at h
    have := (isLittleO_one_iff ℝ).mp h
    exact tendsto_nhds_unique tendsto_const_nhds this
  | succ n ih =>
    intro c h
    -- the constant term
    have hcont : Continuous (fun x : ℝ => ∑ k ∈ range (n + 1 + 1), c k * x ^ k) := by fun_prop
    have hlim0 : Tendsto (fun x : ℝ => ∑ k ∈ range (n + 1 + 1), c k * x ^ k) (𝓝[≠] 0) (𝓝 0) := by
      have hp : Tendsto (fun x : ℝ => x ^ (n + 1)) (𝓝[≠] 0) (𝓝 0) := by
        have : Continuous (fun x : ℝ => x ^ (n + 1)) := by fun_prop
        simpa using (this.tendsto 0).mono_left nhdsWithin_le_nhds
      exact h.trans_tendsto hp
    have hlimc : Tendsto (fun x : ℝ => ∑ k ∈ range (n + 1 + 1), c k * x ^ k) (𝓝[≠] 0) (𝓝 (c 0)) := by
      have := (hcont.tendsto 0).mono_left (nhdsWithin_le_nhds (s := {0}ᶜ))
      convert this using 2
      rw [Finset.sum_range_succ']
      simp
    have hc0 : c 0 = 0 := tendsto_nhds_unique hlimc hlim0
    -- divide by `x`
    have hdiv : (fun x : ℝ => ∑ k ∈ range (n + 1), c (k + 1) * x ^ k) =o[𝓝[≠] 0] (fun x => x ^ n) := by
      have h1 := h.mul_isBigO (isBigO_refl (fun x : ℝ => x⁻¹) (𝓝[≠] 0))
      refine (h1.congr' ?_ ?_)
      · filter_upwards [self_mem_nhdsWithin] with x hx
        have hx' : x ≠ 0 := hx
        rw [Finset.sum_range_succ', hc0, zero_mul, add_zero, Finset.sum_mul]
        refine Finset.sum_congr rfl (fun k _ => ?_)
        rw [pow_succ]; field_simp
      · filter_upwards [self_mem_nhdsWithin] with x hx
        have hx' : x ≠ 0 := hx
        rw [pow_succ]; field_simp
    intro k hk
    rcases k with _ | k
    · exact hc0
    · exact ih (fun k => c (k + 1)) hdiv k (by omega)

/-- **Uniqueness of Taylor coefficients.** If `f` is `Cⁿ` on `(−δ, δ)` and
`f(x) = Σ_{k≤n} a_k x^k + o(xⁿ)`, then `f^{(k)}(0) = k! a_k` for `k ≤ n`. -/
theorem iteratedDeriv_eq_of_isLittleO {f : ℝ → ℝ} {δ : ℝ} (hδ : 0 < δ) (n : ℕ)
    (hf : ContDiffOn ℝ n f (Ioo (-δ) δ)) (a : ℕ → ℝ)
    (h : (fun x => f x - ∑ k ∈ range (n + 1), a k * x ^ k) =o[𝓝 0] (fun x => x ^ n)) :
    ∀ k ≤ n, iteratedDeriv k f 0 = (Nat.factorial k : ℝ) * a k := by
  have h0 : (0 : ℝ) ∈ Ioo (-δ) δ := ⟨by linarith, hδ⟩
  have hT := taylor_isLittleO (convex_Ioo (-δ) δ) h0 hf
  rw [isOpen_Ioo.nhdsWithin_eq h0] at hT
  simp only [sub_zero] at hT
  have hd := (h.sub hT).mono (nhdsWithin_le_nhds (s := {0}ᶜ))
  have hp := poly_isLittleO_zero n
    (fun k => (Nat.factorial k : ℝ)⁻¹ * iteratedDerivWithin k f (Ioo (-δ) δ) 0 - a k) ?_
  · intro k hk
    have := hp k hk
    rw [← iteratedDerivWithin_of_isOpen isOpen_Ioo h0]
    have hkf : (Nat.factorial k : ℝ) ≠ 0 := by positivity
    field_simp at this ⊢
    linarith
  · refine hd.congr' (Eventually.of_forall fun x => ?_) (Eventually.of_forall fun _ => rfl)
    simp only [taylor_within_apply, sub_zero, smul_eq_mul]
    rw [sub_sub_sub_cancel_left, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun k _ => ?_)
    ring

theorem pow_isLittleO {j k : ℕ} (h : j < k) : (fun x : ℝ => x ^ k) =o[𝓝 0] (fun x => x ^ j) := by
  have ht : Tendsto (fun x : ℝ => x ^ (k - j)) (𝓝 0) (𝓝 0) := by
    have : Continuous (fun x : ℝ => x ^ (k - j)) := by fun_prop
    simpa [zero_pow (by omega : k - j ≠ 0)] using this.tendsto 0
  have h1 := (isBigO_refl (fun x : ℝ => x ^ j) (𝓝 0)).mul_isLittleO ((isLittleO_one_iff ℝ).mpr ht)
  simp only [mul_one] at h1
  refine h1.congr_left (fun x => ?_)
  rw [← pow_add]; congr 1; omega

/-! ### Even `n = 2M` -/

open Havelock.Conj13.Config Havelock.Conj13.Even Havelock.Conj13.EvenQuartic in
/-- **Even `n = 2M`: `h₂` and `h₄` as derivatives of `H` along the symmetric curve.** With
`F(τ) = H(curve((1+λ)τ)) − H(ring)` (paper eq. (1); `J₃ = J₃(ring)`, `J₁ = J₂ = 0` on the curve,
`g(u±) = g(1) ± τ`): `F⁽¹⁾(0) = 0`, `F⁽²⁾(0) = 2h₂`, `F⁽³⁾(0) = 0`, `F⁽⁴⁾(0) = 24h₄`. -/
theorem even_iteratedDeriv (M : ℕ) (hM : M ≠ 0) (b : ℝ) (hb : |b| < 1) :
    iteratedDeriv 1 (fun τ => Ham b (curve M b ((1 + b) * τ)) - Ham b (stag M 1 1)) 0 = 0 ∧
    iteratedDeriv 2 (fun τ => Ham b (curve M b ((1 + b) * τ)) - Ham b (stag M 1 1)) 0
      = 2 * h2 M b ∧
    iteratedDeriv 3 (fun τ => Ham b (curve M b ((1 + b) * τ)) - Ham b (stag M 1 1)) 0 = 0 ∧
    iteratedDeriv 4 (fun τ => Ham b (curve M b ((1 + b) * τ)) - Ham b (stag M 1 1)) 0
      = 24 * h4 M b := by
  have hb' := abs_lt.mp hb
  have h1 : 0 < 1 + b := by linarith
  set δ := 1 / (1 + b) with hδ
  have hδpos : 0 < δ := by positivity
  set s := Ioo (-δ) δ with hs_def
  have hx : ∀ τ ∈ s, -1 < (1 + b) * τ ∧ (1 + b) * τ < 1 := by
    intro τ hτ
    obtain ⟨hl, hr⟩ := hτ
    have e : (1 + b) * δ = 1 := by rw [hδ]; field_simp
    constructor <;> nlinarith
  have hbx : ∀ τ ∈ s, -1 < b * ((1 + b) * τ) ∧ b * ((1 + b) * τ) < 1 := by
    intro τ hτ
    obtain ⟨hl, hr⟩ := hx τ hτ
    constructor <;> nlinarith
  have hEq : EqOn (fun τ => Ham b (curve M b ((1 + b) * τ)) - Ham b (stag M 1 1))
      (fun τ => Phi M b ((1 + b) * τ)) s := by
    intro τ hτ
    obtain ⟨hl, hr⟩ := hx τ hτ
    exact Ham_curve M hM b _ hb (abs_lt.mpr ⟨hl, hr⟩)
  have hlog : ∀ g : ℝ → ℝ, ContDiff ℝ 4 g → (∀ τ ∈ s, 0 < g τ) →
      ContDiffOn ℝ 4 (fun τ => Real.log (g τ)) s :=
    fun g hg hpos => hg.contDiffOn.log (fun τ hτ => (hpos τ hτ).ne')
  have hA : ContDiffOn ℝ 4 (fun τ => Real.log (1 - ((1 + b) * τ) ^ 2)) s :=
    hlog _ (by fun_prop) (fun τ hτ => by obtain ⟨hl, hr⟩ := hx τ hτ; nlinarith)
  have hB : ContDiffOn ℝ 4 (fun τ => Real.log (1 - (b * ((1 + b) * τ)) ^ 2)) s :=
    hlog _ (by fun_prop) (fun τ hτ => by obtain ⟨hl, hr⟩ := hbx τ hτ; nlinarith)
  have hath1 : ContDiffOn ℝ 4 (fun τ => ath ((1 + b) * τ)) s := by
    unfold ath
    exact ((hlog _ (by fun_prop) (fun τ hτ => by obtain ⟨hl, hr⟩ := hx τ hτ; linarith)).sub
      (hlog _ (by fun_prop) (fun τ hτ => by obtain ⟨hl, hr⟩ := hx τ hτ; linarith))).div_const 2
  have hath2 : ContDiffOn ℝ 4 (fun τ => ath (b * ((1 + b) * τ))) s := by
    unfold ath
    exact ((hlog _ (by fun_prop) (fun τ hτ => by obtain ⟨hl, hr⟩ := hbx τ hτ; linarith)).sub
      (hlog _ (by fun_prop) (fun τ hτ => by obtain ⟨hl, hr⟩ := hbx τ hτ; linarith))).div_const 2
  have hI : ContDiffOn ℝ 4 (fun τ => (M : ℝ) * (ath ((1 + b) * τ) + ath (b * ((1 + b) * τ))) / 2) s :=
    (contDiffOn_const.mul (hath1.add hath2)).div_const 2
  have hcosh : ContDiffOn ℝ 4
      (fun τ => Real.log (Real.cosh ((M : ℝ) * (ath ((1 + b) * τ) + ath (b * ((1 + b) * τ))) / 2))) s :=
    ((Real.contDiff_cosh (n := 4)).comp_contDiffOn hI).log (fun τ _ => (Real.cosh_pos _).ne')
  have hC : ContDiffOn ℝ 4 (fun τ => Phi M b ((1 + b) * τ)) s := by
    unfold Phi
    exact (contDiffOn_const.mul (hA.add hB)).sub (contDiffOn_const.mul hcosh)
  have hF := hC.congr hEq
  set a : ℕ → ℝ := fun k => if k = 2 then h2 M b else if k = 4 then h4 M b else 0 with ha
  have hjet : (fun τ => (Ham b (curve M b ((1 + b) * τ)) - Ham b (stag M 1 1))
      - ∑ k ∈ range (4 + 1), a k * τ ^ k) =o[𝓝 0] (fun τ => τ ^ 4) := by
    refine ((Ham_curve_jet M hM b hb).trans_isLittleO (pow_isLittleO (by norm_num))).congr_left
      (fun τ => ?_)
    simp only [Finset.sum_range_succ, Finset.range_zero, Finset.sum_empty, ha]
    norm_num
    ring
  have key := iteratedDeriv_eq_of_isLittleO hδpos 4 hF a hjet
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [key 1 (by norm_num)]; simp [ha]
  · rw [key 2 (by norm_num)]; simp [ha, Nat.factorial]
  · rw [key 3 (by norm_num)]; simp [ha]
  · rw [key 4 (by norm_num)]; simp [ha, Nat.factorial]

/-! ### Odd `n = 2M+1` -/

open Havelock.Conj13.Config Havelock.Conj13.OddSums Havelock.Conj13.Odd in
/-- **Odd `n = 2M+1 ≥ 3`: `f₂` as a derivative of `H` along the kernel line.** With
`F(A) = H(oddLine A) − H(ring)` (paper eq. (1); canonical coordinates `w_k = w₀ + A cos(lφ_k)`,
`ϑ_k = φ_k`): `F'(0) = 0` and `F''(0) = 2f₂`, for every `λ ∈ (−1, 1)`. -/
theorem odd_iteratedDeriv (M : ℕ) (hM : 1 ≤ M) (b : ℝ) (hb : -1 < b) (hb1 : b < 1) :
    iteratedDeriv 1 (fun A => Ham b (oddLine M b A) - Ham b (ringConf (2 * M + 1))) 0 = 0 ∧
    iteratedDeriv 2 (fun A => Ham b (oddLine M b A) - Ham b (ringConf (2 * M + 1))) 0
      = 2 * f2 ((2 * M + 1 : ℕ) : ℝ) b := by
  set n := 2 * M + 1 with hn
  set D := fun (i j : ℕ) (A : ℝ) => Dpair b (cc M i) (cc M j) (CC M ((j : ℝ) - i)) A with hD
  have hw0 : 0 < w0 b := by unfold w0; exact div_pos one_pos (by linarith)
  -- a ball on which every `D_ij` (`i ≠ j`) is positive and the line stays in the chart
  have hDcont : ∀ i j : ℕ, Continuous (D i j) := by
    intro i j; simp only [hD]; unfold Dpair Pw; fun_prop
  have hev : ∀ᶠ A in 𝓝 (0 : ℝ), |A| < w0 b ∧
      ∀ i ∈ range n, ∀ j ∈ range n, i ≠ j → 0 < D i j A := by
    refine Filter.Eventually.and ?_ ?_
    · exact continuous_abs.continuousAt.eventually_lt continuousAt_const (by simpa using hw0)
    · rw [Filter.eventually_all_finset]
      intro i hi
      rw [Filter.eventually_all_finset]
      intro j hj
      by_cases hij : i = j
      · exact Eventually.of_forall (fun _ h => absurd hij h)
      · have hC := CC_lt_one M i j (Finset.mem_range.mp hi) (Finset.mem_range.mp hj) hij
        have h0 : 0 < D i j 0 := by
          simp only [hD]; rw [Dpair_zero _ _ _ _ hb]
          have : 0 < 1 + b := by linarith
          have : 0 < 1 - CC M ((j : ℝ) - i) := by linarith
          positivity
        exact ((hDcont i j).tendsto 0 |>.eventually (lt_mem_nhds h0)).mono (fun A hA _ => hA)
  obtain ⟨δ, hδ, hball⟩ := Metric.eventually_nhds_iff_ball.mp hev
  have hs : Metric.ball (0 : ℝ) δ = Ioo (-δ) δ := by rw [Real.ball_eq_Ioo]; simp
  rw [hs] at hball
  set s := Ioo (-δ) δ with hs_def
  have hEq : EqOn (fun A => Ham b (oddLine M b A) - Ham b (ringConf n))
      (fun A => -(1 / 2) * (∑ i ∈ range n, ∑ j ∈ range n,
        (if i = j then 0 else Real.log (D i j A))) - Ham b (ringConf n)) s := by
    intro A hA
    simp only [hD]
    rw [Ham_oddLine M b A hb hb1 (hball A hA).1]
  have hDdiff : ∀ i j : ℕ, i ∈ range n → j ∈ range n → i ≠ j →
      ContDiffOn ℝ 2 (fun A => Real.log (D i j A)) s := by
    intro i j hi hj hij
    have hsq : ContDiffOn ℝ 2 (fun A => √(Pw b (w0 b + A * cc M i) * Pw b (w0 b + A * cc M j))) s := by
      refine ContDiffOn.sqrt ?_ ?_
      · unfold Pw; fun_prop
      · intro A hA
        obtain ⟨h1, h2⟩ := w_pos M b A (i : ℝ) hb hb1 (hball A hA).1
        obtain ⟨h3, h4⟩ := w_pos M b A (j : ℝ) hb hb1 (hball A hA).1
        unfold Pw; positivity
    have hpoly : ContDiffOn ℝ 2 (fun A : ℝ => (w0 b + A * cc M i) + (w0 b + A * cc M j)
        - 2 * b * (w0 b + A * cc M i) * (w0 b + A * cc M j)) s := by fun_prop
    have hDd : ContDiffOn ℝ 2 (D i j) s := by
      simp only [hD]; unfold Dpair
      exact hpoly.sub ((contDiffOn_const.mul hsq).mul contDiffOn_const)
    exact hDd.log (fun A hA => ((hball A hA).2 i hi j hj hij).ne')
  have hG : ContDiffOn ℝ 2 (fun A => -(1 / 2) * (∑ i ∈ range n, ∑ j ∈ range n,
      (if i = j then 0 else Real.log (D i j A))) - Ham b (ringConf n)) s := by
    refine (contDiffOn_const.mul (ContDiffOn.sum (fun i hi => ContDiffOn.sum (fun j hj => ?_)))).sub
      contDiffOn_const
    by_cases hij : i = j
    · simp only [hij, ite_true]; exact contDiffOn_const
    · simp only [hij, ite_false]; exact hDdiff i j hi hj hij
  have hF := hG.congr hEq
  set a : ℕ → ℝ := fun k => if k = 2 then f2 ((2 * M + 1 : ℕ) : ℝ) b else 0 with ha
  have hjet : (fun A => (Ham b (oddLine M b A) - Ham b (ringConf n))
      - ∑ k ∈ range (2 + 1), a k * A ^ k) =o[𝓝 0] (fun A => A ^ 2) := by
    refine ((Ham_oddLine_jet M hM b hb hb1).trans_isLittleO (pow_isLittleO (by norm_num))).congr_left
      (fun A => ?_)
    simp only [Finset.sum_range_succ, Finset.range_zero, Finset.sum_empty, ha]
    norm_num
    rfl
  have key := iteratedDeriv_eq_of_isLittleO hδ 2 hF a hjet
  refine ⟨?_, ?_⟩
  · rw [key 1 (by norm_num)]; simp [ha]
  · rw [key 2 (by norm_num)]; simp [ha, Nat.factorial]

/-! ### At the threshold `b = b_n` -/

open Havelock.Conj13.Config Havelock.Conj13.EvenQuartic in
/-- Even `n = 2M ≥ 4` at `b_n`: the second derivative of `H` along the symmetric curve vanishes
(the kernel direction), and the fourth is `24 h₄(b_n)`. -/
theorem even_iteratedDeriv_at_bn (M : ℕ) (hM : 2 ≤ M) :
    iteratedDeriv 2 (fun τ => Ham (bn (2 * M)) (curve M (bn (2 * M)) ((1 + bn (2 * M)) * τ))
      - Ham (bn (2 * M)) (stag M 1 1)) 0 = 0 ∧
    iteratedDeriv 4 (fun τ => Ham (bn (2 * M)) (curve M (bn (2 * M)) ((1 + bn (2 * M)) * τ))
      - Ham (bn (2 * M)) (stag M 1 1)) 0 = 24 * h4 M (bn (2 * M)) := by
  obtain ⟨h1, h2', -⟩ := bn_spec (2 * M) (by omega)
  obtain ⟨-, hd2, -, hd4⟩ := even_iteratedDeriv M (by omega) (bn (2 * M)) (abs_lt.mpr ⟨h1, h2'⟩)
  refine ⟨?_, hd4⟩
  rw [hd2, (c4_even_at_bn M hM).1, mul_zero]

open Havelock.Conj13.Config Havelock.Conj13.Odd in
/-- Odd `n = 2M+1 ≥ 5` at `b_n`: the second derivative of `H` along the kernel line vanishes. -/
theorem odd_iteratedDeriv_at_bn (M : ℕ) (hM : 2 ≤ M) :
    iteratedDeriv 2 (fun A => Ham (bn (2 * M + 1)) (oddLine M (bn (2 * M + 1)) A)
      - Ham (bn (2 * M + 1)) (ringConf (2 * M + 1))) 0 = 0 := by
  obtain ⟨h1, h2', -⟩ := bn_spec (2 * M + 1) (by omega)
  rw [(odd_iteratedDeriv M (by omega) _ h1 h2').2, (c4_odd_at_bn M hM).1, mul_zero]

end Havelock.Conj13.Taylor
