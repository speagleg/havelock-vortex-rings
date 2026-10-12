/-
Havelock — Conjecture 13, even `n = 2M`: `h₂`, `h₄` and `c₄` derived from the Hamiltonian
(paper `conj13/paper/main.tex`, proof of Thm even).

The curve `Fix(D_M) ∩ J⁻¹(μ)`: two staggered `M`-gons with squared radii `u±`, parametrised by
`g(u±) = g(1) ± τ`, `g(u) = u/(1+λu)`; with `x = (1+λ)τ` this is `u± = (1±x)/(1∓λx)` (`curve`).

* `J3_curve`, `Jp_stag`    the curve lies in the momentum level set: `J₃ = J₃(ring)`, `J₁ + iJ₂ = 0`;
                           `g_uP`, `g_uM`: `g(u±) = g(1) ± τ`.
* `Ham_curve`              `H(curve x) − H(ring) = Phi M λ x` (closed form, `EvenCurve.lean`).
* `Ham_curve_jet`          **`H(curve((1+λ)τ)) − H(ring) = h₂τ² + h₄τ⁴ + O(τ⁶)`** with `h₂ = h2 M λ`,
                           `h₄ = h4 M λ` the definitions of `Conj13Algebra.lean`.
* `epsK_stag`              the component along the unit kernel vector `e_k = n^{-1/2}(−1)^k ẑ_k` is
                           exactly `ε = (√n/2)(ρ₊ − ρ₋)`; `epsK_curve_slope`: `ε(τ)/τ → (√n/2)(1+λ)²`.
* `c4even_limit`           at the threshold (`h₂ = 0`): `(H(curve) − H(ring))/ε⁴ → c4evenOf M λ`
                           (`= h₄(2/(1+λ)²)⁴/n²`, the definition used in `Conj13Algebra.lean`).
-/
import Havelock.Conj13.EvenConfig

open Real Asymptotics Filter Topology

namespace Havelock.Conj13.EvenQuartic

open Havelock.Conj13.Config Havelock.Conj13.Even Havelock.Conj13.Series

/-- `u₊(x) = (1+x)/(1−λx)`. -/
noncomputable def uP (b x : ℝ) : ℝ := (1 + x) / (1 - b * x)
/-- `u₋(x) = (1−x)/(1+λx)`. -/
noncomputable def uM (b x : ℝ) : ℝ := (1 - x) / (1 + b * x)

/-- The symmetric curve: radii `ρ± = √u±(x)`. -/
noncomputable def curve (M : ℕ) (b x : ℝ) : Fin M × Fin 2 → ℂ :=
  stag M (√(uP b x)) (√(uM b x))

/-- `g(u) = u/(1+λu)` (so that `J₃ = Σ g(|z_k|²)`). -/
noncomputable def gfun (b u : ℝ) : ℝ := u / (1 + b * u)

theorem g_uP (b x : ℝ) (hb1 : 1 + b ≠ 0) (h : 1 - b * x ≠ 0) :
    gfun b (uP b x) = gfun b 1 + x / (1 + b) := by
  unfold gfun uP
  have h2 : 1 + b * ((1 + x) / (1 - b * x)) = (1 + b) / (1 - b * x) := by field_simp; ring
  rw [h2, mul_one, div_div_div_cancel_right₀ h]
  field_simp

theorem g_uM (b x : ℝ) (hb1 : 1 + b ≠ 0) (h : 1 + b * x ≠ 0) :
    gfun b (uM b x) = gfun b 1 - x / (1 + b) := by
  unfold gfun uM
  have h2 : 1 + b * ((1 - x) / (1 + b * x)) = (1 + b) / (1 + b * x) := by field_simp; ring
  rw [h2, mul_one, div_div_div_cancel_right₀ h]
  field_simp

/-- `J₃ = Σ_k |z_k|²/σ_k` (paper eq. (1)). -/
noncomputable def J3 {ι : Type*} [Fintype ι] (b : ℝ) (z : ι → ℂ) : ℝ :=
  ∑ k, Complex.normSq (z k) / sig b (z k)

/-- `J₁ + iJ₂ = Σ_k z_k/σ_k` (paper eq. (1)). -/
noncomputable def Jp {ι : Type*} [Fintype ι] (b : ℝ) (z : ι → ℂ) : ℂ :=
  ∑ k, z k / (sig b (z k) : ℂ)

theorem J3_stag (M : ℕ) (b ρp ρm : ℝ) :
    J3 b (stag M ρp ρm) = M * (gfun b (ρp ^ 2) + gfun b (ρm ^ 2)) := by
  unfold J3 gfun
  simp only [Fintype.sum_prod_type, Fin.sum_univ_two, stag, rad, Fin.isValue, ↓reduceIte,
    one_ne_zero, sig_rpow, Complex.normSq_mul, normSq_om_pow, Complex.normSq_ofReal, mul_one]
  simp [sq]
  ring

/-- The curve lies in `J₃⁻¹(J₃(ring))`. -/
theorem J3_curve (M : ℕ) (b τ : ℝ) (hb : -1 < b) (hx : |(1 + b) * τ| < 1) (hb1 : b < 1) :
    J3 b (curve M b ((1 + b) * τ)) = J3 b (stag M 1 1) := by
  set x := (1 + b) * τ with hxdef
  have hbx : |b * x| < 1 := by
    rw [abs_mul]; calc |b| * |x| ≤ 1 * |x| := by gcongr; exact abs_le.mpr ⟨by linarith, by linarith⟩
      _ < 1 := by linarith
  rw [abs_lt] at hx hbx
  have hup : 0 ≤ uP b x := by unfold uP; apply div_nonneg <;> linarith
  have hum : 0 ≤ uM b x := by unfold uM; apply div_nonneg <;> linarith
  unfold curve
  have hb1' : 1 + b ≠ 0 := by linarith
  rw [J3_stag, J3_stag, Real.sq_sqrt hup, Real.sq_sqrt hum,
    g_uP b x hb1' (by linarith), g_uM b x hb1' (by linarith)]
  norm_num

/-- `J₁ = J₂ = 0` on every staggered configuration (`M ≥ 2`). -/
theorem Jp_stag (M : ℕ) (hM : 2 ≤ M) (b ρp ρm : ℝ) : Jp b (stag M ρp ρm) = 0 := by
  have hz := om_sq_prim M (by omega)
  have hgeom : ∑ j : Fin M, (om M ^ 2) ^ (j : ℕ) = 0 := by
    rw [Fin.sum_univ_eq_sum_range (fun j => (om M ^ 2) ^ j) M]
    exact hz.geom_sum_eq_zero (by omega)
  unfold Jp
  have hterm : ∀ (j : Fin M) (s : Fin 2), stag M ρp ρm (j, s) / (sig b (stag M ρp ρm (j, s)) : ℂ)
      = ((rad ρp ρm s : ℂ) / ((1 + b * rad ρp ρm s ^ 2 : ℝ) : ℂ) * om M ^ (s : ℕ))
        * (om M ^ 2) ^ (j : ℕ) := by
    intro j s
    unfold stag
    rw [sig_rpow, pow_add, pow_mul]
    ring
  rw [Fintype.sum_prod_type, Finset.sum_comm]
  simp_rw [hterm, ← Finset.mul_sum, hgeom, mul_zero, Finset.sum_const_zero]

/-- `(e^a + e^b)/2 = e^{(a+b)/2} cosh((a−b)/2)`. -/
theorem exp_add_exp_half (a c : ℝ) :
    (Real.exp a + Real.exp c) / 2 = Real.exp ((a + c) / 2) * Real.cosh ((a - c) / 2) := by
  rw [Real.cosh_eq, mul_div_assoc', mul_add, ← Real.exp_add, ← Real.exp_add]
  congr 2 <;> congr 1 <;> ring

theorem log_half_sum_pow (M : ℕ) (u v : ℝ) (hu : 0 < u) (hv : 0 < v) :
    Real.log ((√u ^ M + √v ^ M) / 2)
      = M * (Real.log u + Real.log v) / 4 + Real.log (Real.cosh (M * (Real.log u - Real.log v) / 4)) := by
  have hs : ∀ w : ℝ, 0 < w → √w ^ M = Real.exp (M * Real.log w / 2) := by
    intro w hw
    rw [← Real.exp_log (Real.sqrt_pos.mpr hw), ← Real.exp_nat_mul, Real.log_sqrt hw.le]
    ring_nf
  rw [hs u hu, hs v hv, exp_add_exp_half, Real.log_mul (Real.exp_pos _).ne'
    (Real.cosh_pos _).ne', Real.log_exp]
  congr 1
  · ring
  · congr 2; ring

/-- **The reduced Hamiltonian on the curve**: `H(curve x) − H(ring) = Phi M λ x` for `|x| < 1`,
`|λ| < 1`. Here `H` is paper eq. (1) (`Config.Ham`), so this derives the closed form `Phi`. -/
theorem Ham_curve (M : ℕ) (hM : M ≠ 0) (b x : ℝ) (hb : |b| < 1) (hx : |x| < 1) :
    Ham b (curve M b x) - Ham b (stag M 1 1) = Phi M b x := by
  have hbx : |b * x| < 1 := by
    rw [abs_mul]; calc |b| * |x| ≤ 1 * |x| := by gcongr
      _ < 1 := by linarith
  rw [abs_lt] at hx hbx hb
  have hxp : 0 < 1 + x := by linarith
  have hxm : 0 < 1 - x := by linarith
  have hbxp : 0 < 1 + b * x := by linarith
  have hbxm : 0 < 1 - b * x := by linarith
  have hb1 : 0 < 1 + b := by linarith
  have hup : 0 < uP b x := div_pos hxp hbxm
  have hum : 0 < uM b x := div_pos hxm hbxp
  have hσp : 1 + b * √(uP b x) ^ 2 = (1 + b) / (1 - b * x) := by
    rw [Real.sq_sqrt hup.le]; unfold uP; field_simp; ring
  have hσm : 1 + b * √(uM b x) ^ 2 = (1 + b) / (1 + b * x) := by
    rw [Real.sq_sqrt hum.le]; unfold uM; field_simp; ring
  unfold curve
  rw [Ham_stag_sub M hM b _ _ (by linarith) (Real.sqrt_pos.mpr hup) (Real.sqrt_pos.mpr hum)
      (by rw [hσp]; positivity) (by rw [hσm]; positivity),
    log_half_sum_pow M _ _ hup hum, hσp, hσm, Real.sq_sqrt hup.le, Real.sq_sqrt hum.le]
  have hD : (M : ℝ) * (Real.log (uP b x) - Real.log (uM b x)) / 4
      = M * (ath x + ath (b * x)) / 2 := by
    unfold uP uM ath
    rw [Real.log_div hxp.ne' hbxm.ne', Real.log_div hxm.ne' hbxp.ne']
    ring
  rw [hD]
  unfold Phi uP uM
  rw [Real.log_div hxp.ne' hbxm.ne', Real.log_div hxm.ne' hbxp.ne',
    Real.log_div hb1.ne' hbxm.ne', Real.log_div hb1.ne' hbxp.ne',
    show 1 - x ^ 2 = (1 - x) * (1 + x) by ring, show 1 - (b * x) ^ 2 = (1 - b * x) * (1 + b * x) by ring,
    Real.log_mul hxm.ne' hxp.ne', Real.log_mul hbxm.ne' hbxp.ne']
  ring

/-- **Even `n = 2M`, Taylor coefficients from the Hamiltonian.** Along the symmetric curve
(`J₃ = J₃(ring)`, `J₁ = J₂ = 0`, parameter `τ` with `g(u±) = g(1) ± τ`),
`H(curve) − H(ring) = h₂ τ² + h₄ τ⁴ + O(τ⁶)`, with `h₂ = h2 M λ`, `h₄ = h4 M λ` exactly the
expressions stated in the paper's proof of Thm even (definitions in `Conj13Algebra.lean`). -/
theorem Ham_curve_jet (M : ℕ) (hM : M ≠ 0) (b : ℝ) (hb : |b| < 1) :
    (fun τ => Ham b (curve M b ((1 + b) * τ)) - Ham b (stag M 1 1)
        - h2 M b * τ ^ 2 - h4 M b * τ ^ 4) =O[𝓝 0] (fun τ => τ ^ 6) := by
  refine (Phi_jet M b).congr' ?_ (Eventually.of_forall fun _ => rfl)
  have hc : Continuous (fun τ : ℝ => |(1 + b) * τ|) := by fun_prop
  have hev : ∀ᶠ τ in 𝓝 (0 : ℝ), |(1 + b) * τ| < 1 :=
    hc.continuousAt.eventually_lt continuousAt_const (by simp)
  filter_upwards [hev] with τ hτ
  rw [Ham_curve M hM b _ hb hτ]

/-! ### The kernel coordinate `ε` and the quartic coefficient -/

/-- The component `⟨q − x₀, e⟩` of a configuration `q` along the Euclidean-unit kernel vector
`e_k = n^{-1/2}(−1)^k ẑ_k` (paper Def. c4, even `n = 2M`; label `k = 2j + s`, so `(−1)^k = (−1)^s`),
with `⟨a, b⟩ = Re(ā b)` on each `ℂ ≅ ℝ²` factor and `ẑ_k = x₀,k` the unit ring positions. -/
noncomputable def epsK (M : ℕ) (z : Fin M × Fin 2 → ℂ) : ℝ :=
  ∑ p, (-1 : ℝ) ^ (p.2 : ℕ) / √(2 * M) * ((starRingEnd ℂ) (stag M 1 1 p) * (z p - stag M 1 1 p)).re

/-- The component along the rotation direction `t_k = i ẑ_k` (paper Def. c4: `⟨q − x₀, t⟩ = 0`). -/
noncomputable def tK (M : ℕ) (z : Fin M × Fin 2 → ℂ) : ℝ :=
  ∑ p, ((starRingEnd ℂ) (Complex.I * stag M 1 1 p) * (z p - stag M 1 1 p)).re

theorem conj_mul_stag (M : ℕ) (ρp ρm : ℝ) (p : Fin M × Fin 2) :
    (starRingEnd ℂ) (stag M 1 1 p) * (stag M ρp ρm p - stag M 1 1 p) = ((rad ρp ρm p.2 - 1 : ℝ) : ℂ) := by
  have hr1 : rad 1 1 p.2 = 1 := by unfold rad; split_ifs <;> rfl
  unfold stag
  rw [hr1]
  simp only [Complex.ofReal_one, one_mul]
  have hn : (starRingEnd ℂ) (om M ^ (2 * (p.1 : ℕ) + (p.2 : ℕ))) * om M ^ (2 * (p.1 : ℕ) + (p.2 : ℕ)) = 1 := by
    rw [mul_comm, Complex.mul_conj, normSq_om_pow]; simp
  generalize om M ^ (2 * (p.1 : ℕ) + (p.2 : ℕ)) = w at hn ⊢
  push_cast
  linear_combination ((rad ρp ρm p.2 : ℂ) - 1) * hn

/-- **`ε` on the staggered configuration**: exactly `ε = (√n/2)(ρ₊ − ρ₋)`, `n = 2M` (the paper's
`ε = (√n/2)(ρ₊−ρ₋) + O(ρ−1)²` holds with no error term). Also `⟨q − x₀, t⟩ = 0`. -/
theorem epsK_stag (M : ℕ) (hM : M ≠ 0) (ρp ρm : ℝ) :
    epsK M (stag M ρp ρm) = √(2 * M) / 2 * (ρp - ρm) ∧ tK M (stag M ρp ρm) = 0 := by
  have hs : (0 : ℝ) < √(2 * M) := Real.sqrt_pos.mpr (by positivity)
  have hsq : √(2 * (M : ℝ)) ^ 2 = 2 * M := Real.sq_sqrt (by positivity)
  constructor
  · unfold epsK
    simp_rw [conj_mul_stag, Complex.ofReal_re]
    simp only [Fintype.sum_prod_type, Fin.sum_univ_two, rad, Fin.isValue, ↓reduceIte, one_ne_zero,
      Fin.val_zero, Fin.val_one, pow_zero, pow_one, Finset.sum_const, Finset.card_univ,
      Fintype.card_fin, nsmul_eq_mul]
    have hsq2 : √((M : ℝ) * 2) ^ 2 = M * 2 := Real.sq_sqrt (by positivity)
    field_simp
    linear_combination (-(ρp - ρm)) * hsq2
  · unfold tK
    have : ∀ p : Fin M × Fin 2, ((starRingEnd ℂ) (Complex.I * stag M 1 1 p)
        * (stag M ρp ρm p - stag M 1 1 p)).re = 0 := by
      intro p
      rw [map_mul, mul_assoc, conj_mul_stag]
      simp
    exact Finset.sum_eq_zero (fun p _ => this p)

/-- `ε` along the curve has slope `κ = (√n/2)(1+λ)²` at `τ = 0`. -/
theorem epsK_curve_hasDerivAt (M : ℕ) (hM : M ≠ 0) (b : ℝ) :
    HasDerivAt (fun τ => epsK M (curve M b ((1 + b) * τ))) (√(2 * M) / 2 * (1 + b) ^ 2) 0 := by
  have hP : HasDerivAt (fun τ : ℝ => uP b ((1 + b) * τ)) ((1 + b) ^ 2) 0 := by
    have hn : HasDerivAt (fun τ : ℝ => 1 + (1 + b) * τ) (1 + b) 0 := by
      simpa using ((hasDerivAt_id (0 : ℝ)).const_mul (1 + b)).const_add 1
    have hd : HasDerivAt (fun τ : ℝ => 1 - b * ((1 + b) * τ)) (-(b * (1 + b))) 0 := by
      simpa using (((hasDerivAt_id (0 : ℝ)).const_mul (1 + b)).const_mul b).const_sub 1
    have := hn.div hd (by simp)
    unfold uP
    convert this using 1
    simp; ring
  have hM' : HasDerivAt (fun τ : ℝ => uM b ((1 + b) * τ)) (-(1 + b) ^ 2) 0 := by
    have hn : HasDerivAt (fun τ : ℝ => 1 - (1 + b) * τ) (-(1 + b)) 0 := by
      simpa using ((hasDerivAt_id (0 : ℝ)).const_mul (1 + b)).const_sub 1
    have hd : HasDerivAt (fun τ : ℝ => 1 + b * ((1 + b) * τ)) (b * (1 + b)) 0 := by
      simpa using (((hasDerivAt_id (0 : ℝ)).const_mul (1 + b)).const_mul b).const_add 1
    have := hn.div hd (by simp)
    unfold uM
    convert this using 1
    simp; ring
  have h1 : uP b ((1 + b) * 0) = 1 := by simp [uP]
  have h2 : uM b ((1 + b) * 0) = 1 := by simp [uM]
  have hsP := hP.sqrt (by rw [h1]; norm_num)
  have hsM := hM'.sqrt (by rw [h2]; norm_num)
  have := (hsP.sub hsM).const_mul (√(2 * M) / 2)
  have e : (fun τ => epsK M (curve M b ((1 + b) * τ)))
      = fun τ => √(2 * M) / 2 * (√(uP b ((1 + b) * τ)) - √(uM b ((1 + b) * τ))) := by
    funext τ; unfold curve; rw [(epsK_stag M hM _ _).1]
  rw [e]
  convert this using 1
  rw [h1, h2, Real.sqrt_one]; ring

/-- **Even `n`: `c₄` from the Hamiltonian.** If `h₂ = 0` (the threshold; `h2_eq_zero_iff`), then
along the symmetric curve `(H − H(ring))/ε⁴ → c4evenOf M λ = h₄(2/(1+λ)²)⁴/n²` as `τ → 0`. -/
theorem c4even_limit (M : ℕ) (hM : M ≠ 0) (b : ℝ) (hb : |b| < 1) (hthr : h2 M b = 0) :
    Tendsto (fun τ => (Ham b (curve M b ((1 + b) * τ)) - Ham b (stag M 1 1))
        / epsK M (curve M b ((1 + b) * τ)) ^ 4) (𝓝[≠] 0) (𝓝 (c4evenOf M b)) := by
  have hb' := abs_lt.mp hb
  have hb1 : 0 < 1 + b := by linarith
  set κ := √(2 * M) / 2 * (1 + b) ^ 2 with hκ
  have hκ0 : κ ≠ 0 := by
    have : 0 < √(2 * (M : ℝ)) := Real.sqrt_pos.mpr (by positivity)
    positivity
  -- `ε(τ)/τ → κ`
  have hslope := (epsK_curve_hasDerivAt M hM b).tendsto_slope_zero
  have he0 : epsK M (curve M b ((1 + b) * 0)) = 0 := by
    unfold curve; rw [(epsK_stag M hM _ _).1]; simp [uP, uM]
  simp only [zero_add, he0, sub_zero] at hslope
  -- `Φ(τ)/τ⁴ → h₄`
  set Φ := fun τ => Ham b (curve M b ((1 + b) * τ)) - Ham b (stag M 1 1) with hΦ
  have hjet := Ham_curve_jet M hM b hb
  rw [hthr] at hjet
  have hq : Tendsto (fun τ => (Φ τ - h4 M b * τ ^ 4) * (τ ^ 4)⁻¹) (𝓝[≠] 0) (𝓝 0) := by
    have h1 := (hjet.mono nhdsWithin_le_nhds).mul (isBigO_refl (fun τ : ℝ => (τ ^ 4)⁻¹) (𝓝[≠] 0))
    have h2 : Tendsto (fun τ : ℝ => τ ^ 6 * (τ ^ 4)⁻¹) (𝓝[≠] 0) (𝓝 0) := by
      have : Tendsto (fun τ : ℝ => τ ^ 2) (𝓝[≠] 0) (𝓝 0) := by
        have hc : Continuous (fun τ : ℝ => τ ^ 2) := by fun_prop
        simpa using (hc.tendsto 0).mono_left nhdsWithin_le_nhds
      refine this.congr' ?_
      filter_upwards [self_mem_nhdsWithin] with τ hτ
      field_simp
    refine (h1.trans_tendsto h2).congr' ?_
    filter_upwards with τ
    simp only [zero_mul, sub_zero]
    rfl
  have hΦ4 : Tendsto (fun τ => Φ τ / τ ^ 4) (𝓝[≠] 0) (𝓝 (h4 M b)) := by
    have := hq.add_const (h4 M b)
    rw [zero_add] at this
    refine this.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with τ hτ
    have : τ ^ 4 ≠ 0 := pow_ne_zero 4 hτ
    rw [sub_mul, mul_inv_cancel_right₀ this, div_eq_mul_inv]
    ring
  have hlim := hΦ4.div (hslope.pow 4) (pow_ne_zero 4 hκ0)
  have hval : h4 M b / κ ^ 4 = c4evenOf M b := by
    unfold c4evenOf
    rw [hκ]
    have hsq : √(2 * (M : ℝ)) ^ 2 = 2 * M := Real.sq_sqrt (by positivity)
    have hM0 : (M : ℝ) ≠ 0 := by exact_mod_cast hM
    have : (√(2 * (M : ℝ)) / 2 * (1 + b) ^ 2) ^ 4 = (2 * M) ^ 2 / 16 * (1 + b) ^ 8 := by
      rw [show (√(2 * (M : ℝ)) / 2 * (1 + b) ^ 2) ^ 4 = (√(2 * (M : ℝ)) ^ 2) ^ 2 / 16 * (1 + b) ^ 8 by ring,
        hsq]
    rw [this]
    field_simp
    ring
  rw [← hval]
  refine hlim.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with τ hτ
  have : τ ≠ 0 := hτ
  simp only [Pi.div_apply, smul_eq_mul, hΦ]
  rw [mul_pow, inv_pow, ← div_eq_inv_mul, div_div_div_cancel_right₀ (pow_ne_zero 4 this)]

/-- **Even `n = 2M ≥ 4` at `b = b_n`**: along the symmetric curve in the momentum level set,
`(H − H(ring))/ε⁴ → (n−2)(m⁴+40m²+48m+16)/(192n(n−1))`, `m = n − 2` (paper Thm even), where `H` is
paper eq. (1) and `ε` the component along the unit kernel vector. -/
theorem c4even_limit_at_bn (M : ℕ) (hM : 2 ≤ M) :
    Tendsto (fun τ => (Ham (bn (2 * M)) (curve M (bn (2 * M)) ((1 + bn (2 * M)) * τ))
          - Ham (bn (2 * M)) (stag M 1 1))
        / epsK M (curve M (bn (2 * M)) ((1 + bn (2 * M)) * τ)) ^ 4) (𝓝[≠] 0)
      (𝓝 (c4even ((2 * M : ℕ) : ℝ))) := by
  obtain ⟨h1, h2, -⟩ := bn_spec (2 * M) (by omega)
  obtain ⟨hz, hc⟩ := c4_even_at_bn M hM
  rw [← hc]
  exact c4even_limit M (by omega) _ (abs_lt.mpr ⟨h1, h2⟩) hz

end Havelock.Conj13.EvenQuartic
