/-
Havelock — Conjecture 13, odd `n = 2M+1`: the kernel coefficient `f₂` derived from the Hamiltonian
(paper `conj13/paper/main.tex`, proof of Thm odd, "Result": the `A²` coefficient of `H` along the kernel).

Canonical coordinates `w = |z|²/σ`, `ϑ = arg z`; the kernel line is `w_k = w₀ + A cos(lφ_k)`,
`ϑ_k = φ_k = 2πk/n`, `w₀ = 1/(1+λ)`, `l = M = (n−1)/2` (`oddLine`). Main results:

* `pairTerm_canonical`  `|z_i−z_j|²/(σ_iσ_j) = w_i + w_j − 2λw_iw_j − 2√(P(w_i)P(w_j)) cos(ϑ_i−ϑ_j)`,
                        `P(w) = w(1−λw)` (paper, proof of Thm odd).
* `oddLine_zero`        `A = 0` is the unit ring; `J3_oddLine`: `J₃` is constant along the line.
* `Ham_oddLine_jet`     **`H(oddLine A) − H(ring) = f₂ A² + O(A³)`**, `f₂ = f2 n λ` the definition of
                        `Conj13Algebra.lean`, for every odd `n ≥ 3` and `λ ∈ (−1,1)` (no threshold assumed).
-/
import Havelock.Conj13.OddSums
import Havelock.Conj13.EvenQuartic

open Real Asymptotics Filter Topology Finset

namespace Havelock.Conj13.Odd

open Havelock.Conj13.Series Havelock.Conj13.Config Havelock.Conj13.OddSums

/-! ### Second-order jets -/

/-- `log D(A) = log d₀ + (d₁/d₀)A + (d₂/d₀ − d₁²/(2d₀²))A² + O(A³)` if `D(A) = d₀ + d₁A + d₂A² + O(A³)`,
`d₀ > 0`. -/
theorem log_jet2 {D : ℝ → ℝ} {d0 d1 d2 : ℝ} (hd0 : 0 < d0)
    (hD : (fun A => D A - (d0 + d1 * A + d2 * A ^ 2)) =O[𝓝 0] (fun A => A ^ 3)) :
    (fun A => Real.log (D A) - Real.log d0 - d1 / d0 * A - (d2 / d0 - d1 ^ 2 / (2 * d0 ^ 2)) * A ^ 2)
      =O[𝓝 0] (fun A => A ^ 3) := by
  set R := fun A => D A - (d0 + d1 * A + d2 * A ^ 2) with hR
  set t := fun A => 1 - D A / d0 with ht_def
  have hpoly1 : (fun A : ℝ => A ^ 1 * (d1 + d2 * A)) =O[𝓝 0] (fun A => A ^ 1) := mul_cont (by fun_prop)
  have ht : t =O[𝓝 0] (fun A => A ^ 1) := by
    refine ((hpoly1.add (bigO_weaken (by norm_num) hD)).const_mul_left (-1 / d0)).congr_left
      (fun A => ?_)
    simp only [ht_def, hR]; field_simp; ring
  have hL := comp_bigO (by norm_num) log_one_sub_jet ht
  have hplus : (fun A => t A + d1 / d0 * A) =O[𝓝 0] (fun A => A ^ 2) := by
    have h1 : (fun A : ℝ => A ^ 2 * d2) =O[𝓝 0] (fun A => A ^ 2) := mul_cont continuousAt_const
    refine ((h1.add (bigO_weaken (by norm_num) hD)).const_mul_left (-1 / d0)).congr_left
      (fun A => ?_)
    simp only [ht_def, hR]; field_simp; ring
  have hminus : (fun A => t A - d1 / d0 * A) =O[𝓝 0] (fun A => A ^ 1) := by
    have h1 : (fun A : ℝ => A ^ 1 * (d1 / d0)) =O[𝓝 0] (fun A => A ^ 1) := mul_cont continuousAt_const
    exact (ht.sub h1).congr_left (fun A => by ring)
  have hprod := mul_pow_bigO hplus hminus
  -- `D > 0` near `0`
  have hDt : Tendsto D (𝓝 0) (𝓝 d0) := by
    have h1 : Tendsto (fun A : ℝ => d0 + d1 * A + d2 * A ^ 2) (𝓝 0) (𝓝 d0) := by
      have hc : Continuous (fun A : ℝ => d0 + d1 * A + d2 * A ^ 2) := by fun_prop
      simpa using hc.tendsto 0
    have := h1.add (tendsto_of_bigO (by norm_num) hD)
    simp only [add_zero] at this
    refine this.congr (fun A => ?_)
    ring
  have hpos : ∀ᶠ A in 𝓝 0, 0 < D A := hDt.eventually (lt_mem_nhds hd0)
  have := ((hL.add ((bigO_weaken (by norm_num) hD).const_mul_left (1 / d0))).sub
    (hprod.const_mul_left (1 / 2)))
  refine this.congr' ?_ (Eventually.of_forall fun _ => rfl)
  filter_upwards [hpos] with A hA
  have hlog : Real.log (D A) = Real.log d0 + Real.log (1 - t A) := by
    simp only [ht_def, sub_sub_cancel]
    rw [← Real.log_mul hd0.ne' (div_pos hA hd0).ne', mul_div_cancel₀ _ hd0.ne']
  rw [hlog]
  simp only [ht_def, hR]
  field_simp
  ring

/-- `√g(A) = s + p₁A + p₂A² + O(A³)` if `g(A) = g₀ + g₁A + g₂A² + O(A³)`, `g₀ > 0`, with `s = √g₀`,
`p₁ = g₁/(2s)`, `p₂ = (g₂ − p₁²)/(2s)`. -/
theorem sqrt_jet2 {g : ℝ → ℝ} {g0 g1 g2 : ℝ} (hg0 : 0 < g0)
    (hg : (fun A => g A - (g0 + g1 * A + g2 * A ^ 2)) =O[𝓝 0] (fun A => A ^ 3)) :
    (fun A => √(g A) - (√g0 + g1 / (2 * √g0) * A
      + (g2 - (g1 / (2 * √g0)) ^ 2) / (2 * √g0) * A ^ 2)) =O[𝓝 0] (fun A => A ^ 3) := by
  set s := √g0 with hs_def
  have hs : 0 < s := Real.sqrt_pos.mpr hg0
  have hs2 : s ^ 2 = g0 := Real.sq_sqrt hg0.le
  set p1 := g1 / (2 * s) with hp1
  set p2 := (g2 - p1 ^ 2) / (2 * s) with hp2
  set p := fun A : ℝ => s + p1 * A + p2 * A ^ 2 with hp
  have hgt : Tendsto g (𝓝 0) (𝓝 g0) := by
    have h1 : Tendsto (fun A : ℝ => g0 + g1 * A + g2 * A ^ 2) (𝓝 0) (𝓝 g0) := by
      have hc : Continuous (fun A : ℝ => g0 + g1 * A + g2 * A ^ 2) := by fun_prop
      simpa using hc.tendsto 0
    have := h1.add (tendsto_of_bigO (by norm_num) hg)
    simp only [add_zero] at this
    exact this.congr (fun A => by ring)
  have hden : Tendsto (fun A => √(g A) + p A) (𝓝 0) (𝓝 (2 * s)) := by
    have h1 := (Real.continuous_sqrt.tendsto g0).comp hgt
    have h2 : Tendsto p (𝓝 0) (𝓝 s) := by
      have hc : Continuous p := by simp only [hp]; fun_prop
      simpa [hp] using hc.tendsto 0
    have := h1.add h2
    rw [show √g0 + s = 2 * s by rw [← hs_def]; ring] at this
    exact this
  have hinv : (fun A => (√(g A) + p A)⁻¹) =O[𝓝 0] (fun _ => (1 : ℝ)) :=
    (hden.inv₀ (by positivity)).isBigO_one ℝ
  -- `g − p² = O(A³)`
  have hgp : (fun A => g A - p A ^ 2) =O[𝓝 0] (fun A => A ^ 3) := by
    have h1 : (fun A : ℝ => A ^ 3 * (-(2 * p1 * p2 + p2 ^ 2 * A))) =O[𝓝 0] (fun A => A ^ 3) :=
      mul_cont (by fun_prop)
    refine (hg.add h1).congr_left (fun A => ?_)
    simp only [hp]
    have e1 : g1 = 2 * s * p1 := by rw [hp1]; field_simp
    have e2 : g2 = 2 * s * p2 + p1 ^ 2 := by rw [hp2]; field_simp; ring
    rw [e1, e2, ← hs2]
    ring
  have := hgp.mul hinv
  simp only [mul_one] at this
  refine this.congr' ?_ (Eventually.of_forall fun _ => rfl)
  have hgpos : ∀ᶠ A in 𝓝 0, 0 < g A := hgt.eventually (lt_mem_nhds hg0)
  have hdpos : ∀ᶠ A in 𝓝 0, 0 < √(g A) + p A :=
    hden.eventually (lt_mem_nhds (by positivity))
  filter_upwards [hgpos, hdpos] with A hA hd
  have hsq : √(g A) ^ 2 = g A := Real.sq_sqrt hA.le
  rw [mul_inv_eq_iff_eq_mul₀ hd.ne']
  simp only [hp]
  linear_combination (-1 : ℝ) * hsq

/-! ### One pair of vortices -/

/-- `w₀ = 1/(1+λ)` (the ring in canonical coordinates). -/
noncomputable def w0 (b : ℝ) : ℝ := 1 / (1 + b)

/-- `P(w) = w(1 − λw)`. -/
noncomputable def Pw (b w : ℝ) : ℝ := w * (1 - b * w)

/-- The pair quantity `w_i + w_j − 2λw_iw_j − 2√(P(w_i)P(w_j)) C` along the kernel line,
`w_k = w₀ + A c_k`. -/
noncomputable def Dpair (b ci cj C A : ℝ) : ℝ :=
  (w0 b + A * ci) + (w0 b + A * cj) - 2 * b * (w0 b + A * ci) * (w0 b + A * cj)
    - 2 * √(Pw b (w0 b + A * ci) * Pw b (w0 b + A * cj)) * C

/-- The `A` coefficient of `log D_ij`. -/
noncomputable def alphaP (b ci cj : ℝ) : ℝ := (1 - b) * (1 + b) * (ci + cj) / 2

/-- The `A²` coefficient of `log D_ij`. -/
noncomputable def betaP (b ci cj C : ℝ) : ℝ :=
  -(1 + b) ^ 4 / 4 * (ci * cj / (1 - C)) + (1 + b) ^ 4 / 8 * ((ci ^ 2 + cj ^ 2) / (1 - C))
    - (1 + b) ^ 2 * (1 + b ^ 2) / 4 * (ci ^ 2 + cj ^ 2)

theorem Dpair_zero (b ci cj C : ℝ) (hb : -1 < b) :
    Dpair b ci cj C 0 = 2 / (1 + b) ^ 2 * (1 - C) := by
  have h1 : 0 < 1 + b := by linarith
  have hP : Pw b (w0 b) = (1 / (1 + b)) ^ 2 := by unfold Pw w0; field_simp; ring
  unfold Dpair
  simp only [zero_mul, add_zero]
  rw [hP, show (1 / (1 + b)) ^ 2 * (1 / (1 + b)) ^ 2 = ((1 / (1 + b)) ^ 2) ^ 2 by ring,
    Real.sqrt_sq (by positivity)]
  unfold w0; field_simp; ring

/-- **One pair**: `log D_ij(A) = log D_ij(0) + α_ij A + β_ij A² + O(A³)`. -/
theorem Dpair_jet (b ci cj C : ℝ) (hb : -1 < b) (hC : C < 1) :
    (fun A => Real.log (Dpair b ci cj C A) - Real.log (Dpair b ci cj C 0)
      - alphaP b ci cj * A - betaP b ci cj C * A ^ 2) =O[𝓝 0] (fun A => A ^ 3) := by
  have h1 : 0 < 1 + b := by linarith
  have hC' : 0 < 1 - C := by linarith
  set P0 := (1 / (1 + b)) ^ 2 with hP0
  set P1 := (1 - b) / (1 + b) with hP1
  have hP0pos : 0 < P0 := by positivity
  have hPw : ∀ h : ℝ, Pw b (w0 b + h) = P0 + P1 * h - b * h ^ 2 := by
    intro h; unfold Pw w0; rw [hP0, hP1]; field_simp; ring
  set g := fun A => Pw b (w0 b + A * ci) * Pw b (w0 b + A * cj) with hg_def
  have hg : (fun A => g A - (P0 ^ 2 + P0 * P1 * (ci + cj) * A
      + (P1 ^ 2 * ci * cj - b * P0 * (ci ^ 2 + cj ^ 2)) * A ^ 2)) =O[𝓝 0] (fun A => A ^ 3) := by
    have : (fun A : ℝ => A ^ 3 * (-(b * P1 * ci * cj * (ci + cj)) + b ^ 2 * ci ^ 2 * cj ^ 2 * A))
        =O[𝓝 0] (fun A => A ^ 3) := mul_cont (by fun_prop)
    refine this.congr_left (fun A => ?_)
    simp only [hg_def, hPw]
    ring
  have hsq := sqrt_jet2 (by positivity) hg
  have hs : √(P0 ^ 2) = P0 := Real.sqrt_sq hP0pos.le
  rw [hs] at hsq
  set p1 := P0 * P1 * (ci + cj) / (2 * P0) with hp1
  set p2 := (P1 ^ 2 * ci * cj - b * P0 * (ci ^ 2 + cj ^ 2) - p1 ^ 2) / (2 * P0) with hp2
  -- `D = a(A) − 2C√g` with `a` exactly quadratic
  set d0 := 2 * w0 b - 2 * b * w0 b ^ 2 - 2 * C * P0 with hd0
  set d1 := (ci + cj) * (1 - 2 * b * w0 b) - 2 * C * p1 with hd1
  set d2 := -2 * b * ci * cj - 2 * C * p2 with hd2
  have hD : (fun A => Dpair b ci cj C A - (d0 + d1 * A + d2 * A ^ 2)) =O[𝓝 0] (fun A => A ^ 3) := by
    refine (hsq.const_mul_left (-2 * C)).congr_left (fun A => ?_)
    simp only [Dpair, hg_def, hd0, hd1, hd2]
    ring
  have hd0v : d0 = 2 * P0 * (1 - C) := by
    rw [hd0, hP0]; unfold w0; field_simp; ring
  have hd0pos : 0 < d0 := by rw [hd0v]; positivity
  have hL := log_jet2 hd0pos hD
  have hz : Dpair b ci cj C 0 = d0 := by
    rw [Dpair_zero b ci cj C hb, hd0v, hP0]; field_simp
  refine hL.congr_left (fun A => ?_)
  rw [hz]
  have ha : d1 / d0 = alphaP b ci cj := by
    rw [hd0v, hd1, hp1, hP0, hP1]; unfold alphaP w0; field_simp; ring
  have hbeta : d2 / d0 - d1 ^ 2 / (2 * d0 ^ 2) = betaP b ci cj C := by
    rw [hd0v, hd1, hd2, hp2, hp1, hP0, hP1]; unfold betaP w0; field_simp; ring
  rw [ha, hbeta]

/-! ### The configuration in canonical coordinates -/

/-- `u = |z|² = w/(1 − λw)` from `w = |z|²/σ`. -/
noncomputable def uOf (b w : ℝ) : ℝ := w / (1 - b * w)

/-- The kernel line in canonical coordinates: `w_k = w₀ + A cos(lφ_k)`, `ϑ_k = φ_k = 2πk/n`,
`n = 2M+1`, `l = M`; i.e. `z_k = √(u(w_k)) e^{iφ_k}`. -/
noncomputable def oddLine (M : ℕ) (b A : ℝ) (k : Fin (2 * M + 1)) : ℂ :=
  (√(uOf b (w0 b + A * cc M ((k : ℕ) : ℝ))) : ℂ)
    * Complex.exp (↑(2 * π * ((k : ℕ) : ℝ) / (2 * M + 1 : ℕ)) * Complex.I)

theorem normSq_polar_sub (r r' θ θ' : ℝ) :
    Complex.normSq ((r : ℂ) * Complex.exp (↑θ * Complex.I) - (r' : ℂ) * Complex.exp (↑θ' * Complex.I))
      = r ^ 2 + r' ^ 2 - 2 * r * r' * Real.cos (θ - θ') := by
  rw [Complex.normSq_apply]
  simp only [Complex.sub_re, Complex.sub_im, Complex.mul_re, Complex.mul_im, Complex.ofReal_re,
    Complex.ofReal_im, Complex.exp_ofReal_mul_I_re, Complex.exp_ofReal_mul_I_im, zero_mul, sub_zero,
    add_zero]
  rw [Real.cos_sub]
  linear_combination r ^ 2 * Real.sin_sq_add_cos_sq θ + r' ^ 2 * Real.sin_sq_add_cos_sq θ'

theorem normSq_polar (r θ : ℝ) :
    Complex.normSq ((r : ℂ) * Complex.exp (↑θ * Complex.I)) = r ^ 2 := by
  rw [Complex.normSq_mul, Complex.normSq_ofReal, Complex.normSq_eq_norm_sq,
    Complex.norm_exp_ofReal_mul_I]
  ring

theorem uOf_rel (b w : ℝ) (h1 : 0 < 1 - b * w) :
    (1 + b * uOf b w) * (1 - b * w) = 1 ∧ uOf b w * (1 - b * w) = w := by
  unfold uOf
  constructor
  · rw [add_mul, one_mul, mul_assoc, div_mul_cancel₀ _ h1.ne']; ring
  · exact div_mul_cancel₀ _ h1.ne'

/-- **Canonical coordinates** (paper, proof of Thm odd): with `u = w/(1−λw)`, `z = √u e^{iϑ}`,
`|z−z'|²/(σσ') = w + w' − 2λww' − 2√(P(w)P(w')) cos(ϑ−ϑ')`. -/
theorem pairTerm_canonical (b w w' θ θ' : ℝ) (hw : 0 < w) (hw' : 0 < w') (h1 : 0 < 1 - b * w)
    (h1' : 0 < 1 - b * w') :
    pairTerm b ((√(uOf b w) : ℂ) * Complex.exp (↑θ * Complex.I))
        ((√(uOf b w') : ℂ) * Complex.exp (↑θ' * Complex.I))
      = Real.log (w + w' - 2 * b * w * w' - 2 * √(Pw b w * Pw b w') * Real.cos (θ - θ')) := by
  have hu : 0 ≤ uOf b w := (div_pos hw h1).le
  have hu' : 0 ≤ uOf b w' := (div_pos hw' h1').le
  have hr : √(uOf b w) * (1 - b * w) = √(Pw b w) := by
    rw [← Real.sqrt_sq h1.le, ← Real.sqrt_mul hu]
    congr 1; unfold uOf Pw; field_simp
  have hr' : √(uOf b w') * (1 - b * w') = √(Pw b w') := by
    rw [← Real.sqrt_sq h1'.le, ← Real.sqrt_mul hu']
    congr 1; unfold uOf Pw; field_simp
  have hP : 0 ≤ Pw b w := by unfold Pw; positivity
  obtain ⟨E1, E2⟩ := uOf_rel b w h1
  obtain ⟨E1', E2'⟩ := uOf_rel b w' h1'
  have hs : 1 + b * uOf b w ≠ 0 := fun h => by rw [h, zero_mul] at E1; exact zero_ne_one E1
  have hs' : 1 + b * uOf b w' ≠ 0 := fun h => by rw [h, zero_mul] at E1'; exact zero_ne_one E1'
  unfold pairTerm sig
  rw [normSq_polar_sub, normSq_polar, normSq_polar, Real.sq_sqrt hu, Real.sq_sqrt hu',
    Real.sqrt_mul hP, ← hr, ← hr']
  congr 1
  rw [div_eq_iff (mul_ne_zero hs hs')]
  set u := uOf b w
  set u' := uOf b w'
  set r := √u
  set r' := √u'
  set C := Real.cos (θ - θ')
  linear_combination (-(u + w' * (1 + b * u') - 2 * C * r * r' * ((1 + b * u') * (1 - b * w')))) * E1
    + (-(w * (1 + b * u) + u' - 2 * C * r * r')) * E1' + (1 + b * u) * E2 + (1 + b * u') * E2'

theorem cc_abs_le (M : ℕ) (x : ℝ) : |cc M x| ≤ 1 := Real.abs_cos_le_one _

/-- On `|A| < w₀` all vortices stay in the chart: `w_k > 0`, `1 − λw_k > 0`. -/
theorem w_pos (M : ℕ) (b A x : ℝ) (hb : -1 < b) (hb1 : b < 1) (hA : |A| < w0 b) :
    0 < w0 b + A * cc M x ∧ 0 < 1 - b * (w0 b + A * cc M x) := by
  have h1 : 0 < 1 + b := by linarith
  have hc := cc_abs_le M x
  have hAc : |A * cc M x| ≤ |A| := by
    rw [abs_mul]; calc |A| * |cc M x| ≤ |A| * 1 := by gcongr
      _ = |A| := mul_one _
  have hb' : |b| ≤ 1 := abs_le.mpr ⟨by linarith, by linarith⟩
  have hbAc : |b * (A * cc M x)| ≤ |A| := by
    rw [abs_mul]; nlinarith [abs_nonneg b, abs_nonneg (A * cc M x)]
  have h1' : 1 + b ≠ 0 := h1.ne'
  have e : 1 - b * (w0 b + A * cc M x) = w0 b - b * (A * cc M x) := by
    unfold w0; rw [mul_add, show b * (1 / (1 + b)) = 1 - 1 / (1 + b) by field_simp; ring]; ring
  rw [e]
  constructor
  · have := neg_abs_le (A * cc M x); linarith
  · have := le_abs_self (b * (A * cc M x)); linarith

theorem CC_lt_one (M i j : ℕ) (hi : i < 2 * M + 1) (hj : j < 2 * M + 1) (hij : i ≠ j) :
    CC M ((j : ℝ) - i) < 1 := by
  refine lt_of_le_of_ne (Real.cos_le_one _) (fun h => hij ?_)
  obtain ⟨k, hk⟩ := (Real.cos_eq_one_iff _).mp h
  have hn : ((2 * M + 1 : ℕ) : ℝ) ≠ 0 := by positivity
  have hk' : (k : ℝ) * (2 * M + 1 : ℕ) = (j : ℝ) - i := by
    field_simp at hk
    nlinarith [Real.pi_pos]
  have hz : (k * (2 * M + 1 : ℕ) : ℤ) = (j : ℤ) - i := by exact_mod_cast hk'
  have hdvd : ((2 * M + 1 : ℕ) : ℤ) ∣ (j : ℤ) - i := ⟨k, by rw [← hz]; ring⟩
  have : (j : ℤ) - i = 0 := by
    apply Int.eq_zero_of_abs_lt_dvd hdvd
    rw [abs_lt]; constructor <;> push_cast <;> omega
  omega

theorem sum_fin_offdiag (n : ℕ) (F : ℕ → ℕ → ℝ) :
    (∑ i : Fin n, ∑ j : Fin n, if i = j then 0 else F i j)
      = ∑ i ∈ range n, ∑ j ∈ range n, if i = j then 0 else F i j := by
  rw [← Fin.sum_univ_eq_sum_range (fun i => ∑ j ∈ range n, if i = j then 0 else F i j)]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [← Fin.sum_univ_eq_sum_range (fun j => if (i : ℕ) = j then 0 else F i j)]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  simp [Fin.ext_iff]

/-- `H` along the kernel line, as a sum over pairs of `log D_ij` (for `|A| < w₀`). -/
theorem Ham_oddLine (M : ℕ) (b A : ℝ) (hb : -1 < b) (hb1 : b < 1) (hA : |A| < w0 b) :
    Ham b (oddLine M b A) = -(1 / 2) * ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else Real.log (Dpair b (cc M i) (cc M j) (CC M ((j : ℝ) - i)) A)) := by
  unfold Ham
  congr 1
  rw [← sum_fin_offdiag]
  refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => ?_))
  split_ifs with h
  · rfl
  · obtain ⟨hwi, hbi⟩ := w_pos M b A ((i : ℕ) : ℝ) hb hb1 hA
    obtain ⟨hwj, hbj⟩ := w_pos M b A ((j : ℕ) : ℝ) hb hb1 hA
    unfold oddLine
    rw [pairTerm_canonical b _ _ _ _ hwi hwj hbi hbj]
    unfold Dpair
    congr 2
    unfold CC
    rw [← Real.cos_neg]
    congr 1
    ring_nf

/-- `A = 0` is the unit ring `z_k = e^{2πik/n}`. -/
theorem oddLine_zero (M : ℕ) (b : ℝ) (hb : -1 < b) : oddLine M b 0 = ringConf (2 * M + 1) := by
  have h1 : 0 < 1 + b := by linarith
  have hu : uOf b (w0 b) = 1 := by unfold uOf w0; field_simp; ring
  funext k
  unfold oddLine ringConf
  rw [zero_mul, add_zero, hu, Real.sqrt_one, Complex.ofReal_one, one_mul]
  congr 1
  push_cast
  ring

/-- `J₃ = Σ w_k` is constant along the kernel line (it lies in `J₃⁻¹(μ₃)`). -/
theorem J3_oddLine (M : ℕ) (hM : 1 ≤ M) (b A : ℝ) (hb : -1 < b) (hb1 : b < 1) (hA : |A| < w0 b) :
    EvenQuartic.J3 b (oddLine M b A) = EvenQuartic.J3 b (oddLine M b 0) := by
  have hterm : ∀ A : ℝ, |A| < w0 b → ∀ k : Fin (2 * M + 1),
      Complex.normSq (oddLine M b A k) / sig b (oddLine M b A k) = w0 b + A * cc M ((k : ℕ) : ℝ) := by
    intro A hA k
    obtain ⟨hw, h1⟩ := w_pos M b A ((k : ℕ) : ℝ) hb hb1 hA
    have hu : 0 ≤ uOf b (w0 b + A * cc M ((k : ℕ) : ℝ)) := (div_pos hw h1).le
    obtain ⟨E1, E2⟩ := uOf_rel b _ h1
    have hs : 1 + b * uOf b (w0 b + A * cc M ((k : ℕ) : ℝ)) ≠ 0 :=
      fun h => by rw [h, zero_mul] at E1; exact zero_ne_one E1
    unfold oddLine sig
    rw [normSq_polar, Real.sq_sqrt hu, div_eq_iff hs]
    linear_combination (-(uOf b (w0 b + A * cc M ((k : ℕ) : ℝ)))) * E1
      + (1 + b * uOf b (w0 b + A * cc M ((k : ℕ) : ℝ))) * E2
  have h0 : |(0 : ℝ)| < w0 b := by
    rw [abs_zero]; unfold w0; exact div_pos one_pos (by linarith)
  unfold EvenQuartic.J3
  rw [Finset.sum_congr rfl (fun k _ => hterm A hA k), Finset.sum_congr rfl (fun k _ => hterm 0 h0 k),
    Finset.sum_add_distrib, Finset.sum_add_distrib,
    ← Finset.mul_sum, ← Finset.mul_sum,
    Fin.sum_univ_eq_sum_range (fun k => cc M (k : ℝ)) (2 * M + 1), sum_cc M hM]
  ring

/-! ### The `A²` coefficient: `f₂` -/

theorem sum_alpha (M : ℕ) (hM : 1 ≤ M) (b : ℝ) :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else alphaP b (cc M i) (cc M j)) = 0 := by
  have h : ∀ i j : ℕ, (if i = j then (0 : ℝ) else alphaP b (cc M i) (cc M j))
      = (1 - b) * (1 + b) / 2 * (if i = j then 0 else cc M i + cc M j) := by
    intro i j; unfold alphaP; split_ifs <;> ring
  simp_rw [h, ← Finset.mul_sum, pairB0 M hM, mul_zero]

theorem sum_beta (M : ℕ) (hM : 1 ≤ M) (b : ℝ) :
    ∑ i ∈ range (2 * M + 1), ∑ j ∈ range (2 * M + 1),
      (if i = j then 0 else betaP b (cc M i) (cc M j) (CC M ((j : ℝ) - i)))
        = -2 * f2 ((2 * M + 1 : ℕ) : ℝ) b := by
  have h : ∀ i j : ℕ, (if i = j then (0 : ℝ) else betaP b (cc M i) (cc M j) (CC M ((j : ℝ) - i)))
      = -(1 + b) ^ 4 / 4 * (if i = j then 0 else cc M i * cc M j / (1 - CC M ((j : ℝ) - i)))
        + (1 + b) ^ 4 / 8 * (if i = j then 0 else (cc M i ^ 2 + cc M j ^ 2) / (1 - CC M ((j : ℝ) - i)))
        - (1 + b) ^ 2 * (1 + b ^ 2) / 4 * (if i = j then 0 else cc M i ^ 2 + cc M j ^ 2) := by
    intro i j; unfold betaP; split_ifs <;> ring
  simp only [h, Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]
  rw [pairA1 M hM, pairA2 M hM, pairB2 M hM]
  unfold f2
  ring

/-- **Odd `n = 2M+1 ≥ 3`, the kernel coefficient from the Hamiltonian.** Along the kernel line in
canonical coordinates (`w_k = w₀ + A cos(lφ_k)`, `ϑ_k = φ_k`; it starts at the ring and lies in
`J₃⁻¹(μ₃)`, `J3_oddLine`), `H(oddLine A) − H(ring) = f₂ A² + O(A³)`, with `f₂ = f2 n λ` exactly the
expression of the paper's proof of Thm odd (definition in `Conj13Algebra.lean`), for every
`λ ∈ (−1, 1)`. (In particular the `A¹` coefficient vanishes, and `f₂ = 0` iff MT (3.6),
`f2_eq_zero_iff`.) -/
theorem Ham_oddLine_jet (M : ℕ) (hM : 1 ≤ M) (b : ℝ) (hb : -1 < b) (hb1 : b < 1) :
    (fun A => Ham b (oddLine M b A) - Ham b (ringConf (2 * M + 1))
        - f2 ((2 * M + 1 : ℕ) : ℝ) b * A ^ 2) =O[𝓝 0] (fun A => A ^ 3) := by
  set n := 2 * M + 1 with hn
  set D := fun (i j : ℕ) (A : ℝ) => Dpair b (cc M i) (cc M j) (CC M ((j : ℝ) - i)) A with hD
  set E := fun A => ∑ i ∈ range n, ∑ j ∈ range n, (if i = j then 0 else
    Real.log (D i j A) - Real.log (D i j 0) - alphaP b (cc M i) (cc M j) * A
      - betaP b (cc M i) (cc M j) (CC M ((j : ℝ) - i)) * A ^ 2) with hE_def
  have hE : E =O[𝓝 0] (fun A => A ^ 3) := by
    rw [hE_def]
    refine IsBigO.fun_sum (fun i hi => ?_)
    refine IsBigO.fun_sum (fun j hj => ?_)
    by_cases hij : i = j
    · simp only [hij, ite_true]; exact isBigO_zero _ _
    · simp only [hij, ite_false]
      exact Dpair_jet b _ _ _ hb (CC_lt_one M i j (Finset.mem_range.mp hi) (Finset.mem_range.mp hj) hij)
  have hw0 : 0 < w0 b := by unfold w0; exact div_pos one_pos (by linarith)
  have hev : ∀ᶠ A in 𝓝 (0 : ℝ), |A| < w0 b := by
    have hc : Continuous (fun A : ℝ => |A|) := continuous_abs
    exact hc.continuousAt.eventually_lt continuousAt_const (by simpa using hw0)
  have := hE.const_mul_left (-(1 / 2))
  refine this.congr' ?_ (Eventually.of_forall fun _ => rfl)
  filter_upwards [hev] with A hA
  have h0 : |(0 : ℝ)| < w0 b := by rw [abs_zero]; exact hw0
  rw [← oddLine_zero M b hb, Ham_oddLine M b A hb hb1 hA, Ham_oddLine M b 0 hb hb1 h0]
  have hsa := sum_alpha M hM b
  have hsb := sum_beta M hM b
  have key : (∑ i ∈ range n, ∑ j ∈ range n, (if i = j then 0 else Real.log (D i j A)))
      = E A + (∑ i ∈ range n, ∑ j ∈ range n, (if i = j then 0 else alphaP b (cc M i) (cc M j))) * A
        + (∑ i ∈ range n, ∑ j ∈ range n,
            (if i = j then 0 else betaP b (cc M i) (cc M j) (CC M ((j : ℝ) - i)))) * A ^ 2
        + ∑ i ∈ range n, ∑ j ∈ range n, (if i = j then 0 else Real.log (D i j 0)) := by
    simp only [hE_def, Finset.sum_mul, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => ?_))
    split_ifs <;> ring
  simp only [hD] at key
  rw [key, hsa, hsb]
  ring

/-! ### The normalisation `A ↔ ε` -/

/-- The component `⟨q − x₀, e⟩` along the Euclidean-unit kernel vector `e_k = (2/n)^{1/2} cos(lφ_k) ẑ_k`
(paper Def. c4, odd `n`), `⟨a, b⟩ = Re(ā b)`, `ẑ_k = e^{2πik/n}` the ring positions. -/
noncomputable def epsOdd (M : ℕ) (z : Fin (2 * M + 1) → ℂ) : ℝ :=
  ∑ k : Fin (2 * M + 1), √(2 / (2 * M + 1 : ℕ)) * cc M ((k : ℕ) : ℝ)
    * ((starRingEnd ℂ) (ringConf (2 * M + 1) k) * (z k - ringConf (2 * M + 1) k)).re

theorem uOf_hasDerivAt (b c : ℝ) (hb : -1 < b) :
    HasDerivAt (fun A => √(uOf b (w0 b + A * c))) (c * (1 + b) ^ 2 / 2) 0 := by
  have h1 : 0 < 1 + b := by linarith
  have hn : HasDerivAt (fun A : ℝ => w0 b + A * c) c 0 := by
    simpa using ((hasDerivAt_id (0 : ℝ)).mul_const c).const_add (w0 b)
  have hd : HasDerivAt (fun A : ℝ => 1 - b * (w0 b + A * c)) (-(b * c)) 0 := by
    simpa using (hn.const_mul b).const_sub 1
  have hden : 1 - b * (w0 b + 0 * c) ≠ 0 := by
    rw [zero_mul, add_zero]; unfold w0
    rw [show 1 - b * (1 / (1 + b)) = 1 / (1 + b) by field_simp; ring]; positivity
  have hq : HasDerivAt (fun A => (w0 b + A * c) / (1 - b * (w0 b + A * c)))
      ((c * (1 - b * (w0 b + 0 * c)) - (w0 b + 0 * c) * -(b * c)) / (1 - b * (w0 b + 0 * c)) ^ 2) 0 :=
    hn.div hd hden
  have hu0 : (w0 b + 0 * c) / (1 - b * (w0 b + 0 * c)) = 1 := by
    rw [zero_mul, add_zero]; unfold w0; field_simp; ring
  have hs := hq.sqrt (by rw [hu0]; norm_num)
  unfold uOf
  convert hs using 1
  rw [hu0, Real.sqrt_one]
  simp only [zero_mul, add_zero]
  unfold w0
  field_simp
  ring

/-- **Odd `n`: the normalisation.** `dε/dA = (n/2)^{1/2}(1+λ)²/2` at the ring, i.e.
`A = ε (2/n)^{1/2} 2/(1+λ)² + O(ε²)` (paper, proof of Thm odd); `oddNorm_pow4` turns this into the
factor `64/(n²(1+λ)⁸)` of `c4oddSet`. -/
theorem epsOdd_hasDerivAt (M : ℕ) (hM : 1 ≤ M) (b : ℝ) (hb : -1 < b) :
    HasDerivAt (fun A => epsOdd M (oddLine M b A))
      (√(2 / (2 * M + 1 : ℕ)) * ((2 * M + 1 : ℕ) / 2) * ((1 + b) ^ 2 / 2)) 0 := by
  have hterm : ∀ (A : ℝ) (k : Fin (2 * M + 1)),
      ((starRingEnd ℂ) (ringConf (2 * M + 1) k) * (oddLine M b A k - ringConf (2 * M + 1) k)).re
        = √(uOf b (w0 b + A * cc M ((k : ℕ) : ℝ))) - 1 := by
    intro A k
    have hz : oddLine M b A k = (√(uOf b (w0 b + A * cc M ((k : ℕ) : ℝ))) : ℂ) * ringConf (2 * M + 1) k := by
      unfold oddLine ringConf; congr 1; congr 1; push_cast; ring
    have hn : (starRingEnd ℂ) (ringConf (2 * M + 1) k) * ringConf (2 * M + 1) k = 1 := by
      rw [mul_comm, Complex.mul_conj, ringConf, Complex.normSq_eq_norm_sq]
      rw [show 2 * ↑π * Complex.I * ((k : ℕ) : ℂ) / ((2 * M + 1 : ℕ) : ℂ)
          = ↑(2 * π * ((k : ℕ) : ℝ) / (2 * M + 1 : ℕ)) * Complex.I by push_cast; ring,
        Complex.norm_exp_ofReal_mul_I]
      simp
    rw [hz]
    generalize ringConf (2 * M + 1) k = w at hn ⊢
    have : (starRingEnd ℂ) w * ((√(uOf b (w0 b + A * cc M ((k : ℕ) : ℝ))) : ℂ) * w - w)
        = ((√(uOf b (w0 b + A * cc M ((k : ℕ) : ℝ))) - 1 : ℝ) : ℂ) := by
      push_cast; linear_combination ((√(uOf b (w0 b + A * cc M ((k : ℕ) : ℝ))) : ℂ) - 1) * hn
    rw [this, Complex.ofReal_re]
  have e : (fun A => epsOdd M (oddLine M b A))
      = fun A => ∑ k : Fin (2 * M + 1), √(2 / (2 * M + 1 : ℕ)) * cc M ((k : ℕ) : ℝ)
          * (√(uOf b (w0 b + A * cc M ((k : ℕ) : ℝ))) - 1) := by
    funext A; unfold epsOdd; simp_rw [hterm]
  rw [e]
  have hs := HasDerivAt.fun_sum (u := Finset.univ) (fun (k : Fin (2 * M + 1)) _ =>
    ((uOf_hasDerivAt b (cc M ((k : ℕ) : ℝ)) hb).sub_const 1).const_mul
      (√(2 / (2 * M + 1 : ℕ)) * cc M ((k : ℕ) : ℝ)))
  convert hs using 1
  have hsq := sum_cc_sq M hM
  rw [← Fin.sum_univ_eq_sum_range (fun k => cc M (k : ℝ) ^ 2) (2 * M + 1)] at hsq
  rw [← hsq, Finset.mul_sum, Finset.sum_mul]
  refine Finset.sum_congr rfl (fun k _ => ?_)
  ring

/-- The paper's factor: `(dε/dA)⁻⁴ = 64/(n²(1+λ)⁸)` (so `c₄ = 64/(n²(1+λ)⁸)·min E`). -/
theorem oddNorm_pow4 (M : ℕ) (b : ℝ) (hb : -1 < b) :
    (√(2 / (2 * M + 1 : ℕ)) * ((2 * M + 1 : ℕ) / 2) * ((1 + b) ^ 2 / 2))⁻¹ ^ 4
      = 64 / (((2 * M + 1 : ℕ) : ℝ) ^ 2 * (1 + b) ^ 8) := by
  have h1 : 0 < 1 + b := by linarith
  have hn : (0 : ℝ) < ((2 * M + 1 : ℕ) : ℝ) := by positivity
  have hsq : √(2 / ((2 * M + 1 : ℕ) : ℝ)) ^ 2 = 2 / ((2 * M + 1 : ℕ) : ℝ) :=
    Real.sq_sqrt (by positivity)
  rw [inv_pow, show (√(2 / ((2 * M + 1 : ℕ) : ℝ)) * (((2 * M + 1 : ℕ) : ℝ) / 2) * ((1 + b) ^ 2 / 2)) ^ 4
      = (√(2 / ((2 * M + 1 : ℕ) : ℝ)) ^ 2) ^ 2 * (((2 * M + 1 : ℕ) : ℝ) / 2) ^ 4 * ((1 + b) ^ 2 / 2) ^ 4
      by ring, hsq]
  field_simp
  ring

end Havelock.Conj13.Odd
