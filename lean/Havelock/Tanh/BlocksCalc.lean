import Mathlib.Analysis.InnerProductSpace.Calculus
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Analysis.Complex.RealDeriv
import Havelock.Tanh.BlocksSums

/-! # Tanh model, Fourier blocks: the Hamiltonian in log-polar coordinates and its derivatives

The model of note §"The model", eq. (H), in the coordinates of note §"Coordinates and the relative
equilibrium":

* `ringZ N ξ s θ j = √ξ w^j e^{s_j + iθ_j}`, `w = e^{2πi/N}` (`s θ : Fin N → ℝ`);
* `ringH = Σ_{j<k} (−ln|z_j − z_k| + ln|1 − z_j z̄_k|)` (note eq. (H)), `ringJ = Σ_k (1+|z_k|²)/(1−|z_k|²)`,
  `ringF = ringH − Ω ringJ` (`F = H − ΩJ`).

**Notion of derivative used.** For a direction `(s, θ) ∈ ℝ^N × ℝ^N` we restrict `F` to the line
`t ↦ F(t s, t θ)` through the ring (`s = θ = 0`). The first derivative of this one-variable function at
`t = 0` is the directional derivative `dF(s,θ)`, and its iterated second derivative `deriv (deriv ·) 0`
is the Hessian quadratic form `Q(s,θ) = d²F(s,θ)[(s,θ),(s,θ)]`. (`Blocks.lean` recovers the Hessian
entries from `Q` by polarisation.)

Main results:
* `grad_eq`: `d/dt F(ts,tθ)|₀ = Σ_{j<k} lPair_{jk} − Ω Σ_j κ₀ s_j`;
* `hess_eq`: `d²/dt² F(ts,tθ)|₀ = Σ_{j<k} qPair_{jk} − Ω Σ_j 8ξ(1+ξ)/(1−ξ)³ s_j²`, where
  `qPair_{jk} = −α_{jk}((s_j−s_k)² − (θ_j−θ_k)²) − Re ρ_{jk} ((s_j+s_k)² − (θ_j−θ_k)²)
  + 2 Im ρ_{jk} (s_j+s_k)(θ_j−θ_k)`, `ζ_{jk} = w^j conj(w^k)`, `α_{jk} = 1/|1−ζ_{jk}|²`,
  `ρ_{jk} = ξζ_{jk}/(1−ξζ_{jk})²`. (The first two pieces come from the pair term `−ln|z_j−z_k|`, the
  others from the image term `ln|1 − z_j z̄_k|`.)
-/

noncomputable section

namespace Havelock.Tanh

open Filter Topology
open scoped ComplexConjugate

/-! ## Second derivatives along a line -/

/-- `f` has, near `0`, derivative `f1`, and `f1` has derivative `c` at `0`; then
`deriv f 0 = f1 0` and `deriv (deriv f) 0 = c` (`HasD2.deriv1`, `HasD2.deriv2`). -/
def HasD2 (f f1 : ℝ → ℝ) (c : ℝ) : Prop :=
  (∀ᶠ t in 𝓝 (0 : ℝ), HasDerivAt f (f1 t) t) ∧ HasDerivAt f1 c 0

namespace HasD2

variable {f g f1 g1 : ℝ → ℝ} {c d : ℝ}

theorem deriv1 (h : HasD2 f f1 c) : deriv f 0 = f1 0 := h.1.self_of_nhds.deriv

theorem deriv2 (h : HasD2 f f1 c) : deriv (deriv f) 0 = c := by
  have h' : deriv f =ᶠ[𝓝 0] f1 := h.1.mono fun t ht => ht.deriv
  rw [h'.deriv_eq]; exact h.2.deriv

theorem add (hf : HasD2 f f1 c) (hg : HasD2 g g1 d) :
    HasD2 (fun t => f t + g t) (fun t => f1 t + g1 t) (c + d) :=
  ⟨(hf.1.and hg.1).mono fun _ ht => ht.1.add ht.2, hf.2.add hg.2⟩

theorem sub (hf : HasD2 f f1 c) (hg : HasD2 g g1 d) :
    HasD2 (fun t => f t - g t) (fun t => f1 t - g1 t) (c - d) :=
  ⟨(hf.1.and hg.1).mono fun _ ht => ht.1.sub ht.2, hf.2.sub hg.2⟩

theorem neg (hf : HasD2 f f1 c) : HasD2 (fun t => -f t) (fun t => -f1 t) (-c) :=
  ⟨hf.1.mono fun _ ht => ht.neg, hf.2.neg⟩

theorem const_mul (a : ℝ) (hf : HasD2 f f1 c) :
    HasD2 (fun t => a * f t) (fun t => a * f1 t) (a * c) :=
  ⟨hf.1.mono fun _ ht => ht.const_mul a, hf.2.const_mul a⟩

theorem zero : HasD2 (fun _ => 0) (fun _ => 0) 0 :=
  ⟨Eventually.of_forall fun t => hasDerivAt_const t 0, hasDerivAt_const 0 0⟩

theorem ite (p : Prop) [Decidable p] (hf : p → HasD2 f f1 c) :
    HasD2 (fun t => if p then f t else 0) (fun t => if p then f1 t else 0) (if p then c else 0) := by
  by_cases hp : p
  · simp only [if_pos hp]; exact hf hp
  · simp only [if_neg hp]; exact zero

theorem sum {ι : Type*} (s : Finset ι) {F F1 : ι → ℝ → ℝ} {C : ι → ℝ}
    (h : ∀ i ∈ s, HasD2 (F i) (F1 i) (C i)) :
    HasD2 (fun t => ∑ i ∈ s, F i t) (fun t => ∑ i ∈ s, F1 i t) (∑ i ∈ s, C i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using zero
  | insert a s ha ih =>
    simp only [Finset.sum_insert ha]
    exact (h a (Finset.mem_insert_self a s)).add (ih fun i hi => h i (Finset.mem_insert_of_mem hi))

theorem congr_c (h : HasD2 f f1 c) (hc : c = d) : HasD2 f f1 d := hc ▸ h

theorem congr_f (h : HasD2 f f1 c) (hf : ∀ t, f t = g t) : HasD2 g f1 c := by
  have : f = g := funext hf
  exact this ▸ h

end HasD2

theorem hasDerivAt_re_comp {h : ℝ → ℂ} {h' : ℂ} {t : ℝ} (H : HasDerivAt h h' t) :
    HasDerivAt (fun u => (h u).re) h'.re t := by
  have := Complex.reCLM.hasFDerivAt.comp_hasDerivAt t H
  exact this

/-- `d/dt ln|g(t)| = Re(g'/g)` where `g ≠ 0`. -/
theorem hasDerivAt_logNorm {g : ℝ → ℂ} {g1 : ℂ} {t : ℝ} (hg : HasDerivAt g g1 t) (h : g t ≠ 0) :
    HasDerivAt (fun u => Real.log ‖g u‖) ((g1 / g t).re) t := by
  have h2 : ‖g t‖ ^ 2 ≠ 0 := pow_ne_zero 2 (norm_ne_zero_iff.mpr h)
  have H := (hg.norm_sq.log h2).div_const 2
  have e : (fun u => Real.log ‖g u‖) = fun u => Real.log (‖g u‖ ^ 2) / 2 := by
    funext u; rw [Real.log_pow]; push_cast; ring
  rw [e]
  convert H using 1
  have hn : Complex.normSq (g t) ≠ 0 := (Complex.normSq_pos.mpr h).ne'
  rw [Complex.inner, Complex.sq_norm, Complex.div_re]
  simp only [Complex.mul_re, Complex.conj_re, Complex.conj_im]
  field_simp
  ring

/-- `ln|g|` has second derivative `Re((g''g − g'²)/g²)` at `0` (if `g 0 ≠ 0`). -/
theorem hasD2_logNorm {g g1 g2 : ℝ → ℂ} (hg : ∀ t, HasDerivAt g (g1 t) t)
    (hg1 : ∀ t, HasDerivAt g1 (g2 t) t) (h0 : g 0 ≠ 0) :
    HasD2 (fun t => Real.log ‖g t‖) (fun t => (g1 t / g t).re)
      (((g2 0 * g 0 - g1 0 * g1 0) / g 0 ^ 2).re) := by
  refine ⟨?_, ?_⟩
  · have hev : ∀ᶠ t in 𝓝 (0 : ℝ), g t ≠ 0 := (hg 0).continuousAt.eventually_ne h0
    exact hev.mono fun t ht => hasDerivAt_logNorm (hg t) ht
  · exact hasDerivAt_re_comp ((hg1 0).div (hg 0) h0)

/-- `t ↦ c e^{tε}`. -/
def expLin (c ε : ℂ) (t : ℝ) : ℂ := c * Complex.exp (t * ε)

theorem hasDerivAt_expLin (c ε : ℂ) (t : ℝ) : HasDerivAt (expLin c ε) (expLin (c * ε) ε t) t := by
  have h1 : HasDerivAt (fun u : ℝ => (u : ℂ) * ε) ε t := by
    simpa using ((hasDerivAt_id t).ofReal_comp).mul_const ε
  have := h1.cexp.const_mul c
  unfold expLin
  convert this using 1; ring

@[simp] theorem expLin_zero (c ε : ℂ) : expLin c ε 0 = c := by simp [expLin]

/-! ## The model -/

/-- `ε_j = s_j + iθ_j`. -/
def epsC {N : ℕ} (s θ : Fin N → ℝ) (j : Fin N) : ℂ := (s j : ℂ) + (θ j : ℂ) * Complex.I

/-- `z_j = √ξ w^j e^{s_j + iθ_j}` (note §"Coordinates and the relative equilibrium"). -/
def ringZ (N : ℕ) (ξ : ℝ) (s θ : Fin N → ℝ) (j : Fin N) : ℂ :=
  (Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ) * Complex.exp ((s j : ℂ) + (θ j : ℂ) * Complex.I)

/-- The tanh Hamiltonian of `N` unit vortices on the Poincaré disk (note eq. (H)):
`H = Σ_{j<k} (−ln|z_j − z_k| + ln|1 − z_j z̄_k|)`. -/
def ringH (N : ℕ) (ξ : ℝ) (s θ : Fin N → ℝ) : ℝ :=
  ∑ j : Fin N, ∑ k : Fin N, if j < k then
    (-Real.log ‖ringZ N ξ s θ j - ringZ N ξ s θ k‖
      + Real.log ‖1 - ringZ N ξ s θ j * conj (ringZ N ξ s θ k)‖) else 0

/-- `J = Σ_k (1+|z_k|²)/(1−|z_k|²)` (note eq. (H)). -/
def ringJ (N : ℕ) (ξ : ℝ) (s θ : Fin N → ℝ) : ℝ :=
  ∑ k : Fin N, (1 + ‖ringZ N ξ s θ k‖ ^ 2) / (1 - ‖ringZ N ξ s θ k‖ ^ 2)

/-- `F = H − ΩJ` as a function of the `2N` real variables `(s, θ)`. -/
def ringF (N : ℕ) (ξ Ω : ℝ) (s θ : Fin N → ℝ) : ℝ := ringH N ξ s θ - Ω * ringJ N ξ s θ

/-- `κ₀ = 4ξ/(1−ξ)²` (note §"Coordinates and the relative equilibrium"). -/
def kappa0 (ξ : ℝ) : ℝ := 4 * ξ / (1 - ξ) ^ 2

/-- `ζ_{jk} = w^j conj(w^k) = w^{j−k}`. -/
def zeta (N : ℕ) (j k : Fin N) : ℂ := ringW N ^ (j : ℕ) * conj (ringW N ^ (k : ℕ))

/-- `α(ζ) = 1/|1−ζ|²`. -/
def alphaF (ζ : ℂ) : ℝ := 1 / Complex.normSq (1 - ζ)

/-- Second derivative of the `(j,k)` pair + image term along the line, at `0`. -/
def qPair (N : ℕ) (ξ : ℝ) (s θ : Fin N → ℝ) (j k : Fin N) : ℝ :=
  -alphaF (zeta N j k) * ((s j - s k) ^ 2 - (θ j - θ k) ^ 2)
    - (rho ξ (zeta N j k)).re * ((s j + s k) ^ 2 - (θ j - θ k) ^ 2)
    + 2 * (rho ξ (zeta N j k)).im * (s j + s k) * (θ j - θ k)

/-- First derivative of the `(j,k)` pair + image term along the line, at `0`. -/
def lPair (N : ℕ) (ξ : ℝ) (s θ : Fin N → ℝ) (j k : Fin N) : ℝ :=
  -(((ringW N ^ (j : ℕ) * epsC s θ j - ringW N ^ (k : ℕ) * epsC s θ k)
      / (ringW N ^ (j : ℕ) - ringW N ^ (k : ℕ))).re)
    - (rho1 ξ (zeta N j k) * (epsC s θ j + conj (epsC s θ k))).re

/-! ## Restriction to a line through the ring -/

theorem ringZ_line {N : ℕ} (ξ : ℝ) (s θ : Fin N → ℝ) (j : Fin N) (t : ℝ) :
    ringZ N ξ (t • s) (t • θ) j = expLin ((Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ)) (epsC s θ j) t := by
  unfold ringZ expLin epsC
  simp only [Pi.smul_apply, smul_eq_mul]
  congr 2
  push_cast; ring

theorem ringZ_mul_conj_line {N : ℕ} {ξ : ℝ} (h0 : 0 ≤ ξ) (s θ : Fin N → ℝ) (j k : Fin N) (t : ℝ) :
    ringZ N ξ (t • s) (t • θ) j * conj (ringZ N ξ (t • s) (t • θ) k)
      = expLin ((ξ : ℂ) * zeta N j k) (epsC s θ j + conj (epsC s θ k)) t := by
  rw [ringZ_line, ringZ_line]
  unfold expLin zeta
  have hs : ((Real.sqrt ξ : ℝ) : ℂ) * ((Real.sqrt ξ : ℝ) : ℂ) = (ξ : ℂ) := by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt h0]
  have he : conj (Complex.exp ((t : ℂ) * epsC s θ k)) = Complex.exp ((t : ℂ) * conj (epsC s θ k)) := by
    rw [← Complex.exp_conj, map_mul, Complex.conj_ofReal]
  simp only [map_mul, Complex.conj_ofReal, he]
  rw [mul_add, Complex.exp_add, ← hs]
  ring

theorem norm_ringZ_sq_line {N : ℕ} {ξ : ℝ} (hN : 0 < N) (h0 : 0 ≤ ξ) (s θ : Fin N → ℝ) (j : Fin N)
    (t : ℝ) : ‖ringZ N ξ (t • s) (t • θ) j‖ ^ 2 = ξ * Real.exp (2 * t * s j) := by
  unfold ringZ
  rw [norm_mul, norm_mul, norm_ringW_pow hN, Complex.norm_exp, Complex.norm_real,
    Real.norm_of_nonneg (Real.sqrt_nonneg _)]
  simp only [Pi.smul_apply, smul_eq_mul, Complex.add_re, Complex.ofReal_re, Complex.mul_re,
    Complex.ofReal_im, Complex.I_re, Complex.I_im, mul_zero, mul_one, sub_zero, add_zero]
  rw [mul_pow, Real.sq_sqrt h0, ← Real.exp_nat_mul]
  push_cast; ring_nf

/-! ## The three kinds of terms -/

/-- The pair term `−ln|z_j − z_k|` along the line. -/
theorem hasD2_pairLog (a b ε δ : ℂ) (hab : a ≠ b) :
    HasD2 (fun t => -Real.log ‖expLin a ε t - expLin b δ t‖)
      (fun t => -((expLin (a * ε) ε t - expLin (b * δ) δ t) / (expLin a ε t - expLin b δ t)).re)
      (-(((a * ε * ε * (a - b) - b * δ * δ * (a - b) - (a * ε - b * δ) * (a * ε - b * δ))
          / (a - b) ^ 2).re)) := by
  have hg : ∀ t, HasDerivAt (fun t => expLin a ε t - expLin b δ t)
      (expLin (a * ε) ε t - expLin (b * δ) δ t) t :=
    fun t => (hasDerivAt_expLin a ε t).sub (hasDerivAt_expLin b δ t)
  have hg1 : ∀ t, HasDerivAt (fun t => expLin (a * ε) ε t - expLin (b * δ) δ t)
      (expLin (a * ε * ε) ε t - expLin (b * δ * δ) δ t) t :=
    fun t => (hasDerivAt_expLin _ ε t).sub (hasDerivAt_expLin _ δ t)
  have h0 : expLin a ε 0 - expLin b δ 0 ≠ 0 := by simpa using sub_ne_zero.mpr hab
  have := (hasD2_logNorm hg hg1 h0).neg
  refine this.congr_c ?_
  simp only [expLin_zero]
  congr 3
  ring

/-- The image term `ln|1 − z_j z̄_k|` along the line. -/
theorem hasD2_imageLog (C μ : ℂ) (hC : C ≠ 1) :
    HasD2 (fun t => Real.log ‖1 - expLin C μ t‖)
      (fun t => ((-expLin (C * μ) μ t) / (1 - expLin C μ t)).re)
      ((-(C / (1 - C) ^ 2) * μ ^ 2).re) := by
  have hg : ∀ t, HasDerivAt (fun t => 1 - expLin C μ t) (-expLin (C * μ) μ t) t := by
    intro t; simpa using (hasDerivAt_expLin C μ t).const_sub 1
  have hg1 : ∀ t, HasDerivAt (fun t => -expLin (C * μ) μ t) (-expLin (C * μ * μ) μ t) t :=
    fun t => (hasDerivAt_expLin _ μ t).neg
  have h0 : 1 - expLin C μ 0 ≠ 0 := by simpa using sub_ne_zero.mpr (Ne.symm hC)
  refine (hasD2_logNorm hg hg1 h0).congr_c ?_
  simp only [expLin_zero]
  have h1 : (1 : ℂ) - C ≠ 0 := sub_ne_zero.mpr (Ne.symm hC)
  congr 1
  field_simp
  ring

/-- The one-body term `(1+u)/(1−u)`, `u = ξ e^{2tσ}`. -/
theorem hasD2_J (ξ σ : ℝ) (h1 : ξ < 1) :
    HasD2 (fun t => (1 + ξ * Real.exp (2 * t * σ)) / (1 - ξ * Real.exp (2 * t * σ)))
      (fun t => 4 * σ * (ξ * Real.exp (2 * t * σ)) / (1 - ξ * Real.exp (2 * t * σ)) ^ 2)
      (8 * σ ^ 2 * ξ * (1 + ξ) / (1 - ξ) ^ 3) := by
  have hu : ∀ t, HasDerivAt (fun t => ξ * Real.exp (2 * t * σ)) (2 * σ * (ξ * Real.exp (2 * t * σ))) t := by
    intro t
    have h2 : HasDerivAt (fun t => 2 * t * σ) (2 * σ) t := by
      simpa using ((hasDerivAt_id t).const_mul 2).mul_const σ
    have := (h2.exp).const_mul ξ
    convert this using 1; ring
  have hcont : ContinuousAt (fun t => ξ * Real.exp (2 * t * σ)) 0 := (hu 0).continuousAt
  have hev : ∀ᶠ t in 𝓝 (0 : ℝ), 1 - ξ * Real.exp (2 * t * σ) ≠ 0 := by
    have : ∀ᶠ t in 𝓝 (0 : ℝ), ξ * Real.exp (2 * t * σ) ≠ 1 :=
      hcont.eventually_ne (by simp; exact h1.ne)
    exact this.mono fun t ht => sub_ne_zero.mpr (Ne.symm ht)
  refine ⟨hev.mono fun t ht => ?_, ?_⟩
  · have := ((hu t).const_add 1).div ((hu t).const_sub 1) ht
    convert this using 1
    field_simp
    ring
  · have hne : (1 : ℝ) - ξ ≠ 0 := by linarith
    have h0 : 1 - ξ * Real.exp (2 * 0 * σ) ≠ 0 := by simpa using hne
    have := ((hu 0).const_mul (4 * σ)).div ((hu 0).const_sub 1 |>.pow 2) (pow_ne_zero 2 h0)
    convert this using 1
    simp only [Pi.pow_apply, mul_zero, zero_mul, Real.exp_zero, mul_one]
    push_cast
    field_simp
    ring

/-! ## Value lemmas -/

theorem re_sq_sub {N : ℕ} (s θ : Fin N → ℝ) (j k : Fin N) :
    ((epsC s θ j - epsC s θ k) ^ 2).re = (s j - s k) ^ 2 - (θ j - θ k) ^ 2 := by
  simp [epsC, sq, Complex.mul_re]

theorem re_mu {N : ℕ} (s θ : Fin N → ℝ) (j k : Fin N) :
    (epsC s θ j + conj (epsC s θ k)).re = s j + s k := by simp [epsC]

theorem im_mu {N : ℕ} (s θ : Fin N → ℝ) (j k : Fin N) :
    (epsC s θ j + conj (epsC s θ k)).im = θ j - θ k := by simp [epsC]; ring

/-- `−ab/(a−b)² = 1/|1−a b̄|²` for `|a| = |b| = 1` (here `a = w^j`, `b = w^k`). -/
theorem alpha_eq {N : ℕ} (j k : Fin N) (hjk : ringW N ^ (j : ℕ) ≠ ringW N ^ (k : ℕ)) (c : ℂ)
    (hc : c ≠ 0) :
    -(c * ringW N ^ (j : ℕ) * (c * ringW N ^ (k : ℕ)))
        / (c * ringW N ^ (j : ℕ) - c * ringW N ^ (k : ℕ)) ^ 2
      = ((alphaF (zeta N j k) : ℝ) : ℂ) := by
  set a := ringW N ^ (j : ℕ)
  set b := ringW N ^ (k : ℕ)
  have ha : a ≠ 0 := pow_ne_zero _ (ringW_ne_zero N)
  have hb : b ≠ 0 := pow_ne_zero _ (ringW_ne_zero N)
  have hab : a - b ≠ 0 := sub_ne_zero.mpr hjk
  have hca : conj a = a⁻¹ := conj_pow_ringW N j
  have hcb : conj b = b⁻¹ := conj_pow_ringW N k
  unfold alphaF zeta
  push_cast
  rw [← Complex.mul_conj]
  simp only [map_sub, map_one, map_mul, Complex.conj_conj]
  rw [hca, hcb]
  have h1 : (1 - a * b⁻¹) * (1 - a⁻¹ * b) = -(a - b) ^ 2 / (a * b) := by
    field_simp; ring
  rw [h1]
  field_simp

theorem rho_conj (ξ : ℝ) (ζ : ℂ) : rho ξ (conj ζ) = conj (rho ξ ζ) := by
  unfold rho; simp [map_div₀, map_mul, map_pow, map_sub, map_one, Complex.conj_ofReal]

theorem rho1_conj (ξ : ℝ) (ζ : ℂ) : rho1 ξ (conj ζ) = conj (rho1 ξ ζ) := by
  unfold rho1; simp [map_div₀, map_mul, map_sub, map_one, Complex.conj_ofReal]

theorem xi_zeta_ne_one {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) (j k : Fin N) :
    (ξ : ℂ) * zeta N j k ≠ 1 := by
  intro h
  have hn : ‖(ξ : ℂ) * zeta N j k‖ = ξ := by
    unfold zeta
    rw [norm_mul, norm_mul, norm_ringW_pow hN, RCLike.norm_conj, norm_ringW_pow hN, Complex.norm_real,
      Real.norm_of_nonneg h0]; ring
  rw [h, norm_one] at hn; linarith

theorem ringW_pow_ne {N : ℕ} (hN : 0 < N) {j k : Fin N} (hjk : j ≠ k) :
    ringW N ^ (j : ℕ) ≠ ringW N ^ (k : ℕ) := by
  intro h
  exact hjk (Fin.ext ((ringW_prim hN).pow_inj j.2 k.2 h))

/-- Per-pair second derivative (pair + image term). -/
theorem hasD2_pairTerm {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (s θ : Fin N → ℝ)
    {j k : Fin N} (hjk : j ≠ k) :
    HasD2 (fun t => -Real.log ‖ringZ N ξ (t • s) (t • θ) j - ringZ N ξ (t • s) (t • θ) k‖
        + Real.log ‖1 - ringZ N ξ (t • s) (t • θ) j * conj (ringZ N ξ (t • s) (t • θ) k)‖)
      (fun t =>
        -((expLin ((Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ) * epsC s θ j) (epsC s θ j) t
            - expLin ((Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ) * epsC s θ k) (epsC s θ k) t)
          / (expLin ((Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ)) (epsC s θ j) t
            - expLin ((Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ)) (epsC s θ k) t)).re
        + ((-expLin ((ξ : ℂ) * zeta N j k * (epsC s θ j + conj (epsC s θ k)))
              (epsC s θ j + conj (epsC s θ k)) t)
            / (1 - expLin ((ξ : ℂ) * zeta N j k) (epsC s θ j + conj (epsC s θ k)) t)).re)
      (qPair N ξ s θ j k) := by
  have hsq : ((Real.sqrt ξ : ℝ) : ℂ) ≠ 0 := by
    rw [Complex.ofReal_ne_zero]; exact (Real.sqrt_pos.mpr h0).ne'
  have hw := ringW_pow_ne hN hjk
  have hab : (Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ) ≠ (Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ) := by
    intro h; exact hw (mul_left_cancel₀ hsq h)
  have HP := hasD2_pairLog _ _ (epsC s θ j) (epsC s θ k) hab
  have HI := hasD2_imageLog ((ξ : ℂ) * zeta N j k) (epsC s θ j + conj (epsC s θ k))
    (xi_zeta_ne_one hN h0.le h1 j k)
  have H := (HP.add HI).congr_f (g := fun t =>
      -Real.log ‖ringZ N ξ (t • s) (t • θ) j - ringZ N ξ (t • s) (t • θ) k‖
        + Real.log ‖1 - ringZ N ξ (t • s) (t • θ) j * conj (ringZ N ξ (t • s) (t • θ) k)‖)
    (fun t => by rw [ringZ_mul_conj_line h0.le, ringZ_line, ringZ_line])
  refine H.congr_c ?_
  -- evaluate the constants
  set a := (Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ)
  set b := (Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ)
  have hab' : a - b ≠ 0 := sub_ne_zero.mpr hab
  have e1 : (a * epsC s θ j * epsC s θ j * (a - b) - b * epsC s θ k * epsC s θ k * (a - b)
      - (a * epsC s θ j - b * epsC s θ k) * (a * epsC s θ j - b * epsC s θ k)) / (a - b) ^ 2
      = (-(a * b) / (a - b) ^ 2) * (epsC s θ j - epsC s θ k) ^ 2 := by
    field_simp; ring
  have e2 : -(a * b) / (a - b) ^ 2 = ((alphaF (zeta N j k) : ℝ) : ℂ) := by
    have := alpha_eq j k hw _ hsq
    simpa [a, b, mul_assoc, mul_comm, mul_left_comm] using this
  rw [e1, e2, Complex.re_ofReal_mul, re_sq_sub]
  have e3 : (ξ : ℂ) * zeta N j k / (1 - (ξ : ℂ) * zeta N j k) ^ 2 = rho ξ (zeta N j k) := rfl
  rw [e3]
  simp only [neg_mul, Complex.neg_re, Complex.mul_re, sq, re_mu, im_mu, Complex.mul_im]
  unfold qPair
  ring

/-! ## The full function along a line -/

/-- `F(ts, tθ)` has, at `t = 0`, first derivative `Σ_{j<k} lPair_{jk} − Ω Σ_k κ₀ s_k` and second
derivative `Σ_{j<k} qPair_{jk} − Ω Σ_k 8ξ(1+ξ)/(1−ξ)³ s_k²`. -/
theorem hasD2_ringF {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ)
    (s θ : Fin N → ℝ) :
    ∃ f1 : ℝ → ℝ, HasD2 (fun t => ringF N ξ Ω (t • s) (t • θ)) f1
        (∑ j : Fin N, ∑ k : Fin N, (if j < k then qPair N ξ s θ j k else 0)
          - Ω * ∑ k : Fin N, 8 * s k ^ 2 * ξ * (1 + ξ) / (1 - ξ) ^ 3)
      ∧ f1 0 = ∑ j : Fin N, ∑ k : Fin N, (if j < k then lPair N ξ s θ j k else 0)
          - Ω * ∑ k : Fin N, kappa0 ξ * s k := by
  have HP := HasD2.sum (Finset.univ : Finset (Fin N)) (fun j _ =>
    HasD2.sum (Finset.univ : Finset (Fin N)) (fun k _ =>
      HasD2.ite (j < k) (fun hjk => hasD2_pairTerm hN h0 h1 s θ (ne_of_lt hjk))))
  have HJ := (HasD2.sum (Finset.univ : Finset (Fin N)) (fun k _ => hasD2_J ξ (s k) h1)).const_mul Ω
  have H := HP.sub HJ
  refine ⟨_, H.congr_f (fun t => ?_), ?_⟩
  · unfold ringF ringH ringJ
    simp only [norm_ringZ_sq_line hN h0.le]
  · have hsq : ((Real.sqrt ξ : ℝ) : ℂ) ≠ 0 := by
      rw [Complex.ofReal_ne_zero]; exact (Real.sqrt_pos.mpr h0).ne'
    simp only [expLin_zero, mul_zero, zero_mul, Real.exp_zero, mul_one]
    congr 1
    · refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => ?_
      split_ifs with hjk
      · have hw := ringW_pow_ne hN (ne_of_lt hjk)
        have hw' : ringW N ^ (j : ℕ) - ringW N ^ (k : ℕ) ≠ 0 := sub_ne_zero.mpr hw
        have e1 : ((Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ) * epsC s θ j
              - (Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ) * epsC s θ k)
            / ((Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ) - (Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ))
            = (ringW N ^ (j : ℕ) * epsC s θ j - ringW N ^ (k : ℕ) * epsC s θ k)
              / (ringW N ^ (j : ℕ) - ringW N ^ (k : ℕ)) := by
          rw [show (Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ) * epsC s θ j
              - (Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ) * epsC s θ k
              = (Real.sqrt ξ : ℂ) * (ringW N ^ (j : ℕ) * epsC s θ j - ringW N ^ (k : ℕ) * epsC s θ k)
              by ring,
            show (Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ) - (Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ)
              = (Real.sqrt ξ : ℂ) * (ringW N ^ (j : ℕ) - ringW N ^ (k : ℕ)) by ring,
            mul_div_mul_left _ _ hsq]
        have e2 : -((ξ : ℂ) * zeta N j k * (epsC s θ j + conj (epsC s θ k)))
              / (1 - (ξ : ℂ) * zeta N j k)
            = -(rho1 ξ (zeta N j k) * (epsC s θ j + conj (epsC s θ k))) := by
          unfold rho1; ring
        rw [e1, e2, Complex.neg_re]
        unfold lPair
        ring
      · rfl
    · congr 1
      refine Finset.sum_congr rfl fun k _ => ?_
      unfold kappa0
      ring

/-- First derivative along a line through the ring (directional derivative `dF(s,θ)`). -/
theorem grad_eq {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ) (s θ : Fin N → ℝ) :
    deriv (fun t : ℝ => ringF N ξ Ω (t • s) (t • θ)) 0
      = ∑ j : Fin N, ∑ k : Fin N, (if j < k then lPair N ξ s θ j k else 0)
          - Ω * ∑ k : Fin N, kappa0 ξ * s k := by
  obtain ⟨f1, H, hf⟩ := hasD2_ringF hN h0 h1 Ω s θ
  rw [H.deriv1, hf]

/-- Second derivative along a line through the ring (Hessian quadratic form `Q(s,θ)`). -/
theorem hess_eq {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ) (s θ : Fin N → ℝ) :
    deriv (deriv (fun t : ℝ => ringF N ξ Ω (t • s) (t • θ))) 0
      = ∑ j : Fin N, ∑ k : Fin N, (if j < k then qPair N ξ s θ j k else 0)
          - Ω * ∑ k : Fin N, 8 * s k ^ 2 * ξ * (1 + ξ) / (1 - ξ) ^ 3 := by
  obtain ⟨f1, H, _⟩ := hasD2_ringF hN h0 h1 Ω s θ
  exact H.deriv2

end Havelock.Tanh
