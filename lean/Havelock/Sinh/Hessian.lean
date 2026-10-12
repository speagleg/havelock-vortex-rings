import Havelock.Tanh.BlocksFourier

/-! # Sinh model: the ring is critical and its Hessian is diagonal on Fourier modes (note Prop. `prop:modes`)

The sinh model of `sinh/note.tex`, eq. (H): `N` unit vortices on the Poincaré disk with
`H = −Σ_{j<k} ln|z_j − z_k| + ((N−1)/2) Σ_k ln(1 − |z_k|²)` (`sinhH`; this is
`−Σ_{j<k} ln sinh(d_{jk}/2)` up to a constant) and `J = Σ_k (1+|z_k|²)/(1−|z_k|²)`, `F = H − ΩJ`
(`sinhF`). The coordinates are those of the note §2 and of the tanh development (`Havelock.Tanh.ringZ`):
`z_j = √ξ w^j e^{s_j + iθ_j}`, `w = e^{2πi/N}`; the ring is `s = θ = 0`.

We reuse the tanh machinery (`Havelock/Tanh/Blocks*.lean`): the pair term is the tanh pair term, and the
tanh image term vanishes identically at `ξ = 0` (`rho 0 = rho1 0 = 0`), so the tanh pair-only lattice sums,
Hessian matrices and Fourier identities are those of the tanh model evaluated at `ξ = 0`.

**Notions used** (as in `Havelock/Tanh/Blocks.lean`): derivatives along lines `t ↦ F(ts, tθ)`; the Hessian
quadratic form `sinhHessQ s θ = d²/dt² F(ts,tθ)|₀`, its polarisation `sinhHessB`, and the entries
`sinhHessSS a b = ∂²F/∂s_a∂s_b`, `sinhHessTT`, `sinhHessST`, `sinhHessTS`.

Main results (all for `0 < ξ < 1`):
* `sinh_grad_formula`, `sinh_ring_critical`, `sinh_ring_critical_iff` (note §2, "One-body terms"): the
  derivative along `(s, θ)` is `−((N−1)(1+ξ)/(2(1−ξ)) + Ωκ₀) Σ s_j`, so the ring is critical exactly for
  `Ω = sinhOmega N ξ = −(N−1)(1+ξ)(1−ξ)/(8ξ)`.
* `sinhD_omega`: the one-body diagonal `W''(0) − ΩJ₁''(0)` equals `C₁(ξ) = (N−1)(1+ξ²)/(1−ξ)²`
  `= (N−1)(1+c²)/2`, `c = cR ξ = (1+ξ)/(1−ξ)`.
* `sinh_blocks_eigen` (note Prop. `prop:modes`): for `0 ≤ m ≤ N` and `v_b = w^{mb}`:
  `Σ_b ∂²F/∂s_a∂s_b v_b = a_m v_a`, `Σ_b ∂²F/∂θ_a∂θ_b v_b = b_m v_a`, and the mixed blocks vanish
  (`g_m = 0`), with `a_m = sinhA N m ξ = (N−1)(1+c²)/2 − m(N−m)/2`, `b_m = sinhB N m = m(N−m)/2`.
  `sinhA_eq_C1`: `a_m = C₁(ξ) − m(N−m)/2` (the note's form).
* `sinh_blocks_hermitian`: the same in Hermitian-form version (cross terms between modes vanish).
-/

noncomputable section

namespace Havelock.Sinh

open Havelock.Tanh
open Filter Topology
open scoped ComplexConjugate

/-! ## The model -/

/-- The sinh Hamiltonian (note eq. (H)) in the coordinates `z_j = √ξ w^j e^{s_j + iθ_j}`:
`H = −Σ_{j<k} ln|z_j − z_k| + ((N−1)/2) Σ_k ln(1 − |z_k|²)`. -/
def sinhH (N : ℕ) (ξ : ℝ) (s θ : Fin N → ℝ) : ℝ :=
  (∑ j : Fin N, ∑ k : Fin N,
      if j < k then -Real.log ‖ringZ N ξ s θ j - ringZ N ξ s θ k‖ else 0)
    + ((N : ℝ) - 1) / 2 * ∑ k : Fin N, Real.log (1 - ‖ringZ N ξ s θ k‖ ^ 2)

/-- `F = H − ΩJ` (`J = ringJ`, the generator of rotations, note §1). -/
def sinhF (N : ℕ) (ξ Ω : ℝ) (s θ : Fin N → ℝ) : ℝ := sinhH N ξ s θ - Ω * ringJ N ξ s θ

/-- The angular velocity of note §2: `Ω = −(N−1)(1+ξ)(1−ξ)/(8ξ)`. -/
def sinhOmega (N : ℕ) (ξ : ℝ) : ℝ := -(((N : ℝ) - 1) * (1 + ξ) * (1 - ξ)) / (8 * ξ)

/-- `C₁(ξ) = (N−1)(1+ξ²)/(1−ξ)²` (note §2). -/
def sinhC1 (N : ℕ) (ξ : ℝ) : ℝ := ((N : ℝ) - 1) * (1 + ξ ^ 2) / (1 - ξ) ^ 2

/-- `a_m = (N−1)(1+c²)/2 − m(N−m)/2`, `c = (1+ξ)/(1−ξ)` (note Prop. `prop:modes`, via `sinhA_eq_C1`). -/
def sinhA (N m : ℕ) (ξ : ℝ) : ℝ := ((N : ℝ) - 1) * (1 + cR ξ ^ 2) / 2 - (m : ℝ) * ((N : ℝ) - m) / 2

/-- `b_m = m(N−m)/2` (note Prop. `prop:modes`). -/
def sinhB (N m : ℕ) : ℝ := (m : ℝ) * ((N : ℝ) - m) / 2

theorem sinhA_eq_C1 (N m : ℕ) {ξ : ℝ} (h1 : ξ ≠ 1) :
    sinhA N m ξ = sinhC1 N ξ - (m : ℝ) * ((N : ℝ) - m) / 2 := by
  have hne : (1 : ℝ) - ξ ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  unfold sinhA sinhC1 cR
  field_simp
  ring

theorem rho_zero (ζ : ℂ) : rho 0 ζ = 0 := by simp [rho]

theorem rho1_zero (ζ : ℂ) : rho1 0 ζ = 0 := by simp [rho1]

/-! ## Derivatives along a line -/

theorem hasDerivAt_u (ξ σ t : ℝ) :
    HasDerivAt (fun t => ξ * Real.exp (2 * t * σ)) (2 * σ * (ξ * Real.exp (2 * t * σ))) t := by
  have h2 : HasDerivAt (fun t => 2 * t * σ) (2 * σ) t := by
    simpa using ((hasDerivAt_id t).const_mul 2).mul_const σ
  have := (h2.exp).const_mul ξ
  convert this using 1; ring

/-- The background one-body term `ln(1 − u)`, `u = ξ e^{2tσ}`. -/
theorem hasD2_logW (ξ σ : ℝ) (h1 : ξ < 1) :
    HasD2 (fun t => Real.log (1 - ξ * Real.exp (2 * t * σ)))
      (fun t => -(2 * σ * (ξ * Real.exp (2 * t * σ))) / (1 - ξ * Real.exp (2 * t * σ)))
      (-(4 * σ ^ 2 * ξ / (1 - ξ) ^ 2)) := by
  have hu := hasDerivAt_u ξ σ
  have hcont : ContinuousAt (fun t => ξ * Real.exp (2 * t * σ)) 0 := (hu 0).continuousAt
  have hev : ∀ᶠ t in 𝓝 (0 : ℝ), 1 - ξ * Real.exp (2 * t * σ) ≠ 0 := by
    have : ∀ᶠ t in 𝓝 (0 : ℝ), ξ * Real.exp (2 * t * σ) ≠ 1 :=
      hcont.eventually_ne (by simp; exact h1.ne)
    exact this.mono fun t ht => sub_ne_zero.mpr (Ne.symm ht)
  refine ⟨hev.mono fun t ht => ?_, ?_⟩
  · exact ((hu t).const_sub 1).log ht
  · have hne : (1 : ℝ) - ξ ≠ 0 := by linarith
    have h0 : 1 - ξ * Real.exp (2 * 0 * σ) ≠ 0 := by simpa using hne
    have := (((hu 0).const_mul (2 * σ)).neg).div ((hu 0).const_sub 1) h0
    convert this using 1
    simp only [Pi.neg_apply, mul_zero, zero_mul, Real.exp_zero, mul_one]
    field_simp
    ring

/-- First derivative of the pair term `−ln|z_j − z_k|` along the line. -/
def pairD1 (N : ℕ) (ξ : ℝ) (s θ : Fin N → ℝ) (j k : Fin N) (t : ℝ) : ℝ :=
  -((expLin ((Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ) * epsC s θ j) (epsC s θ j) t
        - expLin ((Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ) * epsC s θ k) (epsC s θ k) t)
      / (expLin ((Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ)) (epsC s θ j) t
        - expLin ((Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ)) (epsC s θ k) t)).re

/-- The pair term: its second derivative is the tanh `qPair` at `ξ = 0` (no image term). -/
theorem hasD2_pairS {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (s θ : Fin N → ℝ)
    {j k : Fin N} (hjk : j ≠ k) :
    HasD2 (fun t => -Real.log ‖ringZ N ξ (t • s) (t • θ) j - ringZ N ξ (t • s) (t • θ) k‖)
      (pairD1 N ξ s θ j k) (qPair N 0 s θ j k) := by
  have hsq : ((Real.sqrt ξ : ℝ) : ℂ) ≠ 0 := by
    rw [Complex.ofReal_ne_zero]; exact (Real.sqrt_pos.mpr h0).ne'
  have hw := ringW_pow_ne hN hjk
  have hab : (Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ) ≠ (Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ) := by
    intro h; exact hw (mul_left_cancel₀ hsq h)
  have HP := hasD2_pairLog _ _ (epsC s θ j) (epsC s θ k) hab
  have H := HP.congr_f (g := fun t =>
      -Real.log ‖ringZ N ξ (t • s) (t • θ) j - ringZ N ξ (t • s) (t • θ) k‖)
    (fun t => by rw [ringZ_line, ringZ_line])
  refine H.congr_c ?_
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
  unfold qPair
  simp only [rho_zero, Complex.zero_re, Complex.zero_im]
  ring

theorem pairD1_zero {N : ℕ} {ξ : ℝ} (h0 : 0 < ξ) (s θ : Fin N → ℝ) (j k : Fin N) :
    pairD1 N ξ s θ j k 0 = lPair N 0 s θ j k := by
  have hsq : ((Real.sqrt ξ : ℝ) : ℂ) ≠ 0 := by
    rw [Complex.ofReal_ne_zero]; exact (Real.sqrt_pos.mpr h0).ne'
  unfold pairD1 lPair
  simp only [expLin_zero, rho1_zero, zero_mul, Complex.zero_re, sub_zero]
  rw [show (Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ) * epsC s θ j
      - (Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ) * epsC s θ k
      = (Real.sqrt ξ : ℂ) * (ringW N ^ (j : ℕ) * epsC s θ j - ringW N ^ (k : ℕ) * epsC s θ k)
      by ring,
    show (Real.sqrt ξ : ℂ) * ringW N ^ (j : ℕ) - (Real.sqrt ξ : ℂ) * ringW N ^ (k : ℕ)
      = (Real.sqrt ξ : ℂ) * (ringW N ^ (j : ℕ) - ringW N ^ (k : ℕ)) by ring,
    mul_div_mul_left _ _ hsq]

/-- `F(ts, tθ)` has, at `t = 0`, first derivative
`Σ_{j<k} lPair⁰_{jk} + ((N−1)/2) Σ_k (−2ξ s_k/(1−ξ)) − Ω Σ_k κ₀ s_k` and second derivative
`Σ_{j<k} qPair⁰_{jk} + ((N−1)/2) Σ_k (−4ξ s_k²/(1−ξ)²) − Ω Σ_k 8ξ(1+ξ) s_k²/(1−ξ)³`
(`lPair⁰ = lPair N 0`, `qPair⁰ = qPair N 0`: the tanh pair terms without image). -/
theorem hasD2_sinhF {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ)
    (s θ : Fin N → ℝ) :
    ∃ f1 : ℝ → ℝ, HasD2 (fun t => sinhF N ξ Ω (t • s) (t • θ)) f1
        (∑ j : Fin N, ∑ k : Fin N, (if j < k then qPair N 0 s θ j k else 0)
          + ((N : ℝ) - 1) / 2 * ∑ k : Fin N, -(4 * s k ^ 2 * ξ / (1 - ξ) ^ 2)
          - Ω * ∑ k : Fin N, 8 * s k ^ 2 * ξ * (1 + ξ) / (1 - ξ) ^ 3)
      ∧ f1 0 = ∑ j : Fin N, ∑ k : Fin N, (if j < k then lPair N 0 s θ j k else 0)
          + ((N : ℝ) - 1) / 2 * ∑ k : Fin N, -(2 * s k * ξ) / (1 - ξ)
          - Ω * ∑ k : Fin N, kappa0 ξ * s k := by
  have HP := HasD2.sum (Finset.univ : Finset (Fin N)) (fun j _ =>
    HasD2.sum (Finset.univ : Finset (Fin N)) (fun k _ =>
      HasD2.ite (j < k) (fun hjk => hasD2_pairS hN h0 s θ (ne_of_lt hjk))))
  have HW := (HasD2.sum (Finset.univ : Finset (Fin N))
    (fun k _ => hasD2_logW ξ (s k) h1)).const_mul (((N : ℝ) - 1) / 2)
  have HJ := (HasD2.sum (Finset.univ : Finset (Fin N)) (fun k _ => hasD2_J ξ (s k) h1)).const_mul Ω
  have H := (HP.add HW).sub HJ
  refine ⟨_, H.congr_f (fun t => ?_), ?_⟩
  · unfold sinhF sinhH ringJ
    simp only [norm_ringZ_sq_line hN h0.le]
  · simp only [mul_zero, zero_mul, Real.exp_zero, mul_one]
    congr 1
    congr 1
    · refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => ?_
      split_ifs with hjk
      · exact pairD1_zero h0 s θ j k
      · rfl
    · congr 1
      refine Finset.sum_congr rfl fun k _ => ?_
      unfold kappa0
      ring

/-! ## Criticality (note §2, "One-body terms") -/

/-- The directional derivative of `F` at the ring along `(s, θ)`. -/
theorem sinh_grad_formula {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ)
    (s θ : Fin N → ℝ) :
    deriv (fun t : ℝ => sinhF N ξ Ω (t • s) (t • θ)) 0
      = -(((N : ℝ) - 1) * (1 + ξ) / (2 * (1 - ξ)) + Ω * kappa0 ξ) * ∑ j : Fin N, s j := by
  obtain ⟨f1, H, hf⟩ := hasD2_sinhF hN h0 h1 Ω s θ
  rw [H.deriv1, hf, sum_lPair hN (le_refl 0) zero_lt_one s θ, zero_pow hN.ne']
  have e2 : ∑ k : Fin N, -(2 * s k * ξ) / (1 - ξ) = (-(2 * ξ) / (1 - ξ)) * ∑ k : Fin N, s k := by
    rw [Finset.mul_sum]; refine Finset.sum_congr rfl fun k _ => ?_; ring
  rw [e2, ← Finset.mul_sum]
  have hne : (1 : ℝ) - ξ ≠ 0 := by linarith
  field_simp
  ring

theorem sinhOmega_kappa0 (N : ℕ) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) :
    sinhOmega N ξ * kappa0 ξ = -(((N : ℝ) - 1) * (1 + ξ) / (2 * (1 - ξ))) := by
  have hne : (1 : ℝ) - ξ ≠ 0 := by linarith
  unfold sinhOmega kappa0
  field_simp
  ring

/-- **Note §2 (equilibrium condition).** With `Ω = sinhOmega N ξ = −(N−1)(1+ξ)(1−ξ)/(8ξ)` the ring
`s = θ = 0` is a critical point of `F = H − ΩJ`: along every line `t ↦ (ts, tθ)` the derivative at `t = 0`
is `0`. -/
theorem sinh_ring_critical {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1)
    (s θ : Fin N → ℝ) :
    HasDerivAt (fun t : ℝ => sinhF N ξ (sinhOmega N ξ) (t • s) (t • θ)) 0 0 := by
  obtain ⟨f1, H, _⟩ := hasD2_sinhF hN h0 h1 (sinhOmega N ξ) s θ
  have hd := H.1.self_of_nhds
  have := sinh_grad_formula hN h0 h1 (sinhOmega N ξ) s θ
  rw [hd.deriv, sinhOmega_kappa0 N h0 h1, add_neg_cancel, neg_zero, zero_mul] at this
  rw [this] at hd
  exact hd

/-- Uniqueness of `Ω`: the ring is critical iff `Ω = sinhOmega N ξ`. -/
theorem sinh_ring_critical_iff {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ) :
    (∀ s θ : Fin N → ℝ, deriv (fun t : ℝ => sinhF N ξ Ω (t • s) (t • θ)) 0 = 0)
      ↔ Ω = sinhOmega N ξ := by
  have hk : kappa0 ξ ≠ 0 := by
    unfold kappa0
    have : (1 : ℝ) - ξ ≠ 0 := by linarith
    positivity
  have hO := sinhOmega_kappa0 N h0 h1
  constructor
  · intro h
    have := h (fun _ => 1) 0
    rw [sinh_grad_formula hN h0 h1] at this
    simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one] at this
    have hN' : (N : ℝ) ≠ 0 := by exact_mod_cast hN.ne'
    rcases mul_eq_zero.mp this with h' | h'
    · apply mul_right_cancel₀ hk
      linarith
    · exact absurd h' hN'
  · rintro rfl s θ
    rw [sinh_grad_formula hN h0 h1, hO, add_neg_cancel, neg_zero, zero_mul]

/-! ## The Hessian -/

/-- The one-body diagonal of the Hessian: `W''(0) − ΩJ₁''(0)` with
`W(s) = ((N−1)/2) ln(1−ξe^{2s})`, `J₁(s) = (1+ξe^{2s})/(1−ξe^{2s})`. -/
def sinhD (N : ℕ) (ξ Ω : ℝ) : ℝ :=
  -(2 * ((N : ℝ) - 1) * ξ / (1 - ξ) ^ 2) - Ω * (8 * ξ * (1 + ξ) / (1 - ξ) ^ 3)

/-- Note §2: at `Ω = sinhOmega N ξ`, `W''(0) − ΩJ₁''(0) = C₁(ξ) = (N−1)(1+c²)/2`. -/
theorem sinhD_omega (N : ℕ) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) :
    sinhD N ξ (sinhOmega N ξ) = sinhC1 N ξ
      ∧ sinhD N ξ (sinhOmega N ξ) = ((N : ℝ) - 1) * (1 + cR ξ ^ 2) / 2 := by
  have hne : (1 : ℝ) - ξ ≠ 0 := by linarith
  have hξ : ξ ≠ 0 := h0.ne'
  constructor
  · unfold sinhD sinhOmega sinhC1
    field_simp
    ring
  · unfold sinhD sinhOmega cR
    field_simp
    ring

/-- `∂²F/∂s_a∂s_b` in closed form: the tanh pair-only matrix plus the one-body diagonal. -/
def sSS (N : ℕ) (ξ Ω : ℝ) (a b : Fin N) : ℝ := hSS N 0 0 a b + if a = b then sinhD N ξ Ω else 0

theorem sSS_symm (N : ℕ) (ξ Ω : ℝ) (a b : Fin N) : sSS N ξ Ω a b = sSS N ξ Ω b a := by
  unfold sSS
  rw [hSS_symm N 0 0 a b]
  by_cases h : a = b
  · subst h; rfl
  · simp [h, Ne.symm h]

theorem hST_zero (N : ℕ) (a b : Fin N) : hST N 0 a b = 0 := by
  unfold hST gammaM
  simp [rho_zero]

/-- The sinh Hessian quadratic form `d²/dt² F(ts,tθ)|_{t=0}`. -/
def sinhHessQ (N : ℕ) (ξ Ω : ℝ) (s θ : Fin N → ℝ) : ℝ :=
  deriv (deriv (fun t : ℝ => sinhF N ξ Ω (t • s) (t • θ))) 0

theorem qForm_pair (N : ℕ) (s θ : Fin N → ℝ) :
    qForm (hSS N 0 0) (hST N 0) (hTT N 0) s θ
      = ∑ j : Fin N, ∑ k : Fin N, (if j < k then qPair N 0 s θ j k else 0) := by
  simp_rw [qPair_eq]
  rw [sum_lt_add_swap (fun j k => R1 N 0 s θ j k + R2 N 0 s θ j k)]
  have hq : qForm (hSS N 0 0) (hST N 0) (hTT N 0) s θ = ∑ a : Fin N, ∑ b : Fin N,
      (s a * hSS N 0 0 a b * s b + 2 * (s a * hST N 0 a b * θ b) + θ a * hTT N 0 a b * θ b) := by
    unfold qForm bil
    simp only [Finset.mul_sum, ← Finset.sum_add_distrib]
  rw [hq]
  simp_rw [row_eq]
  simp [Vj]

/-- **Hessian quadratic form in matrix form.** -/
theorem sinhHessQ_eq {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ)
    (s θ : Fin N → ℝ) :
    sinhHessQ N ξ Ω s θ = qForm (sSS N ξ Ω) (hST N 0) (hTT N 0) s θ := by
  obtain ⟨f1, H, _⟩ := hasD2_sinhF hN h0 h1 Ω s θ
  unfold sinhHessQ
  rw [H.deriv2]
  have hbil : bil (sSS N ξ Ω) s s = bil (hSS N 0 0) s s + ∑ a : Fin N, sinhD N ξ Ω * s a ^ 2 := by
    unfold bil sSS
    simp only [mul_add, add_mul, Finset.sum_add_distrib]
    congr 1
    refine Finset.sum_congr rfl fun a _ => ?_
    simp only [mul_ite, ite_mul, mul_zero, zero_mul, Finset.sum_ite_eq, Finset.mem_univ, ite_true]
    ring
  have hq : qForm (sSS N ξ Ω) (hST N 0) (hTT N 0) s θ
      = qForm (hSS N 0 0) (hST N 0) (hTT N 0) s θ + ∑ a : Fin N, sinhD N ξ Ω * s a ^ 2 := by
    unfold qForm; rw [hbil]; ring
  rw [hq, qForm_pair]
  have hD : ((N : ℝ) - 1) / 2 * ∑ k : Fin N, -(4 * s k ^ 2 * ξ / (1 - ξ) ^ 2)
      - Ω * ∑ k : Fin N, 8 * s k ^ 2 * ξ * (1 + ξ) / (1 - ξ) ^ 3
      = ∑ a : Fin N, sinhD N ξ Ω * s a ^ 2 := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    unfold sinhD; ring
  rw [← hD]
  ring

/-! ## Polarisation: the Hessian entries -/

/-- Polarisation of `sinhHessQ`: the Hessian bilinear form. -/
def sinhHessB (N : ℕ) (ξ Ω : ℝ) (s θ s' θ' : Fin N → ℝ) : ℝ :=
  (sinhHessQ N ξ Ω (s + s') (θ + θ') - sinhHessQ N ξ Ω (s - s') (θ - θ')) / 4

/-- `∂²F/∂s_a∂s_b` at the ring. -/
def sinhHessSS (N : ℕ) (ξ Ω : ℝ) (a b : Fin N) : ℝ :=
  sinhHessB N ξ Ω (Pi.single a 1) 0 (Pi.single b 1) 0
/-- `∂²F/∂θ_a∂θ_b` at the ring. -/
def sinhHessTT (N : ℕ) (ξ Ω : ℝ) (a b : Fin N) : ℝ :=
  sinhHessB N ξ Ω 0 (Pi.single a 1) 0 (Pi.single b 1)
/-- `∂²F/∂s_a∂θ_b` at the ring. -/
def sinhHessST (N : ℕ) (ξ Ω : ℝ) (a b : Fin N) : ℝ :=
  sinhHessB N ξ Ω (Pi.single a 1) 0 0 (Pi.single b 1)
/-- `∂²F/∂θ_a∂s_b` at the ring. -/
def sinhHessTS (N : ℕ) (ξ Ω : ℝ) (a b : Fin N) : ℝ :=
  sinhHessB N ξ Ω 0 (Pi.single a 1) (Pi.single b 1) 0

section entries
variable {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ) (a b : Fin N)
include hN h0 h1

theorem sinhHessSS_eq : sinhHessSS N ξ Ω a b = sSS N ξ Ω a b := by
  unfold sinhHessSS sinhHessB
  rw [sinhHessQ_eq hN h0 h1, sinhHessQ_eq hN h0 h1]
  unfold qForm
  simp only [add_zero, sub_zero, bil_add_left, bil_add_right, bil_sub_left, bil_sub_right,
    bil_zero_right, bil_single]
  rw [sSS_symm N ξ Ω b a]; ring

theorem sinhHessTT_eq : sinhHessTT N ξ Ω a b = hTT N 0 a b := by
  unfold sinhHessTT sinhHessB
  rw [sinhHessQ_eq hN h0 h1, sinhHessQ_eq hN h0 h1]
  unfold qForm
  simp only [add_zero, zero_add, sub_self, bil_add_left, bil_add_right, bil_sub_left,
    bil_sub_right, bil_zero_left, bil_zero_right, bil_single]
  rw [hTT_symm N 0 b a]; ring

theorem sinhHessST_eq : sinhHessST N ξ Ω a b = 0 := by
  unfold sinhHessST sinhHessB
  rw [sinhHessQ_eq hN h0 h1, sinhHessQ_eq hN h0 h1]
  unfold qForm
  simp only [add_zero, sub_zero, zero_add, bil_sub_left,
    bil_sub_right, bil_zero_left, bil_zero_right, bil_single]
  rw [hST_zero]; ring

theorem sinhHessTS_eq : sinhHessTS N ξ Ω a b = 0 := by
  unfold sinhHessTS sinhHessB
  rw [sinhHessQ_eq hN h0 h1, sinhHessQ_eq hN h0 h1]
  unfold qForm
  simp only [add_zero, sub_zero, zero_add, zero_sub, bil_neg_left, bil_neg_right, bil_single]
  rw [hST_zero]; ring

end entries

/-! ## The Fourier blocks (note Prop. `prop:modes`) -/

theorem sBeta_zero (N m : ℕ) : sBeta N 0 m = 0 := by
  unfold sBeta; simp [rho_zero]

/-- **Note Prop. `prop:modes`.** For `0 < ξ < 1`, `Ω = sinhOmega N ξ`, `0 ≤ m ≤ N` and
`v_b = (w^b)^m = w^{mb}`, the Hessian of `F = H − ΩJ` at the ring satisfies
`Σ_b ∂²F/∂s_a∂s_b v_b = a_m v_a`, `Σ_b ∂²F/∂θ_a∂θ_b v_b = b_m v_a`, and the mixed blocks vanish
(`Σ_b ∂²F/∂s_a∂θ_b v_b = Σ_b ∂²F/∂θ_a∂s_b v_b = 0`, i.e. `g_m = 0`), with
`a_m = (N−1)(1+c²)/2 − m(N−m)/2`, `b_m = m(N−m)/2`, `c = (1+ξ)/(1−ξ)`. -/
theorem sinh_blocks_eigen {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (m : ℕ)
    (hmN : m ≤ N) (a : Fin N) :
    ∑ b : Fin N, ((sinhHessSS N ξ (sinhOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
        = (sinhA N m ξ : ℂ) * (ringW N ^ (a : ℕ)) ^ m
      ∧ ∑ b : Fin N, ((sinhHessTT N ξ (sinhOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
        = (sinhB N m : ℂ) * (ringW N ^ (a : ℕ)) ^ m
      ∧ ∑ b : Fin N, ((sinhHessST N ξ (sinhOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m = 0
      ∧ ∑ b : Fin N, ((sinhHessTS N ξ (sinhOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
        = 0 := by
  have : NeZero N := ⟨hN.ne'⟩
  have hsA := sAlpha_eq (N := N) m hmN
  refine ⟨?_, ?_, ?_, ?_⟩
  · simp_rw [sinhHessSS_eq hN h0 h1]
    unfold sSS
    simp only [Complex.ofReal_add, add_mul, Finset.sum_add_distrib]
    rw [eig_hSS, hsA, sBeta_zero, sBeta_zero]
    have hdiag : ∑ b : Fin N, (((if a = b then sinhD N ξ (sinhOmega N ξ) else 0 : ℝ)) : ℂ)
        * (ringW N ^ (b : ℕ)) ^ m = (sinhD N ξ (sinhOmega N ξ) : ℂ) * (ringW N ^ (a : ℕ)) ^ m := by
      rw [Finset.sum_eq_single a]
      · simp
      · intro b _ hb; simp [Ne.symm hb]
      · intro h; exact absurd (Finset.mem_univ a) h
    rw [hdiag, (sinhD_omega N h0 h1).2]
    unfold sinhA Vj
    push_cast
    ring
  · simp_rw [sinhHessTT_eq hN h0 h1]
    rw [eig_hTT, hsA, sBeta_zero, sBeta_zero]
    unfold sinhB
    push_cast
    ring
  · simp_rw [sinhHessST_eq hN h0 h1]
    simp
  · simp_rw [sinhHessTS_eq hN h0 h1]
    simp

/-- **Note Prop. `prop:modes`, Hermitian form.** For `0 ≤ m, m' ≤ N−1`, `v_b = w^{mb}`, `v'_a = w^{m'a}`:
`Σ_{a,b} conj(v'_a) ∂²F/∂s_a∂s_b v_b = N a_m [m = m']`, `Σ_{a,b} conj(v'_a) ∂²F/∂θ_a∂θ_b v_b = N b_m [m = m']`,
and the mixed `sθ`, `θs` forms vanish: the Hessian is block-diagonal over Fourier modes, `diag(a_m, b_m)`
on mode `m`. -/
theorem sinh_blocks_hermitian {N : ℕ} {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (m m' : ℕ)
    (hmN : m < N) (hm' : m' < N) :
    (∑ a : Fin N, ∑ b : Fin N, conj ((ringW N ^ (a : ℕ)) ^ m')
        * ((sinhHessSS N ξ (sinhOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
      = if m = m' then (N : ℂ) * (sinhA N m ξ : ℂ) else 0)
    ∧ (∑ a : Fin N, ∑ b : Fin N, conj ((ringW N ^ (a : ℕ)) ^ m')
        * ((sinhHessTT N ξ (sinhOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m
      = if m = m' then (N : ℂ) * (sinhB N m : ℂ) else 0)
    ∧ (∑ a : Fin N, ∑ b : Fin N, conj ((ringW N ^ (a : ℕ)) ^ m')
        * ((sinhHessST N ξ (sinhOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m = 0)
    ∧ (∑ a : Fin N, ∑ b : Fin N, conj ((ringW N ^ (a : ℕ)) ^ m')
        * ((sinhHessTS N ξ (sinhOmega N ξ) a b : ℝ) : ℂ) * (ringW N ^ (b : ℕ)) ^ m = 0) := by
  have hN : 0 < N := by omega
  have E := sinh_blocks_eigen hN h0 h1 m hmN.le
  refine ⟨herm_of_eig hN (fun a => (E a).1) m' hmN hm',
    herm_of_eig hN (fun a => (E a).2.1) m' hmN hm', ?_, ?_⟩
  · have := herm_of_eig hN (μ := 0) (fun a => by rw [(E a).2.2.1, zero_mul]) m' hmN hm'
    rw [this]; simp
  · have := herm_of_eig hN (μ := 0) (fun a => by rw [(E a).2.2.2, zero_mul]) m' hmN hm'
    rw [this]; simp

end Havelock.Sinh
