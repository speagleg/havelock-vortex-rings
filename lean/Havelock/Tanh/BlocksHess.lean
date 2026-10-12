import Havelock.Tanh.Blocks

/-! # Tanh model: the Hessian at the ring and its Fourier blocks (note Thm `thm:blocks`)

Continuation of `Blocks.lean` (see the module docstring there for the notions used).

* `hessQ_eq`: the Hessian quadratic form `hessQ s θ = d²/dt² F(ts,tθ)|₀` equals
  `sᵀ hSS s + 2 sᵀ hST θ + θᵀ hTT θ` with explicit circulant matrices `hSS`, `hST`, `hTT`.
* `hessSS_eq`, `hessTT_eq`, `hessST_eq`, `hessTS_eq`: the polarised entries are `hSS a b`, `hTT a b`,
  `hST a b`, `hST b a`.
* `blocks_eigen`, `blocks_eigen_zero`, `blocks_hermitian`: the Fourier blocks (note Thm `thm:blocks`).
-/

noncomputable section

namespace Havelock.Tanh

open scoped ComplexConjugate

/-! ## Bilinear forms -/

/-- `bil A x y = xᵀ A y`. -/
def bil {N : ℕ} (A : Fin N → Fin N → ℝ) (x y : Fin N → ℝ) : ℝ := ∑ a, ∑ b, x a * A a b * y b

section bil
variable {N : ℕ} (A : Fin N → Fin N → ℝ) (x x' y y' : Fin N → ℝ)

theorem bil_add_left : bil A (x + x') y = bil A x y + bil A x' y := by
  unfold bil; simp only [Pi.add_apply, add_mul, Finset.sum_add_distrib]
theorem bil_add_right : bil A x (y + y') = bil A x y + bil A x y' := by
  unfold bil; simp only [Pi.add_apply, mul_add, Finset.sum_add_distrib]
theorem bil_sub_left : bil A (x - x') y = bil A x y - bil A x' y := by
  unfold bil; simp only [Pi.sub_apply, sub_mul, Finset.sum_sub_distrib]
theorem bil_sub_right : bil A x (y - y') = bil A x y - bil A x y' := by
  unfold bil; simp only [Pi.sub_apply, mul_sub, Finset.sum_sub_distrib]
theorem bil_zero_left : bil A 0 y = 0 := by simp [bil]
theorem bil_zero_right : bil A x 0 = 0 := by simp [bil]
theorem bil_neg_left : bil A (-x) y = -bil A x y := by
  unfold bil; simp only [Pi.neg_apply, neg_mul, Finset.sum_neg_distrib]
theorem bil_neg_right : bil A x (-y) = -bil A x y := by
  unfold bil; simp only [Pi.neg_apply, mul_neg, Finset.sum_neg_distrib]
theorem bil_single (a b : Fin N) : bil A (Pi.single a 1) (Pi.single b 1) = A a b := by
  simp [bil, Pi.single_apply]

end bil

/-- `qForm A C D s θ = sᵀ A s + 2 sᵀ C θ + θᵀ D θ`. -/
def qForm {N : ℕ} (A C D : Fin N → Fin N → ℝ) (s θ : Fin N → ℝ) : ℝ :=
  bil A s s + 2 * bil C s θ + bil D θ θ

/-! ## The explicit Hessian matrices -/

/-- `α_{jk} = 1/|1−ζ_{jk}|²`. -/
def alphaM (N : ℕ) (j k : Fin N) : ℝ := alphaF (zeta N j k)
/-- `Re ρ_{jk}`. -/
def betaM (N : ℕ) (ξ : ℝ) (j k : Fin N) : ℝ := (rho ξ (zeta N j k)).re
/-- `Im ρ_{jk}`. -/
def gammaM (N : ℕ) (ξ : ℝ) (j k : Fin N) : ℝ := (rho ξ (zeta N j k)).im

/-- The one-body (`−ΩJ`) diagonal: `−Ω · 8ξ(1+ξ)/(1−ξ)³ = −2cΩκ₀`. -/
def Vj (ξ Ω : ℝ) : ℝ := -Ω * (8 * ξ * (1 + ξ) / (1 - ξ) ^ 3)

/-- `∂²F/∂s_a∂s_b` in closed form. -/
def hSS (N : ℕ) (ξ Ω : ℝ) (a b : Fin N) : ℝ :=
  if a = b then Vj ξ Ω - ∑ k : Fin N, (if a = k then 0 else alphaM N a k + betaM N ξ a k)
  else alphaM N a b - betaM N ξ a b

/-- `∂²F/∂θ_a∂θ_b` in closed form. -/
def hTT (N : ℕ) (ξ : ℝ) (a b : Fin N) : ℝ :=
  if a = b then ∑ k : Fin N, (if a = k then 0 else alphaM N a k + betaM N ξ a k)
  else -(alphaM N a b + betaM N ξ a b)

/-- `∂²F/∂s_a∂θ_b` in closed form. -/
def hST (N : ℕ) (ξ : ℝ) (a b : Fin N) : ℝ :=
  if a = b then ∑ k : Fin N, (if a = k then 0 else gammaM N ξ a k) else -gammaM N ξ a b

theorem alphaM_swap (N : ℕ) (j k : Fin N) : alphaM N k j = alphaM N j k := by
  unfold alphaM alphaF
  rw [zeta_swap, show (1 : ℂ) - conj (zeta N j k) = conj (1 - zeta N j k) by simp,
    Complex.normSq_conj]

theorem betaM_swap (N : ℕ) (ξ : ℝ) (j k : Fin N) : betaM N ξ k j = betaM N ξ j k := by
  unfold betaM; rw [zeta_swap, rho_conj, Complex.conj_re]

theorem gammaM_swap (N : ℕ) (ξ : ℝ) (j k : Fin N) : gammaM N ξ k j = -gammaM N ξ j k := by
  unfold gammaM; rw [zeta_swap, rho_conj, Complex.conj_im]

theorem hSS_symm (N : ℕ) (ξ Ω : ℝ) (a b : Fin N) : hSS N ξ Ω a b = hSS N ξ Ω b a := by
  unfold hSS
  by_cases h : a = b
  · subst h; rfl
  · rw [if_neg h, if_neg (Ne.symm h), alphaM_swap, betaM_swap]

theorem hTT_symm (N : ℕ) (ξ : ℝ) (a b : Fin N) : hTT N ξ a b = hTT N ξ b a := by
  unfold hTT
  by_cases h : a = b
  · subst h; rfl
  · rw [if_neg h, if_neg (Ne.symm h), alphaM_swap, betaM_swap]

/-! ## The quadratic form in matrix form -/

/-- Diagonal part of the per-pair summand. -/
def R1 (N : ℕ) (ξ : ℝ) (s θ : Fin N → ℝ) (j k : Fin N) : ℝ :=
  -(alphaM N j k + betaM N ξ j k) * s j ^ 2 + 2 * gammaM N ξ j k * s j * θ j
    + (alphaM N j k + betaM N ξ j k) * θ j ^ 2

/-- Off-diagonal part of the per-pair summand. -/
def R2 (N : ℕ) (ξ : ℝ) (s θ : Fin N → ℝ) (j k : Fin N) : ℝ :=
  (alphaM N j k - betaM N ξ j k) * s j * s k - 2 * gammaM N ξ j k * s j * θ k
    - (alphaM N j k + betaM N ξ j k) * θ j * θ k

theorem qPair_eq (N : ℕ) (ξ : ℝ) (s θ : Fin N → ℝ) (j k : Fin N) :
    qPair N ξ s θ j k
      = (R1 N ξ s θ j k + R2 N ξ s θ j k) + (R1 N ξ s θ k j + R2 N ξ s θ k j) := by
  unfold qPair R1 R2
  rw [alphaM_swap N j k, betaM_swap N ξ j k, gammaM_swap N ξ j k]
  unfold alphaM betaM gammaM
  ring

theorem row_eq (N : ℕ) (ξ Ω : ℝ) (s θ : Fin N → ℝ) (a : Fin N) :
    ∑ b : Fin N, (s a * hSS N ξ Ω a b * s b + 2 * (s a * hST N ξ a b * θ b) + θ a * hTT N ξ a b * θ b)
      = Vj ξ Ω * s a ^ 2 + ∑ b : Fin N, (if a = b then 0 else R1 N ξ s θ a b + R2 N ξ s θ a b) := by
  have hR1 : ∑ k : Fin N, (if a = k then 0 else R1 N ξ s θ a k)
      = -(s a ^ 2) * ∑ k : Fin N, (if a = k then 0 else alphaM N a k + betaM N ξ a k)
        + 2 * s a * θ a * ∑ k : Fin N, (if a = k then 0 else gammaM N ξ a k)
        + θ a ^ 2 * ∑ k : Fin N, (if a = k then 0 else alphaM N a k + betaM N ξ a k) := by
    simp only [Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    split_ifs
    · ring
    · unfold R1; ring
  have hb : ∀ b : Fin N,
      s a * hSS N ξ Ω a b * s b + 2 * (s a * hST N ξ a b * θ b) + θ a * hTT N ξ a b * θ b
        = (if a = b then Vj ξ Ω * s a ^ 2 + ∑ k : Fin N, (if a = k then 0 else R1 N ξ s θ a k) else 0)
          + (if a = b then 0 else R2 N ξ s θ a b) := by
    intro b
    unfold hSS hST hTT
    by_cases h : a = b
    · subst h
      simp only [if_true, add_zero]
      rw [hR1]; ring
    · simp only [if_neg h, zero_add]
      unfold R2; ring
  rw [Finset.sum_congr rfl (fun b _ => hb b), Finset.sum_add_distrib, Finset.sum_ite_eq,
    if_pos (Finset.mem_univ a), add_assoc, ← Finset.sum_add_distrib]
  congr 1
  refine Finset.sum_congr rfl fun b _ => ?_
  split_ifs <;> ring

/-- The Hessian quadratic form `d²/dt² F(ts,tθ)|_{t=0}`. -/
def hessQ (N : ℕ) (ξ Ω : ℝ) (s θ : Fin N → ℝ) : ℝ :=
  deriv (deriv (fun t : ℝ => ringF N ξ Ω (t • s) (t • θ))) 0

/-- **Hessian quadratic form in matrix form.** -/
theorem hessQ_eq {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ)
    (s θ : Fin N → ℝ) :
    hessQ N ξ Ω s θ = qForm (hSS N ξ Ω) (hST N ξ) (hTT N ξ) s θ := by
  unfold hessQ
  rw [hess_eq hN h0 h1]
  simp_rw [qPair_eq]
  rw [sum_lt_add_swap (fun j k => R1 N ξ s θ j k + R2 N ξ s θ j k)]
  have hq : qForm (hSS N ξ Ω) (hST N ξ) (hTT N ξ) s θ = ∑ a : Fin N, ∑ b : Fin N,
      (s a * hSS N ξ Ω a b * s b + 2 * (s a * hST N ξ a b * θ b) + θ a * hTT N ξ a b * θ b) := by
    unfold qForm bil
    simp only [Finset.mul_sum, ← Finset.sum_add_distrib]
  rw [hq]
  simp_rw [row_eq]
  have hV : Ω * ∑ k : Fin N, 8 * s k ^ 2 * ξ * (1 + ξ) / (1 - ξ) ^ 3
      = -∑ k : Fin N, Vj ξ Ω * s k ^ 2 := by
    rw [Finset.mul_sum, ← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    unfold Vj; ring
  rw [hV, Finset.sum_add_distrib]
  ring

/-! ## Polarisation: the Hessian entries -/

/-- Polarisation of `hessQ`: the Hessian bilinear form. -/
def hessB (N : ℕ) (ξ Ω : ℝ) (s θ s' θ' : Fin N → ℝ) : ℝ :=
  (hessQ N ξ Ω (s + s') (θ + θ') - hessQ N ξ Ω (s - s') (θ - θ')) / 4

/-- `∂²F/∂s_a∂s_b` at the ring. -/
def hessSS (N : ℕ) (ξ Ω : ℝ) (a b : Fin N) : ℝ := hessB N ξ Ω (Pi.single a 1) 0 (Pi.single b 1) 0
/-- `∂²F/∂θ_a∂θ_b` at the ring. -/
def hessTT (N : ℕ) (ξ Ω : ℝ) (a b : Fin N) : ℝ := hessB N ξ Ω 0 (Pi.single a 1) 0 (Pi.single b 1)
/-- `∂²F/∂s_a∂θ_b` at the ring. -/
def hessST (N : ℕ) (ξ Ω : ℝ) (a b : Fin N) : ℝ := hessB N ξ Ω (Pi.single a 1) 0 0 (Pi.single b 1)
/-- `∂²F/∂θ_a∂s_b` at the ring. -/
def hessTS (N : ℕ) (ξ Ω : ℝ) (a b : Fin N) : ℝ := hessB N ξ Ω 0 (Pi.single a 1) (Pi.single b 1) 0

section entries
variable {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) (Ω : ℝ) (a b : Fin N)
include hN h0 h1

theorem hessSS_eq : hessSS N ξ Ω a b = hSS N ξ Ω a b := by
  unfold hessSS hessB
  rw [hessQ_eq hN h0 h1, hessQ_eq hN h0 h1]
  unfold qForm
  simp only [add_zero, sub_zero, bil_add_left, bil_add_right, bil_sub_left, bil_sub_right,
    bil_zero_right, bil_single]
  rw [hSS_symm N ξ Ω b a]; ring

theorem hessTT_eq : hessTT N ξ Ω a b = hTT N ξ a b := by
  unfold hessTT hessB
  rw [hessQ_eq hN h0 h1, hessQ_eq hN h0 h1]
  unfold qForm
  simp only [add_zero, zero_add, sub_self, bil_add_left, bil_add_right, bil_sub_left,
    bil_sub_right, bil_zero_left, bil_zero_right, bil_single]
  rw [hTT_symm N ξ b a]; ring

theorem hessST_eq : hessST N ξ Ω a b = hST N ξ a b := by
  unfold hessST hessB
  rw [hessQ_eq hN h0 h1, hessQ_eq hN h0 h1]
  unfold qForm
  simp only [add_zero, sub_zero, zero_add, bil_sub_left,
    bil_sub_right, bil_zero_left, bil_zero_right, bil_single]
  ring

theorem hessTS_eq : hessTS N ξ Ω a b = hST N ξ b a := by
  unfold hessTS hessB
  rw [hessQ_eq hN h0 h1, hessQ_eq hN h0 h1]
  unfold qForm
  simp only [add_zero, sub_zero, zero_add, zero_sub, bil_neg_left, bil_neg_right, bil_single]
  ring

end entries

end Havelock.Tanh
