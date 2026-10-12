import Havelock.Tanh.BlocksFourier
import Havelock.Tanh.LargeN.Dynamics


/-! # Tanh model: the linearised flow of the full `2N`-dimensional system (note Prop. `prop:lin`)

At the ring (`0 < ξ < 1`, `Ω = ringOmega N ξ`), with the Hessian `H = [[H_ss, H_sθ], [H_θs, H_θθ]]` of
`F = H − ΩJ` in the coordinates `(s, θ)` (entries `hessSS`, `hessST`, `hessTS`, `hessTT` of
`BlocksHess.lean`, i.e. polarised second derivatives of `t ↦ F(ts, tθ)`), the symplectic form
`κ₀ Σ ds_j ∧ dθ_j` (`κ₀ = 4ξ/(1−ξ)²`) and `J = [[0, I_N], [−I_N, 0]]`, the linearised flow is
`linFlow N ξ = (1/κ₀) J H` (complex `2N × 2N` matrix; sign convention `κ₀ (ṡ, θ̇) = J ∇F`, as in the
note's `κ₀ (Ṡ, Ṫ)ᵀ = J₀ M_m (S, T)ᵀ`).

* `linFlow_mode` (`1 ≤ m ≤ N−1`): `L (S v, T v) = (1/κ₀) (J₀ M_m (S,T)ᵀ) ⊗ v`, `v_b = w^{mb}`
  (`ringFlow N m ξ = J₀ M_m = flowBlock a_m b_m g_m`); `linFlow_modeVec` is the same for every
  `0 ≤ m ≤ N−1` with the mode-`0` block `[[c₀, 0], [0, 0]]` (`modeBlock`).
* `linFlow_hasEigenvalue_iff`: `σ` is an eigenvalue of `L` iff it is an eigenvalue of
  `(1/κ₀) J₀ M_m` for some `0 ≤ m ≤ N−1` (mode `0`: `b_0 = g_0 = 0`, eigenvalue `0` only).
  The proof uses the discrete Fourier basis (`fourier_inv`, `eq_zero_of_fhat`) and the block
  diagonalisation `fhat_linFlow`, `fhat_hessMat`.
* `linFlow_spectral_iff`: all eigenvalues of `L` have zero real part iff `specStable N ξ`.
* `hessMat_parseval`: `x* H x = N⁻¹ Σ_m x̂_m* M_m x̂_m`.
* `hessMat_posDef_compl_iff`: `x* H x > 0` for every nonzero `x ∈ ℂ^{2N}` with `Σ_b x_{s,b} = 0 = Σ_b x_{θ,b}`
  (orthogonal complement of mode `0`) iff `emDefinite N ξ`; `hessQ_posDef_compl_iff`: the same for the
  real second derivative `hessQ` on real `(s, θ) ≠ 0` with `Σ s = Σ θ = 0`.

Not formalised: that `L` is the derivative of the vortex vector field (the Hamiltonian dynamics itself);
`L` is defined as `(1/κ₀) J H` from the Hessian.
-/

noncomputable section

namespace Havelock.Tanh.LargeN

open Matrix ComplexOrder
open scoped ComplexConjugate
open Complex (I)

/-! ### Discrete Fourier modes on `ℂ^{2N} = (Fin N ⊕ Fin N → ℂ)` -/

/-- The Fourier vector of mode `m`: `v_b = w^{mb} = (w^b)^m`, `w = e^{2πi/N}`. -/
def fv (N : ℕ) (m b : Fin N) : ℂ := (ringW N ^ (b : ℕ)) ^ (m : ℕ)

/-- The mode vector `(S v, T v)` for `z = (S, T)` (first block `s`, second block `θ`). -/
def modeVec {N : ℕ} (m : Fin N) (z : Fin 2 → ℂ) : Fin N ⊕ Fin N → ℂ :=
  Sum.elim (fun b => z 0 * fv N m b) (fun b => z 1 * fv N m b)

/-- The Fourier coefficient of `x` on mode `m`: `(Σ_b conj(v_b) x_{s,b}, Σ_b conj(v_b) x_{θ,b})`. -/
def fhat {N : ℕ} (m : Fin N) (x : Fin N ⊕ Fin N → ℂ) : Fin 2 → ℂ :=
  ![∑ b, conj (fv N m b) * x (Sum.inl b), ∑ b, conj (fv N m b) * x (Sum.inr b)]

theorem fv_zero_left {N : ℕ} (b : Fin N) (m : Fin N) (hm : (m : ℕ) = 0) : fv N m b = 1 := by
  simp [fv, hm]

theorem fv_zero_right {N : ℕ} (m : Fin N) (b : Fin N) (hb : (b : ℕ) = 0) : fv N m b = 1 := by
  simp [fv, hb]

/-- Orthogonality of the characters (sum over the modes). -/
theorem fv_orth {N : ℕ} (a b : Fin N) :
    ∑ m : Fin N, fv N m a * conj (fv N m b) = if a = b then (N : ℂ) else 0 := by
  have hN : 0 < N := a.pos
  have e : ∀ m : Fin N, fv N m a * conj (fv N m b)
      = (ringW N ^ (m : ℕ)) ^ (a : ℕ) * conj (ringW N ^ (m : ℕ)) ^ (b : ℕ) := by
    intro m
    rw [fv, fv, ← pow_mul, ← pow_mul, mul_comm (a : ℕ), mul_comm (b : ℕ), pow_mul, pow_mul, map_pow]
  rw [Finset.sum_congr rfl (fun m _ => e m), rootSum hN, Nat.mod_eq_of_lt a.2, Nat.mod_eq_of_lt b.2]
  simp only [Fin.val_inj]

/-- Orthogonality of the characters (sum over the sites), diagonal case. -/
theorem fv_norm {N : ℕ} (m : Fin N) : ∑ b : Fin N, fv N m b * conj (fv N m b) = N := by
  have e : ∀ b : Fin N, fv N m b * conj (fv N m b) = 1 := by
    intro b; rw [fv, ← pow_mul]; exact ringW_pow_mul_conj N _
  simp [e]

/-- Characters of nonzero modes sum to zero. -/
theorem fv_sum {N : ℕ} (m : Fin N) (hm : (m : ℕ) ≠ 0) : ∑ b : Fin N, fv N m b = 0 := by
  have hN : 0 < N := m.pos
  have h := rootSum hN (m : ℕ) 0
  rw [if_neg (by rw [Nat.mod_eq_of_lt m.2, Nat.zero_mod]; exact hm)] at h
  rw [← h]
  refine Finset.sum_congr rfl fun b _ => ?_
  simp [fv]

theorem modeVec_zero {N : ℕ} (m : Fin N) : modeVec m 0 = 0 := by
  funext i; rcases i with b | b <;> simp [modeVec]

theorem modeVec_smul {N : ℕ} (m : Fin N) (c : ℂ) (z : Fin 2 → ℂ) :
    modeVec m (c • z) = c • modeVec m z := by
  funext i; rcases i with b | b <;> simp [modeVec, mul_assoc]

theorem modeVec_ne_zero {N : ℕ} (m : Fin N) {z : Fin 2 → ℂ} (hz : z ≠ 0) : modeVec m z ≠ 0 := by
  intro h
  have hN : 0 < N := m.pos
  have h0 := congrFun h (Sum.inl ⟨0, hN⟩)
  have h1 := congrFun h (Sum.inr ⟨0, hN⟩)
  simp only [modeVec, Sum.elim_inl, Sum.elim_inr, Pi.zero_apply,
    fv_zero_right m ⟨0, hN⟩ rfl, mul_one] at h0 h1
  apply hz; ext i; fin_cases i <;> simp [h0, h1]

theorem fhat_smul {N : ℕ} (m : Fin N) (c : ℂ) (x : Fin N ⊕ Fin N → ℂ) :
    fhat m (c • x) = c • fhat m x := by
  ext i; fin_cases i <;> simp [fhat, Finset.mul_sum, mul_left_comm]

/-- Fourier inversion: `x = N⁻¹ Σ_m (fhat_m x) ⊗ v_m`. -/
theorem fourier_inv {N : ℕ} (x : Fin N ⊕ Fin N → ℂ) :
    x = (N : ℂ)⁻¹ • ∑ m : Fin N, modeVec m (fhat m x) := by
  have key : ∀ (y : Fin N → ℂ) (a : Fin N),
      y a = (N : ℂ)⁻¹ * ∑ m : Fin N, (∑ b, conj (fv N m b) * y b) * fv N m a := by
    intro y a
    have hN : (N : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr a.pos.ne'
    have : ∑ m : Fin N, (∑ b, conj (fv N m b) * y b) * fv N m a
        = ∑ b, y b * ∑ m : Fin N, fv N m a * conj (fv N m b) := by
      simp_rw [Finset.sum_mul, Finset.mul_sum]
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl fun b _ => Finset.sum_congr rfl fun m _ => ?_
      ring
    rw [this]
    simp_rw [fv_orth]
    simp only [mul_ite, mul_zero, Finset.sum_ite_eq, Finset.mem_univ, if_true]
    rw [mul_comm (y a), ← mul_assoc, inv_mul_cancel₀ hN, one_mul]
  funext i
  rcases i with a | a
  · simp only [Pi.smul_apply, Finset.sum_apply, modeVec, Sum.elim_inl, smul_eq_mul, fhat]
    simpa using key (fun b => x (Sum.inl b)) a
  · simp only [Pi.smul_apply, Finset.sum_apply, modeVec, Sum.elim_inr, smul_eq_mul, fhat]
    simpa using key (fun b => x (Sum.inr b)) a

/-- A vector with all Fourier coefficients zero is zero. -/
theorem eq_zero_of_fhat {N : ℕ} {x : Fin N ⊕ Fin N → ℂ} (h : ∀ m : Fin N, fhat m x = 0) : x = 0 := by
  rw [fourier_inv x]
  simp [h, modeVec_zero]

theorem star_modeVec_dot {N : ℕ} (m : Fin N) (z : Fin 2 → ℂ) (y : Fin N ⊕ Fin N → ℂ) :
    star (modeVec m z) ⬝ᵥ y = star z ⬝ᵥ fhat m y := by
  simp only [dotProduct, Fintype.sum_sum_type, Fin.sum_univ_two, modeVec, fhat, Pi.star_apply,
    Sum.elim_inl, Sum.elim_inr, star_mul', Finset.mul_sum, Matrix.cons_val_zero,
    Matrix.cons_val_one]
  congr 1 <;> refine Finset.sum_congr rfl fun b _ => ?_ <;> simp only [RCLike.star_def] <;> ring

theorem dot_modeVec {N : ℕ} (m : Fin N) (z : Fin 2 → ℂ) (y : Fin N ⊕ Fin N → ℂ) :
    star y ⬝ᵥ modeVec m z = star (fhat m y) ⬝ᵥ z := by
  simp only [dotProduct, Fintype.sum_sum_type, Fin.sum_univ_two, modeVec, fhat, Pi.star_apply,
    Sum.elim_inl, Sum.elim_inr, Matrix.cons_val_zero,
    Matrix.cons_val_one, RCLike.star_def, map_sum, map_mul, Finset.sum_mul, Complex.conj_conj]
  congr 1 <;> refine Finset.sum_congr rfl fun b _ => ?_ <;> ring

/-- If `Kᴴ` maps mode `m` to itself by `Pᴴ`, then the Fourier coefficient of `K x` on mode `m`
is `P` applied to that of `x`. -/
theorem fhat_mulVec {N : ℕ} (m : Fin N) (K : Matrix (Fin N ⊕ Fin N) (Fin N ⊕ Fin N) ℂ)
    (P : Matrix (Fin 2) (Fin 2) ℂ) (h : ∀ z, Kᴴ *ᵥ modeVec m z = modeVec m (Pᴴ *ᵥ z))
    (x : Fin N ⊕ Fin N → ℂ) : fhat m (K *ᵥ x) = P *ᵥ fhat m x := by
  have key : ∀ z, star z ⬝ᵥ fhat m (K *ᵥ x) = star z ⬝ᵥ (P *ᵥ fhat m x) := by
    intro z
    have e1 : star (modeVec m z) ᵥ* K = star (Kᴴ *ᵥ modeVec m z) := by
      rw [star_mulVec, conjTranspose_conjTranspose]
    have e2 : star z ᵥ* P = star (Pᴴ *ᵥ z) := by
      rw [star_mulVec, conjTranspose_conjTranspose]
    rw [← star_modeVec_dot, dotProduct_mulVec, e1, h, star_modeVec_dot, dotProduct_mulVec, e2]
  ext i
  have := key (Pi.single i 1)
  simpa using this

/-! ### The eigenvalue API -/

theorem hasEigenvalue_iff_exists {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℂ) (σ : ℂ) :
    Module.End.HasEigenvalue (Matrix.toLin' A) σ ↔ ∃ x, x ≠ 0 ∧ A *ᵥ x = σ • x := by
  constructor
  · intro h
    obtain ⟨x, hx⟩ := h.exists_hasEigenvector
    rw [Module.End.hasEigenvector_iff, Module.End.mem_eigenspace_iff, Matrix.toLin'_apply] at hx
    exact ⟨x, hx.2, hx.1⟩
  · rintro ⟨x, hx0, hx⟩
    apply Module.End.hasEigenvalue_of_hasEigenvector (x := x)
    rw [Module.End.hasEigenvector_iff, Module.End.mem_eigenspace_iff, Matrix.toLin'_apply]
    exact ⟨hx, hx0⟩

theorem hasEigenvalue_smul_iff {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℂ) {c : ℂ}
    (hc : c ≠ 0) (σ : ℂ) :
    Module.End.HasEigenvalue (Matrix.toLin' (c • A)) σ ↔
      Module.End.HasEigenvalue (Matrix.toLin' A) (c⁻¹ * σ) := by
  rw [hasEigenvalue_iff_exists, hasEigenvalue_iff_exists]
  refine exists_congr fun x => and_congr_right fun _ => ?_
  rw [smul_mulVec, mul_smul]
  constructor
  · intro h; rw [← h, smul_smul, inv_mul_cancel₀ hc, one_smul]
  · intro h; rw [h, smul_smul, mul_inv_cancel₀ hc, one_smul]

/-! ### The Hessian, the symplectic matrix and the linearised flow of the ring -/

/-- The `2N × 2N` Hessian of `F = H − Ω J` at the ring (`Ω = ringOmega N ξ`) in the coordinates
`(s, θ)`, as a complex matrix: `[[∂²F/∂s∂s, ∂²F/∂s∂θ], [∂²F/∂θ∂s, ∂²F/∂θ∂θ]]` (entries from
`BlocksHess.lean`). -/
def hessMat (N : ℕ) (ξ : ℝ) : Matrix (Fin N ⊕ Fin N) (Fin N ⊕ Fin N) ℂ :=
  fromBlocks (of fun a b => ((hessSS N ξ (ringOmega N ξ) a b : ℝ) : ℂ))
    (of fun a b => ((hessST N ξ (ringOmega N ξ) a b : ℝ) : ℂ))
    (of fun a b => ((hessTS N ξ (ringOmega N ξ) a b : ℝ) : ℂ))
    (of fun a b => ((hessTT N ξ (ringOmega N ξ) a b : ℝ) : ℂ))

/-- The symplectic matrix `J = [[0, I_N], [−I_N, 0]]`. -/
def symplJ (N : ℕ) : Matrix (Fin N ⊕ Fin N) (Fin N ⊕ Fin N) ℂ := fromBlocks 0 1 (-1) 0

/-- The linearised flow of the full `2N`-dimensional system at the ring (note Prop. `prop:lin`):
`L = (1/κ₀) J H`, i.e. `κ₀ (ṡ, θ̇)ᵀ = J H (s, θ)ᵀ` (Hamilton's equations for the symplectic form
`κ₀ Σ ds_j ∧ dθ_j` with the note's convention `κ₀ ṡ = ∂F/∂θ`, `κ₀ θ̇ = −∂F/∂s`, matching
`κ₀ (Ṡ, Ṫ)ᵀ = J₀ M_m (S, T)ᵀ` of note Prop. `prop:lin`). -/
def linFlow (N : ℕ) (ξ : ℝ) : Matrix (Fin N ⊕ Fin N) (Fin N ⊕ Fin N) ℂ :=
  ((kappa0 ξ : ℂ))⁻¹ • (symplJ N * hessMat N ξ)

/-- The `2 × 2` Hessian block on mode `m`: `M_m = [[a_m, i g_m], [−i g_m, b_m]]` for `m ≠ 0`; on mode `0`
it is `[[c₀, 0], [0, 0]]` with `c₀ = Σ_b ∂²F/∂s_0∂s_b` (`b_0 = g_0 = 0`, `BlocksFourier.blocks_eigen_zero`). -/
def modeBlock (N : ℕ) (ξ : ℝ) (m : Fin N) : Matrix (Fin 2) (Fin 2) ℂ :=
  if (m : ℕ) = 0 then hessBlock (∑ b : Fin N, hessSS N ξ (ringOmega N ξ) m b) 0 0
  else ringHess N m ξ

/-- The linearised flow on mode `m`: `(1/κ₀) J₀ M_m`. -/
def modeFlow (N : ℕ) (ξ : ℝ) (m : Fin N) : Matrix (Fin 2) (Fin 2) ℂ :=
  ((kappa0 ξ : ℂ))⁻¹ • (J0 * modeBlock N ξ m)

theorem sum_mul_left_c {N : ℕ} (f g : Fin N → ℂ) (c : ℂ) :
    ∑ b, f b * (c * g b) = c * ∑ b, f b * g b := by
  rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun b _ => by ring

theorem modeBlock_isHermitian {N : ℕ} (ξ : ℝ) (m : Fin N) : (modeBlock N ξ m).IsHermitian := by
  unfold modeBlock ringHess; split_ifs <;> exact hessBlock_isHermitian _ _ _

section ring
variable {N : ℕ} {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1)
include h0 h1

theorem hessSS_rowsum (Ω : ℝ) (a a' : Fin N) :
    ∑ b, hessSS N ξ Ω a b = ∑ b, hessSS N ξ Ω a' b := by
  have hN : 0 < N := a.pos
  have : NeZero N := ⟨hN.ne'⟩
  have e := eig_hSS ξ Ω a 0
  have e' := eig_hSS ξ Ω a' 0
  simp only [pow_zero, mul_one] at e e'
  apply Complex.ofReal_injective
  push_cast
  simp_rw [hessSS_eq hN h0 h1]
  rw [e, e']

theorem modeBlock_eig (m a : Fin N) :
    ∑ b, ((hessSS N ξ (ringOmega N ξ) a b : ℝ) : ℂ) * fv N m b = modeBlock N ξ m 0 0 * fv N m a
    ∧ ∑ b, ((hessST N ξ (ringOmega N ξ) a b : ℝ) : ℂ) * fv N m b = modeBlock N ξ m 0 1 * fv N m a
    ∧ ∑ b, ((hessTS N ξ (ringOmega N ξ) a b : ℝ) : ℂ) * fv N m b = modeBlock N ξ m 1 0 * fv N m a
    ∧ ∑ b, ((hessTT N ξ (ringOmega N ξ) a b : ℝ) : ℂ) * fv N m b
      = modeBlock N ξ m 1 1 * fv N m a := by
  have hN : 0 < N := a.pos
  by_cases hm : (m : ℕ) = 0
  · obtain ⟨z1, z2, z3⟩ := blocks_eigen_zero h0 h1 (ringOmega N ξ) hN a
    have hr := hessSS_rowsum h0 h1 (ringOmega N ξ) a m
    simp only [modeBlock, if_pos hm, hessBlock, fv_zero_left _ m hm, mul_one]
    simp only [← Complex.ofReal_sum, z1, z2, z3, hr]
    simp
  · obtain ⟨e1, e2, e3, e4⟩ := blocks_eigen h0 h1 m (by omega) m.2 a
    refine ⟨?_, ?_, ?_, ?_⟩
    · simpa [modeBlock, hm, ringHess, hessBlock, fv] using e1
    · simpa [modeBlock, hm, ringHess, hessBlock, fv] using e3
    · simpa [modeBlock, hm, ringHess, hessBlock, fv] using e4
    · simpa [modeBlock, hm, ringHess, hessBlock, fv] using e2

theorem hessMat_modeVec (m : Fin N) (z : Fin 2 → ℂ) :
    hessMat N ξ *ᵥ modeVec m z = modeVec m (modeBlock N ξ m *ᵥ z) := by
  funext i
  rcases i with a | a <;>
  · obtain ⟨e1, e2, e3, e4⟩ := modeBlock_eig h0 h1 m a
    rw [hessMat, fromBlocks_mulVec]
    simp only [modeVec, Sum.elim_inl, Sum.elim_inr, Pi.add_apply,
      mulVec, dotProduct, Function.comp_apply, of_apply, Fin.sum_univ_two, sum_mul_left_c]
    simp only [e1, e2, e3, e4]
    ring

theorem hessMat_isHermitian : (hessMat N ξ).IsHermitian := by
  ext i j
  have hN : 0 < N := by rcases i with a | a <;> exact a.pos
  rcases i with a | a <;> rcases j with b | b <;>
    simp [hessMat, conjTranspose_apply, hessSS_eq hN h0 h1, hessST_eq hN h0 h1, hessTS_eq hN h0 h1,
      hessTT_eq hN h0 h1, hSS_symm N ξ _ a b, hTT_symm N ξ a b]

/-- The Fourier coefficient of `H x` on mode `m` is `M_m` times that of `x` (block diagonalisation). -/
theorem fhat_hessMat (m : Fin N) (x : Fin N ⊕ Fin N → ℂ) :
    fhat m (hessMat N ξ *ᵥ x) = modeBlock N ξ m *ᵥ fhat m x :=
  fhat_mulVec m _ _ (fun z => by
    rw [(hessMat_isHermitian h0 h1).eq, (modeBlock_isHermitian ξ m).eq]
    exact hessMat_modeVec h0 h1 m z) x

end ring

theorem J0_conjTranspose : J0ᴴ = !![0, -1; 1, 0] := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [J0, conjTranspose_apply]

theorem symplJ_modeVec {N : ℕ} (m : Fin N) (z : Fin 2 → ℂ) :
    symplJ N *ᵥ modeVec m z = modeVec m (J0 *ᵥ z) := by
  funext i
  rw [symplJ, fromBlocks_mulVec]
  rcases i with a | a <;>
    simp [modeVec, J0, mulVec, dotProduct, Fin.sum_univ_two, Function.comp_def, one_apply]

theorem symplJ_conj_modeVec {N : ℕ} (m : Fin N) (z : Fin 2 → ℂ) :
    (symplJ N)ᴴ *ᵥ modeVec m z = modeVec m (J0ᴴ *ᵥ z) := by
  rw [J0_conjTranspose]
  funext i
  rw [symplJ, fromBlocks_conjTranspose, fromBlocks_mulVec]
  rcases i with a | a <;>
    simp [modeVec, mulVec, dotProduct, Fin.sum_univ_two, Function.comp_def, one_apply]

theorem fhat_symplJ {N : ℕ} (m : Fin N) (x : Fin N ⊕ Fin N → ℂ) :
    fhat m (symplJ N *ᵥ x) = J0 *ᵥ fhat m x :=
  fhat_mulVec m _ _ (symplJ_conj_modeVec m) x

theorem linFlow_mulVec {N : ℕ} (ξ : ℝ) (x : Fin N ⊕ Fin N → ℂ) :
    linFlow N ξ *ᵥ x = ((kappa0 ξ : ℂ))⁻¹ • (symplJ N *ᵥ (hessMat N ξ *ᵥ x)) := by
  rw [linFlow, smul_mulVec, ← mulVec_mulVec]

section ring2
variable {N : ℕ} {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1)
include h0 h1

/-- Fourier modes are invariant: `L (z ⊗ v_m) = ((1/κ₀) J₀ M_m z) ⊗ v_m` (every `0 ≤ m ≤ N−1`). -/
theorem linFlow_modeVec (m : Fin N) (z : Fin 2 → ℂ) :
    linFlow N ξ *ᵥ modeVec m z = modeVec m (modeFlow N ξ m *ᵥ z) := by
  rw [linFlow_mulVec, hessMat_modeVec h0 h1, symplJ_modeVec, ← modeVec_smul, modeFlow, smul_mulVec,
    ← mulVec_mulVec]

/-- The Fourier coefficients of `L x` on mode `m` are `(1/κ₀) J₀ M_m` times those of `x`. -/
theorem fhat_linFlow (m : Fin N) (x : Fin N ⊕ Fin N → ℂ) :
    fhat m (linFlow N ξ *ᵥ x) = modeFlow N ξ m *ᵥ fhat m x := by
  rw [linFlow_mulVec, fhat_smul, fhat_symplJ, fhat_hessMat h0 h1, modeFlow, smul_mulVec,
    ← mulVec_mulVec]

/-- Spectrum of `L` = union of the spectra of the mode blocks `(1/κ₀) J₀ M_m`, `0 ≤ m ≤ N−1`. -/
theorem linFlow_hasEigenvalue_iff_modeFlow (σ : ℂ) :
    Module.End.HasEigenvalue (Matrix.toLin' (linFlow N ξ)) σ ↔
      ∃ m : Fin N, Module.End.HasEigenvalue (Matrix.toLin' (modeFlow N ξ m)) σ := by
  rw [hasEigenvalue_iff_exists]
  simp_rw [hasEigenvalue_iff_exists]
  constructor
  · rintro ⟨x, hx0, hx⟩
    by_contra hcon
    push Not at hcon
    apply hx0
    apply eq_zero_of_fhat
    intro m
    by_contra hm
    apply hcon m (fhat m x) hm
    rw [← fhat_linFlow h0 h1, hx, fhat_smul]
  · rintro ⟨m, y, hy0, hy⟩
    exact ⟨modeVec m y, modeVec_ne_zero m hy0, by rw [linFlow_modeVec h0 h1, hy, modeVec_smul]⟩

end ring2

/-! ### Mode `0`, the user-facing statements -/

theorem bBlock_zero {N : ℕ} {ξ : ℝ} (hX : ξ ^ N ≠ 1) : bBlock N 0 ξ = 0 := by
  unfold bBlock; rw [Nat.sub_zero, W_zero_eq N ξ hX]; push_cast; ring

theorem gBlock_zero {N : ℕ} {ξ : ℝ} (hX : ξ ^ N ≠ 1) : gBlock N 0 ξ = 0 := by
  unfold gBlock; rw [Nat.sub_zero, W_zero_eq N ξ hX]; ring

theorem kappa0_pos {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) : 0 < kappa0 ξ := by
  unfold kappa0
  have : 0 < 1 - ξ := by linarith
  positivity

/-- `J₀ [[a, 0], [0, 0]]` has the single eigenvalue `0`. -/
theorem flowBlock_zero_hasEigenvalue_iff (a : ℝ) (τ : ℂ) :
    Module.End.HasEigenvalue (Matrix.toLin' (flowBlock a 0 0)) τ ↔ τ = 0 := by
  rw [flowBlock_hasEigenvalue_iff]; simp

theorem modeFlow_of_ne {N : ℕ} (ξ : ℝ) {m : Fin N} (hm : (m : ℕ) ≠ 0) :
    modeFlow N ξ m = ((kappa0 ξ : ℂ))⁻¹ • ringFlow N m ξ := by
  simp [modeFlow, modeBlock, hm, ringFlow]

theorem modeFlow_zero_hasEigenvalue_iff {N : ℕ} {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) {m : Fin N}
    (hm : (m : ℕ) = 0) (σ : ℂ) :
    Module.End.HasEigenvalue (Matrix.toLin' (modeFlow N ξ m)) σ ↔ σ = 0 := by
  have hk : ((kappa0 ξ : ℂ))⁻¹ ≠ 0 :=
    inv_ne_zero (Complex.ofReal_ne_zero.mpr (kappa0_pos h0 h1).ne')
  rw [modeFlow, hasEigenvalue_smul_iff _ hk, modeBlock, if_pos hm, ← flowBlock,
    flowBlock_zero_hasEigenvalue_iff]
  simp [(kappa0_pos h0 h1).ne']

theorem ringFlow_zero_smul_hasEigenvalue_iff {N : ℕ} (hN : 0 < N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1)
    (σ : ℂ) :
    Module.End.HasEigenvalue (Matrix.toLin' (((kappa0 ξ : ℂ))⁻¹ • ringFlow N 0 ξ)) σ ↔ σ = 0 := by
  have hk : ((kappa0 ξ : ℂ))⁻¹ ≠ 0 :=
    inv_ne_zero (Complex.ofReal_ne_zero.mpr (kappa0_pos h0 h1).ne')
  have hX : ξ ^ N ≠ 1 := (pow_lt_one₀ h0.le h1 hN.ne').ne
  rw [hasEigenvalue_smul_iff _ hk, ringFlow, ringHess, ← flowBlock, bBlock_zero hX, gBlock_zero hX,
    flowBlock_zero_hasEigenvalue_iff]
  simp [(kappa0_pos h0 h1).ne']

theorem fhat_modeVec_self {N : ℕ} (m : Fin N) (y : Fin 2 → ℂ) : fhat m (modeVec m y) = (N : ℂ) • y := by
  have e : ∀ c : ℂ, ∑ b, conj (fv N m b) * (c * fv N m b) = N * c := by
    intro c; rw [← fv_norm m, Finset.sum_mul]
    exact Finset.sum_congr rfl fun b _ => by ring
  ext i; fin_cases i <;> simp [fhat, modeVec, e]

theorem natCast_inv_pos {N : ℕ} (hN : 0 < N) : (0 : ℂ) < (N : ℂ)⁻¹ := by
  have : (N : ℂ)⁻¹ = (((N : ℝ)⁻¹ : ℝ) : ℂ) := by push_cast; rfl
  rw [this]; exact Complex.zero_lt_real.mpr (by positivity)

section main
variable {N : ℕ} {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1)
include h0 h1

/-- **note Prop. `prop:lin`, mode invariance.** For `1 ≤ m ≤ N−1` and `(S, T) ∈ ℂ²`, with
`v_b = w^{mb}`: `L (S v, T v) = (1/κ₀) · ((J₀ M_m (S, T)ᵀ) ⊗ v)`, where
`ringFlow N m ξ = J₀ M_m = flowBlock a_m b_m g_m`. -/
theorem linFlow_mode {m : ℕ} (hm1 : 1 ≤ m) (hmN : m < N) (S T : ℂ) :
    linFlow N ξ *ᵥ Sum.elim (fun b : Fin N => S * (ringW N ^ (b : ℕ)) ^ m)
        (fun b : Fin N => T * (ringW N ^ (b : ℕ)) ^ m)
      = ((kappa0 ξ : ℂ))⁻¹ • Sum.elim
          (fun b : Fin N => (ringFlow N m ξ *ᵥ ![S, T]) 0 * (ringW N ^ (b : ℕ)) ^ m)
          (fun b : Fin N => (ringFlow N m ξ *ᵥ ![S, T]) 1 * (ringW N ^ (b : ℕ)) ^ m) := by
  have h := linFlow_modeVec h0 h1 ⟨m, hmN⟩ ![S, T]
  rw [modeFlow_of_ne ξ (m := ⟨m, hmN⟩) (by simp; omega), smul_mulVec, modeVec_smul] at h
  simpa [modeVec, fv] using h

/-- **note Prop. `prop:lin`, spectrum.** `σ` is an eigenvalue of the linearised flow `L` (a complex
`2N × 2N` matrix) iff it is an eigenvalue of `(1/κ₀) J₀ M_m` for some `0 ≤ m ≤ N−1`
(`J₀ M_m = flowBlock a_m b_m g_m`; on mode `0`, `b_0 = g_0 = 0` and the only eigenvalue is `0`). -/
theorem linFlow_hasEigenvalue_iff (σ : ℂ) :
    Module.End.HasEigenvalue (Matrix.toLin' (linFlow N ξ)) σ ↔
      ∃ m : ℕ, m < N ∧
        Module.End.HasEigenvalue (Matrix.toLin' (((kappa0 ξ : ℂ))⁻¹ • ringFlow N m ξ)) σ := by
  rw [linFlow_hasEigenvalue_iff_modeFlow h0 h1]
  constructor
  · rintro ⟨m, hm⟩
    refine ⟨m, m.2, ?_⟩
    by_cases h : (m : ℕ) = 0
    · rw [modeFlow_zero_hasEigenvalue_iff h0 h1 h] at hm
      rw [h, ringFlow_zero_smul_hasEigenvalue_iff m.pos h0 h1]
      exact hm
    · rwa [modeFlow_of_ne ξ h] at hm
  · rintro ⟨m, hmN, hm⟩
    refine ⟨⟨m, hmN⟩, ?_⟩
    by_cases h : m = 0
    · subst h
      rw [ringFlow_zero_smul_hasEigenvalue_iff hmN h0 h1] at hm
      rw [modeFlow_zero_hasEigenvalue_iff h0 h1 rfl]
      exact hm
    · rwa [modeFlow_of_ne ξ (m := ⟨m, hmN⟩) h]

/-- `L` always has the eigenvalue `0` (mode `0`, rotation `(0, 1) ⊗ 1`), for `N ≥ 1`. -/
theorem linFlow_hasEigenvalue_zero (hN : 0 < N) :
    Module.End.HasEigenvalue (Matrix.toLin' (linFlow N ξ)) 0 :=
  (linFlow_hasEigenvalue_iff h0 h1 0).mpr
    ⟨0, hN, (ringFlow_zero_smul_hasEigenvalue_iff hN h0 h1 0).mpr rfl⟩

/-- **note Prop. `prop:lin`, spectral stability.** All eigenvalues of `L` have zero real part iff
`specStable N ξ` (`a_m b_m ≥ 0` for `1 ≤ m ≤ N−1`). -/
theorem linFlow_spectral_iff :
    (∀ σ : ℂ, Module.End.HasEigenvalue (Matrix.toLin' (linFlow N ξ)) σ → σ.re = 0) ↔
      specStable N ξ := by
  have hk0 := kappa0_pos h0 h1
  have hk : ((kappa0 ξ : ℂ))⁻¹ ≠ 0 := inv_ne_zero (Complex.ofReal_ne_zero.mpr hk0.ne')
  rw [specStable_iff_spectral]
  constructor
  · intro H m hm1 hmN σ hσ
    have h2 : Module.End.HasEigenvalue (Matrix.toLin' (((kappa0 ξ : ℂ))⁻¹ • ringFlow N m ξ))
        (((kappa0 ξ : ℂ))⁻¹ * σ) := by
      rwa [hasEigenvalue_smul_iff _ hk, inv_mul_cancel_left₀ hk]
    have h3 := H _ ((linFlow_hasEigenvalue_iff h0 h1 _).mpr ⟨m, hmN, h2⟩)
    rw [← Complex.ofReal_inv, Complex.re_ofReal_mul] at h3
    exact (mul_eq_zero.mp h3).resolve_left (inv_ne_zero hk0.ne')
  · intro H σ hσ
    obtain ⟨m, hmN, hm⟩ := (linFlow_hasEigenvalue_iff h0 h1 σ).mp hσ
    by_cases h : m = 0
    · subst h
      rw [ringFlow_zero_smul_hasEigenvalue_iff hmN h0 h1] at hm
      simp [hm]
    · rw [hasEigenvalue_smul_iff _ hk, inv_inv] at hm
      have := H m (by omega) hmN _ hm
      rw [Complex.re_ofReal_mul] at this
      exact (mul_eq_zero.mp this).resolve_left hk0.ne'

/-- Parseval / block diagonalisation of the Hessian form:
`x* H x = N⁻¹ Σ_m (x̂_m)* M_m x̂_m`. -/
theorem hessMat_parseval (x : Fin N ⊕ Fin N → ℂ) :
    star x ⬝ᵥ (hessMat N ξ *ᵥ x)
      = (N : ℂ)⁻¹ * ∑ m : Fin N, star (fhat m x) ⬝ᵥ (modeBlock N ξ m *ᵥ fhat m x) := by
  have hx := fourier_inv x
  calc star x ⬝ᵥ (hessMat N ξ *ᵥ x)
      = star x ⬝ᵥ (hessMat N ξ *ᵥ ((N : ℂ)⁻¹ • ∑ m : Fin N, modeVec m (fhat m x))) := by
        rw [← hx]
    _ = _ := by
        rw [mulVec_smul, mulVec_sum, dotProduct_smul, dotProduct_sum, smul_eq_mul]
        congr 1
        refine Finset.sum_congr rfl fun m _ => ?_
        rw [hessMat_modeVec h0 h1, dot_modeVec]

/-- **note Prop. `prop:lin`, definiteness.** The (complexified) Hessian `H` of `F` at the ring is positive
definite on the complement of mode `0` — every nonzero `x ∈ ℂ^{2N}` with `Σ_b x_{s,b} = 0` and
`Σ_b x_{θ,b} = 0` (i.e. orthogonal to the mode-`0` subspace) has `x* H x > 0` — iff `emDefinite N ξ`
(`a_m > 0`, `a_m b_m − g_m² > 0` for `1 ≤ m ≤ N−1`). -/
theorem hessMat_posDef_compl_iff :
    (∀ x : Fin N ⊕ Fin N → ℂ, x ≠ 0 → ∑ b, x (Sum.inl b) = 0 → ∑ b, x (Sum.inr b) = 0 →
        0 < star x ⬝ᵥ (hessMat N ξ *ᵥ x)) ↔ emDefinite N ξ := by
  rw [emDefinite_iff_posDef]
  constructor
  · intro H m hm1 hmN
    rw [posDef_iff_dotProduct_mulVec]
    refine ⟨hessBlock_isHermitian _ _ _, fun y hy => ?_⟩
    have hm' : ((⟨m, hmN⟩ : Fin N) : ℕ) ≠ 0 := by simp; omega
    have hsum := fv_sum (⟨m, hmN⟩ : Fin N) hm'
    have hx := H (modeVec ⟨m, hmN⟩ y) (modeVec_ne_zero _ hy)
      (by simp [modeVec, ← Finset.mul_sum, hsum]) (by simp [modeVec, ← Finset.mul_sum, hsum])
    rw [hessMat_modeVec h0 h1, dot_modeVec, fhat_modeVec_self, modeBlock, if_neg hm', star_smul,
      smul_dotProduct, star_natCast, smul_eq_mul] at hx
    have e : star y ⬝ᵥ (ringHess N m ξ *ᵥ y)
        = (N : ℂ)⁻¹ * ((N : ℂ) * (star y ⬝ᵥ (ringHess N m ξ *ᵥ y))) := by
      rw [← mul_assoc, inv_mul_cancel₀ (Nat.cast_ne_zero.mpr (by omega)), one_mul]
    rw [e]
    exact mul_pos (natCast_inv_pos (by omega)) hx
  · intro H x hx0 hs ht
    rw [hessMat_parseval h0 h1]
    have hf0 : ∀ m : Fin N, (m : ℕ) = 0 → fhat m x = 0 := by
      intro m hm; ext i; fin_cases i <;> simp [fhat, fv_zero_left _ m hm, hs, ht]
    obtain ⟨m0, hm0⟩ : ∃ m, fhat m x ≠ 0 := by
      by_contra h; push Not at h; exact hx0 (eq_zero_of_fhat h)
    have hm0' : (m0 : ℕ) ≠ 0 := fun h => hm0 (hf0 m0 h)
    have hpos : ∀ m : Fin N, (m : ℕ) ≠ 0 → fhat m x ≠ 0 →
        0 < star (fhat m x) ⬝ᵥ (modeBlock N ξ m *ᵥ fhat m x) := by
      intro m hm hf
      rw [modeBlock, if_neg hm]
      exact (posDef_iff_dotProduct_mulVec.mp (H m (by omega) m.2)).2 hf
    have hterm : ∀ m ∈ Finset.univ, 0 ≤ star (fhat m x) ⬝ᵥ (modeBlock N ξ m *ᵥ fhat m x) := by
      intro m _
      by_cases hm : (m : ℕ) = 0
      · simp [hf0 m hm]
      · by_cases hf : fhat m x = 0
        · simp [hf]
        · exact (hpos m hm hf).le
    exact mul_pos (natCast_inv_pos m0.pos)
      (Finset.sum_pos' hterm ⟨m0, Finset.mem_univ _, hpos m0 hm0' hm0⟩)

end main

/-! ### The real form: the actual second derivative `hessQ` -/

/-- The real `2N × 2N` Hessian (`hessMat` is its complexification). -/
def hessR (N : ℕ) (ξ : ℝ) : Matrix (Fin N ⊕ Fin N) (Fin N ⊕ Fin N) ℝ :=
  fromBlocks (of fun a b => hessSS N ξ (ringOmega N ξ) a b) (of fun a b => hessST N ξ (ringOmega N ξ) a b)
    (of fun a b => hessTS N ξ (ringOmega N ξ) a b) (of fun a b => hessTT N ξ (ringOmega N ξ) a b)

theorem hessMat_eq_map (N : ℕ) (ξ : ℝ) : hessMat N ξ = (hessR N ξ).map (fun r : ℝ => (r : ℂ)) := by
  ext i j; rcases i with a | a <;> rcases j with b | b <;> simp [hessMat, hessR]

/-- For a real symmetric matrix, `x* H x = Σ_{ij} H_ij (Re x_i Re x_j + Im x_i Im x_j)`. -/
theorem quad_split {ι : Type*} [Fintype ι] (Hr : Matrix ι ι ℝ) (hs : ∀ i j, Hr j i = Hr i j)
    (x : ι → ℂ) :
    star x ⬝ᵥ (Hr.map (fun r : ℝ => (r : ℂ)) *ᵥ x)
      = ((∑ i, ∑ j, Hr i j * ((x i).re * (x j).re + (x i).im * (x j).im) : ℝ) : ℂ) := by
  have hsw : ∑ i, ∑ j, Hr i j * ((x i).re * (x j).im)
      = ∑ i, ∑ j, Hr i j * ((x i).im * (x j).re) := by
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    rw [hs]; ring
  apply Complex.ext
  · simp only [dotProduct, mulVec, Pi.star_apply, map_apply, Complex.re_sum, Finset.mul_sum,
      Complex.ofReal_re]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    simp [Complex.mul_re, Complex.mul_im]; ring
  · simp only [dotProduct, mulVec, Pi.star_apply, map_apply, Complex.im_sum, Finset.mul_sum,
      Complex.ofReal_im]
    have : ∑ i, ∑ j, (star (x i) * ((Hr i j : ℂ) * x j)).im
        = ∑ i, ∑ j, Hr i j * ((x i).re * (x j).im) - ∑ i, ∑ j, Hr i j * ((x i).im * (x j).re) := by
      rw [← Finset.sum_sub_distrib]
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [← Finset.sum_sub_distrib]
      refine Finset.sum_congr rfl fun j _ => ?_
      simp [Complex.mul_re, Complex.mul_im]; ring
    rw [this, hsw, sub_self]

section realForm
variable {N : ℕ} {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1)
include h0 h1

theorem hessR_symm (i j : Fin N ⊕ Fin N) : hessR N ξ j i = hessR N ξ i j := by
  have hN : 0 < N := by rcases i with a | a <;> exact a.pos
  rcases i with a | a <;> rcases j with b | b <;>
    simp [hessR, hessSS_eq hN h0 h1, hessST_eq hN h0 h1, hessTS_eq hN h0 h1,
      hessTT_eq hN h0 h1, hSS_symm N ξ _ a b, hTT_symm N ξ a b]

/-- The second derivative `d²/dt² F(t s, t θ)|₀` is the quadratic form of `hessR`. -/
theorem hessQ_eq_sum (hN : 0 < N) (s θ : Fin N → ℝ) :
    hessQ N ξ (ringOmega N ξ) s θ
      = ∑ i, ∑ j, hessR N ξ i j * (Sum.elim s θ i * Sum.elim s θ j) := by
  rw [hessQ_eq hN h0 h1]
  unfold qForm bil
  have hsw : ∑ a, ∑ b, hST N ξ b a * (θ a * s b) = ∑ a, ∑ b, s a * hST N ξ a b * θ b := by
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring
  simp only [Fintype.sum_sum_type, hessR, fromBlocks_apply₁₁, fromBlocks_apply₁₂, fromBlocks_apply₂₁,
    fromBlocks_apply₂₂, of_apply, Sum.elim_inl, Sum.elim_inr, hessSS_eq hN h0 h1, hessST_eq hN h0 h1,
    hessTS_eq hN h0 h1, hessTT_eq hN h0 h1, Finset.sum_add_distrib, hsw]
  have e1 : ∑ a, ∑ b, hSS N ξ (ringOmega N ξ) a b * (s a * s b)
      = ∑ a, ∑ b, s a * hSS N ξ (ringOmega N ξ) a b * s b :=
    Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring
  have e2 : ∑ a, ∑ b, hST N ξ a b * (s a * θ b) = ∑ a, ∑ b, s a * hST N ξ a b * θ b :=
    Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring
  have e3 : ∑ a, ∑ b, hTT N ξ a b * (θ a * θ b) = ∑ a, ∑ b, θ a * hTT N ξ a b * θ b :=
    Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring
  rw [e1, e2, e3]; ring

/-- **note Prop. `prop:lin`, definiteness (real form).** The second derivative of `F = H − ΩJ` at the ring
(`hessQ`, `d²/dt² F(ts, tθ)|₀`) is positive on every nonzero real `(s, θ)` with `Σ s_b = 0` and
`Σ θ_b = 0` (the complement of mode `0`) iff `emDefinite N ξ`. -/
theorem hessQ_posDef_compl_iff (hN : 0 < N) :
    (∀ s θ : Fin N → ℝ, (s ≠ 0 ∨ θ ≠ 0) → ∑ b, s b = 0 → ∑ b, θ b = 0 →
        0 < hessQ N ξ (ringOmega N ξ) s θ) ↔ emDefinite N ξ := by
  rw [← hessMat_posDef_compl_iff h0 h1]
  have key : ∀ x : Fin N ⊕ Fin N → ℂ, star x ⬝ᵥ (hessMat N ξ *ᵥ x)
      = ((hessQ N ξ (ringOmega N ξ) (fun b => (x (Sum.inl b)).re) (fun b => (x (Sum.inr b)).re)
          + hessQ N ξ (ringOmega N ξ) (fun b => (x (Sum.inl b)).im) (fun b => (x (Sum.inr b)).im)
          : ℝ) : ℂ) := by
    intro x
    rw [hessMat_eq_map, quad_split _ (hessR_symm h0 h1), hessQ_eq_sum h0 h1 hN,
      hessQ_eq_sum h0 h1 hN, ← Finset.sum_add_distrib]
    congr 1
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun j _ => ?_
    rcases i with a | a <;> rcases j with b | b <;> simp <;> ring
  have hQ0 : ∀ s θ : Fin N → ℝ, s = 0 → θ = 0 → hessQ N ξ (ringOmega N ξ) s θ = 0 := by
    rintro s θ rfl rfl; rw [hessQ_eq_sum h0 h1 hN]; simp
  constructor
  · intro H x hx0 hs ht
    rw [key, Complex.zero_lt_real]
    have hsr : ∑ b, (x (Sum.inl b)).re = 0 := by rw [← Complex.re_sum, hs]; rfl
    have hsi : ∑ b, (x (Sum.inl b)).im = 0 := by rw [← Complex.im_sum, hs]; rfl
    have htr : ∑ b, (x (Sum.inr b)).re = 0 := by rw [← Complex.re_sum, ht]; rfl
    have hti : ∑ b, (x (Sum.inr b)).im = 0 := by rw [← Complex.im_sum, ht]; rfl
    have hnn : ∀ s θ : Fin N → ℝ, ∑ b, s b = 0 → ∑ b, θ b = 0 →
        0 ≤ hessQ N ξ (ringOmega N ξ) s θ := by
      intro s θ h1' h2'
      by_cases hz : s ≠ 0 ∨ θ ≠ 0
      · exact (H s θ hz h1' h2').le
      · push Not at hz; rw [hQ0 s θ hz.1 hz.2]
    by_cases hre : (fun b => (x (Sum.inl b)).re) ≠ 0 ∨ (fun b => (x (Sum.inr b)).re) ≠ 0
    · exact add_pos_of_pos_of_nonneg (H _ _ hre hsr htr) (hnn _ _ hsi hti)
    · have him : (fun b => (x (Sum.inl b)).im) ≠ 0 ∨ (fun b => (x (Sum.inr b)).im) ≠ 0 := by
        push Not at hre
        by_contra hc
        push Not at hc
        apply hx0
        funext i
        rcases i with b | b
        · exact Complex.ext (congrFun hre.1 b) (congrFun hc.1 b)
        · exact Complex.ext (congrFun hre.2 b) (congrFun hc.2 b)
      exact add_pos_of_nonneg_of_pos (hnn _ _ hsr htr) (H _ _ him hsi hti)
  · intro H s θ hne hs ht
    have hx := H (Sum.elim (fun b => (s b : ℂ)) (fun b => (θ b : ℂ)))
      (by
        intro h0'
        rcases hne with hne | hne
        · exact hne (funext fun b => by simpa using congrFun h0' (Sum.inl b))
        · exact hne (funext fun b => by simpa using congrFun h0' (Sum.inr b)))
      (by simp [← Complex.ofReal_sum, hs]) (by simp [← Complex.ofReal_sum, ht])
    rw [key, Complex.zero_lt_real] at hx
    simp only [Sum.elim_inl, Sum.elim_inr, Complex.ofReal_re, Complex.ofReal_im] at hx
    rwa [hQ0 (fun _ => 0) (fun _ => 0) rfl rfl, add_zero] at hx

end realForm

end Havelock.Tanh.LargeN
