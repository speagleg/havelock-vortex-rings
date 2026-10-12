import Havelock.Tanh.LargeN.Threshold
import Havelock.Tanh.LargeN.RowStability
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.LinearAlgebra.Eigenspace.Charpoly
import Mathlib.LinearAlgebra.Charpoly.ToMatrix

/-! # Tanh model: linearisation and stability criteria on one mode (note Prop. `prop:lin`)

For real `a, b, g` let `M = [[a, i g], [−i g, b]]` (`hessBlock`, Hermitian) and `L = J₀ M` with
`J₀ = [[0, 1], [−1, 0]]` (`flowBlock`; the linearised flow on one Fourier mode up to the positive factor
`1/κ₀`). Eigenvalues are taken in the sense of `Module.End.HasEigenvalue (Matrix.toLin' L)`.

* `flowBlock_hasEigenvalue_iff`: `σ` is an eigenvalue of `L` iff `(σ + i g)² + a b = 0`;
  `flowBlock_hasEigenvalue_iff_eq`: the eigenvalues are `−i g ± s` with `s² = −a b`.
* `flowBlock_re_eigenvalues_iff`: all eigenvalues of `L` have zero real part iff `a b ≥ 0`.
* `hessBlock_posDef_iff`: `M` is positive definite iff `a > 0` and `a b − g² > 0`.
* Ring: `specStable_iff_spectral`, `emDefinite_iff_posDef` restate `specStable` / `emDefinite`
  (Threshold.lean) as these spectral / definiteness statements for the blocks `a_m, b_m, g_m`, `1 ≤ m ≤ N−1`.
* Row: `rowMat ℓ θ = [[M_xx, M_xy], [M_yx, M_yy]]` (RowSymbol.lean); `row_spectral_iff`,
  `row_posDef_iff` restate note Cor. `cor:row` in these terms.

Not formalised: the derivation of the blocks (or the row's Bloch matrix) from the full Hamiltonian,
i.e. that `L` is the linearisation of the vortex dynamics on mode `m`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real
open Complex (I)
open Matrix ComplexOrder

/-- The symplectic matrix `J₀ = [[0, 1], [−1, 0]]` (note Prop. `prop:lin`). -/
def J0 : Matrix (Fin 2) (Fin 2) ℂ := !![0, 1; -1, 0]
/-- The Hermitian `2 × 2` block `M = [[a, i g], [−i g, b]]` of the Hessian on one Fourier mode
(note Thm 2 / Prop. `prop:lin`), for real `a, b, g`. -/
def hessBlock (a b g : ℝ) : Matrix (Fin 2) (Fin 2) ℂ := !![(a : ℂ), I * g; -I * g, b]
/-- The linearised flow on one mode, `L = J₀ M` (note Prop. `prop:lin`: `κ₀ (Ṡ, Ṫ)ᵀ = J₀ M (S, T)ᵀ`;
the positive factor `1/κ₀` is dropped). -/
def flowBlock (a b g : ℝ) : Matrix (Fin 2) (Fin 2) ℂ := J0 * hessBlock a b g

/-- note Prop. `prop:lin`: `σ` is an eigenvalue of `L = J₀ M` iff `(σ + i g)² + a b = 0`
(the characteristic polynomial of `L` vanishes at `σ`). -/
theorem flowBlock_hasEigenvalue_iff (a b g : ℝ) (σ : ℂ) :
    Module.End.HasEigenvalue (Matrix.toLin' (flowBlock a b g)) σ ↔
      (σ + I * g) ^ 2 + a * b = 0 := by
  rw [Module.End.hasEigenvalue_iff_isRoot_charpoly, Matrix.charpoly_toLin',
    Polynomial.IsRoot.def, Matrix.eval_charpoly, Matrix.det_fin_two]
  have hI := Complex.I_sq
  simp [flowBlock, J0, hessBlock, Matrix.scalar_apply]
  constructor <;> intro h <;> linear_combination h + 0 * hI

/-- note Prop. `prop:lin`: the eigenvalues of `L` are exactly `−i g ± s` for any square root `s` of `−a b`. -/
theorem flowBlock_hasEigenvalue_iff_eq (a b g : ℝ) {s : ℂ} (hs : s ^ 2 = -(a * b)) (σ : ℂ) :
    Module.End.HasEigenvalue (Matrix.toLin' (flowBlock a b g)) σ ↔
      σ = -(I * g) + s ∨ σ = -(I * g) - s := by
  have e : (σ + I * g) ^ 2 + a * b = 0 ↔ (σ + I * g) ^ 2 = s ^ 2 := by
    rw [hs]; constructor <;> intro h <;> linear_combination h
  rw [flowBlock_hasEigenvalue_iff, e, sq_eq_sq_iff_eq_or_eq_neg]
  constructor <;> rintro (h | h)
  · left; linear_combination h
  · right; linear_combination h
  · left; linear_combination h
  · right; linear_combination h

/-- note Prop. `prop:lin` (spectral criterion): every eigenvalue of `L` has zero real part iff `a b ≥ 0`.
(If `a b < 0`, `−i g + √(−a b)` is an eigenvalue with positive real part.) -/
theorem flowBlock_re_eigenvalues_iff (a b g : ℝ) :
    (∀ σ, Module.End.HasEigenvalue (Matrix.toLin' (flowBlock a b g)) σ → σ.re = 0) ↔
      0 ≤ a * b := by
  constructor
  · intro h
    by_contra hneg
    replace hneg := not_le.mp hneg
    have hs : ((Real.sqrt (-(a * b)) : ℝ) : ℂ) ^ 2 = -(a * b) := by
      rw [← Complex.ofReal_pow, Real.sq_sqrt (by linarith)]; push_cast; ring
    have he := h _ ((flowBlock_hasEigenvalue_iff_eq a b g hs _).2 (Or.inl rfl))
    simp at he
    have : 0 < Real.sqrt (-(a * b)) := Real.sqrt_pos.2 (by linarith)
    linarith
  · intro hab σ hσ
    rw [flowBlock_hasEigenvalue_iff] at hσ
    have hre := congrArg Complex.re hσ
    have him := congrArg Complex.im hσ
    simp [sq] at hre him
    have h2 : σ.re * (σ.im + g) = 0 := by linarith
    rcases mul_eq_zero.mp h2 with h | h
    · exact h
    · rw [h] at hre
      nlinarith [sq_nonneg σ.re]

/-- `M` is Hermitian. -/
theorem hessBlock_isHermitian (a b g : ℝ) : (hessBlock a b g).IsHermitian := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [hessBlock, Matrix.conjTranspose_apply]

/-- note Prop. `prop:lin` (definiteness): `M` is positive definite iff `a > 0` and `a b − g² > 0`. -/
theorem hessBlock_posDef_iff (a b g : ℝ) :
    (hessBlock a b g).PosDef ↔ 0 < a ∧ 0 < a * b - g ^ 2 := by
  rw [Matrix.posDef_iff_dotProduct_mulVec]
  have key : ∀ x : Fin 2 → ℂ, star x ⬝ᵥ (hessBlock a b g *ᵥ x) =
      ((a * ((x 0).re ^ 2 + (x 0).im ^ 2) + b * ((x 1).re ^ 2 + (x 1).im ^ 2)
        - 2 * g * ((x 0).re * (x 1).im - (x 0).im * (x 1).re) : ℝ) : ℂ) := by
    intro x
    apply Complex.ext <;>
      simp [hessBlock, dotProduct, Matrix.mulVec, Fin.sum_univ_two, sq] <;> ring
  simp only [key, Complex.zero_lt_real]
  constructor
  · rintro ⟨-, h⟩
    have h1 := h (x := ![1, 0]) (by intro h0; simpa using congrFun h0 0)
    simp at h1
    refine ⟨h1, ?_⟩
    have h2 := h (x := ![(g : ℂ), I * a]) (by
      intro h0; have := congrFun h0 1; simp at this; linarith)
    simp at h2
    nlinarith
  · rintro ⟨ha, hd⟩
    refine ⟨hessBlock_isHermitian a b g, fun x hx => ?_⟩
    have hnu := Complex.normSq_nonneg (x 0)
    have hnv := Complex.normSq_nonneg (x 1)
    rw [Complex.normSq_apply] at hnu hnv
    by_cases hv : x 1 = 0
    · have hu : x 0 ≠ 0 := by
        intro hu; apply hx; ext i; fin_cases i <;> simp [hu, hv]
      have hpos := Complex.normSq_pos.mpr hu
      rw [Complex.normSq_apply] at hpos
      simp only [hv, Complex.zero_re, Complex.zero_im]
      nlinarith
    · have hpos := Complex.normSq_pos.mpr hv
      rw [Complex.normSq_apply] at hpos
      have hQ : 0 < a * (a * ((x 0).re ^ 2 + (x 0).im ^ 2) + b * ((x 1).re ^ 2 + (x 1).im ^ 2)
          - 2 * g * ((x 0).re * (x 1).im - (x 0).im * (x 1).re)) := by
        nlinarith [sq_nonneg (a * (x 0).re - g * (x 1).im), sq_nonneg (a * (x 0).im + g * (x 1).re),
          mul_pos hd hpos]
      exact pos_of_mul_pos_right hQ ha.le

/-! ### The ring: blocks `a_m, b_m, g_m` (note Thm 2) -/

/-- The Hessian block of the ring on mode `m` (note Thm 2). -/
def ringHess (N m : ℕ) (ξ : ℝ) : Matrix (Fin 2) (Fin 2) ℂ :=
  hessBlock (aBlock N m ξ) (bBlock N m ξ) (gBlock N m ξ)

/-- The linearised flow of the ring on mode `m`, `J₀ M_m` (note Prop. `prop:lin`). -/
def ringFlow (N m : ℕ) (ξ : ℝ) : Matrix (Fin 2) (Fin 2) ℂ := J0 * ringHess N m ξ

/-- note Prop. `prop:lin`: `specStable N ξ` (all `a_m b_m ≥ 0`) iff, for every `1 ≤ m ≤ N−1`, every
eigenvalue of the mode-`m` flow `J₀ M_m` has zero real part. -/
theorem specStable_iff_spectral (N : ℕ) (ξ : ℝ) :
    specStable N ξ ↔ ∀ m : ℕ, 1 ≤ m → m < N →
      ∀ σ : ℂ, Module.End.HasEigenvalue (Matrix.toLin' (ringFlow N m ξ)) σ → σ.re = 0 := by
  unfold specStable
  refine forall_congr' fun m => forall_congr' fun _ => forall_congr' fun _ => ?_
  exact (flowBlock_re_eigenvalues_iff _ _ _).symm

/-- note Prop. `prop:lin`: `emDefinite N ξ` iff every block `M_m`, `1 ≤ m ≤ N−1`, is positive definite. -/
theorem emDefinite_iff_posDef (N : ℕ) (ξ : ℝ) :
    emDefinite N ξ ↔ ∀ m : ℕ, 1 ≤ m → m < N → (ringHess N m ξ).PosDef := by
  unfold emDefinite
  refine forall_congr' fun m => forall_congr' fun _ => forall_congr' fun _ => ?_
  exact (hessBlock_posDef_iff _ _ _).symm

/-! ### The horocyclic row (note Prop. `prop:row`, Cor. `cor:row`) -/

/-- The Bloch matrix of the row, `[[M_xx, M_xy], [M_yx, M_yy]]` in the basis `(X, Y)` (RowSymbol.lean). -/
def rowMat (ℓ θ : ℝ) : Matrix (Fin 2) (Fin 2) ℂ := !![rowMxx ℓ θ, rowMxy ℓ θ; rowMyx ℓ θ, rowMyy ℓ θ]

theorem rowMyx_eq_ofReal {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    rowMyx ℓ θ = -I * (((2 * π / ℓ) ^ 2 * gamma (1 - 2 * θ) (2 * π / ℓ) : ℝ) : ℂ) :=
  (hasSum_termYX hℓ hθ).tsum_eq

/-- The row's Bloch matrix is a block of the form `hessBlock a b g` with `a = τ²β`,
`b = τ²(E_∞ − 1/(4 sinh²τ) − ψ)`, `g = τ²γ` (`τ = 2π/ℓ`, `x = 1 − 2θ`). -/
theorem rowMat_eq {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    rowMat ℓ θ = hessBlock ((2 * π / ℓ) ^ 2 * beta (1 - 2 * θ) (2 * π / ℓ))
      ((2 * π / ℓ) ^ 2 * (Einf (2 * π / ℓ) - 1 / (4 * sinh (2 * π / ℓ) ^ 2)
        - psi (1 - 2 * θ) (2 * π / ℓ)))
      ((2 * π / ℓ) ^ 2 * gamma (1 - 2 * θ) (2 * π / ℓ)) := by
  rw [rowMat, hessBlock, rowMxx_eq_ofReal hℓ hθ, rowMxy_eq_ofReal hℓ hθ, rowMyx_eq_ofReal hℓ hθ,
    rowMyy_eq_ofReal hℓ hθ]

/-- note Cor. `cor:row` (spectral form): every eigenvalue of `J₀ M(θ)` has zero real part for every
`0 < θ < 1` iff `ℓ ≥ ℓ_* = 2π/τ_∞`. -/
theorem row_spectral_iff {ℓ : ℝ} (hℓ : 0 < ℓ) :
    (∀ θ ∈ Set.Ioo (0 : ℝ) 1, ∀ σ : ℂ,
        Module.End.HasEigenvalue (Matrix.toLin' (J0 * rowMat ℓ θ)) σ → σ.re = 0) ↔
      2 * π / tauInf ≤ ℓ := by
  rw [← row_stable_iff hℓ]
  refine forall_congr' fun θ => forall_congr' fun hθ => ?_
  have hθ' : θ ∈ Set.Icc (0 : ℝ) 1 := Set.Ioo_subset_Icc_self hθ
  rw [rowMat_eq hℓ hθ', ← flowBlock, flowBlock_re_eigenvalues_iff, rowMxx_eq_ofReal hℓ hθ',
    rowMyy_eq_ofReal hℓ hθ', Complex.ofReal_re, Complex.ofReal_re]

/-- note Cor. `cor:row` (definiteness form): `M(θ)` is positive definite for every `0 < θ < 1` iff
`ℓ > ℓ_* = 2π/τ_∞`. -/
theorem row_posDef_iff {ℓ : ℝ} (hℓ : 0 < ℓ) :
    (∀ θ ∈ Set.Ioo (0 : ℝ) 1, (rowMat ℓ θ).PosDef) ↔ 2 * π / tauInf < ℓ := by
  rw [← row_definite_iff hℓ]
  refine forall_congr' fun θ => forall_congr' fun hθ => ?_
  have hθ' : θ ∈ Set.Icc (0 : ℝ) 1 := Set.Ioo_subset_Icc_self hθ
  have hβ : 0 < (2 * π / ℓ) ^ 2 * beta (1 - 2 * θ) (2 * π / ℓ) :=
    mul_pos (by positivity) (beta_pos (abs_one_sub_two_mul_lt hθ) (by positivity))
  rw [rowMat_eq hℓ hθ', hessBlock_posDef_iff, rowMxx_eq_ofReal hℓ hθ', rowMyy_eq_ofReal hℓ hθ',
    rowMxy_eq_ofReal hℓ hθ', Complex.ofReal_re, Complex.ofReal_re, Complex.normSq_mul,
    Complex.normSq_I, Complex.normSq_ofReal, one_mul, and_iff_right hβ, ← sq]

end Havelock.Tanh.LargeN
