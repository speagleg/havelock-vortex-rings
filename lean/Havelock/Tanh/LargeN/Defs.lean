import Mathlib.Analysis.SpecialFunctions.Trigonometric.DerivHyp
import Mathlib.Analysis.SpecialFunctions.Arsinh

/-! # Tanh model, large `N`: the explicit symbol functions

Definitions for `tanh/note_tanh.tex`, section "Large N" (Theorems 11–18). With
`τ = −(N/2) ln ξ` (so `ξ = exp(−2τ/N)`) and `x = 1 − 2m/N`, the Fourier blocks of
`Havelock.Tanh.aBlock`/`bBlock`/`gBlock` are `N²` times the functions below (note Thm 11).
All names follow the note: `Q, R, β, ψ, γ, E_N, Φ_N, Φ, K, 𝓛`. -/

noncomputable section

namespace Havelock.Tanh.LargeN

open Real

/-- `coth t = cosh t / sinh t`. -/
def coth (t : ℝ) : ℝ := cosh t / sinh t

/-- Langevin function `𝓛(η) = coth η − 1/η`. -/
def langevin (η : ℝ) : ℝ := coth η - 1 / η

/-- `ξ = exp(−2τ/N)`, the inverse of `τ = −(N/2) ln ξ`. -/
def xiOf (N : ℕ) (τ : ℝ) : ℝ := exp (-2 * τ / N)

/-- `x = 1 − 2m/N`. -/
def xOf (N m : ℕ) : ℝ := 1 - 2 * (m : ℝ) / N

/-- `Q(x,τ) = cosh τ cosh xτ − x sinh τ sinh xτ`. -/
def Q (x τ : ℝ) : ℝ := cosh τ * cosh (x * τ) - x * sinh τ * sinh (x * τ)

/-- `R(x,τ) = cosh τ sinh xτ − x sinh τ cosh xτ`. -/
def R (x τ : ℝ) : ℝ := cosh τ * sinh (x * τ) - x * sinh τ * cosh (x * τ)

/-- `β(x,τ) = (1−x²)/8 + (1−Q)/(4 sinh²τ)`, the symbol of `b_m` (note Thm 11). -/
def beta (x τ : ℝ) : ℝ := (1 - x ^ 2) / 8 + (1 - Q x τ) / (4 * sinh τ ^ 2)

/-- `ψ(x,τ) = (1−x²)/8 + Q/(4 sinh²τ)`, the `m`-dependent part of `a_m` (note Thm 11). -/
def psi (x τ : ℝ) : ℝ := (1 - x ^ 2) / 8 + Q x τ / (4 * sinh τ ^ 2)

/-- `γ(x,τ) = R/(4 sinh²τ)`, the symbol of `g_m` (note Thm 11 and the paragraph before Lemma 17). -/
def gamma (x τ : ℝ) : ℝ := R x τ / (4 * sinh τ ^ 2)

/-- `k = c/N = coth(τ/N)/N`. -/
def kN (N : ℕ) (τ : ℝ) : ℝ := coth (τ / N) / N

/-- `E_N(τ) = k coth τ − k²/2 − 1/(2N²)`. -/
def EN (N : ℕ) (τ : ℝ) : ℝ := kN N τ * coth τ - kN N τ ^ 2 / 2 - 1 / (2 * (N : ℝ) ^ 2)

/-- `E_∞(τ) = coth τ/τ − 1/(2τ²)`. -/
def Einf (τ : ℝ) : ℝ := coth τ / τ - 1 / (2 * τ ^ 2)

/-- `Φ_N(τ) = E_N − 1/(4 sinh²τ) − ψ(0,τ)`; for even `N` it is `a_{N/2}/N²`. -/
def PhiN (N : ℕ) (τ : ℝ) : ℝ := EN N τ - 1 / (4 * sinh τ ^ 2) - psi 0 τ

/-- `Φ(τ) = E_∞ − 1/(4 sinh²τ) − ψ(0,τ)`, the `N → ∞` limit of `a_{⌊N/2⌋}/N²`. -/
def Phi (τ : ℝ) : ℝ := Einf τ - 1 / (4 * sinh τ ^ 2) - psi 0 τ

/-- `K(x,τ)` of note Lemma 17. -/
def K (x τ : ℝ) : ℝ :=
  cosh τ * (1 - cosh (x * τ)) + x * sinh τ * sinh (x * τ)
    - x ^ 2 * (cosh τ - 1) / 2 * (2 - (1 - x ^ 2) * (cosh τ - 1))

/-- `G = ψ(0) − ψ(x) − γ²/β` of note Lemma 17. -/
def G (x τ : ℝ) : ℝ := psi 0 τ - psi x τ - gamma x τ ^ 2 / beta x τ

end Havelock.Tanh.LargeN
