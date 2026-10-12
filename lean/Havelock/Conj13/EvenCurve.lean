/-
Havelock — Conjecture 13, even `n = 2M`: the Taylor coefficients `h₂`, `h₄` of the reduced
Hamiltonian along the symmetric curve (paper `conj13/paper/main.tex`, proof of Thm even).

On `Fix(D_M) ∩ J⁻¹(μ)` the configuration is two staggered regular `M`-gons with squared radii
`u±`, constrained by `g(u₊) + g(u₋) = 2g(1)`, `g(u) = u/(1+λu)`, parametrised by
`g(u±) = g(1) ± τ`. With `x = (1+λ)τ` this gives `u± = (1±x)/(1∓λx)` (`EvenConfig.lean` proves that
`H(config) − H(ring)` equals `Phi M λ x` below). This file proves

  `Phi M λ ((1+λ)τ) = h₂ τ² + h₄ τ⁴ + O(τ⁶)`,

with `h₂ = Havelock.Conj13.h2`, `h₄ = Havelock.Conj13.h4` exactly the expressions of the paper (the
definitions used in `Conj13Algebra.lean`). Here

  `Phi M λ x = −(M(2M−1)/2)(log(1−x²) + log(1−λ²x²)) − 2M log cosh(M(artanh x + artanh λx)/2)`.
-/
import Havelock.Conj13.Series
import Havelock.Conj13Algebra

open Asymptotics Filter Topology Real

namespace Havelock.Conj13.Even

open Havelock.Conj13.Series

/-- `artanh x`, written with logarithms. -/
noncomputable def ath (x : ℝ) : ℝ := (log (1 + x) - log (1 - x)) / 2

/-- The reduced Hamiltonian along the even symmetric curve, as a function of `x = (1+λ)τ`
(`EvenConfig.H_stag_sub_ring` identifies it with `H(config) − H(ring)`). -/
noncomputable def Phi (M b x : ℝ) : ℝ :=
  -(M * (2 * M - 1) / 2) * (log (1 - x ^ 2) + log (1 - (b * x) ^ 2))
    - 2 * M * log (cosh (M * (ath x + ath (b * x)) / 2))

/-- `artanh x + artanh(λx) = (1+λ)x + (1+λ³)x³/3 + O(x⁵)`. -/
theorem ath_sum_jet (b : ℝ) :
    (fun x => ath x + ath (b * x) - (1 + b) * x - (1 + b ^ 3) / 3 * x ^ 3)
      =O[𝓝 0] (fun x => x ^ 5) := by
  have hy : (fun x : ℝ => b * x) =O[𝓝 0] (fun x => x ^ 1) := by
    have : (fun x : ℝ => x ^ 1 * b) =O[𝓝 0] (fun x => x ^ 1) := mul_cont continuousAt_const
    exact this.congr_left (fun x => by ring)
  have h2 := comp_bigO (by norm_num) artanh_jet hy
  have h2' : (fun x : ℝ => ((log (1 + b * x) - log (1 - b * x)) / 2 - b * x - (b * x) ^ 3 / 3))
      =O[𝓝 0] (fun x => x ^ 5) := by simpa using h2
  refine (artanh_jet.add h2').congr_left (fun x => ?_)
  simp only [ath]
  ring

/-- **Even `n`: the `τ²` and `τ⁴` coefficients.** Along the symmetric curve,
`Phi M λ ((1+λ)τ) − h₂ τ² − h₄ τ⁴ = O(τ⁶)`, where `h₂ = h2 M λ`, `h₄ = h4 M λ` are the
coefficients stated in the paper's proof of Thm even (definitions of `Conj13Algebra`). -/
theorem Phi_jet (M b : ℝ) :
    (fun τ => Phi M b ((1 + b) * τ) - h2 M b * τ ^ 2 - h4 M b * τ ^ 4)
      =O[𝓝 0] (fun τ => τ ^ 6) := by
  -- expansion in `x`
  have hY : (fun x => M * (ath x + ath (b * x)) / 2 - (M * (1 + b) / 2) * x
      - (M * (1 + b ^ 3) / 6) * x ^ 3) =O[𝓝 0] (fun x => x ^ 5) := by
    refine ((ath_sum_jet b).const_mul_left (M / 2)).congr_left (fun x => ?_)
    ring
  have hC := log_cosh_comp_jet hY
  have hL1 : (fun x : ℝ => log (1 - (1 * x) ^ 2) + 1 ^ 2 * x ^ 2 + 1 ^ 4 * x ^ 4 / 2)
      =O[𝓝 0] (fun x => x ^ 6) := log_one_sub_sq_jet 1
  have hL1' : (fun x : ℝ => log (1 - x ^ 2) + x ^ 2 + x ^ 4 / 2) =O[𝓝 0] (fun x => x ^ 6) := by
    simpa using hL1
  have hL2 := log_one_sub_sq_jet b
  have hx : (fun x => Phi M b x
      - (-(M * (2 * M - 1) / 2) * (-(1 + b ^ 2) * x ^ 2 - (1 + b ^ 4) / 2 * x ^ 4)
        - 2 * M * ((M * (1 + b) / 2) ^ 2 / 2 * x ^ 2
          + ((M * (1 + b) / 2) * (M * (1 + b ^ 3) / 6) - (M * (1 + b) / 2) ^ 4 / 12) * x ^ 4)))
      =O[𝓝 0] (fun x => x ^ 6) := by
    have := ((hL1'.add hL2).const_mul_left (-(M * (2 * M - 1) / 2))).sub
      (hC.const_mul_left (2 * M))
    refine this.congr_left (fun x => ?_)
    simp only [Phi]
    ring
  -- substitute `x = (1+λ)τ`
  have hlin : (fun τ : ℝ => (1 + b) * τ) =O[𝓝 0] (fun τ => τ ^ 1) := by
    have : (fun τ : ℝ => τ ^ 1 * (1 + b)) =O[𝓝 0] (fun τ => τ ^ 1) := mul_cont continuousAt_const
    exact this.congr_left (fun τ => by ring)
  have := comp_bigO (by norm_num) hx hlin
  refine (this.trans (by simpa using isBigO_refl (fun x : ℝ => x ^ 6) (𝓝 0))).congr_left (fun τ => ?_)
  simp only [h2, h4]
  ring

end Havelock.Conj13.Even
