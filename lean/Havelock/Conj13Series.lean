/-
Havelock — Conjecture 13 paper: the series inputs of `Conj13Algebra.lean` derived from the Hamiltonian
(paper eq. (1), H = −Σ_{i<j} log(|z_i − z_j|²/(σ_iσ_j)), σ = 1 + λ|z|²).

* `Conj13.Series`: one-variable jets (log, artanh, cosh, log cosh) with `IsBigO` remainders.
* `Conj13.EvenCurve`, `Conj13.EvenConfig`, `Conj13.EvenQuartic`: even n = 2M — the staggered configurations
  on Fix(D_M) ∩ J⁻¹(μ) (`Ham_stag_sub`, `J3_curve`, `Jp_stag`), `Ham_curve_jet` (H − H₀ = h₂τ² + h₄τ⁴ + O(τ⁶)),
  the ε-normalisation (`epsK_stag`), and `c4even_limit_at_bn` ((H − H₀)/ε⁴ → c4even at b = b_n).
* `Conj13.OddSums`, `Conj13.OddQuadratic`, `Conj13.OddConstraint`: odd n — the pair sums, `Ham_oddLine_jet`
  (H − H₀ = f₂A² + O(A³) along the kernel line), `J1_oddAnsatz_jet` (the J₁ constraint `J1con`), and the
  ε-normalisation 64/(n²(1+λ)⁸) (`epsOdd_hasDerivAt`, `oddNorm_pow4`).
* `Conj13.TaylorCoeffs`: the same as iterated derivatives (`even_iteratedDeriv`, `odd_iteratedDeriv`, `_at_bn`).

Not derived here: the odd fourth-order coefficients E00, E10, E01, E20, E02 (still the stated definitions in
`Conj13Algebra`), and the structural steps (symmetric criticality, completeness of the odd ansatz, identifying
c₄ with the reduced Hamiltonian's quartic coefficient).
-/
import Havelock.Conj13.Series
import Havelock.Conj13.EvenCurve
import Havelock.Conj13.EvenConfig
import Havelock.Conj13.EvenQuartic
import Havelock.Conj13.OddSums
import Havelock.Conj13.OddQuadratic
import Havelock.Conj13.OddConstraint
import Havelock.Conj13.TaylorCoeffs
