/-
Havelock — tanh model, large N (tanh/note_tanh.tex, section "Large N", Theorems 11–18, Prop 19, Cor 20).

Modules, in dependency order (CI builds them one at a time in this order):
* `LargeN.Defs`             the explicit symbol functions `Q, R, β, ψ, γ, E_N, Φ_N, Φ, K, G, 𝓛`.
* `LargeN.Symbol`           note Thm 11: the blocks are `N²` times the symbol (`aBlock_xiOf`, ...).
* `LargeN.Stiffness`        note Lemma 12: `beta_pos`, `bBlock_pos` (all `N`, `m`, radii).
* `LargeN.EMCore`           note Lemma 17 (analytic part): `K_pos`, `beta_mul_G`.
* `LargeN.PhiLimit`         eq. (phi), `Φ` strictly decreasing, its root `tauInf`, `kappaInf`, bracket.
* `LargeN.CritMode`         note Lemma 13: `psi_strictAntiOn`; `aBlock_half_le/lt` (critical mode, all `N`).
* `LargeN.EvenN`            note Prop 15: even-`N` closed form, monotonicity, limits, one threshold.
* `LargeN.Threshold`        note Cor 14 and Lemma 17's block identity: `specStable`, `emDefinite`,
                            `xiStar`, `xiEM`, `specStable_iff`, `det_xiOf`, `emDefinite_iff_even`.
* `LargeN.ThresholdLimit`   note Thms 16, 18(c): `cR_xiStar_rate`, `tendsto_cR_xiStar`, `cR_xiEM_rate`.
* `LargeN.ThresholdKrein`   note Thm 18(b): `krein_window_odd`, `xiEM_sub_xiStar_le` (O(N⁻³)).
* `LargeN.ThresholdEven`    note Thm 18(a), Prop 15: `xiEM_eq_xiStar_of_even`, `specStable_all_of_even_le_six`.
* `LargeN.ThresholdSpacing` note Thm 16: `tendsto_spacing`, `spacing_rate` (ℓ_N = ℓ* + O(N⁻²)).
* `LargeN.RowFourier`       the four lattice sums of note Prop 19 (via Fourier coefficients on the circle).
* `LargeN.RowSymbol`        note Prop 19: `row_bloch_symbol` (Bloch symbol of the horocyclic row, velocity `v`).
* `LargeN.RowStability`     note Cor 20: `row_stable_iff`, `row_definite_iff` (critical spacing ℓ* = 2π/τ∞).

Not formalised: the derivation of the Fourier blocks (Thm 2 is the starting point, as in
`TanhAlgebra.lean`); for the row, the per-pair Hessian entries are taken as definitions from the proof of
Prop 19 (the second derivatives of the pair term are not derived in Lean); and the dynamics (stability is
stated through the block/Bloch-matrix criteria of note Prop 3).
-/
import Havelock.Tanh.LargeN.Defs
import Havelock.Tanh.LargeN.Symbol
import Havelock.Tanh.LargeN.Stiffness
import Havelock.Tanh.LargeN.EMCore
import Havelock.Tanh.LargeN.PhiLimit
import Havelock.Tanh.LargeN.CritMode
import Havelock.Tanh.LargeN.EvenN
import Havelock.Tanh.LargeN.Threshold
import Havelock.Tanh.LargeN.ThresholdLimit
import Havelock.Tanh.LargeN.ThresholdKrein
import Havelock.Tanh.LargeN.ThresholdEven
import Havelock.Tanh.LargeN.ThresholdSpacing
import Havelock.Tanh.LargeN.RowFourier
import Havelock.Tanh.LargeN.RowSymbol
import Havelock.Tanh.LargeN.RowStability
