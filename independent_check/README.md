# Independent second-reader check

This folder is a check of both results written from scratch by a separate AI agent (Claude) that had not
seen the derivations. It imports nothing from `havelock`, `conj13/`, `tanh/` or `sinh/`, uses none of
their series coefficients or Fourier blocks, and works directly from the Hamiltonian in Cartesian
coordinates on the Poincaré disk. Its numbers are NUMERICAL (40–50-digit floating point, no rigorous error
bound), so they corroborate the papers but do not replace their proofs or certificates.

Run everything with `./run_all.sh` (needs `mpmath` and `sympy`; about 30 minutes, dominated by the tanh
N = 16 bisection). Outputs are in `out/`.

## Setting

Poincaré disk, curvature −1, area form 4/(1−|z|²)² dx∧dy, N identical unit vortices.

| model | pair energy |
|---|---|
| tanh (Kimura) | −ln tanh²(d/2) = −ln(\|z−w\|²/\|1−z w̄\|²) |
| sinh (Montaldi–Tokieda (1.12), λ = −1) | −ln sinh²(d/2) = −ln(\|z−w\|²/((1−\|z\|²)(1−\|w\|²))) |
| plane (`plane.py`) | −ln\|z−w\|² |

Momenta: J₃ = Σ 2/(1−|z|²) (rotations), J₁ = Σ 4x/(1−|z|²), J₂ = Σ 4y/(1−|z|²). For the ring
z_k = ρ e^{2πik/N}, ξ = tanh²(r/2) = ρ². Gradients and Hessians are analytic (Wirtinger calculus with z and
z̄ treated as independent), not finite differences.

## What each script checks

| script | check | against |
|---|---|---|
| `vort.py` | Library: Hamiltonians, momenta, relative equilibria, the linearisation Π·Hess(H − ΩJ₃). | — |
| `bis.py` | Spectral threshold by bisection on max\|Re λ\| of the full 2N×2N linearisation. The tolerance 10^(−dps/4) sits above the √ε noise from the Jordan block at 0 that every relative equilibrium has. | `docs/results.md` tanh table; sinh q(N) |
| `levelset.py` | At the sinh bifurcation point cosh²r* = q(N): the reduced Hessian on ker dJ₃ (its kernel and the sign of the rest), then the quartic E(s) = min over transverse w of F(x₀ + s k + w) − F(x₀) on the level set J₃ = J₃(x₀), solved by Newton on the Lagrange system, for several kernel directions k and both signs of s. E₄ is extrapolated from s = 2·10⁻³ and 10⁻³. | Theorems even/odd: c₄ closed forms |
| `levelset_full.py` | The same with all three momenta J₁, J₂, J₃ fixed. | `levelset.py` (odd n, where the slaved mode 2m ≡ −1 meets m = 1) |
| `plane.py` | Plane, n = 7 (b₇ = 0), ring radius 1, all momenta fixed. | c₄(7) = 9/14 |
| `krein.py` | tanh N = 9 just above ξ*: spectral stability and the inertia of the reduced Hessian (all three momenta). | Krein-window width 3.220·10⁻⁸ |
| `arith.py` | The stored f_N (`tanh/irreducibility/polys.gp`): palindromy, irreducibility over ℚ (sympy, not PARI), unique root in (0, 1) equal to the table, unit leading coefficient, and whether 1/ξ*(N) is a Salem number. | Table, degrees, unit claim |
| `table.py` | Assembles ξ*²·E₄ for every n and compares it with the closed forms. | Closed forms |

## Results (as recorded in `out/`)

**c₄ closed forms.** With E₄ the level-set quartic for a unit Cartesian kernel vector, ξ*²·E₄ reproduces
the closed forms for every n tested, odd and even, under that one normalisation:

| n | ξ*²·E₄ | closed form | relative difference |
|---|---|---|---|
| 7 (plane, radius 1, E₄ itself) | 0.642857138 | 9/14 | 7·10⁻⁹ |
| 8–21 (sinh) | see `out/table.out` | (n−2)(m⁴+40m²+48m+16)/(192n(n−1)) or (n−1)(n−3)(n²−9n+38)/(128n) | ≤ 7·10⁻⁹ |

The residuals come from the small-s extrapolation. At every bifurcation point tested the reduced Hessian
on ker dJ₃ has exactly a two-dimensional kernel besides the rotation zero, and is positive on the rest
(the paper's Proposition on spectral positivity). E(s)/s⁴ is the same for every kernel direction and both
signs of s (no cubic term). For n = 9 and 13, fixing J₁ and J₂ as well does not change E₄.

**tanh thresholds.** The bisection reproduces the tables in `docs/results.md` to every digit computed
(≥ 20) for N = 8, 9, 10, 11, 12, 16; the sinh runs reproduce cosh²r* = q(N) for N = 8, 9, 13.

**Krein window, N = 9.** For ξ* < ξ < ξ* + 3.0·10⁻⁸ the ring is spectrally stable (max\|Re λ\| ~ 10⁻²⁴)
while the reduced Hessian has two negative directions; at ξ* + 3.4·10⁻⁸ it is definite. This agrees with
the certified width 3.220·10⁻⁸.

**Arithmetic.** All nine f_N are irreducible over ℚ and palindromic with the stated degrees, each has exactly
one root in (0, 1), equal to the table; only f₈ is monic. 1/ξ*(N) is **not** a Salem number for any
N = 8–16 (several conjugates lie off the unit circle).

## What this does not check

- The nonlinear-stability argument of `conj13/paper` §5 (the cross-section, the penalty function, and
  stability for perturbations of any momentum). This check confirms its input, a strict quartic minimum,
  not the argument.
- The sphere cases n = 4, 5, 6 and the case n = 3 (the code is disk-only, plus the plane at n = 7).
- The Lean build, PARI, or any Arb certificate.
- Prior work (Hwang–Kim 2013, Kurakin 2003/2004, Kurakin–Yudovich) and the Boatto figure comparison.
