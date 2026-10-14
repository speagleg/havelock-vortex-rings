# Exploration, 2026-10-14: what stabilises a vortex row on H^2?

Continues `../2026-10-13/`. NUMERICAL and exploratory except where a result is marked as derived; no specialist has
reviewed it. Run each script from its own folder.

## rows/ — rows on curves of constant geodesic curvature (tanh model)

Question: the horocycle row is stable once its spacing exceeds l* = 0.8439 (note_tanh, Cor. 20), while a straight row
in the plane is never stable. Is it the ambient curvature of H^2 that stabilises, or the bending of the curve that
carries the row? H^2 has a one-parameter family of homogeneous curves: geodesics (geodesic curvature kappa = 0),
hypercycles (0 < kappa < 1), horocycles (kappa = 1) and circles (kappa > 1, i.e. the finite rings of the note).
A row equally spaced along any of them is a relative equilibrium.

| Script | What it does |
|---|---|
| `rows.py` | Linear stability of hypercycle rows (log coordinates z = exp(u + i phi), exponentially convergent images) and of the horocycle row (power-law images), supercell of M vortices |
| `asymptotic_check.py` | sympy check of the large-spacing criterion below |
| `rows_boundary_gen.py` | Critical spacing l*(kappa) on the hypercycle branch -> `rows_boundary.txt` |
| `circle_branch.py` | Circle branch from the package's CERTIFIED thresholds (kappa = coth r*, l = 2 pi sinh r*/N) -> `circle_branch.out` |
| `plot_rows.py` | `rows_boundary.png` |
| `rows_kernel.py` | Same analysis for a screened kernel, W = Q_1(cosh d) (Green's function of Laplacian - 2 on H^2, decays like exp(-2d)) |
| `pairing_root.py` | The threshold as the zero of one lattice sum (pairing-mode stiffness); critical Bloch mode |
| `next_order.py` | sympy: exact expansion of that lattice sum to third order in exp(-l sin alpha) |
| `boundary_series.py` | One-, two- and three-term asymptotic boundaries vs the full threshold |
| `collar.py`, `plot_collar.py` | Necklaces of N vortices around the Bolza systole in its collar (hyperbolic-cylinder model) -> `collar.out`, `collar.png` |
| `pair_terms.py`, `pair_sign.py` | sympy: closed form of each pair's contribution to the pairing-mode symbol, and its sign |
| `instability_proof_check.py` | Checks every step of the proof in item 9 -> `instability_proof_check.out` |
| `symbol_q.py` | sympy: pair terms of the full Bloch symbol at every wavenumber q |
| `exact_symbol.py` | The exact symbol and the criterion F(q) >= 0, checked against the supercell -> `exact_symbol.out` |
| `stability_bound.py` | Proofs of item 10 (explicit stability bound, all-mode instability for kappa <= 1/2) -> `stability_bound.out` |

**Results.**

1. *Checks.* Horocycle threshold falls in (0.840, 0.845) vs the note's 0.8439085; log-coordinate velocities agree
   with the upper-half-plane formula to 2e-12.
2. *Geodesic and weakly curved rows are never stable.* kappa = 0 and 0.5: unstable at every spacing tested
   (0.5 to 80 curvature radii). So ambient curvature alone does not stabilise a row; the curve must bend.
3. *Critical curvature kappa_c = 1/sqrt 2.* Bisection at spacing 40 and 80 gives kappa_c = 0.7071064 (1/sqrt 2 =
   0.7071068). DERIVED (leading order, large spacing): only nearest neighbours matter, W ~ 2 sin(phi_i) sin(phi_j)
   exp(-|u_i - u_j|); the linearised Bloch symbol has equal diagonal entries and off-diagonal entries proportional
   to cos^2(alpha) cos q + 2 cos^2(alpha) - sin^2(alpha) and to cos q - 1 <= 0. Stability for all q therefore needs
   the first to be >= 0; the worst mode is the pairing mode q = pi, where it equals cos(2 alpha). Hence stable at
   large spacing iff kappa = cos(alpha) >= 1/sqrt 2 (alpha = 45 deg in the upper half-plane; the hypercycle at
   distance arsinh 1 = 0.8814 from its geodesic). Not proved: error bounds on the asymptotics, and instability for
   kappa < 1/sqrt 2 at every finite spacing (numerical only).
4. *One boundary for all four curve types.* l*(kappa): 4.76 (kappa = 0.72), 3.19 (0.75), 2.26 (0.80), 1.44 (0.90),
   1.01 (0.98), 0.874 (0.999), approaching the horocycle value 0.8439 roughly like sqrt(1 - kappa); the rings
   continue it for kappa > 1: (1.007, 0.839) at N = 64 down to (2.056, 0.437) at N = 8. The finite rings, the
   horocycle row and the hypercycle rows lie on one continuous stability boundary, which diverges at kappa = 1/sqrt 2.

5. *kappa_c does not depend on the kernel.* Redoing the asymptotics for any kernel decaying like exp(-gamma d)
   gives the normal-stiffness factor gamma cos^2(alpha) cos q + gamma cos^2(alpha) + 2 cos^2(alpha) - 1, whose minimum
   over q (at q = pi) is cos(2 alpha) for every gamma: kappa_c = 1/sqrt 2 is universal for exponentially screened
   interactions on H^2. Numerical check with W = Q_1(cosh d): kappa_c = 0.7071061 at spacings 30 and 60. The critical
   spacing itself is kernel-dependent: Q_1 gives l* = 4.09 (kappa = 0.75), 3.01 (0.8), 1.98 (0.9), 1.69 (0.95),
   1.49 (0.99), versus 3.19, 2.26, 1.44, 1.17, 0.95 for the tanh kernel.

6. *The threshold is one lattice sum.* Just below the threshold the growing mode is always the pairing mode
   (Bloch wavenumber 1/2), and at q = pi the Bloch symbol has zero diagonal. So the row is stable iff
   A_UP(pi)(alpha, l') >= 0, an explicit sum over the row. Its zero reproduces the full supercell threshold to 4-5
   digits for kappa <= 0.98 (e.g. 4.7635 vs 4.7635 at kappa = 0.72, 1.4366 vs 1.4366 at 0.90).
7. *Asymptotics near kappa_c.* With e = exp(-l sin alpha), DERIVED (sympy, checked against the exact sum):
   A_UP(pi)/(4 sin^2 alpha) = (2k^2 - 1) e + (8k^4 - 3) e^2 + 4(8k^6 - 8k^4 + 2k^2 - 1) e^3 + O(e^4).
   At kappa_c the e^2 coefficient is exactly -1, so as kappa -> 1/sqrt 2 from above
   l*(kappa) ~ -ln(2 kappa^2 - 1)/sqrt(1 - kappa^2): within 0.1% of the full threshold at kappa = 0.72, 1.3% at 0.75,
   6% at 0.80. The series is only asymptotic (truncations alternate and fail beyond kappa ~ 0.8); the lattice sum
   in item 6 is the usable criterion.

8. *Bolza systole: no stable necklace inside the collar.* Near a closed geodesic of length L the surface is the
   hyperbolic cylinder H^2/<translation by L>, and a necklace of N vortices at distance h is a hypercycle row with
   period L. For the Bolza systole (L = 2 arccosh(1 + sqrt 2) = 3.0571) the collar half-width is w = 0.4407 and the
   largest curvature inside it is exactly tanh w = sqrt 2 - 1 = 0.414 < 1/sqrt 2. Result: for N = 2..12 (odd and even)
   no necklace is stable anywhere in the collar; each first stabilises at h* = 1.06 (N = 2) to 2.11 (N = 12), beyond
   the collar. Even-N thresholds lie on the infinite-row boundary; odd N (no exact pairing mode) sit slightly below
   it (`collar.png`). Caveat: the cylinder keeps only the translation subgroup of the Bolza group, so this is a local
   model of the collar, not a computation on the closed surface.

9. *Instability below kappa_c at every spacing: DERIVED (full argument, sympy-checked, unrefereed).* For
   kappa <= 1/sqrt 2 the row is linearly unstable at every spacing l > 0, through the pairing mode. Write K = kappa^2,
   c_j = cosh(j l'), D = c - cos(2 alpha). Pair by pair (`pair_terms.py`):
   A_PU(pi) = -4 s^4 sum_{j odd} (c_j^2 + c_j - 1 - cos 2alpha)/((c_j - 1) D_j^2) < 0 for every kappa, and
   A_UP(pi) = 2 s^2 sum_{j >= 1} g_j with g_odd = ((2K-1)c - 1)/((c-1)D), g_even = ((4K-1)c - (4K^2-2K+1))/D^2.
   (i) Every g_j is strictly increasing in K, so for K <= 1/2 the sum is at most its value at K = 1/2,
   G(l') = sum_m [-1/((c_{2m-1} - 1)c_{2m-1}) + (c_{2m} - 1)/c_{2m}^2].
   (ii) G < 0 for every l' > 0: for l' >= 0.75 keep the j = 1 term and bound the even terms by sech, giving
   G < -4q/(1+q)^2 + 2q/(1-q) < 0 with q = e^{-2l'} < sqrt 5 - 2; for l' <= 1.2 bound the even terms by
   integral of sech = pi/(4 l') and use (cosh x - 1) cosh x < 4x/pi on (0, 1.2] (convex, negative at both ends).
   So A_UP(pi) < 0 < -A_PU(pi), and the q = pi symbol (zero diagonal) has the real exponent sqrt(A_UP A_PU) > 0.
   The even terms are not individually negative for 1/4 < K < 1/2 (`pair_sign.py`); only their sum against the
   odd terms is. Combined with item 3 this makes kappa_c = 1/sqrt 2 exact for the tanh model: no spacing is stable
   at or below it, and every spacing large enough is stable above it (the latter still leading-order only).
   Consequence for item 8: the Bolza-collar result for even N is now a theorem of the cylinder model for all N
   with an exact pairing mode (any even N, any h inside the collar), not only N <= 12 (superseded by item 10(b)).

10. *Exact symbol, explicit stability bound, all-mode instability: DERIVED (sympy-checked, unrefereed).*
    The full Bloch symbol at wavenumber q is [[i delta, B], [C, i delta]] with equal imaginary diagonal entries,
    B(q) = 2 sum_j [a(c_j) + b(c_j) cos qj], C(q) = 2 sum_j d(c_j)(1 - cos qj), d < 0, with
    a = s^2 [(3K-1)c^2 - (4K^2-K+1)c + 2K^2]/((c-1)D^2), b = s^2 [K c^2 + (1-3K)c + 2K^2-2K+1]/((c-1)D^2).
    So the row is stable iff F(q) = sum_j a_j + b_j cos(qj) >= 0 for all q. This matches the supercell growth
    rates to 7 digits and the thresholds to 6-7 digits up to kappa = 0.999 (the finite-difference lattice sum of
    item 6 gave only 4-5). The critical q is pi in every case computed (NUMERICAL).
    (a) b > 0 always and a - b = s^2((2K-1)c - 1)/((c-1)D), so for kappa > 1/sqrt 2 the row is stable whenever
    cosh l' >= 1/(2 kappa^2 - 1): arc spacing l_A = arccosh(1/(2k^2-1))/sqrt(1-k^2), 1.2-2.3 times the true l*,
    -> 2 as kappa -> 1. With item 9, kappa_c = 1/sqrt 2 is now exact for the tanh model on both sides.
    (b) a + b = s^2((4K-1)c - (4K^2-2K+1))/D^2 < 0 for kappa <= 1/2, so F(q) < 0 for every q: every Bloch mode
    grows. Hence in the collar model no Bolza necklace of any N >= 2, odd or even, is stable inside the systole
    collar (kappa <= sqrt 2 - 1 < 1/2). By (a) a necklace is stable once tanh^2 h >= (1 + 1/cosh(L/N))/2
    (h >= 1.22 for N = 2; computed threshold 1.06).
    `rows_boundary.png` now shows the proved bound l_A (dashed).

**Figures.** `rows/rows_boundary.png` (one boundary for all curve types), `rows/collar.png` (Bolza systole).

**What this says.** On H^2 the stabilising agent is the geodesic curvature of the curve carrying the vortices, with a
sharp minimum kappa_c = 1/sqrt 2, not the ambient curvature by itself. The flat-plane statement "straight rows are
unstable" survives on H^2 for every curve flatter than kappa_c. Written into the tanh note as the subsection "Rows on
curves of constant geodesic curvature" at the end of the large-N section (Prop. 22, Fig. 1, closed-surface remark; existing
theorem numbers unchanged), with entries in `docs/results.md`. The note's earlier heuristic ("screening explains the
stability") is corrected there.

**Open.** A proof that the critical mode is always q = pi (this would make l*(kappa) the exact zero of F(pi));
the same exact-symbol argument for other kernels (Q_1); the analogous family on the sphere (latitude circles only); why 45 degrees (the hypercycle at distance arsinh 1).
