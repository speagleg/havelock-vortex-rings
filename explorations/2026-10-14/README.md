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
| `qpi_proof.py` | Item 11: generating-function identities, the critical mode is q = pi, the threshold series, odd-N necklaces -> `qpi_proof.out` |
| `horo_limit.py` | Item 11: the alpha -> 0 limit of the threshold series is the horocycle criterion of Cor. 20 -> `horo_limit.out` |
| `pair_generic.py` | Item 12: pair functions for any kernel w(cosh d) |
| `kernel_coeffs.py` | Item 12: series coefficients of the pair functions for the tanh and Q_1 kernels (Cauchy integral) |
| `pt_modes.py` | Item 12: the coefficients are Poschl-Teller mode weights (checked for Q_0 and Q_1) |
| `modes_general.py` | Item 12: mode-sum criterion F_nu(pi) for the kernels Q_nu; Q_1 thresholds and horocycle value -> `modes_general.out` |
| `dual_series.py` | Item 13: theta-function form, exact modular dual series, elementary closed-form threshold -> `dual_series.out` |
| `plot_closed_form.py` | `closed_form.png`: exact threshold vs the closed form and the proved bound |
| `single_root_cert.py` | Item 14: Arb certificate that each kappa has exactly one threshold (needs `pip install python-flint`; ~5 s) -> `single_root_cert.out` |
| `q1_instability_cert.py`, `q1_mode_consistency.py` | Item 15: Arb certificate of instability below 1/sqrt 2 at every spacing for the Q_1 kernel (Delta = 2) -> `.out` |

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

11. *The critical mode is always the pairing mode; exact threshold series: DERIVED (exact sympy identities).*
    In z = e^{-x}, c = cosh x, the pair functions are Chebyshev-type generating functions:
    b = 2 s^2 sum_n n cos^2(n alpha) z^n,  d = -2 s^2 sum_n n sin^2(n alpha) z^n,  a = d + sin^2(2 alpha)/(2(c - cos 2alpha))
    (so b - d = s^2/(c - 1)). Every coefficient of b is >= 0, so sum_j b(c_j) cos(qj) is a positive combination of
    Poisson kernels P(r_n, q) = sum_j r_n^j cos(qj), r_n = e^{-n l'}, each strictly decreasing in q on (0, pi).
    Hence F(q) is strictly decreasing: the row is stable iff F(pi) >= 0, and a necklace of N is stable iff
    F(2 pi floor(N/2)/N) >= 0. Summing the geometric series,
    F(pi) = sum_n [(sin 2a sin 2na - 2n s^2 sin^2 na) r_n/(1 - r_n) - 2n s^2 cos^2 na r_n/(1 + r_n)],
    a one-variable series converging like e^{-n l'}. Its root agrees with the supercell thresholds to 1e-6
    (l*(0.8) = 2.2554766253, l*(0.9999) = 0.8530342); F(pi) changes sign exactly once in l' for every kappa
    tested (NUMERICAL). The Bolza necklace thresholds h* = 1.0574 ... 2.1070 (N = 2..12, both parities) from this
    criterion agree with `collar.out` to all printed digits. Horocycle limit: alpha -> 0 at fixed arc spacing turns
    the series into int_0^inf [(2 sin 2t - 2t sin^2 t)/(e^{lt} - 1) - 2t cos^2 t/(e^{lt} + 1)] dt
    = -eta^2 coth^2(eta)/4 + eta coth(2 eta) - 1/4, eta = pi/l, which vanishes exactly when
    (coth eta - 1/eta)^2 = 2 tanh(eta)/eta: the note's horocycle criterion (Cor. 20), l* = 0.8439085312925335.
    So the hypercycle family limits onto the horocycle row by an independent route.

12. *Why it works: a mode expansion, valid for the whole screened family Q_nu. DERIVED (identities checked).*
    In (u, phi) the operator -Delta + nu(nu+1) is sin^2(phi)[-d_u^2 - d_phi^2 + nu(nu+1)/sin^2(phi)] and separates.
    Poschl-Teller modes psi_n = sin^{nu+1}(phi) C^{(nu+1)}_{n-nu-1}(cos phi), eigenvalue n^2, n >= nu+1, so
    W_nu = 2 pi G_nu = 2 pi sum_n e^{-n|u|} psi_n(phi) psi_n(phi0)/(2 n N_n). Then, with w_n = 2 pi e^{-n x}/(2 n N_n),
    b = s^2 sum w_n psi_n'(alpha)^2 >= 0,  d = -s^2 sum n^2 w_n psi_n(alpha)^2 <= 0,
    a = sum w_n [2 s k psi_n psi_n' + (nu(nu+1) - n^2 s^2) psi_n^2].
    (nu = 0, psi_n = sin(n phi), gives item 11.) So for every nu: same symbol structure, critical mode q = pi, threshold =
    zero of a one-variable mode series F_nu(pi) (matches direct lattice sums to 6e-12). The lowest mode
    psi ~ sin^m(phi), m = nu + 1, contributes exactly m sin^{2m}(alpha) cos(2 alpha): the leading term is ~ cos(2 alpha) e^{-m l'}
    for every nu. This answers "why 45 degrees": the hypercycle phi = alpha meets the ideal boundary at angle alpha, and the
    critical one meets it at 45 degrees, where the lowest mode's along-curve (psi'^2) and normal (psi^2) weights balance.
    Q_1 (nu = 1): roots of F_1(pi) give l* = 5.7296, 4.0856, 3.0055, 1.9820, 1.6859, 1.4949 at kappa = 0.72, 0.75, 0.8,
    0.9, 0.95, 0.99 (supercell: 4.09, 3.01, 1.98, 1.69, 1.49), one sign change each. New: the horocyclic row with the
    Q_1 kernel has critical spacing l* = 1.4491260 (zero of the alpha -> 0 integral; kappa = 0.9999 gives 1.4496).
    For nu = 1, 2, F_nu(pi) < 0 at all 3000 grid points with kappa <= 1/sqrt 2 (NUMERICAL; proved only for nu = 0).

13. *Closed-form threshold: DERIVED (theta functions + modular transformation; checked to 1e-36).*
    Using r/(1+r) = r/(1-r) - 2r^2/(1-r^2), F(pi) = sin2a S(a; r) - 2s^2[E(r) - E(r^2)] + 2s^2 C(a; r^2), i.e. log-derivatives
    of Jacobi theta_1 at nomes e^{-l'/2}, e^{-l'} and the Eisenstein sum E(r) = sum n r^n/(1-r^n). The modular transformation
    gives the exact dual series (eta = pi a/l', Q1 = e^{-4pi^2/l'}, Q2 = e^{-2pi^2/l'}):
      F(pi) = F0 - (2pi/l') sum Q1^n/(1-Q1^n)[sin2a sinh(4pi n a/l') - 4pi n s^2/l']
                 - (2pi^2 s^2/l'^2) sum n Q2^n/(1-Q2^n)[1 + cosh(2pi n a/l')],
      F0 = (pi sin2a/(2l')) coth(2 eta) - (pi^2 s^2/(4 l'^2)) coth^2(eta) - a sin2a/l' - cos2a/4.
    Corrections are O(e^{-2pi(pi - a)/l'}) and, for kappa >= 1/sqrt 2, both positive, so F < F0 (F0 < 0 => unstable).
    The ELEMENTARY equation F0 = 0 is the threshold to 1e-4 at kappa = 0.8, 2e-7 at 0.85, 1e-11 at 0.9, 4e-21 at 0.95;
    it degrades only near kappa_c (20% at 0.72, where l' is large and the original series is the right tool). As a -> 0
    at fixed arc spacing F0 is exactly the horocycle criterion. This is the hypercycle generalisation of the note's
    Cor. 20: one closed formula for the whole family between kappa ~ 0.8 and the horocycle. Figure `closed_form.png`.

14. *The threshold is a single curve: CERTIFIED (Arb ball arithmetic, python-flint).* For every kappa in (1/sqrt 2, 1),
    F(pi) has exactly one root in l', crossed transversally (unstable below, stable above). Method: show dF/dl' > 0 at
    every zero with alpha < pi/4 (zeros with alpha >= pi/4 are excluded by item 9), so each F(alpha, .) can only cross
    upwards; F < 0 at small spacing and F > 0 for cosh l' >= 1/(2k^2-1) (item 10). Coverage of all zeros:
    l' >= 6 analytic with an explicit O(r^3) remainder (-P2 in [1.95, 2.01] vs needed 0.32); 0.25 <= l' <= 6 by 15,225
    adaptive boxes (sign of G, or dG/dl' > 0) with Lambert-series tail bounds; l' < 0.25 in (alpha, eta = pi alpha/l')
    via the dual series of item 13: eta >= 7 elementary (F <= F0 < 0), eta <= 2.5 elementary (F0 >= h(eta) - O(eta^2 t)
    beats R), and 113 boxes in between. Consequences: the stable set is exactly {l >= l*(kappa)} with l* real-analytic;
    l'* = l* sin(alpha) is strictly decreasing in kappa (every g_j increases with kappa^2, item 9).

15. *kappa_c = 1/sqrt 2 is exact for Delta = 2 too: CERTIFIED (Arb).* For the Q_1 kernel every hypercycle row with
    kappa <= 1/sqrt 2 is unstable at every spacing (F_1(pi) < 0 on alpha in [pi/4, pi/2], all l'). Three regions:
    - l' >= 4.5: F_1 = (8/3) s^4 cos 2a r^2 + c3 r^3 + E, with c3 = (16/3) s^4 (6s^4 - 9s^2 + 2) <= -4/3 and r Eb = 1.047 < 4/3.
    - 0.05 <= l' <= 4.5: 6,781 Arb boxes on the mode form, with an exact polynomial tail constant.
    - l' <= 0.05: elementary bounds. b is positive and decreasing, and d < 0, both from the mode expansion. The
      remainder e = a - d = 2 Q1(C) + 2 k^2 (C-1)|Q1'(C)| is positive and decreasing, with integral <= 9 pi/4.
      Dropping a log term gives d <= D_up. Then F_1 <= (s^2/l'^2)[P(l') + (9 pi/4) l'/s^2] < 0 (Arb, margin 0.004).
    Closed forms and mode expansions agree to 6e-30 (`q1_mode_consistency.py`). Together with item 12 (leading-order
    stability above 1/sqrt 2 for every nu), the critical curvature is exact for the massless (Delta = 1) and the
    Delta = 2 fields.

**Figures.** `rows/rows_boundary.png` (one boundary for all curve types), `rows/collar.png` (Bolza systole),
`rows/closed_form.png` (exact threshold vs the elementary closed form).

**What this says.** On H^2 the stabilising agent is the geodesic curvature of the curve carrying the vortices, with a
sharp minimum kappa_c = 1/sqrt 2, not the ambient curvature by itself. The flat-plane statement "straight rows are
unstable" survives on H^2 for every curve flatter than kappa_c. Written into the tanh note as the subsection "Rows on
curves of constant geodesic curvature" at the end of the large-N section (Prop. 22, Fig. 1, closed-surface remark; existing
theorem numbers unchanged), with entries in `docs/results.md`. The note's earlier heuristic ("screening explains the
stability") is corrected there.

## sphere/ — the K = +1 member and all three geometries

See `sphere/README.md`. Latitude rings on the sphere reproduce the classical criterion (Polvani & Dritschel 1993):
stable iff (N-1)(1 + cos^2 theta) >= floor(N^2/4), with closed-form Fourier blocks. In geodesic curvature the N = 4
threshold is exactly kappa = 1/sqrt 2 (the square is a face of the inscribed cube), the same value as the hypercycle
criterion on H^2. Tested (sphere/README item 3): both criteria contain the same geometric stiffness 2 kappa^2 + K from the
conformal background, but on H^2 it is the whole criterion (kernel terms cancel at q = pi), while on the sphere it
competes with the flat interaction, giving the mirror criterion 2 kappa^2 - K >= 0. Changing the kernel decides it:
screening moves the sphere square's threshold (0.691 ... 0.871 for mu = -0.5 ... 6), but leaves the H^2 value at 1/sqrt 2.
So the two share an ingredient but not a mechanism. `sphere/three_geometries.png`
puts sphere, plane and H^2 in one (kappa, spacing) picture: positive curvature confines stability to small dense rings
(N <= 6), the plane allows N <= 7, negative curvature makes large spacing stabilising for every N.

## lorentzian/ — Lorentzian spacetimes and gravity

See `lorentzian/README.md`. (1) Exact Euclidean-AdS2 dictionary: the strip coordinates are Euclidean global AdS2, a
hypercycle of curvature kappa is the worldline of a static observer with acceleration a = kappa, a horocycle is the
extremal (T = 0) horizon, a circle is a thermal circle with the Deser-Levin temperature sqrt(a^2 - 1)/(2 pi), and the
kernels Q_nu are AdS2 scalar propagators with Delta = nu + 1. So the 1/sqrt 2 rule is a critical acceleration
a_c = 1/(sqrt 2 L_AdS), and the ring thresholds are temperatures, N T*_N -> 1/l* = 1.18497. (2) Real-time dynamics:
comoving vortices in an open universe follow the H^2 dynamics exactly in tau = int dt/a^2, so all thresholds hold in
comoving units. In the Milne wedge of Minkowski space tau is finite, so instabilities freeze out (growth map
`lorentzian/growth_map.png`).

## channel/ — walls versus curvature

See `channel/README.md`. In strip coordinates the tanh kernel is exactly the Green's function of a flat channel with
walls at phi = 0, pi, so a hypercycle row and a channel row share the pair interaction. They differ only by the walls'
image self-term (1/2) ln sin(phi) = -sigma/2, which the H^2 geometric self-energy cancels, and by the mobility
sin^2(phi). Exact identity: B_H2(q) = sin^2(a) B_ch(q) + sin(2a) V_ch - cos(2a)/2, C_H2 = sin^2(a) C_ch. Result: the channel row
is never stable (25,600-point scan, min_q B_ch <= -1/2; at large spacing B_ch(pi) -> -csc^2(a)/2), while the H^2 row is
stable beyond l*(kappa) for kappa > 1/sqrt 2. Walls destabilise; constant negative curvature stabilises.

## realize/ — physical surfaces and their edges

See `realize/README.md`. The three curve types are the parallels of Minding's three K = -1 surfaces of revolution
(pseudosphere: horocycles; hyperboloid type: hypercycles; cone type: circles). On an embedded hyperboloid-type surface
with waist radius a the 45 degree parallel exists iff a <= R; for a = R it is the rim. On a finite film band, edges add
a self-term that destabilises large spacings, so the rule becomes a stability window l* <= l <= l_max, which needs
bands several curvature radii wide. On embedded Minding surfaces stable necklaces (N = 2, near the rim) were found
only for waists a <~ 0.09 R. A flat-space design (E x B plasma with B ~ csc^2 and a compensating potential) follows
from the channel identity.

**Open.** A pen-and-paper version of item 14 (the certificate is computer-assisted); instability below 1/sqrt 2 at every spacing for nu >= 2; a closed form for the Q_1 horocycle threshold 1.4491;
on the sphere, a configuration whose threshold is set by 2 kappa^2 + K alone (a sphere counterpart of the kernel-universal
H^2 rule).
