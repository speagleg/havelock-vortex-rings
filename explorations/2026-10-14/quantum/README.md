# quantum/ — quantised vortex rows on H^2: phonons, a flat band at kappa_c, zero-point fluctuations

Point-vortex dynamics is first order, so a vortex's two coordinates are canonically conjugate. In the strip coordinates
of ../rows/ the Darboux pair is (u, Y = -cot phi), with {u_i, Y_i} = 1 and area form du dY = dA/R^2. Quantising gives
[u_i, Y_j] = i eps delta_ij. This is the lowest Landau level (LLL) on the hyperbolic plane (Comtet 1987). For vortices
in a superfluid, eps = 1/(2 pi n R^2), with n the areal number density and R the curvature radius: one quantum state
per particle of the fluid (Haldane & Wu 1985). Harmonic (Gaussian) theory throughout; tanh model.

| Script | What it does |
|---|---|
| `symbol.py` | Exact Bloch symbol [[i delta, B], [C, i delta]] of ../rows/exact_symbol.py, vectorised, with its q -> pi and q -> 0 expansions |
| `row_phonons.py` | The quantum phonon branch; energetic vs spectral threshold; velocities at the critical point -> `row_phonons.out` |
| `flat_band.py` | sympy: leading-order quadratic Hamiltonian, its determinant, the flat band and the bond-bending identity at kappa_c -> `flat_band.out` |
| `zero_point.py` | Zero-point fluctuations as Lindemann ratios; the logarithmic divergence at threshold -> `zero_point.out` |
| `plot_dispersion.py` | `dispersion.png` |
| `landau.py` | Landau coefficients of the zigzag transition from the exact lattice energy (fixed and relaxed ensembles) -> `landau.out` |
| `tricritical.py` | The tricritical curvature kappa_t where a4 = 0 on the threshold, and a6 there -> `tricritical.out` |
| `zigzag_dynamics.py` | Exact 2-vortex-cell dynamics just past threshold vs the Landau bound -> `zigzag_dynamics.out`, `zigzag_dynamics.png` |
| `plot_landau.py` | `landau.png` |

**1. One chiral branch** (DERIVED; `row_phonons.py`). Per Bloch wavenumber the quadratic form of H - V I (where
I = sum Y generates translation along u) is h(q) = [[-C/s^2, -i delta], [i delta, B s^2]]. There is one degree of
freedom per vortex, so there is one branch over the whole zone:
    E(q) = eps [ sqrt(-B(q) C(q)) + delta(q) ],   q in (-pi, pi].
It is chiral because delta is odd in q: modes moving one way along the row cost more than modes moving the other way.
The ground-state covariances do not depend on delta, which only shifts energies.

**2. Spectral and energetic thresholds coincide** (NUMERICAL). A quantum ground state needs E(q) >= 0 for all q
(no negative-energy, Krein-type modes). That is stronger than spectral stability. Near q = pi at threshold,
E ~ v|k| + delta'(pi) k, and near q = 0, E ~ c0|q| + delta'(0) q. Both ratios |delta'(pi)|/v and |delta'(0)|/c0 are
below 1 for every kappa tested (0.72-0.995), so E > 0 throughout the spectrally stable region. Every spectrally
stable row is an energetically stable ground state. Unlike odd-N rings, rows have no Krein window.

**3. The critical point is a Doppler-shifted relativistic soft mode** (DERIVED structure, NUMERICAL values). At l*,
E(pi + k) = v|k| + delta'(pi) k: a linearly dispersing (z = 1) mode seen from a frame drifting at speed delta'(pi).
The drift ratio |delta'(pi)|/v grows from 0.03 (kappa = 0.995) to 0.93 (0.72). As kappa -> 1/sqrt 2,
1 - |delta'|/v ~ 2 e^{-l'}, so one of the two movers stops.

**4. At kappa_c the quadratic Hamiltonian is degenerate at every wavenumber: a flat band** (DERIVED, exact sympy;
`flat_band.py`). To leading order in e^{-l'} (nearest neighbours, w_j = 2 sin phi_j sin phi_{j+1} e^{-(u_{j+1}-u_j)}):
    det h(q) = 16 e^{-2l'} s^6 (2 kappa^2 - 1)(1 - cos q).
- *Spectral vs energetic.* The q = pi spectral criterion of Prop. 22 only says the zigzag mode softens first. The
  energetic statement is stronger: h(q) is positive definite for every q when kappa > 1/sqrt 2, and singular for
  every q at kappa = 1/sqrt 2.
- *Mechanical meaning of 45 degrees.* At kappa = 1/sqrt 2 the energy is exactly w * sum_j (delta ln w_j)^2. Every
  deformation that keeps each nearest-neighbour bond weight fixed (fixed hyperbolic bond lengths, to this order)
  costs zero energy at fixed impulse: the chain bends freely. Above kappa_c bending costs energy; below, it releases
  energy.
- *Quantum meaning.* At kappa_c the branch is E(q) = sqrt 8 s^3 e^{-l'} (|sin q| - sin q). It is identically zero on
  0 < q < pi, and the other half propagates one way only (`dispersion.png`). Near kappa_c on the stable side, that half
  of the zone has bandwidth ~ e^{-2l'}, against e^{-l'} for the other half. A flat band means the harmonic ground state
  is nearly degenerate, so anharmonic terms (and, at T > 0, thermal fluctuations) decide what the row does.
- Beyond leading order the degeneracy is lifted with the sign that makes kappa = 1/sqrt 2 itself unstable (../rows/
  item 9).

**5. Zero-point fluctuations: quantum effects need few fluid particles per spacing** (NUMERICAL; `zero_point.py`).
With N_l = n (l R)^2 (fluid particles per spacing squared), the Lindemann ratios are gamma^2 = Lambda / N_l:
- transverse: Lambda_perp = 0.001-0.3;
- neighbour along-row: Lambda_par = 0.007-0.07 (kappa = 0.75-0.98, l/l* = 1.001-3).

The transverse ratio diverges logarithmically at the threshold. The slope d Lambda_perp / d ln(1/(l/l* - 1)) =
sqrt(2|C(pi)|/B''(pi))/(8 pi^2 l*^2) = 0.0218 at kappa = 0.9, matched by the direct integrals to 4 digits. So:
- *Ordinary superfluids:* quantum fluctuations wash out the 45 degree threshold only in a window of width
  ~ exp(-N_l gamma_L^2 / 0.022) around l*. For N_l ~ 10^2-10^4, typical of vortex spacings in BECs and helium films,
  the window is invisible.
- *Strongly quantum vortex matter:* N_l of a few to ~10 particles per vortex, as near the predicted melting of
  rotating-BEC vortex lattices, is where the harmonic picture fails and the row's quantum phase diagram becomes a real
  question.

**6. The zigzag transition is first order over most of the boundary; a tricritical point at kappa_t ~ 0.92-0.96**
(NUMERICAL, 30-digit lattice sums; `landau.py`, `tricritical.py`, `landau.png`).
Setup:
- *Order parameter.* The staggered amplitude eta, with sublattices at Y_A = Y0 + eta and Y_B = Y0 - eta. Swapping the
  sublattices sends eta -> -eta.
- *Landau function.* The exact lattice energy per cell, f(eta) = a2 eta^2 + a4 eta^4 + a6 eta^6.
- *Consistency.* a2 = B(pi) s^2 to 10 digits, so a2 vanishes exactly at l*.
- *Ensembles.* a4 at fixed spacing and impulse (what a closed periodic row conserves), and a4 with the long-wavelength
  spacing and impulse modes relaxed at fixed conjugate forces (the effective quartic of the local field theory). The
  relaxation matrix is positive definite everywhere tested.

On the threshold curve, a4 changes sign:
- *fixed:* a4 < 0 for kappa <= 0.90, > 0 for kappa >= 0.95; kappa_t = 0.92027 (l* = 1.3235).
- *relaxed:* a4 < 0 up to kappa = 0.95, > 0 for kappa >= 0.98; kappa_t = 0.96020 (l* = 1.1178).
- a6 > 0 at both points, so these are genuine tricritical points of the Landau theory.

Below kappa_t (most of the boundary, including the whole approach to 1/sqrt 2) the transition is subcritical. There is
no nearby zigzag minimum: at kappa = 0.8, f(eta) falls monotonically past a small barrier. Just above l* the row is
a local minimum separated by a barrier from lower-energy configurations, in which the two sublattices move apart
towards opposite sides of the ideal boundary. Above kappa_t (near the horocycle end) the transition is an ordinary
supercritical pitchfork to a zigzag row of amplitude eta^2 = -a2/(2 a4).

*Dynamical check* (`zigzag_dynamics.py`). In a 2-vortex cell the exact dynamics is the fixed-ensemble Landau problem:
eta is conjugate to the staggered u-shift xi, and orbits follow contours of H(eta, xi). Just past threshold
(l = 0.98 l*), the orbit's largest |eta| equals the quartic bound sqrt(-a2/a4) to 2e-4 at kappa = 0.98
(0.81304 vs 0.81288) and stays below it at 0.94. At kappa = 0.9 and 0.8 (a4 < 0) the orbits go further (1.43, 1.27)
and close only through higher-order and xi terms.

*Global caveat.* At fixed impulse the straight row is never the global energy minimum: W > 0, so spreading the
sublattices to the boundary drives the energy to 0, below the row's. All statements here are local, about the
neighbourhood of the row.

**7. Quantum reading of section 6** (interpretation, not computed). In the LLL, eta's conjugate is xi, with stiffness
h_uu(pi) = -C(pi)/s^2 > 0. Integrating xi out gives a relativistic (Doppler-tilted, z = 1) eta^4 - eta^6 theory in
1+1 dimensions. Its generic phase structure:
- *above kappa_t:* a continuous transition, naturally in the 1+1D Ising class, as for the quantum zigzag transition of
  ion chains;
- *below kappa_t:* a first-order transition, with the row metastable and decaying by quantum tunnelling (bubble
  nucleation);
- *at kappa_t:* the tricritical Ising point, whose continuum theory has central charge 7/10 and emergent spacetime
  supersymmetry (Friedan, Qiu & Shenker 1984).

So a quantised hypercycle row is a candidate realisation of the supersymmetric tricritical Ising point, tuned by the
geodesic curvature of the curve. Caveats:
- the location kappa_t is a mean-field (Landau) value and depends on the ensemble;
- quantum fluctuations shift it, although only weakly for many fluid particles per vortex (section 5);
- the coupling of eta^2 to the gapless acoustic phonons (the relaxed-vs-fixed difference) could change the
  universality class;
- the chiral (Doppler) term is harmless only while |delta'(pi)| < v (section 2).

**Not yet done.** A field-theory treatment of the eta^2-phonon coupling (does it preserve the tricritical point?). The
anharmonic lifting of the flat band near kappa_c. Finite temperature, where the flat
band makes thermal fluctuations large.

References: Comtet, Ann. Phys. 173, 185 (1987) (Landau levels on H^2); Haldane & Wu, PRL 55, 2887 (1985) (vortex
quantisation in superfluids); Shimshoni, Morigi & Fishman, PRL 106, 010401 (2011) (quantum zigzag transition);
Sinova, Hanna & MacDonald, PRL 89, 030403 (2002) and Cooper, Wilkin & Gunn, PRL 87, 120405 (2001) (quantum vortex
matter, Tkachenko modes, melting); Friedan, Qiu & Shenker, PRL 52, 1575 (1984) and Phys. Lett. B 151, 37 (1985)
(unitarity and supersymmetry of the tricritical Ising model).
