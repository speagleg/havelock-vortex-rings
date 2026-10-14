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
| `phonon_coupling.py` | Couplings of eta^2 to the long-wavelength sector, Larkin-Pikin shift, critical vs sound velocities -> `phonon_coupling.out` |
| `pinned.py` | Tricritical point when the spacing is pinned (long-wavelength sector gapped) -> `pinned.out` |
| `chain_model.py`, `chain_check.py` | The pinned row as a 1D quantum chain, its pair interactions compressed into exponentially decaying separable couplings, and a classical check against the exact lattice energy -> `chain_check.out` |
| `dmrg_chain.py` | iDMRG (TeNPy 1.1) of the chain, two-site unit cell |
| `scan_ising.py`, `scan_tau.py` | Scans in spacing (weak quantum strength) and in quantum strength tau -> `scan_ising.out`, `scan_tau_coarse/fine/window.out` |
| `central_charge.py` | Finite-entanglement scaling S vs ln xi -> `central_charge_0.98_0.95_*.out` |
| `plot_dmrg.py` | `dmrg_scan.png` |
| `flat_nonlinear.py` | Nearest-neighbour model: spacings eliminated exactly, flat band to all orders, thermal factorisation, single-site potential and barrier, checked against the full tanh symbol -> `flat_nonlinear.out` |
| `flat_lifting.py` | sympy: relaxed transverse stiffness of the exact symbol through e^{-2l'} (the lifting of the flat band) -> `flat_lifting.out` |
| `plot_site_potential.py` | `site_potential.png` |

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

Correction (section 8): a free row is NOT a clean candidate for the tricritical Ising point, because the acoustic
phonons couple relevantly there. A row with pinned spacing is.

**8. Coupling to the acoustic phonons; pinned rows** (scaling from the literature; numbers NUMERICAL;
`phonon_coupling.py`, `pinned.py`).

*Scaling.* For a critical point coupled to elastic (acoustic) modes, the coupling is relevant iff
alpha_q = 2 - nu(d + z) > 0 (Larkin-Pikin, Sak; quantum form in Chandra et al. and arXiv:2307.03062).
- *1+1D Ising line* (nu = 1): alpha_q = 0, so the coupling is marginal. Its sign is not determined here.
- *Tricritical Ising point* (energy operator of dimension 1/5, nu = 5/9): alpha_q = 8/9 > 0, so the coupling is
  relevant, with eigenvalue 2 - 2/5 = 8/5. The supersymmetric fixed point is unstable to the eta^2-phonon coupling.

*Numbers* (`phonon_coupling.out`). The Larkin-Pikin shift a4_fixed - a4_relaxed is comparable to a4 itself
(0.3-1 times |a4| for kappa = 0.94-0.995), so the coupling is not weak. The sound is 2.8-3.6 times faster than the
critical mode (c0/v). With the chiral drifts, though, the slowest sound mover is only 1.0-1.6 times faster than the
fastest critical mover, so neither the static (relaxed) nor the instantaneous (fixed) limit is controlled.

*Pinned spacing* (`pinned.out`). A periodic potential along the curve that pins each u_j (for vortices in a BEC, an
optical pinning lattice; Tung, Schweikhard & Cornell 2006):
- It gaps the whole long-wavelength sector: u and Y are conjugate, so Y0 pairs with the pinned u0 into a gapped
  oscillator.
- It leaves a2, the threshold l*, and the fixed-spacing Landau function untouched (the pinning potential vanishes on
  the zigzag configurations).
- Below the gap, Y0 follows adiabatically, which shifts a4 locally: a4_pinned = a4_fixed - g_Y^2/(2 M_YY).
- The tricritical point is then at kappa_t = 0.95991 (l* = 1.1193), almost the relaxed value 0.96020. Nearly all of
  the Larkin-Pikin shift comes from the impulse channel; the strain channel moves kappa_t by only 3e-4.

With no gapless mode left, the eta^4-eta^6 theory is the standard Z2 one. So a pinned quantised row at kappa ~ 0.96 is
a genuine mean-field candidate for the 1+1D tricritical Ising point. A free row has the relevant phonon coupling;
how much, and over what length scales, only a numerical quantum computation (DMRG on the LLL chain) can say.

**9. First quantum many-body numerics: quantum melting of the zigzag in a pinned row** (NUMERICAL, iDMRG;
`dmrg_chain.py`, `dmrg_scan.png`).

*Model.* In the stiff-pinning limit each vortex keeps one quantum coordinate y_j = Y_j - Y0, and the chain is
    H = sum_j (tau/2)(-d^2/dy_j^2) + sum_{i<j} W(|i-j| l'; Y0 + y_i, Y0 + y_j) - V sum_j y_j,
with tau = k eps^2 (pinning stiffness times ħ_eff^2). The -V term is the uniform counterflow that holds the row at
rest in the pinning lattice; without it the pinning cannot stop the drift along u.
- *Compression.* The pair interactions decay like 1/r^2 out to r ~ 1/l' and exponentially beyond, so they are
  compressed into 34 exponentially decaying separable couplings (MPO bond dimension 36) with pair error <= 4e-6.
- *Classical check* (`chain_check.out`, with the DMRG compression settings). The compressed model reproduces the exact
  zigzag energy profile to ~1e-6 per cell for |eta| <= 0.5 and to 8e-6 (1-10% of f) at eta = 0.8. Its curvature error
  corresponds to a threshold shift of ~1e-4 in l/l*.
- *Local space.* A sinc-DVR grid of 12 points.

*Results at kappa = 0.98, l = 0.95 l* (classically zigzag-ordered), bond dimension chi = 16:*
- *Weak quantum strength* (tau = 5e-4, scan in spacing): the quantum zigzag amplitude follows the classical Landau
  value (0.565 vs 0.59 at l = 0.99 l*). Quantum effects are negligible, as section 5 predicts.
- *Raising tau melts the order:* |eta| = 0.72, 0.68, 0.62, 0.55, 0.51, 0.43, 0.07, 0.02, 0, 0 at
  tau = 0.004-0.010. The correlation length peaks sharply (xi = 4 -> 19.5) at tau_c ~ 0.0079. At this resolution it
  looks like a continuous quantum phase transition.
- *Central charge: NOT determined.* At tau = 0.0079 the finite-entanglement slopes 6 dS/d ln xi are 1.03, 0.72 and
  0.28 for chi = 12-16, 16-20 and 20-24 (all-point fit 0.81). Meanwhile eta -> 0 with chi, so that point is slightly on
  the disordered side and the slopes measure a crossover. The test of c = 1/2 here (and c = 7/10 at kappa_t) needs the
  pseudo-critical point tau_c(chi) at each bond dimension and chi >~ 64. At about 2 s per sweep at chi = 16 on 2 cores
  (cost ~ chi^3), that is beyond this machine.

**10. The flat band beyond quadratic order and at finite temperature: independent sites and a conjugate-curvature
barrier** (classical; DERIVED for the nearest-neighbour model, which is the exact leading order in e^{-l'}; checked
against the full tanh model; `flat_nonlinear.py`, `flat_lifting.py`, `site_potential.png`).

*(a) Spacings eliminated exactly.* Take H = sum w_j with w_j = 2 sin(phi_j) sin(phi_{j+1}) e^{-D_j}, and
F = H + V sum cot(phi_j), at fixed period sum D_j = N l'. Minimising over the spacings D_j makes every bond weight
equal, so
    F_rel(phi) = 2N e^{-l'} (prod_j sin^2 phi_j)^{1/N} + V sum_j cot(phi_j)
(numerical minimisation agrees to 1e-15). The first term depends only on the MEAN of ln sin(phi_j), so every
non-uniform deformation sees a sum of independent single-site potentials:
- *Flat to all orders.* The relaxed transverse "band" is exactly flat in q, not just at quadratic order. Its
  stiffness is k_perp = 2P csc^2(alpha) cos(2 alpha), with P = 2 sin^2(alpha) e^{-l'} the bond energy; this matches
  the numerical Hessian to 7 digits.
- *At kappa_c all N - 1 non-uniform modes go soft together.*

*(b) Finite temperature: the transverse coordinates are exactly independent.* Work in the pressure ensemble with
Liouville measure du dY = du dphi/sin^2(phi). Integrating each spacing gives
    int exp(-beta(c e^{-D} + P D)) dD = Gamma(beta P)(beta c)^{-beta P}
(checked), with c_j = 2 sin phi_j sin phi_{j+1}. The joint weight therefore factorises into independent sites:
    p(phi) ~ exp(-beta h_T(phi)),   h_T(phi) = 2(P + T) ln sin(phi) + V cot(phi).
The thermal row is a gas of independent vortices, each in its own potential well. Temperature enters only as
P -> P + T, an entropic stiffening from the measure.

*(c) The 45 degree rule as a saddle-node of the single-site potential.* At T = 0, h'(phi) = 0 <=> sin(2 phi) = V/P.
- The row is the minimum (phi = alpha).
- The barrier sits at phi = pi/2 - alpha, the CONJUGATE hypercycle, with curvature sqrt(1 - kappa^2).
- They merge at the self-conjugate curvature kappa_c = 1/sqrt 2.
- *Barrier per vortex:* Delta E = 2P[ln cot(alpha) - cos(2 alpha)]
  = 4(1 - kappa^2) e^{-l'} [ (1/2) ln(kappa^2/(1 - kappa^2)) - (2 kappa^2 - 1) ]
  = (16/3) P (pi/4 - alpha)^3 + ... ~ (32 sqrt 2/3)(kappa - kappa_c)^3 e^{-l'}.
- *Beyond the barrier* h falls without bound towards the far side of the geodesic. That is the same unzipping
  direction as the subcritical Landau branch (section 6).

So near kappa_c a thermal row loses vortices one at a time, at an Arrhenius rate ~ exp(-Delta E/T) with a barrier
that vanishes cubically. This is a fourth reading of the 45 degree rule, after the boundary angle, the ground-mode
balance and free bending: minimum and barrier sit at conjugate curvatures kappa and sqrt(1 - kappa^2), and 1/sqrt 2
is where they coincide.

*(d) How next-nearest neighbours lift the flat band* (exact series, `flat_lifting.out`). For the full tanh symbol,
the relaxed transverse stiffness (spacings relaxed: the Schur complement Sigma(q) = s^2(B|C| - delta^2)/|C|, which is
>= 0 exactly when E(q) >= 0) is
    Sigma(q) = 4 s^4 [ (2 kappa^2 - 1) e^{-l'} + (8 kappa^4 - 1 + 2 cos q) e^{-2l'} ] + O(e^{-3l'}).
- *Width.* The band acquires a cos q dispersion of width 16 s^4 e^{-2l'} only at second order. Numerically the mean
  tends to 4 s^4 (2 kappa^2 - 1) e^{-l'} (ratio 4.0002) and the width to 16 s^4 e^{-2l'}.
- *At q = pi* it reproduces the threshold series of ../rows/ item 7.
- *At kappa_c,* Sigma = (1 + 2 cos q) e^{-2l'}, negative exactly for |q| > 2 pi/3. That band of modes makes kappa_c
  itself unstable.
- *Away from threshold* the sites are therefore nearly independent: transverse thermal fluctuations are spatially
  white, with <delta Y_j^2> ~ T/Sigma, up to O(e^{-l'}) correlations.

**Not yet done.** The central charge at the pinned-row critical line and at kappa_t (larger chi, or a faster DMRG
code). Which quantum strengths are physical depends on the platform: tau = k eps^2 with eps = 1/(2 pi n R^2). The fate
of the free (unpinned) row at the marginal/relevant phonon coupling (DMRG). A rigorous (non-leading-order) version of
section 10, and the full-model thermal statistics beyond the nearest-neighbour factorisation (Monte Carlo).

References: Comtet, Ann. Phys. 173, 185 (1987) (Landau levels on H^2); Haldane & Wu, PRL 55, 2887 (1985) (vortex
quantisation in superfluids); Shimshoni, Morigi & Fishman, PRL 106, 010401 (2011) (quantum zigzag transition);
Sinova, Hanna & MacDonald, PRL 89, 030403 (2002) and Cooper, Wilkin & Gunn, PRL 87, 120405 (2001) (quantum vortex
matter, Tkachenko modes, melting); Friedan, Qiu & Shenker, PRL 52, 1575 (1984) and Phys. Lett. B 151, 37 (1985)
(unitarity and supersymmetry of the tricritical Ising model); Larkin & Pikin, Sov. Phys. JETP 29, 891 (1969); Sak, PRB 10,
3957 (1974) (critical points on compressible lattices); Chandra, Coleman, Continentino & Lonzarich, PRR 2, 043440 (2020) and
arXiv:2307.03062 (quantum versions); Tung, Schweikhard & Cornell, PRL 97, 240402 (2006) (optical pinning of BEC vortices); Pollmann, Mukerjee, Turner & Moore, PRL 102, 255701 (2009) (finite-entanglement scaling);
Hauschild & Pollmann, SciPost Phys. Lect. Notes 5 (2018) (TeNPy).
