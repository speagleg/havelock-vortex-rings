# realize/ — physical surfaces that carry the H^2 vortex dynamics, and what their edges do

**Minding's surfaces.** The three kinds of curves of constant geodesic curvature in H^2 are exactly the parallels of
the three surfaces of revolution of constant Gaussian curvature K = -1 (Minding, 1830s):

| surface of revolution (K = -1) | quotient of H^2 by | its parallels are | kappa of a parallel |
|---|---|---|---|
| pseudosphere (Beltrami) | a parabolic translation | horocycles | 1 |
| hyperboloid type (has a waist geodesic) | a hyperbolic translation of length 2 pi a | hypercycles at distance h from the waist | tanh h |
| cone type | a rotation | circles | coth r |

So a ring of vortices around a parallel of such a surface is exactly one of the configurations in ../rows/: a
necklace on a hyperboloid-type surface is a hypercycle row (the collar model of ../rows/ item 8). In R^3 these surfaces
exist only as bands. The hyperboloid type with waist radius a (in units of the curvature radius R) reaches distance
He = arsinh(R/a) from the waist, where the parallels have curvature tanh He = 1/sqrt(1 + a^2/R^2). Hence:

- **The 45 degree parallel lies on the surface iff the waist radius is at most R.** For a = R the critical parallel,
  at h = arsinh 1 = 0.8814 R, is exactly the rim of the embedded surface. For a > R every parallel of the embedded
  surface has kappa < 1/sqrt 2, so no stable necklace exists, by the proved H^2 instability theorem (in the idealised
  model with no film edges).

| Script | What it does |
|---|---|
| `band.py` | Film on a finite band phi1 < phi < phi2 of H^2 with impermeable edges: strip Green's function plus the self-term S = (flat Robin term of the strip) + (geometric self-energy) |
| `band_scan.py` | Edge self-term stiffness, and row stability on a grid for bands of half-width Hb -> `band_scan.out` |
| `minding_necklaces.py` | Is any necklace (N >= 2) stable on an embedded hyperboloid-type surface? -> `minding_necklaces.out` |
| `flat_lemma.py` | The local equivalence lemma, checked by direct linearisation in a strip and a disk -> `flat_lemma.out` |
| `ring_blocks.py` | Exact Fourier blocks for rings in the disk (analytic second derivatives), and the exact stability window in S |
| `shear_design.py` | Trap and channel thresholds in background shear vs the H^2 values -> `shear_design.out`, `shear_design.png` |
| `havelock_check.py` | Independent sympy check of the unsheared (Havelock) thresholds, N = 2, 3 -> `havelock_check.out` |
| `nonlinear_sim.py` | Nonlinear runs: trap with passive background, channel with Couette shear, trap with a central column -> `nonlinear_sim.out`, `nonlinear_sim.png` |
| `active_background.py`, `ab_conv.py` | Trap with a self-consistently advected point-vortex background; convergence in N_b -> `ab_conv.out` |
| `ring_center.py` | Full (N+1)-vortex linearisation (mpmath, 50 digits) of ring + mobile central column vs H^2; Havelock g-intervals -> `ring_center.out` |
| `ring_center_n16.py` | N = 16: trap + mobile centre at g* vs H^2 across r0* -> `ring_center_n16.out` |
| `ring_center_df.py`, `ring_center_plot.py` | The Durkin-Fajans geometry (r0 = 0.25); g_lo vs g* -> `ring_center_df.out`, `ring_center.png` |
| `channel_symbol_shim.py` | Sigma* and Sigma_c from ../channel/channel_symbol.py |

**Results (NUMERICAL; edges modelled as impermeable, one choice of net circulation around the waist).**

1. *Edges replace the 45 degree rule by a window.* On a finite band the vortex self-term no longer cancels. Its stiffness
   B_self = d/dphi[sin^2(phi) S'(phi)] is negative at every point tested (e.g. -0.23, -0.04, -0.004 at kappa = 1/sqrt 2
   for half-widths 2, 3, 5). Pair couplings decay like e^{-l'} while the self-term does not, so every row becomes
   unstable again at large spacing. Stability survives only in a window l*(kappa) <= l <= l_max(kappa, Hb). The window
   is empty for half-width 2, opens near kappa = 0.9 at half-width 3, and widens at 5 (l in [3, 5] at kappa = 0.8,
   [2, 12] at kappa = 0.9).
2. *Embedded Minding surfaces: stable necklaces only for thin waists.* Scanning waist radius a, N = 2..12 and every
   parallel on the embedded surface, the first stable necklace appears only for a <~ 0.09 R. It is always N = 2,
   near the flaring rim, at kappa ~ 0.97-0.998 and arc spacing ~ 1.0-1.35 R. For a >= 0.1 R none was found.

## Flat laboratory realisations: one knob, and a Penning-trap experiment

The idea from ../channel/ (an E x B plasma with B ~ csc^2 and a log-sin potential) would reproduce the H^2 dynamics
exactly, but it is not buildable as stated. A z-directed field that varies across the slab needs a volume current
(curl B != 0), and a potential proportional to ln sin needs a background charge density proportional to csc^2. For
LINEAR stability much less is needed.

**1. Local equivalence lemma** (exact; `flat_lemma.py`). Take vortices with pair kernel W, mobility m depending only
on the transverse coordinate y, and self term U(y). A uniformly translating row (or rigidly rotating ring) at y0 has
linearisation
    J = mu(y0) [[ -H_ty, -H_tt ], [ H_yy + S I, H_yt ]],    S = U'' + (mu'/mu)(W_y + U'),
with mu = m (strip) or m/r (disk) and W_y the transverse interaction force. So two systems with the same W have
proportional linear spectra if they have the same S (only this direction is shown). Every other property of m and U is irrelevant. The H^2 tanh
kernel is exactly the Green's function of a flat channel (../channel/) and of a flat disk with a conducting wall,
so:
- *Channel.* The H^2 row is the channel row plus a uniform background vorticity, i.e. a Couette shear of rate
  Sigma*(alpha, l') = 4 cot(a) sum_j W_phi(j) + csc^2(a)/2, of the sign opposite to the vortices. Uniform vorticity
  is invariant under any incompressible flow in a fixed domain, so the set-up is an exact 2D Euler solution (the
  equivalence with H^2 is linear only). Eigenvalues agree
  to 3e-11, while the same row without shear is unstable (growth 0.7-3.0).
- *Penning-Malmberg trap* (uniform B, conducting circular wall: Havelock's problem). The H^2 ring of N at hyperbolic
  radius rho is the trap ring at r0 = tanh(rho/2) plus a background whose shear r dOmega_b/dr at the ring equals
  s*(N, r0). A rigid rotation (uniform background) changes nothing; only the shear counts. Eigenvalues agree to 1e-6
  (finite differences).

**2. Thresholds in the knob** (`shear_design.py`; exact blocks, `ring_blocks.py`). Each Fourier mode's condition is
linear in S. A trap ring is stable iff s >= s_c(N, r0), and a channel row iff Sigma >= Sigma_c = -B_ch(pi). On the
whole grid only lower bounds occur.
- At s = 0 this is Havelock's problem. Stable for r0 < r_H: p = r_H^2 = 0.213740, 0.321281, 0.329840, 0.346101,
  0.299121 (N = 2-6), and never for N >= 7. N = 2 and 3 were checked independently with an exact sympy Jacobian
  (`havelock_check.out`); I could not retrieve Havelock's own table to compare.
- The crossings s* = s_c reproduce the note's H^2 ring thresholds r0* = tanh(rho*/2) to 6 digits: 0.259531, 0.407571,
  0.497411, 0.609657 for N = 8, 10, 12, 16.
- Channel: Sigma* crosses Sigma_c iff kappa > 1/sqrt 2 (the 45 degree rule seen as a shear threshold), and
  Sigma_c >= 0.68 on the grid.

**3. The buildable version: a central column** (`ring_center.py`, `ring_center_df.py`, `ring_center.png`). A central
vortex of strength g (ring column = 1) is, for the ring, a background of shear exactly 2g/r0^2, with no vorticity at the
ring. So the H^2 ring corresponds to Havelock's ring-plus-centre problem with
                    g*(N, r0) = s*(N, r0) r0^2 / 2.
The centre can move, but its displacement couples only to ring modes k = +-1. Results:
- *The full (N+1)-vortex linearisation (50 digits)* has the same stability verdict and the same largest growth rate as
  the H^2 ring in all six cases tested. Only the k = +-1 eigenvalues differ.
- *In every case checked, Havelock's lower limit g_lo(N, r0) is set by the modes the centre does not touch.* It equals s_c r0^2/2 from the
  exact blocks to 4-5 digits (N = 5, 7, 8 tested). The upper limit (12.2 for N = 8) lies far above g*.
- **Hence, in every case checked, the H^2 ring theorem is a statement about Havelock's lower limit: the H^2 ring of N
  is stable iff g*(N, r0) >= g_lo(N, r0).** (Rigorous for a fixed centre, by the lemma; with a mobile centre the k = +-1
  sector is checked case by case, not proved.) The trap ring with g = g*(N, r0) changes from unstable to stable across r0*(N)(1 -+ 0.003)
  for N = 8, 9, 10, 16 (with the mobile centre). In all 9 cases scanned (N = 3-12), the Havelock interval contains g*
  exactly when the H^2 ring is stable.
- *Validation against experiment-era theory.* For N = 6 at r0 = 0.25 the code gives -0.2478 <= g <= 6.2498. Durkin &
  Fajans (Phys. Fluids 12, 289, 2000) quote -0.25 <= g <= 6.25.
- *The experiment exists.* Durkin & Fajans built rings with a central column in a Penning-Malmberg trap (3 T, wall
  radius 2 cm, ring radius 0.25 of the wall), could inject 0 <= g <= 4, and measured the lower limit for N = 7-10.
  At their radius:

  | N | g_lo(N, 0.25) | g*(N, 0.25) | H^2 ring |
  |---|---|---|---|
  | 7 | 0.0014 | 0.3956 | stable (N <= 7 is stable at every radius) |
  | 8 | 0.5002 | 0.4622 | unstable (r0* = 0.2595 lies just above 0.25) |
  | 9 | 1.0001 | 0.5289 | unstable |
  | 10 | 1.7500 | 0.5956 | unstable |

  So their N = 8 lower-limit measurement sits 8% in g from the H^2 equivalence point. A laboratory analogue of the
  H^2 threshold would vary r0 across r0*(N) with g = g*(N, r0) (for N = 8: g* ~ 0.46-0.55 for r0 in [0.25, 0.27]),
  inside their injection range. What such an experiment tests is Havelock's point-vortex linear theory at the
  parameters where it coincides with the H^2 ring problem; the H^2 statement itself is mathematics (the lemma) and
  needs no experiment. Feasibility caveats:
  - g* - g_lo is only about -0.04 to +0.04 over r0 in [0.25, 0.27], comparable to their quoted 10% control of column
    strengths;
  - the columns have finite cores (radius / spacing ~ 0.13 at r0 = 0.25), so corrections of order (a/d)^2 ~ 2% are of
    the same size;
  - growth rates vanish at threshold, so near r0* instability may not show within their ~300-rotation lifetimes.

**4. Nonlinear runs** (`nonlinear_sim.py`, `nonlinear_sim.png`; perturbation 1e-4, 100 rotation periods; at this
amplitude they mainly confirm the linear analysis). N = 8 at
r0 = 0.5 is Havelock-unstable but H^2-stable. Radial spread:
- passive background with s = s*: 4.0e-4, identical to the H^2 run;
- no background: 0.64;
- s = s_c - 2: 0.62;
- s = s_c + 2: 7e-4;
- mobile central column at g* = 2.22: 4.3e-4;
- g_lo - 0.1: 0.11;
- g_lo + 0.1: 9e-4.

Channel (kappa = 0.924, l' = 1): with Sigma* the row tracks H^2 (2.4e-4); without shear, 1.3.

**5. Active smooth background: inconclusive** (`active_background.py`, `ab_conv.out`). I represented the peaked
background zeta_b = A(1 - r^2) by N_b weak blob vortices advected self-consistently. The ring spread after one period
falls with N_b (0.22, 0.081, 0.013 for N_b = 200, 600, 1500 at blob radius 0.03), so graininess dominates at accessible
N_b. Two physical concerns remain.
- The needed background is strong: its total circulation is 4.4 times the columns' for N = 8 at r0 = 0.5.
- Its gradient at the ring is stirred by the columns: the capture radius sqrt(r0/s*) ~ 0.17 is comparable to the
  half-spacing.

The central column avoids both, because it puts no vorticity at the ring.

Literature: Havelock (1931), Phil. Mag. 11, 617; Durkin & Fajans (2000), Phys. Fluids 12, 289; Kurakin (2010), Regul.
Chaotic Dyn. 15(1), doi:10.1134/S1560354710010028 (nonlinear stability for p <= p*_n).

**Original flat-space idea (from ../channel/; superseded by the section above).** The identity there shows what a flat system needs to reproduce the H^2
dynamics: the channel interaction, no wall-image drift, and mobility sin^2(phi). One candidate is a guiding-centre
(E x B) plasma of line charges in a slab with conducting walls. A magnetic field profile B(y) proportional to
csc^2(pi y/w) supplies the mobility, since the drift speed is E/B. An applied transverse potential proportional to
ln sin(pi y/w) cancels the image drift, which is the same for every like charge. This is a design idea; it has not
been simulated here.
