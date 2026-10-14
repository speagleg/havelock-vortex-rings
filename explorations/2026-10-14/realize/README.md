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

**Flat-space realisation (from ../channel/).** The identity there shows what a flat system needs to reproduce the H^2
dynamics: the channel interaction, no wall-image drift, and mobility sin^2(phi). One candidate is a guiding-centre
(E x B) plasma of line charges in a slab with conducting walls. A magnetic field profile B(y) proportional to
csc^2(pi y/w) supplies the mobility, since the drift speed is E/B. An applied transverse potential proportional to
ln sin(pi y/w) cancels the image drift, which is the same for every like charge. This is a design idea; it has not
been simulated here.
