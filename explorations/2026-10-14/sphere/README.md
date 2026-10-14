# sphere/ — the K = +1 member: latitude rings on the sphere, and all three geometries in one picture

On the sphere the only homogeneous curves are the latitude circles, all closed, so there are no rows, only rings. Kernel:
the sphere's Green's function, W = -ln sin(d/2) (with a uniform compensating background).

| Script | What it does |
|---|---|
| `sphere_rings.py` | Direct linear stability of an N-ring at height z = cos(theta) (analytic field, mpmath Jacobian) |
| `sphere_blocks.py` | sympy: pair second derivatives and the exact Fourier blocks per mode m |
| `sphere_formula.py` | Closed form of the blocks and the stability thresholds -> `sphere_formula.out` |
| `plot_three_geometries.py` | `three_geometries.png`: sphere, plane and hyperbolic plane in (kappa, spacing) |

**Results (classical; reproduced here).**

1. *Closed-form blocks.* With z = cos(theta), S = sin(theta), mode m has
   A_m = [(N-1)(1+z^2) - m(N-m)]/(2 S^4), B_m = m(N-m)/2 and no cross term, so lambda^2 = -A_m B_m. The ring is stable
   iff (N-1)(1 + cos^2 theta) >= floor(N^2/4): N = 2, 3 at every latitude; cos^2 theta >= 1/3, 1/2, 4/5 for N = 4, 5, 6;
   N = 7 only at the pole (the planar Thomson heptagon, marginal); N >= 8 never. This is Polvani & Dritschel (1993),
   as stated in Laurent-Polz, Montaldi & Roberts (2004, Thm 5.2); the closed forms are checked against the exact block
   sums to 40 digits and against the direct Jacobian (`sphere_rings.py`, thresholds z^2 = 0.33333, 0.5004, 0.8003).
2. *In geodesic curvature.* With kappa = cot(theta), the criterion reads kappa^2 (2(N-1) - M) >= M - (N-1), M = floor(N^2/4).
   N = 4: 2 kappa^2 - 1 >= 0, i.e. kappa >= 1/sqrt 2. At that threshold the square ring is a face of the inscribed cube
   (cos^2 theta = 1/3). N = 5: kappa >= 1; N = 6: kappa >= 2.
3. *Observation.* The sphere's N = 4 criterion, 2 kappa^2 - 1 >= 0, has the same form as the large-spacing criterion of
   hypercycle rows on H^2 (cos 2 alpha = 2 kappa^2 - 1 >= 0, ../rows/ items 3, 12). Both concern the pairing mode
   (m = N/2 = 2 on the sphere, q = pi on H^2). Whether this is more than a coincidence is not known; on H^2 the
   criterion comes from the lowest Poschl-Teller mode, on the sphere from the two neighbour distances of a square.
4. *Three geometries* (`three_geometries.png`). In (kappa, arc spacing):
   - sphere: stable rings exist only near the poles and only for N <= 6;
   - plane: no length scale, rings stable iff N <= 7;
   - H^2 (tanh model): rings of every N are stable once the spacing is large enough, and rows on hypercycles are
     stable iff kappa > 1/sqrt 2 and l >= l*(kappa).
   Positive curvature confines stability to small, dense rings; negative curvature with its screened Green's function
   makes large spacing stabilising.

**Sources.** L. M. Polvani & D. G. Dritschel, J. Fluid Mech. 255 (1993) 35-64; F. Laurent-Polz, J. Montaldi &
M. Roberts, "Point vortices on the sphere: stability of relative equilibria" (2004), arXiv:math/0402430, Thm 5.2.
