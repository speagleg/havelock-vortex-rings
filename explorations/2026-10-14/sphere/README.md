# sphere/ — the K = +1 member: latitude rings on the sphere, and all three geometries in one picture

On the sphere the only homogeneous curves are the latitude circles, all closed, so there are no rows, only rings. Kernel:
the sphere's Green's function, W = -ln sin(d/2) (with a uniform compensating background).

| Script | What it does |
|---|---|
| `sphere_rings.py` | Direct linear stability of an N-ring at height z = cos(theta) (analytic field, mpmath Jacobian) |
| `sphere_blocks.py` | sympy: pair second derivatives and the exact Fourier blocks per mode m |
| `sphere_formula.py` | Closed form of the blocks and the stability thresholds -> `sphere_formula.out` |
| `plot_three_geometries.py` | `three_geometries.png`: sphere, plane and hyperbolic plane in (kappa, spacing) |
| `n4_generic.py` | Pairing block of the square ring for an arbitrary kernel w(cos d) -> `n4_generic.out` |
| `sphere_screened.py` | Square-ring threshold for the screened sphere kernels (-Delta + mu) -> `sphere_screened.out` |
| `mirror_identity.py` | Both 1/sqrt 2 criteria in conformal language -> `mirror_identity.out` |

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
3. *Is the shared 1/sqrt 2 a coincidence? Tested; answer: a shared ingredient, but not a shared mechanism.*
   - *Shared ingredient.* Write either surface as e^{2 sigma(y)}|dw|^2 with the curve at y = y0, so kappa = -e^{-sigma} sigma'
     and K = -e^{-2 sigma} sigma''. The conformal background contributes the same normal stiffness to both,
     e^{-2 sigma}(2 sigma'^2 - sigma'') = 2 kappa^2 + K. On H^2 this is the whole pairing-mode criterion: at q = pi every
     kernel-dependent term cancels (../rows/ items 3, 12), so rows need 2 kappa^2 + K >= 0. In the sphere's block it is
     the term (N-1)(1 + cos^2 theta) = (N-1) sin^2 theta (2 kappa^2 + K) (`mirror_identity.py`).
   - *Different balance.* On the sphere that stiffness competes with the flat interaction, m(N - m). For the square,
     3 sin^2 theta (2 kappa^2 + 1) >= 4 with sin^2 theta = 1/(1 + kappa^2) gives 2 kappa^2 - K >= 0. This is the mirror
     image of the H^2 criterion 2 kappa^2 + K >= 0 under K -> -K, and both read 2 kappa^2 >= 1 because |K| = 1. At the
     sphere threshold the geometric stiffness is 2, not 0.
   - *Decisive test: change the kernel.* For an arbitrary kernel w(cos d) the square's pairing block is
     A_2 = 2[w'(c2) - w'(c1) + 4 z^2 w''(c2)], with c1 = z^2 and c2 = 2z^2 - 1 (`n4_generic.py`). It depends on the
     kernel's shape at the two neighbour distances, and only the Laplace kernel gives 3 z^2 = 1. With the screened sphere
     kernels (-Delta + mu), the square's critical curvature moves:
     kappa_4 = 0.6909, 0.7071, 0.7229, 0.7679, 0.8710 for mu = -0.5, 0, 0.5, 2, 6. The same numbers follow from the
     closed-form block (`sphere_screened.py`). On H^2 the corresponding screened family Q_nu keeps
     kappa_c = 1/sqrt 2 exactly at leading order for every nu (lowest mode gives m sin^{2m}(a) cos 2a; ../rows/ item 12).
   - *Verdict.* The H^2 value is geometric and kernel-universal. The sphere value is a property of the Laplace kernel
     and of N = 4. The two meet at 1/sqrt 2 because both criteria are built on 2 kappa^2 + K, with opposite signs of
     the curvature term, and |K| = 1. They agree only for the unscreened sphere kernel.
4. *Three geometries* (`three_geometries.png`). In (kappa, arc spacing):
   - sphere: stable rings exist only near the poles and only for N <= 6;
   - plane: no length scale, rings stable iff N <= 7;
   - H^2 (tanh model): rings of every N are stable once the spacing is large enough, and rows on hypercycles are
     stable iff kappa > 1/sqrt 2 and l >= l*(kappa).
   Positive curvature confines stability to small, dense rings; negative curvature with its screened Green's function
   makes large spacing stabilising.

**Sources.** L. M. Polvani & D. G. Dritschel, J. Fluid Mech. 255 (1993) 35-64; F. Laurent-Polz, J. Montaldi &
M. Roberts, "Point vortices on the sphere: stability of relative equilibria" (2004), arXiv:math/0402430, Thm 5.2.
