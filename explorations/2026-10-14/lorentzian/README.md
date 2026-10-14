# lorentzian/ — what the hyperbolic-plane results say about Lorentzian spacetimes and gravity

The hyperbolic plane enters relativity in two exact ways. Both are checked here.

1. **As Euclidean AdS2.** Wick-rotating one coordinate turns it into a Lorentzian spacetime. This is a dictionary for
   static (Euclidean) configurations, and the place where black-hole throats live.
2. **As the spatial slice of an expanding universe.** Point-vortex dynamics carries over exactly here, with real
   Lorentzian time. The cleanest case is the Milne wedge of flat Minkowski space.

| Script | What it does |
|---|---|
| `ads2_dictionary.py` | sympy: Wick rotation of the strip metric to global AdS2, accelerations, Deser-Levin temperature, kernels as AdS2 propagators -> `ads2_dictionary.out` |
| `thermal_rings.py` | Certified ring thresholds rewritten as temperatures -> `thermal_rings.out` |
| `growth_map.py` | Growth rates of every hypercycle row from the exact symbol; Milne freeze-out -> `growth_map.png`, `growth_map.out` |
| `correlator_identity.py` | The pair functions as worldline correlators of the AdS2 field (checked to 1e-30 for Delta = 1, 2) -> `correlator_identity.out` |
| `ground_mode.py` | Why 1/sqrt 2 is the same for every Delta: two ground-mode identities -> `ground_mode.out` |

## 1. Euclidean AdS2 dictionary (exact; `ads2_dictionary.py`)

| Hyperbolic plane (this repo) | AdS2 / gravity |
|---|---|
| strip coordinates z = e^{u+i phi}, metric (du^2 + dphi^2)/sin^2 phi | global AdS2, (-dt^2 + dphi^2)/sin^2 phi, after u -> i t: the row direction u is Euclidean global time |
| hypercycle phi = alpha, geodesic curvature kappa = cos alpha < 1 | Euclidean worldline of the static global observer at sigma = alpha, proper acceleration a = kappa (AdS units), no horizon |
| horocycle, kappa = 1 | Poincare-static observer, a = 1: the extremal (T = 0) horizon |
| circle of radius r, kappa = coth r > 1 | Rindler-AdS2 observer, a > 1: Euclidean thermal circle, T = 1/(2 pi sinh r) = sqrt(a^2 - 1)/(2 pi) (Deser-Levin) |
| kernel Q_nu(cosh d), Green's function of -Delta + nu(nu+1) | EAdS2 bulk-to-bulk propagator of a scalar with m^2 = Delta(Delta-1), Delta = nu + 1 |
| tanh model, -ln tanh(d/2) = Q_0 | massless scalar, Delta = 1 |
| N vortices on a curve | N point sources along the Euclidean worldline, spacing = Euclidean proper time |

The near-horizon geometry of an extremal or near-extremal charged black hole is AdS2 x S2, so the same table describes
the throat:
- BH-static observers near an extremal horizon are the horocycles.
- Near-extremal black holes give the thermal circles, centred on the Euclidean horizon.
- The global observers (hypercycles) are the horizonless ones.

**The results, translated.**

- *Critical acceleration.* A row of sources along the worldline of a static global-AdS2 observer is stable iff
  a > a_c = 1/(sqrt 2 L_AdS) and the Euclidean-time spacing is at least l*(a). This is exact for Delta = 1 (proved, with
  the single-curve certificate), and at leading order for every Delta (the Q_nu family, ../rows/ item 12).
  In Lorentzian terms, a_c is the observer at sigma = pi/4, whose redshift factor 1/sin(sigma) equals sqrt 2.
- *Temperature thresholds for rings* (`thermal_rings.out`). A ring of N sources on the thermal circle is stable iff
  T <= T*_N. The values are T*_N = 0.2860, 0.2101, 0.1628, ... for N = 8, 9, 10, ... (AdS units); N <= 7 is stable at
  every T. The products N T*_N = 2.288, 1.891, ..., 1.192 (N = 64) decrease to 1/l* = 1.18497. Equivalently, a ring of
  N sources is stable iff the Euclidean-time spacing beta/N exceeds a critical value l_N that increases to 0.8439 L_AdS as
  N -> infinity (the extremal limit).
- *Zigzag reading.* The criterion F(pi) >= 0 is a statement about the energy Hessian alone: the stiffness of the
  staggered (zigzag) transverse displacement of the chain. So the 45 degree rule says that, in EAdS2, a chain of like
  scalar charges held on a worldline of acceleration a <= 1/sqrt 2 always gains energy from zigzagging.

**What does not carry over.** Point-vortex motion is a first-order flow on the Euclidean plane, and Wick rotation
does not turn it into a Lorentzian process. The dictionary rows above are statements about static configurations and
their energy, not about real-time evolution.

## 2. Real-time dynamics: expanding universes and the Milne wedge

Consider a thin superfluid (or ideal-fluid) film comoving with an expanding universe whose spatial slices are
hyperbolic, ds^2 = -c^2 dt^2 + a(t)^2 h_{H^2}. With uniform density on each slice, the vortical part of the flow is
incompressible on the comoving slice. By the 2D conformal invariance of the Laplace Green's function, the stream
function of point vortices is Gamma G_h, independent of a(t). The comoving velocity is a(t)^{-2} times the unit-H^2
velocity. Hence:

**The comoving vortex dynamics is exactly the tanh-model dynamics on the unit hyperbolic plane, in the time
tau = int dt / a(t)^2.**

Every result in ../rows/ and ../../tanh/ therefore holds verbatim for comoving vortex arrays in an open universe:
- the 1/sqrt 2 rule;
- the thresholds l*(kappa), in comoving curvature units;
- the ring thresholds;
- the Bolza-collar statements, if the slice is a compact hyperbolic surface.

In comoving units these criteria are frozen: they do not drift as the universe expands.

**Milne: Lorentzian and flat.** The future light cone of an event in Minkowski space, foliated by unit hyperboloids,
is ds^2 = -c^2 dt^2 + (ct)^2 h_{H^2}, so a = ct. Two consequences:
- *Finite total evolution.* tau = int dt/(c t)^2 is finite, so a pattern formed at time t0 has only Delta tau =
  1/(c^2 t0) of evolution left. An unstable mode with rate lambda_max can grow by at most
  exp(lambda_max Gamma / (2 pi c^2 t0)).
  - `growth_map.png` maps lambda_max over all hypercycle rows (Gamma = 2 pi, unit curvature). It is up to ~10^2 for
    dense rows and ~10^-1 at spacing 3; it vanishes at the boundary l*(kappa) and is zero beyond it.
  - In every case computed, the fastest-growing mode is the pairing (zigzag) mode q = pi. That q = pi decides
    stability is proved (../rows/ item 11); that it is also the fastest-growing mode is not.
  - For superfluid 4He, Gamma/c^2 = 1.1e-24 s, so any macroscopic pattern is frozen in: in the Milne wedge an
    instability has too little tau left to act.
  - More generally the total tau is finite whenever a(t) grows faster than t^{1/2}: matter- and curvature-dominated
    open universes freeze vortex patterns, while radiation domination (a ~ t^{1/2}) is the marginal case.
- *Rapidity space.* The unit hyperboloid is the space of 2+1 velocities, and hyperbolic distance is relative rapidity.
  A comoving point is the velocity of a free particle leaving the origin event. So vortex arrays in a Milne-expanding
  film are arrays in velocity space. A hypercycle is the set of velocities at fixed rapidity distance from a boost
  line, and the critical one sits at rapidity arsinh 1 = 0.8814 from it.

**Caveats.**
- This is the incompressible point-vortex limit on the comoving slice: low Mach number relative to the Hubble flow, and
  vortex speeds small compared with c. Vortex inertia and relativistic corrections are not included.
- The film must be two-dimensional and lie in a hyperbolic slice. In 3+1 open FRW that means a totally geodesic H^2
  inside the H^3 slices; in 2+1 it is the slice itself.

## 3. The criterion as field theory on AdS2: correlators, thermal images, and why Delta drops out

**Pair functions are worldline correlators** (`correlator_identity.py`). Take the field phi with Euclidean propagator
W = 2 pi G, two points on the static worldline at Euclidean separation x, proper derivatives d_tau (along) and d_rho
(normal), and a = cos(alpha). Then:

    b = <d_rho phi(x) d_rho phi(0)>,     d = <d_tau phi(x) d_tau phi(0)>,
    a = mu W + d + 2 a <d_rho phi(x) phi(0)>      (the field equation, Delta = d_tau^2 + d_rho^2 - a d_rho).

**The row is an image sum.** A row with spacing l' in Euclidean global time is one insertion in thermal global AdS2
at beta = l'. The static observer then sees the local Tolman temperature T_loc = 1/l, with l the arc spacing. The
Bloch wavenumber q is a twist e^{iq} around the thermal circle, and the pairing mode q = pi is the Z2-twisted
(antiperiodic, "fermionic") sector. With S_+ the periodic and S_- the antiperiodic image sum, excluding the vacuum
(j = 0) term:

    2 F(pi) = S_+[ mu <phi phi> + <d_tau phi d_tau phi> + 2a <d_rho phi phi> ]  +  S_-[ <d_rho phi d_rho phi> ].

These are thermal (respectively twisted) excesses of local field bilinears at the observer. In this language the
stability boundary l*(a) is a critical local temperature, T_c(a) = 1/l*(a). It runs from 0 at a = 1/sqrt 2 to
1/0.8439 = 1.185 L_AdS^{-1} at a = 1.

**Why the critical acceleration does not depend on Delta** (`ground_mode.py`). At low temperature both image sums
are dominated by the global ground mode, psi = sin^Delta(sigma), with energy Delta. For that mode, seen by the static
observer, two identities hold:

    radial stress  P = phi_rho^2 + phi_tau^2 - m^2 phi^2 = Delta phi^2,
    gradient term  a d_rho(phi^2)                         = 2 Delta a^2 phi^2.

The criterion is gradient minus stress, Delta phi^2 (2 a^2 - 1). Both terms scale with Delta, so the sign change sits
at a^2 = 1/2 for every field dimension. This is the field-theory form of the lowest-mode identity in ../rows/ item 12,
and it is why the H^2 rule is kernel-universal while the sphere's N = 4 value is not (../sphere/ item 3).

## Open

- A Lorentzian process whose linearisation is the vortex symbol itself; for example, vortices in a 2+1 relativistic
  superfluid on AdS2 x S1 (dual to charged particles, by particle-vortex duality).
- Section 3 gives a_c a field-theory meaning at the level of thermal and twisted bilinears. Still open: whether it
  marks a change in something directly measurable, such as an Unruh-DeWitt detector with a derivative coupling to
  d_rho phi in the twisted state.
- Rings on the thermal circle with Delta != 1 (the Q_1 ring thresholds).
