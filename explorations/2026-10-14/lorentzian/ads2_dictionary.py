"""Dictionary: the hyperbolic-plane vortex problem as Euclidean AdS2 (sympy checks).
(1) Strip coordinates z = e^{u + i phi}: metric (du^2 + dphi^2)/sin^2 phi. Wick rotation u -> i t gives
    (-dt^2 + dphi^2)/sin^2 phi, global AdS2 (radius 1). The row direction u is Euclidean global time.
(2) The curve phi = alpha (hypercycle) is the Euclidean worldline of the static observer at sigma = alpha; its proper
    acceleration in AdS2 equals the geodesic curvature kappa = cos(alpha) of the hypercycle.
(3) Circles of radius r (kappa = coth r > 1) are Euclidean thermal circles: Euclidean proper period 2 pi sinh r, so the
    local temperature is T = 1/(2 pi sinh r) = sqrt(kappa^2 - 1)/(2 pi): the Deser-Levin temperature of an observer with
    acceleration a = kappa in AdS2. Horocycles (kappa = 1) are T = 0 (Poincare horizon); hypercycles (kappa < 1) see none.
(4) Kernels: Q_nu(cosh d) solves (Delta_{H^2} - nu(nu+1)) Q = 0 away from the source: it is the EAdS2 propagator of a
    scalar with m^2 = nu(nu+1) = Delta(Delta - 1), Delta = nu + 1 (normalisable fall-off e^{-Delta d}). The tanh model,
    -ln tanh(d/2) = Q_0(cosh d), is the massless (Delta = 1) case.
Usage: python ads2_dictionary.py"""
import sympy as sp
u, t, phi, d, nu = sp.symbols('u t phi d nu', positive=True)
# (1) metric under Wick rotation
g_E = sp.diag(1, 1)/sp.sin(phi)**2
J = sp.Matrix([[sp.I, 0], [0, 1]])                 # du = i dt
g_L = sp.simplify(J.T*g_E*J)
print('(1) Wick-rotated strip metric:', g_L.tolist(), ' = (-dt^2 + dphi^2)/sin^2 phi (global AdS2)')
# (2) acceleration of the static worldline phi = alpha in g = (-dt^2 + dphi^2)/sin^2 phi:
#     a = |grad ln sqrt(-g_tt)| = sin(phi) * |d/dphi ln(1/sin phi)| = |cos phi|
a = sp.simplify(sp.sin(phi)*sp.Abs(sp.diff(sp.log(1/sp.sin(phi)), phi)))
print('    check a - |cos phi| on (0, pi/2):', sp.simplify((a - sp.cos(phi)).subs(phi, sp.Rational(7, 10))), sp.simplify(sp.sin(phi)*(sp.cos(phi)/sp.sin(phi)) - sp.cos(phi)))
print('(2) proper acceleration of the static observer at sigma = phi:', a, '(= geodesic curvature cos(alpha) of the hypercycle)')
# (3) circle in H^2: kappa = coth r, circumference 2 pi sinh r
r = sp.symbols('r', positive=True); kap = sp.coth(r)
T = 1/(2*sp.pi*sp.sinh(r))
print('(3) T - sqrt(kappa^2 - 1)/(2 pi) =', sp.simplify(T - sp.sqrt(kap**2 - 1)/(2*sp.pi)))
# (4) radial Laplacian on H^2: f'' + coth(d) f'  ;  Q_0 = -ln tanh(d/2),  Q_1 = C artanh(1/C) - 1, C = cosh d
lap = lambda f: sp.diff(f, d, 2) + sp.cosh(d)/sp.sinh(d)*sp.diff(f, d)
Q0 = -sp.log(sp.tanh(d/2)); C = sp.cosh(d); Q1 = C*sp.atanh(1/C) - 1
print('(4) Delta Q_0 =', sp.simplify(sp.expand(lap(Q0).rewrite(sp.exp))), ';  Delta Q_1 - 2 Q_1 =', sp.simplify((lap(Q1) - 2*Q1).rewrite(sp.log)))
print('    fall-off: Q_0 ~', sp.limit(Q0*sp.exp(d), d, sp.oo), 'e^{-d};  Q_1 ~', sp.limit(Q1*sp.exp(2*d), d, sp.oo), 'e^{-2d}  (Delta = 1, 2)')
