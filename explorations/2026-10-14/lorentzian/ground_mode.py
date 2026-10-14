"""Why the critical acceleration 1/sqrt 2 does not depend on the field's dimension Delta (sympy).
AdS2 scalar of dimension Delta = m (mass^2 = m(m-1)); global ground mode psi = sin^m(sigma), frequency m.
For the static observer at sigma = alpha (acceleration a = cos alpha, lapse 1/sin alpha), with proper derivatives
d_rho = sin(alpha) d_sigma and proper frequency m sin(alpha):
   radial "pressure"   P := phi_rho^2 + phi_tau^2 - m(m-1) phi^2  =  m phi^2            (phi = psi, phi_tau = m s psi)
   gradient term       a d_rho(phi^2)                             =  2 m a^2 phi^2
The pairing-mode bracket of the vortex row (ground-mode term of F(pi)) is  a d_rho(phi^2) - P = m phi^2 (2a^2 - 1):
both pieces scale with m = Delta, so the sign change sits at 2a^2 = 1 for every Delta.
Usage: python ground_mode.py"""
import sympy as sp
sg, m = sp.symbols('sigma m', positive=True)
psi = sp.sin(sg)**m; s, k = sp.sin(sg), sp.cos(sg)
phi_rho = s*sp.diff(psi, sg); phi_tau = m*s*psi
P = sp.simplify(phi_rho**2 + phi_tau**2 - m*(m-1)*psi**2)
grad = sp.simplify(k*s*sp.diff(psi**2, sg))
print('P / phi^2              =', sp.simplify(P/psi**2))
print('a d_rho(phi^2) / phi^2 =', sp.simplify(grad/psi**2), '  (a = cos sigma)')
print('bracket / (m phi^2)    =', sp.simplify((grad - P)/(m*psi**2)), ' = 2 a^2 - 1 = cos 2 sigma')
# cross-check against the bracket used in modes_general.py: 2 s k psi psi' + (mu - s^2 n^2) psi^2 - s^2 psi'^2 at n = m
br = 2*s*k*psi*sp.diff(psi, sg) + (m*(m-1) - s**2*m**2)*psi**2 - s**2*sp.diff(psi, sg)**2
print('matches modes_general bracket:', sp.simplify(br - (grad - P)) == 0)
