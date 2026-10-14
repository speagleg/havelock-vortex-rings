"""The two 1/sqrt 2 thresholds in one language. Write the surface as e^{2 sigma(y)} |dw|^2 with the carrying curve y = y0.
Then the geodesic curvature of the curve is kappa = -e^{-sigma} sigma' and the Gaussian curvature is K = -e^{-2 sigma} sigma''.
  H^2 hypercycle rows (y = phi, sigma = -ln sin phi): large-spacing criterion cos 2a >= 0  <=>  2 kappa^2 + K >= 0
                                                       <=>  2 sigma'^2 - sigma'' >= 0.
  Sphere N = 4 ring (y = Mercator u, sigma = -ln cosh u, z = tanh u): criterion 3 z^2 - 1 >= 0  <=>  2 kappa^2 - K >= 0
                                                       <=>  2 sigma'^2 + sigma'' >= 0.
The two criteria are mirror images under K -> -K (sigma'' -> -sigma''); with |K| = 1 both read 2 kappa^2 >= 1.
Usage: python mirror_identity.py"""
import sympy as sp
y = sp.symbols('y', positive=True)
def kK(sigma):
    s1, s2 = sp.diff(sigma, y), sp.diff(sigma, y, 2)
    return sp.simplify(-sp.exp(-sigma)*s1), sp.simplify(-sp.exp(-2*sigma)*s2), s1, s2
# hyperbolic plane, strip coordinate phi
kap, K, s1, s2 = kK(-sp.log(sp.sin(y)))
print('H^2:    kappa =', kap, ' K =', K)
print('        2 kappa^2 + K =', sp.simplify(2*kap**2 + K), ' = cos(2 phi):', sp.simplify(2*kap**2 + K - sp.cos(2*y)) == 0)
print('        e^{-2 sigma}(2 sigma\'^2 - sigma\'\') - (2 kappa^2 + K) =', sp.simplify(sp.sin(y)**2*(2*s1**2 - s2) - (2*kap**2 + K)))
# sphere, Mercator coordinate u (z = cos theta = tanh u)
kap, K, s1, s2 = kK(-sp.log(sp.cosh(y)))
z = sp.tanh(y)
print('sphere: kappa =', kap, ' K =', K)
print('        (3 z^2 - 1)/(2 kappa^2 - K) =', sp.simplify((3*z**2 - 1)/(2*kap**2 - K)), '(positive factor sech^2 u)')
print('        e^{-2 sigma}(2 sigma\'^2 + sigma\'\') - (2 kappa^2 - K) =', sp.simplify(sp.cosh(y)**2*(2*s1**2 + s2) - (2*kap**2 - K)))
