"""Exact pairing-mode (q = pi) stiffness of a hypercycle row, expanded in eps = exp(-l'), with the full tanh kernel
W = (1/2) ln((C+1)/(C-1)), C = 1 + [cosh(u_i - u_j) - cos(phi_i - phi_j)]/(sin phi_i sin phi_j), neighbours |j| <= 3.
At q = pi the 2x2 Bloch symbol has zero diagonal, so the row is pairing-stable iff A_UP(pi) >= 0 (A_PU(pi) < 0).
Result: A_UP(pi)/(4 sin^2 alpha) = (2k^2-1) eps + (8k^4-3) eps^2 + 4(8k^6-8k^4+2k^2-1) eps^3 + O(eps^4), k = kappa.
Checked against the exact lattice sum (pairing_root.AUP_pi): the residual after two terms is O(eps^3) with the
predicted coefficient (-4.05 measured vs -3.99 at kappa = 0.72).
Usage: python next_order.py   (~2 min)"""
import sympy as sp
a, x = sp.symbols('alpha x', positive=True)          # x = exp(l')
J = 3
ph = sp.symbols('p_m3 p_m2 p_m1 p0 p1 p2 p3')
idx = {j: ph[j+J] for j in range(-J, J+1)}

def C(j):
    cosh = (x**j + x**(-j))/2                      # pair (0, j): delta u = -j l', cosh is even
    return 1 + (cosh - sp.cos(idx[0]-idx[j]))/(sp.sin(idx[0])*sp.sin(idx[j]))

def Wp(Cv):
    return -1/(Cv**2-1)                            # dW/dC

ud0 = sp.sin(idx[0])**2*sum(Wp(C(j))*sp.diff(C(j), idx[0]) for j in range(-J, J+1) if j != 0)
sub = {idx[j]: a for j in range(-J, J+1)}
AUP_pi = sum(((-1)**abs(j))*sp.diff(ud0, idx[j]) for j in range(-J, J+1)).subs(sub)
eps = sp.symbols('epsilon', positive=True)
expr = sp.simplify(AUP_pi.subs(x, 1/eps))
ser = sp.series(expr, eps, 0, 4).removeO()
c1 = sp.simplify(ser.coeff(eps, 1)); c2 = sp.simplify(ser.coeff(eps, 2)); c3 = sp.simplify(ser.coeff(eps, 3))
print('A_UP(pi) =', sp.factor(c1), '* eps  +', sp.factor(c2), '* eps^2  + O(eps^3)')
k = sp.symbols('kappa', positive=True)
tok = lambda e: sp.factor(sp.simplify(sp.expand(sp.expand_trig(e/(4*sp.sin(a)**2))).subs(sp.cos(a), k).subs(sp.sin(a), sp.sqrt(1-k**2))))
print('in kappa:  A_UP(pi) / (4 sin^2 alpha) =', tok(c1), '* eps +', tok(c2), '* eps^2 +', tok(c3), '* eps^3 + O(eps^4)')
