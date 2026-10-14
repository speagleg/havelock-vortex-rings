"""Closed-form pair contributions to the pairing-mode (q = pi) symbol of a hypercycle row (tanh model).
Pair (0, j), with Delta u = j l', c := cosh(j l'), s := sin(alpha), k := cos(alpha).
A_UP(pi) = sum_{j != 0} a(c) + (-1)^j b(c),   A_PU(pi) = sum_{j != 0} (1 - (-1)^j) d_self(c).
Used by pair_sign.py and instability_proof_check.py.
Usage: python pair_terms.py"""
import sympy as sp
p0, pj, du, al = sp.symbols('p0 pj du alpha', real=True)
cc = sp.symbols('c', positive=True)
C = 1 + (sp.cosh(du) - sp.cos(p0-pj))/(sp.sin(p0)*sp.sin(pj))
Wp = -1/(C**2-1)
# u_dot_0 from pair j: sin^2 p0 * W'(C) dC/dp0 ;  phi_dot_0 from pair j: -sin^2 p0 * W'(C) dC/du0  (du = u0 - uj)
T = sp.sin(p0)**2*Wp*sp.diff(C, p0)
S = -sp.sin(p0)**2*Wp*sp.diff(C, du)
eq = {p0: al, pj: al}
def clean(e):
    e = e.subs(eq)
    e = e.subs(sp.cosh(du), cc).subs(sp.sinh(du)**2, cc**2-1)
    return sp.factor(sp.simplify(e))
a_self = clean(sp.diff(T, p0))          # d/dphi_0 of pair term
b_other = clean(sp.diff(T, pj))         # d/dphi_j
print('A_UP(pi) pair term:  a + (-1)^j b')
print('  a =', a_self)
print('  b =', b_other)
# A_PU: d(phi_dot_0)/d(u_j) with weights; phi_dot_0 depends on u_0 - u_j
d_self = clean(sp.diff(S, du))          # d/du_0 = d/d(du)
d_other = clean(-sp.diff(S, du))        # d/du_j = -d/d(du)
print('A_PU(pi) pair term:  d_self + (-1)^j d_other  ->  (1 - (-1)^j) * d_self')
print('  d_self =', d_self)
print('  a + b (even j) =', sp.factor(sp.simplify(a_self + b_other)))
print('  a - b (odd j)  =', sp.factor(sp.simplify(a_self - b_other)))
