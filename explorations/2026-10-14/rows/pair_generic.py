"""sympy: hypercycle-row pair functions for a generic kernel W = w(C), C = cosh d, in terms of W1 = w'(C), W2 = w''(C),
X = cosh(u_0 - u_j) - 1, at the row phi_0 = phi_j = alpha. Used in kernel_coeffs.py.
Usage: python pair_generic.py"""
import sympy as sp
p0, pj, u, al = sp.symbols('p0 pj u alpha', real=True)
W1, W2 = sp.symbols('W1 W2')   # w'(C), w''(C)
C = 1 + (sp.cosh(u) - sp.cos(p0-pj))/(sp.sin(p0)*sp.sin(pj))
def D1(f, v): return sp.diff(f, v)
Cp0, Cpj, Cu = D1(C, p0), D1(C, pj), D1(C, u)
# T = sin^2 p0 w'(C) C_p0 ;  S = -sin^2 p0 w'(C) C_u
b = sp.sin(p0)**2*(W2*Cpj*Cp0 + W1*D1(Cp0, pj))
a = D1(sp.sin(p0)**2, p0)*W1*Cp0 + sp.sin(p0)**2*(W2*Cp0**2 + W1*D1(Cp0, p0))
d = -sp.sin(p0)**2*(W2*Cu**2 + W1*D1(Cu, u))
X = sp.symbols('X', positive=True)
for name, e in [('a', a), ('b', b), ('d', d)]:
    e = sp.simplify(e.subs({p0: al, pj: al}))
    e = sp.simplify(e.subs(sp.sinh(u)**2, sp.cosh(u)**2 - 1).subs(sp.cosh(u), X + 1))
    print(name, '=', sp.factor(e))
