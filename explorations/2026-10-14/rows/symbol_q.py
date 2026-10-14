"""Pair terms of the full Bloch symbol of a hypercycle row (tanh model), K = kappa^2, c = cosh(u_0 - u_j).
uu_self, pp_self, pp_other: diagonal-entry pair terms (odd in u_0 - u_j, so the diagonal is i*delta(q), equal in both slots);
a, b: d(u_dot_0)/d(phi_0), d(u_dot_0)/d(phi_j);  d: d(phi_dot_0)/d(u_0).  Used by exact_symbol.py.
Usage: python symbol_q.py"""
import sympy as sp
p0, pj, du, al = sp.symbols('p0 pj du alpha', real=True)
c, K = sp.symbols('c K', positive=True)
C = 1 + (sp.cosh(du) - sp.cos(p0-pj))/(sp.sin(p0)*sp.sin(pj))
Wp = -1/(C**2-1)
T = sp.sin(p0)**2*Wp*sp.diff(C, p0)     # u_dot_0 pair term
S = -sp.sin(p0)**2*Wp*sp.diff(C, du)    # phi_dot_0 pair term
def clean(e):
    e = sp.simplify(e.subs({p0: al, pj: al}))
    e = e.subs(sp.sinh(du)**2, sp.cosh(du)**2-1)
    e = sp.simplify(e).subs(sp.cosh(du), c)
    e = e.subs(sp.cos(2*al), 2*sp.cos(al)**2-1)
    e = e.subs(sp.sin(al), sp.sqrt(1-K)).subs(sp.cos(al), sp.sqrt(K))
    return sp.factor(sp.simplify(e))
for name, e in [('uu_self', sp.diff(T, du)), ('pp_self', sp.diff(S, p0)), ('pp_other', sp.diff(S, pj)),
                ('a', sp.diff(T, p0)), ('b', sp.diff(T, pj)), ('d', sp.diff(S, du))]:
    r = clean(e); print(name, '=', r)
