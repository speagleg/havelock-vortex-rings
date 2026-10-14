"""Horocycle limit of the hypercycle threshold series (tanh model).
As alpha -> 0 with arc spacing l = l'/sin(alpha) fixed, n alpha -> t and the Lambert series F(pi)/alpha... becomes
  I(l) = int_0^inf [ (2 sin 2t - 2t sin^2 t)/(e^{l t} - 1) - 2t cos^2 t/(e^{l t} + 1) ] dt,
whose zero should be the horocycle threshold l* = pi/eta_inf = 0.8439085 (note, Cor. 20).
Closed form via int_0^inf sin(at)/(e^{bt} -+ 1) dt and its a-derivative (sympy below).
Usage: python horo_limit.py"""
import numpy as np, sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq

def I(l):
    f = lambda t: (2*np.sin(2*t) - 2*t*np.sin(t)**2)/np.expm1(l*t) - 2*t*np.cos(t)**2/(np.exp(l*t) + 1)
    return quad(f, 0, 60/l, limit=400)[0]
print('numerical zero of I(l):', brentq(I, 0.5, 1.5, xtol=1e-13))

# closed form: use t sin^2 t = (t - t cos 2t)/2, t cos^2 t = (t + t cos 2t)/2, and
#   int sin(at)/(e^{bt}-1) = pi/(2b) coth(pi a/b) - 1/(2a),   int sin(at)/(e^{bt}+1) = 1/(2a) - pi/(2b) csch(pi a/b),
#   int t/(e^{bt}-1) = pi^2/(6b^2),   int t/(e^{bt}+1) = pi^2/(12 b^2),   int t cos(at)/(...) = d/da int sin(at)/(...)
a, b = sp.symbols('a b', positive=True)
Sm = sp.pi/(2*b)*sp.coth(sp.pi*a/b) - 1/(2*a)
Sp = 1/(2*a) - sp.pi/(2*b)*sp.csch(sp.pi*a/b)
Cm, Cp = sp.diff(Sm, a), sp.diff(Sp, a)          # int t cos(at)/(e^{bt} -+ 1)
l = sp.symbols('l', positive=True)
expr = (2*Sm - (sp.pi**2/(6*b**2) - Cm) - (sp.pi**2/(12*b**2) + Cp)).subs(a, 2).subs(b, l)
eta = sp.symbols('eta', positive=True)
expr_eta = sp.simplify((expr.subs(l, sp.pi/eta)).rewrite(sp.exp))
print('closed form I(pi/eta) =', sp.simplify(expr.subs(l, sp.pi/eta)))
f = sp.lambdify(eta, expr.subs(l, sp.pi/eta), 'mpmath')
import mpmath as mp; mp.mp.dps = 30
root = mp.findroot(f, 3.72)
print('eta root =', root, '-> l* = pi/eta =', mp.pi/root)
g = lambda e: (mp.coth(e) - 1/e)**2 - 2*mp.tanh(e)/e
print("note's eta_inf (root of (coth e - 1/e)^2 = 2 tanh(e)/e) =", mp.findroot(g, 3.72))
