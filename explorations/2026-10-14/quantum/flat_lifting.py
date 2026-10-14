"""Next order: the relaxed transverse stiffness Sigma(q) = s^2 (B|C| - delta^2)/|C| of the exact tanh symbol, expanded
in eps = e^{-l'} through eps^2. Pair functions are expanded in eps first; cos(jq) = T_j(X), sin(jq) = sin(q) U_{j-1}(X).
Usage: python flat_lifting.py -> flat_lifting.out"""
import sympy as sp
K, X, e = sp.symbols('K X epsilon', positive=True); s2 = 1 - K; k = sp.sqrt(K)
ORDER = 4
def ser(expr): return sp.series(expr, e, 0, ORDER).removeO()
def pair(j):
    c = (e**(-j) + e**j)/2; sh = (e**(-j) - e**j)/2; D = c + 1 - 2*K; den = (c - 1)*D**2
    a = s2*((3*K - 1)*c**2 - (4*K**2 - K + 1)*c + 2*K**2)/den
    b = s2*(K*c**2 + (1 - 3*K)*c + 2*K**2 - 2*K + 1)/den
    d = -s2**2*(c**2 + c - 2*K)/den
    f = -k*s2**sp.Rational(3, 2)*sh/D**2
    return [sp.expand(ser(x)) for x in (a, b, d, f)]
S2 = 1 - X**2                                   # sin^2 q
B = 0; C = 0; dl_over_sin = 0
for j in (1, 2, 3):
    a, b, d, f = pair(j)
    B += 2*(a + b*sp.chebyshevt(j, X)); C += 2*d*(1 - sp.chebyshevt(j, X)); dl_over_sin += 2*f*sp.chebyshevu(j - 1, X)
# Sigma = s^2 (B - delta^2/|C|) with |C| = -C ; expand the ratio in eps
ratio = sp.expand(ser(sp.expand(dl_over_sin**2*S2)/sp.expand(-C)))
Sig = sp.expand(ser(s2*(B - ratio)))
c1 = sp.factor(Sig.coeff(e, 1)); c2 = sp.factor(sp.simplify(Sig.coeff(e, 2)))
lines = [f'Sigma = eps * [{c1}] + eps^2 * [{c2}] + O(eps^3),   X = cos q']
lines.append(f'q-dependence at eps^2: d/dX = {sp.factor(sp.diff(c2, X))};  spread (X = 1 minus X = -1) = {sp.factor(c2.subs(X, 1) - c2.subs(X, -1))}')
lines.append(f'at K = 1/2: Sigma = eps^2 * [{sp.factor(c2.subs(K, sp.Rational(1, 2)))}] + O(eps^3)')
open('flat_lifting.out', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
