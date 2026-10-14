"""Symbolic check of the large-spacing criterion for hypercycle rows (sympy).
At large spacing only nearest neighbours matter and W ~ 1/cosh d ~ 2 sin(phi_i) sin(phi_j) exp(-|u_i - u_j|).
Linearising the chain about u_k = k l', phi_k = alpha and Fourier-transforming (wavenumber q) gives a 2x2 symbol
with equal diagonal entries, so the row is stable iff the product of the off-diagonal entries is <= 0.
Result: off-diagonal (U <- P) is proportional to  cos^2(alpha) cos(q) + 2 cos^2(alpha) - sin^2(alpha),
off-diagonal (P <- U) is proportional to  4 cos(q) - 4 <= 0.
Worst mode q = pi (pairing): cos(2 alpha) >= 0  <=>  kappa = cos(alpha) >= 1/sqrt(2).
Usage: python asymptotic_check.py"""
import sympy as sp

K = 5
u = sp.symbols('u0:%d' % K); ph = sp.symbols('p0:%d' % K)
H = sum(2*sp.sin(ph[k])*sp.sin(ph[k+1])*sp.exp(-(u[k+1]-u[k])) for k in range(K-1))
k = 2
ud = sp.sin(ph[k])**2*sp.diff(H, ph[k]); pd = -sp.sin(ph[k])**2*sp.diff(H, u[k])
al, lp, q = sp.symbols('alpha lp q', positive=True)
sub = {ph[j]: al for j in range(K)}; sub.update({u[j]: j*lp for j in range(K)})

def symbol(expr):
    cU = sum(sp.diff(expr, u[j]).subs(sub)*sp.exp(sp.I*q*(j-k)) for j in range(K))
    cP = sum(sp.diff(expr, ph[j]).subs(sub)*sp.exp(sp.I*q*(j-k)) for j in range(K))
    return sp.simplify(cU), sp.simplify(cP)

(AUU, AUP), (APU, APP) = symbol(ud), symbol(pd)
print('diagonal entries equal:', sp.simplify(AUU-APP) == 0)
crit = sp.simplify(AUP/(4*sp.sin(al)**2*sp.exp(-lp)))
print('normal-stiffness factor:', sp.simplify(sp.expand(crit.rewrite(sp.cos))))
print('at q = pi:', sp.simplify(crit.subs(q, sp.pi)), '-> zero at alpha =', sp.solve(sp.Eq(sp.simplify(crit.subs(q, sp.pi)), 0), al))
print('other factor:', sp.simplify(APU/(sp.sin(al)**4*sp.exp(-lp))))
