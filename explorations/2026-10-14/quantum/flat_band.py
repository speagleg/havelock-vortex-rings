"""Leading order (nearest neighbours, large spacing) structure of the quadratic Hamiltonian of a hypercycle row, and the
flat band at kappa = 1/sqrt 2.
Nearest-neighbour model (exact leading term of -ln tanh(d/2) ~ 2 e^{-d}, e^{-d} ~ sin(phi_i) sin(phi_j) e^{-|du|}/... ):
    H = sum_j w_j,  w_j = 2 sin(phi_j) sin(phi_{j+1}) e^{-(u_{j+1}-u_j)},   I = sum_j Y_j,  Y = -cot(phi)  ({u, Y} = 1).
A row u_j = j l', phi_j = alpha is a critical point of H - V I for the right V. We compute the Hessian of H - V I in the
canonical coordinates (u, Y), Bloch-transform it, and check:
 (i)  the symbol agrees with the leading terms of the exact pair functions (a, b, d, f of symbol.py);
 (ii) det h(q) = 16 eps^2 s^6 (2K - 1)(1 - cos q), eps = e^{-l'}/... (normalisation below);
 (iii) at K = 1/2 the Hessian equals (w/2) sum_j (delta ln w_j)^2: bending at fixed bond weights costs nothing;
 (iv) the chiral dispersion E(q) = sqrt(-BC) + delta vanishes identically on 0 < q < pi at K = 1/2.
Usage: python flat_band.py -> flat_band.out"""
import sympy as sp
K = sp.symbols('K', positive=True); s = sp.sqrt(1 - K); k = sp.sqrt(K)
q, lp, V = sp.symbols('q lp V', real=True)
# three consecutive vortices j-1, j, j+1 suffice for the Hessian rows of vortex j (nearest neighbours)
u = sp.symbols('u0:3', real=True); Y = sp.symbols('Y0:3', real=True)
phi = [sp.acot(-y) if False else sp.atan2(1, -y) for y in Y]     # Y = -cot(phi), phi in (0, pi)
w = lambda i: 2*sp.sin(phi[i])*sp.sin(phi[i + 1])*sp.exp(-(u[i + 1] - u[i]))
Hloc = w(0) + w(1) - V*Y[1]                                     # terms of H - V I involving vortex 1
Y0 = -sp.cot(sp.acos(k))                                          # = -k/s
pt = {u[0]: -lp, u[1]: 0, u[2]: lp, Y[0]: Y0, Y[1]: Y0, Y[2]: Y0}
Vsol = sp.solve(sp.Eq(sp.diff(Hloc, Y[1]).subs(pt), 0), V)[0]    # equilibrium: dH/dY_1 = V
def d2(x, y): return sp.simplify(sp.diff(Hloc, x, y).subs(pt).subs(V, Vsol))
# Bloch symbol of the Hessian: h(q) = sum_m Hess(1, 1+m) e^{i q m}
h_uu = d2(u[1], u[1]) + 2*d2(u[1], u[2])*sp.cos(q)
h_YY = d2(Y[1], Y[1]) + 2*d2(Y[1], Y[2])*sp.cos(q)
h_uY = d2(u[1], Y[1]) + d2(u[1], Y[2])*sp.exp(sp.I*q) + d2(u[1], Y[0])*sp.exp(-sp.I*q)
h_uY = sp.simplify(h_uY)
eps = sp.exp(-lp)
print('V =', sp.simplify(Vsol))
print('h_uu(q) =', sp.factor(sp.simplify(h_uu/eps)), '* e^{-l\'}')
print('h_YY(q) =', sp.factor(sp.simplify(h_YY/eps)), '* e^{-l\'}')
print('h_uY(q) =', sp.simplify(h_uY/eps), '* e^{-l\'}')
Kr = sp.symbols('Kr', positive=True)
det = sp.expand((h_uu*h_YY - h_uY*sp.conjugate(h_uY)).subs(K, Kr))
det = sp.simplify(sp.expand_complex(det)/eps**2)
print('(ii) det h(q) / e^{-2l\'} - 16 s^6 (2K-1)(1-cos q) =', sp.simplify(sp.trigsimp(det - 16*(1 - Kr)**3*(2*Kr - 1)*(1 - sp.cos(q))).subs(sp.Abs(Kr - 1), 1 - Kr)))
# (i) compare with the exact symbol's leading terms: -C/s^2 = h_uu, B s^2 = h_YY, -i delta = h_uY  (symbol.py conventions,
# with B, C, delta for the 2 pi-normalised tanh kernel whose leading pair term is 2 e^{-d}: the j = 1 pair functions at large c)
c = sp.symbols('c', positive=True); D = c + 1 - 2*K
a_ = s**2*((3*K - 1)*c**2 - (4*K**2 - K + 1)*c + 2*K**2)/((c - 1)*D**2)
b_ = s**2*(K*c**2 + (1 - 3*K)*c + 2*K**2 - 2*K + 1)/((c - 1)*D**2)
d_ = -s**4*(c**2 + c - 2*K)/((c - 1)*D**2)
f_ = -k*s**3*sp.sqrt(c**2 - 1)/D**2
lead = lambda e: 2*sp.limit(e*c, c, sp.oo)                        # coefficient of e^{-l'} (1/c ~ 2 e^{-l'})
Bq = 2*(lead(a_) + lead(b_)*sp.cos(q)); Cq = 2*lead(d_)*(1 - sp.cos(q)); dq = 2*lead(f_)*sp.sin(q)
chk = [sp.simplify(-Cq/s**2 - h_uu/eps), sp.simplify(Bq*s**2 - h_YY/eps), sp.simplify((-sp.I*dq - h_uY/eps).rewrite(sp.sin))]
print('(i) exact-symbol leading terms minus nearest-neighbour Hessian (uu, YY, uY):', [sp.simplify(sp.expand_complex(x)) for x in chk])
# (iii) K = 1/2: Hessian of H - V I equals (w/2) * sum (delta ln w_j)^2 on the three-vortex patch (rows of vortex 1)
dv = sp.symbols('du0:3 dY0:3', real=True); du, dY = dv[:3], dv[3:]
lnw = lambda i: sp.log(w(i))
def lin(expr):
    return sum(sp.diff(expr, x).subs(pt)*dx for x, dx in zip(list(u) + list(Y), list(du) + list(dY)))
quad_full = sum(sp.diff(Hloc, x, y).subs(pt).subs(V, Vsol)*dx*dy for x, dx in zip(list(u) + list(Y), list(du) + list(dY))
                for y, dy in zip(list(u) + list(Y), list(du) + list(dY)))
w0 = w(0).subs(pt)
target = w0*(lin(lnw(0))**2 + lin(lnw(1))**2)
diff_ = sp.simplify((quad_full - target).subs(K, sp.Rational(1, 2)))
# only the terms that belong to vortex 1's rows matter; compare the coefficients that involve du1 or dY1
coeffs = [sp.simplify(sp.diff(diff_, a, b)) for a in (du[1], dY[1]) for b in list(du) + list(dY)]
print('(iii) K = 1/2: Hessian - w * sum (d ln w)^2, coefficients in vortex 1 rows:', set(coeffs))
# (iv) chiral dispersion at K = 1/2
E = sp.sqrt(sp.simplify(h_uu*h_YY)) - sp.I*h_uY      # sqrt(-BC) + delta with -i delta = h_uY -> delta = i h_uY
E12 = sp.simplify((sp.sqrt(sp.factor(sp.simplify(h_uu*h_YY))) + sp.I*h_uY).subs(K, sp.Rational(1, 2))/eps)
import numpy as np
f = sp.lambdify(q, E12, 'numpy'); qs = np.linspace(-3.1, 3.1, 13)
print('(iv) K = 1/2: E(q)/e^{-l\'} at q =', np.round(qs, 2), '->', np.round(np.real(f(qs)), 12))
