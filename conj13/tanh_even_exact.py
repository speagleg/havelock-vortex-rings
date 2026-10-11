"""Conj13, tanh (image / 'Kimura') model on H^2, even n = 2M: exact reduced function on the
D_{M}-fixed family (two staggered M-gons, radii^2 u_+, u_-) restricted to the momentum level set,
then threshold (h2 = 0) and the decisive quartic h4, with an Arb-certified sign.

Model (unit disk, curvature -1), 2x the item's H_tanh:
  H = -sum_{i<j} log|z_i-z_j|^2 + sum_{i<j} log|1 - z_i conj z_j|^2,   J3 = sum u/(1-u)  (u=|z|^2).
Product formulas (eta = M-th roots of 1, zeta = M-th roots of -1):
  prod_{k=1}^{M-1} (1 - u eta^k) = (1-u^M)/(1-u),   prod_zeta (1 - a zeta) = 1 + a^M,
  prod_zeta (x - y zeta) = x^M + y^M.
Hence (up to additive constants)
  sum log|z_i-z_j|^2   = M(M-1)/2 (log u+ + log u-) + 2M log(u+^{M/2} + u-^{M/2})
  sum log|1-z_i z_j*|^2 = M log[(1-u+^M)/(1-u+)] + M log[(1-u-^M)/(1-u-)] + 2M log(1 + (u+ u-)^{M/2}).
Level set: g(u_pm) = g0 +- tau, g(u) = u/(1-u), u = g/(1+g).
"""
import sympy as sp, mpmath as mp
from flint import arb, fmpz_poly, ctx

tau, u0 = sp.symbols('tau u0', positive=True)

def reduced(Mv):
    g0 = u0/(1-u0)
    gp, gm = g0 + tau, g0 - tau
    up, um = gp/(1+gp), gm/(1+gm)
    h = -(sp.Rational(Mv*(Mv-1), 2)*(sp.log(up)+sp.log(um)) + 2*Mv*sp.log(up**sp.Rational(Mv, 2) + um**sp.Rational(Mv, 2))) \
        + Mv*sp.log((1-up**Mv)/(1-up)) + Mv*sp.log((1-um**Mv)/(1-um)) + 2*Mv*sp.log(1 + (up*um)**sp.Rational(Mv, 2))
    s = sp.series(h, tau, 0, 5).removeO()
    c = [sp.factor(sp.simplify(s.coeff(tau, k))) for k in range(5)]
    return c

# numeric sanity check of the product formulas against a direct sum
mp.mp.dps = 30
def direct(Mv, up_, um_):
    n = 2*Mv; w = mp.exp(2j*mp.pi/n)
    z = [(mp.sqrt(up_) if k % 2 == 0 else mp.sqrt(um_))*w**k for k in range(n)]
    return sum(-mp.log(abs(z[i]-z[j])**2) + mp.log(abs(1-z[i]*mp.conj(z[j]))**2) for i in range(n) for j in range(i+1, n))
def closed(Mv, up_, um_):
    return -(Mv*(Mv-1)/mp.mpf(2)*(mp.log(up_)+mp.log(um_)) + 2*Mv*mp.log(up_**(mp.mpf(Mv)/2)+um_**(mp.mpf(Mv)/2))) \
        + Mv*mp.log((1-up_**Mv)/(1-up_)) + Mv*mp.log((1-um_**Mv)/(1-um_)) + 2*Mv*mp.log(1+(up_*um_)**(mp.mpf(Mv)/2))
for Mv in (2, 3, 4, 5):
    d = (direct(Mv, mp.mpf('0.3'), mp.mpf('0.1')) - direct(Mv, mp.mpf('0.2'), mp.mpf('0.25'))) \
        - (closed(Mv, mp.mpf('0.3'), mp.mpf('0.1')) - closed(Mv, mp.mpf('0.2'), mp.mpf('0.25')))
    print(f"product formula check M={Mv}: {mp.nstr(d, 5)}")

for Mv in (4, 5, 6):
    n = 2*Mv
    c = reduced(Mv)
    assert c[1] == 0 and c[3] == 0
    h2num = sp.factor(sp.numer(sp.together(c[2])))
    print(f"\nn={n}: h2 = {c[2]}")
    print(f"n={n}: h4 = {c[4]}")
    # threshold: root of h2 numerator in (0,1)
    P = sp.Poly(sp.numer(sp.together(c[2])), u0)
    facs = sp.factor_list(P.as_expr())[1]
    roots = []
    for f, e in facs:
        fp = sp.Poly(f, u0)
        if fp.degree() < 1:
            continue
        for iv in fp.intervals(eps=sp.Rational(1, 10**20)):
            (a, bnd), mult = iv
            if 0 < a and bnd < 1:
                roots.append((fp, a, bnd))
    print(f"n={n}: roots of h2 in (0,1): {[(r[0].degree(), float(r[1])) for r in roots]}")
    for fp, a, bnd in roots:
        # Arb: refine the root by bisection in balls, then evaluate h4 there
        ctx.prec = 400
        coeffs = [int(cc) for cc in reversed(fp.all_coeffs())]
        fz = fmpz_poly(coeffs)
        print(f"   factor (deg {fp.degree()}): {fp.as_expr()}")
        lo, hi = sp.Rational(a), sp.Rational(bnd)
        for _ in range(400):
            mid = (lo+hi)/2
            if sp.sign(fp.eval(mid)) == sp.sign(fp.eval(lo)):
                lo = mid
            else:
                hi = mid
        # exact rational enclosure [lo,hi] of xi*, width 2^-400 * initial
        X = arb(sp.Rational(lo+hi, 2).p) / arb(sp.Rational(lo+hi, 2).q)
        X = X.union(arb(lo.p)/arb(lo.q)).union(arb(hi.p)/arb(hi.q))
        h4f = sp.lambdify(u0, c[4], 'mpmath')
        mp.mp.dps = 120
        print(f"   xi* = {mp.nstr(mp.mpf(lo.p)/lo.q, 50)}  (exact isolating interval width {float(hi-lo):.2e})")
        # certified evaluation of h4 on the ball X: rational function in u0 -> evaluate numer/denom with arb
        num, den = sp.fraction(sp.together(c[4]))
        pn = sp.Poly(sp.expand(num), u0); pd = sp.Poly(sp.expand(den), u0)
        def ev(poly):
            acc = arb(0)
            for (k,), cf in poly.terms():
                acc += arb(int(cf)) * X**k
            return acc
        val = ev(pn)/ev(pd)
        print(f"   h4(xi*) in Arb ball: {val.str(30, radius=True)}   positive: {val > 0}")
        print(f"   h4(xi*) mpmath 50 digits: {mp.nstr(h4f(mp.mpf(lo.p)/lo.q), 30)}")
