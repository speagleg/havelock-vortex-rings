"""Computer-assisted proof (Arb, python-flint) that for the screened kernel Q_1 (Delta = 2) every hypercycle row with
kappa = cos(alpha) <= 1/sqrt 2 is linearly unstable at EVERY spacing, so kappa_c = 1/sqrt 2 is exact for Delta = 2 too.
Criterion (../rows/modes_general.py): unstable iff F(pi) < 0, with the mode (Lambert) form, r = e^{-l'},
   F = sum_{n>=2} w_n [A_n r^n/(1-r^n) - B_n r^n/(1+r^n)],  w_n = 2/(n(n^2-1)),
   psi_n = n cos(n a) - cot(a) sin(n a),  A_n = 2 s k psi psi' + (2 - s^2 n^2) psi^2,  B_n = s^2 psi'^2,
and the pair (j) form F = sum_{j>=1} [d(jl') + e(jl') + (-1)^j b(jl')], a = d + e, from the closed-form kernel Q_1.
Domain a in [pi/4, pi/2] (s^2 in [1/2, 1], |cot| <= 1, csc^2 <= 2), so |psi_n| <= n+1, |psi_n'| <= n^2+n+2.
Regions:
 L  l' >= L1: F = c2 r^2 + c3 r^3 + E with c2 = (8/3) s^4 cos 2a <= 0, c3 = (16/3) s^4 (6 s^4 - 9 s^2 + 2) <= -4/3
    (6y^2 - 9y + 2 < 0 on [1/2, 1]) and |E| <= Eb r^4; F < 0 if r1 Eb(r1) < 4/3.
 M  l0 <= l' <= L1: adaptive Arb boxes on the mode form (with tail bound) show F < 0.
 S  l' <= l0: lemmas (each proved below or in the docstring of a check):
    b = s^2 sum w_n psi_n'^2 e^{-nx} > 0 decreasing  =>  sum_j (-1)^j b_j <= -b_1 + b_2 <= 0;
    d = -s^2 sum n^2 w_n psi_n^2 e^{-nx} < 0           =>  sum_j d_j <= d(l');
    e = a - d = 2 Q1(C) + 2 k^2 f(C), f(C) = (C-1)|Q1'(C)|, C = 1 + (cosh x - 1)/s^2: positive and decreasing in x
        (Q1' < 0 since atanh t < t/(1-t^2); f' < 0 since atanh t > t/(1+t) + t^2/(1+t)^2, t = 1/C), so
        sum_j e_j <= (1/l') int_0^inf e <= (1/l') [int ln(1 + 4/x^2) + k^2 int sech] = (2 pi + pi k^2/2)/l' <= (9 pi/4)/l',
        using Q1 <= Q0 = atanh(1/C) (atanh t <= t/(1-t)), f(C) <= 1/(2C) (polynomial check), C - 1 >= x^2/2, C >= cosh x;
    d(x) <= D_up(x) = s^2 [c (X+s^2)(X+2s^2) - 2 (X+2) s^4]/(X (X+2s^2)^2), X = cosh x - 1 (drop the atanh term).
    Hence F <= (s^2/l'^2) [P(l') + (9 pi/4) l'/s^2], P = x^2 D_up/s^2, and P(l') + (9 pi/4) l'/s^2 < 0 is checked in Arb.
Usage: python q1_instability_cert.py"""
import math, time
from flint import arb, ctx
import sympy as sp
ctx.prec = 96
PI = arb.pi()
def iv(a, b):
    a, b = float(a), float(b); return arb((a + b)/2, (b - a)/2 + 4e-16*max(abs(a), abs(b), 1e-300))
def hi(x): return float(x.mid()) + float(x.rad())
A_LO, A_HI = math.pi/4, math.pi/2
def Apoly(n): return 2*(n + 1)*(n*n + n + 2) + (n*n + 2)*(n + 1)**2
def Bpoly(n): return (n*n + n + 2)**2
def wt(n): return 2/(n*(n*n - 1))

# ---------- lemma checks (exact / symbolic) ----------
def lemmas():
    t, y = sp.symbols('t y', positive=True)
    ok = {}
    # f(C) <= 1/(2C):  f = 1/(1+t) - ((1-t)/t) atanh t <= 1/(1+t) - (1-t)(1 + t^2/3) (atanh t >= t + t^3/3); need <= t/2
    h = sp.together(t/2 - 1/(1 + t) + (1 - t)*(1 + t**2/3))
    num = sp.Poly(sp.numer(h)*1, t)         # denominator 6(1+t) > 0
    ok['f(C) <= 1/(2C)'] = all(v >= 0 for v in [sp.N(num.eval(sp.Rational(i, 1000))) for i in range(0, 1001)]) and \
        len([r for r in sp.real_roots(num) if 0 < r < 1]) == 0
    # f' < 0: g(t) = atanh t - t/(1+t) - t^2/(1+t)^2 has g(0) = 0, g'(t) = 4t^2/((1-t)(1+t)^3) > 0
    g = sp.atanh(t) - t/(1 + t) - t**2/(1 + t)**2
    ok["f' < 0 (g' = 4t^2/((1-t)(1+t)^3))"] = sp.simplify(sp.diff(g, t) - 4*t**2/((1 - t)*(1 + t)**3)) == 0
    # c3 sign: 6y^2 - 9y + 2 < 0 on [1/2, 1]
    ok['6y^2-9y+2 < 0 on [1/2,1]'] = all(r < sp.Rational(1, 2) or r > 1 for r in sp.real_roots(6*y**2 - 9*y + 2))
    # c3 <= -4/3: q(y) = y^2(6y^2 - 9y + 2) <= -1/4 on [1/2, 1] (the maximum is at y = 1/2)
    qpoly = sp.Poly(sp.expand(-sp.Rational(1, 4) - y**2*(6*y**2 - 9*y + 2)), y)
    inside = [r for r in sp.real_roots(qpoly) if sp.Rational(1, 2) < r < 1]
    ok['c3 <= -4/3 on the domain'] = len(inside) == 0 and qpoly.eval(sp.Rational(3, 4)) >= 0
    # c2, c3 closed forms against the mode definition
    ph = sp.symbols('phi', positive=True)
    def AB(n):
        psi = n*sp.cos(n*ph) - sp.cot(ph)*sp.sin(n*ph); dpsi = sp.diff(psi, ph); s_, k_ = sp.sin(ph), sp.cos(ph)
        return 2*s_*k_*psi*dpsi + (2 - s_**2*n**2)*psi**2, s_**2*dpsi**2
    A2, B2 = AB(2); A3, B3 = AB(3)
    # identities of trigonometric polynomials (degree <= 16 in e^{i phi}): checked at 60 points with 50 digits
    import mpmath as mpm
    mpm.mp.dps = 50
    e2 = sp.lambdify(ph, sp.Rational(1, 3)*(A2 - B2) - sp.Rational(8, 3)*sp.sin(ph)**4*sp.cos(2*ph), 'mpmath')
    e3 = sp.lambdify(ph, sp.Rational(1, 12)*(A3 - B3) - sp.Rational(16, 3)*sp.sin(ph)**4*(6*sp.sin(ph)**4 - 9*sp.sin(ph)**2 + 2), 'mpmath')
    pts = [mpm.mpf(i)/61*mpm.pi for i in range(1, 61)]
    ok['c2 formula'] = max(abs(e2(p_)) for p_ in pts) < mpm.mpf(10)**-40
    ok['c3 formula'] = max(abs(e3(p_)) for p_ in pts) < mpm.mpf(10)**-40
    return ok

# ---------- region L ----------
def region_L(L1):
    """r1 * Eb(r1), Eb bounding |E|/r^4 for every r <= r1 and every alpha in the domain (Eb increases with r)."""
    r = math.exp(-L1)
    Eb = wt(2)*(Apoly(2)/(1 - r*r) + Bpoly(2)) + wt(3)*(Apoly(3)/(1 - r**3) + Bpoly(3))*r*r
    for n in range(4, 61):
        Eb += wt(n)*((Apoly(n) + Bpoly(n))*r**(n - 4) + (Apoly(n)/(1 - r**n) + Bpoly(n))*r**(2*n - 4))
    # n >= 61: g(n) = w_n (Apoly + Bpoly) <= 5n, and r^(2n-4)/(1-r^n) <= r^(n-4)/(1-r):
    qq = r*(1 + 1/61); Eb += 5*(1 + 1/(1 - r))*61*r**(61 - 4)/(1 - qq)
    return r*Eb

# ---------- region M ----------
def F_box(al, lp, N):
    s, k = al.sin(), al.cos(); cot, csc2 = k/s, 1/(s*s)
    r = (-lp).exp(); rmax = math.exp(-(float(lp.mid()) - float(lp.rad())))
    F = arb(0)
    for n in range(2, N + 1):
        sn, cn = (n*al).sin(), (n*al).cos()
        psi = n*cn - cot*sn; dpsi = -n*n*sn + csc2*sn - n*cot*cn
        A = 2*s*k*psi*dpsi + (2 - s*s*n*n)*psi*psi; B = s*s*dpsi*dpsi
        rn = r**n
        F += wt(n)*(A*rn/(1 - rn) - B*rn/(1 + rn))
    # tail: sum_{n>N} w_n (Apoly + Bpoly) r^n/(1-r), ratio bound
    # tail: g(n) = w_n (Apoly + Bpoly)(n) <= 5 n for n >= 30 (exact polynomial check in tail_constant), so
    # sum_{n>N} g(n) r^n/(1-r^n) <= 5/(1-rmax) sum_{n>=N+1} n rmax^n  (ratio bound)
    m = N + 1; qq = rmax*(1 + 1/m); assert m >= 30 and qq < 1
    T = 5*m*rmax**m/(1 - qq)/(1 - rmax)
    return F + arb(0, T)
def region_M(l0, L1):
    stack = [(A_LO, A_HI, l0, L1, 0)]; nb = 0
    while stack:
        a0, a1, x0, x1, dpt = stack.pop(); nb += 1
        assert dpt < 40, ('depth', a0, a1, x0, x1)
        N = max(int(45/x0) + 6, 30)
        F = F_box(iv(a0, a1), iv(x0, x1), N)
        if F < 0: continue
        if (a1 - a0)*6 > (x1 - x0): m = (a0 + a1)/2; stack += [(a0, m, x0, x1, dpt+1), (m, a1, x0, x1, dpt+1)]
        else: m = (x0 + x1)/2; stack += [(a0, a1, x0, m, dpt+1), (a0, a1, m, x1, dpt+1)]
    return nb

def tail_constant():
    """exact: 5 n * n (n^2 - 1)/2 - (Apoly + Bpoly)(n) >= 0 for every integer n >= 30."""
    n = sp.symbols('n', positive=True)
    poly = sp.expand(sp.Rational(5, 2)*n*n*(n*n - 1) - (2*(n + 1)*(n*n + n + 2) + (n*n + 2)*(n + 1)**2 + (n*n + n + 2)**2))
    roots = [r for r in sp.real_roots(sp.Poly(poly, n)) if r >= 30]
    return poly, len(roots) == 0 and poly.subs(n, 30) > 0

# ---------- region S ----------
def region_S(l0):
    """P(x) + (9 pi/4) x/s^2 < 0 for x in (0, l0], a in [pi/4, pi/2]; x^2/X enclosed in [2/cosh x, 2]."""
    worst = -1e9; stack = [(A_LO, A_HI, 0.0, l0)]; nb = 0
    while stack:
        a0, a1, x0, x1 = stack.pop(); nb += 1
        al, x = iv(a0, a1), iv(x0, x1)
        s2 = al.sin()**2; c = x.cosh(); X = c - 1
        q_lo = 2/arb(x1).cosh(); ratio = arb((float(q_lo.mid()) + 2)/2, (2 - float(q_lo.mid()))/2 + 1e-12)   # x^2/X in [2/cosh x1, 2]
        Xs = X if x0 > 0 else arb((hi(X))/2, hi(X)/2 + 1e-30)                                              # X in [0, cosh x1 - 1]
        P = ratio*(c*(Xs + s2)*(Xs + 2*s2) - 2*(Xs + 2)*s2*s2)/((Xs + 2*s2)**2)
        val = P + (9*PI/4)*x/s2
        if val < 0: worst = max(worst, hi(val)); continue
        if (a1 - a0) > (x1 - x0)*4: m = (a0 + a1)/2; stack += [(a0, m, x0, x1), (m, a1, x0, x1)]
        else: m = (x0 + x1)/2; stack += [(a0, a1, x0, m), (a0, a1, m, x1)]
        assert nb < 200000
    return nb, worst

if __name__ == '__main__':
    t0 = time.time()
    print('lemmas:', lemmas())
    print('tail constant: g(n) <= 5n for n >= 30:', tail_constant())
    L1, l0 = 4.5, 0.05
    rEb = region_L(L1)
    print(f'L (l\' >= {L1}): r1*Eb = {rEb:.4f} < 4/3: {rEb < 4/3}')
    nbS, w = region_S(l0)
    print(f'S (l\' <= {l0}): {nbS} boxes, max of P + (9pi/4) l\'/s^2 <= {w:.4f} < 0')
    nbM = region_M(l0, L1)
    print(f'M ({l0} <= l\' <= {L1}): {nbM} boxes, F(pi) < 0 on all')
    print(f'total {time.time() - t0:.0f} s')
