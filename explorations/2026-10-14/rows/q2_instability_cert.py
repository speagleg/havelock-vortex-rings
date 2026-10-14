"""Computer-assisted proof (Arb, python-flint): for the screened kernel Q_2 (Delta = 3) every hypercycle row with
kappa = cos(alpha) <= 1/sqrt 2 is linearly unstable at every spacing (F_2(pi) < 0 on alpha in [pi/4, pi/2], all l').
Modes (Darboux ladder, nu = 2): psi_n = (1 - n^2 + 3q^2) sin(n a) - 3 n q cos(n a), q = cot(a), n >= 3,
  psi_n' = (-n^3 + 6 n q^2 + 4n) cos(n a) + (3 n^2 q - 6 q^3 - 6 q) sin(n a),  N_n = pi (n^2-1)(n^2-4)/2,
  w_n = 2/(n (n^2-1)(n^2-4)),  A_n = 2 s k psi psi' + (6 - s^2 n^2) psi^2,  B_n = s^2 psi'^2,
  F = sum_n w_n [A_n r^n/(1-r^n) - B_n r^n/(1+r^n)],  r = e^{-l'}.   On the domain |q| <= 1, so
  |psi_n| <= (n+1)(n+2), |psi_n'| <= n^3 + 3n^2 + 10n + 12.
Regions:
 L  l' >= L1: F = c3 r^3 + c4 r^4 + E, c3 = (16/5) s^6 cos 2a <= 0, c4 <= c4max (Arb over alpha; c4max = -0.760, true max -4/5), |E| <= Eb r^5;
    need r1 Eb < -c4max.
 M  l0 <= l' <= L1: Arb boxes on the mode form with an exact polynomial tail constant.
 S  l' <= l0: F <= sum d_j + sum e_j + sum (-1)^j b_j with b > 0 decreasing and d < 0 (mode expansion), so
    F <= d(l') + sum_j e(j l');  e = a - d = 6 Q2(C) - 2 k^2 (C-1) Q2'(C). From Laplace's integral
    Q_nu(C) = int_0^inf (C + sqrt(C^2-1) cosh t)^{-nu-1} dt and C + sqrt(C^2-1) cosh t >= C: Q2 <= Q0/C^2. And
    -Q2' <= 3/(C^2 (C^2-1)) (in fact <= 1/(C^2(C^2-1))): with t = 1/C, C^2(C^2-1)(-Q2')/3 = sum_{k>=2} 2 t^(2k-4)/(4k^2-1),
    whose coefficients are positive, so its sup is the value at t = 1, sum_{k>=2} 2/(4k^2-1) = 1/3 (proof: q2_lemmas.py), so e <= 3 ln(1 + 4/x^2) sech^2 x + 3 sech^2(x)/(cosh x + 1) =: e_up(x),
    which is positive and decreasing, giving sum_j e_j <= (1/l') int_0^inf e_up = K/l' (K by a rigorous Riemann sum).
    Dropping the (negative) Q0 terms of d gives d <= D_up(x) (rational in c = cosh x, X = c - 1, s^2), and
    F <= (s^2/l'^2)[P(l') + K l'/s^2],  P = x^2 D_up/s^2 -> -1, checked < 0 in Arb.
Usage: python q2_instability_cert.py"""
import math, time
from flint import arb, ctx
import mpmath as mp
ctx.prec = 96
PI = arb.pi()
def iv(a, b):
    a, b = float(a), float(b); return arb((a + b)/2, (b - a)/2 + 4e-16*max(abs(a), abs(b), 1e-300))
def hi(x): return float(x.mid()) + float(x.rad())
A_LO, A_HI = math.pi/4, math.pi/2
def wt(n): return 2/(n*(n*n - 1)*(n*n - 4))
def psib(n): return (n + 1)*(n + 2)
def dpsib(n): return n**3 + 3*n*n + 10*n + 12
def Apoly(n): return psib(n)*dpsib(n) + (n*n + 6)*psib(n)**2
def Bpoly(n): return dpsib(n)**2
def g(n): return wt(n)*(Apoly(n) + Bpoly(n))

def AB_arb(n, al):
    s, k = al.sin(), al.cos(); q = k/s
    S, Cn = (n*al).sin(), (n*al).cos()
    psi = (1 - n*n + 3*q*q)*S - 3*n*q*Cn
    dpsi = (-n**3 + 6*n*q*q + 4*n)*Cn + (3*n*n*q - 6*q**3 - 6*q)*S
    return 2*s*k*psi*dpsi + (6 - s*s*n*n)*psi*psi, s*s*dpsi*dpsi

def checks():
    mp.mp.dps = 50; ok = {}
    # mode equation residual and c3 formula at 40 points (trigonometric-rational identities)
    def psi_mp(n, a):
        q = mp.cot(a); return (1 - n*n + 3*q*q)*mp.sin(n*a) - 3*n*q*mp.cos(n*a)
    res = max(abs(-mp.diff(lambda t: psi_mp(5, t), a, 2) + 6/mp.sin(a)**2*psi_mp(5, a) - 25*psi_mp(5, a)) for a in [mp.mpf(i)/41*mp.pi for i in range(1, 41)])
    ok['mode equation (n=5)'] = res < mp.mpf(10)**-30
    def c3_mp(a):
        s, k = mp.sin(a), mp.cos(a); q = k/s; n = 3
        psi = (1 - n*n + 3*q*q)*mp.sin(n*a) - 3*n*q*mp.cos(n*a)
        dpsi = (-n**3 + 6*n*q*q + 4*n)*mp.cos(n*a) + (3*n*n*q - 6*q**3 - 6*q)*mp.sin(n*a)
        return mp.mpf(2)/(3*8*5)*(2*s*k*psi*dpsi + (6 - s*s*9)*psi**2 - s*s*dpsi**2) - mp.mpf(16)/5*s**6*mp.cos(2*a)
    ok['c3 = (16/5) s^6 cos 2a'] = max(abs(c3_mp(mp.mpf(i)/41*mp.pi)) for i in range(1, 41)) < mp.mpf(10)**-35
    # norm N_4 by quadrature (Darboux: N_n = pi (n^2-1)(n^2-4)/2)
    ok['norm N_4'] = abs(mp.quad(lambda t: psi_mp(4, t)**2, [0, mp.pi/2, mp.pi]) - mp.pi*15*12/2) < 1e-20
    # Arb upper bound c4max of c4 over the domain (true max -4/5 at a = pi/4)
    worst = -1e9
    for i in range(400):
        al = iv(A_LO + (A_HI - A_LO)*i/400, A_LO + (A_HI - A_LO)*(i + 1)/400)
        A4, B4 = AB_arb(4, al); worst = max(worst, hi(wt(4)*(A4 - B4)))
    ok['c4max'] = worst                     # Arb upper bound of c4 over the domain (true max -4/5 at a = pi/4)
    # tail constant g(n) <= 7n for n >= 30 (polynomial in n, checked exactly with sympy)
    import sympy as sp
    nn = sp.symbols('n', positive=True)
    poly = sp.expand(7*nn*nn*(nn*nn - 1)*(nn*nn - 4)/2 - ((nn + 1)*(nn + 2)*(nn**3 + 3*nn**2 + 10*nn + 12) + (nn**2 + 6)*((nn + 1)*(nn + 2))**2 + (nn**3 + 3*nn**2 + 10*nn + 12)**2))
    ok['g(n) <= 7n for n >= 30'] = poly.subs(nn, 30) > 0 and all(r < 30 for r in sp.real_roots(sp.Poly(poly, nn)))
    return ok

def region_L(L1):
    r = math.exp(-L1)
    Eb = g(3)*r/(1 - r**3) + g(4)*r**3/(1 - r**4)
    for n in range(5, 61):
        Eb += g(n)*(r**(n - 5) + r**(2*n - 5)/(1 - r**n))
    qq = r*(1 + 1/61); Eb += 7*(1 + 1/(1 - r))*61*r**(61 - 5)/(1 - qq)
    return r*Eb

def F_box(al, lp, N):
    r = (-lp).exp(); rmax = math.exp(-(float(lp.mid()) - float(lp.rad())))
    F = arb(0)
    for n in range(3, N + 1):
        A, B = AB_arb(n, al); rn = r**n
        F += wt(n)*(A*rn/(1 - rn) - B*rn/(1 + rn))
    m = N + 1; qq = rmax*(1 + 1/m); assert m >= 30 and qq < 1
    return F + arb(0, 7*m*rmax**m/(1 - qq)/(1 - rmax))
def region_M(l0, L1):
    stack = [(A_LO, A_HI, l0, L1, 0)]; nb = 0
    while stack:
        a0, a1, x0, x1, dpt = stack.pop(); nb += 1
        assert dpt < 45, ('depth', a0, a1, x0, x1)
        F = F_box(iv(a0, a1), iv(x0, x1), max(int(45/x0) + 6, 30))
        if F < 0: continue
        if (a1 - a0)*6 > (x1 - x0): m = (a0 + a1)/2; stack += [(a0, m, x0, x1, dpt+1), (m, a1, x0, x1, dpt+1)]
        else: m = (x0 + x1)/2; stack += [(a0, a1, x0, m, dpt+1), (a0, a1, m, x1, dpt+1)]
    return nb

def K_integral(h=1e-3, xmax=30.0):
    """rigorous upper bound for int_0^inf e_up, e_up(x) = 3 ln(1+4/x^2) sech^2 x + 3 sech^2 x/(cosh x + 1) (decreasing)."""
    H = arb(h)
    first = 3*(H*(1 + 4/H**2).log() + 4*(H/2).atan()) + 3*H/2          # sech^2 <= 1, 1/(cosh+1) <= 1/2 on [0, h]
    tot = first; x = h
    while x < xmax:
        step = 1e-3 if x < 2 else 1e-2
        X = arb(x); f = 3*(1 + 4/X**2).log()/X.cosh()**2 + 3/(X.cosh()**2*(X.cosh() + 1))
        tot += f*step; x += step
    # tail: e_up <= (12/x^2 + 6) 4 e^{-2x} for x >= xmax
    tail = arb(4)*(arb(12)/xmax**2 + 6)*arb(-2*xmax).exp()/2
    return hi(tot + tail)

def region_S(l0, K):
    worst = -1e9; stack = [(A_LO, A_HI, 0.0, l0)]; nb = 0
    while stack:
        a0, a1, x0, x1 = stack.pop(); nb += 1; assert nb < 400000
        al, x = iv(a0, a1), iv(x0, x1)
        s2 = al.sin()**2; c = x.cosh(); X = c - 1
        q_lo = 2/arb(x1).cosh(); ratio = arb((float(q_lo.mid()) + 2)/2, (2 - float(q_lo.mid()))/2 + 1e-12)
        Xs = X if x0 > 0 else arb(hi(X)/2, hi(X)/2 + 1e-30)
        T = 3*(Xs + s2)**2 - s2*s2
        sing = c*T/(2*(Xs + 2*s2)) - (Xs + 2)*(Xs + s2)*T/(Xs + 2*s2)**2          # multiplies 1/X
        reg = 3*c/2 + 6*(Xs + 2)*(Xs + s2)/(Xs + 2*s2)                             # regular part
        P = (ratio*sing + x*x*reg)/s2
        val = P + K*x/s2
        if val < 0: worst = max(worst, hi(val)); continue
        if (a1 - a0) > (x1 - x0)*4: m = (a0 + a1)/2; stack += [(a0, m, x0, x1), (m, a1, x0, x1)]
        else: m = (x0 + x1)/2; stack += [(a0, a1, x0, m), (a0, a1, m, x1)]
    return nb, worst

if __name__ == '__main__':
    t0 = time.time()
    ck = checks(); print('checks:', ck)
    L1 = 5.0; rEb = region_L(L1)
    print(f'L (l\' >= {L1}): r1*Eb = {rEb:.4f} < -c4max = {-ck["c4max"]:.4f}: {rEb < -ck["c4max"]}')
    # sanity: D_up >= d (closed forms) at sample points
    mp.mp.dps = 30
    def d_true(a, x):
        s2 = mp.sin(a)**2; c = mp.cosh(x); X = c - 1; C = 1 + X/s2; Q0 = mp.atanh(1/C)
        W1 = 3*C*Q0 - (3*C**2 - 1)/(2*(C**2 - 1)) - mp.mpf(3)/2
        W2 = 3*Q0 - 6*C/(C**2 - 1) + C*(3*C**2 - 1)/(C**2 - 1)**2
        return -W1*c - W2*(c**2 - 1)/s2
    def d_up(a, x):
        s2 = mp.sin(a)**2; c = mp.cosh(x); X = c - 1; T = 3*(X + s2)**2 - s2**2
        return (c*T/(2*(X + 2*s2)) - (X + 2)*(X + s2)*T/(X + 2*s2)**2)/X + 3*c/2 + 6*(X + 2)*(X + s2)/(X + 2*s2)
    print('sanity D_up >= d:', all(d_up(a, x) >= d_true(a, x) for a in [mp.pi/4, mp.mpf(1), mp.pi/2] for x in [mp.mpf(v) for v in ['0.001', '0.01', '0.1', '1', '3']]))
    K = K_integral(); print(f'K = int e_up <= {K:.4f}')
    l0 = 0.02
    nbS, w = region_S(l0, K); print(f'S (l\' <= {l0}): {nbS} boxes, max of P + K l\'/s^2 <= {w:.4f} < 0')
    nbM = region_M(l0, L1); print(f'M ({l0} <= l\' <= {L1}): {nbM} boxes, F_2(pi) < 0 on all')
    print(f'total {time.time() - t0:.0f} s')
