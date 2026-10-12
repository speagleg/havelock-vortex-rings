# Large-N limit of the tanh-model ring threshold: checks behind note_tanh.tex, Section "Large N".
# Everything here is a numerical check of exact statements proved in the note (PROVED there), except the parts
# marked NUMERICAL, which are evidence only.
#   python large_n.py > large_n.out      (about 2 minutes)
import random
import mpmath as mp
import numpy as np
import sympy as sp
from blocks import blocks

mp.mp.dps = 60
L = lambda y: mp.coth(y) - 1/y                         # Langevin function


def Q(x, s):
    return mp.cosh(s)*mp.cosh(x*s) - x*mp.sinh(s)*mp.sinh(x*s)


def symbol(N, m, xi):
    """The exact symbol form of Theorem (exact symbol): returns a_m/N^2, b_m/N^2, g_m/N^2."""
    s = -mp.mpf(N)/2*mp.log(xi); x = 1 - mp.mpf(2*m)/N; k = mp.coth(s/N)/N; sh2 = mp.sinh(s)**2
    a = k*mp.coth(s) - k**2/2 - 1/(2*mp.mpf(N)**2) - (1-x**2)/8 - (1+Q(x, s))/(4*sh2)
    b = (1-x**2)/8 + (1-Q(x, s))/(4*sh2)
    th = (1-x)/2
    g = ((th*mp.exp(x*s) - (1-th)*mp.exp(-x*s))/(2*mp.sinh(s)) + mp.exp(-s)*mp.sinh(x*s)/(2*sh2))/2
    return a, b, g


def limit_symbol(x, s):
    sh2 = mp.sinh(s)**2
    a = mp.coth(s)/s - 1/(2*s**2) - (1-x**2)/8 - (1+Q(x, s))/(4*sh2)
    b = (1-x**2)/8 + (1-Q(x, s))/(4*sh2)
    th = (1-x)/2
    g = ((th*mp.exp(x*s) - (1-th)*mp.exp(-x*s))/(2*mp.sinh(s)) + mp.exp(-s)*mp.sinh(x*s)/(2*sh2))/2
    return a, b, g


Phi = lambda s: (2*mp.tanh(s/2)/(s/2) - L(s/2)**2)/8

print("1. Exact symbol form (Theorem: exact symbol), 400 random (N, m, xi), N <= 300")
random.seed(20261012); worst = 0
for _ in range(400):
    N = random.randint(2, 300); m = random.randint(1, N-1)
    xi = mp.mpf(random.randint(1, 10**6 - 1))/10**6
    A, B, G, _ = blocks(N, xi, m); a, b, g = symbol(N, m, xi)
    worst = max(worst, abs(A/N**2 - a)/(1 + abs(a)), abs(B/N**2 - b)/(1 + abs(b)), abs(G/N**2 - g)/(1 + abs(g)))
print("   max relative difference from blocks.py:", mp.nstr(worst, 3))

print("\n2. The constant")
yinf = mp.findroot(lambda y: L(y)**2 - 2*mp.tanh(y)/y, 3.7); sinf = 2*yinf
kinf = 1/sinf
print("   y_inf   (L(y)^2 = 2 tanh(y)/y)     =", mp.nstr(yinf, 40))
print("   kappa_inf = 1/(2 y_inf)             =", mp.nstr(kinf, 40))
print("   ell* = 2 pi kappa_inf = pi/y_inf    =", mp.nstr(mp.pi/yinf, 40))
u = mp.e**(-1/kinf)
old = kinf - kinf**2/2 + 2*kinf*u**2/(1-u**2) - mp.mpf(1)/8 - u/(2*(1-u)**2)
print("   residual of the note's earlier eq. (kappa) at kappa_inf:", mp.nstr(old, 3))
print("   algebraic approximation y = 2+sqrt3: kappa = 1-sqrt3/2 =", mp.nstr(1-mp.sqrt(3)/2, 12),
      ", ell ~ pi(2-sqrt3) =", mp.nstr(mp.pi*(2-mp.sqrt(3)), 12))
print("   Phi(1e-8) =", mp.nstr(Phi(mp.mpf('1e-8')), 12), "(limit 1/4);  Phi(200) =", mp.nstr(Phi(mp.mpf(200)), 12), "(limit -1/8);",
      " Phi'(s_inf) =", mp.nstr(mp.diff(Phi, sinf), 12))

print("\n3. Even-N closed form 8 a_{N/2}/N^2 = 2 mu tanh y - (coth y - mu)^2 - 4/N^2, mu = (2/N) coth(2y/N)")
worst = 0
for N in [2, 4, 8, 10, 16, 50, 400]:
    for xi in ['0.01', '0.2', '0.5', '0.9', '0.999']:
        xi = mp.mpf(xi); y = -mp.mpf(N)/4*mp.log(xi); mu = 2/mp.mpf(N)*mp.coth(2*y/N)
        worst = max(worst, abs(8*blocks(N, xi, N//2)[0]/N**2 - (2*mu*mp.tanh(y) - (mp.coth(y)-mu)**2 - 4/mp.mpf(N)**2)))
print("   max difference from blocks.py:", mp.nstr(worst, 3))
for N in [6, 7, 8, 9]:
    m = N//2; lim0 = 8*blocks(N, mp.mpf('1e-30'), m)[0]/N**2
    exact = -mp.mpf(N*N-8*N+8)/N**2 if N % 2 == 0 else -mp.mpf((N-1)*(N-7))/N**2
    print(f"   planar limit xi -> 0, N = {N}: 8 a_crit/N^2 = {mp.nstr(lim0, 12)}, predicted {mp.nstr(exact, 12)}")

print("\n4. Thresholds against the second-order expansion cosh r* = kappa_inf N + C1/N + O(N^-3)")
dPhi = mp.diff(Phi, sinf)
D_even = (sinf*mp.coth(sinf) - 1)/3 - mp.mpf(1)/2
psixx = -mp.mpf(1)/4 + (sinf**2*mp.cosh(sinf) - 2*sinf*mp.sinh(sinf))/(4*mp.sinh(sinf)**2)
D_odd = D_even - psixx/2
C1 = {0: sinf/3 + D_even/(dPhi*sinf**2), 1: sinf/3 + D_odd/(dPhi*sinf**2)}
print("   C1(even) =", mp.nstr(C1[0], 15), "  C1(odd) =", mp.nstr(C1[1], 15))
print("      N   (cosh r* - kappa_inf N) N    C1(parity)")
for N in [8, 9, 16, 17, 50, 51, 100, 101, 200, 201, 400, 401, 800, 801]:
    m = N//2; s0 = sinf + (D_even if N % 2 == 0 else D_odd)/(-dPhi)/N**2
    xs = mp.findroot(lambda x: blocks(N, x, m)[0], mp.e**(-2*s0/N))
    c = (1+xs)/(1-xs)
    print(f"   {N:4d}   {mp.nstr((c - kinf*N)*N, 15):>20}   {mp.nstr(C1[N % 2], 15)}")

print("\n5. Odd N: the extra term P3 of 8a/N^2 against the strictly decreasing part P12 (NUMERICAL margin)")
mp.mp.dps = 80


def parts(N, v):
    y = N*v/2; mu = 2/mp.mpf(N)*mp.coth(v)
    P12 = 2*mu*mp.tanh(y) - (mp.coth(y)-mu)**2
    P3 = (2*mp.sinh(v)*mp.sinh(N*v)/N - 2*mp.cosh(N*v)*(mp.cosh(v)-1))/mp.sinh(N*v)**2
    return P12, P3


for N in [3, 5, 7, 9, 11, 15, 21, 51, 101, 1001]:
    r = 0
    for k in range(-60, 61):
        v = mp.mpf(10)**(mp.mpf(k)/20)/N*4
        if N*v > 60:
            continue
        d12 = mp.diff(lambda t: parts(N, t)[0], v); d3 = mp.diff(lambda t: parts(N, t)[1], v)
        assert d12 < 0
        if d3 > 0:
            r = max(r, d3/(-d12))
    print(f"   N={N:5d}: max_v P3'/|P12'| (where P3' > 0) = {mp.nstr(r, 4)}")
mp.mp.dps = 60

print("\n6. Infinite vortex row on the horocycle y = 1 (upper half-plane), spacing ell (NUMERICAL, K = 2e5 terms)")
x0, y0, x1, y1 = sp.symbols('x0 y0 x1 y1', real=True)
gpair = -(sp.log(sp.sqrt((x0-x1)**2+(y0-y1)**2)) - sp.log(sp.sqrt((x0-x1)**2+(y0+y1)**2)))
qv = [x0, y0, x1, y1]
grad = sp.lambdify(qv, [sp.diff(gpair, v) for v in qv], 'numpy')
hess = sp.lambdify(qv, [[sp.diff(gpair, a, b) for b in qv] for a in qv], 'numpy')


def row_symbol(ell, theta, K=200000):
    k = np.concatenate([np.arange(1, K+1), -np.arange(1, K+1)]).astype(float)
    X1, Y1 = k*ell, np.ones_like(k)
    v = np.sum(grad(0.0, 1.0, X1, Y1)[1]) + 4.0/(K*ell**2)     # force balance, plus the tail of sum 2/((k ell)^2+4)
    Hh = [[np.broadcast_to(np.asarray(h, float), k.shape) for h in row] for row in hess(0.0, 1.0, X1, Y1)]
    onsite = np.array([[np.sum(Hh[a][b]) for b in range(2)] for a in range(2)])
    onsite[1, 1] += 2.0*v                                       # -v J with J = -sum 1/y
    ph = np.exp(2j*np.pi*theta*k)
    return onsite + np.array([[np.sum(ph*Hh[a][b+2]) for b in range(2)] for a in range(2)])


worst = 0
for ell in [0.5, float(mp.pi/yinf), 1.5, 3.0]:
    for theta in [0.5, 0.3, 0.1]:
        M = row_symbol(ell, theta)*(ell/(2*np.pi))**2
        a, b, g = limit_symbol(mp.mpf(1-2*theta), 2*mp.pi/ell)
        worst = max(worst, abs(M[1, 1].real - float(a)), abs(M[0, 0].real - float(b)), abs(M[0, 1].imag - float(g)))
print("   max |(ell/2pi)^2 Hess(row) - limit symbol (alpha, beta, gamma)| over ell in {0.5, ell*, 1.5, 3},"
      " theta in {0.5, 0.3, 0.1}:", f"{worst:.1e}")
print("   radial entry at ell = ell*, theta = 1/2:", f"{row_symbol(float(mp.pi/yinf), 0.5)[1, 1].real:.1e}")

print("\n7. Limit symbol off the critical mode, s <= s_inf: min (alpha beta - gamma^2)/beta (NUMERICAL)")
worst = mp.inf
for s in [sinf*mp.mpf(k)/40 for k in range(1, 41)]:
    for i in range(1, 100):
        a, b, g = limit_symbol(mp.mpf(i)/100, s); worst = min(worst, (a*b - g*g)/b)
print("   min over s = s_inf k/40, x = i/100:", mp.nstr(worst, 5))
h = mp.mpf('1e-8'); a0, b0, _ = limit_symbol(0, sinf); a1, _, g1 = limit_symbol(h, sinf)
W = 2*(g1/h)**2/(b0*(-dPhi))
print("   predicted odd-N Krein window xi_EM - xi_* ~ W/N^3 with W =", mp.nstr(W, 6))
for N in [9, 15, 31, 61, 121]:
    m = (N-1)//2; s0 = sinf + D_odd/(-dPhi)/N**2
    xs = mp.findroot(lambda x: blocks(N, x, m)[0], mp.e**(-2*s0/N))

    def Dm(x):
        a, b, g, _ = blocks(N, x, m); return a*b - g*g
    xe = mp.findroot(Dm, xs + W/N**3)
    print(f"   N={N:4d}: N^3 (xi_EM - xi_*) = {mp.nstr(N**3*(xe-xs), 6)}")

print("\n8. Energy-momentum threshold (Lemma 17, Theorem 18)")
mp.mp.dps = 60


def psi(x, t): return (1-x**2)/8 + Q(x, t)/(4*mp.sinh(t)**2)


def Kf(x, t):
    c = mp.cosh(t)
    return c*(1-mp.cosh(x*t)) + x*mp.sinh(t)*mp.sinh(x*t) - x**2*(c-1)/2*(2-(1-x**2)*(c-1))


worst = 0
for x, t in [(mp.mpf('0.3'), mp.mpf(2)), (mp.mpf('0.7'), mp.mpf(9)), (mp.mpf('0.05'), mp.mpf('0.4')), (mp.mpf('0.99'), mp.mpf(30))]:
    a, b, g = limit_symbol(x, t); G = psi(0, t) - psi(x, t) - g*g/b
    worst = max(worst, abs(b*G - (mp.cosh(t)+1)**2*Kf(x, t)/(32*mp.sinh(t)**4))/abs(b*G))
print("   identity beta G = (cosh+1)^2 K/(32 sinh^4), max rel err:", mp.nstr(worst, 3))
ts = [mp.mpf(10)**(mp.mpf(e)/10) for e in range(-30, 21)]; xg = [mp.mpf(i)/100 for i in range(1, 100)]
print("   K > 0 on the grid tau = 10^(-3..2), x = 0.01..0.99:", all(Kf(x, t) > 0 for t in ts for x in xg))
print("   even N: min_m (a_m b_m - g_m^2) just above / below xi* (EM threshold = spectral threshold)")
for N in [8, 20, 50, 100]:
    s0 = sinf + D_even/(-dPhi)/N**2
    xs = mp.findroot(lambda x: blocks(N, x, N//2)[0], mp.e**(-2*s0/N))

    def minD(xi):
        return min(blocks(N, xi, m)[0]*blocks(N, xi, m)[1] - blocks(N, xi, m)[2]**2 for m in range(1, N//2+1))
    e = mp.mpf('1e-30')
    print(f"   N={N:4d}: above {mp.nstr(minD(xs+e), 3)}   below {mp.nstr(minD(xs-e), 3)}")
print("   odd N: xi* < xi_EM <= xi(first zero of Phi_N), the bound of Theorem 18(b)")
for N in [9, 15, 31, 61, 121, 241]:
    m = (N-1)//2

    def PhiN(t):
        k = mp.coth(t/N)/N
        return k*mp.coth(t) - k**2/2 - 1/(2*mp.mpf(N)**2) - 1/(4*mp.sinh(t)**2) - psi(0, t)
    tPhi = mp.findroot(PhiN, sinf)
    tspec = mp.findroot(lambda t: PhiN(t) + psi(0, t) - psi(mp.mpf(1)/N, t), sinf)
    xs, xb = mp.e**(-2*tspec/N), mp.e**(-2*tPhi/N)

    def Dm(xi):
        a, b, g, _ = blocks(N, xi, m); return a*b - g*g
    xe = mp.findroot(Dm, (xs+xb)/2)
    print(f"   N={N:4d}: N^3 (xi_EM - xi*) = {mp.nstr(N**3*(xe-xs), 6)}   N^3 (bound - xi*) = {mp.nstr(N**3*(xb-xs), 6)}"
          f"   xi* < xi_EM <= bound: {xs < xe <= xb}")
print("   odd N: on 9 interior points of (xi*, xi_EM): critical block indefinite, all other blocks definite (NUMERICAL)")
for N in [17, 21, 31, 61, 121, 241]:
    m0 = (N-1)//2

    def PhiN(t):
        k = mp.coth(t/N)/N
        return k*mp.coth(t) - k**2/2 - 1/(2*mp.mpf(N)**2) - 1/(4*mp.sinh(t)**2) - psi(0, t)
    tPhi = mp.findroot(PhiN, sinf)
    tspec = mp.findroot(lambda t: PhiN(t) + psi(0, t) - psi(mp.mpf(1)/N, t), sinf)
    xs, xb = mp.e**(-2*tspec/N), mp.e**(-2*tPhi/N)

    def Dm(xi, m):
        a, b, g, _ = blocks(N, xi, m); return a*b - g*g
    xe = mp.findroot(lambda xi: Dm(xi, m0), (xs+xb)/2)
    ok = True
    for j in range(1, 10):
        xi = xs + (xe-xs)*j/10
        ok &= Dm(xi, m0) < 0 and all(Dm(xi, m) > 0 and blocks(N, xi, m)[0] > 0 for m in range(1, N//2+1) if m != m0)
    print(f"   N={N:4d}: {ok}")
