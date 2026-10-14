"""Hypercycle rows (tanh model): the critical Bloch mode is always the pairing mode q = pi, and the exact threshold.
With z = e^{-x}, c = cosh x, s = sin(alpha), K = cos^2(alpha), the pair functions of exact_symbol.py expand as
  b(c) =  2 s^2 sum_{n>=1} n cos^2(n alpha) z^n          (all coefficients >= 0)
  d(c) = -2 s^2 sum_{n>=1} n sin^2(n alpha) z^n
  a(c) = d(c) + sin^2(2 alpha) / (2 (c - cos 2alpha)) = sum_n [sin(2 alpha) sin(2n alpha) - 2 n s^2 sin^2(n alpha)] z^n,
equivalently b - d = s^2/(c - 1).  Hence with r_n = e^{-n l'}:
  sum_j b(c_j) cos(qj) = 2 s^2 sum_n n cos^2(n alpha) P(r_n, q),  P(r, q) = sum_j r^j cos(qj) = (r cos q - r^2)/(1 - 2 r cos q + r^2),
and dP/dq = -r (1 - r^2) sin q / (1 - 2 r cos q + r^2)^2 < 0 on (0, pi).  So F(q) is strictly decreasing on (0, pi):
  * the row is stable iff F(pi) >= 0 (critical mode = pairing mode, for every kappa and spacing);
  * a periodic necklace of N vortices (allowed q = 2 pi m/N) is stable iff F(2 pi floor(N/2)/N) >= 0.
Threshold function (Lambert-type series, converges like e^{-n l'}):
  F(pi) = sum_n [ (sin 2a sin 2na - 2n s^2 sin^2 na) r_n/(1 - r_n) - 2n s^2 cos^2(na) r_n/(1 + r_n) ].
Usage: python qpi_proof.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings('ignore')
import numpy as np, sympy as sp
from scipy.optimize import brentq
from exact_symbol import ab, Fmin

# 1. generating-function identities, exact in E = e^{i alpha} (all functions rational in z and E)
z, E = sp.symbols('z E', nonzero=True)
ca, sa = (E + 1/E)/2, (E - 1/E)/(2*sp.I)
K = ca**2; s2 = sa**2; c = (z + 1/z)/2; D = c + 1 - 2*K; cos2a = 2*K - 1; sin2a = 2*sa*ca
a = s2*((3*K-1)*c**2 - (4*K**2-K+1)*c + 2*K**2)/((c-1)*D**2)
b = s2*(K*c**2 + (1-3*K)*c + 2*K**2-2*K+1)/((c-1)*D**2)
d = -s2**2*(c**2 + c - 2*K)/((c-1)*D**2)
w1, w2 = z*E**2, z/E**2
S_n = z/(1-z)**2                                    # sum n z^n
S_cos = (w1/(1-w1)**2 + w2/(1-w2)**2)/2             # sum n cos(2 n alpha) z^n
S_sin = sin2a*z/(1 - 2*z*cos2a + z**2)              # sum sin(2 n alpha) z^n
chk = lambda e: sp.cancel(sp.together(e))
print('1. b - 2s^2 sum n cos^2(na) z^n      =', chk(b - s2*(S_n + S_cos)))
print('   d + 2s^2 sum n sin^2(na) z^n      =', chk(d + s2*(S_n - S_cos)))
print('   a - d - sin(2a) sum sin(2na) z^n  =', chk(a - d - sin2a*S_sin))
print('   b - d - s^2/(c-1)                 =', chk(b - d - s2/(c-1)))
print('   a - d - sin^2(2a)/(2(c - cos2a))  =', chk(a - d - sin2a**2/(2*(c - cos2a))))
q, r = sp.symbols('q r', positive=True)
P = (r*sp.cos(q) - r**2)/(1 - 2*r*sp.cos(q) + r**2)
print('   dP/dq =', sp.factor(sp.simplify(sp.diff(P, q))))

# 2. the Lambert series for F(pi) vs the direct lattice sum
def Fpi_series(kappa, lp, nmax=None):
    a_ = np.arccos(kappa); s2_ = np.sin(a_)**2
    nmax = nmax or int(np.ceil(80/lp)) + 5
    n = np.arange(1, nmax+1); rn = np.exp(-n*lp)
    A = np.sin(2*a_)*np.sin(2*n*a_) - 2*n*s2_*np.sin(n*a_)**2
    B = 2*n*s2_*np.cos(n*a_)**2
    return np.sum(A*rn/(1-rn) - B*rn/(1+rn))
def Fq_direct(kappa, lp, q):
    J = max(4, int(np.ceil(60/lp))); j = np.arange(1, J+1); cj = np.cosh(np.minimum(j*lp, 700.0))
    A, B = ab(kappa**2, cj); return A.sum() + (B*np.cos(q*j)).sum()
worst = 0
for kappa in [0.3, 0.71, 0.8, 0.95, 0.999]:
    for lp in [0.05, 0.3, 1.0, 3.0]:
        x, y = Fpi_series(kappa, lp), Fq_direct(kappa, lp, np.pi); worst = max(worst, abs(x-y)/max(abs(y), 1e-300))
print(f'2. Lambert series vs direct lattice sum for F(pi): max relative difference {worst:.1e}')

# 3. monotonicity in q, numerically, on a grid
bad = 0
for kappa in np.linspace(0.05, 0.999, 12):
    for lp in np.geomspace(0.02, 10, 15):
        qs = np.linspace(0.01, np.pi, 100); F = [Fq_direct(kappa, lp, x) for x in qs]
        bad += np.sum(np.diff(F) > 1e-12*max(1, abs(F[0])))
print(f'3. grid points where F(q) increases: {bad} (of 12 x 15 x 99)')

# 4. exact threshold l*(kappa) = root of F(pi) (arc length), and uniqueness of the root on a fine grid
print('4. l*(kappa) from F(pi) = 0; number of sign changes of F(pi) on l\' in [0.005, 40]')
for kappa in [0.72, 0.75, 0.8, 0.9, 0.95, 0.99, 0.999, 0.9999]:
    s = np.sqrt(1-kappa**2); g = np.geomspace(0.005, 40, 1500); v = np.array([Fpi_series(kappa, x) for x in g])
    ch = np.where(np.diff(np.sign(v)) != 0)[0]
    root = brentq(lambda x: Fpi_series(kappa, x), g[ch[0]], g[ch[0]+1], xtol=1e-14)/s if len(ch) else np.nan
    print(f'   kappa {kappa:7.4f}: l* = {root:.10f}   sign changes {len(ch)}')

# 5. odd-N necklaces: critical mode q = pi (N-1)/N; Bolza collar thresholds vs collar.out
L = 2*np.arccosh(1+np.sqrt(2))
print('5. Bolza necklaces: threshold h* from F(2 pi floor(N/2)/N) = 0 at l\' = L/N')
for N in range(2, 13):
    qN = 2*np.pi*(N//2)/N; lp = L/N
    f = lambda h: Fq_direct(np.tanh(h), lp, qN)
    print(f'   N = {N:2d}: h* = {brentq(f, 0.5, 4.0, xtol=1e-12):.4f}')
