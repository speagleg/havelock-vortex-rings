"""Checks for the proof that hypercycle rows with kappa <= 1/sqrt 2 are unstable at EVERY spacing (tanh model).
Notation: K = kappa^2, s^2 = 1 - K, c_j = cosh(j l'), D(c) = c - (2K - 1) = c - cos(2 alpha) > 0.
Pairing mode (q = pi): the 2x2 Bloch symbol has zero diagonal, eigenvalues +-sqrt(A_UP(pi) A_PU(pi)), with
  A_PU(pi) = 4 sum_{j odd >= 1} d(c_j),  d(c) = -(c^2 + c - 1 - cos 2a) s^4 / ((c - 1) D^2) < 0,
  A_UP(pi) = 2 s^2 sum_{j >= 1} g_j,   g_odd(c) = ((2K-1)c - 1)/((c-1) D),  g_even(c) = ((4K-1)c - (4K^2-2K+1))/D^2.
Step 1: each g is strictly increasing in K (c > 1)  =>  sum g(K) < sum g(1/2) for K < 1/2.
Step 2: G(l') := sum g(1/2) = sum_m [-1/((c_{2m-1}-1) c_{2m-1}) + (c_{2m}-1)/c_{2m}^2] < 0 for every l' > 0:
  (a) l' >= 0.75:  G < -1/c_1^2 + sum_m sech(2 m l') <= -4q/(1+q)^2 + 2q/(1-q) < 0,  q = e^{-2l'} < sqrt5 - 2;
  (b) l' <= 1.2:   G < -1/((c_1-1)c_1) + pi/(4 l'),  and (cosh x - 1) cosh x < 4x/pi on (0, 1.2] (convexity).
Hence A_UP(pi) < 0 < -A_PU(pi), a real positive exponent: the row is linearly unstable.
Usage: python instability_proof_check.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings('ignore')
import numpy as np, sympy as sp
from pairing_root import AUP_pi

c, K, x = sp.symbols('c K x', positive=True)
D = c - (2*K - 1)
g_odd = ((2*K-1)*c - 1)/((c-1)*D)
g_even = ((4*K-1)*c - (4*K**2-2*K+1))/D**2

# 0. closed forms vs the independent finite-difference lattice sum
def A_closed(Kv, lp, J=60):
    s2 = 1 - Kv; tot = 0.0
    for j in range(1, J+1):
        cj = np.cosh(j*lp)
        if not np.isfinite(cj): break
        Dj = cj - (2*Kv - 1)
        tot += ((2*Kv-1)*cj - 1)/((cj-1)*Dj) if j % 2 else ((4*Kv-1)*cj - (4*Kv**2-2*Kv+1))/Dj**2
    return 2*s2*tot
worst = 0
for Kv in [0.1, 0.3, 0.45, 0.5, 0.6, 0.8]:
    for lp in [0.3, 0.7, 1.5, 3.0]:
        a = np.arccos(np.sqrt(Kv)); ex = AUP_pi(a, lp); cl = A_closed(Kv, lp)
        worst = max(worst, abs(ex-cl)/max(abs(ex), 1e-12))
print(f'0. closed form vs finite-difference lattice sum: max relative difference {worst:.1e}')

# 1. monotonicity in K
dgo = sp.factor(sp.simplify(sp.diff(g_odd, K)))
dge = sp.factor(sp.simplify(sp.diff(g_even, K)))
print('1. d g_odd/dK  =', dgo, '   (> 0 for c > 1)')
print('   d g_even/dK =', dge, '   (> 0 for c > 1, K < 1: (2c-1)(c+1) > 2 > 2K)')

# 2a. large spacing: -4q/(1+q)^2 + 2q/(1-q) < 0  <=>  q^2 + 4q - 1 < 0
q = sp.symbols('q', positive=True)
print('2a. 2(1-q) - (1+q)^2 =', sp.expand(2*(1-q) - (1+q)**2), '; root q =', sp.nsimplify(sp.solve(q**2+4*q-1, q)[0]),
      '-> l\' >', float(-sp.log(sp.sqrt(5)-2)/2))
# 2b. small spacing: phi(x) = (cosh x - 1) cosh x - 4x/pi convex with phi(0) = 0, phi(1.2) < 0
phi = (sp.cosh(x)-1)*sp.cosh(x) - 4*x/sp.pi
print('2b. phi(1.2) =', float(phi.subs(x, sp.Rational(6, 5))), '; phi\'\' =', sp.simplify(sp.diff(phi, x, 2)), '(> 0)')

# 3. direct numerical confirmation of G < 0 and of the reduction
lps = np.concatenate([np.geomspace(1e-3, 0.1, 50), np.linspace(0.1, 12, 600)])
G = [A_closed(0.5, lp)/(2*0.5) for lp in lps]
print(f'3. max over {len(lps)} spacings of G(l\') = {max(G):.3e} (must be < 0); G*e^(2l\') at l\'=12: {G[-1]*np.exp(24):.4f} (asymptote -2)')
