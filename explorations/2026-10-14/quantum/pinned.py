"""Pinned spacing: a periodic potential along the hypercycle pins u_j, so the long-wavelength sector (u, Y conjugate)
is gapped and couples to eta^2 only locally. At frequencies below the gap Y0 follows adiabatically at fixed spacing:
    a4_pinned = a4_fixed - g_Y^2 / (2 M_YY),   g_Y = d a2/d Y0,  M_YY = d^2 e0/d Y0^2  (per cell).
kappa_t(pinned) = zero of a4_pinned on the threshold curve. With no gapless mode left, the eta^4-eta^6 theory is the
standard Z2 one and its tricritical point is the tricritical Ising point (if the transition there is continuous).
Usage: python pinned.py -> pinned.out"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows'))
import numpy as np, mpmath as mp
from scipy.optimize import brentq
from landau import coeffs, e_cell
from exact_symbol import lstar
mp.mp.dps = 30
def a4_pinned(kappa, full=False):
    ls = lstar(kappa); s = np.sqrt(1 - kappa**2); lp, Y0 = mp.mpf(ls*s), mp.mpf(-kappa/s)
    a2, a4 = coeffs(lp, Y0)
    gY = mp.diff(lambda y: coeffs(lp, y)[0], Y0); MYY = mp.diff(lambda y: e_cell(y, y, lp), Y0, 2)
    val = a4 - gY**2/(2*MYY)
    return (float(a4), float(gY), float(MYY), float(val)) if full else float(val)
out = []; P = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
P(' kappa    a4_fixed     g_Y        M_YY      a4_pinned')
for k in [0.90, 0.92, 0.93, 0.95, 0.98]:
    a4, gY, M, v = a4_pinned(k, True); P(f' {k:5.3f}  {a4:10.3e}  {gY:10.3e}  {M:9.4f}  {v:10.3e}')
kt = brentq(a4_pinned, 0.95, 0.98, xtol=1e-7)
P(f'kappa_t (pinned spacing) = {kt:.6f}, l*(kappa_t) = {lstar(kt):.6f}; compare fixed 0.920274, relaxed 0.960198')
open('pinned.out', 'w').write('\n'.join(out) + '\n')
