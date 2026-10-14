"""Tricritical point of the zigzag transition: the curvature kappa_t where the quartic a4 vanishes on the threshold
curve l*(kappa), for the fixed (l', Y0) and the relaxed (fixed forces) ensembles; and a6 there (> 0 needed).
Usage: python tricritical.py -> tricritical.out"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows'))
import numpy as np, mpmath as mp
from scipy.optimize import brentq
from landau import relaxed, e_cell
from exact_symbol import lstar
mp.mp.dps = 30
def a4s(kappa):
    ls = lstar(kappa); s = np.sqrt(1 - kappa**2)
    a2, a4, a4r, M = relaxed(ls*s, -kappa/s); return float(a4), float(a4r), ls
def a6(kappa):
    ls = lstar(kappa); s = np.sqrt(1 - kappa**2); lp = mp.mpf(ls*s); Y0 = mp.mpf(-kappa/s)
    f = lambda eta: e_cell(Y0 + eta, Y0 - eta, lp); return mp.diff(f, 0, 6)/720
out = []; P = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
kf = brentq(lambda k: a4s(k)[0], 0.90, 0.95, xtol=1e-7)
kr = brentq(lambda k: a4s(k)[1], 0.95, 0.98, xtol=1e-7)
for name, k in [('fixed (l\', Y0)', kf), ('relaxed (fixed forces)', kr)]:
    P(f'{name:24s}: kappa_t = {k:.6f},  l*(kappa_t) = {lstar(k):.6f},  a6(kappa_t, fixed) = {mp.nstr(a6(k), 6)},'
      f'  hypercycle distance artanh(kappa_t) = {np.arctanh(k):.5f}')
open('tricritical.out', 'w').write('\n'.join(out) + '\n')
