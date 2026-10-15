"""Robustness of the 'integrable-looking' verdict: finite-time Lyapunov exponents (T = 800) of the periodic Toda-limit
vortex row for many random initial conditions (N = 3 and 4), with a wider spread than before, vs periodic Toda.
Usage: python many_ics.py -> many_ics.out"""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'quantum'))
import numpy as np
from toda_test import nn_rhs, toda_rhs, lyap
out = open('many_ics.out', 'w')
for N in (3, 4):
    L = N*1.0; rng = np.random.default_rng(100 + N); vals = []
    for k in range(10):
        u0 = np.arange(N)*L/N + 0.35*rng.standard_normal(N); Y0 = rng.uniform(-3, 1.5, N)
        vals.append(lyap(nn_rhs(L), np.concatenate([u0, Y0]), T=800.0, seed=k))
        line = f'N {N} ic {k}: lambda(800) = {vals[-1]:.3e}'; print(line, flush=True); out.write(line + '\n'); out.flush()
    tv = [lyap(toda_rhs(L), np.concatenate([np.arange(N) + 0.4*rng.standard_normal(N), rng.standard_normal(N)]), T=800.0, seed=k) for k in range(3)]
    line = f'N {N}: row max {max(vals):.3e}, median {np.median(vals):.3e};  Toda controls {np.round(tv, 5)}'
    print(line, flush=True); out.write(line + '\n'); out.flush()
