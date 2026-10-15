"""Long-time check: finite-time Lyapunov estimates vs T for the NN vortex row and the periodic Toda lattice
(integrable: decays ~ ln T / T) and the full tanh row (chaotic: plateaus). Usage: python toda_long.py -> toda_long.out"""
import numpy as np
from toda_test import nn_rhs, toda_rhs, full_rhs, lyap
N = 4; L = N*1.0; rng = np.random.default_rng(7)
u0 = np.arange(N)*L/N + 0.15*rng.standard_normal(N); Y0 = -1.2 + 0.6*rng.standard_normal(N); z0 = np.concatenate([u0, Y0])
q0 = np.arange(N)*1.0 + 0.4*rng.standard_normal(N); p0 = 0.8*rng.standard_normal(N); w0 = np.concatenate([q0, p0])
out = open('toda_long.out', 'w')
for T in [300, 1200, 4800]:
    line = f'T = {T:5d}: NN vortex row {lyap(nn_rhs(L), z0, T=T):.3e}   periodic Toda {lyap(toda_rhs(L), w0, T=T):.3e}'
    print(line, flush=True); out.write(line + '\n'); out.flush()
line = f'T =  1200: full tanh row {lyap(full_rhs(L, N), z0, T=1200):.3e}'; print(line, flush=True); out.write(line + '\n')
