"""Long-time and second-initial-condition check for the nearest-neighbour -ln cosh(d/2) chain and the Toda-limit row.
Usage: python cosh_long.py -> cosh_long.out"""
import numpy as np
from kernel_integrability import KERNELS, rhs_of
from toda_test import lyap, nn_rhs
N = 4; L = 4.0; out = open('cosh_long.out', 'w')
for seed in (7, 11):
    rng = np.random.default_rng(seed)
    u0 = np.arange(N)*L/N + 0.15*rng.standard_normal(N); Y0 = -1.2 + 0.6*rng.standard_normal(N); z0 = np.concatenate([u0, Y0])
    for T in (300.0, 1200.0, 4800.0):
        lc = lyap(rhs_of(KERNELS['cosh  -ln cosh(d/2)'], L), z0, T=T)
        ln_ = lyap(nn_rhs(L), z0, T=T) if seed == 11 else float('nan')
        line = f'seed {seed} T {T:6.0f}: -ln cosh(d/2) chain {lc:.3e}' + (f'   Toda-limit row {ln_:.3e}' if seed == 11 else '')
        print(line, flush=True); out.write(line + '\n'); out.flush()
