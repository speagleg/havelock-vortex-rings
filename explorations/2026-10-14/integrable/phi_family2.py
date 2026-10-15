"""Characterise the integrable phi in H = sum phi(Y_j) e^{-(u_{j+1}-u_j)}: is phi = 1/Q(Y), Q quadratic, integrable for
every Q (complex, double, real roots)? Control 1/(1+Y^2)^2. Usage: python phi_family2.py -> phi_family2.out"""
import numpy as np
from phi_family import rhs_phi
from toda_test import lyap
PHIS = {'1/(Y^2 + 0.6 Y + 2)':   (lambda Y: 1/(Y**2 + 0.6*Y + 2), lambda Y: -(2*Y + 0.6)/(Y**2 + 0.6*Y + 2)**2, (-1.5, 1.5)),
        '1/Y^2 (Y > 0)':         (lambda Y: 1/Y**2, lambda Y: -2/Y**3, (1.5, 3.0)),
        '1/(Y^2 - 1) (Y > 1)':   (lambda Y: 1/(Y**2 - 1), lambda Y: -2*Y/(Y**2 - 1)**2, (2.0, 3.5)),
        '1/(1+Y^2)^2 (control)': (lambda Y: 1/(1 + Y**2)**2, lambda Y: -4*Y/(1 + Y**2)**3, (-1.5, 1.5)),
        'Y^2 (Y > 0)':           (lambda Y: Y**2, lambda Y: 2*Y, (1.5, 3.0))}
N = 4; L = 4.0; out = open('phi_family2.out', 'w')
for name, (phi, dphi, rngY) in PHIS.items():
    res = []
    for seed in (1, 2):
        rng = np.random.default_rng(seed)
        z0 = np.concatenate([np.arange(N)*1.0 + 0.2*rng.standard_normal(N), rng.uniform(*rngY, N)])
        res.append((lyap(rhs_phi(phi, dphi, L), z0, T=300.0), lyap(rhs_phi(phi, dphi, L), z0, T=1200.0)))
    line = f'{name:24s}: ' + '; '.join(f'{a:.3e} -> {b:.3e} (ratio {a/b:4.2f})' for a, b in res)
    print(line, flush=True); out.write(line + '\n'); out.flush()
