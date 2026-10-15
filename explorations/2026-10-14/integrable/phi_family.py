"""Is H = sum_j phi(Y_j) e^{-(u_{j+1}-u_j)} (periodic ring, canonical u, Y) integrable for arbitrary phi?
The Toda-limit vortex row is gauge-equivalent (canonical shift u_j -> u_j + g(Y_j)) to phi = 2/(1+Y^2).
Finite-time Lyapunov exponents at T = 300 and 1200 for several phi. Usage: python phi_family.py -> phi_family.out"""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'quantum'))
import numpy as np
from toda_test import lyap
PHIS = {'2/(1+Y^2) (vortex row)': (lambda Y: 2/(1 + Y**2), lambda Y: -4*Y/(1 + Y**2)**2),
        'exp(Y)': (np.exp, np.exp),
        '1 + Y^2': (lambda Y: 1 + Y**2, lambda Y: 2*Y),
        '1 + 0.5 sin(2Y) + 0.1 Y^2': (lambda Y: 1 + 0.5*np.sin(2*Y) + 0.1*Y**2, lambda Y: np.cos(2*Y) + 0.2*Y),
        'cosh(Y)': (np.cosh, np.sinh)}
def rhs_phi(phi, dphi, L):
    def f(t, z):
        N = len(z)//2; u, Y = z[:N], z[N:]; D = np.roll(u, -1) - u; D[-1] += L; e = np.exp(-D)
        w = phi(Y)*e
        dHdu = w - np.roll(w, 1); dHdY = dphi(Y)*e
        return np.concatenate([dHdY, -dHdu])
    return f
if __name__ == '__main__':
    N = 4; L = 4.0; out = open('phi_family.out', 'w')
    for name, (phi, dphi) in PHIS.items():
        res = []
        for seed in (1, 2):
            rng = np.random.default_rng(seed)
            z0 = np.concatenate([np.arange(N)*1.0 + 0.2*rng.standard_normal(N), rng.uniform(-1.5, 1.5, N)])
            res.append((lyap(rhs_phi(phi, dphi, L), z0, T=300.0), lyap(rhs_phi(phi, dphi, L), z0, T=1200.0)))
        line = f'{name:28s}: ' + '; '.join(f'lambda(300) = {a:.3e}, lambda(1200) = {b:.3e} (ratio {a/b:4.2f})' for a, b in res)
        print(line, flush=True); out.write(line + '\n'); out.flush()
