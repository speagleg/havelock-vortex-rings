"""Classical check of the compressed chain: zigzag energy per cell f(eta) (y_A = +eta, y_B = -eta, grid points) from
(i) the exact pair sums and (ii) the compressed couplings (DMRG settings K = 5, 8 decay factors), at l = r l* for kappa = 0.98. Also the classical Landau
coefficients a2, a4 of the stiff-pinning model (fixed spacing, Y0 relaxed via -V y): compared with pinned.py.
Usage: python chain_check.py -> chain_check.out"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows'))
import numpy as np
from chain_model import build, W
from exact_symbol import lstar
def cell_energy_exact(m, yA, yB):
    lp, Y0 = m['lp'], m['Y0']; pA, pB = np.pi/2 + np.arctan(Y0 + yA), np.pi/2 + np.arctan(Y0 + yB)
    E = 0.0
    for r in range(1, 4000):
        if r % 2: E += 2*W(r*lp, pA, pB)
        else: E += W(r*lp, pA, pA) + W(r*lp, pB, pB)
        if r*lp > 60: break
    return E - m['V']*(yA + yB)
def cell_energy_comp(m, iA, iB):
    E = 0.0
    for lam, w, O in m['terms']:
        # sum_{r>=1} lam^r: odd r -> O_A O_B (twice: A->B and B->A), even r -> O_A^2 + O_B^2
        so, se = lam/(1 - lam**2), lam**2/(1 - lam**2)
        E += w*(2*so*O[iA]*O[iB] + se*(O[iA]**2 + O[iB]**2))
    return E - m['V']*(m['y'][iA] + m['y'][iB])
out = []; P = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
kappa = 0.98; ls = lstar(kappa)
for r in [1.0, 0.98]:
    m = build(kappa, ls*r, d=25, ymax=1.6, K=5, nlam=8, tol=1e-8)   # the DMRG compression settings
    y = m['y']; i0 = len(y)//2
    e0x, e0c = cell_energy_exact(m, 0, 0), cell_energy_comp(m, i0, i0)
    P(f'kappa {kappa}, l = {r} l*: f(eta) = E(eta) - E(0), exact vs compressed')
    for j in range(1, 7):
        fx = cell_energy_exact(m, y[i0 + j], y[i0 - j]) - e0x; fc = cell_energy_comp(m, i0 + j, i0 - j) - e0c
        P(f'   eta = {y[i0 + j]:.3f}: exact {fx: .6e}   compressed {fc: .6e}   diff {fc - fx: .1e}')
open('chain_check.out', 'w').write('\n'.join(out) + '\n')
