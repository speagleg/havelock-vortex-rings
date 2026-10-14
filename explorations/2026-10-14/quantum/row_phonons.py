"""Quantised phonons of a vortex row on a hypercycle (tanh model).
Darboux coordinates (u, Y = -cot phi) are canonical ({u_i, Y_i} = 1, area form du dY = dA/R^2), so quantisation is
[u_i, Y_j] = i eps delta_ij with eps = 1/(2 pi n R^2) for superfluid vortices (n = areal number density, R = curvature
radius; Haldane-Wu). Per Bloch wavenumber the quadratic form of H - V I (I = sum Y generates u-translation) is
    h(q) = [[ -C/s^2, -i delta ], [ i delta, B s^2 ]],   J h(q) = [[i delta, B s^2], [C/s^2, i delta]],
so there is ONE chiral branch (one degree of freedom per vortex):  E(q) = eps [ sqrt(-B C) + delta(q) ],  q in (-pi, pi].
Spectral stability: B >= 0 for all q.  Energetic stability (h(q) >= 0, i.e. a quantum ground state): E(q) >= 0 for all q,
i.e. sqrt(-B C) >= |delta|.  Near q = pi at the spectral threshold: E ~ (v |k| + delta'(pi) k), v = sqrt(|C(pi)| B''(pi)/2).
Near q = 0: E ~ (c0 |q| + delta'(0) q), c0 = sqrt(B(0) |C''(0)/2|).
Usage: python row_phonons.py -> row_phonons.out"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows'))
import numpy as np
from scipy.optimize import brentq
from symbol import symbol, derivs_at_pi, derivs_at_0
from exact_symbol import lstar

def Emin(kappa, lp, nq=2001):
    q = np.linspace(-np.pi, np.pi, nq)
    B, C, dl = symbol(kappa, lp, q)
    root = np.sqrt(np.maximum(-B*C, 0)); E = root + dl
    E[B < 0] = np.nan
    return np.nanmin(E) if np.isfinite(E).any() else np.nan, q, E

if __name__ == '__main__':
    out = []; P = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
    P('1. Checks: symbol from the closed forms vs the q -> pi and q -> 0 expansions')
    for kappa, l in [(0.8, 3.0), (0.95, 1.5)]:
        s = np.sqrt(1 - kappa**2); lp = l*s
        Bpp, Cpi, dpi = derivs_at_pi(kappa, lp); k = 1e-3
        B, C, dl = symbol(kappa, lp, np.array([np.pi - k, np.pi, np.pi + k]))
        P(f'   kappa {kappa}, l {l}: B\'\'(pi) fd {(B[0] - 2*B[1] + B[2])/k**2:.6f} vs {Bpp:.6f};  delta\'(pi) fd {(dl[2] - dl[0])/(2*k):.6f} vs {dpi:.6f}')
    P('\n2. At the spectral threshold l*(kappa): critical-mode velocity v vs the chiral shift delta\'(pi), and the q -> 0 sound')
    P('   kappa     l*        v = sqrt(|C| B\'\'/2)   delta\'(pi)    |delta\'|/v    c0 (q->0)   delta\'(0)   |delta\'(0)|/c0')
    rows = []
    for kappa in [0.72, 0.75, 0.8, 0.85, 0.9, 0.95, 0.98, 0.995]:
        ls = lstar(kappa); s = np.sqrt(1 - kappa**2); lp = ls*s
        Bpp, Cpi, dpi = derivs_at_pi(kappa, lp); v = np.sqrt(abs(Cpi)*Bpp/2)
        B0, c2, d0 = derivs_at_0(kappa, lp); c0 = np.sqrt(B0*abs(c2))
        rows.append((kappa, ls, v, dpi, c0, d0))
        P(f'   {kappa:5.3f}  {ls:8.5f}   {v:14.6e}     {dpi:12.5e}   {abs(dpi)/v:8.4f}   {c0:10.4e}  {d0:11.4e}   {abs(d0)/c0:8.4f}')
    P('\n3. Energetic (ground-state) threshold l_E(kappa) = smallest spacing with E(q) >= 0 for all q, vs spectral l*(kappa)')
    for kappa in [0.72, 0.75, 0.8, 0.85, 0.9, 0.95, 0.98, 0.995]:
        ls = lstar(kappa); s = np.sqrt(1 - kappa**2)
        g = lambda l: Emin(kappa, l*s)[0]
        ls_grid = np.linspace(ls*1.0005, ls*6, 300); vals = np.array([g(x) for x in ls_grid])
        ok = np.where(vals >= 0)[0]
        if len(ok) == 0:
            P(f'   kappa {kappa:5.3f}: l* = {ls:.5f}; E(q) has negative values for all l up to {ls_grid[-1]:.3f}'); continue
        i = ok[0]
        lE = ls_grid[0] if i == 0 else brentq(g, ls_grid[i - 1], ls_grid[i])
        allok = np.all(vals[i:] >= 0)
        P(f'   kappa {kappa:5.3f}: spectral l* = {ls:.6f}, energetic l_E = {lE:.6f}  (ratio {lE/ls:.4f}; E >= 0 for all larger l on grid: {allok})')
    open('row_phonons.out', 'w').write('\n'.join(out) + '\n')
