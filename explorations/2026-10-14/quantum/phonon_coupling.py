"""Coupling of the zigzag order parameter to the long-wavelength (acoustic) sector, at the two mean-field tricritical
points and along the threshold.
 - couplings g = d a2 / d(l', Y0) of eta^2 to uniform strain and uniform impulse density;
 - long-wavelength stiffness matrix M (Hessian of the row energy per cell in (l', Y0));
 - Larkin-Pikin shift Delta = g^T M^{-1} g / 2 = a4_fixed - a4_relaxed, and its size relative to |a4_fixed|;
 - velocities: critical (soft) mode v = sqrt(|C(pi)| B''(pi)/2) with drift delta'(pi); sound c0 = sqrt(B(0)|C''(0)/2|)
   with drift delta'(0). The induced eta^4 interaction is retarded: it equals the relaxed value for frequencies
   << c k and vanishes for frequencies >> c k. If the sound is much faster than the critical mode, the relaxed
   ensemble is the right mean-field description of the transition.
Scaling (from the literature, not computed here): alpha_q = 2 - nu (d + z); 1+1D Ising nu = 1 -> marginal;
tricritical Ising nu = 5/9 -> alpha_q = 8/9 > 0, coupling relevant (eigenvalue 2 - 2 Delta_eps = 8/5).
Usage: python phonon_coupling.py -> phonon_coupling.out"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows'))
import numpy as np, mpmath as mp
from landau import coeffs, e_cell
from symbol import derivs_at_pi, derivs_at_0
from exact_symbol import lstar
mp.mp.dps = 30

def analyse(kappa):
    ls = lstar(kappa); s = np.sqrt(1 - kappa**2); lp, Y0 = mp.mpf(ls*s), mp.mpf(-kappa/s)
    a2, a4 = coeffs(lp, Y0)
    g = mp.matrix([mp.diff(lambda x: coeffs(x, Y0)[0], lp), mp.diff(lambda y: coeffs(lp, y)[0], Y0)])
    e0 = lambda x, y: e_cell(y, y, x)
    M = mp.matrix([[mp.diff(e0, (lp, Y0), (2, 0)), mp.diff(e0, (lp, Y0), (1, 1))],
                   [mp.diff(e0, (lp, Y0), (1, 1)), mp.diff(e0, (lp, Y0), (0, 2))]])
    shift = (g.T*mp.inverse(M)*g)[0]/2
    Bpp, Cpi, dpi = derivs_at_pi(kappa, float(lp)); v = np.sqrt(abs(Cpi)*Bpp/2)
    B0, c2, d0 = derivs_at_0(kappa, float(lp)); c0 = np.sqrt(B0*abs(c2))
    return dict(ls=ls, a4=float(a4), g=(float(g[0]), float(g[1])), shift=float(shift), v=v, dpi=dpi, c0=c0, d0=d0)

if __name__ == '__main__':
    out = []; P = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
    P(' kappa     l*      a4_fixed    LP shift   shift/|a4|   g_strain    g_impulse   v (crit)  c0 (sound)  c0/v   slowest sound / fastest crit')
    for kappa in [0.9203, 0.94, 0.9602, 0.98, 0.995]:
        r = analyse(kappa)
        slow_c = r['c0'] - abs(r['d0']); fast_v = r['v'] + abs(r['dpi'])
        P(f' {kappa:6.4f}  {r["ls"]:7.4f}   {r["a4"]:10.3e}  {r["shift"]:10.3e}   {r["shift"]/abs(r["a4"]) if r["a4"] else float("inf"):8.3f}'
          f'   {r["g"][0]:10.3e}  {r["g"][1]:10.3e}   {r["v"]:7.4f}   {r["c0"]:8.4f}   {r["c0"]/r["v"]:5.2f}   {slow_c/fast_v:6.2f}')
    open('phonon_coupling.out', 'w').write('\n'.join(out) + '\n')
