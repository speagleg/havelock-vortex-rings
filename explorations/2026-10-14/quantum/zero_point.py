"""Zero-point fluctuations of a quantised vortex row on a hypercycle (harmonic theory), as Lindemann ratios.
Ground-state covariances per Bloch mode (independent of delta, which only shifts energies):
    <|Y_q|^2> = (eps/2) sqrt(h_uu/h_YY) = (eps/(2 s^2)) sqrt(|C|/B),   <|u_q|^2> = (eps s^2/2) sqrt(B/|C|).
Proper displacements (metric (du^2 + dphi^2)/sin^2 phi, curvature radius R): normal d rho = s dY, along-curve ds = du/s.
    <d rho^2>/R^2          = (eps/(4 pi)) int_{-pi}^{pi} sqrt(|C|/B) dq            (transverse, per vortex)
    <(ds_{j+1}-ds_j)^2>/R^2 = (eps/(2 pi)) int_{-pi}^{pi} (1 - cos q) sqrt(B/|C|) dq   (relative, neighbours)
For superfluid vortices eps = 1/(2 pi n R^2) (Haldane & Wu 1985: one state per particle of the fluid), so with
N_l = n (l R)^2 (fluid particles per spacing squared) the Lindemann ratios are gamma^2 = Lambda / N_l with
    Lambda_perp = (1/(8 pi^2 l^2)) int sqrt(|C|/B) dq,   Lambda_par = (1/(4 pi^2 l^2)) int (1 - cos q) sqrt(B/|C|) dq.
Near the threshold B(pi) -> 0 and Lambda_perp diverges logarithmically (the harmonic theory then fails).
Usage: python zero_point.py -> zero_point.out"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows'))
import numpy as np
from symbol import symbol, derivs_at_pi
from exact_symbol import lstar

def Lambdas(kappa, l, nq=40000):
    s = np.sqrt(1 - kappa**2); lp = l*s
    h = 2*np.pi/nq; q = -np.pi + h*(np.arange(nq) + 0.5)     # midpoint rule, never hits q = 0
    B, C, dl = symbol(kappa, lp, q)
    assert (B > 0).all()
    w = np.full(nq, h)
    Lp = np.sum(np.sqrt(np.abs(C)/B)*w)/(8*np.pi**2*l*l)
    La = np.sum((1 - np.cos(q))*np.sqrt(B/np.abs(C))*w)/(4*np.pi**2*l*l)
    return Lp, La

if __name__ == '__main__':
    out = []; P = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
    P('Lindemann ratios gamma^2 = Lambda / N_l  (N_l = fluid particles per spacing^2):')
    P('  kappa   l/l*    Lambda_perp   Lambda_par    N_l needed for gamma_perp < 0.15')
    for kappa in [0.75, 0.8, 0.9, 0.98]:
        ls = lstar(kappa)
        for r in [1.001, 1.01, 1.05, 1.2, 1.5, 2.0, 3.0]:
            Lp, La = Lambdas(kappa, ls*r)
            P(f'  {kappa:5.2f}  {r:6.3f}   {Lp:11.4e}   {La:11.4e}   {Lp/0.15**2:10.3g}')
    # logarithmic divergence: Lambda_perp vs ln(1/(l/l* - 1))
    P('\nNear threshold (kappa = 0.9): Lambda_perp against ln(1/(l/l* - 1)) -> slope (log divergence)')
    ls = lstar(0.9); xs, ys = [], []
    for r in [1 + 10**(-k) for k in [2, 3, 4, 5]]:
        Lp, _ = Lambdas(0.9, ls*r, nq=400000); xs.append(np.log(1/(r - 1))); ys.append(Lp)
        P(f'  l/l* - 1 = {r - 1:.0e}: Lambda_perp = {Lp:.5f}')
    P(f'  successive slopes d Lambda / d ln: {np.round(np.diff(ys)/np.diff(xs), 5)}')
    s = np.sqrt(1 - 0.81); Bpp, Cpi, _ = derivs_at_pi(0.9, ls*s)
    P(f'  predicted slope (harmonic, B(pi) linear in l - l*): sqrt(2|C(pi)|/B\'\'(pi))/(8 pi^2 l*^2) = '
      f'{np.sqrt(2*abs(Cpi)/Bpp)/(8*np.pi**2*ls**2):.5f}')
    open('zero_point.out', 'w').write('\n'.join(out) + '\n')
