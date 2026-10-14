"""Vortex necklaces around a closed geodesic: the collar model of the Bolza systole (tanh model).
Near a closed geodesic of length L a hyperbolic surface is isometric to the hyperbolic cylinder H^2/<translation by L>.
N equal vortices equally spaced on the hypercycle at distance h from the geodesic form a relative equilibrium
(rotating around the collar); this is the hypercycle row of rows.py with period L and N vortices per period:
kappa = tanh h, spacing l = L cosh h / N. Bolza systole: L = 2 arccosh(1 + sqrt 2); collar half-width
w = arsinh(1/sinh(L/2)) (collar lemma), so inside the collar kappa <= tanh w = sqrt 2 - 1 (exact for Bolza).
Caveat: the cylinder keeps only the translation subgroup of the Bolza group; the rest of the surface is not modelled
(and on the closed surface the net circulation must vanish).
Usage: python collar.py > collar.out   (~6 min)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import rows

L_SYS = 2*np.arccosh(1+np.sqrt(2))
W_COLLAR = np.arcsinh(1/np.sinh(L_SYS/2))

def necklace(N, h, L=L_SYS):
    kappa = np.tanh(h); l = L*np.cosh(h)/N
    if h == 0:
        kappa = 1e-12
    g, res, s = rows.stab_hyper(N, l, kappa)
    return g, l, kappa

def h_threshold(N, L=L_SYS, hmax=6.0, tol=1e-4):
    """Smallest h at which the necklace becomes stable (scan then bisect); None if never stable up to hmax."""
    hs = np.linspace(0.0, hmax, 241)
    st = [necklace(N, h, L)[0] < tol for h in hs]
    if st[0]: return 0.0
    for i in range(1, len(hs)):
        if st[i]:
            lo, hi = hs[i-1], hs[i]
            for _ in range(40):
                m = (lo+hi)/2
                if necklace(N, m, L)[0] < tol: hi = m
                else: lo = m
            return hi
    return None

if __name__ == '__main__':
    print(f'Bolza systole L = {L_SYS:.6f}, collar half-width w = {W_COLLAR:.6f}, max kappa in collar = {np.tanh(W_COLLAR):.4f}')
    print('N   stable inside the collar?   first stable distance h* (kappa*, spacing*)')
    for N in range(2, 13):
        inside = [necklace(N, h)[0] < 1e-4 for h in np.linspace(0, W_COLLAR, 12)]
        hs = h_threshold(N)
        txt = 'never (h <= 6)' if hs is None else f'h* = {hs:.4f}  (kappa* = {np.tanh(hs):.4f}, spacing {L_SYS*np.cosh(hs)/N:.4f})'
        print(f'{N:2d}  {"yes" if all(inside) else ("partly" if any(inside) else "no"):>6}                    {txt}', flush=True)
