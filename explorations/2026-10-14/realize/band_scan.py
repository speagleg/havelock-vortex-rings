"""Finite film band (band.py): how the free edges change the H^2 row picture.
Band = parallels at hyperbolic distance |h| <= Hb from the waist geodesic (kappa = tanh h on a parallel).
(1) Edge self-term at large spacing: B_self(alpha) = d/dphi[sin^2(phi) S'(phi)] at phi = alpha; with C < 0 the row is
    stable at large spacing iff B_self >= 0 (pair couplings decay like e^{-l'}, the self-term does not).
(2) Direct supercell stability on a grid of (kappa, arc spacing) for several band widths.
Usage: python band_scan.py"""
import numpy as np, sys
sys.path.insert(0, '.')
from band import growth
def edges(Hb):
    p1 = np.arccos(np.tanh(Hb)); return p1, np.pi - p1
def Bself(alpha, p1, p2, eps=1e-6):
    w = p2 - p1; k = np.pi/w
    Sp = lambda p: 0.5*k/np.tan((p - p1)*k) - 0.5/np.tan(p)
    f = lambda p: np.sin(p)**2*Sp(p)
    return (f(alpha + eps) - f(alpha - eps))/(2*eps)
if __name__ == '__main__':
    print('(1) edge self-term stiffness B_self at the parallel of curvature kappa (sign decides large-spacing stability)')
    for Hb in [1.0, 1.5, 2.0, 3.0, 5.0]:
        p1, p2 = edges(Hb)
        row = []
        for kap in [0.3, 0.6, 0.7071, 0.8, 0.9]:
            h = np.arctanh(kap)
            if h >= Hb: row.append('   (off)'); continue
            row.append(f'{Bself(np.arccos(kap), p1, p2):+.2e}')
        print(f'   band half-width Hb = {Hb}: ' + '  '.join(row) + '   (kappa = 0.3, 0.6, 0.7071, 0.8, 0.9)')
    print('(2) supercell stability (S = stable) for arc spacings 1, 2, 3, 5, 8, 12')
    for Hb in [2.0, 3.0, 5.0]:
        p1, p2 = edges(Hb)
        print(f'   Hb = {Hb}')
        for kap in [0.5, 0.7, 0.75, 0.8, 0.9]:
            if np.arctanh(kap) >= Hb: continue
            al = np.arccos(kap); s = np.sin(al)
            print(f'      kappa {kap}: ' + ' '.join('S' if growth(al, l*s, p1, p2) < 1e-4 else 'u' for l in [1, 2, 3, 5, 8, 12]), flush=True)
