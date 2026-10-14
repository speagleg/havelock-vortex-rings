"""Strip roll-up vs ring survival.  Usage: python h2strip.py h2   |   python h2strip.py plane
A uniform-vorticity strip (outer radius R, thickness D, in curvature radii) is assumed to roll up into its
fastest-growing m* vortices; the resulting ring is judged by Havelock (plane: m* <= 7) or by the tanh-model
threshold (H^2: cosh r > kappa_inf N + C1/N, the note's large-N formula, adequate as a screen).
Slow at R >= 3 (dense eigenproblems); see h2check.py for targeted m-ranges there."""
import sys; sys.path.insert(0, '.')
from kuo import *
import numpy as np
kinf = 0.1343122142726; C1e = 0.5412188; C1o = 0.4010382

def ring_stable(N, r):
    C1 = C1e if N % 2 == 0 else C1o
    return np.cosh(r) > kinf*N + C1/N

if __name__ == '__main__':
    geom = sys.argv[1]
    print(f'{geom}: strip of uniform vorticity, outer radius R, thickness D (curvature radii); fastest m*, roll-up spacing, survival of an m*-ring')
    for R in [1.0, 2.0, 3.0, 4.0]:
        for D in [0.05, 0.1, 0.2, 0.4]:
            if D >= R: continue
            n = int(1600*(R+3)/3); rmax = R+3.0
            gs = {}
            for m in range(1, 400):
                g = growth(geom, rmax, n, tophat(R-D, R), m)
                gs[m] = g
                if m > 6 and g <= 1e-9 and max(gs.values()) > 0 and m > 2*max(gs, key=gs.get): break
            mstar = max(gs, key=gs.get)
            rmid = R-D/2
            circ = 2*np.pi*(np.sinh(rmid) if geom == 'h2' else rmid)
            sp = circ/mstar
            if geom == 'h2': verdict = 'stable ring' if ring_stable(mstar, rmid) else 'ring unstable -> merge'
            else: verdict = 'stable ring' if mstar <= 7 else 'ring unstable -> merge'
            print(f' R={R} D={D}: m*={mstar:3d} growth={gs[mstar]:.4f} spacing={sp:.3f} -> {verdict}', flush=True)
