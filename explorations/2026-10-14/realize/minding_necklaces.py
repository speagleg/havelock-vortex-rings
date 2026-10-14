"""Necklaces of N vortices on an EMBEDDED constant-curvature (K = -1) surface of revolution of hyperboloid type (Minding):
profile radius r(h) = a cosh h around the waist geodesic (length 2 pi a), parallels = hypercycles at distance h, kappa = tanh h.
In R^3 the surface exists only while |dr/ds| = a sinh|h| <= 1, so the film band is |h| <= He = arsinh(1/a).
A necklace of N on the parallel h has coordinate spacing l' = 2 pi a/N (arc spacing l = 2 pi a cosh h/N).
Film edges impermeable (band.py: strip Green's function with images in u of period 2 pi a; one choice of the net
circulation around the waist). Question: does ANY stable necklace (N >= 2) exist on such a surface?
Usage: python minding_necklaces.py"""
import numpy as np, sys
sys.path.insert(0, '.')
from band import growth
res = []
for a in [0.02, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5]:
    He = np.arcsinh(1/a); p1 = np.arccos(np.tanh(He)); p2 = np.pi - p1
    best = None
    for N in range(2, 13):
        lp = 2*np.pi*a/N
        for h in np.linspace(0.05, He*0.995, 80):
            al = np.arccos(np.tanh(h))
            g = growth(al, lp, p1, p2, M=N)
            if g < 1e-6:
                best = (N, h, np.tanh(h), lp*np.cosh(h)); break
        if best: break
    print(f'waist radius a = {a:4}: edge He = {He:.3f}, max kappa on surface = {np.tanh(He):.3f}; first stable necklace: {best}', flush=True)
