"""Acoustic spacetime of a ring of quantized vortices in a uniform 2D condensate.
Units: lengths in l_e = kappa/(2 pi c) (radius of a single vortex's ergoregion, ~ sqrt(2) healing lengths), speeds in c.
A unit vortex at w gives v = i (z - w)/|z - w|^2 (speed 1/r); a central vortex of charge q gives q/r.
Acoustic metric ds^2 = -(c^2 - v^2) dt^2 - 2 v.dx dt + dx^2.
(1) Lab-frame ergoregion: |v| > 1.
(2) The ring rotates rigidly at Omega = (q + (N-1)/2)/R^2, so the metric is stationary only along the helical field
    chi = d_t + Omega d_phi; g(chi, chi) > 0 (co-rotating ergoregion) where |v - Omega z^ x x| > 1.
Closed forms (azimuthal average, verified on the grid):
    acoustic light cylinder  r_LC = [1 + sqrt(1 + 4 Omega Q)] / (2 Omega),  Q = q + N;
    inner co-rotating ergoregion (q > 0)  r_in = [-1 + sqrt(1 + 4 Omega q)] / (2 Omega).
Validity: uniform-density model; needs points of interest >~ 1-2 l_e from vortex centres.
Usage: python ergo.py   (prints the radius below which the ring's ergoregions merge into one annulus)"""
import numpy as np
from scipy import ndimage

def velocity(X, Y, N, R, q=0.0):
    Z = X + 1j*Y
    V = np.zeros_like(Z)
    for k in range(N):
        w = R*np.exp(2j*np.pi*k/N)
        d = Z - w
        V += 1j*d/np.abs(d)**2
    if q:
        V += 1j*q*Z/np.abs(Z)**2
    return V

def omega_ring(N, R, q=0.0):
    return (q + (N-1)/2)/R**2

def analyse(N, R, q=0.0, frame='lab', L=None, n=1200):
    """Returns (number of ergoregion components, whether one component encircles or covers the centre, mask, axis)."""
    L = L or 1.6*max(R, N, q+1)
    x = np.linspace(-L, L, n); X, Y = np.meshgrid(x, x)
    V = velocity(X, Y, N, R, q)
    if frame == 'corot':
        V = V - 1j*omega_ring(N, R, q)*(X+1j*Y)
    E = np.abs(V) > 1
    lab, nc = ndimage.label(E)
    c = n//2
    encl = False
    for k in range(1, nc+1):
        holes = ndimage.label(~(lab == k))[0]
        if holes[c, c] != holes[0, 0]:
            encl = True; break
    return nc, encl, E, x

def merge_radius(N, q=0.0):
    """Largest ring radius R (in l_e) at which one lab-frame ergoregion component encircles the centre."""
    lo, hi = 0.3, 4.0*max(N, q+1)
    if not analyse(N, lo, q)[1]: return None
    for _ in range(30):
        m = (lo+hi)/2
        if analyse(N, m, q, n=800)[1]: lo = m
        else: hi = m
    return lo

if __name__ == '__main__':
    print('N  R_annulus (l_e)  R_annulus/N   neighbour spacing at R_annulus (l_e)')
    for N in [3, 4, 5, 6, 7, 8, 10, 12, 16, 19, 24]:
        Rm = merge_radius(N)
        sp = 2*Rm*np.sin(np.pi/N)
        print(f'{N:2d}  {Rm:8.3f}        {Rm/N:.3f}         {sp:.3f}')
