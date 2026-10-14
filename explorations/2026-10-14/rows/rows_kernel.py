"""Kernel dependence of the hypercycle-row boundary. Screened vortex models on H^2: W = Q_nu(cosh d), the Green's
function of (Laplacian - nu(nu+1)) on H^2 (Legendre Q). nu = 0 is the tanh model (Q_0 = artanh(1/C) = -ln tanh(d/2));
nu = 1: Q_1(C) = C artanh(1/C) - 1, decaying like exp(-2 d).
Large-spacing prediction (asymptotic_check.py with W ~ exp(-gamma d)): kappa_c = 1/sqrt 2 for every gamma.
Usage: python rows_kernel.py   (~6 min)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import rows

def dQ1(C):
    """d/dC [C artanh(1/C) - 1] = artanh(1/C) - C/(C^2 - 1), with a series for large C (cancellation)."""
    C = np.asarray(C, float); out = np.empty_like(C)
    big = C > 50
    x = 1/C[big]
    out[big] = -sum(2*k/(2*k+1)*x**(2*k+1) for k in range(1, 8))
    Cs = C[~big]
    out[~big] = np.arctanh(1/Cs) - Cs/(Cs**2-1)
    return out

def vel_hyper_Q1(u, phi, L, nimg):
    du = u[:, None, None] - (u[None, :, None] + L*np.arange(-nimg, nimg+1)[None, None, :])
    dp = phi[:, None, None] - phi[None, :, None]
    si = np.sin(phi)[:, None, None]; sj = np.sin(phi)[None, :, None]
    ok = np.abs(du) < 300
    dus = np.where(ok, du, 0.0)
    C = 1 + (np.cosh(dus) - np.cos(dp))/(si*sj)
    M = len(u); ok[np.arange(M), np.arange(M), nimg] = False
    Cs = np.where(ok, C, 2.0)
    Wp = np.where(ok, dQ1(Cs), 0.0)
    dCdu = np.sinh(dus)/(si*sj)
    dCdp = np.sin(dp)/(si*sj) - (Cs-1)*np.cos(phi)[:, None, None]/si
    Hu = (Wp*dCdu).sum(axis=(1, 2)); Hp = (Wp*dCdp).sum(axis=(1, 2))
    s2 = np.sin(phi)**2
    return s2*Hp, -s2*Hu

if __name__ == '__main__':
    rows.vel_hyper = vel_hyper_Q1          # stab_hyper looks the velocity up at call time
    for L in [30, 60]:
        lo, hi = 0.6, 0.8
        for _ in range(30):
            m = (lo+hi)/2
            if rows.stab_hyper(40, L, m)[0] < 1e-3: hi = m
            else: lo = m
        print(f'Q_1 kernel, spacing {L}: kappa_c = {hi:.7f}  (1/sqrt 2 = {1/np.sqrt(2):.7f})', flush=True)
    for k in [0.75, 0.8, 0.9, 0.95, 0.99]:
        print(f'Q_1 kernel: kappa = {k}: l* = {rows.critical(lambda l: rows.stab_hyper(40, l, k), 0.1, 60, tol=1e-3):.4f}', flush=True)
