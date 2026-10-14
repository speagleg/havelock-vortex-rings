"""Exact pairing-mode (q = pi) criterion vs the full numerical threshold, and the critical Bloch mode.
The row's threshold is the zero of one explicit lattice sum, A_UP(pi)(alpha, l') = 0 (normal stiffness of the
pairing mode). Note: AUP_pi truncates at |j| <= 40, which is too short near kappa -> 1 (interactions become
power-law along the horocycle); trust it for kappa <= 0.98.
Usage: python pairing_root.py   (~4 min)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import brentq
import rows

def AUP_pi(alpha, lp, J=None, h=1e-6):
    """d(u_dot_0)/d(phi) for the alternating perturbation phi_j = alpha + h(-1)^j, exact kernel, |j| <= J."""
    J = J or max(2, min(40, int(250/lp)))
    def ud0(phis):
        u = lp*np.arange(-J, J+1); i0 = J
        s0 = np.sin(phis[i0]); tot = 0.0
        for j in range(2*J+1):
            if j == i0: continue
            du = u[i0]-u[j]; dp = phis[i0]-phis[j]; sj = np.sin(phis[j])
            C = 1 + (np.cosh(du) - np.cos(dp))/(s0*sj)
            dCdp = np.sin(dp)/(s0*sj) - (C-1)*np.cos(phis[i0])/s0
            tot += -1/(C**2-1)*dCdp
        return s0**2*tot
    sgn = np.array([(-1)**abs(j) for j in range(-J, J+1)], float)
    base = np.full(2*J+1, alpha)
    return (ud0(base + h*sgn) - ud0(base - h*sgn))/(2*h)

def mode_spectrum(M, l, kappa):
    """Largest growth rate and its Bloch wavenumber theta = j/M (eigenvector FFT) from the supercell Jacobian."""
    alpha = np.arccos(kappa); lp = l*np.sin(alpha); L = M*lp; nimg = max(1, int(60/L)+1)
    u = np.arange(M)*lp; phi = np.full(M, alpha)
    ud, pd = rows.vel_hyper(u, phi, L, nimg); s = ud[0]
    def F(q):
        a, b = rows.vel_hyper(q[:M], q[M:], L, nimg); return np.concatenate([a - s, b])
    q0 = np.concatenate([u, phi]); Jm = np.zeros((2*M, 2*M)); h = 1e-7
    for c in range(2*M):
        qp = q0.copy(); qm = q0.copy(); qp[c] += h; qm[c] -= h
        Jm[:, c] = (F(qp)-F(qm))/(2*h)
    ev, V = np.linalg.eig(Jm)
    k = np.argmax(ev.real); v = V[:, k]
    m = int(np.argmax(np.abs(np.fft.fft(v[M:])) + np.abs(np.fft.fft(v[:M]))))
    return ev[k].real, min(m, M-m)/M

if __name__ == '__main__':
    import warnings; warnings.filterwarnings('ignore')
    data = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rows_boundary.txt'))
    print('kappa   l* (full)   l* (pairing q=1/2 root)   critical theta just below l*   asymptotic -ln(cos 2a)/sin a')
    for k, lfull in data:
        a = np.arccos(k)
        try:
            lpair = brentq(lambda lp: AUP_pi(a, lp), 1e-3, 25)/np.sin(a)
        except ValueError:
            lpair = np.nan
        g, th = mode_spectrum(40, 0.97*lfull, k)
        asym = -np.log(np.cos(2*a))/np.sin(a)
        print(f'{k:6.3f}  {lfull:9.4f}   {lpair:9.4f}                  {th:.3f}                        {asym:.4f}', flush=True)
