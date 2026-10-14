"""Negative result kept for the record: point-vortex 'jets' (a ring + an opposite-signed ring) on the sphere.
Every m up to N/2 is unstable and the grid scale m = N/2 grows fastest (vortex-sheet Kelvin–Helmholtz),
so point-vortex sheets cannot select a polygon number; finite-width profiles (kuo.py) are needed.
Note: the two rings co-rotate only if the second ring's strength is tuned (see README)."""
import numpy as np

def vel(X, G):
    D = 1 - X @ X.T; np.fill_diagonal(D, np.inf)
    S = (G[None, :]/D) @ X / (4*np.pi)
    return np.cross(X, S)

def rings(spec):
    """spec: list of (N, colatitude, Gamma, phase)."""
    X = []; G = []
    for N, th, Gm, ph0 in spec:
        for k in range(N):
            ph = 2*np.pi*k/N + ph0
            X.append([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)]); G.append(Gm)
    return np.array(X), np.array(G, float)

def jacobian(X, G, h=1e-6):
    n = len(G); z = np.array([0, 0, 1.0])
    V0 = vel(X, G)
    w = V0[0] @ np.cross(z, X[0]) / np.linalg.norm(np.cross(z, X[0]))**2
    Vr = lambda Y: vel(Y, G) - w*np.cross(z, Y)
    th = np.arccos(np.clip(X[:, 2], -1, 1)); ph = np.arctan2(X[:, 1], X[:, 0])
    eth = np.stack([np.cos(th)*np.cos(ph), np.cos(th)*np.sin(ph), -np.sin(th)], 1)
    eph = np.stack([-np.sin(ph), np.cos(ph), np.zeros(n)], 1)
    B = [eth, eph]; M = np.zeros((2*n, 2*n))
    for k in range(n):
        for a in range(2):
            Yp = X.copy(); Ym = X.copy()
            Yp[k] = (X[k]+h*B[a][k]); Yp[k] /= np.linalg.norm(Yp[k])
            Ym[k] = (X[k]-h*B[a][k]); Ym[k] /= np.linalg.norm(Ym[k])
            dV = (Vr(Yp)-Vr(Ym))/(2*h)
            M[0::2, 2*k+a] = np.einsum('ij,ij->i', dV, eth)
            M[1::2, 2*k+a] = np.einsum('ij,ij->i', dV, eph)
    return M, w, np.max(np.abs(Vr(X)))

def mode_growth(spec):
    """Growth rate and frequency per azimuthal wavenumber m (all rings must share N)."""
    X, G = rings(spec); M, w, res = jacobian(X, G)
    ev, V = np.linalg.eig(M)
    N = spec[0][0]; out = {}
    for lam, v in zip(ev, V.T):
        a = v[0:2*N:2]
        if np.linalg.norm(a) < 1e-8: a = v[1:2*N:2]
        m = int(np.argmax(np.abs(np.fft.fft(a)))); m = min(m, N-m)
        if m not in out or lam.real > out[m][0]:
            out[m] = (lam.real, abs(lam.imag))
    return out, w, res

def corotating_strength(N, t1, t2):
    """Strength of the second ring (at colatitude t2) that makes both rings co-rotate with a unit ring at t1."""
    z = np.array([0, 0, 1.0])
    def domega(g):
        X, G = rings([(N, t1, 1.0, 0), (N, t2, g, 0)]); V = vel(X, G)
        o = [V[i] @ np.cross(z, X[i]) / np.linalg.norm(np.cross(z, X[i]))**2 for i in (0, N)]
        return o[0]-o[1]
    a, b = -1.5, -0.5
    for _ in range(80):
        m = (a+b)/2
        if np.sign(domega(m)) == np.sign(domega(a)): a = m
        else: b = m
    return (a+b)/2

if __name__ == '__main__':
    N = 80
    for th0 in [15, 30, 45]:
        for d in [2, 6]:
            t1, t2 = np.radians(th0-d/2), np.radians(th0+d/2)
            g2 = corotating_strength(N, t1, t2)
            out, w, res = mode_growth([(N, t1, 1.0, 0), (N, t2, g2, 0)])
            mstar = max(out, key=lambda mm: out[mm][0])
            print(f'theta0={th0} d={d}: Gamma2={g2:.4f} residual={res:.0e} fastest m*={mstar} (N/2={N//2})')
