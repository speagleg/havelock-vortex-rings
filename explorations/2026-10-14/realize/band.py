"""A superfluid film on a FINITE band of a constant-negative-curvature surface of revolution (hyperboloid type, Minding).
Model: the band between the hypercycles phi1 < phi < phi2 of H^2 (strip coordinates z = u + i phi, metric |dz|^2/sin^2 phi),
with impermeable film edges (stream function constant on both edges). The parallels of the surface are the hypercycles
phi = const; the waist geodesic is phi = pi/2. Point-vortex dynamics on a curved film (Kirchhoff-Routh):
   u_dot = sin^2(phi) dH/dphi,  phi_dot = -sin^2(phi) dH/du,
   H = sum_{i<j} W_band(z_i, z_j) + sum_i S(phi_i),
   W_band = flat Green's function (x 2 pi) of the strip of width w = phi2 - phi1 (conformal invariance),
   S(phi) = (1/2) ln sin(pi (phi - phi1)/w) - (1/2) ln sin(phi)   (flat Robin term of the strip + geometric self-energy).
For phi1 -> 0, phi2 -> pi this is the complete H^2 (tanh model): S = const, W_band = -ln tanh(d/2).
Usage: imported by band_scan.py"""
import numpy as np

def grad(u, phi, L, nimg, p1, w):
    k = np.pi/w
    du = u[:, None, None] - (u[None, :, None] + L*np.arange(-nimg, nimg+1)[None, None, :])
    yi, yj = (phi[:, None, None] - p1)*k, (phi[None, :, None] - p1)*k
    M = len(u); mask = np.ones(du.shape, bool); mask[np.arange(M), np.arange(M), nimg] = False
    ok = mask & (np.abs(du*k) < 300); X = np.where(ok, du*k, 0.0)
    # W = ln|(e^{Z} - e^{conj Zj})/(e^{Z} - e^{Zj})| = (1/2) ln[(cosh X - cos(yi + yj))/(cosh X - cos(yi - yj))]
    A = np.cosh(X) - np.cos(yi + yj); B = np.cosh(X) - np.cos(yi - yj)
    B = np.where(ok, B, 1.0); A = np.where(ok, A, 1.0)
    dWdX = 0.5*np.sinh(X)*(1/A - 1/B)
    dWdy = 0.5*(np.sin(yi + yj)/A - np.sin(yi - yj)/B)
    Hu = (np.where(ok, dWdX, 0)*k).sum(axis=(1, 2)); Hp = (np.where(ok, dWdy, 0)*k).sum(axis=(1, 2))
    Hp = Hp + 0.5*k/np.tan((phi - p1)*k) - 0.5/np.tan(phi)          # dS/dphi
    return Hu, Hp
def vel(u, phi, L, nimg, p1, w):
    Hu, Hp = grad(u, phi, L, nimg, p1, w); s2 = np.sin(phi)**2
    return s2*Hp, -s2*Hu
def growth(alpha, lp, p1, p2, M=24, h=1e-7):
    w = p2 - p1; L = M*lp; nimg = max(1, int(60*w/(np.pi*L)) + 1)
    u = np.arange(M)*lp; phi = np.full(M, alpha)
    ud, pd = vel(u, phi, L, nimg, p1, w); s = ud[0]
    F = lambda q: (lambda v: np.concatenate([v[0] - s, v[1]]))(vel(q[:M], q[M:], L, nimg, p1, w))
    q0 = np.concatenate([u, phi]); J = np.zeros((2*M, 2*M))
    for a in range(2*M):
        qp = q0.copy(); qm = q0.copy(); qp[a] += h; qm[a] -= h; J[:, a] = (F(qp) - F(qm))/(2*h)
    ev = np.linalg.eigvals(J); return np.max(ev.real)/np.max(np.abs(ev))
