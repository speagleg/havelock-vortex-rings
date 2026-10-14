"""Planar point vortices with an isotropic pair kernel W(r): relative-equilibrium check, linear stability.
H = sum_{i<j} G_i G_j W(r_ij);  G_k xdot_k = dF/dy_k, G_k ydot_k = -dF/dx_k,  F = H - Omega * sum G_k |z_k|^2/2.
Kernels: 'log' (W = -ln r, plain 2D Euler) and 'k0' (W = K0(r/lam), quasi-geostrophic screening)."""
import numpy as np
from scipy.special import k0, k1

def kernel(kind, lam=1.0):
    if kind == 'log':
        return (lambda r: -1.0/r, lambda r: 1.0/r**2)          # W', W''  for W = -ln r
    if kind == 'k0':
        # W = K0(r/lam): W' = -K1/lam, W'' = (K0 + K1/(r/lam))/lam^2
        return (lambda r: -k1(r/lam)/lam, lambda r: (k0(r/lam) + k1(r/lam)/(r/lam))/lam**2)
    raise ValueError(kind)

def grad_hess(P, G, kind, lam=1.0):
    Wp, Wpp = kernel(kind, lam)
    n = len(G); g = np.zeros(2*n); Hm = np.zeros((2*n, 2*n))
    for i in range(n):
        for j in range(i+1, n):
            d = P[i]-P[j]; r = np.hypot(*d); u = d/r; c = G[i]*G[j]
            gi = c*Wp(r)*u
            g[2*i:2*i+2] += gi; g[2*j:2*j+2] -= gi
            B = c*(Wpp(r)*np.outer(u, u) + Wp(r)/r*(np.eye(2)-np.outer(u, u)))
            Hm[2*i:2*i+2, 2*i:2*i+2] += B; Hm[2*j:2*j+2, 2*j:2*j+2] += B
            Hm[2*i:2*i+2, 2*j:2*j+2] -= B; Hm[2*j:2*j+2, 2*i:2*i+2] -= B
    return g, Hm

def analyse(P, G, kind, lam=1.0):
    """Returns (Omega, relative-equilibrium residual, max Re(eigenvalue)/max|eigenvalue|, eigenvalues)."""
    P = np.asarray(P, float); G = np.asarray(G, float); n = len(G)
    g, Hm = grad_hess(P, G, kind, lam)
    gL = np.concatenate([G[k]*P[k] for k in range(n)])
    Om = g @ gL / (gL @ gL)
    res = np.max(np.abs(g - Om*gL)) / max(1e-300, np.max(np.abs(g)))
    HF = Hm - Om*np.diag(np.repeat(G, 2))
    Jinv = np.zeros((2*n, 2*n))
    for k in range(n):
        Jinv[2*k, 2*k+1] = 1/G[k]; Jinv[2*k+1, 2*k] = -1/G[k]
    ev = np.linalg.eigvals(Jinv @ HF)
    scale = max(1e-300, np.max(np.abs(ev)))
    return Om, res, np.max(ev.real)/scale, ev

def ring(N, R, q=None):
    """N unit vortices on a circle of radius R, plus an optional central vortex of strength q (q != 0)."""
    P = [(R*np.cos(2*np.pi*k/N), R*np.sin(2*np.pi*k/N)) for k in range(N)]
    G = [1.0]*N
    if q is not None:
        P.append((0.0, 0.0)); G.append(q)
    return P, G
