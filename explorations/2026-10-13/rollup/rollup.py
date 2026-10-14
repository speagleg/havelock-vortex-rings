"""Nonlinear roll-up of a thin vorticity strip: H^2 (Poincare disk, tanh kernel) vs the flat plane.
Vortex-blob discretisation. Strip between geodesic radii R-D and R, uniform vorticity 1.
H = sum_{i<j} G_i G_j W_ij,  W = -ln|z-w| (+ ln|1 - z conj(w)| on H^2), first term regularised with a blob.
Equations: xdot = lam^{-1} dH/dy / G, ydot = -lam^{-1} dH/dx / G, lam = 4/(1-|z|^2)^2 on H^2, 1 in the plane.
Usage: python rollup.py GEOM R D T_END [K M SEED]      e.g.  python rollup.py h2 2.0 0.2 100 3 480
Saves roll_<geom>_R<R>_D<D>_s<seed>.npz (positions at 121 times) and prints a quick azimuthal diagnostic;
use clumps.py for the clump count. Runtime: ~10 min (plane, R=1) to ~40 min (H^2, R=2) on one core."""
import sys, numpy as np
from scipy.integrate import solve_ivp

def setup(geom, R, D, K, M, seed):
    rng = np.random.default_rng(seed)
    edges = np.linspace(R-D, R, K+1)
    z = []; G = []
    for i in range(K):
        r1, r2 = edges[i], edges[i+1]
        if geom == 'h2':
            area = 2*np.pi*(np.cosh(r2)-np.cosh(r1)); rm = np.arccosh((np.cosh(r1)+np.cosh(r2))/2)
            rho = np.tanh(rm/2)
        else:
            area = np.pi*(r2*r2-r1*r1); rm = np.sqrt((r1*r1+r2*r2)/2); rho = rm
        for k in range(M):
            ph = 2*np.pi*(k+0.5*(i % 2))/M + 1e-3*rng.standard_normal()
            z.append(rho*np.exp(1j*ph)); G.append(area/M)
    return np.array(z), np.array(G)

def make_rhs(geom, G, delta):
    def rhs(t, y):
        n = len(G); z = y[:n] + 1j*y[n:]
        d = z[:, None]-z[None, :]
        if geom == 'h2':
            lam = 4/(1-np.abs(z)**2)**2
            dE = delta*(1-np.abs(z)**2)/2           # geodesic blob size -> Euclidean
            r2 = np.abs(d)**2 + dE[:, None]**2
            np.fill_diagonal(r2, np.inf)
            img = -z[None, :]/(1-np.conj(z)[:, None]*z[None, :])
            np.fill_diagonal(img, 0)
            Gk = (G[None, :]*(-d/r2 + img)).sum(1)
        else:
            lam = 1.0
            r2 = np.abs(d)**2 + delta**2; np.fill_diagonal(r2, np.inf)
            Gk = (G[None, :]*(-d/r2)).sum(1)
        v = -1j*Gk/lam
        return np.concatenate([v.real, v.imag])
    return rhs

def clumps(z, G, nbins=720, smooth_deg=3.0):
    ang = np.angle(z) % (2*np.pi)
    h, _ = np.histogram(ang, bins=nbins, range=(0, 2*np.pi), weights=G)
    s = smooth_deg/360*nbins
    kx = np.arange(-int(4*s), int(4*s)+1); ker = np.exp(-0.5*(kx/s)**2); ker /= ker.sum()
    hp = np.concatenate([h[-len(kx):], h, h[:len(kx)]]); hs = np.convolve(hp, ker, 'same')[len(kx):-len(kx)]
    thr = 0.5*hs.max()
    peaks = [i for i in range(nbins) if hs[i] > thr and hs[i] >= hs[i-1] and hs[i] > hs[(i+1) % nbins]]
    F = np.abs(np.fft.rfft(h))[1:40]
    return len(peaks), int(np.argmax(F))+1, hs.max()/max(hs.mean(), 1e-30)

if __name__ == '__main__':
    geom, R, D, T = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    K = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    M = int(sys.argv[6]) if len(sys.argv) > 6 else 240
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 0
    z, G = setup(geom, R, D, K, M, seed)
    delta = D/K                                       # blob size ~ layer spacing
    rhs = make_rhs(geom, G, delta)
    y0 = np.concatenate([z.real, z.imag])
    ts = np.linspace(0, T, 121)
    sol = solve_ivp(rhs, (0, T), y0, t_eval=ts, rtol=1e-7, atol=1e-10, method='RK45')
    n = len(G)
    np.savez(f'roll_{geom}_R{R}_D{D}_s{seed}.npz', t=sol.t, y=sol.y, G=G, geom=geom, R=R, D=D)
    print(f'{geom} R={R} D={D} K={K} M={M} seed={seed} blob={delta}  (status {sol.status})')
    for t, y in zip(sol.t, sol.y.T):
        zz = y[:n]+1j*y[n:]
        c, mdom, contrast = clumps(zz, G)
        rad = np.median(np.abs(zz)); spread = np.std(np.abs(zz))
        print(f'  t={t:7.1f}  clumps={c:3d}  dominant m={mdom:3d}  peak/mean={contrast:5.2f}  median rho={rad:.3f} spread={spread:.3f}', flush=True)
