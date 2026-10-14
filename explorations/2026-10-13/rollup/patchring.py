"""Ring of N equal vortex patches (or point vortices, a = 0) on H^2 at geodesic radius r: does it survive?
Each patch = hexagonal blob cluster (NRING rings: 19 blobs for NRING=2, 37 for 3) of geodesic radius a,
carried to its centre by a disk isometry. Small random perturbation (1e-3) of each centre.
Diagnostic: max over patches of |geodesic radius of patch centroid - r| at t = 0, T/2, T, and its max.
Usage: python patchring.py r a T N1 N2 ...        (a = 0 for point vortices)
       set NRING=3 in the environment for finer patches.
Small patches are stiff: a = 0.06 takes ~30 min per N at T = 40."""
import sys, os, numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rollup import make_rhs

NRING = int(os.environ.get('NRING', 2))

def hexpatch(a):
    pts = [0j]
    for ring in range(1, NRING+1):
        cnt = 6*ring
        for k in range(cnt):
            pts.append(ring/NRING*np.exp(2j*np.pi*k/cnt))
    return np.array(pts)*np.tanh(a/2)

def build(N, r, a, seed):
    rng = np.random.default_rng(seed)
    rho0 = np.tanh(r/2); z = []; owner = []
    base = hexpatch(a) if a > 0 else np.array([0j])
    for k in range(N):
        c = rho0*np.exp(2j*np.pi*k/N)*(1+1e-3*rng.standard_normal()) * np.exp(1e-3j*rng.standard_normal())
        for w in base:
            z.append((w + c)/(1 + np.conj(c)*w)); owner.append(k)
    z = np.array(z); owner = np.array(owner)
    G = np.full(len(z), 1.0/len(base))
    return z, G, owner

def radial_dev(z, owner, N, r):
    return max(abs(2*np.arctanh(min(abs(np.mean(z[owner == k])), 1-1e-15)) - r) for k in range(N))

if __name__ == '__main__':
    r, a, T = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    for N in [int(x) for x in sys.argv[4:]]:
        z, G, owner = build(N, r, a, 0)
        delta = a/NRING if a > 0 else 1e-6
        rhs = make_rhs('h2', G, delta)
        ts = np.linspace(0, T, 41)
        sol = solve_ivp(rhs, (0, T), np.concatenate([z.real, z.imag]), t_eval=ts, rtol=1e-8, atol=1e-11)
        n = len(G)
        dev = [radial_dev(y[:n]+1j*y[n:], owner, N, r) for y in sol.y.T]
        i_half = len(dev)//2
        print(f'r={r} a={a} NRING={NRING} N={N:2d}: radial deviation  t=0: {dev[0]:.1e}  t={ts[i_half]:.0f}: '
              f'{dev[i_half]:.1e}  t={T:.0f}: {dev[-1]:.1e}  max: {max(dev):.1e}', flush=True)
