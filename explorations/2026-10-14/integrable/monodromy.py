"""Numerical Painleve (monodromy) test for the reduced family
    db_j/dt = b_j [g(Y_j) b_{j-1} - g(Y_{j+1}) b_{j+1}],   dY_j/dt = b_{j-1} - b_j      (ring of N sites),
with g = (ln phi)'. Members: phi = (1+Y^2)^{lam/2}, g = lam Y/(1+Y^2)  (lam = -2: vortex row, -4: earlier control,
+2: phi = 1+Y^2, chaotic); reference g = 1 (phi = e^Y, the Volterra lattice, which has the Painleve property).
Paths that come close to a pole (|z| > 50) or whose result depends on the tolerance are discarded.
Test: from a real initial condition integrate in complex time to t1 = T + iA along two paths (real-then-imaginary and
imaginary-then-real). Single-valued solutions give the same endpoint; branching around enclosed singularities does not.
Usage: python monodromy.py -> monodromy.out"""
import numpy as np
from scipy.integrate import solve_ivp
N = 3
def make_rhs(g):
    def F(z):
        b, Y = z[:N], z[N:]
        bm, bp = np.roll(b, 1), np.roll(b, -1); gY = g(Y); gYp = np.roll(gY, -1)
        return np.concatenate([b*(gY*bm - gYp*bp), bm - b])
    return F
def path(F, z0, t0, t1, zmax=50.0):
    """integrate along the segment t0 -> t1; reject if |z| exceeds zmax anywhere (close to a pole) or if two
    tolerances disagree (unreliable)."""
    dt = t1 - t0; res = []
    for rtol in (1e-10, 1e-13):
        sol = solve_ivp(lambda s, z: F(z)*dt, (0, 1), z0, rtol=rtol, atol=rtol*1e-2, method='DOP853', dense_output=True)
        if sol.status != 0: return None
        zz = sol.sol(np.linspace(0, 1, 400))
        if not np.all(np.isfinite(zz)) or np.abs(zz).max() > zmax: return None
        res.append(sol.y[:, -1])
    if np.abs(res[0] - res[1]).max() > 1e-7*(1 + np.abs(res[1]).max()): return None
    return res[1]
def test(g, ntrials=150, seed=0):
    F = make_rhs(g); rng = np.random.default_rng(seed); mism = []; fails = 0
    for k in range(ntrials):
        z0 = np.concatenate([rng.uniform(0.2, 1.0, N), rng.uniform(-1.5, 1.5, N)]).astype(complex)
        T, A = rng.uniform(0.3, 3.0), rng.uniform(0.3, 3.0)
        za = path(F, z0, 0, T); za = None if za is None else path(F, za, T, T + 1j*A)
        zb = path(F, z0, 0, 1j*A); zb = None if zb is None else path(F, zb, 1j*A, T + 1j*A)
        if za is None or zb is None: fails += 1; continue
        mb = np.abs(za[:N] - zb[:N]).max()/(1 + np.abs(za[:N]).max())
        dY = za[N:] - zb[N:]
        mism.append((mb, np.abs(dY).max(), np.abs(np.exp(dY) - 1).max(), np.abs((dY/(2j*np.pi)).imag).max() if True else 0))
    mism = np.array(mism)
    return mism, fails
if __name__ == '__main__':
    cases = {'Volterra (g = 1)': lambda Y: np.ones_like(Y),
             'lam = -2 (vortex row)': lambda Y: -2*Y/(1 + Y**2),
             'lam = -4 (control)': lambda Y: -4*Y/(1 + Y**2),
             'lam = +2 (phi = 1+Y^2)': lambda Y: 2*Y/(1 + Y**2),
             'lam = -1': lambda Y: -1*Y/(1 + Y**2),
             'lam = +1': lambda Y: 1*Y/(1 + Y**2)}
    out = open('monodromy.out', 'w')
    for name, g in cases.items():
        m, f = test(g)
        mb, mY = m[:, 0], m[:, 1]
        # Y mismatches that are integer multiples of 2 pi i would be logarithmic monodromy (harmless if Y enters via e^Y)
        k = m[:, 1]/(2*np.pi)
        line = (f'{name:24s}: {len(m)} loops ({f} failed); b-mismatch: fraction > 1e-6 = {np.mean(mb > 1e-6):.2f} (median {np.median(mb):.1e});'
                f'  Y-mismatch: fraction > 1e-6 = {np.mean(mY > 1e-6):.2f}, of which |dY|/2pi near an integer: '
                f'{np.mean(np.abs(k[mY > 1e-6] - np.round(k[mY > 1e-6])) < 1e-4) if np.any(mY > 1e-6) else float("nan"):.2f}')
        print(line, flush=True); out.write(line + '\n'); out.flush()
