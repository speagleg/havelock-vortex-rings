"""Dynamical check of the Landau picture: integrate the exact H^2 (tanh model) dynamics of a periodic supercell of
M vortices on a hypercycle, just past the threshold (l = r l*, r < 1), from a tiny zigzag perturbation.
M = 2 keeps only the q = 0 and q = pi modes, so the dynamics is exactly the fixed-(l', Y0) Landau problem
(period and impulse conserved). Supercritical (a4 > 0): the staggered amplitude eta = (1/M) sum (-1)^j (Y_j - <Y>) stays within sqrt(-a2/a4)
(energy conservation in the double well). Subcritical (a4 < 0): it runs away.
Strip coordinates, H^2 mobility sin^2(phi), no self term; images along u summed (period M l').
Usage: python zigzag_dynamics.py -> zigzag_dynamics.out, zigzag_dynamics.png"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows'))
import numpy as np, mpmath as mp
from scipy.integrate import solve_ivp
from exact_symbol import lstar
from landau import e_cell
from symbol import symbol

def rhs_factory(M, period, J=30):
    ks = np.arange(-J, J + 1)
    def f(t, Z):
        X, Ph = Z[:M], Z[M:]
        dx = np.clip(X[:, None, None] - X[None, :, None] - ks[None, None, :]*period, -600, 600)
        p1, p2 = Ph[:, None, None], Ph[None, :, None]
        ch, sh = np.cosh(dx), np.sinh(dx); A, B = ch - np.cos(p1 + p2), ch - np.cos(p1 - p2)
        with np.errstate(divide='ignore', invalid='ignore'):
            Wx = 0.5*(sh/A - sh/B); Wp = 0.5*(np.sin(p1 + p2)/A - np.sin(p1 - p2)/B)
        mask = (np.eye(M)[:, :, None] == 1) & (ks[None, None, :] == 0)
        Wx = np.where(mask, 0, Wx); Wp = np.where(mask, 0, Wp)
        Wx = np.where(np.eye(M)[:, :, None] == 1, 0, Wx)
        m = np.sin(Ph)**2
        return np.concatenate([m*Wp.sum(axis=(1, 2)), -m*Wx.sum(axis=(1, 2))])
    return f

def run(kappa, r, M=2, eta0=1e-3, seed=0, efolds=40):
    s = np.sqrt(1 - kappa**2); alpha = np.arccos(kappa); lp = lstar(kappa)*s*r
    rng = np.random.default_rng(seed); sgn = (-1.0)**np.arange(M)
    Y = -kappa/s + eta0*sgn + 1e-6*rng.standard_normal(M)
    Z0 = np.concatenate([np.arange(M)*lp + 1e-6*rng.standard_normal(M), np.pi/2 + np.arctan(Y)])
    B, C, _ = symbol(kappa, lp, np.array([np.pi])); lam = np.sqrt(max(B[0]*C[0], 1e-30))
    T = efolds/lam
    sol = solve_ivp(rhs_factory(M, M*lp), (0, T), Z0, rtol=1e-10, atol=1e-12, method='DOP853', dense_output=True)
    ts = np.linspace(0, sol.t[-1], 3000); Ph = sol.sol(ts)[M:]; Yt = -1/np.tan(Ph)
    eta = (sgn[:, None]*(Yt - Yt.mean(axis=0))).mean(axis=0)
    a = [float(mp.diff(lambda e: e_cell(mp.mpf(-kappa/s) + e, mp.mpf(-kappa/s) - e, mp.mpf(lp)), 0, n)/mp.factorial(n)) for n in (2, 4)]
    bound = np.sqrt(-a[0]/a[1]) if a[1] > 0 and a[0] < 0 else np.nan
    return ts*lam, eta, a, bound, sol.status

if __name__ == '__main__':
    out = []; P = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
    res = {}
    for kappa, r in [(0.98, 0.98), (0.94, 0.98), (0.9, 0.98), (0.8, 0.98)]:
        t, eta, a, bound, st = run(kappa, r)
        res[(kappa, r)] = (t, eta, bound)
        P(f'kappa {kappa}, l = {r} l*: a2 = {a[0]:.4e}, a4 = {a[1]:.4e};  max |eta| = {np.abs(eta).max():.4e}'
          f'  (Landau bound sqrt(-a2/a4) = {bound:.4e});  final |eta| = {abs(eta[-1]):.3e}; integrator status {st}')
    open('zigzag_dynamics.out', 'w').write('\n'.join(out) + '\n')
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    COL = ['#2a78d6', '#1baf7a', '#eb6834', '#e34948']
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for c, ((kappa, r), (t, eta, bound)) in zip(COL, res.items()):
        ax.semilogy(t, np.abs(eta) + 1e-12, color=c, lw=1.6, label=f'$\\kappa$ = {kappa}' + ('' if np.isnan(bound) else ' (Landau bound dotted)'))
        if not np.isnan(bound): ax.axhline(bound, color=c, lw=1, ls=':')
    ax.set_xlabel('time x linear growth rate'); ax.set_ylabel('zigzag amplitude |$\\eta$|')
    ax.set_title('2-vortex cell at $\\ell = 0.98\\,\\ell^*$: orbits stop at the Landau bound (dotted)\nfor $\\kappa > \\kappa_t = 0.920$ and overshoot below it', fontsize=9)
    ax.spines[['top', 'right']].set_visible(False); ax.grid(color='#e6e5e0', lw=0.6); ax.legend(fontsize=7, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2)
    fig.tight_layout(); fig.savefig('zigzag_dynamics.png', dpi=150)
