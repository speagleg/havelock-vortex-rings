"""Nonlinear design simulations (passive background flow), and the realisable version with a central column.

TRAP. N = 8 columns (circulation 2 pi each) on a ring of radius r0 = 0.5 in a unit-radius trap. Havelock: unstable
(N >= 7). H^2: stable (r0 > r0* = 0.2595). Background: a centrally peaked, monotone (Rayleigh-stable) same-sign vorticity
    zeta_b(r) = A (1 - r^2),  A = 2 s / r0^2,  which has  zeta_mean_inside(r0) - zeta_b(r0) = s = r dOmega_b/dr at r0.
Runs: (i) H^2 itself (disk model, mobility ((1 - r^2)/2)^2); (ii) trap, no background; (iii) trap, s = s*;
(iv) trap, s = s_c - 2 (just below threshold); (v) trap, s = s_c + 2.
CHANNEL. Row at alpha = pi/8 (kappa = 0.924) with spacing l' = 1 (H^2-stable: threshold 0.529), supercell of M = 6,
width pi: (i) H^2 (strip model), (ii) channel, (iii) channel + uniform background vorticity -Sigma* (Couette shear Sigma*).
The background is held fixed (axisymmetric / uniform). For the channel, uniform vorticity is exactly invariant under
any incompressible flow, so (iii) is an exact solution of 2D Euler with point vortices. For the trap a non-uniform
background is not invariant; see active_background.py.
Usage: python nonlinear_sim.py  -> nonlinear_sim.out, nonlinear_sim.png"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
from ring_blocks import S_window, S_trap, S_H2

# ---------------- trap ----------------
def disk_rhs(mode, N, A=0.0):
    def f(t, Z):
        z = Z[:N] + 1j*Z[N:]
        zi, zj = z[:, None], z[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            g = -zj/(1 - np.conj(zi)*zj) - 1/np.conj(zi - zj)
        np.fill_diagonal(g, 0); g = g.sum(axis=1)
        r2 = np.abs(z)**2
        if mode == 'H2':
            m = ((1 - r2)/2)**2
        else:
            m = 1.0
            g += -z/(1 - r2)                         # Robin: grad (1/2) ln(1 - r^2)
            g += -A*(0.5 - r2/4)*z                   # background: U_b' = -A (r/2 - r^3/4), grad = U_b' z/r
        return np.concatenate([m*g.imag, -m*g.real])
    return f

def ring_run(mode, N, r0, A, T, seed=0, eps=1e-4):
    rng = np.random.default_rng(seed); th = 2*np.pi*np.arange(N)/N
    z0 = r0*np.exp(1j*th) + eps*(rng.standard_normal(N) + 1j*rng.standard_normal(N))
    sol = solve_ivp(disk_rhs(mode, N, A), (0, T), np.concatenate([z0.real, z0.imag]), rtol=1e-10, atol=1e-12,
                    dense_output=True, method='DOP853')
    ts = np.linspace(0, T, 2000); Z = sol.sol(ts); z = Z[:N] + 1j*Z[N:]
    # shape deviation: spread of radii and of consecutive angular gaps
    r = np.abs(z); dev = np.maximum(r.max(axis=0) - r.min(axis=0), 0)
    return ts, dev

def trap_rhs_gam(gam):
    gam = np.asarray(gam, float); n = len(gam)
    def f(t, Z):
        z = Z[:n] + 1j*Z[n:]; zi, zj = z[:, None], z[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            G = (-zj/(1 - np.conj(zi)*zj) - 1/np.conj(zi - zj))*gam[None, :]
        np.fill_diagonal(G, 0); g = G.sum(axis=1) - gam*z/(1 - np.abs(z)**2)
        return np.concatenate([g.imag, -g.real])
    return f

def center_run(N, r0, g, T, seed=0, eps=1e-4):
    rng = np.random.default_rng(seed); th = 2*np.pi*np.arange(N)/N
    z0 = np.append(r0*np.exp(1j*th), 0) + eps*(rng.standard_normal(N + 1) + 1j*rng.standard_normal(N + 1))
    sol = solve_ivp(trap_rhs_gam(np.append(np.ones(N), g)), (0, T), np.concatenate([z0.real, z0.imag]),
                    rtol=1e-10, atol=1e-12, dense_output=True, method='DOP853')
    ts = np.linspace(0, T, 2000); Z = sol.sol(ts); z = Z[:N] + 1j*Z[N + 1:2*N + 1]
    r = np.abs(z); return ts, r.max(axis=0) - r.min(axis=0)

# ---------------- channel ----------------
def strip_rhs(mode, M, period, Sigma=0.0, J=30):
    ks = np.arange(-J, J + 1)
    def f(t, Z):
        X, Y = Z[:M], Z[M:]
        dx = X[:, None, None] - X[None, :, None] - ks[None, None, :]*period
        dx = np.clip(dx, -600, 600)
        y1, y2 = Y[:, None, None], Y[None, :, None]
        ch, sh = np.cosh(dx), np.sinh(dx)
        P, Q = ch - np.cos(y1 + y2), ch - np.cos(y1 - y2)
        with np.errstate(divide='ignore', invalid='ignore'):
            Wx = 0.5*(sh/P - sh/Q); Wy = 0.5*(np.sin(y1 + y2)/P - np.sin(y1 - y2)/Q)
        self_mask = (np.eye(M)[:, :, None] == 1) & (ks[None, None, :] == 0)
        Wx = np.where(self_mask, 0, Wx); Wy = np.where(self_mask, 0, Wy)
        Wx = np.where(np.eye(M)[:, :, None] == 1, 0, Wx)   # self images exert no x-force
        Hx, Hy = Wx.sum(axis=(1, 2)), Wy.sum(axis=(1, 2))
        if mode == 'H2':
            m = np.sin(Y)**2
        else:
            m = 1.0; Hy = Hy + 0.5/np.tan(Y) + Sigma*Y          # Robin + Couette (U_b' = Sigma y up to a constant)
        return np.concatenate([m*Hy, -m*Hx])
    return f

def row_run(mode, alpha, lp, M, Sigma, T, seed=1, eps=1e-4):
    rng = np.random.default_rng(seed)
    X0 = np.arange(M)*lp + eps*rng.standard_normal(M); Y0 = alpha + eps*rng.standard_normal(M)
    sol = solve_ivp(strip_rhs(mode, M, M*lp, Sigma), (0, T), np.concatenate([X0, Y0]), rtol=1e-9, atol=1e-12,
                    dense_output=True, method='DOP853')
    ts = np.linspace(0, sol.t[-1], 1500); Z = sol.sol(ts)
    return ts, Z[M:].max(axis=0) - Z[M:].min(axis=0)

if __name__ == '__main__':
    out = []; P = lambda *a: (print(*a), out.append(' '.join(str(x) for x in a)))
    N, r0 = 8, 0.5
    Strap0 = S_trap(N, r0, 0.0); s_star = S_H2(N, r0) - Strap0; s_c = S_window(N, r0)[0] - Strap0
    P(f'TRAP N = {N}, r0 = {r0}: s* = {s_star:.4f}, s_c = {s_c:.4f}; background circulation for s*: pi s*/r0^2 = {np.pi*s_star/r0**2:.1f} '
      f'(columns: {2*np.pi*N:.1f})')
    # time unit: one ring rotation period in the trap with s = s* (Omega from the rhs)
    f = disk_rhs('trap', N, 2*s_star/r0**2); z = r0*np.exp(2j*np.pi*np.arange(N)/N)
    Om = f(0, np.concatenate([z.real, z.imag]))[N]/r0; Tp = 2*np.pi/abs(Om)
    scale = ((1 - r0**2)/2)**2
    runs = {'H^2 (time rescaled)': ('H2', 0.0, 1/scale), 'trap, no background': ('trap', 0.0, 1.0),
            'trap, s = s*': ('trap', 2*s_star/r0**2, 1.0),
            'trap, s = s_c - 2': ('trap', 2*(s_c - 2)/r0**2, 1.0), 'trap, s = s_c + 2': ('trap', 2*(s_c + 2)/r0**2, 1.0)}
    res = {}
    for name, (mode, A, tf) in runs.items():
        ts, dev = ring_run(mode, N, r0, A, 100*Tp*tf)
        res[name] = (ts/tf/Tp, dev)
        P(f'  {name:22s}: max radial spread over 100 periods = {dev.max():.2e} (initial {dev[0]:.1e})')
    P(f'  (period Tp = {Tp:.4f} in trap units at s = s*)')
    # channel
    al, lp, M = np.pi/8, 1.0, 6
    from channel_symbol_shim import Sigmas
    Sstar, Sc = Sigmas(al, lp)
    P(f'CHANNEL alpha = pi/8, l\' = {lp}: Sigma* = {Sstar:.4f}, Sigma_c = {Sc:.4f}')
    s2 = np.sin(al)**2; cres = {}
    for name, (mode, Sg, tf) in {'H^2 (time rescaled)': ('H2', 0, 1/s2), 'channel, no shear': ('ch', 0, 1.0),
                                 'channel, Sigma = Sigma*': ('ch', Sstar, 1.0)}.items():
        ts, dev = row_run(mode, al, lp, M, Sg, 60*tf)
        cres[name] = (ts/tf, dev)
        P(f'  {name:24s}: max transverse spread up to t = 60 = {dev.max():.2e} (initial {dev[0]:.1e})')
    # trap with a mobile central column (the realisable design)
    g_star = s_star*r0**2/2; g_lo = s_c*r0**2/2
    P(f'TRAP + CENTRAL COLUMN, N = {N}, r0 = {r0}: g* = {g_star:.4f}, Havelock lower limit g_lo = {g_lo:.4f}')
    f2 = trap_rhs_gam(np.append(np.ones(N), g_star)); z = r0*np.exp(2j*np.pi*np.arange(N)/N)
    Zc = np.concatenate([np.append(z, 0).real, np.append(z, 0).imag]); Om2 = f2(0, Zc)[N + 1]/r0; Tp2 = 2*np.pi/abs(Om2)
    gres = {}
    for name, g in [('centre g = g*', g_star), ('centre g = g_lo - 0.1', g_lo - 0.1), ('centre g = g_lo + 0.1', g_lo + 0.1)]:
        ts, dev = center_run(N, r0, g, 100*Tp2)
        gres[name] = (ts/Tp2, dev)
        P(f'  {name:22s}: max radial spread over 100 periods = {dev.max():.2e} (initial {dev[0]:.1e})')
    open('nonlinear_sim.out', 'w').write('\n'.join(out) + '\n')

    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.2))
    for name, (t, d) in res.items():
        if name.startswith('H^2'): ax[0].semilogy(t, np.maximum(d, 1e-12), 'k--', label=name + ' (on top of s = s*)', lw=0.8, zorder=5)
        else: ax[0].semilogy(t, np.maximum(d, 1e-12), label=name, lw=1)
    ax[0].set_xlabel('ring rotation periods'); ax[0].set_ylabel('radial spread of the ring')
    ax[0].set_title(f'Trap, N = {N} at $r_0$ = {r0} (Havelock-unstable, $H^2$-stable)', fontsize=9); ax[0].legend(fontsize=7)
    for name, (t, d) in cres.items():
        if name.startswith('H^2'): ax[1].semilogy(t, np.maximum(d, 1e-12), 'k--', label=name + ' (on top of Sigma*)', lw=0.8, zorder=5)
        else: ax[1].semilogy(t, np.maximum(d, 1e-12), label=name, lw=1)
    for name, (t, d) in gres.items():
        ax[2].semilogy(t, np.maximum(d, 1e-12), label=name, lw=1)
    ax[2].set_xlabel('ring rotation periods'); ax[2].set_ylabel('radial spread of the ring')
    ax[2].set_title(f'Trap, N = {N} at $r_0$ = {r0} + mobile central column', fontsize=9); ax[2].legend(fontsize=7)
    ax[1].set_xlabel('time (channel units)'); ax[1].set_ylabel('transverse spread of the row')
    ax[1].set_title(r"Channel row, $\kappa$ = 0.924, $\ell'$ = 1 (supercell of 6)", fontsize=9); ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig('nonlinear_sim.png', dpi=150)
