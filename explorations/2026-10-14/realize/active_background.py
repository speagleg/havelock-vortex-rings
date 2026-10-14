"""Trap design with an ACTIVE background: the peaked background zeta_b = A(1 - r^2) is represented by Nb weak point
vortices (desingularised with a small blob radius delta for every interaction that involves a background vortex),
advected self-consistently with the 8 columns. Tests whether the equivalence survives background dynamics
(the passive version is nonlinear_sim.py). Background vortices are laid on concentric rings, equal circulation each,
ring populations proportional to the annulus circulation, random phases.
Usage: python active_background.py [Nb] [periods]  -> active_background.out, active_background.png"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
from ring_blocks import S_window, S_trap, S_H2

def make_background(A, Nb, rng, nrings=24):
    edges = np.linspace(0, 0.98, nrings + 1); circ = []
    for a, b in zip(edges[:-1], edges[1:]):   # circulation of zeta = A(1 - r^2) in the annulus
        F = lambda r: 2*np.pi*A*(r**2/2 - r**4/4); circ.append(F(b) - F(a))
    circ = np.array(circ); tot = circ.sum(); n = np.maximum(1, np.round(Nb*circ/tot)).astype(int)
    zs = []
    for k, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        rk = np.sqrt((a*a + b*b)/2); ph = rng.uniform(0, 2*np.pi)
        zs.append(rk*np.exp(1j*(ph + 2*np.pi*np.arange(n[k])/n[k])))
    zb = np.concatenate(zs); g = tot/(2*np.pi*len(zb))   # gam units: a column (circulation 2 pi) has gam = 1
    return zb, g, tot

def rhs_factory(Nc, gam, delta):
    def f(t, Z):
        n = len(gam); z = Z[:n] + 1j*Z[n:]
        zi, zj = z[:, None], z[None, :]
        d = zi - zj
        bg = (gam[:, None] < 1) | (gam[None, :] < 1)          # pairs involving a background vortex get a blob
        den = np.abs(d)**2 + np.where(bg, delta**2, 0.0)
        np.fill_diagonal(den, 1.0)
        free = -d/den                                          # grad_i of -(1/2) ln(|d|^2 + delta^2)
        img = -zj/(1 - np.conj(zi)*zj)                         # grad_i of ln|1 - z_i conj z_j|
        G = (free + img)*gam[None, :]; np.fill_diagonal(G, 0)
        g = G.sum(axis=1)
        r2 = np.abs(z)**2
        g += -gam*z/(1 - r2)                                   # self image (Robin), weighted by own circulation
        # velocity of vortex i: u + i v = -i * conj(... ) ; with H weighted by circulations: x_i' = dH/dy_i / gam_i
        return np.concatenate([g.imag, -g.real])
    return f

def run(s, Nb, periods, seed=3, delta=0.03, eps=1e-4, N=8, r0=0.5):
    rng = np.random.default_rng(seed)
    th = 2*np.pi*np.arange(N)/N
    zc = r0*np.exp(1j*th) + eps*(rng.standard_normal(N) + 1j*rng.standard_normal(N))
    if s > 0:
        zb, g, tot = make_background(2*s/r0**2, Nb, rng)
    else:
        zb, g, tot = np.zeros(0, complex), 0.0, 0.0
    z = np.concatenate([zc, zb]); gam = np.concatenate([np.ones(N), np.full(len(zb), g)])
    # in these units a column has circulation 2 pi (gam = 1), background vortices gam = g < 1
    f = rhs_factory(N, gam, delta)
    Om = f(0, np.concatenate([(r0*np.exp(1j*th)).real, np.zeros(len(zb)) + zb.real, (r0*np.exp(1j*th)).imag, zb.imag]))[len(gam)]/r0
    T = periods*2*np.pi/abs(Om)
    sol = solve_ivp(f, (0, T), np.concatenate([z.real, z.imag]), rtol=1e-7, atol=1e-9, method='RK45', dense_output=True)
    ts = np.linspace(0, sol.t[-1], 600); Z = sol.sol(ts); n = len(gam)
    zc_t = Z[:N] + 1j*Z[n:n + N]; r = np.abs(zc_t)
    return ts/(2*np.pi/abs(Om)), r.max(axis=0) - r.min(axis=0), len(zb), g, tot, sol.nfev

if __name__ == '__main__':
    Nb = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    periods = float(sys.argv[2]) if len(sys.argv) > 2 else 40
    N, r0 = 8, 0.5
    S0 = S_trap(N, r0, 0.0); s_star = S_H2(N, r0) - S0; s_c = S_window(N, r0)[0] - S0
    out = []; P = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
    P(f'N = {N}, r0 = {r0}: s* = {s_star:.3f}, s_c = {s_c:.3f}; Nb ~ {Nb}, blob delta = 0.03, {periods} periods')
    res = {}
    for name, s in [('no background', 0.0), ('active, s = s*', s_star), ('active, s = s_c - 2', s_c - 2), ('active, s = s_c + 4', s_c + 4)]:
        t, dev, nb, g, tot, nfev = run(s, Nb, periods)
        res[name] = (t, dev)
        i1 = np.searchsorted(t, 10)
        P(f'  {name:22s}: Nb = {nb}, background circulation {tot:.1f} (columns {2*np.pi*N:.1f}); radial spread: '
          f'initial {dev[0]:.1e}, max to 10 periods {dev[:i1].max():.2e}, max overall {dev.max():.2e}, final {dev[-1]:.2e}  [nfev {nfev}]')
    open('active_background.out', 'w').write('\n'.join(out) + '\n')
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    for name, (t, d) in res.items(): ax.semilogy(t, np.maximum(d, 1e-12), label=name, lw=1)
    ax.set_xlabel('ring rotation periods'); ax.set_ylabel('radial spread of the 8 columns')
    ax.set_title(f'Trap, N = 8, $r_0$ = 0.5, active point-vortex background ($N_b$ ~ {Nb})', fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig('active_background.png', dpi=150)
