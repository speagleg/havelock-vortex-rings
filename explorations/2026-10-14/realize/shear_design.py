"""Background shear as the single design knob (uses the lemma in flat_lemma.py).

Trap (Penning-Malmberg, unit wall radius, uniform B, columns of circulation 2 pi):
  a ring of N at radius r0 with background shear s = r dOmega_b/dr at the ring is stable iff s >= s_c(N, r0);
  the H^2 ring (hyperbolic radius rho, r0 = tanh(rho/2)) is the trap ring with s = s*(N, r0).
  At s = 0 this is Havelock's problem (1931): stable iff r0 < r_H(N) for N <= 6, never for N >= 7.
Channel (width pi, uniform Couette background of shear Sigma):
  stable iff Sigma >= Sigma_c(alpha, l') = -B_ch(pi) (critical mode q = pi, from ../rows/ item 11);
  the H^2 hypercycle row is Sigma = Sigma*(alpha, l') = 4 cot(a) sum_j W_phi(j) + csc^2(a)/2.
Usage: python shear_design.py  -> shear_design.out, shear_design.png"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'channel'))
import numpy as np
from channel_symbol import symbols, pieces

from ring_blocks import S_window, S_trap, S_H2

def s_crit(N, r0):
    """trap threshold: stable iff s >= s_c (exact Fourier blocks; an upper window edge, if any, is reported by S_window)."""
    lo, hi = S_window(N, r0); assert hi == np.inf
    return lo - S_trap(N, r0, 0.0)

def s_star(N, r0):
    return S_H2(N, r0) - S_trap(N, r0, 0.0)

def crossings(f, a, b, n=400):
    xs = np.linspace(a, b, n); v = np.array([f(x) for x in xs]); roots = []
    for i in range(n - 1):
        if v[i]*v[i+1] < 0:
            lo, hi = xs[i], xs[i+1]
            for _ in range(60):
                mid = 0.5*(lo + hi)
                if f(lo)*f(mid) <= 0: hi = mid
                else: lo = mid
            roots.append(0.5*(lo + hi))
    return roots, v

def channel_maps(alpha, lps):
    Sc, Ss = [], []
    for lp in lps:
        BH, CH, Bc, Cc, V = symbols(alpha, lp, np.array([np.pi]))
        j, a, b, d, Wp = pieces(alpha, lp)
        Sc.append(-Bc[0]); Ss.append(4/np.tan(alpha)*Wp.sum() + 0.5/np.sin(alpha)**2)
    return np.array(Sc), np.array(Ss)

if __name__ == '__main__':
    out = []
    P = lambda *a: (print(*a), out.append(' '.join(str(x) for x in a)))
    P('Havelock check (s = 0): the unsheared trap ring is stable iff s_c(N, r0) <= 0; sign changes of s_c in (0, 0.99):')
    for N in range(2, 9):
        roots, v = crossings(lambda r: s_crit(N, r), 0.005, 0.99)
        P(f'  N = {N}: s_c < 0 at r0 = 0.005: {v[0] < 0};  sign changes at r0 = {[round(x, 6) for x in roots]}'
          + (f'  (p = r0^2 = {roots[0]**2:.6f})' if roots else ''))
    P('\nTrap: H^2 ring <-> trap ring with background shear s* = s*(N, r0); H^2 stable iff s* >= s_c')
    cir = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows', 'circle_branch.out'), skiprows=1)
    r0s = np.linspace(0.02, 0.97, 200); ring = {}
    for N in [3, 5, 7, 8, 10, 12, 16]:
        sc = np.array([s_crit(N, r) for r in r0s]); ss = np.array([s_star(N, r) for r in r0s]); ring[N] = (sc, ss)
        roots, _ = crossings(lambda r: s_star(N, r) - s_crit(N, r), 0.02, 0.97)
        note = ''
        for NN, k, l in cir:
            if int(NN) == N: note = f';  note threshold r0* = tanh(arccoth(a*)/2) = {np.tanh(np.arctanh(1/k)/2):.6f}'
        P(f'  N = {N:2d}: s* - s_c changes sign at r0 = {[round(x, 6) for x in roots]}{note};  s*(0.5) = {s_star(N, 0.5):.3f}, s_c(0.5) = {s_crit(N, 0.5):.3f}')
    P('  physical scale: required background vorticity contrast zeta_mean_inside - zeta_at_ring = s* Gamma/(2 pi R^2);')
    P('  relative to a column of radius a: (s*/2)(a/R)^2. Durkin-Fajans 2000 trap (a/R = 0.025): factor 3.1e-4.')
    for N in [8, 10, 12]:
        for r0 in [0.3, 0.5, 0.7]:
            P(f'    N = {N}, r0 = {r0}: s* = {s_star(N, r0):7.3f},  contrast / column vorticity = {s_star(N, r0)/2*0.025**2:.2e}')
    P('\nChannel (width pi): Couette shear Sigma_c(alpha, l\') needed, and the H^2 value Sigma*')
    lps = np.geomspace(0.1, 6, 60); chan = {}
    for al in [np.pi/8, np.pi/6, np.pi/4, np.pi/3]:
        Sc, Ss = channel_maps(al, lps); chan[al] = (Sc, Ss)
        P(f'  alpha = {al:.4f} (kappa {np.cos(al):.3f}): Sigma_c at l\' = 0.1, 1, 6: {np.interp([0.1, 1, 6], lps, Sc).round(3)};'
          f'  Sigma* - Sigma_c >= 0 for l\' >= {lps[np.argmax(Ss - Sc >= 0)] if np.any(Ss - Sc >= 0) else None}')
    P('  minimum over the grid of Sigma_c:', round(min(c[0].min() for c in chan.values()), 4))
    open('shear_design.out', 'w').write('\n'.join(out) + '\n')

    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(ring)))
    for c, (N, (sc, ss)) in zip(cols, ring.items()):
        ax[0].plot(r0s, sc, color=c, label=f'N={N}: $s_c$'); ax[0].plot(r0s, ss, '--', color=c)
    ax[0].axhline(0, color='k', lw=0.6); ax[0].set_ylim(-10, 80)
    ax[0].set_xlabel('ring radius $r_0$ / wall radius'); ax[0].set_ylabel(r'background shear $r\,d\Omega_b/dr$  [$\Gamma/2\pi R^2$]')
    ax[0].set_title('Penning trap: solid = threshold $s_c$ (stable above),\ndashed = shear $s^*$ that reproduces the $H^2$ ring', fontsize=9)
    ax[0].legend(fontsize=7, ncol=2)
    for c, (al, (Sc, Ss)) in zip(plt.cm.plasma(np.linspace(0, 0.8, 4)), chan.items()):
        ax[1].loglog(lps, Sc, color=c, label=fr'$\kappa$={np.cos(al):.3f}: $\Sigma_c$'); ax[1].loglog(lps, Ss, '--', color=c)
    ax[1].set_xlabel("row spacing $\\ell'$ (channel width $\\pi$)"); ax[1].set_ylabel(r'Couette shear $\Sigma$  [$\Gamma/2\pi$ per unit$^2$]')
    ax[1].set_title('Channel: solid = threshold $\\Sigma_c$ (stable above),\ndashed = $\\Sigma^*$ of the $H^2$ row (crosses iff $\\kappa>1/\\sqrt{2}$)', fontsize=9)
    ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig('shear_design.png', dpi=150)
