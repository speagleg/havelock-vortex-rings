"""H^2 rings realised in a Penning-Malmberg trap by a CENTRAL column (Havelock's ring-plus-centre problem).

A central vortex of strength g (in units of a ring column) is, for the ring, a background with stream function -g ln r,
whose shear at the ring is r (U'/r)' = 2 g / r0^2. By the lemma (flat_lemma.py) the trap ring with a FIXED central
vortex of strength
        g*(N, r0) = s*(N, r0) r0^2 / 2
has the H^2 ring's linearisation (time rescaled), mode by mode. The central vortex can move; by symmetry its
displacement couples only to ring modes k = +-1 (mod N). This script checks the full (N+1)-vortex linearisation
(exact, mpmath finite differences at 50 digits) against the H^2 ring, and computes Havelock's stability interval
[g_lo, g_hi] for the trap ring with a mobile centre, to compare with g*.
Usage: python ring_center.py -> ring_center.out"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, mpmath as mp
from ring_blocks import S_window, S_trap, S_H2, ring_spectrum
mp.mp.dps = 50

def vel(Z, gam, mode):
    """Z = [x..., y...] (mpmath), circulations gam (column = 1); trap (Robin, m = 1) or H2 (m = ((1-r^2)/2)^2, no Robin)."""
    n = len(gam); z = [mp.mpc(Z[i], Z[n + i]) for i in range(n)]; out_x, out_y = [], []
    for i in range(n):
        g = mp.mpc(0)
        for j in range(n):
            if i == j: continue
            g += gam[j]*(-z[j]/(1 - mp.conj(z[i])*z[j]) - 1/mp.conj(z[i] - z[j]))
        r2 = abs(z[i])**2
        if mode == 'trap': g += -gam[i]*z[i]/(1 - r2); m = 1
        else: m = ((1 - r2)/2)**2
        out_x.append(m*g.imag); out_y.append(-m*g.real)
    return out_x + out_y

def spectrum(N, r0, g, mode, center=True):
    th = [2*mp.pi*k/N for k in range(N)]
    xs = [r0*mp.cos(t) for t in th]; ys = [r0*mp.sin(t) for t in th]; gam = [1]*N
    if center: xs.append(mp.mpf(0)); ys.append(mp.mpf(0)); gam.append(g)
    n = len(gam); Z = xs + ys
    v = vel(Z, gam, mode); Om = v[n]/r0                       # y-velocity of column 0 at (r0, 0)
    def vrot(Z):
        v = vel(Z, gam, mode); return [v[i] + Om*Z[n + i] for i in range(n)] + [v[n + i] - Om*Z[i] for i in range(n)]
    assert max(abs(x) for x in vrot(Z)) < mp.mpf(10)**-40
    h = mp.mpf(10)**-20; J = mp.matrix(2*n, 2*n)
    for k in range(2*n):
        Zp = list(Z); Zm = list(Z); Zp[k] += h; Zm[k] -= h
        a, b = vrot(Zp), vrot(Zm)
        for i in range(2*n): J[i, k] = (a[i] - b[i])/(2*h)
    ev = mp.eig(J, left=False, right=False)
    return np.array([complex(e) for e in ev])

def maxre(ev): return np.abs(ev.real).max()

def g_interval(N, r0, gmin=-5.0, gmax=40.0, n=181):
    gs = np.linspace(gmin, gmax, n); st = np.array([maxre(spectrum(N, mp.mpf(r0), mp.mpf(g), 'trap')) < 1e-8 for g in gs])
    if not st.any(): return None
    idx = np.where(st)[0]
    def refine(a, b, a_stable):
        for _ in range(40):
            c = 0.5*(a + b); sc = maxre(spectrum(N, mp.mpf(r0), mp.mpf(c), 'trap')) < 1e-8
            if sc == a_stable: a = c
            else: b = c
        return 0.5*(a + b)
    runs = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
    out = []
    for r in runs:
        lo = gs[r[0]] if r[0] == 0 else refine(gs[r[0]], gs[r[0] - 1], True)
        hi = gs[r[-1]] if r[-1] == n - 1 else refine(gs[r[-1]], gs[r[-1] + 1], True)
        out.append((lo, hi))
    return out

if __name__ == '__main__':
    out = []; P = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
    P('1. Full (N+1)-vortex trap spectrum with g = g* vs the H^2 ring (rescaled by ((1-r0^2)/2)^2):')
    for N, r0 in [(5, 0.4), (8, 0.3), (8, 0.5), (10, 0.6), (12, 0.45), (12, 0.55)]:
        s_star = S_H2(N, r0) - S_trap(N, r0, 0.0); gs = s_star*r0**2/2
        eT = spectrum(N, mp.mpf(r0), mp.mpf(gs), 'trap')*((1 - r0**2)/2)**2
        eH = spectrum(N, mp.mpf(r0), 0, 'H2', center=False)
        D = np.abs(eH[:, None] - eT[None, :]); unmatched = np.sort(D.min(axis=0))[-4:]
        P(f'  N = {N:2d}, r0 = {r0}: g* = {gs:.5f};  max Re: H^2 {maxre(eH):.2e}, trap+centre {maxre(eT):.2e};'
          f'  H^2 eigenvalues found in trap spectrum to {D.min(axis=1).max():.1e}; 4 largest trap-side mismatches {np.round(unmatched, 4)}  (k = +-1 sector, where the centre moves)')
    P('\n2. Validation: N = 6 at r0 = 0.25 (the ring radius of Durkin & Fajans 2000, who quote -0.25 <= g <= 6.25):')
    iv = g_interval(6, 0.25, -1.0, 8.0, 91); P(f'  stable g-interval(s): {[(round(a, 4), round(b, 4)) for a, b in iv]}')
    P('\n3. Havelock interval: trap ring with a mobile central column is stable iff g in [g_lo, g_hi]; vs g*:')
    for N, r0 in [(3, 0.4), (5, 0.4), (8, 0.25), (8, 0.27), (8, 0.5), (10, 0.39), (10, 0.42), (12, 0.48), (12, 0.52)]:   # N = 16: ring_center_n16.py
        s_star = S_H2(N, r0) - S_trap(N, r0, 0.0); gs = s_star*r0**2/2
        iv = g_interval(N, r0, -1.0, 40.0, 83)
        P(f'  N = {N:2d}, r0 = {r0}: g* = {gs:.4f}; stable g-intervals {[(round(a, 4), round(b, 4)) for a, b in iv] if iv else None}; g* inside: {any(a <= gs <= b for a, b in iv) if iv else False};'
          f'  H^2 ring stable: {maxre(spectrum(N, mp.mpf(r0), 0, "H2", center=False)) < 1e-8}')
    open('ring_center.out', 'w').write('\n'.join(out) + '\n')
