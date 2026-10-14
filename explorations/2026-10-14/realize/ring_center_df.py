"""The Durkin-Fajans geometry (normalized ring radius 0.25, N = 7-10, central column strength g in [0, 4] accessible):
Havelock's lower limit g_lo(N, 0.25) vs the H^2-equivalent strength g*(N, 0.25), and a check that the trap ring with g = g*(N, r0)
changes stability at the H^2 threshold r0*(N) of the note. Usage: python ring_center_df.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, mpmath as mp
from ring_center import spectrum, maxre
from ring_blocks import S_trap, S_H2
mp.mp.dps = 50
stable = lambda N, r0, g: maxre(spectrum(N, mp.mpf(r0), mp.mpf(g), 'trap')) < 1e-8
gstar = lambda N, r0: (S_H2(N, r0) - S_trap(N, r0, 0.0))*r0**2/2

def g_lo(N, r0, a=-1.0, b=6.0):
    """lowest stable g (bisection between an unstable a and a stable b)."""
    assert not stable(N, r0, a) and stable(N, r0, b), (N, r0)
    for _ in range(40):
        c = 0.5*(a + b)
        if stable(N, r0, c): b = c
        else: a = c
    return b

cir = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows', 'circle_branch.out'), skiprows=1)
r0star = {int(N): float(np.tanh(np.arctanh(1/k)/2)) for N, k, l in cir}
out = []; P = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
P(' N   g_lo(N, 0.25)   g*(N, 0.25)   H^2 ring at r0 = 0.25 | trap with g = g*(N, r0) at r0 = r0*(N)(1 -+ 0.003)')
for N in [7, 8, 9, 10]:
    gl = g_lo(N, 0.25); gs = gstar(N, 0.25)
    if N in r0star:
        rs = r0star[N]
        chk = ' / '.join('stable' if stable(N, r, gstar(N, r)) else 'unstable' for r in (rs*0.997, rs*1.003))
        chk = f'r0* = {rs:.6f}: {chk}'
    else:
        chk = 'N = 7: H^2 ring stable at every radius'
    P(f'{N:2d}   {gl:12.6f}   {gs:11.6f}   {"stable" if gs >= gl else "unstable":>8s}            | {chk}')
open('ring_center_df.out', 'w').write('\n'.join(out) + '\n')
