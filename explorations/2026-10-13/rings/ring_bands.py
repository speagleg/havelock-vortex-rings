"""Threads 1 and 2 of the 2026-10-13 exploration, as one script (about 40 min on 1 core).
(1) Screened (K0) rings in the plane: which N stay stable, and the N=6 threshold radius.
(2) Ring of N around a central vortex of strength q: the stable q-window, log kernel vs K0 kernel.
Usage: python ring_bands.py"""
import sys; sys.path.insert(0, '.')
import numpy as np
from screen import analyse, ring

TOL = 1e-6
stable = lambda N, R, q, k: analyse(*ring(N, R, q), k)[2] < TOL

def check_havelock():
    for N in range(2, 10):
        print(f'  log ring N={N}:', 'stable' if analyse(*ring(N, 1.0), 'log')[2] < TOL else 'unstable')

def screened_rings():
    for N in range(3, 9):
        Rs = np.linspace(0.05, 40, 400)
        fl = [analyse(*ring(N, R), 'k0')[2] > TOL for R in Rs]
        ch = [round(Rs[i], 2) for i in range(len(Rs)-1) if fl[i] != fl[i+1]]
        print(f'  K0 ring N={N}: {"unstable" if fl[0] else "stable"} at small R; changes near R/lambda = {ch}')
    lo, hi = 1.0, 2.0
    for _ in range(50):
        m = (lo+hi)/2
        if analyse(*ring(6, m), 'k0')[2] > TOL: hi = m
        else: lo = m
    print(f'  N=6 K0 threshold R/lambda = {lo:.6f}')

def window(N, R, kind):
    qs = np.linspace(0.01, 60, 1200); s = [stable(N, R, q, kind) for q in qs]
    if not any(s): return None
    i = s.index(True); j = len(s)-1-s[::-1].index(True)
    def bis(a, b):
        for _ in range(40):
            m = (a+b)/2
            if stable(N, R, m, kind): b = m
            else: a = m
        return b
    lo = bis(qs[i-1], qs[i]) if i > 0 else 0.0
    hi = bis(qs[j+1], qs[j]) if j < len(qs)-1 else np.inf
    return lo, hi

if __name__ == '__main__':
    print('Havelock check (log kernel, plane; N=7 is marginal):'); check_havelock()
    print('Screened rings:'); screened_rings()
    print('Central vortex window, log kernel; closed forms: q_lo = (N^2-8N+8)/16 (even), (N-1)(N-7)/16 (odd); q_hi = (N-1)^2/4')
    for N in range(5, 13):
        lo, hi = window(N, 1.0, 'log')
        f = (N*N-8*N+8)/16 if N % 2 == 0 else (N-1)*(N-7)/16
        print(f'  N={N}: [{lo:.4f}, {hi:.4f}]  closed form [{max(f, 0):.4f}, {(N-1)**2/4:.4f}] (scan starts at q=0.01)')
    for R in [0.5, 1, 2, 3]:
        print(f'Central vortex window, K0 kernel, R/lambda = {R}')
        for N in range(5, 13):
            w = window(N, R, 'k0')
            print(f'  N={N}:', 'none' if w is None else f'[{w[0]:.3f}, {w[1]:.3f}]')
    for N in [7, 8]:
        lo, hi = 0.05, 6
        for _ in range(40):
            m = (lo+hi)/2
            if stable(N, m, 1.0, 'k0'): lo = m
            else: hi = m
        print(f'q=1, N={N}: stable for R/lambda < {lo:.4f}' + (' (search cap)' if lo > 5.9 else ''))
