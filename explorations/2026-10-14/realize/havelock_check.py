"""Independent check of the unsheared trap (Havelock) thresholds for N = 2, 3: exact sympy Jacobian in Cartesian
coordinates of the full N-vortex system in the unit disk (with image/Robin self terms), rotating frame, then the
characteristic polynomial; threshold = where a mode leaves the imaginary axis. Compares with ring_blocks.S_window.
Usage: python havelock_check.py"""
import sympy as sp, numpy as np
from ring_blocks import S_window, S_trap

def jac_numeric(N, r0v):
    xs = sp.symbols(f'x0:{N}', real=True); ys = sp.symbols(f'y0:{N}', real=True)
    H = 0
    for i in range(N):
        H += sp.log(1 - xs[i]**2 - ys[i]**2)/2
        for j in range(i+1, N):
            H += sp.log((1 - (xs[i]*xs[j] + ys[i]*ys[j]))**2 + (xs[i]*ys[j] - ys[i]*xs[j])**2)/2 \
                 - sp.log((xs[i] - xs[j])**2 + (ys[i] - ys[j])**2)/2
    Om = sp.Symbol('Om')
    V = [sp.diff(H, ys[i]) + Om*ys[i] for i in range(N)] + [-sp.diff(H, xs[i]) - Om*xs[i] for i in range(N)]
    Z = list(xs) + list(ys)
    Jm = sp.Matrix(V).jacobian(Z)
    th = [2*sp.pi*k/N for k in range(N)]
    sub = {xs[k]: r0v*sp.cos(th[k]) for k in range(N)}; sub.update({ys[k]: r0v*sp.sin(th[k]) for k in range(N)})
    # Om from equilibrium: x-velocity of vortex 0 (at (r0,0)) in the rotating frame is 0 automatically; y-velocity:
    om = sp.solve(V[N].subs(sub), Om)[0]
    J = np.array(Jm.subs(sub).subs(Om, om).evalf(30).tolist(), dtype=float)
    return np.linalg.eigvals(J)

for N in (2, 3):
    lo, hi = 0.05, 0.95
    f = lambda r: S_trap(N, r, 0.0) - S_window(N, r)[0]     # >= 0 stable (blocks)
    g = lambda r: np.abs(jac_numeric(N, sp.Rational(str(round(r, 12)))).real).max()
    for _ in range(45):
        mid = 0.5*(lo + hi)
        if f(mid) >= 0: lo = mid
        else: hi = mid
    rb = lo
    print(f'N = {N}: blocks threshold r_H = {rb:.9f} (p = {rb*rb:.9f});  sympy max|Re eig| at 0.999 r_H: {g(0.999*rb):.2e},  at 1.001 r_H: {g(1.001*rb):.2e}')
