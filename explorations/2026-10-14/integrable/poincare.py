"""Poincare sections for N = 3 periodic chains at fixed impulse (sum Y), section D_1 = L/3 crossed upwards; each
trajectory has its own energy, so a regular one traces a curve and a chaotic one scatters,
plotting (D_2, Y_1 - Y_2). Integrable -> smooth closed curves; chaotic -> scattered points.
rows near Y0 = -2 (kappa = 2/sqrt 5, stable side), energy slightly above the row's.
(a) the Toda-limit vortex row, (b) the full tanh-kernel row (periodic images). Usage: python poincare.py -> poincare.png"""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'quantum'))
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from toda_test import nn_rhs, full_rhs
N = 3; L = 3.0; Y0 = -2.0; I0 = 3*Y0
def H_nn(z):
    u, Y = z[:N], z[N:]; s = 1/np.sqrt(1 + Y**2); D = np.roll(u, -1) - u; D[-1] += L
    return np.sum(2*s*np.roll(s, -1)*np.exp(-D))
def H_full(z, J=40):
    u, Y = z[:N], z[N:]; ph = np.pi/2 + np.arctan(Y); E = 0.0
    for i in range(N):
        for j in range(N):
            for k in range(-J, J + 1):
                if i == j and k == 0: continue
                du = u[i] - u[j] - k*L
                if abs(du) > 600: continue
                ch = np.cosh(du); w = 0.5*np.log((ch - np.cos(ph[i] + ph[j]))/(ch - np.cos(ph[i] - ph[j])))
                E += 0.5*w
    return E
def section(rhs, Hf, E, ics, T):
    pts = []
    for (Y1, Y2, d2) in ics:
        Y3 = I0 - Y1 - Y2
        z0 = np.array([0, L/3, L/3 + d2, Y1, Y2, Y3])
        ev = lambda t, z: (z[1] - z[0]) - L/3; ev.direction = 1
        sol = solve_ivp(rhs, (0, T), z0, events=ev, rtol=1e-10, atol=1e-12, method='DOP853')
        for z in sol.y_events[0]:
            pts.append((z[2] - z[1], z[3] - z[4], len(pts)))
        pts.append(None)
    return pts
if __name__ == '__main__':
    rng = np.random.default_rng(3)
    ics = [(Y0 + rng.uniform(-0.4, 0.4), Y0 + rng.uniform(-0.4, 0.4), L/3 + rng.uniform(-0.25, 0.25)) for _ in range(12)]
    z_ref = np.array([0, 1.12, 2.0, Y0 + 0.35, Y0 - 0.2, Y0 - 0.15])
    res = {}
    for name, rhs, Hf, T in [('Toda-limit row', nn_rhs(L), H_nn, 3000.0), ('full tanh row', full_rhs(L, N), H_full, 600.0)]:
        E = Hf(z_ref); res[name] = section(rhs, Hf, E, ics, T)
        print(name, 'E =', round(E, 5), 'points', sum(1 for p in res[name] if p), flush=True)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    cols = plt.cm.viridis(np.linspace(0, 0.95, len(ics)))
    for a, (name, pts) in zip(ax, res.items()):
        k = 0; xs, ys = [], []
        for p in pts:
            if p is None:
                if xs: a.plot(xs, ys, '.', ms=1.2, color=cols[k % len(cols)]); k += 1
                xs, ys = [], []
            else: xs.append(p[0]); ys.append(p[1])
        a.set_title(f'{name} (N = 3, fixed impulse; one colour per trajectory)', fontsize=9); a.set_xlabel('$D_2$ at $D_1 = L/3$'); a.set_ylabel('$Y_1 - Y_2$')
        a.spines[['top', 'right']].set_visible(False); a.set_xlim(0.6, 1.25)
    fig.tight_layout(); fig.savefig('poincare.png', dpi=150); print('saved')
