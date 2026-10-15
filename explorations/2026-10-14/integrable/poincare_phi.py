"""Poincare sections (N = 3 ring, section D_1 = L/3 upward, plot (D_2, Y_1 - Y_2)) for H = sum phi(Y_j) e^{-D_j}:
phi = 2/(1+Y^2) (vortex row), 1/(1+Y^2)^2 (control), 1 + Y^2, and a decaying non-rational phi = exp(-Y^2/2) (Gaussian).
Usage: python poincare_phi.py -> poincare_phi.png"""
import numpy as np
from scipy.integrate import solve_ivp
from phi_family import rhs_phi
N = 3; L = 3.0
PHIS = {'2/(1+Y^2)  (vortex row)': (lambda Y: 2/(1 + Y**2), lambda Y: -4*Y/(1 + Y**2)**2),
        '1/(1+Y^2)^2': (lambda Y: 1/(1 + Y**2)**2, lambda Y: -4*Y/(1 + Y**2)**3),
        'exp(-Y^2/2)': (lambda Y: np.exp(-Y**2/2), lambda Y: -Y*np.exp(-Y**2/2)),
        '1 + Y^2': (lambda Y: 1 + Y**2, lambda Y: 2*Y)}
def section(rhs, ics, T):
    out = []
    for (Y1, Y2, Y3, d2) in ics:
        z0 = np.array([0, L/3, L/3 + d2, Y1, Y2, Y3])
        ev = lambda t, z: (z[1] - z[0]) - L/3; ev.direction = 1
        sol = solve_ivp(rhs, (0, T), z0, events=ev, rtol=1e-10, atol=1e-12, method='DOP853')
        out.append(np.array([[z[2] - z[1], z[3] - z[4]] for z in sol.y_events[0]]))
    return out
if __name__ == '__main__':
    rng = np.random.default_rng(3); Y0 = -1.2
    ics = [(Y0 + rng.uniform(-0.5, 0.5), Y0 + rng.uniform(-0.5, 0.5), Y0 + rng.uniform(-0.5, 0.5), L/3 + rng.uniform(-0.25, 0.25)) for _ in range(12)]
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 4, figsize=(16, 4)); cols = plt.cm.viridis(np.linspace(0, 0.95, len(ics)))
    for a, (name, (phi, dphi)) in zip(ax, PHIS.items()):
        secs = section(rhs_phi(phi, dphi, L), ics, 3000.0)
        for c, sct in zip(cols, secs):
            if len(sct): a.plot(sct[:, 0], sct[:, 1], '.', ms=1.0, color=c)
        a.set_title(f'phi = {name}', fontsize=9); a.set_xlabel('$D_2$ at $D_1 = L/3$'); a.set_ylabel('$Y_1 - Y_2$')
        a.spines[['top', 'right']].set_visible(False)
        print(name, [len(s) for s in secs], flush=True)
    fig.tight_layout(); fig.savefig('poincare_phi.png', dpi=140); print('saved')
