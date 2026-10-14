"""Snapshot figure: plane vs H^2 roll-up of the R=2, D=0.2 strip (needs matplotlib and the two .npz files).
Usage: python plot_rollup.py   ->  rollup_R2.png"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 4, figsize=(14, 7.2))
for row, (f, title) in enumerate([('roll_plane_R2.0_D0.2_s0.npz', 'Flat plane'),
                                  ('roll_h2_R2.0_D0.2_s0.npz', 'Hyperbolic plane (Poincare disk)')]):
    d = np.load(f, allow_pickle=True); n = len(d['G'])
    for col, ti in enumerate([0, 8, 40, 120]):
        y = d['y'][:, ti]; z = y[:n]+1j*y[n:]
        ax = axs[row, col]; ax.scatter(z.real, z.imag, s=0.6, c='#1f4e79')
        if row == 1:
            th = np.linspace(0, 2*np.pi, 400); ax.plot(np.cos(th), np.sin(th), lw=0.6, c='0.6'); lim = 1.02
        else:
            lim = 2.6
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f'{title}\nt = {d["t"][ti]:.0f}' if col == 0 else f't = {d["t"][ti]:.0f}', fontsize=10)
fig.suptitle('Strip of uniform vorticity, outer radius R = 2, thickness D = 0.2 (curvature radii): roll-up and merging', fontsize=11)
plt.tight_layout(); plt.savefig('rollup_R2.png', dpi=130)
