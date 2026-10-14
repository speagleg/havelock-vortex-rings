"""Ergoregion maps: lab frame (|v| > c) vs co-rotating frame (|v - Omega x r| > c). Usage: python plot_ergo.py -> ergo.png"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from ergo import velocity, omega_ring

cases = [(6, 3.0, 0, 'Ring of 6, R = 3'), (6, 5.0, 0, 'Ring of 6, R = 5'), (19, 15.0, 22, 'Charge-22 core + ring of 19, R = 15')]
fig, axs = plt.subplots(2, 3, figsize=(13.5, 9))
for col, (N, R, q, title) in enumerate(cases):
    L = 1.25*max(R+N, q+N, 1/omega_ring(N, R, q)) if q else 2.4*R + 4
    x = np.linspace(-L, L, 900); X, Y = np.meshgrid(x, x)
    V = velocity(X, Y, N, R, q)
    for row, frame in enumerate(['lab', 'corot']):
        W = V if frame == 'lab' else V - 1j*omega_ring(N, R, q)*(X+1j*Y)
        ax = axs[row, col]
        ax.contourf(X, Y, (np.abs(W) > 1).astype(float), levels=[0.5, 1.5], colors=['#c9732c'], alpha=0.55)
        ax.contour(X, Y, np.abs(W), levels=[1.0], colors='#7a3f12', linewidths=0.8)
        w = R*np.exp(2j*np.pi*np.arange(N)/N); ax.plot(w.real, w.imag, 'o', ms=3, c='#1f4e79')
        if q: ax.plot(0, 0, '*', ms=10, c='#1f4e79')
        ax.set_aspect('equal'); ax.set_xlim(-L, L); ax.set_ylim(-L, L)
        lab = 'lab frame: |v| > c' if frame == 'lab' else f'co-rotating frame: |v - Omega x r| > c   (Omega = {omega_ring(N, R, q):.3f})'
        ax.set_title(f'{title}\n{lab}', fontsize=9)
        ax.tick_params(labelsize=7)
fig.suptitle('Acoustic ergoregions (shaded) of vortex rings in a uniform 2D condensate. Lengths in l_e = kappa/(2 pi c) ~ 1.4 healing lengths.', fontsize=10)
plt.tight_layout(); plt.savefig('ergo.png', dpi=120)
