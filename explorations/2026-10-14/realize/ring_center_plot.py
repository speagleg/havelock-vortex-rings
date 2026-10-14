"""Figure: Havelock's lower central-column limit g_lo(N, r0) (= s_c r0^2/2 from the exact Fourier blocks; equal to the
mobile-centre value, ring_center.py) and the H^2-equivalent strength g*(N, r0) = s* r0^2/2. The H^2 ring is stable iff
g* >= g_lo; the curves cross at the note's thresholds r0*(N). Usage: python ring_center_plot.py -> ring_center.png"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ring_blocks import S_window, S_trap, S_H2
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
r0s = np.linspace(0.03, 0.8, 160)
fig, ax = plt.subplots(figsize=(6.4, 4.6))
for c, N in zip(plt.cm.viridis(np.linspace(0, 0.9, 5)), [7, 8, 10, 12, 16]):
    gl = [(S_window(N, r)[0] - S_trap(N, r, 0.0))*r*r/2 for r in r0s]
    gs = [(S_H2(N, r) - S_trap(N, r, 0.0))*r*r/2 for r in r0s]
    ax.plot(r0s, gl, color=c, label=f'N = {N}: Havelock $g_{{lo}}$'); ax.plot(r0s, gs, '--', color=c)
ax.axvline(0.25, color='grey', lw=0.7, ls=':'); ax.axhspan(0, 4, color='grey', alpha=0.08)
ax.text(0.255, 10.0, 'Durkin & Fajans:\n$r_0$ = 0.25, $0\\leq g\\leq 4$', fontsize=7, color='grey')
ax.set_ylim(-0.5, 11); ax.set_xlabel('ring radius $r_0$ / wall radius')
ax.set_ylabel('central column strength $g$ (ring column = 1)')
ax.set_title('Trap ring + central column: solid = lower stability limit, dashed = $g^*$\n'
             '($H^2$ ring stable iff dashed $\\geq$ solid; crossings = $H^2$ thresholds)', fontsize=9)
ax.legend(fontsize=7, loc='center left'); fig.tight_layout(); fig.savefig('ring_center.png', dpi=150)
print('ok')
