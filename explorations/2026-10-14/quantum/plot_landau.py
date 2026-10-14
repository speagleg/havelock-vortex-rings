"""Figure: quartic Landau coefficient a4 of the zigzag transition along the threshold curve l*(kappa), fixed and relaxed
ensembles, with the tricritical points (landau.out, tricritical.out). Usage: python plot_landau.py -> landau.png"""
import numpy as np, re, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
rows = [l.split() for l in open('landau.out') if re.match(r'^\s+0\.\d+\s+\d', l)]
k = np.array([float(r[0]) for r in rows]); af = np.array([float(r[3]) for r in rows]); ar = np.array([float(r[4]) for r in rows])
kt = [float(re.search(r'kappa_t = ([\d.]+)', l).group(1)) for l in open('tricritical.out')]
plt.rcParams.update({'font.size': 9})
fig, ax = plt.subplots(figsize=(6.4, 3.9))
ax.plot(k, af, color='#2a78d6', lw=2, marker='o', ms=4, label='fixed spacing and impulse')
ax.plot(k, ar, color='#eb6834', lw=2, marker='o', ms=4, label='relaxed (fixed forces)')
ax.axhline(0, color='#6b6a64', lw=0.8)
for x, c, name in zip(kt, ['#2a78d6', '#eb6834'], ['fixed', 'relaxed']):
    ax.axvline(x, color=c, lw=1, ls=':'); ax.annotate(f'$\\kappa_t$ = {x:.3f} ({name})', (x, -0.06 if name == 'fixed' else -0.08),
                                                      xytext=(4, 0), textcoords='offset points', fontsize=7, color='#1f1f1e')
ax.text(0.905, 0.012, '$a_4>0$: continuous (Ising-type)', fontsize=7, color='#6b6a64')
ax.text(0.73, -0.068, '$a_4<0$: first order (subcritical)', fontsize=7, color='#6b6a64')
ax.set_xlabel('geodesic curvature $\\kappa$ of the hypercycle (on the threshold $\\ell = \\ell^*(\\kappa)$)')
ax.set_ylabel('quartic coefficient $a_4$ (per cell)')
ax.set_title('Zigzag transition: the Landau quartic changes sign at a tricritical point', fontsize=9)
ax.spines[['top', 'right']].set_visible(False); ax.grid(color='#e6e5e0', lw=0.6); ax.legend(fontsize=7, frameon=False, loc='lower right')
fig.tight_layout(); fig.savefig('landau.png', dpi=150); print('ok')
