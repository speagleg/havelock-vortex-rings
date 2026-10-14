"""Growth rates of hypercycle rows (tanh model) from the exact Bloch symbol, and their meaning in an expanding universe.
Symbol at wavenumber q (../rows/exact_symbol.py): [[i delta, B], [C, i delta]], B = 2 sum_j (a_j + b_j cos qj),
C = 2 sum_j d_j (1 - cos qj) < 0, so the growth rate is lambda(q) = sqrt(max(B C, 0)); lambda_max = max_q lambda(q).
Units: H = sum_{i<j} W, W = -ln tanh(d/2) = 2 pi G, i.e. circulation Gamma = 2 pi on the unit-curvature H^2.
Expanding universe: on a spatial slice a(t)^2 h_{H^2}, the stream function of point vortices is Gamma G_h (2D conformal
invariance) and the comoving velocity is a^{-2} times the unit-H^2 velocity, so the comoving dynamics is exactly the
H^2 problem in tau = int dt/a(t)^2. For Milne (flat Minkowski space, slices = unit hyperboloids, a = t) the remaining
tau after t0 is int dt/(c t)^2 = 1/(c^2 t0): a mode with rate lambda_max (for Gamma = 2 pi) grows at most by
exp(lambda_max Gamma/(2 pi c^2 t0)).
Usage: python growth_map.py -> growth_map.png, growth_map.out"""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from scipy.optimize import brentq

def pair(K, c):
    s2 = 1 - K; D = c + 1 - 2*K; den = (c - 1)*D**2
    a = s2*((3*K-1)*c**2 - (4*K**2-K+1)*c + 2*K**2)/den
    b = s2*(K*c**2 + (1-3*K)*c + 2*K**2-2*K+1)/den
    d = -s2**2*(c**2 + c - 2*K)/den
    return a, b, d
def lam_max(kappa, l, nq=181):
    K = kappa**2; lp = l*np.sqrt(1 - K)
    J = max(4, int(np.ceil(60/lp))); j = np.arange(1, J+1); c = np.cosh(np.minimum(j*lp, 700.0))
    a, b, d = pair(K, c)
    q = np.linspace(1e-3, np.pi, nq)
    cq = np.cos(np.outer(q, j))
    B = 2*(a.sum() + (cq*b).sum(axis=1)); C = 2*((1 - cq)*d).sum(axis=1)
    g = np.sqrt(np.maximum(B*C, 0)); i = np.argmax(g)
    return g[i], q[i]
def Fpi(kappa, lp):
    a = np.arccos(kappa); s2 = np.sin(a)**2; n = np.arange(1, int(np.ceil(90/lp)) + 6); r = np.exp(-n*lp)
    return np.sum((np.sin(2*a)*np.sin(2*n*a) - 2*n*s2*np.sin(n*a)**2)*r/(1-r) - 2*n*s2*np.cos(n*a)**2*r/(1+r))
def lstar(kappa):
    g = np.geomspace(0.005, 40, 600); v = np.array([Fpi(kappa, x) for x in g]); i = np.where(np.diff(np.sign(v)) > 0)[0][-1]
    return brentq(lambda x: Fpi(kappa, x), g[i], g[i+1])/np.sqrt(1 - kappa**2)

if __name__ == '__main__':
    ks = np.linspace(0.0, 0.995, 120); ls = np.geomspace(0.3, 10, 120)
    Lm = np.array([[lam_max(k, l)[0] for k in ks] for l in ls])
    with open('growth_map.out', 'w') as f:
        f.write('max growth rate lambda_max (Gamma = 2 pi, unit curvature) and critical q/pi\n')
        for k, l in [(0.0, 1.0), (0.0, 3.0), (0.5, 1.0), (0.5, 3.0), (0.7, 3.0), (0.8, 1.5), (0.9, 1.0), (0.95, 0.8)]:
            g, q = lam_max(k, l); f.write(f'kappa {k:4}, l {l:4}: lambda_max = {g:.6f}  at q/pi = {q/np.pi:.3f}\n')
        f.write('Milne freeze-out: total amplification after t0 is exp(lambda_max * Gamma/(2 pi c^2 t0));\n')
        f.write('superfluid 4He (Gamma = h/m4 = 9.97e-8 m^2/s): Gamma/c^2 = 1.1e-24 s, so for any t0 >> 1e-24 s the pattern is frozen\n')
    print(open('growth_map.out').read())
    INK, MUTED, GRID = '#2b2b29', '#6b6a64', '#e4e3df'
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    Z = np.ma.masked_where(Lm < 1e-9, np.log10(np.maximum(Lm, 1e-12)))     # stable points left white
    top = np.ceil(np.log10(Lm.max())*2)/2
    cs = ax.contourf(ks, ls, Z, levels=np.arange(-5, top + 0.01, 0.5), cmap='Blues', extend='min')
    cb = fig.colorbar(cs, ax=ax); cb.set_label('log10 of the largest growth rate (Gamma = 2 pi, curvature radius 1)', fontsize=8.5, color=INK)
    kk = np.linspace(0.709, 0.995, 40); ax.plot(kk, [lstar(k) for k in kk], '-', color='#eb6834', lw=2.4, label='stability boundary l*(kappa)')
    ax.axvline(1/np.sqrt(2), color=MUTED, ls='--', lw=1)
    ax.set_yscale('log'); ax.set_ylim(0.3, 10)
    ax.set_yticks([0.3, 1, 3, 10]); ax.set_yticklabels(['0.3', '1', '3', '10']); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel('kappa = geodesic curvature = proper acceleration of the matching AdS2 static observer', color=INK, fontsize=9)
    ax.set_ylabel('arc spacing l (curvature radii)', color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8.5, loc='upper left', labelcolor=INK)
    ax.set_title('Hypercycle rows: how fast the unstable ones fall apart (white = stable)', fontsize=10, color=INK, loc='left')
    ax.tick_params(colors=MUTED, labelsize=8.5)
    plt.tight_layout(); plt.savefig('growth_map.png', dpi=140)
