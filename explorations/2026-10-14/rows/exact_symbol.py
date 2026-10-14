"""Exact Bloch symbol of a hypercycle row (tanh model) and an explicit stability criterion.
Pair terms (symbol_q.py), K = kappa^2, s^2 = 1-K, c = cosh(j l'), D = c + 1 - 2K:
  a(c) = s^2 [(3K-1)c^2 - (4K^2-K+1)c + 2K^2] / ((c-1) D^2)
  b(c) = s^2 [K c^2 + (1-3K)c + 2K^2-2K+1]   / ((c-1) D^2)
  d(c) = -s^4 (c^2 + c - 2K) / ((c-1) D^2)            (< 0)
  diagonal entries: both equal 2i sum_j f_j sin(qj), f_j = -sqrt(K) s^3 sinh(j l')/D^2.
Symbol at Bloch wavenumber q: [[i delta, B],[C, i delta]], B = 2 sum_j (a_j + b_j cos qj), C = 2 sum_j d_j (1 - cos qj) < 0.
Exponents: i delta +- sqrt(B C).  So the row is spectrally stable  <=>  F(q) := sum_j a_j + b_j cos(qj) >= 0 for all q.
Usage: python exact_symbol.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings('ignore')
import numpy as np
from scipy.optimize import brentq

def ab(K, c):
    s2 = 1-K; D = c+1-2*K; den = (c-1)*D**2
    return (s2*((3*K-1)*c**2 - (4*K**2-K+1)*c + 2*K**2)/den, s2*(K*c**2 + (1-3*K)*c + 2*K**2-2*K+1)/den)

def Fmin(K, lp, nq=721):
    J = max(4, int(np.ceil(60/lp)))
    j = np.arange(1, J+1); cj = np.cosh(np.minimum(j*lp, 700.0))
    a, b = ab(K, cj)
    q = np.linspace(1e-3, np.pi, nq)
    F = a.sum() + (b[None, :]*np.cos(np.outer(q, j))).sum(axis=1)
    i = np.argmin(F); return F[i], q[i]

def lstar(kappa):
    """Arc-length critical spacing from the exact criterion (sign change of min_q F in l')."""
    K = kappa**2; s = np.sqrt(1-K)
    f = lambda lp: Fmin(K, lp)[0]
    lps = np.geomspace(0.02, 40, 400); v = [f(x) for x in lps]
    for i in range(len(lps)-1, 0, -1):
        if v[i-1] < 0 <= v[i]:
            return brentq(f, lps[i-1], lps[i], xtol=1e-12)/s
    return np.nan

if __name__ == '__main__':
    import rows
    # 1. symbol vs the full supercell linearisation (growth rates)
    print('1. relative growth max Re(lambda)/max|lambda|: exact symbol vs supercell (rows.stab_hyper, M = 40)')
    for kappa, l in [(0.5, 2.0), (0.7, 3.0), (0.8, 1.5), (0.95, 0.9), (0.9, 2.0)]:
        K = kappa**2; s = np.sqrt(1-K); lp = l*s
        J = int(np.ceil(60/lp)); j = np.arange(1, J+1); cj = np.cosh(j*lp); a, b = ab(K, cj)
        D = cj+1-2*K; d = -(1-K)**2*(cj**2+cj-2*K)/((cj-1)*D**2)
        f = -np.sqrt(K)*(1-K)**1.5*np.sinh(j*lp)/D**2
        lam = []
        for q in 2*np.pi*np.arange(0, 40)/40:
            B = 2*(a + b*np.cos(q*j)).sum(); C = 2*(d*(1-np.cos(q*j))).sum(); dl = 2j*(f*np.sin(q*j)).sum()
            r = np.sqrt(complex(B*C)); lam += [dl + r, dl - r]
        lam = np.array(lam)
        print(f'   kappa {kappa}, l {l}: symbol {lam.real.max()/abs(lam).max():.6e}   supercell {rows.stab_hyper(40, l, kappa)[0]:.6e}')
    # 2. thresholds from the exact criterion vs rows_boundary.txt
    print('2. critical arc spacing: exact criterion vs supercell bisection; critical q')
    for kappa, ls in np.loadtxt('rows_boundary.txt'):
        le = lstar(kappa); qc = Fmin(kappa**2, le*np.sqrt(1-kappa**2)*0.999)[1]
        print(f'   kappa {kappa:6.3f}: {le:.6f}  vs {ls:.6f}   q_crit/pi = {qc/np.pi:.3f}')
