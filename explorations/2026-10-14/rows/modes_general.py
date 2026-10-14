"""Hypercycle rows for the screened kernels W_nu = 2 pi G_nu, G_nu the Green's function of -Delta + nu(nu+1) on H^2
(nu = 0: tanh model, -ln tanh(d/2); nu = 1: Q_1(cosh d) = C artanh(1/C) - 1).
Mode expansion. In z = e^{u + i phi} the operator separates: -(d_u^2 + d_phi^2) + mu/sin^2(phi), mu = nu(nu+1). The
Poschl-Teller modes on (0, pi) are psi_n(phi) = sin^{nu+1}(phi) C^{(nu+1)}_{n-nu-1}(cos phi), eigenvalue n^2, n >= nu+1, so
  W(u, phi, phi0) = 2 pi sum_n e^{-n|u|}/(2 n N_n) psi_n(phi) psi_n(phi0),   N_n = int_0^pi psi_n^2.
With s = sin alpha, k = cos alpha, w_n(x) = 2 pi e^{-n x}/(2 n N_n), the row's pair functions are
  b = s^2 sum w_n psi_n'(alpha)^2 >= 0,   d = -s^2 sum n^2 w_n psi_n(alpha)^2 <= 0,
  a = sum w_n [2 s k psi_n psi_n' + (mu - s^2 n^2) psi_n^2].
Consequences for every nu: the critical Bloch mode is q = pi (b is a positive combination of e^{-n x}); the threshold is
  F(pi) = sum_n 2 pi/(2 n N_n) [ (2 s k psi psi' + (mu - s^2 n^2) psi^2) r_n/(1 - r_n) - s^2 psi'^2 r_n/(1 + r_n) ] = 0,
r_n = e^{-n l'}; and the lowest mode (psi ~ sin^{m} phi, m = nu+1) gives the leading term m s^{2m} cos(2 alpha) r_m (times
2 pi/(2 m N_m)): kappa_c = 1/sqrt 2 for every nu.
Usage: python modes_general.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings('ignore')
import numpy as np, sympy as sp
from scipy.special import eval_gegenbauer, gammaln
from scipy.optimize import brentq
from kernel_coeffs import funcs

def modes(nu, alpha, nmax):
    lam = nu + 1; n = np.arange(lam, nmax+1); j = n - lam; x = np.cos(alpha); s = np.sin(alpha)
    C = eval_gegenbauer(j, lam, x)
    dC = np.where(j > 0, 2*lam*eval_gegenbauer(np.maximum(j-1, 0), lam+1, x), 0.0)   # d/dx C_j^lam = 2 lam C_{j-1}^{lam+1}
    psi = s**lam*C
    dpsi = lam*s**(lam-1)*np.cos(alpha)*C - s**(lam+1)*dC                     # d/dphi, x = cos phi
    logN = np.log(np.pi) + (1-2*lam)*np.log(2) + gammaln(j+2*lam) - gammaln(j+1) - np.log(j+lam) - 2*gammaln(lam)
    return n, psi, dpsi, np.exp(logN)

def Fpi_modes(nu, kappa, lp, nmax=None):
    alpha = np.arccos(kappa); s, k = np.sin(alpha), np.cos(alpha); mu = nu*(nu+1)
    nmax = nmax or int(np.ceil(70/lp)) + nu + 5
    n, psi, dpsi, N = modes(nu, alpha, nmax); r = np.exp(-n*lp); wgt = 2*np.pi/(2*n*N)
    return np.sum(wgt*((2*s*k*psi*dpsi + (mu - s**2*n**2)*psi**2)*r/(1-r) - s**2*dpsi**2*r/(1+r)))

def Fq_direct(ker, kappa, lp, q):
    J = max(4, int(np.ceil(60/lp))); j = np.arange(1, J+1); z = np.exp(-j*lp)
    a, b, d = funcs(z, np.arccos(kappa), ker); return a.sum() + (b*np.cos(q*j)).sum()

if __name__ == '__main__':
    # 0. lowest-mode leading term
    al, m = sp.symbols('alpha m', positive=True)
    s_, k_ = sp.sin(al), sp.cos(al); psi = s_**m; dpsi = sp.diff(psi, al)
    lead = 2*s_*k_*psi*dpsi + (m*(m-1) - s_**2*m**2)*psi**2 - s_**2*dpsi**2
    S, Kc = sp.symbols('S Kc', positive=True)       # S = sin a, Kc = cos a
    ex = sp.expand(sp.powsimp(sp.expand(lead/(m*s_**(2*m))), force=True)).subs({sp.sin(al): S, sp.cos(al): Kc})
    ex = sp.simplify(sp.expand(ex).subs(Kc**2, 1 - S**2))
    print('0. lowest-mode coefficient / (m sin^{2m} a) =', ex, ' (= 1 - 2 sin^2 a = cos 2a)')
    # 1. mode sums vs direct lattice sums (pair functions from kernel_coeffs.funcs)
    worst = 0
    for nu, ker in [(0, 'Q0'), (1, 'Q1')]:
        for kappa in [0.4, 0.75, 0.95]:
            for lp in [0.2, 0.8, 2.5]:
                x, y = Fpi_modes(nu, kappa, lp), Fq_direct(ker, kappa, lp, np.pi); worst = max(worst, abs(x-y)/abs(y))
    print(f'1. F(pi): mode sum vs direct lattice sum, max relative difference {worst:.1e}')
    # 2. q = pi is critical for Q_1: F(q) decreasing on a grid
    bad = 0
    for kappa in np.linspace(0.1, 0.995, 10):
        for lp in np.geomspace(0.05, 6, 10):
            F = [Fq_direct('Q1', kappa, lp, q) for q in np.linspace(0.02, np.pi, 80)]
            bad += np.sum(np.diff(F) > 1e-12*np.abs(F).max())
    print(f'2. Q_1: grid points where F(q) increases: {bad} (of 10 x 10 x 79)')
    # 3. Q_1 thresholds (arc length) from F(pi) = 0, vs the supercell values of rows_kernel.py
    print('3. Q_1 critical spacing l*(kappa) = root of F(pi):  [supercell, rows_kernel.py]')
    ref = {0.75: 4.09, 0.8: 3.01, 0.9: 1.98, 0.95: 1.69, 0.99: 1.49}
    for kappa in [0.72, 0.75, 0.8, 0.9, 0.95, 0.99, 0.999, 0.9999]:
        s = np.sqrt(1-kappa**2); g = np.geomspace(0.01, 40, 500); v = np.array([Fpi_modes(1, kappa, x) for x in g])
        ch = np.where(np.diff(np.sign(v)) != 0)[0]
        root = brentq(lambda x: Fpi_modes(1, kappa, x), g[ch[-1]], g[ch[-1]+1], xtol=1e-14)/s
        print(f'   kappa {kappa:7.4f}: l* = {root:.6f}   sign changes {len(ch)}   ' + (f'[{ref[kappa]}]' if kappa in ref else ''))
    # 4. kappa_c for Q_1 and nu = 2 (large spacing): sign of F(pi) just below / above 1/sqrt 2 at l' = 30
    for nu in [1, 2]:
        print(f'4. nu = {nu}: F(pi) at l\'=30, kappa = 0.7070 / 0.7072:', Fpi_modes(nu, 0.7070, 30.0), Fpi_modes(nu, 0.7072, 30.0))
    # 5. nu = 1, 2: F(pi) < 0 for kappa <= 1/sqrt 2 at every spacing tested
    for nu in [1, 2]:
        npos = 0; tot = 0
        for kappa in np.linspace(0.05, 1/np.sqrt(2), 25):
            for lp in np.geomspace(0.02, 25, 120):
                npos += Fpi_modes(nu, kappa, lp) >= 0; tot += 1
        print(f'5. nu = {nu}: points with F(pi) >= 0 for kappa <= 1/sqrt2, l\' in [0.02, 25]: {npos} of {tot}')
    # 6. horocycle limit for nu = 1 (alpha -> 0, n alpha -> t, l = l'/alpha):
    #    I_1(l) = int 2/t^3 { [2(tc - s)(s - t^2 s - tc) + (2 - t^2)(tc - s)^2]/(e^{lt} - 1) - (s - t^2 s - tc)^2/(e^{lt} + 1) } dt
    from scipy.integrate import quad
    def I1(l):
        def f(t):
            c, s = np.cos(t), np.sin(t); A = t*c - s; B = s - t**2*s - t*c
            if t < 1e-3:   # series: A ~ -t^3/3, B ~ -t^3/2 ... use the exact expression anyway with care
                pass
            return 2/t**3*((2*A*B + (2 - t**2)*A**2)/np.expm1(l*t) - B**2/(np.exp(l*t) + 1))
        return quad(f, 1e-6, 80/l, limit=500)[0]
    lh = brentq(I1, 0.8, 3.0, xtol=1e-12)
    print(f'6. nu = 1 horocycle row (limit kappa -> 1): critical arc spacing l* = {lh:.8f}  (cf. kappa = 0.9999 above)')
