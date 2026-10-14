"""Hypercycle rows (tanh model): proofs of an explicit stability bound and of all-mode instability for kappa <= 1/2.
Uses the exact symbol of exact_symbol.py: stable <=> F(q) = sum_{j>=1} a(c_j) + b(c_j) cos(qj) >= 0 for all q, c_j = cosh(j l').
(I)   b(c) > 0 for c > 1 and every K in (0,1).
(II)  a - b = s^2 ((2K-1)c - 1)/((c-1)D),  a + b = s^2 ((4K-1)c - (4K^2-2K+1))/D^2,  D = c + 1 - 2K.
Theorem A (stability): if kappa > 1/sqrt 2 and cosh l' >= 1/(2 kappa^2 - 1), then F(q) >= sum_j (a_j - b_j) > 0 for all q.
Theorem B (all-mode instability): if kappa <= 1/2, then F(q) <= sum_j (a_j + b_j) < 0 for all q, so every Bloch mode
  q != 0 grows; in particular every necklace of N >= 2 vortices (either parity) is unstable at every spacing.
Usage: python stability_bound.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings('ignore')
import numpy as np, sympy as sp
from exact_symbol import ab, lstar

c, K, t = sp.symbols('c K t', positive=True)
s2 = 1-K; D = c+1-2*K; den = (c-1)*D**2
a = s2*((3*K-1)*c**2 - (4*K**2-K+1)*c + 2*K**2)/den
b = s2*(K*c**2 + (1-3*K)*c + 2*K**2-2*K+1)/den
nb = sp.expand(K*(1+t)**2 + (1-3*K)*(1+t) + 2*K**2-2*K+1)
print('(I)  numerator of b in t = c - 1:', sp.Poly(nb, t).all_coeffs(), '-> coefficients K, 1-K, 2(1-K)^2, all > 0')
print('(II) a - b =', sp.factor(sp.simplify(a-b)))
print('     a + b =', sp.factor(sp.simplify(a+b)))
print('     check vs stated forms:', sp.simplify(a-b - s2*((2*K-1)*c-1)/((c-1)*D)), sp.simplify(a+b - s2*((4*K-1)*c-(4*K**2-2*K+1))/D**2))

print('\nTheorem A: sufficient arc spacing l_A = arccosh(1/(2k^2-1))/sqrt(1-k^2) vs exact threshold l*')
for kappa in [0.72, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999]:
    k2 = kappa**2; lA = np.arccosh(1/(2*k2-1))/np.sqrt(1-k2); ls = lstar(kappa)
    print(f'   kappa {kappa:5.3f}: l_A = {lA:.4f}   l* = {ls:.4f}   ratio {lA/ls:.3f}')
print('   kappa -> 1 limit of l_A: 2 (horocycle exact value 0.8439);  kappa -> 1/sqrt2: l_A - l* -> ln2/sqrt(1-k^2) = ', np.log(2)*np.sqrt(2))

print('\nTheorem B check: max_q F over spacings, kappa <= 1/2 (must be < 0)')
worst = -np.inf
for kappa in [0.0, 0.2, 0.4142, 0.5]:
    K_ = kappa**2
    for lp in np.geomspace(0.01, 30, 200):
        J = max(4, int(np.ceil(60/lp))); j = np.arange(1, J+1); cj = np.cosh(np.minimum(j*lp, 700.0)); A, B = ab(K_, cj)
        q = np.linspace(0, np.pi, 361); F = A.sum() + (B[None, :]*np.cos(np.outer(q, j))).sum(axis=1)
        worst = max(worst, F.max()/abs(A.sum()))
print(f'   max over kappa, l\', q of F/|sum a| = {worst:.3e}')

print('\nBolza systole necklaces (cylinder model): Theorem A gives a sufficient distance h_A from the geodesic')
L = 2*np.arccosh(1+np.sqrt(2))
for N in range(2, 13):
    cl = np.cosh(L/N); kA2 = (1+1/cl)/2; hA = np.arctanh(np.sqrt(kA2))
    print(f'   N = {N:2d}: stable for all h >= h_A = {hA:.4f}  (kappa_A = {np.sqrt(kA2):.4f})')
