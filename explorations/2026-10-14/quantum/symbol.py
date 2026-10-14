"""Exact Bloch symbol of a hypercycle row (tanh model), from ../rows/exact_symbol.py, vectorised over q.
Strip coordinates (u, phi), dynamics u' = sin^2(phi) dH/dphi, phi' = -sin^2(phi) dH/du, H = sum W (circulation 2 pi).
Symbol [[i delta, B], [C, i delta]] acting on (du, dphi)."""
import numpy as np
def pieces(kappa, lp):
    K = kappa**2; s2 = 1 - K
    J = max(6, int(np.ceil(80/lp))); j = np.arange(1, J + 1); c = np.cosh(np.minimum(j*lp, 700.0))
    D = c + 1 - 2*K; den = (c - 1)*D**2
    a = s2*((3*K - 1)*c**2 - (4*K**2 - K + 1)*c + 2*K**2)/den
    b = s2*(K*c**2 + (1 - 3*K)*c + 2*K**2 - 2*K + 1)/den
    d = -s2**2*(c**2 + c - 2*K)/den
    f = -np.sqrt(K)*s2**1.5*np.sinh(np.minimum(j*lp, 700.0))/D**2
    return j, a, b, d, f
def symbol(kappa, lp, q):
    j, a, b, d, f = pieces(kappa, lp); q = np.atleast_1d(q)
    cq, sq = np.cos(np.outer(q, j)), np.sin(np.outer(q, j))
    B = 2*(a.sum() + cq @ b); C = 2*((1 - cq) @ d); delta = 2*(sq @ f)
    return B, C, delta
def derivs_at_pi(kappa, lp):
    """B''(pi), C(pi), delta'(pi) analytically from the lattice sums."""
    j, a, b, d, f = pieces(kappa, lp); sgn = (-1.0)**j
    return 2*(-(j**2)*b*sgn).sum(), 2*(d*(1 - sgn)).sum(), 2*(j*f*sgn).sum()
def derivs_at_0(kappa, lp):
    """B(0), C''(0)/2 (C ~ c2 q^2), delta'(0)."""
    j, a, b, d, f = pieces(kappa, lp)
    return 2*(a + b).sum(), (d*j**2).sum(), 2*(j*f).sum()
