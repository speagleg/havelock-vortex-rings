"""Landau theory of the zigzag transition of a hypercycle row (tanh model), from the exact lattice energy.
Strip coordinates, canonical pair (u, Y = -cot phi). Kernel (exact): W = (1/2) ln[(cosh du - cos(p1+p2))/(cosh du - cos(p1-p2))].
Zigzag cell: A at u = 2 m l', B at u = (2m+1) l' (no relative shift: forbidden at this order by the reflection u -> -u),
Y_A = Y0 + eta, Y_B = Y0 - eta. Energy per cell
    e(Y_A, Y_B; l') = sum_{m>=1} [W_AA(2 m l') + W_BB(2 m l')] + sum_{m in Z} W_AB((2m+1) l').
Relative equilibria at fixed (l', Y0): de/d eta = 0 (equal u-velocities). Landau function f(eta) = e - e(0) = a2 eta^2 + a4 eta^4.
Relaxed quartic: let (l', Y0) relax at fixed conjugate forces (P, V) of the straight row:
    a4_rel = a4 - (1/2) g^T M^{-1} g,  g = d a2/d(l', Y0),  M = Hessian of e(0) - V (2 Y0) + 2 P l' in (l', Y0)
(M must be positive definite: the long-wavelength sector of the row is stable).
Usage: python landau.py -> landau.out"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows'))
import numpy as np, mpmath as mp
from symbol import symbol
from exact_symbol import lstar
mp.mp.dps = 30

def W(du, YA, YB):
    pA, pB = mp.acot(-YA) if YA < 0 else mp.pi - mp.acot(YA), mp.acot(-YB) if YB < 0 else mp.pi - mp.acot(YB)
    ch = mp.cosh(du)
    return mp.log((ch - mp.cos(pA + pB))/(ch - mp.cos(pA - pB)))/2

def phi_of(Y):            # Y = -cot(phi), phi in (0, pi)
    return mp.pi/2 + mp.atan(Y)

def Wp(du, pA, pB):
    ch = mp.cosh(du); return mp.log((ch - mp.cos(pA + pB))/(ch - mp.cos(pA - pB)))/2

def e_cell(YA, YB, lp):
    pA, pB = phi_of(YA), phi_of(YB)
    M = int(mp.ceil(80/lp)) + 2; tot = mp.mpf(0)
    for m in range(1, M):
        tot += Wp(2*m*lp, pA, pA) + Wp(2*m*lp, pB, pB)
    for m in range(-M, M):
        tot += Wp((2*m + 1)*lp, pA, pB)
    return tot

def coeffs(lp, Y0):
    f = lambda eta: e_cell(Y0 + eta, Y0 - eta, lp)
    a2 = mp.diff(f, 0, 2)/2; a4 = mp.diff(f, 0, 4)/24
    return a2, a4

def relaxed(lp, Y0):
    lp, Y0 = mp.mpf(lp), mp.mpf(Y0)
    a2, a4 = coeffs(lp, Y0)
    # couplings of eta^2 to the uniform modes
    g_l = mp.diff(lambda x: coeffs(x, Y0)[0], lp); g_Y = mp.diff(lambda y: coeffs(lp, y)[0], Y0)
    e0 = lambda x, y: e_cell(y, y, x)
    # Hessian of e0 - V*2Y0 + 2P*l' : the linear terms drop out of the Hessian
    M = mp.matrix([[mp.diff(e0, (lp, Y0), (2, 0)), mp.diff(e0, (lp, Y0), (1, 1))],
                   [mp.diff(e0, (lp, Y0), (1, 1)), mp.diff(e0, (lp, Y0), (0, 2))]])
    g = mp.matrix([g_l, g_Y])
    corr = (g.T*mp.inverse(M)*g)[0]/2
    return a2, a4, a4 - corr, M

if __name__ == '__main__':
    out = []; P = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
    P('1. Check: a2 is proportional to B(pi) (the q = pi transverse stiffness), so it vanishes at l*:')
    for kappa, l in [(0.8, 2.0), (0.8, 2.6), (0.9, 1.3), (0.9, 1.7)]:
        s = np.sqrt(1 - kappa**2); lp = l*s; Y0 = -kappa/s
        a2, a4 = coeffs(mp.mpf(lp), mp.mpf(Y0)); B, C, dl = symbol(kappa, lp, np.array([np.pi]))
        P(f'   kappa {kappa}, l {l}: a2 = {mp.nstr(a2, 10)},  B(pi) s^2 = {B[0]*s*s:.10f},  ratio a2/(B s^2) = {mp.nstr(a2/(B[0]*s*s), 10)}')
    P('\n2. At the threshold l*(kappa): quartic coefficients (fixed l\', Y0) and (relaxed at fixed forces)')
    P('   kappa     l*        a2(l*)        a4_fixed       a4_relaxed     M positive definite')
    for kappa in [0.72, 0.75, 0.8, 0.85, 0.9, 0.95, 0.98, 0.995]:
        ls = lstar(kappa); s = np.sqrt(1 - kappa**2)
        a2, a4, a4r, M = relaxed(ls*s, -kappa/s)
        pd = M[0, 0] > 0 and M[0, 0]*M[1, 1] - M[0, 1]**2 > 0
        P(f'   {kappa:5.3f}  {ls:8.5f}   {mp.nstr(a2, 3):>12s}   {mp.nstr(a4, 8):>12s}   {mp.nstr(a4r, 8):>12s}     {pd}')
    open('landau.out', 'w').write('\n'.join(out) + '\n')
