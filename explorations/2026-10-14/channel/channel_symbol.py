"""Exact Bloch symbols of a vortex row in a flat channel (0 < y < pi) and on an H^2 hypercycle, from the same pair kernel.
Per pair (x = j l', c = cosh x, K = cos^2 a, s = sin a, k = cos a), the tanh/H^2 pair functions a, b, d of ../rows/ and the
one-sided derivative W_phi = s k/(c - cos 2a) (from a = d + 2 k s W_phi, ../lorentzian/correlator_identity.py) give the
channel pair functions a_f = (a - 2 s k W_phi)/s^2, b_f = b/s^2, d_f = d/s^2, plus the Robin self term (1/2) ln sin(phi):
   B_ch(q) = 2 sum (a_f + b_f cos qj) - csc^2(a)/2,   C_ch(q) = 2 sum d_f (1 - cos qj) < 0,
   V_ch = 2 sum W_phi + cot(a)/2  (row speed).
Exact identity:  B_H2(q) = sin^2(a) B_ch(q) + sin(2a) V_ch - cos(2a)/2,   C_H2(q) = sin^2(a) C_ch(q).
Stability: channel row stable iff B_ch(q) >= 0 for all q; H^2 row iff B_H2(q) >= 0 for all q.
Usage: python channel_symbol.py"""
import numpy as np
def pieces(alpha, lp):
    K = np.cos(alpha)**2; s, k = np.sin(alpha), np.cos(alpha)
    J = max(4, int(np.ceil(60/lp))); j = np.arange(1, J+1); c = np.cosh(np.minimum(j*lp, 700.0))
    s2 = 1 - K; D = c + 1 - 2*K; den = (c - 1)*D**2
    a = s2*((3*K-1)*c**2 - (4*K**2-K+1)*c + 2*K**2)/den
    b = s2*(K*c**2 + (1-3*K)*c + 2*K**2-2*K+1)/den
    d = -s2**2*(c**2 + c - 2*K)/den
    Wp = s*k/(c - np.cos(2*alpha))
    return j, a, b, d, Wp
def symbols(alpha, lp, q):
    j, a, b, d, Wp = pieces(alpha, lp); s, k = np.sin(alpha), np.cos(alpha)
    cq = np.cos(np.outer(q, j))
    BH = 2*(a.sum() + (cq*b).sum(axis=1)); CH = 2*((1 - cq)*d).sum(axis=1)
    af, bf, df = (a - 2*s*k*Wp)/s**2, b/s**2, d/s**2
    Bc = 2*(af.sum() + (cq*bf).sum(axis=1)) - 0.5/s**2; Cc = 2*((1 - cq)*df).sum(axis=1)
    V = 2*Wp.sum() + 0.5*k/s
    return BH, CH, Bc, Cc, V
if __name__ == '__main__':
    q = np.linspace(0.01, np.pi, 200)
    worst = 0
    for al in [0.2, 0.6, 1.0, 1.4]:
        for lp in [0.3, 1.0, 3.0]:
            BH, CH, Bc, Cc, V = symbols(al, lp, q)
            worst = max(worst, np.abs(BH - (np.sin(al)**2*Bc + np.sin(2*al)*V - np.cos(2*al)/2)).max(), np.abs(CH - np.sin(al)**2*Cc).max())
    print(f'identity B_H2 = s^2 B_ch + sin(2a) V_ch - cos(2a)/2, C_H2 = s^2 C_ch: max deviation {worst:.1e}')
    nst = 0; tot = 0; maxB = -np.inf
    for al in np.linspace(0.02, np.pi/2, 160):
        for lp in np.geomspace(0.02, 30, 160):
            BH, CH, Bc, Cc, V = symbols(al, lp, q)
            tot += 1; nst += np.all(Bc*Cc <= 1e-12*np.abs(Bc*Cc).max()); maxB = max(maxB, Bc.min())
    print(f'channel: stable at {nst} of {tot} grid points (alpha in (0, pi/2], l\' in [0.02, 30]); max over grid of min_q B_ch = {maxB:.3e}')
    # the destabilising mechanism: B_ch(pi) -> -csc^2(a)/2 at large spacing (wall image self-term), vs B_H2(pi) -> 0
    for al in [0.3, 0.785, 1.3]:
        BH, CH, Bc, Cc, V = symbols(al, 12.0, np.array([np.pi]))
        print(f'  alpha {al}: l\' = 12: B_ch(pi) = {Bc[0]:.6f} (-csc^2/2 = {-0.5/np.sin(al)**2:.6f}),  B_H2(pi) = {BH[0]:.2e}')
