"""Exact Fourier blocks for a ring of N in the unit disk, kernel W = ln|1 - z conj(w)| - ln|z - w| (the H^2 tanh kernel),
any radial mobility and self term, via the lemma of flat_lemma.py:
    J = mu [[ -H_tr , -H_tt ], [ H_rr + S I , H_rt ]],  S = U'' + (mu'/mu)(W_r + U').
Second derivatives of W are analytic (sympy, lambdified), so the k = 0 zero modes are exact.
Mode k block: h_xy(k) = sum_j H_{x_0 y_j} e^{i k theta_j}. Stability of mode k: both eigenvalues of the 2x2 block imaginary.
Usage: imported by shear_design.py; python ring_blocks.py runs the checks vs flat_lemma (finite differences)."""
import numpy as np, sympy as sp
r1, r2, t = sp.symbols('r1 r2 t', positive=True)
Wsym = sp.log(1 - 2*r1*r2*sp.cos(t) + r1**2*r2**2)/2 - sp.log(r1**2 + r2**2 - 2*r1*r2*sp.cos(t))/2   # t = th1 - th2
d = lambda e, *v: sp.lambdify((r1, r2, t), sp.diff(e, *v), 'numpy')
W_r1 = d(Wsym, r1)
W_r1r1, W_r1r2 = d(Wsym, r1, r1), d(Wsym, r1, r2)
W_r1t, W_r2t, W_tt = d(Wsym, r1, t), d(Wsym, r2, t), d(Wsym, t, t)

def blocks(N, r0):
    """returns W_r (radial force on vortex 0) and, for k = 0..N-1, (h_rr, h_tt, h_rt, h_tr) of the interaction alone."""
    th = 2*np.pi*np.arange(1, N)/N; tt = -th   # t = th_0 - th_j
    Wr = W_r1(r0, r0, tt).sum()
    out = []
    for k in range(N):
        ph = np.exp(1j*k*th)
        # diagonal (j = 0) entries: sums over partners; off-diagonal: pair second derivatives with the partner's variable
        hrr = W_r1r1(r0, r0, tt).sum() + (W_r1r2(r0, r0, tt)*ph).sum()
        htt = W_tt(r0, r0, tt).sum() - (W_tt(r0, r0, tt)*ph).sum()          # d/dth_j = -d/dt
        hrt = W_r1t(r0, r0, tt).sum() - (W_r1t(r0, r0, tt)*ph).sum()        # d^2/dr_0 dth_j = -W_r1t
        htr = W_r1t(r0, r0, tt).sum() + (W_r2t(r0, r0, tt)*ph).sum()        # d^2/dth_0 dr_j = +W_r2t
        out.append((hrr, htt, hrt, htr))
    return Wr, out

def ring_spectrum(N, r0, S):
    Wr, bl = blocks(N, r0); ev = []
    for hrr, htt, hrt, htr in bl:
        ev.extend(np.linalg.eigvals(np.array([[-htr, -htt], [hrr + S, hrt]])))
    return np.array(ev)

def S_trap(N, r0, s):          # m = 1, U = (1/2) ln(1 - r^2) + U_b with r (U_b'/r)' = s
    Wr, _ = blocks(N, r0); r = r0
    Rp, Rpp = -r/(1 - r*r), -(1 + r*r)/(1 - r*r)**2
    return Rpp + s - (Rp + Wr)/r            # note r(U_b'/r)' = U_b'' - U_b'/r, so U_b' drops out except via this

def S_H2(N, r0):
    Wr, _ = blocks(N, r0); r = r0
    return -(4*r/(1 - r*r) + 1/r)*Wr

def growth(N, r0, S):
    return np.abs(ring_spectrum(N, r0, S).real).max()

if __name__ == '__main__':
    from flat_lemma import ring_eigs, match
    one = lambda r: (np.ones_like(r), 0*r)
    for N, r0 in [(3, 0.3), (6, 0.5), (8, 0.7)]:
        e_fd, _ = ring_eigs(N, r0, one, lambda r: -r/(1 - r*r))
        e_ex = ring_spectrum(N, r0, S_trap(N, r0, 0.0))/r0   # mu = m/r = 1/r0
        print(f'N {N} r0 {r0}: exact blocks vs finite-difference trap spectrum: {match(e_ex, e_fd):.1e}')

def S_window(N, r0):
    """exact stability window in S: each block [[a, b], [c + S, d]] has eigenvalues (a+d)/2 +- sqrt(((a-d)/2)^2 + b(c+S));
    with a + d and a - d imaginary, mode k is stable iff ((a-d)/2)^2 + b(c + S) <= 0. Returns (S_lo, S_hi)."""
    Wr, bl = blocks(N, r0); lo, hi = -np.inf, np.inf
    for k, (hrr, htt, hrt, htr) in enumerate(bl):
        a, b, c, dd = -htr, -htt, hrr, hrt
        assert abs((a + dd).real) < 1e-9*(1 + abs(a) + abs(dd)) and abs((a - dd).real) < 1e-9*(1 + abs(a) + abs(dd))
        h2 = (((a - dd)/2)**2).real; b = b.real; c = c.real
        if abs(b) < 1e-12:
            continue   # k = 0: b = 0, eigenvalues (a+d)/2 +- (a-d)/2 imaginary
        Sk = -c - h2/b
        if b < 0: lo = max(lo, Sk)
        else: hi = min(hi, Sk)
    return lo, hi
