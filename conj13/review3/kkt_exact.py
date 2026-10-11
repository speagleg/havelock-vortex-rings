"""Third reader, route B: order-by-order KKT expansion of the level-set minimisation in CARTESIAN
coordinates, with no Fourier ansatz, no mode selection and no residue sums.  All 2n coordinates are
unknowns; the field is either GF(p) (exact; p = 1 mod 4n so zeta_n, i exist; sqrt(n-3) or sqrt(7)
exists) or mpmath (to select the root b_n in (-1,1) and cross-check).

f(eps) = H(q(eps)) where q(eps) is the KKT point of
  H = -sum_{i<j} log|zi-zj|^2 + (n-1) sum log sigma_k,  sigma = 1+lam|z|^2
  s.t.  J1 = sum x/sigma = 0, J2 = sum y/sigma = 0, J3 = sum |z|^2/sigma = n/(1+lam),
        t.(q-q0) = 0 (SO(2) slice), e.(q-q0) = eps, e2.(q-q0) = 0 (odd n),
e_k = cos(l phi_k)*(radial unit) (unnormalised), e2 = sin(l phi_k)*(radial unit), l = floor(n/2).
Envelope: f'(eps) = nu_e(eps), so f_4 = nu_e[3]/4; c4 (unit-norm e) = f_4 * |e|^4 = f_4 (n/2)^2 (odd),
n^2 f_4 (even, |e|^2 = n).  Also checks nu_e[1] = 0 (threshold) and nu_e[2] = 0 (no cubic).
"""
import sys, mpmath as mp
from fractions import Fraction

class Fp:
    p = None
    __slots__ = ('v',)
    def __init__(s, v): s.v = v % Fp.p
    def _c(o): return o if isinstance(o, Fp) else Fp(o)
    def __add__(s, o): return Fp(s.v + Fp._c(o).v)
    __radd__ = __add__
    def __sub__(s, o): return Fp(s.v - Fp._c(o).v)
    def __rsub__(s, o): return Fp(Fp._c(o).v - s.v)
    def __mul__(s, o): return Fp(s.v * Fp._c(o).v)
    __rmul__ = __mul__
    def __neg__(s): return Fp(-s.v)
    def inv(s):
        assert s.v != 0, 'division by zero mod p'
        return Fp(pow(s.v, Fp.p - 2, Fp.p))
    def __truediv__(s, o): return s * Fp._c(o).inv()
    def __rtruediv__(s, o): return Fp._c(o) * s.inv()
    def iszero(s): return s.v == 0
    def __pow__(s, k): return Fp(pow(s.v, k % (Fp.p - 1), Fp.p)) if s.v else Fp(0)
    def __repr__(s): return 'Fp(%d)' % s.v

# ---------- truncated power series in eps (length K) ----------
K = 5
def S(c0): return [c0] + [c0 * 0] * (K - 1)
def sadd(a, b): return [x + y for x, y in zip(a, b)]
def ssub(a, b): return [x - y for x, y in zip(a, b)]
def smul(a, b):
    out = [a[0] * 0] * K
    for i in range(K):
        if is0(a[i]): continue
        for j in range(K - i): out[i + j] = out[i + j] + a[i] * b[j]
    return out
def sinv(a):
    b = [a[0] * 0] * K; b[0] = 1 / a[0]
    for k in range(1, K):
        acc = a[0] * 0
        for j in range(1, k + 1): acc = acc + a[j] * b[k - j]
        b[k] = -acc * b[0]
    return b
def sscal(c, a): return [c * x for x in a]
def is0(x): return x.iszero() if isinstance(x, Fp) else (abs(x) < mp.mpf(10) ** (-(mp.mp.dps - 10)))

def solve(M, b):
    n = len(M); A = [row[:] + [bb] for row, bb in zip(M, b)]
    for c in range(n):
        piv = None
        best = None
        for r in range(c, n):
            if isinstance(A[r][c], Fp):
                if not A[r][c].iszero(): piv = r; break
            else:
                if best is None or abs(A[r][c]) > best: best = abs(A[r][c]); piv = r
        assert piv is not None and (not (best is not None) or best > mp.mpf(10) ** (-mp.mp.dps // 2)), 'singular'
        A[c], A[piv] = A[piv], A[c]
        iv = 1 / A[c][c]
        A[c] = [x * iv for x in A[c]]
        for r in range(n):
            if r != c and not is0(A[r][c]):
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [A[r][n] for r in range(n)]

def run(n, lam, cs, sn, cl, sl, one):
    """cs[k], sn[k] = cos, sin phi_k ; cl[k], sl[k] = cos, sin (l phi_k); 'one' = field unit."""
    zero = one * 0
    odd = n % 2 == 1
    N = 2 * n
    q0 = []
    for k in range(n): q0 += [cs[k], sn[k]]
    e = []; e2 = []; t = []
    for k in range(n):
        e += [cl[k] * cs[k], cl[k] * sn[k]]; e2 += [sl[k] * cs[k], sl[k] * sn[k]]; t += [-sn[k], cs[k]]
    lin = [t, e] + ([e2] if odd else [])
    nc = 3 + len(lin)            # constraint order: J1 J2 J3 t e (e2)
    mu3 = n * one / (1 + lam)
    def E(qs, nus):
        """KKT residual as series: qs = list of N series, nus = list of nc series."""
        g = [S(zero) for _ in range(N)]
        gc = [[S(zero) for _ in range(N)] for _ in range(3)]
        cv = [S(zero) for _ in range(3)]
        for i in range(n):
            xi, yi = qs[2 * i], qs[2 * i + 1]
            for j in range(i + 1, n):
                dx = ssub(xi, qs[2 * j]); dy = ssub(yi, qs[2 * j + 1])
                d2i = sinv(sadd(smul(dx, dx), smul(dy, dy)))
                gx = sscal(-2 * one, smul(dx, d2i)); gy = sscal(-2 * one, smul(dy, d2i))
                g[2 * i] = sadd(g[2 * i], gx); g[2 * i + 1] = sadd(g[2 * i + 1], gy)
                g[2 * j] = ssub(g[2 * j], gx); g[2 * j + 1] = ssub(g[2 * j + 1], gy)
            u = sadd(smul(xi, xi), smul(yi, yi))
            sg = sadd(S(one), sscal(lam, u)); si = sinv(sg); si2 = smul(si, si)
            g[2 * i] = sadd(g[2 * i], sscal(2 * (n - 1) * lam, smul(xi, si)))
            g[2 * i + 1] = sadd(g[2 * i + 1], sscal(2 * (n - 1) * lam, smul(yi, si)))
            cv[0] = sadd(cv[0], smul(xi, si)); cv[1] = sadd(cv[1], smul(yi, si)); cv[2] = sadd(cv[2], smul(u, si))
            # grads: d(x/s)/dx = 1/s - 2 lam x^2/s^2 ; d(x/s)/dy = -2 lam x y/s^2 ; d(u/s) = 2(x,y)/s^2
            xy = smul(xi, yi)
            gc[0][2 * i] = ssub(si, sscal(2 * lam, smul(smul(xi, xi), si2)))
            gc[0][2 * i + 1] = sscal(-2 * lam, smul(xy, si2))
            gc[1][2 * i] = sscal(-2 * lam, smul(xy, si2))
            gc[1][2 * i + 1] = ssub(si, sscal(2 * lam, smul(smul(yi, yi), si2)))
            gc[2][2 * i] = sscal(2 * one, smul(xi, si2)); gc[2][2 * i + 1] = sscal(2 * one, smul(yi, si2))
        res = []
        for a in range(N):
            r = g[a]
            for c in range(3): r = ssub(r, smul(nus[c], gc[c][a]))
            for c, v in enumerate(lin): r = ssub(r, sscal(v[a], nus[3 + c]))
            res.append(r)
        targets = [zero, zero, mu3]
        for c in range(3): res.append(ssub(cv[c], S(targets[c])))
        for c, v in enumerate(lin):
            acc = S(zero)
            for a in range(N): acc = sadd(acc, sscal(v[a], ssub(qs[a], S(q0[a]))))
            if c == 1: acc = ssub(acc, [zero, one] + [zero] * (K - 2))   # e.(q-q0) = eps
            res.append(acc)
        return res
    # order 0: multiplier omega on J3 from vortex 0 (y0 = 0): dH/dx0 = omega dJ3/dx0
    Z0q = [S(x) for x in q0]
    nus0 = [S(zero) for _ in range(nc)]
    r0 = E(Z0q, nus0)
    omega = r0[0][0] / (2 * one / (1 + lam) ** 2)
    nus0[2] = S(omega)
    r0 = E(Z0q, nus0)
    assert all(is0(r[0]) for r in r0), 'ring is not a KKT point'
    # Jacobian M (order-1 coefficient of directional evaluation)
    NZ = N + nc
    M = [[zero] * NZ for _ in range(NZ)]
    for col in range(NZ):
        qs = [s[:] for s in Z0q]; ns = [s[:] for s in nus0]
        if col < N: qs[col][1] = qs[col][1] + one
        else: ns[col - N][1] = ns[col - N][1] + one
        r = E(qs, ns)
        for row in range(NZ): M[row][col] = r[row][1] - (r0[row][1])
    Zq = [s[:] for s in Z0q]; Zn = [s[:] for s in nus0]
    for k in range(1, K - 1 + 1 - 1):     # orders 1..3
        r = E(Zq, Zn)
        rhs = [-r[row][k] for row in range(NZ)]
        sol = solve(M, rhs)
        for a in range(N): Zq[a][k] = sol[a]
        for c in range(nc): Zn[c][k] = sol[N + c]
    r = E(Zq, Zn)
    assert all(is0(r[row][k]) for row in range(NZ) for k in range(4)), 'KKT not solved to order 3'
    nue = Zn[4]
    norm4 = (n * one / 2) ** 2 if odd else (n * one) ** 2
    return dict(omega=omega, nue1=nue[1], nue2=nue[2], c4=nue[3] / 4 * norm4, M=M, Zq=Zq)
