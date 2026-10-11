"""Reimplementation of Boatto (2008) Sec. 4-5 from her own equations, independent of ../blocks.py.

Canonical variables (her eq. 11/21), Gamma = 1:
    q_a = phi_a,  p_a = (cosh(sqrt|K| r_a) - 1)/K   (K<0),   (cos(sqrt K r_a) - 1)/K  (K>0)
Pair quantity (her eq. 22):
    d_ab = (1+K p_a)(1+K p_b) + K sqrt(p_a p_b (2+K p_a)(2+K p_b)) cos(q_a - q_b)
Hamiltonians (her eqs. 16-18, c = -1/(4 pi)):
    'boatto' : K<0: c sum ln((d-1)/(d+1))   K>0: c sum ln((1-d)/(1+d))   [= tanh / tan model]
    'sinh'   : K<0: c sum ln(d-1)           K>0: c sum ln(1-d)           [= sinh / whole-sphere S^2 model, her eq. 18]
Ring (her eq. 23 + Remark 2 of Sec 5.1): q_a = 2 pi (a-1)/N, p_a = p_o with geodesic radius r_o = 1.
Matrix monitored (her eq. 38, 44-45, Sec 5 text): Hessian S~ of H~ = H + nu M in (q,p) at the ring.
M_z = (1/sqrt|K|) sum (1 + K p_a) is LINEAR in p, so Hess(H~) = Hess(H) in (q,p).
We report the minimum eigenvalue over the 2N-1 eigenvalues other than the rotation zero (Sec 5 ii).

Second derivatives: pair function f(phi, pa, pb) differentiated symbolically (sympy) and evaluated with mpmath.
"""
import sys
import sympy as sp
import mpmath as mp

mp.mp.dps = 40
phi, pa, pb, Ks = sp.symbols('phi pa pb K', real=True)
c = -1 / (4 * sp.pi)
d = (1 + Ks * pa) * (1 + Ks * pb) + Ks * sp.sqrt(pa * pb * (2 + Ks * pa) * (2 + Ks * pb)) * sp.cos(phi)
FORMS = {
    ('boatto', -1): c * sp.log((d - 1) / (d + 1)),
    ('boatto', +1): c * sp.log((1 - d) / (1 + d)),
    ('sinh', -1): c * sp.log(d - 1),
    ('sinh', +1): c * sp.log(1 - d),
}
V = [phi, pa, pb]
_cache = {}


def pairfuns(model, sgn):
    key = (model, sgn)
    if key not in _cache:
        f = FORMS[key]
        H = [[sp.diff(f, u, v) for v in V] for u in V]
        G = [sp.diff(f, u) for u in V]
        _cache[key] = (sp.lambdify([phi, pa, pb, Ks], H, 'mpmath'),
                       sp.lambdify([phi, pa, pb, Ks], G, 'mpmath'))
    return _cache[key]


def p0(K, r0=1):
    K = mp.mpf(K)
    if K < 0:
        return (mp.cosh(mp.sqrt(-K) * r0) - 1) / K
    return (mp.cos(mp.sqrt(K) * r0) - 1) / K


def hessian(N, K, model='boatto', r0=1):
    """Full 2N x 2N Hessian of H in (q_1..q_N, p_1..p_N) at the ring; also returns gradient."""
    K = mp.mpf(K)
    sgn = -1 if K < 0 else 1
    Hf, Gf = pairfuns(model, sgn)
    po = p0(K, r0)
    S = mp.zeros(2 * N, 2 * N)
    grad = [mp.mpf(0)] * (2 * N)
    for a in range(N):
        for b in range(a + 1, N):
            ph = 2 * mp.pi * (a - b) / N
            h = Hf(ph, po, po, K)
            g = Gf(ph, po, po, K)
            # variables: phi = q_a - q_b, pa, pb.  Chain rule: dphi/dq_a = 1, dphi/dq_b = -1.
            J = {(a, 'q'): [1, 0, 0], (b, 'q'): [-1, 0, 0], (a, 'p'): [0, 1, 0], (b, 'p'): [0, 0, 1]}
            idx = {(a, 'q'): a, (b, 'q'): b, (a, 'p'): N + a, (b, 'p'): N + b}
            for k1, j1 in J.items():
                grad[idx[k1]] += sum(j1[u] * g[u] for u in range(3))
                for k2, j2 in J.items():
                    S[idx[k1], idx[k2]] += sum(j1[u] * h[u][v] * j2[v] for u in range(3) for v in range(3))
    return S, grad


def eigs(N, K, model='boatto', r0=1):
    S, grad = hessian(N, K, model, r0)
    ev = sorted(mp.eigsy(S, eigvals_only=True) if hasattr(mp, 'eigsy') else mp.eighe(S)[0])
    return ev, grad


def min_nontrivial(N, K, model='boatto', r0=1):
    S, grad = hessian(N, K, model, r0)
    E, Q = mp.eigsy(S)
    ev = [E[i] for i in range(2 * N)]
    # drop the eigenvalue closest to zero whose eigenvector is the rotation (1,..,1,0..0)/sqrt N
    rot = [Q[i, :] for i in range(2 * N)]
    best, bi = None, None
    for j in range(2 * N):
        ov = abs(sum(Q[i, j] for i in range(N))) / mp.sqrt(N)
        if best is None or ov > best:
            best, bi = ov, j
    rest = [ev[j] for j in range(2 * N) if j != bi]
    return min(rest), max(rest), ev[bi], best, grad


if __name__ == '__main__':
    model = sys.argv[1] if len(sys.argv) > 1 else 'boatto'
    for N in range(5, 10):
        S, grad = hessian(N, -0.25, model)
        gq = max(abs(grad[i]) for i in range(N))
        gp = max(abs(grad[N + i] - grad[N]) for i in range(N))
        print('N=%d K=-0.25 grad_q max %s, grad_p spread %s, dH/dp = %s' % (N, mp.nstr(gq, 3), mp.nstr(gp, 3), mp.nstr(grad[N], 8)))
