"""Quartic of F = H - w0*J3 on the LINEAR slice {dJ(x0)(q - x0) = 0} (no level-set curvature),
for comparison with the level-set c4 of my_c4.py. Usage: python slice_c4.py N"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from my_c4 import *

def solve_slice(n, eps, lam, z0, e, ep, t, q_init=None):
    N = 2*n; x0 = []
    for w in z0: x0 += [w.real, w.imag]
    gH0 = gradH(x0, n, lam); gJ0 = gradJ(x0, n, lam)
    w0 = sum(a*b for a, b in zip(gH0, gJ0[2]))/sum(b*b for b in gJ0[2])
    lin = [gJ0[0], gJ0[1], gJ0[2], e, t] + ([ep] if ep is not None else [])
    tg = [0, 0, 0, eps, 0] + ([0] if ep is not None else [])
    m = len(lin); q = list(q_init) if q_init else [x0[i]+eps*e[i] for i in range(N)]; nu = [mpf(0)]*m
    h = mpf(10)**-25
    gF = lambda qq: [a-w0*b for a, b in zip(gradH(qq, n, lam), gradJ(qq, n, lam)[2])]
    for it in range(60):
        g = gF(q)
        r1 = [g[a]-sum(nu[i]*lin[i][a] for i in range(m)) for a in range(N)]
        r2 = [sum(lin[i][a]*(q[a]-x0[a]) for a in range(N))-tg[i] for i in range(m)]
        if max(abs(v) for v in r1+r2) < mpf(10)**-45: break
        cols = fd_jac(gF, q, h); M = matrix(N+m, N+m)
        for a in range(N):
            for b in range(N): M[a, b] = cols[b][a]
            for i in range(m): M[a, N+i] = -lin[i][a]; M[N+i, a] = lin[i][a]
        d = lu_solve(M, matrix([-v for v in r1+r2]))
        q = [q[a]+d[a] for a in range(N)]; nu = [nu[i]+d[N+i] for i in range(m)]
    F = lambda qq: H(qq, n, lam)-w0*J(qq, n, lam)[2]
    return F(q)-F(x0), q

if __name__ == '__main__':
    n = int(sys.argv[1]); lam, z0 = setup(n); e, ep, t = kernel_vectors(n, z0)
    out = []; qi = None
    for s in ('1e-3', '2e-3', '4e-3'):
        eps = mpf(s); f, qi = solve_slice(n, eps, lam, z0, e, ep, t, qi); out.append((eps, f/eps**4))
    A = matrix([[1, a, a**2] for a, _ in out]); y = matrix([r for _, r in out]); sol = lu_solve(A, y)
    print(n, 'slice quartic (no curvature):', mp.nstr(sol[0], 12))
