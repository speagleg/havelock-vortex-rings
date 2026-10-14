"""Independent from-scratch check of the level-set quartic c4 (MT family).
Original chart, Cartesian coords, KKT Newton, mpmath. Shares nothing with the repo.
Usage: python my_c4.py N   (prints c4 extrapolated from eps = 1e-3, 2e-3, 4e-3)"""
import sys
from mpmath import mp, mpf, mpc, sqrt, cos, sin, pi, matrix, lu_solve, log
mp.dps = 60

def bn(n):
    # root in (-1,1] of (1+b^2)/(1+b)^2 = floor(n^2/4)/(2(n-1))
    R = mpf(n*n//4)/(2*(n-1))
    a, bb, c = 1-R, -2*R, 1-R
    if abs(a) < mpf(10)**-50:
        return -c/bb
    d = sqrt(bb*bb-4*a*c)
    roots = [r for r in [(-bb+d)/(2*a), (-bb-d)/(2*a)] if -1 < r <= 1]
    assert len(roots) == 1, roots
    return roots[0]

def setup(n):
    return bn(n), [mpc(cos(2*pi*k/n), sin(2*pi*k/n)) for k in range(n)]

def unpack(q, n):
    return [mpc(q[2*k], q[2*k+1]) for k in range(n)]

def H(q, n, lam):
    z = unpack(q, n); s = mpf(0)
    for i in range(n):
        for j in range(i+1, n):
            s -= log(abs(z[i]-z[j])**2)
    for k in range(n):
        s += (n-1)*log(1+lam*abs(z[k])**2)
    return s

def gradH(q, n, lam):
    z = unpack(q, n); g = []
    for k in range(n):
        c = mpc(0)
        for j in range(n):
            if j != k:
                d = z[k]-z[j]; c -= 2*d/abs(d)**2
        c += (n-1)*2*lam*z[k]/(1+lam*abs(z[k])**2)
        g += [c.real, c.imag]
    return g

def J(q, n, lam):
    z = unpack(q, n)
    return [sum(w.real/(1+lam*abs(w)**2) for w in z), sum(w.imag/(1+lam*abs(w)**2) for w in z),
            sum(abs(w)**2/(1+lam*abs(w)**2) for w in z)]

def gradJ(q, n, lam):
    z = unpack(q, n); G = [[], [], []]
    for w in z:
        x, y = w.real, w.imag; s = 1+lam*(x*x+y*y)
        G[0] += [1/s - x*2*lam*x/s**2, -x*2*lam*y/s**2]
        G[1] += [-y*2*lam*x/s**2, 1/s - y*2*lam*y/s**2]
        G[2] += [2*x/s**2, 2*y/s**2]
    return G

def fd_jac(f, q, h):
    cols = []
    for i in range(len(q)):
        qp = list(q); qm = list(q); qp[i] += h; qm[i] -= h
        fp = f(qp); fm = f(qm)
        cols.append([(a-b)/(2*h) for a, b in zip(fp, fm)])
    return cols

def solve(n, eps, lam, z0, e, ep, t, q_init=None, nu_init=None):
    N = 2*n; x0 = []
    for w in z0: x0 += [w.real, w.imag]
    mu = J(x0, n, lam)
    lin = [e, t] + ([ep] if ep is not None else [])
    targets_lin = [eps, mpf(0)] + ([mpf(0)] if ep is not None else [])
    m = 3 + len(lin)
    q = list(q_init) if q_init else [x0[i] + eps*e[i] for i in range(N)]
    if nu_init: nu = list(nu_init)
    else:
        gH0 = gradH(x0, n, lam); gJ0 = gradJ(x0, n, lam)
        w = sum(a*b for a, b in zip(gH0, gJ0[2]))/sum(b*b for b in gJ0[2])
        nu = [mpf(0), mpf(0), w] + [mpf(0)]*(m-3)
    h = mpf(10)**-25
    for it in range(60):
        gH = gradH(q, n, lam); gJ = gradJ(q, n, lam); cg = gJ + lin
        r1 = [gH[a] - sum(nu[i]*cg[i][a] for i in range(m)) for a in range(N)]
        Jv = J(q, n, lam)
        r2 = [Jv[i]-mu[i] for i in range(3)] + [sum(lin[i][a]*(q[a]-x0[a]) for a in range(N)) - targets_lin[i] for i in range(len(lin))]
        res = max(abs(v) for v in r1+r2)
        if res < mpf(10)**-45: break
        def gL(qq):
            gh = gradH(qq, n, lam); gj = gradJ(qq, n, lam)
            return [gh[a] - sum(nu[i]*gj[i][a] for i in range(3)) for a in range(N)]
        cols = fd_jac(gL, q, h)
        M = matrix(N+m, N+m)
        for a in range(N):
            for b in range(N): M[a, b] = cols[b][a]
            for i in range(m): M[a, N+i] = -cg[i][a]; M[N+i, a] = cg[i][a]
        d = lu_solve(M, matrix([-v for v in r1+r2]))
        q = [q[a]+d[a] for a in range(N)]; nu = [nu[i]+d[N+i] for i in range(m)]
    else:
        raise RuntimeError('no convergence %s' % res)
    return H(q, n, lam) - H(x0, n, lam), q, nu, it

def kernel_vectors(n, z0):
    if n % 2 == 0:
        e = []
        for k in range(n):
            u = z0[k]; e += [((-1)**k)*u.real/sqrt(n), ((-1)**k)*u.imag/sqrt(n)]
        ep = None
    else:
        l = (n-1)//2; e = []; ep = []
        for k in range(n):
            u = z0[k]; ph = 2*pi*k/n
            e += [sqrt(mpf(2)/n)*cos(l*ph)*u.real, sqrt(mpf(2)/n)*cos(l*ph)*u.imag]
            ep += [sqrt(mpf(2)/n)*sin(l*ph)*u.real, sqrt(mpf(2)/n)*sin(l*ph)*u.imag]
    t = []
    for u in z0: t += [-u.imag, u.real]
    return e, ep, t

def c4(n, epss):
    lam, z0 = setup(n); e, ep, t = kernel_vectors(n, z0)
    out = []; qi = nui = None
    for eps in epss:
        f, qi, nui, it = solve(n, eps, lam, z0, e, ep, t, qi, nui)
        out.append((eps, f, f/eps**4))
    return lam, out

if __name__ == '__main__':
    n = int(sys.argv[1])
    lam, out = c4(n, [mpf(x) for x in ('1e-3', '2e-3', '4e-3')])
    A = matrix([[1, e, e**2] for e, _, _ in out]); y = matrix([r for _, _, r in out])
    sol = lu_solve(A, y)
    print(n, mp.nstr(lam, 12), [mp.nstr(r, 12) for _, _, r in out], 'c4_extrap=', mp.nstr(sol[0], 15))
