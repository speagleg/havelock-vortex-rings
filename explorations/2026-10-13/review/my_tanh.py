"""From-scratch check of tanh-model ring thresholds (Poincare disk, Cartesian, no block formulas).
Usage: python my_tanh.py N xi1 xi2 ...   (prints criticality residual, max |Re eigenvalue|,
and the smallest Hessian eigenvalue on the non-rotational modes)"""
import sys
from mpmath import mp, mpf, mpc, sqrt, cos, sin, pi, matrix, eig, eigsy
mp.dps = 50

def ring(N, xi):
    r = sqrt(xi)
    return [mpc(r*cos(2*pi*k/N), r*sin(2*pi*k/N)) for k in range(N)]

def gradH(z):
    N = len(z); g = []
    for k in range(N):
        c = mpc(0)
        for j in range(N):
            if j == k: continue
            d = z[k]-z[j]
            c += -d/abs(d)**2
            c += -z[j]/(1 - z[k].conjugate()*z[j])
        g.append(c)
    return g

def gradJ(z):
    return [4*w/(1-abs(w)**2)**2 for w in z]

def to_real(cl):
    out = []
    for c in cl: out += [c.real, c.imag]
    return out

def from_real(q):
    return [mpc(q[2*k], q[2*k+1]) for k in range(len(q)//2)]

def omega(z):
    gH = to_real(gradH(z)); gJ = to_real(gradJ(z))
    return sum(a*b for a, b in zip(gH, gJ))/sum(b*b for b in gJ)

def hessF(z, Om):
    q = to_real(z); n = len(q); h = mpf(10)**-20
    def gF(qq):
        zz = from_real(qq)
        return [a - Om*b for a, b in zip(to_real(gradH(zz)), to_real(gradJ(zz)))]
    M = matrix(n, n)
    for i in range(n):
        qp = list(q); qm = list(q); qp[i] += h; qm[i] -= h
        fp = gF(qp); fm = gF(qm)
        for r in range(n):
            M[r, i] = (fp[r]-fm[r])/(2*h)
    for i in range(n):
        for j in range(i+1, n):
            v = (M[i, j]+M[j, i])/2; M[i, j] = v; M[j, i] = v
    return M

def analyse(N, xi):
    z = ring(N, xi); Om = omega(z)
    res = max(abs(a - Om*b) for a, b in zip(to_real(gradH(z)), to_real(gradJ(z))))
    Hs = hessF(z, Om); n = 2*N
    lam = 4/(1-xi)**2
    L = matrix(n, n)
    for k in range(N):
        for c in range(n):
            L[2*k, c] = Hs[2*k+1, c]/lam
            L[2*k+1, c] = -Hs[2*k, c]/lam
    ev = eig(L, left=False, right=False)
    maxre = max(abs(e.real) for e in ev)
    rad = []; tan = []
    for w in z:
        u = w/abs(w); rad += [u.real, u.imag]; tan += [-u.imag, u.real]
    vecs = [rad, tan]
    vecs = [[x/sqrt(sum(y*y for y in v)) for x in v] for v in vecs]
    basis = []
    for i in range(n):
        v = [mpf(1) if j == i else mpf(0) for j in range(n)]
        for b in vecs+basis:
            d = sum(x*y for x, y in zip(v, b)); v = [x-d*y for x, y in zip(v, b)]
        s = sqrt(sum(x*x for x in v))
        if s > mpf('1e-10'):
            basis.append([x/s for x in v])
        if len(basis) == n-2: break
    P = matrix(basis)
    evs = eigsy(P*Hs*P.T, eigvals_only=True)
    return res, maxre, min(evs)

if __name__ == '__main__':
    N = int(sys.argv[1])
    for x in [mpf(x) for x in sys.argv[2:]]:
        res, mr, mn = analyse(N, x)
        print(N, mp.nstr(x, 20), 'crit_res', mp.nstr(res, 3), 'max|Re ev|', mp.nstr(mr, 5), 'min Hess eig (modes!=0)', mp.nstr(mn, 8))
