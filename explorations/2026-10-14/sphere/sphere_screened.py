"""Sphere rings with the SCREENED kernel W_mu = Green's function of (-Delta + mu) on the unit sphere,
W_mu(c) propto sum_l (2l+1) P_l(c)/(l(l+1) + mu) = -pi P_nu(-c)/sin(pi nu),  nu(nu+1) = -mu, c = cos d.
mu -> 0 recovers -ln sin(d/2) (plus a constant). Question: does the N = 4 threshold stay at kappa = cot(theta) = 1/sqrt 2?
Usage: python sphere_screened.py"""
import mpmath as mp
mp.mp.dps = 30

def make_wp(mu):
    if mu == 0:
        return lambda c: 1/(2*(1 - c))                 # d/dc of -(1/2) ln(1 - c)
    nu = -mp.mpf(1)/2 + mp.sqrt(mp.mpf(1)/4 - mu + 0j)
    W = lambda c: mp.re(-mp.pi*mp.legenp(nu, 0, -c)/mp.sin(mp.pi*nu))
    return lambda c: mp.diff(W, c)
def field(z, p, wp):
    N = len(z); S = [mp.sqrt(1 - x**2) for x in z]; dz = []; dp = []
    for i in range(N):
        hz = 0; hp = 0
        for j in range(N):
            if j == i: continue
            cs, sn = mp.cos(p[i]-p[j]), mp.sin(p[i]-p[j])
            c = z[i]*z[j] + S[i]*S[j]*cs; w1 = wp(c)
            hz += w1*(z[j] - z[i]*S[j]/S[i]*cs); hp += w1*(-S[i]*S[j]*sn)
        dz.append(-hp); dp.append(hz)
    return dz, dp
def growth(N, z0, wp):
    z = [mp.mpf(z0)]*N; p = [2*mp.pi*k/N for k in range(N)]
    dz0, dp0 = field(z, p, wp); om = dp0[0]
    F = lambda v: (lambda a: a[0] + [x - om for x in a[1]])(field(v[:N], v[N:], wp))
    h = mp.mpf('1e-12'); x0 = z + p; J = mp.matrix(2*N, 2*N)
    for a in range(2*N):
        xp = list(x0); xm = list(x0); xp[a] += h; xm[a] -= h; fp, fm = F(xp), F(xm)
        for b in range(2*N): J[b, a] = (fp[b] - fm[b])/(2*h)
    ev = mp.eig(J, left=False, right=False)
    return max(mp.re(e) for e in ev)/max(abs(e) for e in ev)
def threshold(N, wp, lo=0.05, hi=0.99, tol=1e-9):
    """smallest z0 above which the ring is stable (bisection on a sign change found by scanning)."""
    zs = [lo + (hi - lo)*i/60 for i in range(61)]; st = [growth(N, x, wp) < 1e-9 for x in zs]
    i = next(k for k in range(60) if not st[k] and st[k+1])
    a, b = mp.mpf(zs[i]), mp.mpf(zs[i+1])
    while b - a > tol:
        m = (a + b)/2
        if growth(N, m, wp) < 1e-9: b = m
        else: a = m
    return (a + b)/2
if __name__ == '__main__':
    for mu in [-0.5, 0, 0.5, 2, 6]:
        wp = make_wp(mp.mpf(mu))
        z = threshold(4, wp); k = z/mp.sqrt(1 - z**2)
        print(f'mu = {mu:5}: N = 4 threshold cos^2(theta) = {mp.nstr(z**2, 10)},  kappa = cot(theta) = {mp.nstr(k, 10)}  (1/sqrt2 = 0.7071067812)', flush=True)
