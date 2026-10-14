"""Rings of N equal point vortices on a latitude circle of the unit sphere (the K = +1 member of the family).
Kernel: Green's function of the sphere, W = -ln sin(d/2) = -(1/2) ln(1 - cos d) + const (uniform compensating background).
Canonical coordinates (z = cos theta, phi): dz_i/dt = -dH/dphi_i, dphi_i/dt = dH/dz_i (unit circulations; the overall sign
does not affect stability). Linear stability from the Jacobian at the uniformly rotating ring (finite differences, mpmath).
Usage: python sphere_rings.py"""
import mpmath as mp
mp.mp.dps = 30

def H(z, p):
    N = len(z); h = 0
    for i in range(N):
        for j in range(i+1, N):
            c = z[i]*z[j] + mp.sqrt((1-z[i]**2)*(1-z[j]**2))*mp.cos(p[i]-p[j])
            h += -mp.log(1 - c)/2
    return h
def field(z, p):
    """analytic: dz_i/dt = -dH/dphi_i, dphi_i/dt = dH/dz_i, H = -(1/2) sum ln(1 - c_ij)."""
    N = len(z); dz = []; dp = []
    S = [mp.sqrt(1 - x**2) for x in z]
    for i in range(N):
        hz = 0; hp = 0
        for j in range(N):
            if j == i: continue
            cs, sn = mp.cos(p[i]-p[j]), mp.sin(p[i]-p[j])
            c = z[i]*z[j] + S[i]*S[j]*cs
            hz += (z[j] - z[i]*S[j]/S[i]*cs)/(2*(1 - c))
            hp += (-S[i]*S[j]*sn)/(2*(1 - c))
        dz.append(-hp); dp.append(hz)
    return dz, dp
def growth(N, z0):
    z = [mp.mpf(z0)]*N; p = [2*mp.pi*k/N for k in range(N)]
    dz0, dp0 = field(z, p); om = dp0[0]
    def F(v):
        dz, dp = field(v[:N], v[N:]); return dz + [x - om for x in dp]
    h = mp.mpf('1e-14'); x0 = z + p; J = mp.matrix(2*N, 2*N)
    for a in range(2*N):
        xp = list(x0); xm = list(x0); xp[a] += h; xm[a] -= h
        fp, fm = F(xp), F(xm)
        for b in range(2*N): J[b, a] = (fp[b] - fm[b])/(2*h)
    ev = mp.eig(J, left=False, right=False)
    return max(mp.re(e) for e in ev)/max(abs(e) for e in ev)
if __name__ == '__main__':
    for N in range(2, 9):
        row = []
        for z0 in [0.05, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]:
            g = growth(N, z0); row.append('S' if g < 1e-10 else 'u')
        print(N, ' '.join(row))
