"""Is the nearest-neighbour vortex row integrable (Toda-like)? Largest Lyapunov exponent (Benettin) of periodic chains:
 (a) nearest-neighbour row H = sum_j 2 s_j s_{j+1} e^{-(u_{j+1}-u_j)}, s = (1+Y^2)^{-1/2}, u_{N+1} = u_1 + L;
 (b) the full tanh-kernel row (all pairs, periodic images), the same N vortices;
 (c) the periodic Toda lattice (integrable control) and a Toda lattice with an added cubic on-site term (chaotic control).
Canonical pairs (u, Y) and (q, p). Usage: python toda_test.py -> toda_test.out"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
from zigzag_dynamics import rhs_factory           # full tanh kernel, periodic images, (u, phi) coordinates
def nn_rhs(L):
    def f(t, z):
        N = len(z)//2; u, Y = z[:N], z[N:]; s = 1/np.sqrt(1 + Y**2); dlns = -Y/(1 + Y**2)
        D = np.roll(u, -1) - u; D[-1] += L
        w = 2*s*np.roll(s, -1)*np.exp(-D)
        dHdu = w - np.roll(w, 1)                    # d/du_j of w_j (+) and of w_{j-1} (-)
        dHdY = (w + np.roll(w, 1))*dlns
        return np.concatenate([dHdY, -dHdu])
    return f
def toda_rhs(L, eps=0.0):
    def f(t, z):
        N = len(z)//2; q, p = z[:N], z[N:]; D = np.roll(q, -1) - q; D[-1] += L
        e = np.exp(-D); dV = e - np.roll(e, 1)    # V = sum e^{-(q_{j+1}-q_j)}; dV/dq_j = e_j - e_{j-1}
        dV += eps*3*(q - q.mean())**2*0 + eps*np.sin(q)  # chaotic control: on-site sin potential
        return np.concatenate([p, -dV])
    return f
def full_rhs(L, N):
    g = rhs_factory(N, L)
    def f(t, z):
        u, Y = z[:N], z[N:]; ph = np.pi/2 + np.arctan(Y)
        r = g(t, np.concatenate([u, ph])); du, dph = r[:N], r[N:]
        return np.concatenate([du, dph/np.sin(ph)**2])   # dY/dt = dphi/dt / sin^2 phi
    return f
def lyap(rhs, z0, T=600.0, dt=2.0, d0=1e-8, seed=0):
    rng = np.random.default_rng(seed); v = rng.standard_normal(len(z0)); v /= np.linalg.norm(v)
    z, zp = z0.copy(), z0 + d0*v; acc = 0.0; t = 0.0
    while t < T:
        z = solve_ivp(rhs, (0, dt), z, rtol=1e-11, atol=1e-13, method='DOP853').y[:, -1]
        zp = solve_ivp(rhs, (0, dt), zp, rtol=1e-11, atol=1e-13, method='DOP853').y[:, -1]
        d = np.linalg.norm(zp - z); acc += np.log(d/d0); zp = z + (zp - z)*d0/d; t += dt
    return acc/t
if __name__ == "__main__":
    out = []; P = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
    N = 4; L = N*1.0; rng = np.random.default_rng(7)
    for trial in range(2):
        u0 = np.arange(N)*L/N + 0.15*rng.standard_normal(N); Y0 = -1.2 + 0.6*rng.standard_normal(N)
        z0 = np.concatenate([u0, Y0])
        lnn = lyap(nn_rhs(L), z0); lfull = lyap(full_rhs(L, N), z0)
        q0 = np.arange(N)*1.0 + 0.4*rng.standard_normal(N); p0 = 0.8*rng.standard_normal(N)
        lt = lyap(toda_rhs(N*1.0), np.concatenate([q0, p0])); ltc = lyap(toda_rhs(N*1.0, eps=1.0), np.concatenate([q0, p0]))
        P(f'trial {trial}: max Lyapunov (T = 600): NN vortex row {lnn:.3e} | full tanh row {lfull:.3e} | periodic Toda {lt:.3e} | Toda + sin {ltc:.3e}')
    open('toda_test.out', 'w').write('\n'.join(out) + '\n')
