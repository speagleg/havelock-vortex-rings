"""Which nearest-neighbour vortex chains on H^2 look integrable? Periodic chain of N vortices, canonical (u, Y = -cot phi),
H = sum_j f(C_j), C_j = cosh d(z_j, z_{j+1}) = 1 + (cosh(du) - cos(dphi))/(sin phi_j sin phi_{j+1}), u_{N+1} = u_1 + L.
Kernels: sinh model f = -ln sinh(d/2); tanh model f = -ln tanh(d/2); its large-d asymptote f = 2 e^{-d};
control f = -ln cosh(d/2). Gradients by complex-step differentiation (exact to rounding).
Finite-time largest Lyapunov exponent at T and 4T: integrable chains decay ~ 1/T, chaotic ones plateau.
Hyperboloid picture: a vortex is a classical SU(1,1) spin S with S.S' = -cosh d; the integrable lattice Landau-Lifshitz
chain (Ishimori 1982; Haldane 1982) has H = sum ln(1 + S_j.S_{j+1}), whose hyperbolic analogue is ln(cosh d - 1)
= 2 ln sinh(d/2) + const, i.e. the sinh kernel.  Usage: python kernel_integrability.py -> kernel_integrability.out"""
import numpy as np
from scipy.integrate import solve_ivp
from toda_test import lyap
KERNELS = {
    'sinh  -ln sinh(d/2)': lambda C: -0.5*np.log((C - 1)/2),
    'tanh  -ln tanh(d/2)': lambda C: 0.5*np.log((C + 1)/(C - 1)),
    'asym  2 e^{-d}':      lambda C: 2/(C + np.sqrt(C*C - 1)),
    'cosh  -ln cosh(d/2)': lambda C: -0.5*np.log((C + 1)/2),
}
def H_of(f, L):
    def H(z):
        N = len(z)//2; u, Y = z[:N], z[N:]; ph = np.pi/2 + np.arctan(Y)
        du = np.roll(u, -1) - u; du[-1] = du[-1] + L
        C = 1 + (np.cosh(du) - np.cos(np.roll(ph, -1) - ph))/(np.sin(ph)*np.sin(np.roll(ph, -1)))
        return np.sum(f(C))
    return H
def rhs_of(f, L):
    H = H_of(f, L)
    def rhs(t, z):
        n = len(z); g = np.empty(n); h = 1e-30
        for i in range(n):
            zc = z.astype(complex); zc[i] += 1j*h; g[i] = H(zc).imag/h
        N = n//2; return np.concatenate([g[N:], -g[:N]])
    return rhs
if __name__ == '__main__':
    N = 4; L = N*1.0; rng = np.random.default_rng(7)
    u0 = np.arange(N)*L/N + 0.15*rng.standard_normal(N); Y0 = -1.2 + 0.6*rng.standard_normal(N); z0 = np.concatenate([u0, Y0])
    out = open('kernel_integrability.out', 'w')
    for name, f in KERNELS.items():
        r = rhs_of(f, L); H = H_of(f, L)
        l1 = lyap(r, z0, T=300.0); l2 = lyap(r, z0, T=1200.0)
        line = f'{name:22s}: lambda(T=300) = {l1:.3e}, lambda(T=1200) = {l2:.3e}, ratio {l1/l2:5.2f}  (integrable ~4, chaotic ~1)'
        print(line, flush=True); out.write(line + '\n'); out.flush()
