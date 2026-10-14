"""Infinite vortex rows on curves of constant geodesic curvature kappa in H^2 (tanh model).
kappa = 0 geodesic, 0 < kappa < 1 hypercycle, kappa = 1 horocycle (kappa > 1: circles = the finite rings).
Kernel W(d) = -ln tanh(d/2) (the tanh model's Green's function).

Hypercycles: upper half-plane z = exp(u + i phi); the hypercycle phi = alpha has kappa = cos(alpha) and the
hyperbolic translation along it is u -> u + const. Metric |dw|^2/sin^2(phi), area du dphi / sin^2(phi),
cosh d_ij = 1 + [cosh(u_i - u_j) - cos(phi_i - phi_j)]/(sin phi_i sin phi_j),
du_i/dt = sin^2(phi_i) dH/dphi_i,  dphi_i/dt = -sin^2(phi_i) dH/du_i.
Row: u_k = k l', phi_k = alpha; arc spacing along the curve l = l'/sin(alpha).
Horocycles: z = x + i y, row z_k = k l + i, dz_i/dt = -i y_i^2 sum_j [1/(conj z_i - z_j) - 1/(conj z_i - conj z_j)]
(power-law interactions along a horocycle: many images needed).
A supercell of M vortices with periodic images gives Bloch modes theta = j/M; stability from the eigenvalues
of the Jacobian of the co-moving velocity (relative growth < 1e-4 counts as stable; the zero-mode Jordan block
leaves ~3e-5).
Validated: horocycle threshold falls in (0.840, 0.845) vs the note's l* = 0.8439085; the log-coordinate velocity
agrees with the upper-half-plane formula to 2e-12.
Usage: python rows.py"""
import numpy as np

# ---------- hypercycles, log coordinates ----------
def vel_hyper(u, phi, L, nimg):
    du = u[:, None, None] - (u[None, :, None] + L*np.arange(-nimg, nimg+1)[None, None, :])
    dp = phi[:, None, None] - phi[None, :, None]
    si = np.sin(phi)[:, None, None]; sj = np.sin(phi)[None, :, None]
    ok = np.abs(du) < 300
    dus = np.where(ok, du, 0.0)
    C = 1 + (np.cosh(dus) - np.cos(dp))/(si*sj)
    M = len(u); n0 = nimg
    ok[np.arange(M), np.arange(M), n0] = False
    Cs = np.where(ok, C, 2.0)
    Wp = np.where(ok, -1/(Cs**2 - 1), 0.0)                 # dW/dC
    dCdu = np.sinh(dus)/(si*sj)
    dCdp = np.sin(dp)/(si*sj) - (Cs-1)*np.cos(phi)[:, None, None]/si
    Hu = (Wp*dCdu).sum(axis=(1, 2)); Hp = (Wp*dCdp).sum(axis=(1, 2))
    s2 = np.sin(phi)**2
    return s2*Hp, -s2*Hu

def stab_hyper(M, l, kappa, h=1e-7, nimg=None):
    alpha = np.arccos(kappa); lp = l*np.sin(alpha); L = M*lp
    nimg = nimg or max(1, int(60/L)+1)
    u = np.arange(M)*lp; phi = np.full(M, alpha)
    ud, pd = vel_hyper(u, phi, L, nimg)
    s = ud[0]; res = max(np.max(np.abs(ud-s)), np.max(np.abs(pd)))/max(abs(s), 1e-300)
    def F(q):
        a, b = vel_hyper(q[:M], q[M:], L, nimg)
        return np.concatenate([a - s, b])
    q0 = np.concatenate([u, phi]); J = np.zeros((2*M, 2*M))
    for a in range(2*M):
        qp = q0.copy(); qm = q0.copy(); qp[a] += h; qm[a] -= h
        J[:, a] = (F(qp)-F(qm))/(2*h)
    ev = np.linalg.eigvals(J)
    return np.max(ev.real)/np.max(np.abs(ev)), res, s

# ---------- horocycle ----------
def vel_horo(z, L, nimg):
    W = z[None, :] + (L*np.arange(-nimg, nimg+1))[:, None]
    zb = np.conj(z)[:, None, None]
    with np.errstate(divide='ignore', invalid='ignore'):
        T = 1/(zb - W[None]) - 1/(zb - np.conj(W)[None])
    M = len(z); T[np.arange(M), nimg, np.arange(M)] = 0
    return -1j*np.imag(z)**2*T.sum(axis=(1, 2))

def stab_horo(M, l, h=1e-7, nimg=400):
    z = np.arange(M)*l + 1j; L = M*l
    v = vel_horo(z, L, nimg); s = v[0].real
    q0 = np.concatenate([z.real, z.imag]); J = np.zeros((2*M, 2*M))
    F = lambda q: vel_horo(q[:M]+1j*q[M:], L, nimg) - s
    for a in range(2*M):
        qp = q0.copy(); qm = q0.copy(); qp[a] += h; qm[a] -= h
        d = (F(qp)-F(qm))/(2*h); J[:M, a] = d.real; J[M:, a] = d.imag
    ev = np.linalg.eigvals(J)
    return np.max(ev.real)/np.max(np.abs(ev)), None, s

def critical(stab, lo, hi, tol=1e-4):
    st = lambda l: stab(l)[0] < tol
    if st(lo): return 0.0
    if not st(hi): return np.inf
    for _ in range(40):
        m = (lo+hi)/2
        if st(m): hi = m
        else: lo = m
    return hi

if __name__ == '__main__':
    print('horocycle critical spacing (note: l* = 0.8439085):', critical(lambda l: stab_horo(20, l), 0.5, 1.5))
    for kap in [0.0, 0.5, 0.9, 0.99, 0.999]:
        print(f'kappa={kap}:', ' '.join(f'l={l}: {stab_hyper(40, l, kap)[0]:.1e}' for l in [0.5, 1, 2, 4, 8, 16]), flush=True)
    lo, hi = 0.70, 0.75
    for _ in range(30):
        m = (lo+hi)/2
        if stab_hyper(40, 80, m)[0] < 1e-3: hi = m
        else: lo = m
    print('critical curvature at l = 80:', hi, ' (1/sqrt 2 =', 1/np.sqrt(2), ')')
