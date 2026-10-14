"""Local equivalence lemma for flat realisations of the H^2 row/ring dynamics, checked by direct linearisation.

Setting. N like vortices, pair kernel W (the same in both systems), dynamics x' = m dH/dy, y' = -m dH/dx with a mobility
m that depends only on the transverse coordinate, and H = sum_{i<j} W_ij + sum_i U(transverse_i).
  * channel/strip (coordinates x along, y across): transverse = y, mu = m(y);
  * disk (polar r, theta):                         transverse = r, mu = m(r)/r   (theta' = mu dH/dr, r' = -mu dH/dtheta).
For a uniformly translating row (or rigidly rotating ring) at transverse position y0, the linearisation is
      J = mu(y0) * [[ -H_tx , -H_tt ], [ H_xx + S*I , H_xt ]],     S = U''(y0) + (mu'/mu)(y0) * (W_y + U'(y0)),
with t = along-row/angle, x = transverse, W_y the transverse interaction force on one vortex (the same for all).
So two systems with the same W have proportional linearisations (same stability, time rescaled by mu ratio) if they
have the same S. Everything else about m and U is irrelevant to linear stability.

Instances:
  H^2, strip model:  m = sin^2 y, U = 0                    -> S = 2 cot(a) W_y
  channel:           m = 1, U = (1/2) ln(2 sin y) + U_b(y) -> S = -csc^2(a)/2 + U_b''
  H^2, disk model:   m = ((1 - r^2)/2)^2, U = 0           -> S = -(4 r/(1 - r^2) + 1/r) W_r
  trap (Havelock):   m = 1, U = (1/2) ln(1 - r^2) + U_b(r) -> S = R'' - (R' + W_r)/r + r (U_b'/r)'
A background flow enters only through its SHEAR at the row: channel du_b/dy = U_b''; disk r dOmega_b/dr = r (U_b'/r)'.
A rigid rotation (uniform background vorticity in the disk) changes nothing; uniform vorticity in a channel is a Couette
shear and changes S by its shear rate.
Usage: python flat_lemma.py"""
import numpy as np

# ---------------- channel / strip: periodic supercell -----------------
def strip_W_grad(x1, y1, x2, y2):
    """W = (1/2) ln[(cosh dx - cos(y1+y2))/(cosh dx - cos(y1-y2))] and its gradient in (x1, y1)."""
    dx = np.clip(x1 - x2, -700, 700); ch, sh = np.cosh(dx), np.sinh(dx)
    P, M = ch - np.cos(y1 + y2), ch - np.cos(y1 - y2)
    Wx = 0.5*(sh/P - sh/M)
    Wy = 0.5*(np.sin(y1 + y2)/P - np.sin(y1 - y2)/M)
    return Wx, Wy

def strip_velocity(X, Y, period, mob, Up, J=40):
    """velocities of M vortices in a periodic supercell (period along x). mob(y) -> (m, m'), Up(y) -> U'(y)."""
    M = len(X); Hx = np.zeros(M); Hy = np.zeros(M)
    for i in range(M):
        for j in range(M):
            for k in range(-J, J + 1):
                if i == j and k == 0: continue
                wx, wy = strip_W_grad(X[i], Y[i], X[j] + k*period, Y[j])
                if i == j:   # self-images: only the y-dependence counts, and each image pair is counted once in H
                    wx = 0.0
                    # d/dy1 of (1/2) W(kP; y1, y1): differentiate both slots -> wy(first slot) + wy(second slot) = 2 wy
                    # times 1/2 -> wy. strip_W_grad gives the first-slot derivative only; by symmetry the second equals it.
                Hx[i] += wx; Hy[i] += wy
        Hy[i] += Up(Y[i])
    m = mob(Y)[0]
    return m*Hy, -m*Hx

def jacobian(vel, Z, h=1e-6):
    n = len(Z); Jm = np.zeros((n, n))
    for k in range(n):
        e = np.zeros(n); e[k] = h
        Jm[:, k] = (vel(Z + e) - vel(Z - e))/(2*h)
    return Jm

def strip_check(alpha, lp, M=6):
    P = M*lp; X0 = np.arange(M)*lp; Y0 = np.full(M, alpha)
    def make(mob, Up):
        def vel(Z):
            u, v = strip_velocity(Z[:M], Z[M:], P, mob, Up); return np.concatenate([u, v])
        return vel
    vH = make(lambda y: (np.sin(y)**2, np.sin(2*y)), lambda y: 0.0*y)
    # S for H^2 at the row: 2 cot(a) W_y, with W_y from the row itself
    _, Wy = strip_velocity(X0, Y0, P, lambda y: (np.ones_like(y), 0*y), lambda y: 0.0*y)
    Wy = Wy  # placeholder (velocity form); recompute force directly below
    Hy = 0.0
    for k in range(1, 400):
        Hy += 2*strip_W_grad(0.0, alpha, k*lp, alpha)[1]  # +-k images, first-slot derivative
    S_H = 2/np.tan(alpha)*Hy
    Sigma = S_H + 0.5/np.sin(alpha)**2          # Couette shear that matches
    vC = make(lambda y: (np.ones_like(y), 0*y), lambda y: 0.5/np.tan(y) + Sigma*(y - alpha))
    Z0 = np.concatenate([X0, Y0])
    eH = np.linalg.eigvals(jacobian(vH, Z0)); eC = np.linalg.eigvals(jacobian(vC, Z0))
    v0 = make(lambda y: (np.ones_like(y), 0*y), lambda y: 0.5/np.tan(y))      # control: channel without shear
    e0 = np.linalg.eigvals(jacobian(v0, Z0))
    return match(eH, eC*np.sin(alpha)**2), np.abs(eH.real).max(), Sigma, np.abs(e0.real).max()*np.sin(alpha)**2

# ---------------- disk: ring of N -----------------
def disk_velocity(Z, mob, Up):
    """Z = [x..., y...]; W = ln|1 - z conj(w)| - ln|z - w|, self term U(r) (includes Robin when wanted)."""
    N = len(Z)//2; z = Z[:N] + 1j*Z[N:]
    g = np.zeros(N, complex)   # g = dH/dx + i dH/dy
    for i in range(N):
        for j in range(N):
            if i == j: continue
            # d/dz_i-bar of W: grad = 2 dW/dzbar ; W = Re ln(1 - z wb) - Re ln(z - w)
            g[i] += -z[j]/(1 - np.conj(z[i])*z[j]) - 1/np.conj(z[i] - z[j])   # grad W = 2 dW/dzbar
    r = np.abs(z); g += Up(r)*z/r
    m = mob(r)[0]
    return np.concatenate([m*g.imag, -m*g.real])

def _grad_check():
    rng = np.random.default_rng(1); Z = rng.uniform(-0.4, 0.4, 6)
    def H(Z):
        z = Z[:3] + 1j*Z[3:]; s = 0
        for i in range(3):
            for j in range(i+1, 3): s += np.log(abs(1 - z[i]*np.conj(z[j]))) - np.log(abs(z[i] - z[j]))
        return s
    num = np.array([(H(Z + e) - H(Z - e))/2e-7 for e in 1e-7*np.eye(6)])
    v = disk_velocity(Z, lambda r: (np.ones_like(r), 0*r), lambda r: 0*r)
    return max(abs(v[:3] - num[3:]).max(), abs(v[3:] + num[:3]).max())

def ring_eigs(N, r0, mob, Up):
    th = 2*np.pi*np.arange(N)/N; Z0 = np.concatenate([r0*np.cos(th), r0*np.sin(th)])
    v0 = disk_velocity(Z0, mob, Up); Om = v0[N]/r0     # vortex 0 at (r0, 0) moves in +y
    def vrot(Z):
        v = disk_velocity(Z, mob, Up); return v - Om*np.concatenate([-Z[N:], Z[:N]])
    assert abs(vrot(Z0)).max() < 1e-9*max(1, abs(Om))
    return np.linalg.eigvals(jacobian(vrot, Z0, 1e-6)), Om

def ring_force(N, r0):
    """W_r: radial derivative of sum_j W(z_0, z_j) at vortex 0."""
    th = 2*np.pi*np.arange(1, N)/N; w = r0*np.exp(1j*th); h = 1e-6
    f = lambda r: np.sum(np.log(abs(1 - r*np.conj(w))) - np.log(abs(r - w)))
    return (f(r0 + h) - f(r0 - h))/(2*h)

def trap_background_shear(N, r0):
    """the background shear r dOmega_b/dr (m = 1 units, vortex circulation 2 pi) that makes the trap ring equivalent."""
    Wr = ring_force(N, r0); r = r0
    S_H = -(4*r/(1 - r*r) + 1/r)*Wr
    Rp, Rpp = -r/(1 - r*r), -(1 + r*r)/(1 - r*r)**2
    return S_H - (Rpp - (Rp + Wr)/r)

def disk_check(N, r0):
    sh = trap_background_shear(N, r0)
    mobH = lambda r: (((1 - r*r)/2)**2, -r*(1 - r*r))
    eH, _ = ring_eigs(N, r0, mobH, lambda r: 0*r)
    # background with U_b'/r = (sh/(2 r0^2)) r^2  ->  r (U_b'/r)' = sh at r0 (quartic stream function)
    UpT = lambda r: -r/(1 - r*r) + sh/(2*r0**2)*r**3
    eT, _ = ring_eigs(N, r0, lambda r: (np.ones_like(r), 0*r), UpT)
    scale = ((1 - r0*r0)/2)**2
    return match(eH, eT*scale), np.abs(eH.real).max(), sh

def match(a, b):
    """max distance from each eigenvalue of a to the nearest of b (and back): robust to sort order."""
    D = np.abs(a[:, None] - b[None, :]); return max(D.min(axis=1).max(), D.min(axis=0).max())

if __name__ == '__main__':
    print(f'disk gradient check: {_grad_check():.1e}')
    print('strip: H^2 row vs channel row + Couette shear Sigma* (eigenvalues, channel rescaled by sin^2 a):')
    for al, lp in [(0.4, 0.8), (0.6, 1.5), (0.9, 1.0), (1.2, 2.0)]:
        d, mx, Sg, m0 = strip_check(al, lp)
        print(f'  alpha {al}, l\' {lp}: max |eig diff| = {d:.1e}  (max |Re eig| H^2 = {mx:.3f}, channel without shear = {m0:.3f}),  Sigma* = {Sg:.4f}')
    print('disk: H^2 ring vs trap ring + background shear (eigenvalues, trap rescaled by ((1-r^2)/2)^2):')
    for N, r0 in [(3, 0.3), (5, 0.5), (8, 0.6), (8, 0.85), (10, 0.9)]:
        d, mx, sh = disk_check(N, r0)
        print(f'  N {N}, r0 {r0}: max |eig diff| = {d:.1e}  (max |Re eig| = {mx:.3f}),  r dOmega_b/dr = {sh:.4f}')
