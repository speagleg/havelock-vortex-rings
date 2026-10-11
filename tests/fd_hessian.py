"""Brute-force finite-difference Hessian of the full 2N-dimensional system (independent of havelock.blocks).

F(u) = H(z(u)) - Omega J3(z(u)), z_k = ring_k exp(s_k + i theta_k), u = (s, theta). Omega is the momentum
multiplier, solved from the brute-force gradient (dH = Omega dJ3 along s_0); the test then checks that dF = 0
in all 2N directions, so the augmented Hamiltonian is critical at the ring. The Hessian is formed by central
differences with step h in mpmath, and its eigenvalues (w.r.t. sum ds^2 + dtheta^2) are returned sorted."""
import mpmath as mp


def fd_spectrum(model, N, v, h, omega_shift=0):
    v = mp.mpf(v)
    rho = model.chart_radius(v)
    ring = [rho * mp.expjpi(mp.mpf(2 * k) / N) for k in range(N)]
    dim = 2 * N

    def Z(u):
        return [ring[k] * mp.exp(u[k] + 1j * u[N + k]) for k in range(N)]

    def shift(u, i, s):
        u = list(u)
        u[i] += s
        return u

    zero = [mp.mpf(0)] * dim
    H = lambda u: model.H_mp(Z(u), v)
    J3 = lambda u: model.J_mp(Z(u), v)[2]
    omega = (H(shift(zero, 0, h)) - H(shift(zero, 0, -h))) / (J3(shift(zero, 0, h)) - J3(shift(zero, 0, -h)))
    omega += omega_shift
    F = lambda u: H(u) - omega * J3(u)
    grad = max(abs(F(shift(zero, i, h)) - F(shift(zero, i, -h))) / (2 * h) for i in range(dim))
    M = mp.matrix(dim, dim)
    f0 = F(zero)
    for i in range(dim):
        M[i, i] = (F(shift(zero, i, h)) - 2 * f0 + F(shift(zero, i, -h))) / h**2
        for j in range(i + 1, dim):
            pp = F(shift(shift(zero, i, h), j, h))
            pm = F(shift(shift(zero, i, h), j, -h))
            mp_ = F(shift(shift(zero, i, -h), j, h))
            mm = F(shift(shift(zero, i, -h), j, -h))
            M[i, j] = M[j, i] = (pp - pm - mp_ + mm) / (4 * h * h)
    E = mp.eigsy(M, eigvals_only=True)
    return sorted(E[k] for k in range(dim)), omega, grad


def block_spectrum(family, N, v):
    """Eigenvalues of all Fourier blocks m = 0..N-1 from havelock.blocks, sorted."""
    from havelock.blocks import blocks_mp

    out = []
    for m in range(N):
        a, b, g = blocks_mp(family, N, m, v)
        s = mp.sqrt(((a - b) / 2) ** 2 + g * g)
        out += [(a + b) / 2 + s, (a + b) / 2 - s]
    return sorted(out)
