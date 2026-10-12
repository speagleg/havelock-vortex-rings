"""Independent point-vortex ring code (written from scratch, no havelock imports).

Poincare disk, curvature -1: area form sigma = 4/(1-|z|^2)^2 dx^dy.
Models (identical unit vortices, overall constants irrelevant to stability):
  tanh : H = -sum_{i<j} [ ln|z_i-z_j|^2 - ln|1-z_i conj(z_j)|^2 ]          (tanh^2(d/2))
  sinh : H = -sum_{i<j} [ ln|z_i-z_j|^2 - ln(1-|z_i|^2) - ln(1-|z_j|^2) ]   (sinh^2(d/2))
Momentum J = sum 2/(1-|z_i|^2)  (= sum cosh r_i + const), the rotation moment map.
Ring: z_k = rho * exp(2 pi i k/N);  xi = tanh^2(r/2) = rho^2.
Hessians are analytic (Wirtinger, treating z and zbar as independent).
"""
import mpmath as mp


def complex_derivs(Z, model, coefH=1, coefJ=0):
    """Gradient and Hessian of coefH*H + coefJ*J w.r.t. V=(z_1..z_N, zb_1..zb_N)."""
    N = len(Z)
    n2 = 2 * N
    g = [mp.mpc(0)] * n2
    Hs = mp.matrix(n2, n2)
    Zb = [mp.conj(z) for z in Z]
    val = {}

    def var(k):
        return Z[k] if k < N else Zb[k - N]

    def add_log(c, idx, Lval, grad, hess):
        # c * ln L, L depends on variables idx with gradient grad (list) and hessian hess (dict)
        for a, ga in zip(idx, grad):
            g[a] += c * ga / Lval
        for i1, a in enumerate(idx):
            for i2, b in enumerate(idx):
                h = hess.get((i1, i2), 0)
                Hs[a, b] += c * (h / Lval - grad[i1] * grad[i2] / Lval ** 2)

    def add_inv(c, idx, Lval, grad, hess):
        # c / L
        for a, ga in zip(idx, grad):
            g[a] += -c * ga / Lval ** 2
        for i1, a in enumerate(idx):
            for i2, b in enumerate(idx):
                h = hess.get((i1, i2), 0)
                Hs[a, b] += c * (-h / Lval ** 2 + 2 * grad[i1] * grad[i2] / Lval ** 3)

    for i in range(N):
        for j in range(i + 1, N):
            # -ln(z_i - z_j) - ln(zb_i - zb_j)
            add_log(-coefH, [i, j], Z[i] - Z[j], [1, -1], {})
            add_log(-coefH, [N + i, N + j], Zb[i] - Zb[j], [1, -1], {})
            if model == "tanh":
                # + ln(1 - z_i zb_j) + ln(1 - zb_i z_j)
                add_log(coefH, [i, N + j], 1 - Z[i] * Zb[j], [-Zb[j], -Z[i]], {(0, 1): -1, (1, 0): -1})
                add_log(coefH, [N + i, j], 1 - Zb[i] * Z[j], [-Z[j], -Zb[i]], {(0, 1): -1, (1, 0): -1})
    if model == "sinh":
        # each vortex appears in N-1 pairs:  + (N-1) ln(1 - z zb)
        for i in range(N):
            add_log(coefH * (N - 1), [i, N + i], 1 - Z[i] * Zb[i], [-Zb[i], -Z[i]], {(0, 1): -1, (1, 0): -1})
    if coefJ:
        for i in range(N):
            add_inv(2 * coefJ, [i, N + i], 1 - Z[i] * Zb[i], [-Zb[i], -Z[i]], {(0, 1): -1, (1, 0): -1})
    return g, Hs


def real_derivs(Z, model, coefH=1, coefJ=0):
    """Real gradient/Hessian in coordinates (x_1,y_1,...,x_N,y_N)."""
    N = len(Z)
    g, Hs = complex_derivs(Z, model, coefH, coefJ)
    # d/dx_k = d/dz_k + d/dzb_k ; d/dy_k = i (d/dz_k - d/dzb_k)
    T = mp.matrix(2 * N, 2 * N)  # rows: complex var, cols: real var
    for k in range(N):
        T[k, 2 * k] = 1
        T[N + k, 2 * k] = 1
        T[k, 2 * k + 1] = 1j
        T[N + k, 2 * k + 1] = -1j
    gr = [mp.re(sum(T[a, c] * g[a] for a in range(2 * N))) for c in range(2 * N)]
    Hr = T.T * Hs * T
    Hr = mp.matrix([[mp.re(Hr[a, b]) for b in range(2 * N)] for a in range(2 * N)])
    return gr, Hr


def ring(N, rho):
    return [rho * mp.expjpi(mp.mpf(2 * k) / N) for k in range(N)]


def relative_equilibrium(N, rho, model):
    Z = ring(N, rho)
    gH, HH = real_derivs(Z, model, 1, 0)
    gJ, HJ = real_derivs(Z, model, 0, 1)
    Om = gH[0] / gJ[0]  # vortex 0 at (rho, 0): radial derivative = d/dx
    # sanity: full gradient of H - Om J vanishes
    res = max(abs(gH[c] - Om * gJ[c]) for c in range(2 * N))
    return Z, Om, HH - Om * HJ, gJ, res


def linearization(N, rho, model):
    Z, Om, Q, gJ, res = relative_equilibrium(N, rho, model)
    P = mp.matrix(2 * N, 2 * N)
    for k in range(N):
        s = 4 / (1 - abs(Z[k]) ** 2) ** 2
        # Poisson tensor = inverse of sigma*[[0,1],[-1,0]]
        P[2 * k, 2 * k + 1] = 1 / s
        P[2 * k + 1, 2 * k] = -1 / s
    return P * Q, Q, gJ, res


def max_real_part(N, xi, model):
    rho = mp.sqrt(xi)
    L, Q, gJ, res = linearization(N, rho, model)
    ev = mp.eig(L, left=False, right=False)
    return max(abs(mp.re(e)) for e in ev), res


def em_inertia(N, xi, model):
    """Inertia of Q = Hess(H - Om J) restricted to ker dJ (rotation direction gives one zero)."""
    rho = mp.sqrt(xi)
    L, Q, gJ, res = linearization(N, rho, model)
    n = 2 * N
    # basis of ker dJ: eliminate the coordinate with largest |gJ|
    p = max(range(n), key=lambda c: abs(gJ[c]))
    B = mp.matrix(n, n - 1)
    col = 0
    for c in range(n):
        if c == p:
            continue
        B[c, col] = 1
        B[p, col] = -gJ[c] / gJ[p]
        col += 1
    R = B.T * Q * B
    ev = mp.eigsy(R, eigvals_only=True)
    ev = sorted(ev, key=lambda e: e)
    return ev
