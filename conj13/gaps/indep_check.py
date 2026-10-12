"""Independent check of the spectral and nearby-momentum arguments (paper Prop. spectral, Lemmas B2-B3). Does not import the other code.

Model: MT (1.12) in the uniform chart, sigma = 1 + lam|z|^2, Omega = sum 2 dx^dy/sigma^2 (MT (1.8)),
H = -sum_{i<j} log(|zi-zj|^2/(sigma_i sigma_j)),  J = sum (x/sigma, y/sigma, |z|^2/sigma).
Poisson bracket {f,g} = sum sigma^2/2 (f_x g_y - f_y g_x).  Ring z_k = e^{2 pi i k/n}, lam = b_n.

Part 1 (mpmath, dps 50; analytic gradients, Hessian by central differences of the analytic gradient):
  (1a) {J1,J2} = n/2 - lam*J3 at RANDOM configurations (so for lam = 0 it is the constant n/2: the
       SE(2) momentum map is equivariant only for an affine (centrally extended) action; for lam != 0 the
       shift J3 -> J3 - n/(2 lam) makes it linear-equivariant).  Also {J3,J1} = -J2? {J3,J2} = J1?  printed.
  (1b) the 3x3 matrix dJ_i(X_{J_j})(ring): rank 2, kernel spanned by X_{J3} only (orbit cap ker dJ = R t).
  (1c) Cartesian (NOT log-polar) compression of d^2(H - w0 J3) onto ker dJ cap (X_J3)^perp: dimension,
       inertia (#neg, #zero), smallest positive eigenvalue (metric-dependent, inertia is not).
Part 2 (float64 + scipy SLSQP; NUMERICAL evidence only) Lemma B3 directly: minimise
       V - V0 = H - w0 J3 + Kc (J3 - mu3)^2 - V(ring)  over {J1 = J2 = 0, gauge t.(x-x0) = 0, |x-x0| = eps}
  from many random starts, for Kc = 0 and Kc = 5; print min and min/eps^4.  Also Lemma B2 directly:
  min over the same sphere of H subject to J3 = mu3 + j, compared with H(x_j) (concentric ring).
"""
import sys, random
import mpmath as mp
import numpy as np
from scipy.optimize import minimize

def b_n(n):
    T = mp.mpf((n * n) // 4) / (2 * (n - 1))
    # (1+b^2) = T (1+b)^2 ; pick root in (-1, 1]
    if n == 7:
        return mp.mpf(0)
    r = mp.polyroots([1 - T, -2 * T, 1 - T])
    r = [mp.re(x) for x in r if abs(mp.im(x)) < 1e-30 and -1 < mp.re(x) <= 1]
    assert len(r) == 1
    return r[0]

# ---------- generic (works for mp and float) ----------
def gradH(q, lam, n):
    g = [0] * (2 * n)
    for i in range(n):
        xi, yi = q[2 * i], q[2 * i + 1]
        si = 1 + lam * (xi * xi + yi * yi)
        gx = (n - 1) * 2 * lam * xi / si; gy = (n - 1) * 2 * lam * yi / si
        for j in range(n):
            if j == i: continue
            dx, dy = xi - q[2 * j], yi - q[2 * j + 1]; d2 = dx * dx + dy * dy
            gx -= 2 * dx / d2; gy -= 2 * dy / d2
        g[2 * i], g[2 * i + 1] = gx, gy
    return g

def Hval(q, lam, n, log):
    s = 0
    for i in range(n):
        for j in range(i + 1, n):
            xi, yi, xj, yj = q[2 * i], q[2 * i + 1], q[2 * j], q[2 * j + 1]
            si = 1 + lam * (xi * xi + yi * yi); sj = 1 + lam * (xj * xj + yj * yj)
            s -= log(((xi - xj) ** 2 + (yi - yj) ** 2) / (si * sj))
    return s

def Jval(q, lam, n):
    J = [0, 0, 0]
    for k in range(n):
        x, y = q[2 * k], q[2 * k + 1]; s = 1 + lam * (x * x + y * y)
        J[0] += x / s; J[1] += y / s; J[2] += (x * x + y * y) / s
    return J

def gradJ(q, lam, n):
    G = [[0] * (2 * n) for _ in range(3)]
    for k in range(n):
        x, y = q[2 * k], q[2 * k + 1]; s = 1 + lam * (x * x + y * y)
        G[0][2 * k] = (s - 2 * lam * x * x) / s ** 2; G[0][2 * k + 1] = -2 * lam * x * y / s ** 2
        G[1][2 * k] = -2 * lam * x * y / s ** 2;      G[1][2 * k + 1] = (s - 2 * lam * y * y) / s ** 2
        G[2][2 * k] = 2 * x / s ** 2;                  G[2][2 * k + 1] = 2 * y / s ** 2
    return G

def pb(gf, gg, q, lam, n):
    s = 0
    for k in range(n):
        x, y = q[2 * k], q[2 * k + 1]; sg = 1 + lam * (x * x + y * y)
        s += sg ** 2 / 2 * (gf[2 * k] * gg[2 * k + 1] - gf[2 * k + 1] * gg[2 * k])
    return s

def ham_vf(gf, q, lam, n):
    # X_f with Omega(X_f, .) = df :  Xx = (sigma^2/2) f_y, Xy = -(sigma^2/2) f_x  (sign irrelevant here)
    X = [0] * (2 * n)
    for k in range(n):
        x, y = q[2 * k], q[2 * k + 1]; sg = 1 + lam * (x * x + y * y)
        X[2 * k] = sg ** 2 / 2 * gf[2 * k + 1]; X[2 * k + 1] = -sg ** 2 / 2 * gf[2 * k]
    return X

def ring(n, one, cos, sin, pi):
    q = []
    for k in range(n):
        q += [cos(2 * pi * k / n) * one, sin(2 * pi * k / n) * one]
    return q

# ---------- Part 1 ----------
def part1(n):
    mp.mp.dps = 50
    lam = b_n(n); N = 2 * n
    random.seed(1000 + n)
    worst = 0
    for trial in range(3):
        qr = [mp.mpf(random.uniform(-0.5, 0.5)) for _ in range(N)]
        G = gradJ(qr, lam, n); Jr = Jval(qr, lam, n)
        e12 = pb(G[0], G[1], qr, lam, n) - (mp.mpf(n) / 2 - lam * Jr[2])
        e31 = pb(G[2], G[0], qr, lam, n); e32 = pb(G[2], G[1], qr, lam, n)
        # expect {J3,J1} = +-J2, {J3,J2} = -+J1
        e3 = min(abs(e31 - Jr[1]) + abs(e32 + Jr[0]), abs(e31 + Jr[1]) + abs(e32 - Jr[0]))
        eH = max(abs(pb(gradH(qr, lam, n), G[i], qr, lam, n)) for i in range(3))
        worst = max(worst, abs(e12), e3, eH)
    q0 = ring(n, mp.mpf(1), mp.cos, mp.sin, mp.pi)
    G0 = gradJ(q0, lam, n); gH0 = gradH(q0, lam, n)
    J12 = pb(G0[0], G0[1], q0, lam, n)
    X = [ham_vf(G0[i], q0, lam, n) for i in range(3)]
    M = mp.matrix(3, 3)
    for i in range(3):
        for j in range(3):
            M[i, j] = mp.fdot(G0[i], X[j])
    sv = mp.svd_r(M, compute_uv=False)
    # w0 from dH = w0 dJ3 (least squares), residual
    w0 = mp.fdot(gH0, G0[2]) / mp.fdot(G0[2], G0[2])
    res = max(abs(gH0[i] - w0 * G0[2][i]) for i in range(N))
    # Hessian of F by central differences of analytic gradient
    h = mp.mpf(10) ** -18
    def gF(q):
        gh = gradH(q, lam, n); gj = gradJ(q, lam, n)[2]
        return [gh[i] - w0 * gj[i] for i in range(N)]
    Hs = mp.matrix(N, N)
    for c in range(N):
        qp = q0[:]; qm = q0[:]; qp[c] += h; qm[c] -= h
        gp, gm = gF(qp), gF(qm)
        for r in range(N):
            Hs[r, c] = (gp[r] - gm[r]) / (2 * h)
    Hs = (Hs + Hs.T) / 2
    C = mp.matrix(4, N)
    for c in range(N):
        for i in range(3): C[i, c] = G0[i][c]
        C[3, c] = X[2][c]
    P = mp.eye(N) - C.T * mp.inverse(C * C.T) * C
    pe, pv = mp.eigsy(P)
    cols = [c for c in range(N) if abs(pe[c] - 1) < 1e-20]
    W = mp.matrix(N, len(cols))
    for j, c in enumerate(cols):
        for r in range(N): W[r, j] = pv[r, c]
    R = W.T * Hs * W; R = (R + R.T) / 2
    ev = sorted(mp.eigsy(R)[0])
    tol = mp.mpf(10) ** -25
    nz = sum(1 for e in ev if abs(e) < tol); nn = sum(1 for e in ev if e < -tol)
    mpos = min(e for e in ev if e > tol)
    pred12 = mp.mpf(n) * (1 - lam) / (2 * (1 + lam))
    print(f"n={n:2d} b={mp.nstr(lam,10):>14s} | (1a) max err[{{J1,J2}}=n/2-lam J3, so(3)/se(2)/sl2 brackets, {{H,J}}=0] = {mp.nstr(worst,3)}"
          f" | {{J1,J2}}(ring)={mp.nstr(J12,10)} pred {mp.nstr(pred12,10)}"
          f" | (1b) sv[dJ(X_J)]={[mp.nstr(s,6) for s in sv]} | |dF|={mp.nstr(res,3)}"
          f" | (1c) dim={len(cols)} (2n-4={N-4}) neg={nn} zero={nz} min+={mp.nstr(mpos,8)}")
    ok = (worst < 1e-35 and abs(J12 - pred12) < 1e-35 and len(cols) == N - 4 and nn == 0
          and nz == (1 if n % 2 == 0 else 2) and res < 1e-35
          and sum(1 for s in sv if abs(s) < 1e-30) == 1)
    return ok

# ---------- Part 2 ----------
def part2(n, eps_list=(0.08, 0.04, 0.02), Kc_list=(0.0, 5.0), starts=24):
    lam = float(b_n(n)); N = 2 * n
    q0 = np.array(ring(n, 1.0, np.cos, np.sin, np.pi))
    g0 = np.array(gradH(list(q0), lam, n)); GJ0 = np.array(gradJ(list(q0), lam, n))
    w0 = g0 @ GJ0[2] / (GJ0[2] @ GJ0[2])
    mu3 = Jval(list(q0), lam, n)[2]
    t = np.array(ham_vf(list(GJ0[2]), list(q0), lam, n))
    Hf = lambda q: Hval(list(q), lam, n, np.log)
    V0 = Hf(q0) - w0 * mu3
    rng = np.random.default_rng(n)
    out = []
    for Kc in Kc_list:
        for eps in eps_list:
            def V(q):
                J3 = Jval(list(q), lam, n)[2]
                return Hf(q) - w0 * J3 + Kc * (J3 - mu3) ** 2 - V0
            def dV(q):
                J3 = Jval(list(q), lam, n)[2]
                return np.array(gradH(list(q), lam, n)) + (-w0 + 2 * Kc * (J3 - mu3)) * np.array(gradJ(list(q), lam, n)[2])
            cons = [{'type': 'eq', 'fun': lambda q: Jval(list(q), lam, n)[0], 'jac': lambda q: np.array(gradJ(list(q), lam, n)[0])},
                    {'type': 'eq', 'fun': lambda q: Jval(list(q), lam, n)[1], 'jac': lambda q: np.array(gradJ(list(q), lam, n)[1])},
                    {'type': 'eq', 'fun': lambda q: t @ (q - q0), 'jac': lambda q: t},
                    {'type': 'eq', 'fun': lambda q: (q - q0) @ (q - q0) - eps ** 2, 'jac': lambda q: 2 * (q - q0)}]
            best = np.inf
            for s in range(starts):
                d = rng.standard_normal(N); d /= np.linalg.norm(d)
                r = minimize(V, q0 + eps * d, jac=dV, constraints=cons, method='SLSQP',
                             options={'ftol': 1e-16, 'maxiter': 500})
                cviol = max(abs(c['fun'](r.x)) for c in cons)
                if cviol < 1e-10:
                    best = min(best, r.fun)
            out.append((Kc, eps, best))
            print(f"   n={n} Kc={Kc:4.1f} eps={eps:5.3f}: min_(S, gauge, |dx|=eps) (V-V0) = {best: .4e}   /eps^4 = {best/eps**4: .4f}")
    return out

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('1', 'all'):
        ok = True
        for n in range(4, 12):
            ok = part1(n) and ok
        print("PART 1:", "PASS" if ok else "FAIL")
    if which in ('2', 'all'):
        for n in (4, 5, 6, 7, 8, 9):
            part2(n)
