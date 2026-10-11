"""Gap (i) numerical check: full reduced Hessian at the bifurcation point b_n, MT family (1.12).

Independent of engine.py / c4_levelset.py: sympy differentiation, mpmath evaluation (dps 40).
Model: chart C, sigma = 1 + lam|z|^2, Omega_lam = sum 2 dx^dy / sigma^2 (MT (1.8)),
H = -sum_{i<j} log(|zi-zj|^2/(sigma_i sigma_j))   (MT H_lam without kappa^2/4pi),
J = sum (x/sigma, y/sigma, |z|^2/sigma)            (MT (1.9)).
Ring z_k = r0 e^{2 pi i k/n}, r0 = 1, lam = b_n (root of (3.6) in (-1,1]).

Checks, for n = 4..15:
 (a) {H, J_i} = 0 at a random configuration (J really is the momentum map for this omega).
 (b) {J1, J2}(ring) = n(1-b)/(2(1+b)) != 0  (slice S = {J1=J2=0} is symplectic at the ring).
 (c) dF(ring) = 0 for F = H - w0 J3 (w0 solved from the mode-0 radial direction).
 (d) Hessian of F in log-polar coords (s_k, theta_k) [= Jac^T Hess_cart Jac, valid at a critical point],
     compressed to W = ker dJ(ring) cap t^perp (Euclidean in (s,theta)); its spectrum is compared with
     the Fourier prediction
        radial    a_m = 2(n-1)(1+b^2)/(1+b)^2 - m(n-m),   m = 2..n-2
        tangential b_m = m(n-m),                           m = 2..n-2
        mode-1 piece (2-dim, ker dJ_+ cap V_1): (a_1 + b_1 rho^2)/(1+rho^2), rho = (1-b)/(1+b), double.
     and the inertia (zero count / negative count) is reported.  Inertia is metric-independent.
"""
import sympy as sp, mpmath as mp, random
mp.mp.dps = 40

def b_n(n):
    T = mp.mpf((n*n)//4)/(2*(n-1))          # (1+b^2)/(1+b)^2 = T
    # (1+b^2) = T(1+b)^2 -> (1-T) b^2 - 2T b + (1-T) = 0
    if n == 7: return mp.mpf(0)
    A, B, C = 1-T, -2*T, 1-T
    roots = [(-B + s*mp.sqrt(B*B-4*A*C))/(2*A) for s in (1, -1)]
    r = [x for x in roots if -1 < x <= 1]
    assert len(r) == 1, (n, roots)
    return r[0]

def build(n):
    X = sp.symbols('x0:%d' % n, real=True); Y = sp.symbols('y0:%d' % n, real=True); lam = sp.Symbol('lam')
    sig = [1 + lam*(X[k]**2 + Y[k]**2) for k in range(n)]
    H = 0
    for i in range(n):
        for j in range(i+1, n):
            H += -sp.log(((X[i]-X[j])**2 + (Y[i]-Y[j])**2)/(sig[i]*sig[j]))
    J = [sum(X[k]/sig[k] for k in range(n)), sum(Y[k]/sig[k] for k in range(n)),
         sum((X[k]**2+Y[k]**2)/sig[k] for k in range(n))]
    V = []
    for k in range(n): V += [X[k], Y[k]]
    gH = [sp.diff(H, v) for v in V]
    hH = [[sp.diff(g, v) for v in V] for g in gH]
    gJ = [[sp.diff(Ji, v) for v in V] for Ji in J]
    hJ3 = [[sp.diff(g, v) for v in V] for g in gJ[2]]
    f = lambda e: sp.lambdify([V, lam], e, 'mpmath')
    return dict(V=V, sig=sig, lam=lam, gH=f(gH), hH=f(hH), gJ=f(gJ), hJ3=f(hJ3), J=f(J))

def poisson(gF, gG, q, lam, n):
    s = mp.mpf(0)
    for k in range(n):
        x, y = q[2*k], q[2*k+1]; sg = 1 + lam*(x*x+y*y)
        s += sg**2/2*(gF[2*k]*gG[2*k+1] - gF[2*k+1]*gG[2*k])
    return s

def run(n, B):
    lam = b_n(n); N = 2*n
    q0 = []
    for k in range(n):
        ph = 2*mp.pi*k/n; q0 += [mp.cos(ph), mp.sin(ph)]
    # (a) momentum check at random config
    random.seed(n); qr = [mp.mpf(random.uniform(-0.6, 0.6)) for _ in range(N)]
    gHr = B['gH'](qr, lam); gJr = B['gJ'](qr, lam)
    pa = max(abs(poisson(gHr, gJr[i], qr, lam, n)) for i in range(3))
    # (b) {J1,J2} at ring
    gJ0 = [[mp.mpf(v) for v in row] for row in B['gJ'](q0, lam)]
    pb = poisson(gJ0[0], gJ0[1], q0, lam, n); pb_pred = n*(1-lam)/(2*(1+lam))
    # (c) critical point
    gH0 = [mp.mpf(v) for v in B['gH'](q0, lam)]
    rad = q0[:]                                 # d/ds of all radii at once = position vector
    w0 = mp.fdot(gH0, rad)/mp.fdot(gJ0[2], rad)
    gF = [gH0[i] - w0*gJ0[2][i] for i in range(N)]
    pc = max(abs(v) for v in gF)
    # (d) Hessian of F in (s,theta): x = e^{s} cos(phi+th), y = e^{s} sin(phi+th) (r0 = 1)
    HF = mp.matrix(B['hH'](q0, lam)) - w0*mp.matrix(B['hJ3'](q0, lam))
    Jac = mp.zeros(N, N)                         # columns: s_0..s_{n-1}, th_0..th_{n-1}
    for k in range(n):
        x, y = q0[2*k], q0[2*k+1]
        Jac[2*k, k] = x; Jac[2*k+1, k] = y       # d/ds
        Jac[2*k, n+k] = -y; Jac[2*k+1, n+k] = x  # d/dtheta
    Hst = Jac.T*HF*Jac
    dJst = mp.matrix([[mp.fdot([gJ0[i][r] for r in range(N)], [Jac[r, c] for r in range(N)]) for c in range(N)] for i in range(3)])
    t = mp.matrix([0]*n + [1]*n)                # rotation generator in (s,theta)
    # orthonormal basis of W = ker dJst cap t^perp
    C = mp.zeros(4, N)
    for i in range(3):
        for c in range(N): C[i, c] = dJst[i, c]
    for c in range(N): C[3, c] = t[c]
    P = mp.eye(N) - C.T*mp.inverse(C*C.T)*C      # orthogonal projector onto W
    pe, pv = mp.eigsy((P+P.T)/2)
    cols = [c for c in range(N) if abs(pe[c]-1) < mp.mpf(10)**(-25)]
    assert len(cols) == N-4, len(cols)
    Wb = mp.matrix(N, N-4)
    for j, c in enumerate(cols):
        for r in range(N): Wb[r, j] = pv[r, c]
    R = Wb.T*Hst*Wb
    R = (R + R.T)/2
    ev = sorted(mp.eigsy(R)[0])
    # prediction
    C1 = 2*(n-1)*(1+lam**2)/(1+lam)**2
    pred = [C1 - m*(n-m) for m in range(2, n-1)] + [mp.mpf(m*(n-m)) for m in range(2, n-1)]
    rho = (1-lam)/(1+lam); a1, b1 = C1-(n-1), mp.mpf(n-1)
    pred += [(a1 + b1*rho**2)/(1+rho**2)]*2
    pred = sorted(pred)
    dev = max(abs(ev[i]-pred[i]) for i in range(len(ev)))
    tol = mp.mpf(10)**(-20)
    nzero = sum(1 for e in ev if abs(e) < tol); nneg = sum(1 for e in ev if e < -tol)
    pos = [e for e in ev if e > tol]
    print(f"n={n:2d} b={mp.nstr(lam,12):>16s} |{{H,J}}|={mp.nstr(pa,3):>8s} {{J1,J2}}={mp.nstr(pb,10)} "
          f"(pred {mp.nstr(pb_pred,10)}) |dF|={mp.nstr(pc,3)} dim W={N-4} zero={nzero} neg={nneg} "
          f"min+={mp.nstr(min(pos),10)} max|ev-pred|={mp.nstr(dev,3)} radial mode-1 a_1={mp.nstr(a1,8)}")
    expect_zero = 1 if n % 2 == 0 else 2
    return nzero == expect_zero and nneg == 0 and dev < tol and pa < tol and pc < tol and abs(pb-pb_pred) < tol

if __name__ == '__main__':
    ok = True
    for n in range(4, 16):
        B = build(n)
        ok = run(n, B) and ok
    print("ALL CHECKS", "PASS" if ok else "FAIL")
