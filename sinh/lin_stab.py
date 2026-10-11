"""Item 01: independent linear-stability computation for the sinh model (NUMERICAL).

Fresh script (does not reuse checks/).  Model, written out:
  phase space: N points z_k = x_k + i y_k in the unit disk, symplectic form
      w = sum_k  lam(z_k) dx_k ^ dy_k,   lam(z) = 4/(1-|z|^2)^2
  H = -sum_{j<k} ln|z_j - z_k| + (N-1)/2 sum_k ln(1-|z_k|^2)     (= -sum ln sinh(d_jk/2) + const)
  J = sum_k (1+|z_k|^2)/(1-|z_k|^2)   (generator of rotations about 0)
  Hamilton's equations: xdot_k = (1/lam) dF/dy_k,  ydot_k = -(1/lam) dF/dx_k,  F = H - Omega J.
Gradient of H and J is coded analytically in Cartesian coordinates; nothing below uses the
Fourier/log-polar derivation.  Omega is solved from the equilibrium condition numerically.
The Jacobian of the vector field at the ring is formed by central differences at 60-digit
precision (h = 1e-25), and its eigenvalues are computed with mpmath.eig.

Checks:
 (A) spectrum of the Jacobian vs. the predicted set {0,0} U {+-sqrt(-a_m b_m)/kappa : m=1..N-1},
     a_m = C1 - m(N-m)/2, b_m = m(N-m)/2, C1 = (N-1)(1+xi^2)/(1-xi)^2, kappa = 4 xi/(1-xi)^2;
 (B) threshold: bisection in xi on 'max Re(eigenvalue) > tol', compared with xi*(N).
"""
import mpmath as mp
mp.mp.dps = 60

def grad(z, N, Om):
    """dF/dx_k, dF/dy_k for F = H - Om*J (analytic)."""
    gx, gy = [], []
    for k in range(N):
        xk, yk = mp.re(z[k]), mp.im(z[k])
        sx = sy = mp.mpf(0)
        for j in range(N):
            if j == k: continue
            dx, dy = xk - mp.re(z[j]), yk - mp.im(z[j])
            d2 = dx*dx + dy*dy
            sx += -dx/d2; sy += -dy/d2
        rho2 = xk*xk + yk*yk
        # W = (N-1)/2 ln(1-rho2):  dW/dx = -(N-1) x/(1-rho2)
        sx += -(N-1)*xk/(1-rho2); sy += -(N-1)*yk/(1-rho2)
        # J_k = (1+rho2)/(1-rho2): dJ/dx = 4x/(1-rho2)^2
        sx -= Om*4*xk/(1-rho2)**2; sy -= Om*4*yk/(1-rho2)**2
        gx.append(sx); gy.append(sy)
    return gx, gy

def field(v, N, Om):
    z = [mp.mpc(v[k], v[N+k]) for k in range(N)]
    gx, gy = grad(z, N, Om)
    out = [None]*(2*N)
    for k in range(N):
        lam = 4/(1-abs(z[k])**2)**2
        out[k] = gy[k]/lam
        out[N+k] = -gx[k]/lam
    return out

def ring(N, xi):
    r = mp.sqrt(xi)
    z = [r*mp.expjpi(2*mp.mpf(k)/N) for k in range(N)]
    return [mp.re(w) for w in z] + [mp.im(w) for w in z]

def omega_eq(N, xi):
    """Solve dF/dx_0 = 0 at vortex 0 (on the positive real axis) for Omega (F linear in Omega)."""
    v = ring(N, xi); z = [mp.mpc(v[k], v[N+k]) for k in range(N)]
    g0 = grad(z, N, 0)[0][0]; g1 = grad(z, N, 1)[0][0]
    return g0/(g0 - g1)

def jac_eigs(N, xi):
    Om = omega_eq(N, xi)
    v = ring(N, xi)
    res = max(abs(c) for c in field(v, N, Om))
    h = mp.mpf(10)**-25
    M = mp.matrix(2*N, 2*N)
    for i in range(2*N):
        vp = list(v); vm = list(v); vp[i] += h; vm[i] -= h
        fp, fm = field(vp, N, Om), field(vm, N, Om)
        for r in range(2*N):
            M[r, i] = (fp[r] - fm[r])/(2*h)
    ev = mp.eig(M, left=False, right=False)
    return Om, res, ev

def predicted(N, xi):
    C1 = (N-1)*(1+xi**2)/(1-xi)**2
    kap = 4*xi/(1-xi)**2
    out = [mp.mpc(0), mp.mpc(0)]
    for m in range(1, N):
        a = C1 - mp.mpf(m*(N-m))/2; b = mp.mpf(m*(N-m))/2
        s2 = -a*b/kap**2
        s = mp.sqrt(s2) if s2 >= 0 else mp.mpc(0, mp.sqrt(-s2))
        out += [s, -s]
    return out

def match(ev, pr):
    """max distance in an optimal greedy matching of two multisets of complex numbers."""
    ev = list(ev); err = mp.mpf(0)
    for p in pr:
        i = min(range(len(ev)), key=lambda j: abs(ev[j]-p))
        err = max(err, abs(ev[i]-p)); ev.pop(i)
    return err

def xistar(N):
    q = mp.mpf(N-3)/4 if N % 2 else mp.mpf((N-2)**2)/(4*(N-1))
    return (mp.sqrt(q)-1)/(mp.sqrt(q)+1)

if __name__ == '__main__':
    print('(A) spectrum check: N, xi, Omega numeric vs formula, residual |X(ring)|, max |eig - predicted|')
    worst = 0
    for N in range(8, 16):
        for xi in [mp.mpf('0.02'), xistar(N)*mp.mpf('0.9'), xistar(N)*mp.mpf('1.1'), mp.mpf('0.5')]:
            Om, res, ev = jac_eigs(N, xi)
            Omf = -(N-1)*(1+xi)*(1-xi)/(8*xi)
            e = match(ev, predicted(N, xi))
            worst = max(worst, e)
            print(f'{N:3d} {mp.nstr(xi,8):>12} {mp.nstr(Om-Omf,3):>10} {mp.nstr(res,3):>10} {mp.nstr(e,3):>10}   maxRe={mp.nstr(max(mp.re(x) for x in ev),5)}')
    print('worst spectral mismatch:', mp.nstr(worst, 5))
    print()
    print('(B) threshold by bisection on max Re(eig) > 1e-15 (pure dynamics, no formula), vs xi*(N)')
    for N in range(8, 16):
        lo, hi = xistar(N)*mp.mpf('0.5'), xistar(N)*mp.mpf('1.5')
        def unstable(x):
            return max(mp.re(e) for e in jac_eigs(N, x)[2]) > mp.mpf(10)**-15
        assert unstable(lo) and not unstable(hi)
        for _ in range(60):
            mid = (lo+hi)/2
            if unstable(mid): lo = mid
            else: hi = mid
        print(f'{N:3d}  bisection xi in [{mp.nstr(lo,22)}, {mp.nstr(hi,22)}]   xi*(N) = {mp.nstr(xistar(N),22)}   diff = {mp.nstr((lo+hi)/2-xistar(N),3)}')
