"""Conj13 numerical engine (independent of the Fourier formulas).

Cartesian coordinates q = (x_1..x_n, y_1..y_n) on the uniform chart C of Montaldi-Tokieda (MT).
Models (overall factor 1/(4 pi) dropped everywhere; it does not affect signs):
  'MT'  : H = -sum_{i<j} log|z_i-z_j|^2 + (n-1) sum_j log(1+lam|z_j|^2)      MT (1.12) family;
          lam<0: sinh model on H^2 (curvature 4 lam), lam=0: plane, lam>0: sphere.
  'tanh': lam = -1 only. H = -sum_{i<j} [log|z_i-z_j|^2 - log|1 - z_i conj z_j|^2]  (= -2 sum log tanh(d/2)).
Momentum (both models): J = (sum x/s, sum y/s, sum |z|^2/s), s = 1+lam|z|^2 (MT Sec.1; constants dropped).
Symplectic form: sum 2 dx^dy / s^2 (MT); for lam=-1 this is (1/2) of 4dx^dy/(1-|z|^2)^2 -> time rescaling only.
Augmented Hamiltonian Hhat = H - omega*J3, omega from the ring equilibrium.

Effective (reduced, centre-direction) function:
  f_eff(eps) = min { H(q) - H(q0) : q = q0 + eps*e + v, v in span(F) + span(Nn), J(q) = J(q0) }
where e spans (part of) the kernel of the constrained Hessian, F = the other constrained-Hessian eigvectors,
Nn = momentum-normal directions; rotation is gauge-fixed by v _|_ t (t = rotation generator).
"""
import sympy as sp
import mpmath as mp

xi_, yi_, xj_, yj_, lam_s = sp.symbols('xi yi xj yj lam')
_vars = [xi_, yi_, xj_, yj_]

def _mk_pair(expr):
    g = [sp.diff(expr, v) for v in _vars]
    h = [[sp.diff(gi, v) for v in _vars] for gi in g]
    f0 = sp.lambdify(_vars + [lam_s], expr, 'mpmath')
    f1 = sp.lambdify(_vars + [lam_s], g, 'mpmath')
    f2 = sp.lambdify(_vars + [lam_s], h, 'mpmath')
    return f0, f1, f2

def _mk_one(expr):
    v2 = [xi_, yi_]
    g = [sp.diff(expr, v) for v in v2]
    h = [[sp.diff(gi, v) for v in v2] for gi in g]
    return (sp.lambdify(v2 + [lam_s], expr, 'mpmath'), sp.lambdify(v2 + [lam_s], g, 'mpmath'),
            sp.lambdify(v2 + [lam_s], h, 'mpmath'))

d2 = (xi_-xj_)**2 + (yi_-yj_)**2
PAIR_LOG = _mk_pair(-sp.log(d2))
# |1 - z_i conj z_j|^2 with lam: |1 + lam z_i conj z_j|^2 (MT 4.2 form); for lam=-1 this is the tanh image term
re = xi_*xj_ + yi_*yj_; im = yi_*xj_ - xi_*yj_
PAIR_IMG = _mk_pair(sp.log((1 + lam_s*re)**2 + (lam_s*im)**2))
r2 = xi_**2 + yi_**2
s_ = 1 + lam_s*r2
ONE_LOG = _mk_one(sp.log(s_))
J1 = _mk_one(xi_/s_); J2 = _mk_one(yi_/s_); J3 = _mk_one(r2/s_)


class Ring:
    def __init__(self, n, lam, r0, model='MT'):
        self.n, self.lam, self.r0, self.model = n, mp.mpf(lam), mp.mpf(r0), model
        self.q0 = self.ring_q()
        self.omega = mp.mpf(0)
        g = self.gradH(self.q0)
        gJ = self.gradJ(self.q0)[2]
        # stationarity along radial direction of vortex 0 (x_0 component)
        self.omega = g[0] / gJ[0]
        self.Jmu = self.J(self.q0)

    def ring_q(self):
        n = self.n
        z = [self.r0 * mp.expjpi(2*mp.mpf(k)/n) for k in range(n)]
        return [mp.re(w) for w in z] + [mp.im(w) for w in z]

    def _pairs(self):
        n = self.n
        for i in range(n):
            for j in range(i+1, n):
                yield i, j

    def H(self, q):
        n, lam = self.n, self.lam
        x, y = q[:n], q[n:]
        val = mp.mpf(0)
        for i, j in self._pairs():
            val += PAIR_LOG[0](x[i], y[i], x[j], y[j], lam)
            if self.model == 'tanh':
                val += PAIR_IMG[0](x[i], y[i], x[j], y[j], lam)
        if self.model == 'MT':
            for i in range(n):
                val += (n-1)*ONE_LOG[0](x[i], y[i], lam)
        return val

    def gradH(self, q, with_omega=True):
        n, lam = self.n, self.lam
        x, y = q[:n], q[n:]
        g = [mp.mpf(0)]*(2*n)
        def add_pair(F, i, j):
            d = F[1](x[i], y[i], x[j], y[j], lam)
            g[i] += d[0]; g[n+i] += d[1]; g[j] += d[2]; g[n+j] += d[3]
        for i, j in self._pairs():
            add_pair(PAIR_LOG, i, j)
            if self.model == 'tanh':
                add_pair(PAIR_IMG, i, j)
        for i in range(n):
            if self.model == 'MT':
                d = ONE_LOG[1](x[i], y[i], lam); g[i] += (n-1)*d[0]; g[n+i] += (n-1)*d[1]
            if with_omega:
                d = J3[1](x[i], y[i], lam); g[i] -= self.omega*d[0]; g[n+i] -= self.omega*d[1]
        return g

    def hessH(self, q):
        n, lam = self.n, self.lam
        x, y = q[:n], q[n:]
        Hm = mp.zeros(2*n, 2*n)
        idx = lambda i: (i, n+i)
        def add_pair(F, i, j):
            h = F[2](x[i], y[i], x[j], y[j], lam)
            ids = [i, n+i, j, n+j]
            for a in range(4):
                for c in range(4):
                    Hm[ids[a], ids[c]] += h[a][c]
        for i, j in self._pairs():
            add_pair(PAIR_LOG, i, j)
            if self.model == 'tanh':
                add_pair(PAIR_IMG, i, j)
        for i in range(n):
            ids = [i, n+i]
            if self.model == 'MT':
                h = ONE_LOG[2](x[i], y[i], lam)
                for a in range(2):
                    for c in range(2):
                        Hm[ids[a], ids[c]] += (n-1)*h[a][c]
            h = J3[2](x[i], y[i], lam)
            for a in range(2):
                for c in range(2):
                    Hm[ids[a], ids[c]] -= self.omega*h[a][c]
        return Hm

    def J(self, q):
        n, lam = self.n, self.lam
        x, y = q[:n], q[n:]
        return [sum(F[0](x[i], y[i], lam) for i in range(n)) for F in (J1, J2, J3)]

    def gradJ(self, q):
        n, lam = self.n, self.lam
        x, y = q[:n], q[n:]
        out = []
        for F in (J1, J2, J3):
            g = [mp.mpf(0)]*(2*n)
            for i in range(n):
                d = F[1](x[i], y[i], lam); g[i] = d[0]; g[n+i] = d[1]
            out.append(g)
        return out

    def rot_gen(self, q):
        n = self.n
        return [-q[n+i] for i in range(n)] + [q[i] for i in range(n)]

    # ---------- linear algebra helpers
    @staticmethod
    def _gs(vecs, basis=()):
        out = list(basis)
        new = []
        for v in vecs:
            w = mp.matrix(v)
            for u in out:
                w = w - (u.T*w)[0]*u
            nr = mp.norm(w)
            if nr > mp.mpf(10)**(-mp.mp.dps//2):
                w = w/nr; out.append(w); new.append(w)
        return new

    def constrained_hessian(self):
        """Hessian of Hhat restricted to T = ker dJ(q0) cap t^perp; returns (eigvals, eigvecs-in-R^2n, normals)."""
        n = self.n
        q0 = self.q0
        G = self.gradJ(q0)
        t = self.rot_gen(q0)
        # Normal space: span(G rows) and t; T = orthogonal complement
        Nbasis = self._gs([t] + G)
        tvec = Nbasis[0]
        normals = Nbasis[1:]
        std = [[mp.mpf(1) if k == i else mp.mpf(0) for k in range(2*n)] for i in range(2*n)]
        Tb = self._gs(std, basis=Nbasis)
        assert len(Tb) == 2*n - len(Nbasis)
        P = mp.matrix(2*n, len(Tb))
        for c, v in enumerate(Tb):
            for r in range(2*n):
                P[r, c] = v[r]
        # constrained (Lagrangian) Hessian: Hess Hhat + sum nu_k Hess J_k with nu from stationarity;
        # at the ring grad Hhat = 0 so nu = 0 except J3 already in Hhat.  J1,J2 multipliers vanish by symmetry.
        Hm = self.hessH(q0)
        B = P.T*Hm*P
        B = (B + B.T)/2
        E, Q = mp.eigsy(B)
        vecs = [P*Q[:, k] for k in range(len(Tb))]
        return [E[k] for k in range(len(Tb))], vecs, normals, tvec

    def f_eff(self, eps, e, F, normals, guess=None, tol=None):
        """Minimise H over q0 + eps e + span(F) + span(normals) subject to the momentum constraint."""
        n = self.n
        m, k = len(F), len(normals)
        G0 = self.gradJ(self.q0)
        tol = tol or mp.mpf(10)**(-mp.mp.dps + 8)
        X = guess if guess is not None else [mp.mpf(0)]*(m + 2*k)
        if guess is None:
            X[m+k+2] = self.omega   # multiplier of J3 at the ring
        def build(X):
            v = mp.matrix(self.q0) + eps*e
            for i in range(m):
                v += X[i]*F[i]
            for i in range(k):
                v += X[m+i]*normals[i]
            return [v[r] for r in range(2*n)]
        def resid(X):
            q = build(X)
            g = mp.matrix(self.gradH(q, with_omega=False))
            Jq = self.J(q)
            cons = [Jq[a] - self.Jmu[a] for a in range(3)]
            GJ = [mp.matrix(gg) for gg in self.gradJ(q)]
            nu = X[m+k:]
            # Lagrange: grad H - sum nu_a grad J_a  projected on the free directions
            L = g
            for a in range(3):
                L = L - nu[a]*GJ[a]
            R = [(F[i].T*L)[0] for i in range(m)] + [(normals[i].T*L)[0] for i in range(k)] + cons
            return R
        # Newton with numerical Jacobian (central differences)
        for it in range(40):
            R = resid(X)
            nr = max(abs(r) for r in R)
            if nr < tol:
                break
            h = mp.mpf(10)**(-mp.mp.dps//3)
            Jm = mp.matrix(len(R), len(X))
            for c in range(len(X)):
                Xp = list(X); Xm = list(X); Xp[c] += h; Xm[c] -= h
                Rp, Rm = resid(Xp), resid(Xm)
                for r in range(len(R)):
                    Jm[r, c] = (Rp[r]-Rm[r])/(2*h)
            dX = mp.lu_solve(Jm, mp.matrix(R))
            X = [X[i] - dX[i] for i in range(len(X))]
        else:
            raise RuntimeError(f"Newton failed, resid {nr}")
        q = build(X)
        val = self.H(q) - self.H(self.q0)
        return val, X


def fit_even_poly(eps_list, vals, powers):
    """Solve vals = sum_p c_p eps^p (square system), with eps rescaled by max|eps| for conditioning."""
    sc = max(abs(e) for e in eps_list)
    A = mp.matrix(len(eps_list), len(powers))
    for i, e in enumerate(eps_list):
        for j, p in enumerate(powers):
            A[i, j] = (e/sc)**p
    c = mp.lu_solve(A, mp.matrix(vals))
    return [c[j]/sc**p for j, p in enumerate(powers)]
