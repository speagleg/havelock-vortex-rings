"""Conj13: direct integrations near threshold (NUMERICAL).  Independent of engine.py and the Fourier formulas.
Poincare disk, curvature -1, area form 4 dx^dy/(1-|z|^2)^2.  Item normalisation:
  sinh: H = -sum_{i<j} ln|z_i-z_j| + (N-1)/2 sum ln(1-|z|^2)
  tanh: H = -sum_{i<j} [ln|z_i-z_j| - ln|1 - z_i conj z_j|]
Equations: zdot_k = (2i/lam_k) dH/dconj(z_k) with lam = 4/(1-|z|^2)^2, i.e. xdot = H_y/lam, ydot = -H_x/lam.
Initial data: ring at xi = r_E^2 plus a perturbation along the critical mode N/2 (radial, alternating
sign) of relative amplitude delta:  |z_k| = sqrt(xi) (1 + delta (-1)^k).
Diagnostic: shape deviation D(t) = max_{i<j} |d_ij(t) - d_ij(ring)| (hyperbolic distances; invariant under
the whole symmetry group, so no frame/drift issues).  We report max_t D / D(0).
Conservation of H and of the momentum (J0,J1,J2) is monitored.
"""
import sys, numpy as np
from scipy.integrate import solve_ivp

def make(N, model):
    iu = np.triu_indices(N, 1)
    def grad(z):
        dz = z[:, None] - z[None, :]
        np.fill_diagonal(dz, 1.0)
        inv = 1.0/np.conj(dz)          # d/dconj z_k of ln|z_k - z_j| = 1/(2 conj(z_k - z_j))
        np.fill_diagonal(inv, 0.0)
        G = -0.5*inv.sum(axis=1)       # dH/dconj z_k of -sum ln|z_k-z_j|
        r2 = np.abs(z)**2
        if model == 'sinh':
            G += (N-1)/2 * (-z/(1-r2))          # d/dconj z of (N-1)/2 ln(1-|z|^2)
        else:
            # + sum_{j} ln|1 - z_k conj z_j| : d/dconj z_k of ln|1 - z_k conj z_j|  (k as first slot)
            # ln|1 - a conj b| = (1/2)[ln(1-a conj b) + ln(1 - conj a b)];  d/dconj a = (1/2)(-b)/(1-conj a b)
            P = 1 - np.conj(z)[:, None]*z[None, :]
            T = -0.5*z[None, :]/P
            # pair (k,j) counted once in H but the term is symmetric; d/dconj z_k gets the k-slot derivative
            # of ln|1-z_k conj z_j| for j != k (the j-slot term is the same function).
            np.fill_diagonal(T, 0.0)
            G += T.sum(axis=1)
        return G
    def rhs(t, y):
        z = y[:N] + 1j*y[N:]
        lam = 4/(1-np.abs(z)**2)**2
        zd = 2j*grad(z)/lam
        return np.concatenate([zd.real, zd.imag])
    def H(z):
        dz = np.abs(z[:, None]-z[None, :])[iu]
        h = -np.log(dz).sum()
        if model == 'sinh':
            h += (N-1)/2*np.log(1-np.abs(z)**2).sum()
        else:
            h += np.log(np.abs(1 - z[:, None]*np.conj(z)[None, :])[iu]).sum()
        return h
    def J(z):
        r2 = np.abs(z)**2
        return np.array([((1+r2)/(1-r2)).sum(), (2*z.real/(1-r2)).sum(), (2*z.imag/(1-r2)).sum()])
    def dists(z):
        num = 2*np.abs(z[:, None]-z[None, :])**2
        den = (1-np.abs(z)**2)[:, None]*(1-np.abs(z)**2)[None, :]
        return np.arccosh(1 + num/den)[iu]
    return rhs, H, J, dists

def run(N, model, xi, delta, T, rtol=1e-12):
    rhs, H, J, dists = make(N, model)
    k = np.arange(N)
    ring = np.sqrt(xi)*np.exp(2j*np.pi*k/N)
    z0 = ring*(1 + delta*(-1.0)**k)
    y0 = np.concatenate([z0.real, z0.imag])
    d_ring = dists(ring)
    ts = np.linspace(0, T, 4001)
    sol = solve_ivp(rhs, (0, T), y0, method='DOP853', rtol=rtol, atol=rtol*1e-2, t_eval=ts)
    Z = sol.y[:N].T + 1j*sol.y[N:].T
    D = np.array([np.max(np.abs(dists(z) - d_ring)) for z in Z])
    Hs = np.array([H(z) for z in Z]); Js = np.array([J(z) for z in Z])
    # sanity: the unperturbed ring is a relative equilibrium (D stays ~0)
    return D, Hs, Js, sol

if __name__ == '__main__':
    # usage: integrate.py N model T rel delta     (xi = xi*(1+rel); rel>0 is the linearly stable side)
    import mpmath as mp
    mp.mp.dps = 30
    N = int(sys.argv[1]); model = sys.argv[2]; T = float(sys.argv[3])
    rel = float(sys.argv[4]); delta = float(sys.argv[5])
    XS = {('sinh', 8): 8-3*mp.sqrt(7), ('tanh', 8): mp.mpf('0.06735622300408126226843747656623557439402511756974')}
    xs = float(XS[(model, N)])
    Dr, _, _, _ = run(N, model, xs, 0.0, 20.0)
    xi = xs*(1+rel)
    D, Hs, Js, sol = run(N, model, xi, delta, T)
    D0 = D[0]
    q = len(D)//4
    print(f"N={N} {model} xi*={xs:.16g} rel={rel:+.0e} delta={delta:.0e} T={T:g}: unperturbed-ring drift {Dr.max():.1e}; "
          f"max D/D0 = {D.max()/D0:.3f}; per quarter: {[round(float(D[k*q:(k+1)*q+1].max()/D0), 3) for k in range(4)]}; "
          f"|dH|={np.abs(Hs-Hs[0]).max():.1e} |dJ|={np.abs(Js-Js[0]).max():.1e} nfev={sol.nfev}")
    np.save(f"out/D_{model}{N}_{rel:+.0e}_{delta:.0e}.npy", np.vstack([sol.t, D/D0]))
