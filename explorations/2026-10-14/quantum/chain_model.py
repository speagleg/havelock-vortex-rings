"""Pinned quantised hypercycle row as a 1D quantum chain (stiff-pinning limit), and its compression to a compact MPO.
Each site j carries one quantum coordinate y_j = Y_j - Y0 (Y = -cot phi), conjugate to the pinned along-row
displacement x_j ([x_j, Y_j] = i eps). With pinning stiffness k, the stiff-pinning Hamiltonian is
    H = sum_j (k eps^2/2)(-d^2/dy_j^2) + sum_{i<j} W(|i-j| l'; Y0+y_i, Y0+y_j) - V sum_j y_j,
where V is the row velocity (the counter-term that keeps the row at rest in the pinning lattice; it is the Lagrange
multiplier of the impulse in the Landau analysis). Local space: sinc-DVR grid of d points on [-ymax, ymax].
Pair interactions: common orthonormal basis A_s(y) (K functions, from the stacked G_r, r = 1..R), coefficient
matrices c(r) = A^T G_r A, fitted as c(r) = sum_k lambda_k^r M^(k) on a fixed grid of decay factors; each M^(k) is
eigendecomposed, so the MPO has one exponentially decaying coupling per (k, eigenvector).
Usage: python chain_model.py (prints compression errors)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows'))
import numpy as np

def W(du, pa, pb):
    ch = np.cosh(du); return 0.5*np.log((ch - np.cos(pa + pb))/(ch - np.cos(pa - pb)))

def build(kappa, l, d=16, ymax=1.2, R=None, K=6, nlam=12, tol=1e-9):
    s = np.sqrt(1 - kappa**2); lp = l*s; Y0 = -kappa/s
    y = np.linspace(-ymax, ymax, d); ph = np.pi/2 + np.arctan(Y0 + y)
    if R is None: R = int(np.ceil(36/lp))
    rs = np.arange(1, R + 1)
    G = np.array([W(r*lp, ph[:, None], ph[None, :]) for r in rs])          # (R, d, d)
    # velocity V: d/dY of the row energy per vortex = sum_r 2 * dW/dY_i at the row (both neighbours)
    h = 1e-6; ph_p = np.pi/2 + np.arctan(Y0 + h); ph_m = np.pi/2 + np.arctan(Y0 - h); ph0 = np.pi/2 + np.arctan(Y0)
    V = sum(2*(W(r*lp, ph_p, ph0) - W(r*lp, ph_m, ph0))/(2*h) for r in rs)
    # common basis: left singular vectors of the stacked [G_1 | G_2 | ...] weighted to emphasise accuracy
    stack = np.concatenate(list(G), axis=1)
    U, S, _ = np.linalg.svd(stack, full_matrices=False)
    A = U[:, :K]
    c = np.einsum('gs,rgh,ht->rst', A, G, A)                                 # (R, K, K)
    recon = np.einsum('gs,rst,ht->rgh', A, c, A)
    basis_err = np.abs(recon - G).max()
    # exponential fit of every entry of c(r): c(r) ~ sum_k lambda_k^r M_k   (least squares on a fixed lambda grid)
    lam = np.exp(-lp*np.geomspace(0.98, 40, nlam))                            # slowest ~ e^{-l'} (the far tail)
    X = lam[None, :]**rs[:, None]                                             # (R, nlam)
    coef, *_ = np.linalg.lstsq(X, c.reshape(R, -1), rcond=None)
    M = coef.reshape(nlam, K, K); M = 0.5*(M + M.transpose(0, 2, 1))
    cfit = np.einsum('rk,kst->rst', X, M)
    fit_err = np.abs(np.einsum('gs,rst,ht->rgh', A, cfit, A) - G).max()
    terms = []
    for k in range(nlam):
        w, v = np.linalg.eigh(M[k])
        for m in range(K):
            if abs(w[m]) > tol: terms.append((lam[k], w[m], A @ v[:, m]))       # coupling w * O(y_i) O(y_j) lam^(j-i)
    # total-energy check: classical energy per site of the uniform and zigzag configurations from the compressed form
    return dict(y=y, lp=lp, Y0=Y0, V=V, terms=terms, R=R, basis_err=basis_err, fit_err=fit_err, G=G, rs=rs, lam=lam)

if __name__ == '__main__':
    for kappa in [0.98, 0.97]:
        from exact_symbol import lstar
        m = build(kappa, lstar(kappa))
        print(f'kappa {kappa}: R = {m["R"]}, basis error {m["basis_err"]:.2e}, exponential-fit error {m["fit_err"]:.2e},'
              f' MPO couplings {len(m["terms"])}, V = {m["V"]:.6f}')
