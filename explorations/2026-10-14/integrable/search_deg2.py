"""Linear-algebra search for a local degree-2 integral of the open Toda-limit vortex chain (n sites)
    H = sum_j 2 s_j s_{j+1} e^{-D_j},  D_j = u_{j+1} - u_j,  s = (1+Y^2)^{-1/2},  {u_i, Y_k} = delta_ik.
Ansatz: K = sum_j s_j s_{j+1}^2 s_{j+2} P(Y_j, Y_{j+1}, Y_{j+2}) e^{-D_j - D_{j+1}} + sum_j s_j^2 s_{j+1}^2 Q(Y_j, Y_{j+1}) e^{-2 D_j},
P, Q polynomials of degree <= dmax in each variable. Rows: {K_basis, H} at random points (complex-step gradients).
Null vectors = local integrals of degree 2. Usage: python search_deg2.py [n] [dmax]"""
import sys, itertools, numpy as np
n = int(sys.argv[1]) if len(sys.argv) > 1 else 6; dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 2
rng = np.random.default_rng(1)
npts = 900
U = rng.uniform(-1, 1, (npts, n)) + 0.8*np.arange(n); Ys = rng.uniform(-2, 2, (npts, n))
def s(Y): return 1/np.sqrt(1 + Y**2)
def gradH(u, Y):
    S = s(Y); D = u[:, 1:] - u[:, :-1]; w = 2*S[:, :-1]*S[:, 1:]*np.exp(-D)
    gu = np.zeros_like(u); gu[:, :-1] += w; gu[:, 1:] -= w
    dl = -Y/(1 + Y**2); gY = np.zeros_like(Y); gY[:, :-1] += w*dl[:, :-1]; gY[:, 1:] += w*dl[:, 1:]
    return gu, gY
# basis
monoP = list(itertools.product(range(dmax + 1), repeat=3)); monoQ = list(itertools.product(range(dmax + 1), repeat=2))
def term(kind, mono, u, Y):
    S = s(Y); tot = 0
    if kind == 'P':
        i, j, k = mono
        for a in range(n - 2):
            tot = tot + S[:, a]*S[:, a+1]**2*S[:, a+2]*Y[:, a]**i*Y[:, a+1]**j*Y[:, a+2]**k*np.exp(-(u[:, a+2] - u[:, a]))
    else:
        i, j = mono
        for a in range(n - 1):
            tot = tot + S[:, a]**2*S[:, a+1]**2*Y[:, a]**i*Y[:, a+1]**j*np.exp(-2*(u[:, a+1] - u[:, a]))
    return tot
basis = [('P', m) for m in monoP] + [('Q', m) for m in monoQ]
gu, gY = gradH(U, Ys)
cols = []
h = 1e-30
for kind, mono in basis:
    col = np.zeros(npts)
    for c in range(n):
        Uc = U.astype(complex); Uc[:, c] += 1j*h; dFdu = term(kind, mono, Uc, Ys.astype(complex)).imag/h
        Yc = Ys.astype(complex); Yc[:, c] += 1j*h; dFdY = term(kind, mono, U.astype(complex), Yc).imag/h
        col += dFdu*gY[:, c] - dFdY*gu[:, c]
    cols.append(col)
A = np.array(cols).T
A /= np.linalg.norm(A, axis=0, keepdims=True) + 1e-300
sv = np.linalg.svd(A, compute_uv=False)
print(f'n = {n}, dmax = {dmax}: {len(basis)} unknowns; smallest singular values: {np.array2string(sv[-6:], precision=2)}; largest {sv[0]:.2e}')
_, _, Vt = np.linalg.svd(A)
for r in range(1, 4):
    v = Vt[-r]
    if sv[-r] < 1e-8*sv[0]:
        big = np.argsort(-abs(v))[:12]
        print('  null vector', r, ':', [(basis[i], round(float(v[i]/v[big[0]]), 6)) for i in big if abs(v[i]) > 1e-6*abs(v[big[0]])])
