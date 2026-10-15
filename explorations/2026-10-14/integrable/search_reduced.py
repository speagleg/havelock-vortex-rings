"""Search for a degree-2 integral of the reduced Toda-limit dynamics on a ring of n sites:
   dc_j/dt = 4 c_j (f_{j+1} c_{j+1} - f_j c_{j-1}),  dY_j/dt = 2 (c_{j-1} - c_j),  f = Y/(1+Y^2).
Ansatz K = sum_j c_j c_{j+1} P(Y_j,Y_{j+1},Y_{j+2}) + sum_j c_j^2 Q(Y_j,Y_{j+1}), with P, Q = monomials / products of
(1+Y^2)^a. Basis functions that are dependent as functions are removed first (column-pivoted QR on their values).
Usage: python search_reduced.py [dmax] [amax]"""
import sys, itertools, numpy as np
from scipy.linalg import qr
dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 2; amax = int(sys.argv[2]) if len(sys.argv) > 2 else 1
n = 7; npts = 2500; rng = np.random.default_rng(5)
C = rng.uniform(0.05, 1.0, (npts, n)); Yv = rng.uniform(-2.5, 2.5, (npts, n))
def rhs(c, Y):
    f = Y/(1 + Y**2); cp, cm = np.roll(c, -1, axis=1), np.roll(c, 1, axis=1); fp = np.roll(f, -1, axis=1)
    return 4*c*(fp*cp - f*cm), 2*(cm - c)
basis = []
for m in itertools.product(range(dmax + 1), repeat=3):
    for a in itertools.product(range(amax + 1), repeat=3): basis.append(('P', m, a))
for m in itertools.product(range(dmax + 1), repeat=2):
    for a in itertools.product(range(amax + 2), repeat=2): basis.append(('Q', m, a))
def K(b, c, Y):
    kind, m, a = b; tot = 0
    for j in range(n):
        if kind == 'P':
            Y0, Y1, Y2 = Y[:, j], Y[:, (j+1) % n], Y[:, (j+2) % n]
            tot = tot + c[:, j]*c[:, (j+1) % n]*Y0**m[0]*Y1**m[1]*Y2**m[2]/((1+Y0**2)**a[0]*(1+Y1**2)**a[1]*(1+Y2**2)**a[2])
        else:
            Y0, Y1 = Y[:, j], Y[:, (j+1) % n]
            tot = tot + c[:, j]**2*Y0**m[0]*Y1**m[1]/((1+Y0**2)**a[0]*(1+Y1**2)**a[1])
    return tot
dc, dY = rhs(C, Yv); h = 1e-30; cols = []; vals = []
for b in basis:
    col = np.zeros(npts)
    for k in range(n):
        Cc = C.astype(complex); Cc[:, k] += 1j*h; col += K(b, Cc, Yv.astype(complex)).imag/h*dc[:, k]
        Yc = Yv.astype(complex); Yc[:, k] += 1j*h; col += K(b, C.astype(complex), Yc).imag/h*dY[:, k]
    cols.append(col); vals.append(K(b, C, Yv))
A = np.array(cols).T; F = np.array(vals).T
Fn = F/np.linalg.norm(F, axis=0); _, R, piv = qr(Fn, pivoting=True, mode='economic')
rank = int(np.sum(np.abs(np.diag(R)) > 1e-10*abs(R[0, 0]))); keep = np.sort(piv[:rank])
A = A[:, keep]; B = [basis[i] for i in keep]; A = A/np.linalg.norm(A, axis=0)
U_, sv, Vt = np.linalg.svd(A, full_matrices=False)
print(f'dmax {dmax}, amax {amax}: basis {len(basis)} -> independent {rank}; smallest sv {np.array2string(sv[-5:], precision=2)}; largest {sv[0]:.2f}')
if sv[-1] < 1e-9*sv[0]:
    v = Vt[-1]; idx = np.argsort(-abs(v))
    print('  null vector:', [(B[i], round(float(v[i]/v[idx[0]]), 8)) for i in idx if abs(v[i]) > 1e-7*abs(v[idx[0]])])
