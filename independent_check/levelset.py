"""Independent level-set quartic at the sinh-model bifurcation (MT family, lambda<0).

At xi* the reduced Hessian Q|ker dJ has, besides the rotation zero, a kernel (critical mode).
For x = x0 + s*k + w on the level set J = J0 we minimise F = H - Om0*J over w (complement of
span{k, rotation}) by Newton on the Lagrange system, and report E(s) = F(x(s)) - F(x0) and E/s^4.
c4 > 0  <=>  E(s)/s^4 -> positive constant for every kernel direction k.
"""
import os
from fractions import Fraction
import sys
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vort import real_derivs, relative_equilibrium, ring

model = sys.argv[1]
N = int(sys.argv[2])
_q = Fraction(sys.argv[3])  # cosh^2 r*, exact rational, e.g. 9/7
q = mp.mpf(_q.numerator) / _q.denominator  # sinh model: (N-3)/4 for odd N, (N-2)^2/(4(N-1)) for even N
mp.mp.dps = 40
ch = mp.sqrt(q)
xi = (ch - 1) / (ch + 1)
rho = mp.sqrt(xi)
n2 = 2 * N

Z0, Om, Q, gJ, res = relative_equilibrium(N, rho, model)
x0 = []
for z in Z0:
    x0 += [mp.re(z), mp.im(z)]


def to_Z(x):
    return [mp.mpc(x[2 * k], x[2 * k + 1]) for k in range(N)]


def energy(x, which):
    Z = to_Z(x)
    if which == "J":
        return sum(2 / (1 - abs(z) ** 2) for z in Z)
    H = mp.mpf(0)
    for i in range(N):
        for j in range(i + 1, N):
            if model == "tanh":
                H -= mp.log(abs(Z[i] - Z[j]) ** 2 / abs(1 - Z[i] * mp.conj(Z[j])) ** 2)
            else:
                H -= mp.log(abs(Z[i] - Z[j]) ** 2 / ((1 - abs(Z[i]) ** 2) * (1 - abs(Z[j]) ** 2)))
    return H


# orthonormal basis of ker dJ, reduced Hessian, its near-zero eigenvectors
gJv = mp.matrix(gJ)
nJ = gJv / mp.norm(gJv)
# Gram-Schmidt: start from identity, remove nJ
basis = []
for c in range(n2):
    v = mp.matrix(n2, 1)
    v[c] = 1
    v = v - (nJ.T * v)[0] * nJ
    for b in basis:
        v = v - (b.T * v)[0] * b
    nv = mp.norm(v)
    if nv > mp.mpf("1e-10"):
        basis.append(v / nv)
    if len(basis) == n2 - 1:
        break
B = mp.matrix(n2, n2 - 1)
for j, b in enumerate(basis):
    for i in range(n2):
        B[i, j] = b[i]
R = B.T * Q * B
E, V = mp.eigsy(R)
order = sorted(range(n2 - 1), key=lambda i: abs(E[i]))
print(f"{model} N={N} xi*={mp.nstr(xi, 20)} RE residual {mp.nstr(res, 3)}")
print("smallest |eig| of reduced Hessian:", [mp.nstr(E[i], 5) for i in order[:5]])
print("signs of the rest: min", mp.nstr(min(E[i] for i in order[3:]), 5), "max", mp.nstr(max(E[i] for i in order[3:]), 5))
vrot = mp.matrix(n2, 1)
for k in range(N):
    vrot[2 * k] = -x0[2 * k + 1]
    vrot[2 * k + 1] = x0[2 * k]
vrot = vrot / mp.norm(vrot)
# kernel vectors: near-zero eigvecs orthogonal to rotation
ker = []
for i in order[:3]:
    v = B * V[:, i]
    v = v - (vrot.T * v)[0] * vrot
    for kv in ker:
        v = v - (kv.T * v)[0] * kv
    if mp.norm(v) > mp.mpf("0.5"):
        ker.append(v / mp.norm(v))
print("kernel dimension (excluding rotation):", len(ker))

F0 = energy(x0, "H") - Om * energy(x0, "J")
J0 = energy(x0, "J")


def solve(k, s):
    # complement W of span{k, vrot}
    W = []
    for c in range(n2):
        v = mp.matrix(n2, 1)
        v[c] = 1
        for b in [k, vrot] + W:
            v = v - (b.T * v)[0] * b
        if mp.norm(v) > mp.mpf("1e-8"):
            W.append(v / mp.norm(v))
        if len(W) == n2 - 2:
            break
    Wm = mp.matrix(n2, n2 - 2)
    for j, w in enumerate(W):
        for i in range(n2):
            Wm[i, j] = w[i]
    y = mp.matrix(n2 - 2, 1)
    mu = mp.mpf(0)
    for it in range(30):
        x = mp.matrix(x0) + s * k + Wm * y
        xl = [x[i] for i in range(n2)]
        Z = to_Z(xl)
        gH, HH = real_derivs(Z, model, 1, 0)
        gJx, HJ = real_derivs(Z, model, 0, 1)
        gF = mp.matrix([gH[i] - (Om + mu) * gJx[i] for i in range(n2)])
        r1 = Wm.T * gF
        r2 = energy(xl, "J") - J0
        Hm = HH - (Om + mu) * HJ
        A = mp.matrix(n2 - 1, n2 - 1)
        top = Wm.T * Hm * Wm
        cJ = Wm.T * mp.matrix(gJx)
        for i in range(n2 - 2):
            for j in range(n2 - 2):
                A[i, j] = top[i, j]
            A[i, n2 - 2] = -cJ[i]
            A[n2 - 2, i] = cJ[i]
        rhs = mp.matrix(n2 - 1, 1)
        for i in range(n2 - 2):
            rhs[i] = -r1[i]
        rhs[n2 - 2] = -r2
        d = mp.lu_solve(A, rhs)
        for i in range(n2 - 2):
            y[i] += d[i]
        mu += d[n2 - 2]
        if mp.norm(d) < mp.mpf(10) ** (-32):
            break
    return energy(xl, "H") - Om * energy(xl, "J") - F0, mp.norm(y)


for ang in [0, 1, 2, 3]:
    th = mp.pi * ang / 8
    k = ker[0] * mp.cos(th) + (ker[1] * mp.sin(th) if len(ker) > 1 else 0 * ker[0])
    k = k / mp.norm(k)
    out = []
    for s in ["2e-3", "-2e-3", "1e-3"]:
        s = mp.mpf(s)
        Es, ny = solve(k, s)
        out.append(mp.nstr(Es / s ** 4, 10))
    print(f"kernel angle {ang}/8 pi: E(s)/s^4 at s=+2e-3,-2e-3,+1e-3:", out, flush=True)
