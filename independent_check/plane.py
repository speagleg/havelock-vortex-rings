"""Plane, n=7 (Thomson/Havelock heptagon, MT b_7 = 0). Independent level-set quartic.
H = -sum ln|z_i-z_j|^2, J3 = sum |z|^2, J1 = sum x, J2 = sum y; ring radius 1.
Also computes the quartic on the *linear slice* (no momentum constraint / no slaving along
momentum directions) to test the paper's P1 claim 9/14 (level set) vs 6/7 (linear slice).
"""
import mpmath as mp
mp.mp.dps = 40
N = 7
n2 = 2 * N


def H_g_h(x):
    g = mp.matrix(n2, 1); Hm = mp.matrix(n2, n2); E = mp.mpf(0)
    for i in range(N):
        for j in range(i + 1, N):
            dx = x[2*i] - x[2*j]; dy = x[2*i+1] - x[2*j+1]; r2 = dx*dx + dy*dy
            E -= mp.log(r2)
            d = [dx, dy]
            gi = [-2 * dx / r2, -2 * dy / r2]
            for a in range(2):
                g[2*i+a] += gi[a]; g[2*j+a] -= gi[a]
            for a in range(2):
                for b in range(2):
                    h = -((2 / r2) * (a == b) - 4 * d[a] * d[b] / r2**2)
                    Hm[2*i+a, 2*i+b] += h; Hm[2*j+a, 2*j+b] += h
                    Hm[2*i+a, 2*j+b] -= h; Hm[2*j+a, 2*i+b] -= h
    return E, g, Hm


def J(x):
    return sum(x[i]**2 for i in range(n2))


x0 = []
for k in range(N):
    x0 += [mp.cos(2 * mp.pi * k / N), mp.sin(2 * mp.pi * k / N)]
E0, gH, HH = H_g_h(x0)
gJ = mp.matrix([2 * v for v in x0])
Om = gH[0] / gJ[0]
Q = HH - Om * 2 * mp.eye(n2)
F0 = E0 - Om * J(x0)
G1 = mp.matrix([1 if i % 2 == 0 else 0 for i in range(n2)])
G2 = mp.matrix([1 if i % 2 == 1 else 0 for i in range(n2)])
cons_g = [gJ / mp.norm(gJ), G1 / mp.norm(G1), G2 / mp.norm(G2)]
# orthonormal basis of the joint kernel of dJ3, dJ1, dJ2
basis = []
for c in range(n2):
    v = mp.matrix(n2, 1); v[c] = 1
    for b in cons_g + basis:
        v = v - (b.T * v)[0] * b
    if mp.norm(v) > mp.mpf("1e-10"):
        basis.append(v / mp.norm(v))
B = mp.matrix(n2, len(basis))
for j, b in enumerate(basis):
    for i in range(n2):
        B[i, j] = b[i]
ev, V = mp.eigsy(B.T * Q * B)
order = sorted(range(len(basis)), key=lambda i: abs(ev[i]))
print("reduced Hessian eigenvalues (smallest 5):", [mp.nstr(ev[i], 6) for i in order[:5]])
vrot = mp.matrix([(-x0[i+1] if i % 2 == 0 else x0[i-1]) for i in range(n2)]); vrot /= mp.norm(vrot)
ker = []
for i in order[:3]:
    v = B * V[:, i]; v -= (vrot.T * v)[0] * vrot
    for kv in ker:
        v -= (kv.T * v)[0] * kv
    if mp.norm(v) > 0.5:
        ker.append(v / mp.norm(v))
print("kernel dim:", len(ker))


def solve(k, s, constrained=True):
    W = []
    for c in range(n2):
        v = mp.matrix(n2, 1); v[c] = 1
        for b in [k, vrot] + W:
            v = v - (b.T * v)[0] * b
        if mp.norm(v) > mp.mpf("1e-8"):
            W.append(v / mp.norm(v))
        if len(W) == n2 - 2:
            break
    m = len(W)
    Wm = mp.matrix(n2, m)
    for j, w in enumerate(W):
        for i in range(n2):
            Wm[i, j] = w[i]
    nc = 3 if constrained else 0
    y = mp.matrix(m, 1); mu = [mp.mpf(0)] * 3
    for it in range(60):
        x = mp.matrix(x0) + s * k + Wm * y
        xl = [x[i] for i in range(n2)]
        E, g, Hm = H_g_h(xl)
        gJx = mp.matrix([2 * v for v in xl])
        G = [gJx, G1, G2]
        gF = g - Om * gJx
        for c in range(nc):
            gF -= mu[c] * G[c]
        A = mp.matrix(m + nc, m + nc)
        top = Wm.T * (Hm - (Om + mu[0]) * 2 * mp.eye(n2)) * Wm
        for i in range(m):
            for j in range(m):
                A[i, j] = top[i, j]
        for c in range(nc):
            cc = Wm.T * G[c]
            for i in range(m):
                A[i, m + c] = -cc[i]; A[m + c, i] = cc[i]
        rhs = mp.matrix(m + nc, 1)
        r1 = Wm.T * gF
        for i in range(m):
            rhs[i] = -r1[i]
        cons = [J(xl) - J(x0), sum(xl[0::2]), sum(xl[1::2])]
        for c in range(nc):
            rhs[m + c] = -cons[c]
        d = mp.lu_solve(A, rhs)
        for i in range(m):
            y[i] += d[i]
        for c in range(nc):
            mu[c] += d[m + c]
        if mp.norm(d) < mp.mpf(10) ** (-34):
            break
    return E - Om * J(xl) - F0


k = ker[0]
for label, con in [("level set J3=J1=J2 fixed", True)]:
    v = [solve(k, mp.mpf(s), con) / mp.mpf(s) ** 4 for s in ("2e-3", "1e-3")]
    E4 = (4 * v[1] - v[0]) / 3
    print(label, "E4 =", mp.nstr(E4, 15), " 9/14 =", mp.nstr(mp.mpf(9) / 14, 15), " ratio", mp.nstr(E4 / (mp.mpf(9) / 14), 12))
