"""Pairing mode (m = 2) of a square ring (N = 4) on the sphere, for a GENERIC kernel w(c), c = cos d.
Neighbours sit at c1 = z^2 (Delta = pi/2, two of them) and c2 = 2z^2 - 1 (Delta = pi, one).
Block: A_2 = sum_j [H_zizi + H_zizj e^{i m theta_j}], B_2 likewise in phi, cross term C_2.
Shows which kernel data the threshold depends on (for -ln sin(d/2) it reduces to 3 z^2 - 1).
Usage: python n4_generic.py"""
import sympy as sp
zi, zj, pi_, pj, z = sp.symbols('z_i z_j phi_i phi_j z', real=True)
w1a, w2a, w1b, w2b = sp.symbols('wp1 wpp1 wp2 wpp2')   # w'(c1), w''(c1), w'(c2), w''(c2)
c = zi*zj + sp.sqrt(1-zi**2)*sp.sqrt(1-zj**2)*sp.cos(pi_-pj)
def d2(a, b, Delta, W1, W2):
    e = W2*sp.diff(c, a)*sp.diff(c, b) + W1*sp.diff(c, a, b)
    return sp.simplify(e.subs({zi: z, zj: z}).subs({pi_: Delta, pj: 0}))
A = B = C = 0
for j, (Dl, W1, W2) in enumerate([(-sp.pi/2, w1a, w2a), (-sp.pi, w1b, w2b), (sp.pi/2, w1a, w2a)], start=1):
    e = sp.cos(2*sp.pi*2*j/4)        # e^{i m theta_j} with m = 2: (-1)^j
    A += d2(zi, zi, Dl, W1, W2) + d2(zi, zj, Dl, W1, W2)*e
    B += d2(pi_, pi_, Dl, W1, W2) + d2(pi_, pj, Dl, W1, W2)*e
    C += d2(zi, pi_, Dl, W1, W2) + d2(zi, pj, Dl, W1, W2)*e
A, B, C = [sp.factor(sp.simplify(x)) for x in (A, B, C)]
print('A_2 =', A); print('B_2 =', B); print('C_2 =', C)
# Laplace kernel w = -(1/2) ln(1 - c): w' = 1/(2(1-c)), w'' = 1/(2(1-c)^2)
lap = {w1a: 1/(2*(1-z**2)), w2a: 1/(2*(1-z**2)**2), w1b: 1/(2*(2-2*z**2)), w2b: 1/(2*(2-2*z**2)**2)}
print('Laplace kernel: A_2 B_2 - C_2^2 =', sp.factor(sp.simplify((A*B - C**2).subs(lap))))
