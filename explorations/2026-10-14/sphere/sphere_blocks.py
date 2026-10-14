"""Fourier blocks of a ring of N vortices on the sphere (kernel -ln sin(d/2)), exact in z0 = cos(theta), via sympy.
Mode m: Hessian of H in (z, phi) restricted to e^{i m k 2pi/N} perturbations = [[A_m, i C_m], [-i C_m, B_m]];
the linearised flow has lambda^2 = -(A_m B_m - C_m^2), so mode m is stable iff D_m = A_m B_m - C_m^2 >= 0.
Usage: python sphere_blocks.py"""
import sympy as sp
zi, zj, pi_, pj = sp.symbols('z_i z_j phi_i phi_j', real=True)
z0 = sp.symbols('z0', positive=True)
c = zi*zj + sp.sqrt(1-zi**2)*sp.sqrt(1-zj**2)*sp.cos(pi_-pj)
Hp = -sp.log(1-c)/2
D = sp.symbols('Delta', real=True)
def at(e):
    return sp.simplify(e.subs({zi: z0, zj: z0}).subs({pi_: D, pj: 0}))
d2 = {k: at(sp.diff(Hp, *v)) for k, v in {'zizi': (zi, zi), 'zizj': (zi, zj), 'pipi': (pi_, pi_), 'pipj': (pi_, pj),
                                            'zipi': (zi, pi_), 'zipj': (zi, pj)}.items()}
def blocks(N, m):
    A = B = C = 0
    for j in range(1, N):
        th = 2*sp.pi*j/N; e = sp.exp(sp.I*m*th)
        sub = lambda x: x.subs(D, -th)    # Delta = phi_0 - phi_j
        A += sub(d2['zizi']) + sub(d2['zizj'])*e
        B += sub(d2['pipi']) + sub(d2['pipj'])*e
        C += sub(d2['zipi']) + sub(d2['zipj'])*e
    return [sp.nsimplify(sp.simplify(sp.expand_complex(x)), rational=False) for x in (A, B, C)]
if __name__ == '__main__':
    for k, v in d2.items(): print(k, '=', sp.factor(v))
    for N in [4, 5, 6, 7, 8]:
        out = []
        for m in range(1, N//2 + 1):
            A, B, C = blocks(N, m)
            Dm = sp.factor(sp.simplify(A*B + C**2 if sp.im(C) == 0 else A*B - (C/sp.I)**2))
            out.append((m, Dm))
        print(f'N = {N}')
        for m, Dm in out: print(f'   m = {m}: D_m = {Dm}')
