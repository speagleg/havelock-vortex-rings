"""Conj13: reduced (effective) quartic coefficient of the Hamiltonian along the critical mode,
computed numerically (NUMERICAL, mpmath 40 digits) by constrained minimisation (engine.py).
For each case: spectrum check of the constrained Hessian (kernel dimension, min positive eigenvalue),
then c2, c4 (cubic vanishes by symmetry; checked via odd fit) of f_eff(eps) both
  - constrained : on the true momentum level set J = mu (correct reduced Hamiltonian),
  - linear      : on MT's linear symplectic slice N (MT Prop. 9 / Lemma 8 setting).
Even n: compares with the exact closed form h4 of even_exact.py.
usage: python reduced_quartic.py CASE...   (CASE like  MT:-1:8  sphere/plane/hyperbolic via lam sign, or tanh:8)
"""
import sys, mpmath as mp
mp.mp.dps = 40
from engine import Ring, fit_even_poly

def b_threshold(n):
    # MT (3.6): (1+b^2)/(1+b)^2 = floor(n^2/4)/(2(n-1)), root in (-1,1]
    k = mp.mpf(n*n//4)/(2*(n-1))
    # (1-k) b^2 - 2k b + (1-k) = 0
    A, B, C = 1-k, -2*k, 1-k
    if abs(A) < mp.mpf(10)**-30:
        return mp.mpf(0)
    roots = [(-B + s*mp.sqrt(B*B-4*A*C))/(2*A) for s in (1, -1)]
    return [r for r in roots if -1 < r <= 1][0]

TANH8 = mp.mpf('0.06735622300408126226843747656623557439402511756974')

def run(model, n, lam_sign=None):
    if model == 'MT':
        b = b_threshold(n)
        lam, r0 = (b, 1) if b != 0 else (0, 1)
    else:
        lam, r0 = -1, mp.sqrt({8: TANH8}[n])
        b = -r0**2
    R = Ring(n, lam, r0, model)
    E, V, normals, t = R.constrained_hessian()
    order = sorted(range(len(E)), key=lambda i: abs(E[i]))
    nk = 1 if n % 2 == 0 else 2
    pos = [E[i] for i in order[nk:]]
    print(f"== model={model} n={n} lam*r0^2={mp.nstr(b, 20)}  omega={mp.nstr(R.omega, 15)}")
    print(f"   kernel eigs: {[mp.nstr(E[i], 3) for i in order[:nk]]};  min of the other {len(pos)} eigs: {mp.nstr(min(pos), 10)}")
    ker = [V[i] for i in order[:nk]]
    F = [V[i] for i in order[nk:]]
    results = {}
    dirs = [('e1', ker[0], ker[1:])]
    if nk == 2:
        dirs.append(('(e1+e2)/sqrt2', (ker[0]+ker[1])/mp.sqrt(2), [(ker[0]-ker[1])/mp.sqrt(2)]))
    for linear in (False, True):
        for name, e, rest in dirs:
            eps = [mp.mpf(k)/2000 for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5)]
            vals, guess = [], None
            for ep in eps:
                v, _ = R.f_eff(ep, e, F[:]+rest, normals, linear=linear)
                vals.append(v)
            c = fit_even_poly(eps, vals, list(range(2, 12)))
            tag = 'linear-slice (MT)' if linear else 'constrained J=mu '
            print(f"   {tag} dir {name}: c2={mp.nstr(c[0], 5)}  c3={mp.nstr(c[1], 5)}  c4={mp.nstr(c[2], 15)}")
            results[(linear, name)] = c[2]
    if model == 'MT' and n % 2 == 0 and r0 == 1:
        # exact h4 (even_exact.py) converted: tau = g(u+)-g0, u+ = (1+eps/sqrt n)^2, g'(1)=1/(1+b)^2
        nn = mp.mpf(n); m = nn-2
        h4 = (nn-1)*m*(1+b)**4*(1+b*b)**2*(m**4+40*m**2+48*m+16)/(48*nn**3)
        # f_eff(eps) = h2 tau^2 + h4 tau^4, tau = a1 eps + a2 eps^2 + ...; with h2 = 0: c4 = h4 a1^4
        a1 = 2/(mp.sqrt(nn)*(1+b)**2)
        print(f"   exact closed form h4*a1^4 = {mp.nstr(h4*a1**4, 15)}")
    return results

if __name__ == '__main__':
    for arg in sys.argv[1:]:
        model, n = arg.split(':')
        run(model, int(n))
        sys.stdout.flush()
