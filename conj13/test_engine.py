import mpmath as mp, sys
mp.mp.dps = 50
from engine import Ring, fit_even_poly
def run(n, lam, r0, model='MT', kdir=0, mix=False):
    R = Ring(n, lam, r0, model)
    E, V, normals, t = R.constrained_hessian()
    order = sorted(range(len(E)), key=lambda i: abs(E[i]))
    print(f"n={n} lam={lam} r0^2={mp.nstr(R.r0**2,12)} model={model}: smallest |eig|:",
          [mp.nstr(E[i], 6) for i in order[:4]], " min eig:", mp.nstr(min(E), 6))
    nk = 1 if n % 2 == 0 else 2
    ker = [V[i] for i in order[:nk]]
    F = [V[i] for i in order[nk:]]
    e = ker[0] if not mix else (ker[0]+ker[1])/mp.sqrt(2)
    if mix: F = F  # second kernel vector excluded from F: keep it free? no -> include
    if nk == 2:
        other = ker[1] if not mix else (ker[0]-ker[1])/mp.sqrt(2)
        F = F + [other]
    eps = [mp.mpf(k)/400 for k in (1, 2, 3, 4)]
    vals = []
    guess = None
    for ep in eps:
        v, guess = R.f_eff(ep, e, F, normals, guess=guess)
        vals.append(v)
    c = fit_even_poly(eps, vals, [2, 4, 6, 8])
    print("   constrained", "c2,c4,c6 =", [mp.nstr(x, 12) for x in c[:3]])
    return c
if __name__ == '__main__':
    run(7, 0, 1)
