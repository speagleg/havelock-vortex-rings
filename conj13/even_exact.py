"""Conj13, even n = 2M: exact reduced quartic on the D_{n/2}-fixed (two staggered n/2-gons) family,
restricted to the momentum level set J3 = const, for the Montaldi-Tokieda family (their (1.12) Green's
function, curvature 4*lambda; lambda<0 is the sinh model on H^2, lambda>0 the sphere, lambda=0 the plane).

H = -sum_{i<j} log|z_i-z_j|^2 + (n-1) sum_j log(1+lam|z_j|^2)        (MT f(t) times 4 pi, omega-term dropped)
J3 = sum_j g(|z_j|^2),  g(u) = u/(1+lam u)                              (MT momentum, rotations)
Two staggered M-gons: z_{2k} = rho_+ w^{2k}, z_{2k+1} = rho_- w^{2k+1}, w = e^{2 pi i/n}.
Closed form (product formulas, see the proof of Thm. even in paper/main.tex):
  sum_{pairs} log|z_i-z_j|^2 = M(M-1)/2 (log u_+ + log u_-) + 2M log(u_+^{M/2} + u_-^{M/2}) + const.
Constrained curve: g(u_pm) = g0 +- tau (J3 exactly conserved). h(tau) = H on the curve; h is even.
"""
import sympy as sp

M, lam, g0, tau, b = sp.symbols('M lambda g0 tau b', real=True)
n = 2*M

def series_even(expr, var, order=4):
    s = sp.series(expr, var, 0, order+1).removeO()
    return [sp.simplify(s.coeff(var, k)) for k in range(order+1)]

# --- check the product formula numerically first
import mpmath as mp
def pair_sum_direct(Mv, rp, rm):
    nv = 2*Mv; w = mp.exp(2j*mp.pi/nv)
    z = [(rp if k % 2 == 0 else rm)*w**k for k in range(nv)]
    return sum(mp.log(abs(z[i]-z[j])**2) for i in range(nv) for j in range(i+1, nv))
def pair_sum_closed(Mv, rp, rm):
    return Mv*(Mv-1)/mp.mpf(2)*(mp.log(rp**2)+mp.log(rm**2)) + 2*Mv*mp.log(rp**Mv + rm**Mv)
mp.mp.dps = 30
for Mv in (2, 3, 4, 5):
    # constant: the closed form omits sum log|1-eta^k|^2 terms, which are constants; compare differences
    d1 = pair_sum_direct(Mv, mp.mpf('0.7'), mp.mpf('0.4')) - pair_sum_direct(Mv, mp.mpf('0.3'), mp.mpf('0.5'))
    d2 = pair_sum_closed(Mv, mp.mpf('0.7'), mp.mpf('0.4')) - pair_sum_closed(Mv, mp.mpf('0.3'), mp.mpf('0.5'))
    print(f"product-formula check M={Mv}: diff = {mp.nstr(d1-d2, 5)}")

# --- constrained reduced function, exact in M, lambda, g0
gp, gm = g0 + tau, g0 - tau
up, um = gp/(1 - lam*gp), gm/(1 - lam*gm)          # inverse of g
Lsig = sp.log(up) + sp.log(um)
Ldel = sp.log(up) - sp.log(um)
h = -(M**2 - M/2)*Lsig - 2*M*sp.log(sp.cosh(M*Ldel/4)) \
    - (2*M - 1)*M*(sp.log(1 - lam*gp) + sp.log(1 - lam*gm))
c = series_even(h, tau, 4)
assert sp.simplify(c[1]) == 0 and sp.simplify(c[3]) == 0
# express in b = lam*u0 with u0 = 1  ->  g0 = 1/(1+b), lam = b
sub = {lam: b, g0: 1/(1+b)}
c2 = sp.factor(sp.simplify(c[2].subs(sub)))
c4 = sp.factor(sp.simplify(c[4].subs(sub)))
print("tau^2 coefficient h2(M,b) =", c2)
print("tau^4 coefficient h4(M,b) =", c4)

# threshold: MT (3.7)  (1+b^2)/(1+b)^2 = n^2/(8(n-1))
thr = sp.numer(sp.together((1+b**2)*8*(n-1) - n**2*(1+b)**2))
print("threshold polynomial:", sp.factor(sp.expand(thr)))
# h2 must vanish on the threshold
h2num = sp.numer(sp.together(c2))
print("h2 numerator:", sp.factor(h2num))
print("rem(h2num, thr) =", sp.simplify(sp.rem(sp.expand(h2num), sp.expand(thr), b)))

# --- evaluate at thresholds
def report(Mv):
    nv = 2*Mv
    roots = [r for r in sp.solve(thr.subs(M, Mv), b) if r.is_real and -1 < r <= 1]
    for r in roots:
        r = sp.nsimplify(r)
        h2v = sp.simplify(c2.subs({M: Mv, b: r}))
        h4v = sp.radsimp(sp.simplify(c4.subs({M: Mv, b: r})))
        print(f"n={nv}: b*={r} = {sp.N(r, 20)}; h2={h2v}; h4 = {h4v} = {sp.N(h4v, 25)}")
for Mv in (2, 3, 4, 5, 6, 7, 8, 10, 20, 50):
    report(Mv)
