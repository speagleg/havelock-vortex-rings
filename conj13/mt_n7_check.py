"""Check that engine's 'linear-slice' quantity reproduces MT's own n=7 numbers (beta, gamma, delta):
normalisation-free ratio (delta - gamma^2/beta)/delta = (441/4)/(1071/8) = 14/17 from MT p.23."""
import mpmath as mp
mp.mp.dps = 40
from engine import Ring, fit_even_poly
R = Ring(7, 0, 1)
E, V, normals, t = R.constrained_hessian()
o = sorted(range(len(E)), key=lambda i: abs(E[i]))
e = V[o[0]]
eps = [mp.mpf(k)/2000 for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5)]
vals = []
for ep in eps:
    q = mp.matrix(R.q0) + ep*e
    q = [q[i] for i in range(14)]
    vals.append(R.H(q) - R.H(R.q0) - R.omega*(R.J(q)[2] - R.Jmu[2]))
c = fit_even_poly(eps, vals, list(range(2, 12)))
d_straight = c[2]
print("straight-line (no slaving) quartic delta:", mp.nstr(d_straight, 15))
print("linear-slice minimised (from reduced_quartic): 6/7 ; ratio =", mp.nstr(mp.mpf(6)/7/d_straight, 15),
      " MT ratio 14/17 =", mp.nstr(mp.mpf(14)/17, 15))
print("constrained (true reduced) 9/14 ; ratio to straight =", mp.nstr(mp.mpf(9)/14/d_straight, 15))
