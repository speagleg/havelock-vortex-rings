"""Beyond quadratic order and at finite temperature: the flat band of a hypercycle row, in the nearest-neighbour model
(the exact leading order in e^{-l'} of the tanh model).
Model: H = sum_j w_j, w_j = 2 sin(phi_j) sin(phi_{j+1}) e^{-D_j} (D_j = u_{j+1} - u_j), F = H - V I with I = -sum cot(phi_j)
(so F = H + V sum cot phi_j), period L = N l' (sum D_j = L).
(A) Spacings relaxed at fixed L: minimising over D_j makes all bond weights equal (Lagrange), so
        F_rel(phi) = 2 N e^{-l'} (prod_j sin^2 phi_j)^{1/N} + V sum_j cot phi_j.
    The first term depends only on the mean of x_j = ln sin phi_j, so every non-uniform direction sees independent
    single-site potentials: the relaxed transverse stiffness is the same for all non-uniform modes (an exactly flat band,
    to all orders), k_perp = 2 P csc^2(alpha) cos(2 alpha), P = 2 sin^2(alpha) e^{-l'} the bond energy.
(B) Pressure ensemble at temperature T (Liouville measure du dY = du dphi / sin^2 phi): integrating each D_j over R
    against exp(-beta (w_j + P D_j)) gives Gamma(beta P) (beta c_j)^{-beta P}, c_j = 2 sin phi_j sin phi_{j+1}, so the
    phi_j are EXACTLY independent with weight exp(-beta h_T(phi)),  h_T(phi) = 2 (P + T) ln sin(phi) + V cot(phi).
    At T = 0: h' = 0 <=> sin(2 phi) = V/P. The row (phi = alpha) is the minimum and the barrier sits at phi = pi/2 - alpha,
    the conjugate hypercycle (curvature sqrt(1 - kappa^2)); they merge at kappa_c = 1/sqrt 2 (saddle-node).
    Barrier per vortex: Delta E = 2 P [ln cot(alpha) - cos(2 alpha)] = (16/3) P (pi/4 - alpha)^3 + O((pi/4 - alpha)^5).
(C) In the full tanh model the relaxed transverse stiffness is Sigma(q) = s^2 (B |C| - delta^2)/|C| (Schur complement
    of the exact symbol); E(q) >= 0 <=> Sigma(q) >= 0, Sigma(pi) = s^2 B(pi). Its q-spread (the lifting of the flat band by
    next-nearest neighbours) is checked to be O(e^{-2l'}) in absolute terms, against a mean ~ 4 s^4 (2K-1) e^{-l'}.
Usage: python flat_nonlinear.py -> flat_nonlinear.out"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows'))
import numpy as np, sympy as sp
from scipy.optimize import minimize
from symbol import symbol
from exact_symbol import lstar
out = []; P_ = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))

# (A) relaxation of the spacings: numerical minimisation vs the closed form
rng = np.random.default_rng(4); N = 7; lp = 2.0; L = N*lp
for trial in range(3):
    phi = np.pi/4 + 0.3*rng.standard_normal(N)*0.5
    c = 2*np.sin(phi)*np.sin(np.roll(phi, -1))
    Dfree = rng.standard_normal(N - 1)
    def Hb(Df):
        D = np.append(Df, L - Df.sum()); return np.sum(c*np.exp(-D))
    res = minimize(Hb, Dfree, method='BFGS', options={'gtol': 1e-13})
    closed = 2*N*np.exp(-lp)*np.prod(np.sin(phi)**2)**(1/N)
    P_(f'(A) trial {trial}: min over spacings {res.fun:.12f}  closed form {closed:.12f}  diff {abs(res.fun - closed):.1e}')
# Hessian of F_rel at the uniform row: all non-uniform eigenvalues equal k_perp
for kappa in [0.6, 1/np.sqrt(2), 0.8, 0.95]:
    a = np.arccos(kappa); Pb = 2*np.sin(a)**2*np.exp(-lp); V = Pb*np.sin(2*a)
    Frel = lambda ph: 2*N*np.exp(-lp)*np.prod(np.sin(ph)**2)**(1/N) + V*np.sum(1/np.tan(ph))
    h = 1e-4; x0 = np.full(N, a); Hm = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            e1, e2 = np.eye(N)[i]*h, np.eye(N)[j]*h
            Hm[i, j] = (Frel(x0 + e1 + e2) - Frel(x0 + e1 - e2) - Frel(x0 - e1 + e2) + Frel(x0 - e1 - e2))/(4*h*h)
    ev = np.sort(np.linalg.eigvalsh(Hm)); kperp = 2*Pb*np.cos(2*a)/np.sin(a)**2
    grad = np.array([(Frel(x0 + np.eye(N)[i]*h) - Frel(x0 - np.eye(N)[i]*h))/(2*h) for i in range(N)])
    P_(f'(A) kappa {kappa:.4f}: |grad F_rel| at the row = {np.abs(grad).max():.1e}; Hessian eigenvalues {np.round(ev, 7)};'
       f' predicted non-uniform k_perp = {kperp:.7f}')

# (B) the single-site potential, its extrema and the barrier
phi, al, Pp, T = sp.symbols('phi alpha P T', positive=True)
V = Pp*sp.sin(2*al)
hT = 2*(Pp + T)*sp.log(sp.sin(phi)) + V*sp.cot(phi)
h0 = hT.subs(T, 0)
P_('(B) h0\'(alpha) =', sp.simplify(sp.diff(h0, phi).subs(phi, al)), '; h0\'(pi/2 - alpha) =', sp.simplify(sp.diff(h0, phi).subs(phi, sp.pi/2 - al)))
P_('    h0\'\'(alpha) =', sp.simplify(sp.diff(h0, phi, 2).subs(phi, al)), '(= 2 P csc^2 cos 2 alpha)')
dE = sp.simplify(h0.subs(phi, sp.pi/2 - al) - h0.subs(phi, al))
P_('    barrier h0(pi/2 - alpha) - h0(alpha) =', dE, '-> 2P[ln cot(alpha) - cos(2 alpha)] check:',
   sp.simplify(sp.expand_log(dE - 2*Pp*(sp.log(sp.cot(al)) - sp.cos(2*al)), force=True).rewrite(sp.sin)) == 0 or
   abs(float((dE - 2*Pp*(sp.log(sp.cot(al)) - sp.cos(2*al))).subs({Pp: 1.3, al: 0.5}))) < 1e-12)
e = sp.symbols('e', positive=True)
P_('    near kappa_c (alpha = pi/4 - e): barrier/(2P) =', sp.series((sp.log(sp.cot(al)) - sp.cos(2*al)).subs(al, sp.pi/4 - e), e, 0, 7))
# Liouville measure + pressure integral: weight exp(-beta h_T); check the Gamma-function step numerically
beta, Pv, cv = 3.0, 0.7, 1.9
from scipy.integrate import quad
from scipy.special import gamma
lhs = quad(lambda D: np.exp(-beta*(cv*np.exp(-D) + Pv*D)), -30, 60, limit=500)[0]
P_(f'(B) int_R exp(-beta(c e^-D + P D)) dD = {lhs:.10f}  vs Gamma(beta P)(beta c)^(-beta P) = {gamma(beta*Pv)*(beta*cv)**(-beta*Pv):.10f}')

# (C) flatness of the relaxed transverse stiffness in the full tanh model
P_('(C) full tanh model: relaxed transverse stiffness Sigma(q) = s^2 (B|C| - delta^2)/|C| over q in (0, pi]:')
P_('    kappa   l/l*    l\'     mean Sigma     (max-min) over q     e^{-2l\'}     (max-min)/(s^2 e^{-2l\'})   mean/(s^4 e^{-l\'} (2K-1))')
for kappa in [0.72, 0.75, 0.8, 0.9]:
    s = np.sqrt(1 - kappa**2); ls = lstar(kappa)
    for r in [1.5, 2.0, 3.0, 4.0]:
        lpp = ls*r*s; q = np.linspace(0.05, np.pi, 400)
        B, C, dl = symbol(kappa, lpp, q); Sig = s*s*(B*np.abs(C) - dl**2)/np.abs(C)
        sp_ = Sig.max() - Sig.min()
        P_(f'    {kappa:5.2f}   {r:4.1f}   {lpp:6.3f}   {Sig.mean():11.4e}   {sp_:11.4e}        {np.exp(-2*lpp):10.3e}    {sp_/(s*s*np.exp(-2*lpp)):9.4f}'
           f'              {Sig.mean()/(s**4*np.exp(-lpp)*(2*kappa**2 - 1)):8.4f}')
open('flat_nonlinear.out', 'w').write('\n'.join(out) + '\n')
