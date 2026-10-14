"""Symbolic proofs of the inputs of q2_instability_cert.py that the certificate itself checks only at sample points.
(1) -Q2'(C) <= 1/(C^2 (C^2-1)) (the certificate uses the weaker 3/(C^2(C^2-1))). With t = 1/C,
      3 t^5 h = 3t - 2t^3 - 3(1-t^2) atanh t,   h := C^2 (C^2-1)(-Q2')/3,
    and (1-t^2) atanh t = t - sum_{k>=1} 2 t^(2k+1)/(4k^2-1), so h = sum_{k>=2} 2 t^(2k-4)/(4k^2-1): positive coefficients,
    increasing in t, sup at t -> 1 equal to sum_{k>=2} 2/(4k^2-1) = 1/3.
(2) psi_n = (d - 2 cot)(d - cot) sin(n a) (Darboux ladder), and it solves -psi'' + 6 csc^2(a) psi = n^2 psi (symbolic n).
(3) Norms: with A_nu = d - nu cot, A_nu^+ A_nu g = (n^2 - nu^2) g on the ladder and (A g)^2 = g A^+ A g + (g A g)',
    the boundary term vanishing at 0 and pi, so ||psi_n||^2 = (n^2-1)(n^2-4) ||sin n a||^2 = pi (n^2-1)(n^2-4)/2.
(4) c3 = w_3 (A_3 - B_3) = (16/5) s^6 cos 2a exactly.
Usage: python q2_lemmas.py"""
import sympy as sp
C, t, a = sp.symbols('C t a', positive=True); n = sp.symbols('n', integer=True, positive=True)
Q0 = sp.atanh(1/C); Q2 = (3*C**2 - 1)/2*Q0 - sp.Rational(3, 2)*C
h = sp.simplify((C**2*(C**2 - 1)*(-sp.diff(Q2, C))/3).subs(C, 1/t))
print('(1) 3 t^5 h - (3t - 2t^3 - 3(1-t^2) atanh t) =', sp.simplify(3*t**5*h - (3*t - 2*t**3 - 3*(1 - t**2)*sp.atanh(t))))
kk = sp.Symbol('k', integer=True, positive=True)
ser = sp.series((1 - t**2)*sp.atanh(t), t, 0, 16).removeO()
print('    (1-t^2) atanh t = t - sum 2 t^(2k+1)/(4k^2-1) through t^15:',
      sp.expand(ser - (t - sum(sp.Rational(2, 4*k*k - 1)*t**(2*k + 1) for k in range(1, 8)))) == 0)
print('    coefficient identity (general k): (1/(2k+1) - 1/(2k-1)) =', sp.simplify(1/(2*kk + 1) - 1/(2*kk - 1) + 2/(4*kk**2 - 1)), '+ (-2/(4k^2-1))')
print('    sum_{k>=2} 2/(4k^2-1) =', sp.summation(2/(4*kk**2 - 1), (kk, 2, sp.oo)))
q = sp.cot(a)
psi = (1 - n**2 + 3*q**2)*sp.sin(n*a) - 3*n*q*sp.cos(n*a)
f = sp.sin(n*a); f1 = sp.diff(f, a) - q*f; f2 = sp.diff(f1, a) - 2*q*f1
print('(2) ladder image minus psi_n:', sp.simplify(sp.expand_trig(f2 - psi)),
      '; mode-equation residual:', sp.simplify(sp.expand_trig(-sp.diff(psi, a, 2) + 6/sp.sin(a)**2*psi - n**2*psi)))
for nu, g in [(1, f), (2, f1)]:
    Ag = sp.diff(g, a) - nu*q*g; AdA = -sp.diff(Ag, a) - nu*q*Ag
    print(f'(3) nu = {nu}: A^+A g - (n^2-nu^2) g =', sp.simplify(sp.expand_trig(AdA - (n**2 - nu**2)*g)),
          '; (Ag)^2 - g A^+A g - (g Ag)\' =', sp.simplify(sp.expand_trig(Ag**2 - g*AdA - sp.diff(g*Ag, a))),
          '; g Ag at a -> 0+:', sp.limit((g*Ag).subs(n, 7), a, 0, '+'))
s, k = sp.sin(a), sp.cos(a); p3 = psi.subs(n, 3); dp3 = sp.diff(p3, a)
c3 = sp.Rational(2, 3*8*5)*(2*s*k*p3*dp3 + (6 - 9*s*s)*p3**2 - s*s*dp3**2)
res = sp.together(sp.expand_trig(c3 - sp.Rational(16, 5)*s**6*sp.cos(2*a)).subs(sp.cot(a), sp.cos(a)/sp.sin(a)))
S_, K_ = sp.symbols('S K'); num = sp.expand(sp.numer(res).subs({sp.sin(a): S_, sp.cos(a): K_}))
print('(4) c3 - (16/5) s^6 cos 2a: numerator reduced modulo S^2 + K^2 - 1 =', sp.reduced(num, [S_**2 + K_**2 - 1], K_, S_)[1])
