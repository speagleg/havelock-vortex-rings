"""Sign of each pair term of A_UP(pi) as a polynomial in c = cosh(j l') > 1 and K = kappa^2 in (0, 1).
odd j:  f_odd(c) = (a - b),  even j: f_even(c) = (a + b); common positive factor sin^2 a /((c-1)^2 (c - cos 2a)^2).
Result: odd numerator (c-1)(c+1-2K)(2Kc-c-1) < 0 for K <= 1/2; even numerator (c-1)^2((4K-1)c - (4K^2-2K+1)),
negative for K <= 1/4 but positive at large c for 1/4 < K. instability_proof_check.py handles the sum.
Usage: python pair_sign.py"""
import sympy as sp
c, K, t = sp.symbols('c K t', positive=True)   # K = kappa^2, c = 1 + t
s2 = 1 - K
k2 = K
# numerators from pair_terms.py (sign included), with sin^2 = 1-K, cos^2 = K
odd = -(2*c**3*s2 - c**3 + 4*c**2*s2**2 - 6*c**2*s2 + 3*c**2 - c*s2**3 - c*s2**2 + 3*c*s2 - c*k2**3 - 2*c - 2*s2 + 1)
even = -(4*c**3*s2 - 3*c**3 + 4*c**2*s2**2 - 14*c**2*s2 + 9*c**2 + c*s2**3 - 11*c*s2**2 + 19*c*s2 + c*k2**3 - 10*c + 4*s2**2 - 6*s2 + 3)
for name, P in [('odd', odd), ('even', even)]:
    P = sp.expand(P)
    print(name, 'numerator =', sp.factor(P))
    Pt = sp.Poly(sp.expand(P.subs(c, 1+t)), t)
    print('   in t = c - 1:', [sp.factor(co) for co in Pt.all_coeffs()])
