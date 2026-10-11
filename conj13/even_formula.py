"""Conj13 even n: symbolic proof that on the threshold the reduced quartic h4 (even_exact.py) equals
  h4 = (n-1)(n-2)(1+b)^4 (1+b^2)^2 (n^4 - 8n^3 + 64n^2 - 144n + 96) / (48 n^3),
and that n^4-8n^3+64n^2-144n+96 = m^4+40m^2+48m+16 (m=n-2) > 0.  Exact (sympy polynomial remainder)."""
import sympy as sp
M, b, n, m = sp.symbols('M b n m')
h4 = M*(b + 1)**4*(M**4*b**4 + 4*M**4*b**3 + 6*M**4*b**2 + 4*M**4*b + M**4 - 16*M**2*b**4 - 16*M**2*b**3
      - 16*M**2*b - 16*M**2 + 48*M*b**4 + 48*M - 24*b**4 - 24)/96          # printed by even_exact.py
thr = M**2*b**2 + 2*M**2*b + M**2 - 4*M*b**2 - 4*M + 2*b**2 + 2              # MT (3.7), n = 2M
claim = (2*M-1)*(2*M-2)*(1+b)**4*(1+b**2)**2*((2*M)**4 - 8*(2*M)**3 + 64*(2*M)**2 - 144*(2*M) + 96)/(48*(2*M)**3)
diff = sp.together(h4 - claim)
num = sp.expand(sp.numer(diff))
r = sp.rem(sp.Poly(num, b), sp.Poly(thr, b))
print("remainder of numer(h4 - claim) modulo threshold polynomial (in b):", sp.simplify(r.as_expr()))
print("n^4-8n^3+64n^2-144n+96 at n=m+2:", sp.expand((m+2)**4 - 8*(m+2)**3 + 64*(m+2)**2 - 144*(m+2) + 96))
# normalised coefficient used by reduced_quartic.py (unit Cartesian kernel vector, r0 = 1):
c4 = sp.factor(sp.simplify((n-1)*(n-2)*(n**4-8*n**3+64*n**2-144*n+96)/(48*n**3) * (2/sp.sqrt(n))**4
                           * (n**2/(8*(n-1)))**2))
print("c4(n) [unit kernel vector, r0=1] =", c4, " values:", [c4.subs(n, k) for k in (4, 6, 8, 10, 12)])
