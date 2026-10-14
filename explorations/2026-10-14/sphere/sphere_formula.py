"""Closed form for the Fourier blocks of a latitude ring of N vortices on the unit sphere (kernel -ln sin(d/2)),
from the pair second derivatives of sphere_blocks.py and the classical sums
   sum_{j=1}^{N-1} 1/(1 - cos(2 pi j/N)) = (N^2-1)/6,   sum_j (1 - cos(2 pi m j/N))/(1 - cos(2 pi j/N)) = m(N-m):
   A_m = [(N-1)(1 + z^2) - m(N-m)]/(2 S^4),  B_m = m(N-m)/2,  C_m = 0   (z = cos theta, S = sin theta),
   D_m = A_m B_m = m(N-m)[(N-1)(1+z^2) - m(N-m)]/(4 S^4).
Stable iff (N-1)(1 + cos^2 theta) >= floor(N^2/4), i.e. cos^2 theta >= 1/3, 1/2, 4/5 for N = 4, 5, 6, never for N >= 8,
only at the pole for N = 7 (Polvani & Dritschel 1993; Laurent-Polz, Montaldi & Roberts 2004, Thm 5.2).
In the geodesic curvature kappa = cot(theta) of the circle: kappa^2 (2(N-1) - M) >= M - (N-1), M = floor(N^2/4):
   N = 4: 2 kappa^2 - 1 >= 0 (kappa >= 1/sqrt 2: the ring is a face of the inscribed cube), N = 5: kappa >= 1, N = 6: kappa >= 2.
Usage: python sphere_formula.py"""
import sympy as sp
from sphere_blocks import blocks
z = sp.symbols('z0', positive=True)
ok = True
for N in range(3, 10):
    for m in range(1, N):
        A, B, C = blocks(N, m)
        Af = ((N-1)*(1+z**2) - m*(N-m))/(2*(1-z**2)**2); Bf = sp.Rational(m*(N-m), 2)
        for zv in [sp.Rational(1, 7), sp.Rational(1, 2), sp.Rational(9, 10)]:      # numerical check at 40 digits
            ok &= all(abs(sp.N((X - Y).subs(z, zv), 40)) < 1e-30 for X, Y in [(A, Af), (B, Bf), (C, 0)])
print('closed forms A_m, B_m, C_m = 0 agree with the exact block sums (N = 3..9, all m, three latitudes, 40 digits):', ok)
k = sp.symbols('kappa', positive=True)
for N in range(2, 9):
    M = N*N//4
    p = sp.Rational(M, N-1) - 1
    if p <= 0: print(f'N = {N}: stable at every latitude'); continue
    if p >= 1: print(f'N = {N}: ' + ('stable only at the pole' if p == 1 else 'never stable')); continue
    print(f'N = {N}: stable iff cos^2(theta) >= {p},  i.e. kappa = cot(theta) >= {sp.sqrt(p/(1-p))}')
