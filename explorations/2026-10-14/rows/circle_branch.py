"""Circle branch (kappa > 1) of the stability boundary from the havelock package's CERTIFIED tanh thresholds:
a ring of N at its threshold radius r* has geodesic curvature coth r* and spacing 2 pi sinh r* / N.
Usage (repo root, package installed): python explorations/2026-10-14/rows/circle_branch.py > circle_branch.out"""
import math
import havelock as hv

print('N kappa=coth r*  spacing=2 pi sinh r*/N')
for N in [8, 9, 10, 11, 12, 14, 16, 20, 24, 32, 48, 64]:
    r = float(hv.threshold('tanh', N).r)
    print(N, 1/math.tanh(r), 2*math.pi*math.sinh(r)/N)
