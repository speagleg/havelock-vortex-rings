"""Asymptotic stability boundary of hypercycle rows (tanh model) vs the full numerical threshold.
A_UP(pi)/(4 sin^2 alpha) = (2k^2-1) e + (8k^4-3) e^2 + 4(8k^6-8k^4+2k^2-1) e^3 + O(e^4),  e = exp(-l sin alpha), k = kappa
(next_order.py). The row is pairing-stable iff this is >= 0. The series is asymptotic: useful only near kappa_c.
Usage: python boundary_series.py"""
import os, numpy as np
d = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rows_boundary.txt'))
print('kappa   l* full    1-term    2-term    3-term')
for k, lf in d:
    s = np.sqrt(1-k*k); c = [2*k*k-1, 8*k**4-3, 4*(8*k**6-8*k**4+2*k*k-1)]
    out = []
    for n in (1, 2, 3):
        if n == 1:
            e = c[0]
        else:
            r = np.roots(c[:n][::-1])
            r = r[np.isreal(r)].real; r = r[r > 0]
            e = r.min() if len(r) else np.nan
        out.append(-np.log(e)/s if e > 0 else np.nan)
    print(f'{k:6.3f}  {lf:8.4f}  ' + '  '.join(f'{x:8.4f}' for x in out))
