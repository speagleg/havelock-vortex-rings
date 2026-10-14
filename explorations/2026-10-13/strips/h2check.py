"""Grid-convergence check (H^2, R=2, D=0.2) and targeted scans at R=3 (D=0.2 and 0.1).  Usage: python h2check.py"""
import sys; sys.path.insert(0, '.')
from kuo import *
import numpy as np

if __name__ == '__main__':
    for n in [2666, 5000]:
        gs = {m: growth('h2', 5.0, n, tophat(1.8, 2.0), m) for m in [15, 16, 17, 18, 19]}
        print('conv n', n, {m: round(g, 5) for m, g in gs.items()}, 'm*', max(gs, key=gs.get), flush=True)
    for (r1, r2, ms_) in [(2.8, 3.0, range(36, 60, 2)), (2.9, 3.0, range(76, 116, 4))]:
        gs = {m: growth('h2', 5.0, 3000, tophat(r1, r2), m) for m in ms_}
        ms = max(gs, key=gs.get); r = (r1+r2)/2
        print(f'R=3 D={r2-r1:.1f}: m*={ms} (edge of scanned range? {ms in (min(ms_), max(ms_))}); cosh r={np.cosh(r):.3f} '
              f'vs ring threshold {0.1343122*ms+0.54/ms:.3f}; spacing {2*np.pi*np.sinh(r)/ms:.3f}', flush=True)
