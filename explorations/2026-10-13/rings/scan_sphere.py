"""Shell-BEC configuration scan: ring + mirrored antiring (+ pole charges +-q) on the sphere.
Usage: python scan_sphere.py STAGGER q1 q2 ...   (STAGGER = 0 aligned antiring, 1 staggered by pi/N)
Prints the colatitudes (deg) of the +ring at which the configuration is linearly stable, for N = 2..10."""
import sys; sys.path.insert(0, '.')
from sphere import *
import numpy as np
ths = np.radians(np.arange(3, 88.1, 1.5))

def intervals(st):
    out = []; i = 0
    while i < len(ths):
        if st[i]:
            j = i
            while j+1 < len(ths) and st[j+1]: j += 1
            out.append(f'{np.degrees(ths[i]):.0f}-{np.degrees(ths[j]):.0f}'); i = j+1
        else: i += 1
    return ','.join(out) if out else 'none'

if __name__ == '__main__':
    stag = int(sys.argv[1]); qs = [float(x) for x in sys.argv[2:]]
    for q in qs:
        print(f'== q={q} {"staggered" if stag else "aligned"} antiring; stable colatitudes (deg) of the +ring')
        for N in range(2, 11):
            st = [analyse(*config(N, t, q, stag))[2] < 1e-4 for t in ths]
            print('  N', N, intervals(st), flush=True)
