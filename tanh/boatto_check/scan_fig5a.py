"""Reproduce Boatto Fig. 5a: min nontrivial eigenvalue of S~ vs K, r_o = 1, N = 5..9, both pair laws.
Also locate zero crossings (bisection on det-sign of the per-mode problem is avoided: we bisect the
min eigenvalue directly) and print values at a few K for comparison with the figure."""
import sys
import mpmath as mp
from boatto_hessian import min_nontrivial

mp.mp.dps = 30
models = sys.argv[1:] or ['boatto', 'sinh']


def f(N, K, model):
    return min_nontrivial(N, K, model)[0]


for model in models:
    print('==== model', model)
    for N in range(5, 10):
        Ks = [mp.mpf(k) / 100 for k in range(-60, 61, 5) if k != 0] + [mp.mpf('1e-8'), -mp.mpf('1e-8')]
        Ks = sorted(Ks)
        vals = [(K, f(N, K, model)) for K in Ks]
        line = ' '.join('%s:%s' % (mp.nstr(K, 2), mp.nstr(v, 3)) for K, v in vals)
        print('N=%d  %s' % (N, line))
        for (K1, v1), (K2, v2) in zip(vals, vals[1:]):
            if v1 * v2 < 0:
                lo, hi = K1, K2
                if abs(K1) < 1e-6 or abs(K2) < 1e-6:
                    print("   zero crossing N=%d at K = 0 (planar Thomson heptagon)" % N); continue
                for _ in range(60):
                    mid = (lo + hi) / 2
                    if f(N, mid, model) * v1 > 0:
                        lo = mid
                    else:
                        hi = mid
                print('   zero crossing N=%d: K = %s' % (N, mp.nstr((lo + hi) / 2, 15)))
        sys.stdout.flush()
