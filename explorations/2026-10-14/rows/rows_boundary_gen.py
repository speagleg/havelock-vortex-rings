"""Hypercycle branch of the stability boundary: critical spacing l*(kappa). ~8 min.
Usage: python rows_boundary_gen.py -> rows_boundary.txt"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rows import stab_hyper, critical

ks = [0.72, 0.73, 0.74, 0.75, 0.77, 0.8, 0.85, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999]
with open('rows_boundary.txt', 'w') as f:
    for k in ks:
        l = critical(lambda l: stab_hyper(40, l, k), 0.3, 80, tol=1e-3)
        f.write(f'{k} {l}\n'); print(k, l, flush=True)
