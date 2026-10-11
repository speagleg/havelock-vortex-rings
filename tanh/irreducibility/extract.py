"""Extract the threshold polynomials f_N (N=8..16) from two independent sources and check they agree.

Source 1: lean/Havelock/TanhAlgebra.lean, `def fN : List Z := [...]` (constant term first).
Source 2: tanh/exact.py -- factor A_{floor(N/2)} over Q and take the unique irreducible factor
          with a root in (0,1), normalised to positive leading coefficient and content 1.
Writes polys.gp (PARI input; Pol() takes coefficients highest degree first) and prints a comparison.
Run from the repository root:  python3 tanh/irreducibility/extract.py > tanh/irreducibility/extract.out
"""
import os, re, sys
import sympy as sp
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'tanh'))
from exact import polys, x  # noqa: E402

lean = open(os.path.join(ROOT, 'lean', 'Havelock', 'TanhAlgebra.lean')).read()
ok_all = True
gp_lines = []
for N in range(8, 17):
    m = re.search(r'^def f%d : List ℤ := \[([^\]]*)\]' % N, lean, re.M)
    lean_coeffs = [int(t) for t in m.group(1).split(',')]          # constant term first
    A, _, _ = polys(N, N // 2)
    _, facs = A.factor_list()
    cands = []
    for f, e in facs:
        roots = [r for r in sp.Poly(f, x).real_roots() if 0 < r < 1]   # exact isolation
        if roots:
            cands.append(f)
    assert len(cands) == 1, (N, cands)
    f = sp.Poly(cands[0], x).primitive()[1]
    if f.LC() < 0:
        f = -f
    exact_coeffs = [int(c) for c in reversed(f.all_coeffs())]       # constant term first
    same = exact_coeffs == lean_coeffs
    ok_all &= same
    print('N=%2d deg=%2d lead=%3d  lean==exact.py: %s' % (N, len(lean_coeffs) - 1, lean_coeffs[-1], same))
    gp_lines.append('f[%d] = Pol(%s);' % (N, list(reversed(lean_coeffs))))
print('ALL AGREE' if ok_all else 'MISMATCH')
with open(os.path.join(os.path.dirname(__file__), 'polys.gp'), 'w') as fh:
    fh.write('\\\\ Threshold polynomials f_N (from TanhAlgebra.lean; checked equal to exact.py by extract.py)\n')
    fh.write('f = vector(16);\n' + '\n'.join(gp_lines) + '\n')
