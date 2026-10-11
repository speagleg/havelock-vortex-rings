"""Gap (i) spot check + sign meaning: inertia of the reduced Hessian of H - omega J3 on
ker dJ(q0) cap t^perp at lam = b_n (float, mpmath), via the bordered matrix
[[Hs, C^T],[C, 0]] with C = (dJ1, dJ2, dJ3, t): inertia = inertia(Hs|kerC) + (4 pos, 4 neg).
Expect (odd n): exactly 2 zero eigenvalues (the mode-l kernel) and all other restricted eigenvalues > 0,
so the ring is a degenerate MINIMUM direction-wise and c4 > 0 means strict local min (stability).
Also an off-threshold control: at lam = b_n + 1e-3 the e-multiplier nu_e1 must be nonzero."""
import sys, mpmath as mp
from run_exact import lam_roots
from kkt_exact import run
mp.mp.dps = 30
for n in [int(a) for a in sys.argv[1:]]:
    l = n // 2; ph = [2 * mp.pi * k / n for k in range(n)]
    lam = [x for x in lam_roots(n, mp.sqrt(n - 3 if n % 2 else 7)) if -1 < x <= 1][0]
    args = ([mp.cos(x) for x in ph], [mp.sin(x) for x in ph], [mp.cos(l * x) for x in ph], [mp.sin(l * x) for x in ph], mp.mpf(1))
    r = run(n, lam, *args)
    N = 2 * n; M = r['M']; keep = list(range(N + 4))
    B = mp.matrix(N + 4, N + 4)
    for i, a in enumerate(keep):
        for j, b in enumerate(keep):
            B[i, j] = M[a][b] * (-1 if b >= N else 1)
    asym = max(abs(B[i, j] - B[j, i]) for i in range(N + 4) for j in range(N + 4))
    ev = mp.eigsy(B, eigvals_only=True)
    tol = mp.mpf(10) ** -20
    pos = [x for x in ev if x > tol]; neg = [x for x in ev if x < -tol]; zer = [x for x in ev if abs(x) <= tol]
    r2 = run(n, lam + mp.mpf('1e-3'), *args)
    print('n=%d b_n=%s asym=%.1e  bordered inertia (+,-,0)=(%d,%d,%d) -> restricted: %d pos, %d neg, %d zero; '
          'min pos bordered eig=%s; c4=%s; off-threshold nu_e1=%s' % (
          n, mp.nstr(lam, 10), float(asym), len(pos), len(neg), len(zer), len(pos) - 4, len(neg) - 4, len(zer),
          mp.nstr(min(pos), 5), mp.nstr(r['c4'], 12), mp.nstr(r2['nue1'], 5)), flush=True)
