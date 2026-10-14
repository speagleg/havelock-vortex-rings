import sys, numpy as np
from active_background import run
from ring_blocks import S_window, S_trap, S_H2
N, r0 = 8, 0.5; S0 = S_trap(N, r0, 0.0); s_star = S_H2(N, r0) - S0
for Nb in map(int, sys.argv[1:]):
    for delta in (0.03, 0.06):
        t, dev, nb, g, tot, nfev = run(s_star, Nb, 4, delta=delta)
        k = np.argmax(dev > 1e-2) if np.any(dev > 1e-2) else -1
        print(f'Nb {nb:5d} delta {delta}: spread at 1,2,3,4 periods = {[f"{np.interp(p, t, dev):.1e}" for p in (1,2,3,4)]}; '
              f'first > 1e-2 at {t[k] if k >= 0 else None}  [nfev {nfev}]', flush=True)
