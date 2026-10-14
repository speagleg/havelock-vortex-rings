"""Finite-entanglement scaling at a quantum critical point of the pinned row (Pollmann et al. 2009): for increasing bond
dimension chi, S(chi) = (c/6) ln xi(chi) + const. Usage: python central_charge.py kappa r tau -> central_charge_<kappa>_<r>.out"""
import sys, numpy as np
from dmrg_chain import run_point
k, r, tau = map(float, sys.argv[1:4]); chis = tuple(int(x) for x in (sys.argv[4].split(',') if len(sys.argv) > 4 else ['8', '12', '16', '20', '24']))
res, psi, m = run_point(k, r, tau, chis=chis, ymax=1.6, verbose=True)
xs = np.log([x['xi'] for x in res]); Ss = np.array([x['S'] for x in res])
c = 6*np.polyfit(xs, Ss, 1)[0]
lines = [f'kappa {k}, r {r}, tau {tau}'] + [f'  chi {x["chi"]:3d}: xi {x["xi"]:9.3f}  S {x["S"]:.5f}  eta {x["eta"]: .2e}' for x in res]
lines.append(f'  fit S = (c/6) ln xi + const: c = {c:.4f}' + (f';  last-pair slope c = {6*(Ss[-1]-Ss[-2])/(xs[-1]-xs[-2]):.4f}' if len(res) > 1 else ''))
open(f'central_charge_{k}_{r}_{tau}.out', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
