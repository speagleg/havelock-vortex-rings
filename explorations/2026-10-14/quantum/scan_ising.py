"""Locate the quantum critical spacing of the pinned row at kappa = 0.98 (continuous side) for two quantum strengths:
zigzag order parameter eta and correlation length xi vs r = l/l* at bond dimension chi = 16.
Usage: python scan_ising.py -> scan_ising.out"""
import sys; from dmrg_chain import run_point
out = open('scan_ising.out', 'w')
for tau in [5e-4, 2e-3]:
    for r in [0.90, 0.93, 0.95, 0.97, 0.99, 1.01]:
        res, psi, m = run_point(0.98, r, tau, chis=(16,))
        x = res[-1]; line = f'tau {tau:.0e}  r {r:.2f}: eta {x["eta"]: .4f}  <y> {x["ymean"]: .4f}  xi {x["xi"]:8.2f}  S {x["S"]:.4f}'
        print(line, flush=True); out.write(line + '\n'); out.flush()
