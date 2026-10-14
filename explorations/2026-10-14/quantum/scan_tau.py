"""Quantum melting of the zigzag order of the pinned row: kappa = 0.98, l = 0.95 l* (classically ordered),
scan the quantum strength tau = k eps^2 at chi = 16. Run in three passes (tau lists edited between runs):
[0.004, 0.01] -> scan_tau_coarse.out, [0.005..0.009] -> scan_tau_fine.out, [0.0073, 0.0076, 0.0079] -> scan_tau_window.out. Usage: python scan_tau.py -> scan_tau_window.out"""
from dmrg_chain import run_point
out = open('scan_tau_window.out', 'w')
for tau in [0.0076, 0.0073, 0.0079]:
    res, psi, m = run_point(0.98, 0.95, tau, chis=(16,), ymax=1.6)
    x = res[-1]; line = f'tau {tau:.4f}  r 0.95: eta {x["eta"]: .4f}  <y> {x["ymean"]: .4f}  xi {x["xi"]:8.2f}  S {x["S"]:.4f}'
    print(line, flush=True); out.write(line + '\n'); out.flush()
