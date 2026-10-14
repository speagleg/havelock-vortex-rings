"""N = 16 (17 vortices, 50 digits): stability of the trap ring with a mobile central column at g = g*(16, r0), the H^2
ring, and the fixed-centre lower limit g_lo = s_c r0^2/2, just below and above the note's threshold r0* = 0.609657.
(Replaces the slow interval scans of ring_center.py for N = 16.) Usage: python ring_center_n16.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from ring_center import spectrum, maxre
from ring_blocks import S_window, S_trap, S_H2
mp.mp.dps = 50
N = 16; out = []
for r0 in (0.6096575*0.997, 0.6096575*1.003):
    S0 = S_trap(N, r0, 0.0); gs = (S_H2(N, r0) - S0)*r0**2/2; gl = (S_window(N, r0)[0] - S0)*r0**2/2
    eT = maxre(spectrum(N, mp.mpf(r0), mp.mpf(gs), 'trap')); eH = maxre(spectrum(N, mp.mpf(r0), 0, 'H2', center=False))
    line = (f'N = 16, r0 = {r0:.6f}: g* = {gs:.5f}, g_lo (fixed centre) = {gl:.5f};  max Re: trap + mobile centre at g* '
            f'{eT:.2e}, H^2 ring {eH:.2e}')
    print(line, flush=True); out.append(line)
open('ring_center_n16.out', 'w').write('\n'.join(out) + '\n')
