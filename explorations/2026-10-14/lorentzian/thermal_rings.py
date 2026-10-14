"""Ring thresholds of the tanh model read as temperatures (Euclidean AdS2 dictionary, ads2_dictionary.py).
A ring of radius r is the Euclidean thermal circle of an AdS2 observer with acceleration a = coth r (AdS radius 1),
temperature T = 1/(2 pi sinh r) = sqrt(a^2 - 1)/(2 pi); the N vortices sit at Euclidean-time spacing beta/N = 2 pi sinh r/N.
The certified thresholds (rings stable for r >= r*(N), note Thm 7/16) become: stable iff T <= T*_N, with
N T*_N = 1/l_N -> 1/l* = 1/0.8439085 = 1.18497 (the horocycle / extremal limit). N <= 7: stable at every T.
Usage: python thermal_rings.py   (reads ../rows/circle_branch.out)"""
import numpy as np
cir = np.loadtxt('../rows/circle_branch.out', skiprows=1)
print(' N   a* = coth r*   T*_N (AdS units)   beta/N = l_N   N T*_N')
for N, k, l in cir:
    sinh_r = 1/np.sqrt(k**2 - 1); T = 1/(2*np.pi*sinh_r)
    print(f'{int(N):3d}   {k:10.6f}    {T:12.6f}        {l:9.6f}    {N*T:8.5f}')
print(f'limit N -> infinity: N T* -> 1/l* = {1/0.8439085312925335:.5f}')
