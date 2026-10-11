import mpmath as mp
from indep_lin import maxre
mp.mp.dps=40
for N in range(3,8):
    w=max(maxre(N,mp.mpf(i)/100) for i in list(range(1,100,3))+[99])
    print('N=%d max|Re eig| over xi grid (0.01..0.99, 34 pts): %s'%(N,mp.nstr(w,3)))
