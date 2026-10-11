"""NUMERICAL: for N=8..120, find the largest xi in (0,1) where any a_m<0 (m=1..N-1), check it is mode floor(N/2),
and check b_m>0 on a grid. Uses candidate block formulas (validated against my Hessian in compare_blocks.py)."""
import sys, mpmath as mp
sys.path.insert(0,'..'); from blocks import blocks
mp.mp.dps=30
bad=[]
for N in list(range(8,41))+[50,64,81,100,120]:
    grid=[mp.mpf(i)/400 for i in range(1,400)]
    last={}
    minb=mp.inf
    for xi in grid:
        for m in range(1,N//2+1):
            a,b,g,E=blocks(N,xi,m); minb=min(minb,b)
            if a<0: last[m]=xi
    top=max(last.values()); mtop=[m for m in last if last[m]==top]
    if N//2 not in mtop or minb<=0: bad.append((N,mtop,minb))
    if N in (8,9,40,120): print(N,'last unstable grid pt',mp.nstr(top,4),'modes',mtop,'min b',mp.nstr(minb,4))
print('exceptions:',bad)
