import os
import mpmath as mp, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vort import *
model=sys.argv[1]; N=int(sys.argv[2]); lo=mp.mpf(sys.argv[3]); hi=mp.mpf(sys.argv[4])
mp.mp.dps=int(sys.argv[5]) if len(sys.argv)>5 else 60
tol=mp.mpf(10)**(-(mp.mp.dps//4))
def unst(x): return max_real_part(N,x,model)[0]>tol
a,b=unst(lo),unst(hi)
assert a!=b,(a,b)
for _ in range(int(sys.argv[6]) if len(sys.argv)>6 else 80):
  m=(lo+hi)/2
  if unst(m)==a: lo=m
  else: hi=m
print(model,N,"spectral threshold xi in [",mp.nstr(lo,25),",",mp.nstr(hi,25),"]", "unstable below" if a else "unstable above", flush=True)
