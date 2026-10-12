import re,glob
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps=30
def c4(n):
    if n%2==0:
        m=n-2; return Fr((n-2)*(m**4+40*m*m+48*m+16),192*n*(n-1))
    return Fr((n-1)*(n-3)*(n*n-9*n+38),128*n)
print(" n   xi*            E4 (extrap)        xi^2*E4          closed form      rel.diff")
for n in range(8,22):
    try: txt=open(f"out/ls_{n}.out").read()
    except: continue
    m=re.search(r"angle 0/8 pi.*\['([^']+)', '([^']+)', '([^']+)'\]",txt)
    if not m: continue
    v2,_,v1=[mp.mpf(x) for x in m.groups()]
    E4=(4*v1-v2)/3
    q=mp.mpf((n-3))/4 if n%2 else mp.mpf((n-2)**2)/(4*(n-1))
    ch=mp.sqrt(q); xi=(ch-1)/(ch+1)
    val=xi**2*E4; cf=c4(n)
    print(f"{n:2d}  {mp.nstr(xi,12):14s} {mp.nstr(E4,14):18s} {mp.nstr(val,12):16s} {str(cf):16s} {mp.nstr(abs(val/ (mp.mpf(cf.numerator)/cf.denominator)-1),2)}")
