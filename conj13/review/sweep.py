"""Level-set c4 (independent KKT code) vs candidate closed forms, all n given on argv."""
import mpmath as mp, sys
from c4_levelset import c4
mp.mp.dps=45
def bn(n):
    k=mp.mpf((n*n)//4)/(2*(n-1)); rts=mp.polyroots([1-k,-2*k,1-k]) if k!=1 else [mp.mpf(0)]
    r=[mp.re(x) for x in rts if -1<mp.re(x)<=1 and abs(mp.im(x))<1e-30]
    return r[0]
for n in map(int,sys.argv[1:]):
    lam=bn(n) if n!=7 else mp.mpf(0)
    v,_=c4(n,lam,'level',[mp.mpf(j)/200 for j in range(1,7)])
    m=n-2
    F=mp.mpf((n-2)*(m**4+40*m**2+48*m+16))/(192*n*(n-1)) if n%2==0 else mp.mpf((n-1)*(n-3)*(n*n-9*n+38))/(128*n)
    print(n,'b_n=',mp.nstr(lam,15),' c4=',mp.nstr(v,22),' formula=',mp.nstr(F,22),' diff=',mp.nstr(v-F,3),flush=True)
