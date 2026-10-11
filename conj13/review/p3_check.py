"""P3: MT's printed even-n t^4 coefficient  n/(768 pi r0^4 sigma^4) * T  vs independent straight-line quartic.
My units: f = 4pi * Hhat_MT, unit vector e = (1/sqrt n)(-1)^j radial, so eps = sqrt(n) t:
  q4_mine * eps^4 = 4 pi * n T/(768 pi sigma^4) t^4  =>  T_pred = q4_mine * n * 192 * sigma^4 / ... (see code)."""
import mpmath as mp
from c4_levelset import c4
mp.mp.dps=50
for n,lam in ((8,-(8-3*mp.sqrt(7))),(10,None),(6,None),(4,None)):
    if lam is None:
        # solve (1+b^2)/(1+b)^2 = n^2/(8(n-1)) root in (-1,1]: (1-k)b^2 -2k b + (1-k)=0
        k=mp.mpf(n*n)/(8*(n-1)); rts=mp.polyroots([1-k,-2*k,1-k]); lam=[r for r in rts if -1<mp.re(r)<=1][0].real
    ql,_=c4(n,lam,'linear',[mp.mpf(j)/200 for j in range(1,6)])
    # 4pi*f_MT coefficient of t^4 = n T/(192 sigma^4);  ql*eps^4 = ql n^2 t^4
    for name,sig in (('1+lam',1+lam),('1-lam',1-lam)):
        Tpred=ql*n*n*192*sig**4/n
        T=(n-2)*(n**3+2*n**2-12*n+24)*sig**4+24*(n-1)*lam*(19*lam**3-54*lam**2+43*lam-4)
        print(n,'lam=',mp.nstr(lam,12),'sigma=',name,' T_MT=',mp.nstr(T,12),' T_indep=',mp.nstr(Tpred,12),' rel diff',mp.nstr((T-Tpred)/Tpred,4))
