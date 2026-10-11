"""P2 ratio test (2026-10-11): normalisation-free comparison of MT (3.9) printed beta,gamma,delta (n=5,7,9,11)
with our level-set and linear-slice c4, each divided by the unslaved quartic q4 along the kernel vector.
n=11 linear-slice value 5.795598029371687691946575 is from `python c4_levelset.py 11 "-3+2*sqrt(2)" 40`."""
import sys, mpmath as mp
sys.path.insert(0,'.')
from c4_levelset import H_grad_hess, J_grad_hess, setup, vectors
mp.mp.dps=50
MT={5:(37.1,-9.04,15.8),7:(21/2,63/4,1071/8),9:(6.9,13.8,182.4),11:(12.0,29.7,555.1)}
lev={5:mp.mpf('0.225'),7:mp.mpf(9)/14,9:mp.mpf('1.58333333333333333'),11:mp.mpf('3.40909090909090909')}
lin={11:mp.mpf('5.795598029371687691946575'),5:mp.mpf('0.204657287525381'),7:mp.mpf(6)/7,9:mp.mpf('2.57814629687882')}
for n in (5,7,9,11):
    # threshold (n-7)(1+l^2)+2(n+1)l=0, root with |l|<1
    if n==7: lam=mp.mpf(0)
    else:
        a=n-7; b=2*(n+1); r=[(-b+s*mp.sqrt(b*b-4*a*a))/(2*a) for s in (1,-1)]; lam=[x for x in r if abs(x)<1][0]
    q0=setup(n,lam); e,e2,t=vectors(n,q0)
    _,gH0,_=H_grad_hess(q0,n,lam); J0,G0,_=J_grad_hess(q0,n,lam)
    om=mp.fsum(gH0[i]*G0[2][i] for i in range(2*n))/mp.fsum(G0[2][i]**2 for i in range(2*n))
    def Hh(s):
        q=[q0[i]+s*e[i] for i in range(2*n)]
        v,_,_=H_grad_hess(q,n,lam); J,_,_=J_grad_hess(q,n,lam); return v-om*J[2]
    h0=Hh(0)
    fs=[(Hh(mp.mpf(k)/400)+Hh(-mp.mpf(k)/400)-2*h0)/2/(mp.mpf(k)/400)**4 for k in range(1,7)]
    xs=[(mp.mpf(k)/400)**2 for k in range(1,7)]; P=list(fs)
    for k in range(1,6):
        for i in range(6-k): P[i]=(xs[i+k]*P[i]-xs[i]*P[i+1])/(xs[i+k]-xs[i])
    q4=P[0]; b,g,d=MT[n]
    print(f"n={n} lam={mp.nstr(lam,8)} q4={mp.nstr(q4,10)} | ours: level/q4={mp.nstr(lev[n]/q4,5)}"
          + (f" lin/q4={mp.nstr(lin[n]/q4,5)}" if n in lin else "")
          + f" | MT: (d-g^2/b)/d={(d-g*g/b)/d:.4f} (d-3g^2/b)/d={(d-3*g*g/b)/d:.4f}")
