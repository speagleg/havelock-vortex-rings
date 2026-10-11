"""P2: compare MT (3.9) numbers for n=7 (beta=21/2pi, gamma=63/4pi, delta=1071/8pi) with an independent
computation.  Normalisation-free invariant: (Lemma-8 slaved quartic)/(unslaved quartic D a^4).
D = 24 * [eps^4 coeff of Hhat(q0+eps a)] (a = unit kernel vector, lies in Ker dJ);
the linear-slice slaved quartic is c4_lin = 6/7 (from c4_levelset.py, 'linear').
"""
import mpmath as mp
from c4_levelset import H_grad_hess, J_grad_hess, setup, vectors
mp.mp.dps=60
n=7; lam=mp.mpf(0); q0=setup(n,lam); e,e2,t=vectors(n,q0)
_,gH0,_=H_grad_hess(q0,n,lam); J0,G0,_=J_grad_hess(q0,n,lam)
om=mp.fsum(gH0[i]*G0[2][i] for i in range(2*n))/mp.fsum(G0[2][i]**2 for i in range(2*n))
def Hh(s):
    q=[q0[i]+s*e[i] for i in range(2*n)]
    v,_,_=H_grad_hess(q,n,lam); J,_,_=J_grad_hess(q,n,lam); return v-om*J[2]
h0=Hh(0)
fs=[ (Hh(mp.mpf(k)/400)+Hh(-mp.mpf(k)/400)-2*h0)/2/(mp.mpf(k)/400)**4 for k in range(1,7)]
xs=[(mp.mpf(k)/400)**2 for k in range(1,7)]; P=list(fs)
for k in range(1,6):
    for i in range(6-k): P[i]=(xs[i+k]*P[i]-xs[i]*P[i+1])/(xs[i+k]-xs[i])
q4=P[0]; print('unslaved quartic coeff q4 =',mp.nstr(q4,20), mp.identify(q4))
lin=mp.mpf(6)/7; lev=mp.mpf(9)/14
print('ratio lin/q4 =',mp.nstr(lin/q4,20),' ratio level/q4 =',mp.nstr(lev/q4,20))
b,g,d=mp.mpf(21)/2,mp.mpf(63)/4,mp.mpf(1071)/8
for k in (1,3,12): print(' MT (delta - %d gamma^2/beta)/delta ='%k, mp.nstr((d-k*g*g/b)/d,20))
print(' 8/17 =',mp.nstr(mp.mpf(8)/17,20),' 6/17 =',mp.nstr(mp.mpf(6)/17,20))
