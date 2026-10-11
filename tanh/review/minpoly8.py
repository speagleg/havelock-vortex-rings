"""N=8: refine xi* from MY Cartesian Hessian (radial Fourier entry of mode 4, no candidate formula),
recover an integer relation with PSLQ, compare to the candidate minimal polynomial; irreducibility via flint."""
import mpmath as mp, flint, sympy as sp
import indep_lin
from indep_lin import linop
mp.mp.dps=120
N=8
def a_crit(xi):
    L,HF,Om,res=linop(N,xi,True); R=mp.sqrt(xi); m=4; s=0
    for j in range(N):
        for k in range(N):
            tj=2*mp.pi*j/N; tk=2*mp.pi*k/N
            vj=[R*mp.cos(tj),R*mp.sin(tj)]; vk=[R*mp.cos(tk),R*mp.sin(tk)]
            s+=mp.cos(m*(tj-tk))*sum(vj[p]*HF[2*j+p,2*k+q]*vk[q] for p in range(2) for q in range(2))
    return s/N
xs=mp.findroot(a_crit, mp.mpf('0.0673562230040812622684374765662355743940251172'), tol=mp.mpf(10)**-110)
print('xi*(8) from my Hessian, 110 digits:', mp.nstr(xs,105))
P=[1,-12,-37,-72,-84,-88,-84,-72,-37,-12,1]
print('|candidate minpoly(xi*)| =', mp.nstr(abs(mp.polyval(P,xs)),5))
rel=mp.findpoly(xs,10,maxcoeff=10**4)
print('PSLQ degree<=10 relation:',rel)
fp=flint.fmpz_poly(list(reversed(P)))
print('flint factor:',fp.factor())
print('sympy irreducible over Q:', sp.Poly(sp.Poly(list(P),sp.Symbol('x'))).is_irreducible)
print('leading/constant coeff:',P[0],P[-1],'-> unit (monic, norm = +1)')
rel=mp.pslq([xs**k for k in range(11)], maxcoeff=10**6, maxsteps=10**6)
print('pslq on powers 0..10:', rel)
