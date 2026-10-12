import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from vort import *
mp.mp.dps=50
N=9; xs=mp.mpf('0.1128710632366027404397')
def neg_count(xi):
    rho=mp.sqrt(xi); L,Q,gJ,res=linearization(N,rho,'tanh'); n2=2*N
    Z=ring(N,rho); g1=[];g2=[]
    for z in Z:
        x,y=mp.re(z),mp.im(z); D=1-x*x-y*y
        g1+=[4/D+8*x*x/D**2, 8*x*y/D**2]; g2+=[8*x*y/D**2, 4/D+8*y*y/D**2]
    cons=[mp.matrix(gJ),mp.matrix(g1),mp.matrix(g2)]
    basis=[]
    for c in cons:
        v=c
        for b in basis: v=v-(b.T*v)[0]*b
        basis.append(v/mp.norm(v))
    nb=[]
    for c in range(n2):
        v=mp.matrix(n2,1); v[c]=1
        for b in basis+nb: v=v-(b.T*v)[0]*b
        if mp.norm(v)>1e-10: nb.append(v/mp.norm(v))
    B=mp.matrix(n2,len(nb))
    for j,b in enumerate(nb):
        for i in range(n2): B[i,j]=b[i]
    ev=mp.eigsy(B.T*Q*B,eigvals_only=True)
    ev=sorted(ev)
    lam=mp.eig(L,left=False,right=False)
    return sum(1 for e in ev if e< -mp.mpf('1e-20')), mp.nstr(ev[0],5), mp.nstr(max(abs(mp.re(l)) for l in lam),3)
for d in ['-1e-8','1e-9','1.6e-8','3.0e-8','3.4e-8','1e-7']:
    print('xi*+',d,' (#neg reduced-Hessian eigs, min eig, max|Re lambda|):',neg_count(xs+mp.mpf(d)))
