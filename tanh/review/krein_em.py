"""Odd-N EM/spectral split on MY full Cartesian Hessian: in the window (xi*, xi_EM) count negative
eigenvalues of Hess F (coordinate-invariant signature), check spectrum imaginary, and Krein signs
sign(v^H HessF v) for each eigenvector of L. Also identify the Fourier mode of the unstable eigenvector below xi*."""
import mpmath as mp
from indep_lin import linop
mp.mp.dps=50
def analyse(N, xi):
    xi=mp.mpf(xi); L,HF,Om,res=linop(N,xi,True)
    hev=mp.eigsy(HF)[0]
    neg=sum(1 for e in hev if e< -mp.mpf(10)**-30); zero=sum(1 for e in hev if abs(e)<=mp.mpf(10)**-30)
    E,V=mp.eig(L)
    mr=max(abs(mp.re(e)) for e in E)
    kre=[]
    for i,e in enumerate(E):
        v=V[:,i]; q=mp.re((v.H*HF*v)[0])
        if abs(mp.im(e))>1e-20: kre.append((mp.nstr(mp.im(e),8),'+' if q>0 else '-'))
    # dominant Fourier mode of each eigenvector
    def mode(v):
        best=None
        for m in range(N):
            s=0
            for comp in (0,1):
                # complex coordinate z_j perturbation u_j = v_x + i v_y, rotate to local frame
                pass
        return None
    return neg, zero, mr, sorted(kre)
def unstable_mode(N, xi):
    xi=mp.mpf(xi); L,HF,Om,res=linop(N,xi,True)
    E,V=mp.eig(L); i=max(range(len(E)),key=lambda k: mp.re(E[k])); v=V[:,i]
    u=[(v[2*j]+1j*v[2*j+1])*mp.expj(-2*mp.pi*j/N) for j in range(N)]  # local-frame complex perturbation
    ub=[(v[2*j]-1j*v[2*j+1])*mp.expj(2*mp.pi*j/N) for j in range(N)]
    w=[]
    for m in range(N):
        w.append(abs(sum(u[j]*mp.expj(-2*mp.pi*j*m/N) for j in range(N)))**2+abs(sum(ub[j]*mp.expj(-2*mp.pi*j*m/N) for j in range(N)))**2)
    tot=sum(w); return mp.re(E[i]), {m:mp.nstr(w[m]/tot,3) for m in range(N) if w[m]/tot>1e-6}
xs={9:mp.mpf('0.11287106323660274043969401905650541554651262119976'),11:mp.mpf('0.20460117688266929063553583315480299815006661809463'),
    8:mp.mpf('0.06735622300408126226843747656623557439402511756974'),10:mp.mpf('0.16611424152623911808175869978489718411277213101352')}
em={9:mp.mpf('0.112871095435130867829042255257'),11:mp.mpf('0.204601301313693124209019442065')}
for N in (9,11):
    for lab,xi in [('below xi*',xs[N]-mp.mpf('1e-6')),('in window',(xs[N]+em[N])/2),('above xi_EM',em[N]+mp.mpf('1e-6'))]:
        neg,zero,mr,kre=analyse(N,xi)
        print('N=%d %-11s neg eig(HessF)=%d zero=%d max|Re|=%s'%(N,lab,neg,zero,mp.nstr(mr,3)))
        if lab=='in window': print('   Krein (Im sigma, sign):',kre)
for N in (8,9,10,11):
    print('N=%d unstable eigvec at xi*-1e-3: Re=%s, Fourier weights'%(N,mp.nstr(unstable_mode(N,xs[N]-mp.mpf('1e-3'))[0],5)),unstable_mode(N,xs[N]-mp.mpf('1e-3'))[1])
