# Independent check of the tanh Fourier blocks (blocks.py) against the exact Cartesian Hessian / Jacobian.
# Model: H = -sum_{j<k} ln tanh(d_jk/2) = sum_{j<k} [-ln|z_j-z_k| + ln|1 - z_j conj(z_k)|],
# area form 4 dx^dy/(1-|z|^2)^2, moment map J = sum (1+|z|^2)/(1-|z|^2), F = H - Om*J.
# Pair/one-body second derivatives exact (sympy), evaluated in mpmath; Om from dF=0 (linear in Om).
import sympy as sp, mpmath as mp, sys
from blocks import blocks
mp.mp.dps = 50
x1,y1,x2,y2,x,y,Om = sp.symbols('x1 y1 x2 y2 x y Om', real=True)
pair = -sp.Rational(1,2)*sp.log((x1-x2)**2+(y1-y2)**2) + sp.Rational(1,2)*sp.log((1-x1*x2-y1*y2)**2+(x1*y2-y1*x2)**2)
one = -Om*(1+x**2+y**2)/(1-x**2-y**2)
pv=[x1,y1,x2,y2]; ov=[x,y]
pg=sp.lambdify(pv,[sp.diff(pair,v) for v in pv],'mpmath')
ph=sp.lambdify(pv,[[sp.diff(pair,a,b) for b in pv] for a in pv],'mpmath')
og=sp.lambdify(ov+[Om],[sp.diff(one,v) for v in ov],'mpmath')
oh=sp.lambdify(ov+[Om],[[sp.diff(one,a,b) for b in ov] for a in ov],'mpmath')
def ring(N,xi):
    r=mp.sqrt(xi); return [(r*mp.cos(2*mp.pi*k/N), r*mp.sin(2*mp.pi*k/N)) for k in range(N)]
def grad(N,P,O):
    g=[mp.mpf(0)]*(2*N)
    for j in range(N):
        for k in range(j+1,N):
            d=pg(*P[j],*P[k]); g[2*j]+=d[0];g[2*j+1]+=d[1];g[2*k]+=d[2];g[2*k+1]+=d[3]
        d=og(*P[j],O); g[2*j]+=d[0]; g[2*j+1]+=d[1]
    return g
def hess(N,P,O):
    Hm=mp.zeros(2*N)
    for j in range(N):
        for k in range(j+1,N):
            h=ph(*P[j],*P[k]); idx=[2*j,2*j+1,2*k,2*k+1]
            for a in range(4):
                for b in range(4): Hm[idx[a],idx[b]]+=h[a][b]
        h=oh(*P[j],O)
        for a in range(2):
            for b in range(2): Hm[2*j+a,2*j+b]+=h[a][b]
    return Hm
def run(N,xi):
    xi=mp.mpf(xi); P=ring(N,xi)
    f=lambda O: grad(N,P,O)[0]; a0=f(0); O=-a0/(f(1)-a0)
    res=max(abs(v) for v in grad(N,P,O))
    Hm=hess(N,P,O)
    # Omega check: E = Om*4xi/(1-xi)^2 ; here F=H-Om J with J'=4u/(1-u)^2
    E=blocks(N,xi,1)[3]; Oerr=abs(O*4*xi/(1-xi)**2-E)
    # log-polar Hessian: T_j = [[x,-y],[y,x]]
    T=mp.zeros(2*N)
    for j in range(N):
        X,Y=P[j]; T[2*j,2*j]=X;T[2*j+1,2*j]=Y;T[2*j,2*j+1]=-Y;T[2*j+1,2*j+1]=X
    L=T.T*Hm*T
    w=mp.exp(2j*mp.pi/N); err=0
    for m in range(1,N+1):
        M=[[sum(w**(m*k)*L[a,2*k+b] for k in range(N)) for b in range(2)] for a in range(2)]
        a,b,g,_=blocks(N,xi,m)
        err=max(err,abs(M[0][0]-a),abs(M[1][1]-b),abs(M[0][1]-1j*g),abs(M[0][1]-M[1][0].conjugate()),abs(mp.re(M[0][1])))
    # full Jacobian spectrum
    lam=4/(1-xi)**2; A=mp.zeros(2*N)
    for k in range(N):
        for c in range(2*N):
            A[2*k,c]=Hm[2*k+1,c]/lam; A[2*k+1,c]=-Hm[2*k,c]/lam
    ev=mp.eig(A,left=False,right=False)
    kap=4*xi/(1-xi)**2; pred=[mp.mpc(0),mp.mpc(0)]
    for m in range(1,N):
        a,b,g,_=blocks(N,xi,m); s=mp.sqrt(mp.mpc(-a*b))
        pred+=[(-1j*g+s)/kap,(-1j*g-s)/kap]
    used=[False]*len(pred); serr=0
    for e in ev:
        i=min((i for i in range(len(pred)) if not used[i]),key=lambda i:abs(e-pred[i]))
        used[i]=True; serr=max(serr,abs(e-pred[i]))
    return res,Oerr,err,serr,max(mp.re(e) for e in ev)
if __name__=='__main__':
    print('N xi  |gradF|  |Omega-formula|  max|symbol-blocks|  max|spectrum-pred|  maxRe(sigma)')
    for N in [3,5,7,8,9,10,12,13,16]:
        for xi in ['0.02','0.09','0.3','0.75']:
            r=run(N,xi)
            print(N,xi,*['%.1e'%float(v) for v in r[:4]],mp.nstr(r[4],8),flush=True)
