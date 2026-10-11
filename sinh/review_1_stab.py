# Independent reviewer linear-stability check for the sinh model on the Poincare disk.
# H = -sum_{j<k} ln|z_j-z_k| + (N-1)/2 sum ln(1-|z_k|^2); omega = sum 4/(1-|z|^2)^2 dx^dy.
# Rotating-frame F = H - Omega*J, J = sum (1+|z|^2)/(1-|z|^2). Exact symbolic Hessian pieces (sympy),
# evaluated in mpmath at 50 digits. Omega solved from dF=0 numerically (not from the formula).
import sympy as sp, mpmath as mp
mp.mp.dps = 50
x1,y1,x2,y2,x,y,Om,Nn = sp.symbols('x1 y1 x2 y2 x y Om Nn', real=True)
pair = -sp.Rational(1,2)*sp.log((x1-x2)**2+(y1-y2)**2)
one  = (Nn-1)/2*sp.log(1-x**2-y**2) - Om*(1+x**2+y**2)/(1-x**2-y**2)
pv=[x1,y1,x2,y2]; ov=[x,y]
pg = sp.lambdify(pv,[sp.diff(pair,v) for v in pv],'mpmath')
ph = sp.lambdify(pv,[[sp.diff(pair,a,b) for b in pv] for a in pv],'mpmath')
og = sp.lambdify(ov+[Om,Nn],[sp.diff(one,v) for v in ov],'mpmath')
oh = sp.lambdify(ov+[Om,Nn],[[sp.diff(one,a,b) for b in ov] for a in ov],'mpmath')

def ring(N,xi):
    rE=mp.sqrt(xi); return [(rE*mp.cos(2*mp.pi*k/N), rE*mp.sin(2*mp.pi*k/N)) for k in range(N)]
def grad(N,pts,O):
    g=[mp.mpf(0)]*(2*N)
    for j in range(N):
        for k in range(j+1,N):
            d=pg(*pts[j],*pts[k]); g[2*j]+=d[0];g[2*j+1]+=d[1];g[2*k]+=d[2];g[2*k+1]+=d[3]
        d=og(*pts[j],O,N); g[2*j]+=d[0]; g[2*j+1]+=d[1]
    return g
def hess(N,pts,O):
    Hm=mp.zeros(2*N)
    for j in range(N):
        for k in range(j+1,N):
            h=ph(*pts[j],*pts[k]); idx=[2*j,2*j+1,2*k,2*k+1]
            for a in range(4):
                for b in range(4): Hm[idx[a],idx[b]]+=h[a][b]
        h=oh(*pts[j],O,N)
        for a in range(2):
            for b in range(2): Hm[2*j+a,2*j+b]+=h[a][b]
    return Hm
def omega_num(N,xi):
    pts=ring(N,xi)
    # radial derivative at vortex 0 is linear in Omega: solve g0(O)=0
    f=lambda O: grad(N,pts,O)[0]
    a=f(0); b=f(1); return -a/(b-a)
def lin(N,xi):
    O=omega_num(N,xi); pts=ring(N,xi); Hm=hess(N,pts,O)
    res=max(abs(v) for v in grad(N,pts,O))
    lam=4/(1-xi)**2
    # xdot = lam^-1 dF/dy, ydot = -lam^-1 dF/dx
    A=mp.zeros(2*N)
    for k in range(N):
        for c in range(2*N):
            A[2*k,c]=Hm[2*k+1,c]/lam; A[2*k+1,c]=-Hm[2*k,c]/lam
    return O,res,A
def g(N,xi):
    _,_,A=lin(N,xi); ev=mp.eig(A,left=False,right=False)
    ev=[mp.re(e*e) for e in ev if abs(e)>mp.mpf(10)**-20]
    return max(ev)
def xistar(N):
    q=mp.mpf((N-3))/4 if N%2 else mp.mpf((N-2)**2)/(4*(N-1))
    return (mp.sqrt(q)-1)/(mp.sqrt(q)+1)
def C1(N,xi): return (N-1)*(1+xi**2)/(1-xi)**2
if __name__=='__main__':
    for N in range(8,16):
        xs=xistar(N)
        # 1) Omega and spectrum vs Fourier prediction at a few radii
        maxerr=0; omerr=0
        for xi in [mp.mpf('0.03'),xs*mp.mpf('0.8'),xs*mp.mpf('1.2'),mp.mpf('0.6')]:
            O,res,A=lin(N,xi)
            omerr=max(omerr,abs(O+(N-1)*(1+xi)*(1-xi)/(8*xi)), res)
            ev=sorted([mp.re(e*e) for e in mp.eig(A,left=False,right=False)])
            k0=4*xi/(1-xi)**2
            pred=[mp.mpf(0),mp.mpf(0)]
            for m in range(1,N):
                b=mp.mpf(m*(N-m))/2; a=C1(N,xi)-b; pred+=[-a*b/k0**2]*2
            pred=sorted(pred)
            maxerr=max(maxerr,max(abs(e-p) for e,p in zip(ev,pred)))
        # 2) root of g (largest nonzero sigma^2) by bisection, no use of formula except bracket
        lo,hi=xs*mp.mpf('0.7'),xs*mp.mpf('1.3')
        assert g(N,lo)>0 and g(N,hi)<0
        for _ in range(60):
            mid=(lo+hi)/2
            if g(N,mid)>0: lo=mid
            else: hi=mid
        print(N, 'spec_err=%.2e'%float(maxerr), 'Omega/grad_err=%.2e'%float(omerr),
              'root=',mp.nstr((lo+hi)/2,20),'xi*=',mp.nstr(xs,20),'diff=%.1e'%float(abs((lo+hi)/2-xs)))
