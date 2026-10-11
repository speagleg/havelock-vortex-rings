"""Independent second-reader computation of the effective quartic c4 on J^{-1}(mu)/SO(2).
Model (MT (1.12)): chart C, sigma=1+lam|z|^2, omega=sum 2dx^dy/sigma^2,
H = -sum_{i<j} log(|zi-zj|^2/(sigma_i sigma_j))  (no 1/4pi), J = sum (x/sigma, y/sigma, |z|^2/sigma).
Ring z_k = exp(2 pi i k/n), r0 = 1, lam = b_n.
f(eps) = min{ H(q)-H(q0) : e.(q-q0)=eps, t.(q-q0)=0, e2.(q-q0)=0 (odd n), J(q)=J(q0) }
e = unit Euclidean vector, radial Fourier mode l=floor(n/2): dz_k = c cos(l phi_k) e^{i phi_k}.
Solved by Newton on the KKT system at high precision; c4 by polynomial extrapolation in eps.
Also computes the LINEAR-SLICE quantity (MT Lemma 8 on the affine slice q0 + ker dJ cap t^perp)
for comparison:  f_lin(eps) = min{ Hhat(q) : q-q0 in ker dJ(q0), e.(q-q0)=eps, t.(q-q0)=0 },
Hhat = H - omega.J (omega the multiplier at the ring).
"""
import mpmath as mp, sys
mp.mp.dps = int(sys.argv[3]) if len(sys.argv)>3 else 60

def setup(n, lam):
    q0=[]
    for k in range(n):
        ph=2*mp.pi*k/n; q0 += [mp.cos(ph), mp.sin(ph)]
    return q0

def H_grad_hess(q, n, lam):
    N=2*n; g=[mp.mpf(0)]*N; Hs=mp.zeros(N,N); val=mp.mpf(0)
    for i in range(n):
        xi,yi=q[2*i],q[2*i+1]
        for j in range(i+1,n):
            xj,yj=q[2*j],q[2*j+1]
            dx,dy=xi-xj,yi-yj; d2=dx*dx+dy*dy
            val -= mp.log(d2)
            gx,gy = -2*dx/d2, -2*dy/d2       # d(-log d2)/d(xi,yi)
            for (a,b) in ((i,i),(j,j)): pass
            g[2*i]+=gx; g[2*i+1]+=gy; g[2*j]-=gx; g[2*j+1]-=gy
            # hessian of -log d2 wrt (dx,dy): -2/d2 I + 4 d d^T/d2^2
            hxx=-2/d2+4*dx*dx/d2**2; hyy=-2/d2+4*dy*dy/d2**2; hxy=4*dx*dy/d2**2
            blk=[[hxx,hxy],[hxy,hyy]]
            for a in range(2):
                for b in range(2):
                    Hs[2*i+a,2*i+b]+=blk[a][b]; Hs[2*j+a,2*j+b]+=blk[a][b]
                    Hs[2*i+a,2*j+b]-=blk[a][b]; Hs[2*j+a,2*i+b]-=blk[a][b]
    for i in range(n):  # +(n-1) log sigma_i
        x,y=q[2*i],q[2*i+1]; s=1+lam*(x*x+y*y); c=n-1
        val += c*mp.log(s)
        g[2*i]+=c*2*lam*x/s; g[2*i+1]+=c*2*lam*y/s
        Hs[2*i,2*i]+=c*(2*lam/s-4*lam**2*x*x/s**2)
        Hs[2*i+1,2*i+1]+=c*(2*lam/s-4*lam**2*y*y/s**2)
        Hs[2*i,2*i+1]+=c*(-4*lam**2*x*y/s**2); Hs[2*i+1,2*i]+=c*(-4*lam**2*x*y/s**2)
    return val,g,Hs

def J_grad_hess(q,n,lam):
    """3 components; each a sum of single-vortex functions."""
    N=2*n; vals=[mp.mpf(0)]*3; grads=[[mp.mpf(0)]*N for _ in range(3)]; hess=[mp.zeros(N,N) for _ in range(3)]
    for i in range(n):
        x,y=q[2*i],q[2*i+1]; u=x*x+y*y; s=1+lam*u
        # F1=x/s, F2=y/s, F3=u/s
        vals[0]+=x/s; vals[1]+=y/s; vals[2]+=u/s
        sx,sy=2*lam*x,2*lam*y
        grads[0][2*i]+=1/s - x*sx/s**2; grads[0][2*i+1]+= -x*sy/s**2
        grads[1][2*i]+= -y*sx/s**2; grads[1][2*i+1]+=1/s - y*sy/s**2
        grads[2][2*i]+=2*x/s**2; grads[2][2*i+1]+=2*y/s**2
        # hessians (by sympy-free hand formulas)
        # F1 = x/s
        F1xx=-2*sx/s**2 - x*2*lam/s**2 + 2*x*sx*sx/s**3
        F1xy=-sy/s**2 + 2*x*sx*sy/s**3
        F1yy=-x*2*lam/s**2 + 2*x*sy*sy/s**3
        F2yy=-2*sy/s**2 - y*2*lam/s**2 + 2*y*sy*sy/s**3
        F2xy=-sx/s**2 + 2*y*sx*sy/s**3
        F2xx=-y*2*lam/s**2 + 2*y*sx*sx/s**3
        F3xx=2/s**2 - 4*x*sx/s**3; F3yy=2/s**2-4*y*sy/s**3; F3xy=-4*x*sy/s**3
        for m,(a,b,c) in enumerate(((F1xx,F1xy,F1yy),(F2xx,F2xy,F2yy),(F3xx,F3xy,F3yy))):
            hess[m][2*i,2*i]+=a; hess[m][2*i,2*i+1]+=b; hess[m][2*i+1,2*i]+=b; hess[m][2*i+1,2*i+1]+=c
    return vals,grads,hess

def vectors(n,q0):
    l=n//2; N=2*n
    e=[mp.mpf(0)]*N; e2=[mp.mpf(0)]*N; t=[mp.mpf(0)]*N
    for k in range(n):
        ph=2*mp.pi*k/n
        e[2*k]=mp.cos(l*ph)*mp.cos(ph); e[2*k+1]=mp.cos(l*ph)*mp.sin(ph)
        e2[2*k]=mp.sin(l*ph)*mp.cos(ph); e2[2*k+1]=mp.sin(l*ph)*mp.sin(ph)
        t[2*k]=-mp.sin(ph); t[2*k+1]=mp.cos(ph)
    def nz(v):
        s=mp.sqrt(mp.fsum(a*a for a in v)); return [a/s for a in v]
    return nz(e), nz(e2), nz(t)

def solve(n,lam,eps,mode,q_init=None):
    N=2*n; q0=setup(n,lam); e,e2,t=vectors(n,q0)
    J0,G0,_=J_grad_hess(q0,n,lam)
    lin=[e,t]+([e2] if n%2 else [])
    rhs_lin=[eps,0]+([0] if n%2 else [])
    # omega for Hhat: grad H = omega * grad J3 at ring
    _,gH0,_=H_grad_hess(q0,n,lam)
    omega=mp.fsum(gH0[i]*G0[2][i] for i in range(N))/mp.fsum(G0[2][i]**2 for i in range(N))
    q=list(q_init) if q_init else [q0[i]+eps*e[i] for i in range(N)]
    nc=len(lin)+3; mult=[mp.mpf(0)]*nc
    if mode=='level': mult[len(lin)+2]=omega
    for it in range(60):
        val,g,Hs=H_grad_hess(q,n,lam)
        if mode=='level':
            Jv,Jg,Jh=J_grad_hess(q,n,lam)
            cons=[mp.fsum(lin[a][i]*(q[i]-q0[i]) for i in range(N))-rhs_lin[a] for a in range(len(lin))]+[Jv[m]-J0[m] for m in range(3)]
            C=lin+Jg
            L=Hs.copy()
            for m in range(3): L-= mult[len(lin)+m]*Jh[m]
            obj_g=g
        else:  # linear slice, Hhat
            Jv,Jg,Jh=J_grad_hess(q,n,lam)
            obj_g=[g[i]-omega*Jg[2][i] for i in range(N)]
            L=Hs-omega*Jh[2]
            cons=[mp.fsum(lin[a][i]*(q[i]-q0[i]) for i in range(N))-rhs_lin[a] for a in range(len(lin))]+[mp.fsum(G0[m][i]*(q[i]-q0[i]) for i in range(N)) for m in range(3)]
            C=lin+G0
        K=mp.zeros(N+nc,N+nc); r=mp.zeros(N+nc,1)
        for i in range(N):
            for j in range(N): K[i,j]=L[i,j]
            r[i]=obj_g[i]-mp.fsum(mult[a]*C[a][i] for a in range(nc))
        for a in range(nc):
            for i in range(N): K[i,N+a]=-C[a][i]; K[N+a,i]=C[a][i]
            r[N+a]=cons[a]
        d=mp.lu_solve(K,-r)
        q=[q[i]+d[i] for i in range(N)]; mult=[mult[a]+d[N+a] for a in range(nc)]
        if mp.norm(d)<mp.mpf(10)**(-mp.mp.dps+8): break
    val,_,_=H_grad_hess(q,n,lam); val0,_,_=H_grad_hess(q0,n,lam)
    if mode=='level': return val-val0, q
    Jv,_,_=J_grad_hess(q,n,lam)
    return (val-omega*Jv[2])-(val0-omega*J0[2]), q

def c4(n,lam,mode,epss):
    fs=[]
    for s in epss:
        fp,_=solve(n,lam,s,mode); fm,_=solve(n,lam,-s,mode)
        fs.append((fp+fm)/2/s**4)
    # fit fs = c4 + c6 s^2 + c8 s^4 ... (Neville at s^2 -> 0)
    xs=[s**2 for s in epss]; P=list(fs); m=len(xs)
    for k in range(1,m):
        for i in range(m-k):
            P[i]=(xs[i+k]*P[i]-xs[i]*P[i+1])/(xs[i+k]-xs[i])
    return P[0], fs

if __name__=='__main__':
    n=int(sys.argv[1]); lamstr=sys.argv[2]
    lam = mp.mpf(0) if lamstr=='0' else mp.mpf(eval(lamstr,{'sqrt':mp.sqrt}))
    epss=[mp.mpf(k)/200 for k in range(1,7)]
    for mode in ('level','linear'):
        v,fs=c4(n,lam,mode,epss)
        print(n, lamstr, mode, 'c4 =', mp.nstr(v,25), ' raw f/eps^4:', [mp.nstr(f,10) for f in fs[:3]])
        if mode=='level': print('   identify:', mp.identify(v), ' pslq w/ sqrt7:', mp.identify(v,['sqrt(7)']))
