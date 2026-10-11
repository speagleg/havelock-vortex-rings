"""Odd n: symbolic-in-n derivation of the level-set quartic c4 (second reader, independent route).

Coordinates: canonical (theta_k, w_k), w = u/(1+lam u) (u=|z|^2): J3 = sum w_k is LINEAR,
J1 + i J2 = sum sqrt(P(w_k)) e^{i theta_k}, P(w) = w(1-lam w) (momentum map: momentum_check.py).
Pair term: |zi-zj|^2/(sig_i sig_j) = w_i + w_j - 2 lam w_i w_j - 2 sqrt(P_i P_j) cos(theta_i-theta_j).
Ring r0 = 1: w0 = 1/(1+lam).  l = (n-1)/2,  H = -sum_{i<j} log(...).
Ansatz on J^{-1}(mu)/SO(2) to the order that affects A^4:
  delta_k = A cos(l phi_k) + A^2 p cos(phi_k),  psi_k = A^2 s sin(phi_k).
Reasons: J3 linear => mode 0 of delta is exactly 0; mode-0 psi = SO(2) gauge; for n >= 5 only modes 0, +-1
meet a^2 at cubic order (2l = n-1 = -1 mod n); reflection k -> -k fixes a, so only cos-delta/sin-psi in mode 1.
A^3, A^4 corrections do not enter f at O(A^4): dH vanishes on Ker dJ3 (dH = omega dJ3, J3 linear) and
a spans a null direction of the full Hessian (radial/tangential decouple at equal radii).
J1 = 0 at O(A^2) is one linear equation in (p,s); J2 = 0 by reflection.
Offset sums: t = e^{i pi m/n}, sigma = (-1)^m = t^n, e^{i l phi_m} = sigma/t, e^{i phi_m} = t^2,
  sum_{m=1}^{n-1} g = 1/2 sum_{t^{2n}=1, t != +-1} g   (g invariant under (t,sigma)->(-t,-sigma), asserted),
  sum_{t^n=s, R finite} R(t) = -Res_{t=1} - Res_{t=-1} - Res_0 - Res_inf of R(t) n t^{n-1}/(t^n - s);
  Res_0 = 0 and Res_inf = -n [t^0 coeff of R at inf] once n exceeds R's Laurent degree at 0 and inf (printed).
"""
import sympy as sp, time, sys, pickle
T0=time.time()
def log(*a): print('[%.0fs]'%(time.time()-T0),*a,flush=True)
A,lam,p,s,X,Y,U,V,t,sig,n,tau=sp.symbols('A lam p s X Y U V t sigma n tau')
c,sn,x,y,z=sp.symbols('c sn x y z')
Lp=sp.Symbol('Lp',positive=True)   # Lp = 1+lam > 0 (lam in (-1,1]); lets sqrt(P0^2) = P0 simplify
lamL=Lp-1
w0=1/Lp
P=lambda w: w*(1-lamL*w)
def sqrtP(w): return sp.sqrt(P(w))
D=(w0+x)+(w0+y)-2*lamL*(w0+x)*(w0+y)-2*sqrtP(w0+x)*sqrtP(w0+y)*(c*sp.cos(z)-sn*sp.sin(z))
G=-sp.log(D)
# Taylor coefficients needed: x^i y^j z^k with i+j+2k <= 4
coef={}
for k in range(3):
    for i in range(5):
        for j in range(5):
            if i+j+2*k<=4 and i+j+k>=1:
                d=sp.diff(G,x,i,y,j,z,k).subs({x:0,y:0,z:0})
                d=sp.factor(sp.simplify(d))
                assert not d.has(sp.sqrt(Lp)) , d
                d=sp.factor(sp.cancel(d.subs(Lp,1+lam)))
                coef[(i,j,k)]=d/(sp.factorial(i)*sp.factorial(j)*sp.factorial(k))
log('Taylor coefficients done:',len(coef))
xa_i,xb_i,xa_j,xb_j,zz=sp.symbols('xa_i xb_i xa_j xb_j zz')
Fpoly=sum(co*(A*xa_i+A**2*xb_i)**i*(A*xa_j+A**2*xb_j)**j*(A**2*zz)**k for (i,j,k),co in coef.items())
Fpoly=sp.expand(Fpoly)
Fk={k:Fpoly.coeff(A,k) for k in (1,2,3,4)}
subsd={xa_i:(X+1/X)/2, xa_j:(X*U+1/(X*U))/2, xb_i:p*(Y+1/Y)/2, xb_j:p*(Y*V+1/(Y*V))/2,
       zz:s*((Y-1/Y)-(Y*V-1/(Y*V)))/(2*sp.I)}
monos=[xa_i,xa_j,xb_i,xb_j,zz]
avg_cache={}
def avg_mono(m):
    if m in avg_cache: return avg_cache[m]
    e=sp.expand(m.subs(subsd)); out=0
    for term in sp.Add.make_args(e):
        # handle negative powers
        pw=term.as_powers_dict(); a=pw.get(X,0); b=pw.get(Y,0)
        if a==2*b: out+=term.subs({X:1,Y:1})
    avg_cache[m]=sp.expand(out); return avg_cache[m]
def average(F):
    Pl=sp.Poly(F,*monos); out=0
    for exps,co in Pl.terms():
        m=sp.Mul(*[v**e for v,e in zip(monos,exps)])
        out+=co*avg_mono(m)
    return out
def to_t(e):
    return e.subs({U:sig/t, V:t**2, c:(t**2+t**-2)/2, sn:-(t**2-t**-2)/(2*sp.I)})   # pair angle theta_i-theta_j = -phi_m + z, so c = cos(phi_m), sn = sin(-phi_m)
NS=14
def ser_inv(poly_coeffs,N):
    # inverse of power series a0 + a1 tau + ... (list), truncated to N terms
    a=poly_coeffs+[0]*(N-len(poly_coeffs)); b=[sp.Integer(0)]*N; b[0]=1/a[0]
    for k in range(1,N):
        b[k]=sp.expand(-sum(a[j]*b[k-j] for j in range(1,k+1))/a[0])
    return b
def ser_mul(a,b,N):
    return [sp.expand(sum(a[j]*b[k-j] for j in range(k+1) if j<len(a) and k-j<len(b))) for k in range(N)]
def kernel_series(t0,sval,N):
    """Laurent series of n t^{n-1}/(t^n - s) at t = t0+tau (n odd, so t0^n = t0). Returns (shift, coeffs)."""
    B=[sp.expand(sp.expand_func(sp.binomial(n,k))*sp.Integer(t0)**(-k)) for k in range(N+2)]
    inv_t=[sp.Integer(1)/t0*(-sp.Integer(1)/t0)**k for k in range(N+2)]   # 1/(t0+tau)
    if t0==sval:   # pole: K = n B / ((t0+tau)(B-1)),  (B-1) = tau*C
        C=B[1:]; Ci=ser_inv(C,N+1); K=ser_mul(ser_mul(B,Ci,N+1),inv_t,N+1)
        return -1,[sp.expand(n*k) for k in K]
    else:          # K = n B / ((t0+tau)(B+1))
        Bp=[B[0]+1]+B[1:]; Bi=ser_inv(Bp,N+1); K=ser_mul(ser_mul(B,Bi,N+1),inv_t,N+1)
        return 0,[sp.expand(n*k) for k in K]
def laurent_at(Rnum,Rden,t0,N):
    """Laurent coefficients of Rnum/Rden (polys in t, coeffs may contain lam) at t0: returns (shift, list)."""
    num=sp.Poly(sp.expand(Rnum.subs(t,t0+tau)),tau); den=sp.Poly(sp.expand(Rden.subs(t,t0+tau)),tau)
    # strip tau^k from den
    dc=den.all_coeffs()[::-1]; k=0
    while dc[k]==0: k+=1
    dc=dc[k:]; nc=num.all_coeffs()[::-1]
    nl=0
    while nl<len(nc) and nc[nl]==0: nl+=1
    return -k, ser_mul(nc,ser_inv(dc,N),N)
def root_sum(R,sval):
    Rt=sp.together(R); Rnum,Rden=sp.fraction(Rt)
    Rnum=sp.expand(Rnum); Rden=sp.expand(Rden)
    total=0; degs=[]
    for t0 in (1,-1):
        sh,L=laurent_at(Rnum,Rden,t0,NS)
        ksh,K=kernel_series(t0,sval,NS)
        # residue: coefficient of tau^{-1} in tau^{sh+ksh} * L*K
        need=-1-(sh+ksh)
        if need>=0:
            res=sp.expand(sum(L[j]*K[need-j] for j in range(need+1)))
            total-=res
        assert -sh < NS-2, 'series too short'
    # infinity / zero degrees
    pn=sp.Poly(Rnum,t); pd=sp.Poly(Rden,t)
    deg_inf=pn.degree()-pd.degree()
    low_n=min(m[0] for m in pn.monoms()); low_d=min(m[0] for m in pd.monoms())
    pole0=low_d-low_n
    degs=(deg_inf,pole0)
    if deg_inf>=0:
        q=sp.Symbol('q'); r0=sp.limit(sp.together(Rt.subs(t,1/q)) - 0, q, 0) if deg_inf==0 else None
        if deg_inf>0:
            ser=sp.series(Rt.subs(t,1/q),q,0,1).removeO(); r0=sp.expand(ser).coeff(q,0)
        total+=n*r0
    return sp.simplify(total),degs
def offset_sum(g):
    g=sp.together(g)
    chk=sp.simplify(sp.together(g.subs({t:-t,sig:-sig},simultaneous=True)-g))
    assert chk==0,'not invariant'
    out=0; maxdeg=0
    for sval in (1,-1):
        r,degs=root_sum(g.subs(sig,sval),sval); out+=r; maxdeg=max(maxdeg,*degs)
    return sp.simplify(out/2),maxdeg
if __name__=='__main__':
    res={}
    for k in (2,3,4):
        Ek=average(Fk[k]); log('A^%d averaged'%k)
        if k==3:
            log('  A^3 average =',sp.simplify(Ek)); continue
        Ek=sp.expand(to_t(Ek))
        Pl=sp.Poly(Ek,p,s)
        for mon,co in Pl.terms():
            val,md=offset_sum(co)
            res[(k,)+mon]=sp.factor(sp.simplify(n/2*val))
            log('  A^%d p^%d s^%d :'%((k,)+mon), res[(k,)+mon],'  [valid for n >',md,']')
    pickle.dump(res,open('odd_general.pkl','wb'))
