"""P2: direct evaluation of Montaldi-Tokieda's beta, gamma, delta (arXiv:1210.5662v1, p.23) in THEIR normalisation.
Hhat = (H - omega0 J3)/(4 pi)  (MT p.8 with kappa = 1, r0 = 1);  alpha = alpha_r^(l) (MT (2.5), l=(n-1)/2, Euclidean
length sqrt(n/2));  v1 = sigma alpha_r^(1) - sigmatilde beta_theta^(1) (MT (2.8)).
beta = d2Hhat(v1,v1), delta = d4Hhat(alpha^4), gamma = d3Hhat(alpha,alpha,v1) = d^3Hhat/dx^2du.
Result: printed beta,gamma,delta at n=7 (exact), 9, 11 are reproduced; n=5 printed values are not.
(delta - 3 gamma^2/beta)/delta equals our linear-slice c4 / q4 (see p2_ratio_test.out)."""
import mpmath as mp, sys
mp.mp.dps=60
def run(n,lam):
    l=(n-1)//2; phi=[2*mp.pi*k/n for k in range(n)]
    # Hhat = (1/4pi)[ -sum log(|zi-zj|^2/(si sj)) ] - omega*sum |z|^2/s ; r0=1
    def H(z):
        s=[1+lam*abs(w)**2 for w in z]; v=0
        for i in range(n):
            for j in range(i+1,n): v-=mp.log(abs(z[i]-z[j])**2/(s[i]*s[j]))
        return v
    def J3(z): return mp.fsum(abs(w)**2/(1+lam*abs(w)**2) for w in z)
    om=-(n-1)*(1-lam**2)/2
    F=lambda z:(H(z)-om*J3(z))/(4*mp.pi)
    z0=[mp.expj(p) for p in phi]
    sig=1+lam; sigt=1-lam
    al=[mp.cos(l*p)*mp.expj(p) for p in phi]              # alpha_r^(l), radial unit vectors
    ar1=[mp.cos(p)*mp.expj(p) for p in phi]               # alpha_r^(1)
    bt1=[-mp.sin(p)*1j*mp.expj(p) for p in phi]           # beta_theta^(1) (tangential unit)
    br1=[-mp.sin(p)*mp.expj(p) for p in phi]
    at1=[mp.cos(p)*1j*mp.expj(p) for p in phi]
    v1=[sig*a-sigt*b for a,b in zip(ar1,bt1)]; v2=[sig*a+sigt*b for a,b in zip(br1,at1)]
    def f(*cs):
        def g(*ts):
            return F([z0[k]+sum(t*c[k] for t,c in zip(ts,cs)) for k in range(n)])
        return g
    h=mp.mpf('1e-12')
    # beta: second derivative along v1
    g=f(v1); beta=mp.diff(lambda t:g(t),0,2)
    g=f(al); delta=mp.diff(lambda t:g(t),0,4)
    out=[]
    for v in (v1,v2,[-x for x in v1]):
        g2=f(al,v); out.append(mp.diff(lambda a,b:g2(a,b),(0,0),(2,1)))
    return beta,delta,out
for n,lam in ((7,mp.mpf(0)),):
    b,d,c=run(n,lam)
    print(n,'beta*pi',mp.nstr(b*mp.pi,15),'delta*pi',mp.nstr(d*mp.pi,15),'Hxxu*pi',[mp.nstr(x*mp.pi,15) for x in c])
for n in (9,11):
    a=n-7; bb=2*(n+1); r=[(-bb+s*mp.sqrt(bb*bb-4*a*a))/(2*a) for s in (1,-1)]; lam=[x for x in r if abs(x)<1][0]
    b,d,c=run(n,lam); print(n,'beta',mp.nstr(b,8),'delta',mp.nstr(d,8),'Hxxu',[mp.nstr(x,8) for x in c], 'd-3g2/b', mp.nstr(d-3*c[0]**2/b,8))
n=5; a=n-7; bb=2*(n+1); r=[(-bb+s*mp.sqrt(bb*bb-4*a*a))/(2*a) for s in (1,-1)]; lam=[x for x in r if abs(x)<1][0]
b,d,c=run(n,lam); print(n,'lam',lam,'beta',mp.nstr(b,8),'delta',mp.nstr(d,8),'Hxxu',mp.nstr(c[0],8),'(d-3g2/b)/d',mp.nstr((d-3*c[0]**2/b)/d,8))
