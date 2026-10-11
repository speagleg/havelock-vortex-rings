"""Even n = 2M: independent symbolic check of c4 on Fix(D_M) cap J3^{-1}(mu) (two staggered M-gons).
Pair sum via prod_{zeta^M=-1}(x - y zeta) = x^M + y^M and prod_{k=1}^{M-1}(1-eta^k u)... (checked below for M=2..6).
H = -[M(M-1)/2 (log u+ + log u-) + 2M log(u+^{M/2}+u-^{M/2})] + (n-1) M (log s+ + log s-),  s = 1+lam u.
Constraint: g(u+) + g(u-) = 2 g(1), g(u) = u/(1+lam u); parametrise g(u+-) = g0 +- tau.
eps (unit Cartesian e = n^{-1/2}(-1)^k radial) = sqrt(n) * drho;  tau = 2/(1+lam)^2 * drho + O(drho^2).
"""
import sympy as sp
M,lam,tau,x,y,u=sp.symbols('M lam tau x y u',positive=True)
# product formula checks
for m in range(2,7):
    # exact: prod over roots of T^m+1 of (x - y T) = y^m * (x/y)^m + y^m ; check via resultant
    T=sp.Symbol('T')
    assert sp.expand(sp.resultant(T**m+1, x-y*T, T)*(-1)**0 - (x**m+y**m))==0 or sp.expand(sp.resultant(T**m+1,x-y*T,T)+(x**m+y**m))==0
    q=sp.cancel((T**m-1)/(T-1))   # roots eta^k, k=1..m-1
    r=sp.resultant(q,1-u*T,T)
    assert sp.expand(r-sp.cancel((1-u**m)/(1-u)))==0 or sp.expand(r+sp.cancel((1-u**m)/(1-u)))==0
print('product formulas OK for M=2..6')
# direct check of the pair-sum formula numerically for M=3, random radii
import mpmath as mp
for Mv in (2,3,4,5):
    rp,rm=mp.mpf('1.13'),mp.mpf('0.87'); pts=[rp*mp.exp(2j*mp.pi*k/Mv) for k in range(Mv)]+[rm*mp.exp(1j*mp.pi*(2*k+1)/Mv) for k in range(Mv)]
    S=mp.fsum(mp.log(abs(pts[i]-pts[j])**2) for i in range(2*Mv) for j in range(i+1,2*Mv))
    C=mp.fsum(mp.log(abs(mp.exp(2j*mp.pi*i/Mv)-mp.exp(2j*mp.pi*j/Mv))**2) for i in range(Mv) for j in range(i+1,Mv))
    F=Mv*(Mv-1)/2*(mp.log(rp**2)+mp.log(rm**2))+2*C+2*Mv*mp.log(rp**Mv+rm**Mv)
    print(' pair-sum formula M=%d diff'%Mv, mp.nstr(S-F,5))
g=lambda w: w/(1+lam*w); ginv=lambda G: G/(1-lam*G)
g0=g(1)
up=ginv(g0+tau); um=ginv(g0-tau)
n=2*M
Hf=-(M*(M-1)/2*(sp.log(up)+sp.log(um))+2*M*sp.log(sp.exp(M/2*sp.log(up))+sp.exp(M/2*sp.log(um))))+(n-1)*M*(sp.log(1+lam*up)+sp.log(1+lam*um))
ser=sp.series(Hf,tau,0,5).removeO()
h2=sp.factor(sp.simplify(ser.coeff(tau,2))); h4=sp.factor(sp.simplify(ser.coeff(tau,4)))
print('h2 =',h2); print('h4 =',h4)
# impose b_n: (1+lam^2)/(1+lam)^2 = n^2/(8(n-1))
crit=sp.expand((1+lam**2)*8*(n-1)-n**2*(1+lam)**2)
print('h2 numerator rem crit:', sp.simplify(sp.rem(sp.numer(sp.together(h2)),crit,lam)))
c4=sp.together(h4*(2/(1+lam)**2)**4/n**2)
num,den=sp.fraction(c4)
# reduce modulo crit (quadratic in lam): c4 = (a + b lam)/(c + d lam) -> rationalise
r_num=sp.rem(sp.expand(num),crit,lam); r_den=sp.rem(sp.expand(den),crit,lam)
mm=sp.Symbol('m'); target=(n-2)*((n-2)**4+40*(n-2)**2+48*(n-2)+16)/(192*n*(n-1))
diff=sp.rem(sp.expand(r_num*sp.denom(sp.together(target))-r_den*sp.numer(sp.together(target))),crit,lam)
print('c4 - formula (mod b_n relation) =', sp.simplify(diff))
