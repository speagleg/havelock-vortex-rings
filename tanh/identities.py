# Exact identity checks for the closed forms stated in derivation.txt (sympy).
import sympy as sp
x,p,n,k=sp.symbols('x p n k',positive=True)
N_,m_=sp.symbols('N m',positive=True,integer=True)
# W_m = sum_{p>=1, p=m mod N} p x^p : compare exact power series to order 80
okW=True
for N in range(2,13):
    for m in range(0,N):
        mm=m if m else N
        target=mm*x**mm/(1-x**N)+N*x**(mm+N)/(1-x**N)**2
        ser=sum(p*x**p for p in range(1,80) if p%N==m)
        okW&= sp.expand(sp.series(target,x,0,80).removeO()-ser)==0
print('W_m closed form (series to x^79, N=2..12, all residues):',okW)
from exact import polys
def a_raw(N,m):
    X=x**N; R=sp.Rational
    W=lambda j: j*x**j/(1-X)+N*x**(j+N)/(1-X)**2
    W0=N*X/(1-X)**2; G0=X/(1-X)
    E=-R(N-1,2)+x/(1-x)-N*G0
    Vpp=2*x/(1-x)**2-2*E*(1+x)/(1-x)
    return Vpp-R(m*(N-m),2)-N*W0-R(N,2)*(W(m)+W(N-m))
c=(1+x)/(1-x)
ok=True
for N in range(3,17):
    for m in range(1,N):
        X=x**N; R=sp.Rational
        compact=(N-1)*c-(c-1)**2/2-R(m*(N-m),2)+2*N*c*X/(1-X)-N**2*X/(1-X)**2 \
            -R(N,2)*(m*x**m+(N-m)*x**(N-m)+(N-m)*x**(N+m)+m*x**(2*N-m))/(1-X)**2
        ok&= sp.cancel(compact-a_raw(N,m))==0
    if N%2==0:
        y=x**(N//2)
        even=(N-1)*c-(c-1)**2/2-sp.Rational(N*N,8)+2*N*c*y**2/(1-y**2)-N**2*y/(2*(1-y)**2)
        ok&= sp.cancel(even-a_raw(N,N//2))==0
print('compact a_m form and even-N critical form agree with block formula, N=3..16, all m:',ok)
# palindromy origin: x -> 1/x symmetry of the critical numerator
for N in range(8,17):
    A=polys(N,N//2)[0]; co=A.all_coeffs()
    print(N,'A_crit palindromic up to sign:',co==co[::-1] or co==[-v for v in co[::-1]])
