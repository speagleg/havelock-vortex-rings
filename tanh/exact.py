# Exact (sympy) threshold polynomials for the tanh model and their algebraic data.
import sympy as sp
x=sp.symbols('x')  # x = xi = r_E^2 = tanh^2(r/2)
def polys(N,m):
    X=x**N; R=sp.Rational
    a=(N-1)*(1+x)/(1-x)-2*x**2/(1-x)**2+2*N*(1+x)*X/((1-x)*(1-X))-R(m*(N-m),2)-N**2*X/(1-X)**2 \
      -R(N,2)*((m*x**m+(N-m)*x**(N-m))/(1-X)+N*X*(x**m+x**(N-m))/(1-X)**2)
    b=R(m*(N-m),2)+N**2*X/(1-X)**2-R(N,2)*((m*x**m+(N-m)*x**(N-m))/(1-X)+N*X*(x**m+x**(N-m))/(1-X)**2)
    g=R(N,2)*((m*x**m-(N-m)*x**(N-m))/(1-X)+N*X*(x**m-x**(N-m))/(1-X)**2)
    A=sp.Poly(sp.cancel(2*(1-x)**2*(1-X)**2*a),x)
    B=sp.Poly(sp.cancel(2*(1-X)**2*b),x)
    G=sp.Poly(sp.cancel(2*(1-X)**2*g),x)
    return A,B,G
def minpoly_xistar(N):
    m=N//2; A,_,_=polys(N,m)
    c,facs=A.factor_list()
    return A,facs
if __name__=='__main__':
    for N in range(8,17):
        A,facs=minpoly_xistar(N)
        print('N=%d  A_{%d} degree %d  factorisation degrees:'%(N,N//2,A.degree()),[(f.degree(),e) for f,e in facs])
