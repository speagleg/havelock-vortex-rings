"""Symbolic proof (all N, m) that a_m(1/x)=a_m(x), and behaviour of b_m, g_m, under x->1/x.
Treat x, X=x^N, t=x^m as independent symbols; x^{N-m}=X/t, x^{N+m}=X t, x^{2N-m}=X^2/t.
Inversion x->1/x sends (x,X,t)->(1/x,1/X,1/t). An identity in the free symbols implies it for all N,m."""
import sympy as sp
x,X,t,N,m = sp.symbols('x X t N m')
c=(1+x)/(1-x)
a=(N-1)*c-(c-1)**2/2-m*(N-m)/2+2*N*c*X/(1-X)-N**2*X/(1-X)**2 \
  -N/2*(m*t+(N-m)*X/t+(N-m)*X*t+m*X**2/t)/(1-X)**2
Wm = m*t/(1-X)+N*X*t/(1-X)**2
WNm = (N-m)*(X/t)/(1-X)+N*X*(X/t)/(1-X)**2
WN = N*X/(1-X)+N*X**2/(1-X)**2   # p=N mod N: N x^N/(1-X)+N x^{2N}/(1-X)^2 ; candidate uses N*X/(1-X)^2 (same)
print('W_N forms agree:', sp.simplify(WN - N*X/(1-X)**2)==0)
b = m*(N-m)/2 + N*N*X/(1-X)**2 - N/2*(Wm+WNm)
g = N/2*(Wm-WNm)
inv={x:1/x,X:1/X,t:1/t}
for name,e in [('a',a),('b',b),('g',g)]:
    d=sp.factor(sp.cancel(e.subs(inv,simultaneous=True)-e)); s=sp.factor(sp.cancel(e.subs(inv,simultaneous=True)+e))
    print(name,'(e(1/x)-e(x)) =',d if len(str(d))<200 else '...', '| e(1/x)+e(x) zero?', s==0)
