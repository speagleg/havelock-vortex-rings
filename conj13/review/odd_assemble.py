"""Assemble the odd-n level-set quartic from odd_general.pkl (symbolic in n and lam), impose
J1 = 0 at O(A^2) and the threshold relation, compare with (n-1)(n-3)(n^2-9n+38)/(128 n)."""
import pickle, sympy as sp
n,lam,p,s=sp.symbols('n lam p s')
res=pickle.load(open('odd_general.pkl','rb'))
f2=res.get((2,0,0),0)
crit=sp.expand(8*(1+lam**2)-(n+1)*(1+lam)**2)      # (1+lam^2)/(1+lam)^2 = (n+1)/8  (odd n)
print('f2 (kernel eigenvalue, A^2 coefficient) =',sp.factor(f2))
print('  f2 numerator mod threshold relation:',sp.rem(sp.expand(sp.numer(sp.together(f2))),crit,lam))
E=sum(res.get((4,i,j),0)*p**i*s**j for i in range(3) for j in range(3))
# J1 constraint: Q'(w0) p - Q0 s + Q''(w0)/4 = 0, Q = sqrt(w(1-lam w))
w=sp.Symbol('w'); Q=sp.sqrt(w*(1-lam*w)); w0=1/(1+lam)
Q0=sp.simplify(Q.subs(w,w0)); Q1=sp.simplify(sp.diff(Q,w).subs(w,w0)); Q2=sp.simplify(sp.diff(Q,w,2).subs(w,w0))
Q0=sp.simplify(sp.sqrt(sp.factor((w*(1-lam*w)).subs(w,w0))));
print('Q0,Q1,Q2 =',Q0,Q1,Q2)
s_of_p=sp.solve(sp.Eq(Q1*p-Q0*s+Q2/4,0),s)[0]
Ep=sp.expand(sp.together(E.subs(s,s_of_p)))
a2=sp.together(Ep.coeff(p,2)); a1=sp.together(Ep.coeff(p,1)); a0=sp.together(Ep.coeff(p,0))
print('curvature in the free mode-1 direction a2 =',sp.factor(a2))
fmin=sp.together(a0-a1**2/(4*a2))
c4=sp.together(fmin*64/(n**2*(1+lam)**8))
num,den=sp.fraction(sp.factor(c4))
target=(n-1)*(n-3)*(n**2-9*n+38)/(128*n)
tn,td=sp.fraction(sp.together(target))
diff=sp.rem(sp.expand(num*td-den*tn),crit,lam)
print('c4*den - formula*den  mod threshold relation =',sp.factor(diff))
dr=sp.rem(sp.expand(den),crit,lam); print('den mod relation (must be nonzero):',sp.factor(dr))
# numeric sanity at n=9,11,13 (b_n from relation)
for nv in (5,7,9,11,13,15,17):
    for r in sp.Poly(crit.subs(n,nv),lam).nroots(n=40):
        if -1<r<=1:
            print(nv,'b_n=',sp.N(r,15),' c4 =',sp.N(c4.subs({n:nv,lam:r}),25),' formula=',sp.N(target.subs(n,nv),25),
                  ' a2=',sp.N(a2.subs({n:nv,lam:r}),8))
# n = 7 separately (the reduction mod the threshold relation divides by (n-7)); exact, lam = 0:
print('n=7 exact: c4 =', sp.nsimplify(sp.simplify(c4.subs({n:7,lam:0}))), ' formula =', target.subs(n,7))
# denominator vanishing check: den = 0 mod relation iff 2 lam (n+1) + n - 7 = 0; with the relation this forces lam = 0, n = 7
print('a2 > 0 on (-1,1):', sp.factor(a2))
