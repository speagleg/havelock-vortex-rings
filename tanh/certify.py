# Rigorous certification for the tanh model, N = 2..16.
#  - exact integer polynomials A_m, B_m, D_m = A_m B_m - (1-x)^2 G_m^2 (sympy, from exact.py)
#  - root counting / isolation on (0,1) with exact rational arithmetic (sympy Sturm/Vincent-type isolation)
#  - Arb ball enclosure of xi*(N) (python-flint acb roots of the irreducible factor, 400 bits)
# Sign facts used:  a_m = A_m / (2(1-x)^2(1-x^N)^2),  b_m = B_m / (2(1-x^N)^2),
#   a_m b_m - g_m^2 = D_m / (4 (1-x)^2 (1-x^N)^4)  -- all denominators > 0 on (0,1).
# Spectral stability <=> a_m b_m > 0 for all m=1..N-1 ; energy-momentum (Lyapunov, sufficient) <=> a_m>0, a_m b_m-g_m^2>0.
import sympy as sp, flint, sys
from fractions import Fraction as Fr
from exact import polys, x
flint.ctx.prec = 400
def strip(P):  # remove factors x and (x-1)
    P=sp.Poly(P,x)
    while P.eval(0)==0: P=sp.quo(P,sp.Poly(x,x))
    while P.eval(1)==0: P=sp.quo(P,sp.Poly(x-1,x))
    return P
def roots01(P):
    """exact isolating intervals (rational) of the distinct real roots of P in (0,1)"""
    P=strip(P); S=sp.Poly(sp.quo(P,sp.gcd(P,P.diff(x))),x)
    return [(Fr(int(sp.numer(a)),int(sp.denom(a))),Fr(int(sp.numer(b)),int(sp.denom(b)))) for (a,b),k in S.intervals(inf=0,sup=1,eps=sp.Rational(1,10**15)) if b>0 and a<1], P
def sign_at(P,q): return sp.sign(P.eval(sp.Rational(q.numerator,q.denominator)))
def palin_reduce(p):  # p palindromic degree 2k -> Q(c) with p(x)=x^k Q(x+1/x)
    c=sp.symbols('c'); k=p.degree()//2
    # do it by repeated subtraction of Chebyshev-like (x+1/x)^j
    e=sp.expand(p.as_expr()/x**k); Q=0
    for j in range(k,-1,-1):
        co=sp.Poly(sp.expand(e*x**k),x).coeff_monomial(x**(k+j))
        Q+=co*c**j; e=sp.expand(e-co*(x+1/x)**j)
    assert sp.simplify(e)==0
    return sp.Poly(Q,c)
out=[]
def pr(*a):
    s=' '.join(str(v) for v in a); print(s,flush=True); out.append(s)
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 16
for N in range(2,NMAX+1):
    data={}
    for m in range(1,N//2+1):
        A,B,G=polys(N,m)
        D=A*B-sp.Poly((1-x)**2,x)*G*G
        rA,As=roots01(A); rB,Bs=roots01(B); rD,Ds=roots01(D)
        assert rB==[] and sign_at(Bs,Fr(1,2))>0, ('b_m not positive',N,m)
        data[m]=(rA,As,rD,Ds)
    mc=N//2
    if all(data[m][0]==[] for m in data):
        assert all(sign_at(data[m][1],Fr(1,2))>0 and data[m][2]==[] and sign_at(data[m][3],Fr(1,2))>0 for m in data)
        pr('N=%d: b_m>0, a_m>0, det>0 on all of (0,1) for every m  => spectrally AND energy-momentum stable for all r (certified, exact Sturm)'%N)
        continue
    # largest root of each A_m
    tops={m:(data[m][0][-1] if data[m][0] else None) for m in data}
    lo,hi=tops[mc]
    for m in data:
        if m!=mc and tops[m] is not None: assert tops[m][1]<lo, ('other mode crosses later',N,m)
        # above its largest root each A_m is positive:
        t=tops[m][1] if tops[m] else Fr(0)
        assert sign_at(data[m][1],(t+1)/2)>0
    assert len(data[mc][0])==1, 'critical mode has several roots'
    unstable_modes=[m for m in data if tops[m] is not None]
    # Arb enclosure from the irreducible factor
    A=polys(N,mc)[0]; facs=A.factor_list()[1]
    fac=[f for f,e in facs if roots01(f)[0]]
    assert len(fac)==1; f=fac[0]
    fp=flint.fmpz_poly([int(c) for c in reversed(f.all_coeffs())])
    ball=[r for r,mult in fp.complex_roots() if abs(r.imag)<1e-100 and lo<= float(r.real.mid()) <=hi]
    assert len(ball)==1; z=ball[0].real
    rad=float(z.rad())
    # independent exact refinement of the isolating interval
    S=f; I=[iv for iv,k in sp.Poly(S,x).intervals(eps=sp.Rational(1,10**45),inf=sp.Rational(lo.numerator,lo.denominator),sup=sp.Rational(hi.numerator,hi.denominator))]
    assert len(I)==1; a,b=I[0]
    assert flint.arb(sp.Rational(a).p)/sp.Rational(a).q <= z.mid() <= flint.arb(sp.Rational(b).p)/sp.Rational(b).q
    pal = f.all_coeffs()==f.all_coeffs()[::-1]
    lc=f.LC(); Q=palin_reduce(f) if pal else None
    # energy-momentum threshold
    em=max((data[m][2][-1] for m in data if data[m][2]),key=lambda iv:iv[0])
    em_mode=[m for m in data if data[m][2] and data[m][2][-1]==em]
    Dm=data[em_mode[0]][3]; Dfac=[g for g,e in Dm.factor_list()[1] if roots01(g)[0] and roots01(g)[0][-1][0]<=em[0]<=roots01(g)[0][-1][1]]
    gp=flint.fmpz_poly([int(c) for c in reversed(Dfac[0].all_coeffs())])
    emb=[r.real for r,mu in gp.complex_roots() if abs(r.imag)<1e-100 and em[0]<=float(r.real.mid())<=em[1]]
    pr('N=%d critical mode m=%d; modes with a_m<0 somewhere in (0,1): %s; b_m>0 all m; other modes cross strictly earlier: certified'%(N,mc,unstable_modes))
    pr('   xi*  = %s   (Arb radius %.1e; exact rational isolating interval width %.1e agrees)'%(z.str(55,radius=False),rad,float(b-a)))
    pr('   r_E* = %s   r* = 2 artanh sqrt(xi*) = %s'%(z.sqrt().str(30,radius=False),(2*z.sqrt().atanh()).str(30,radius=False)))
    pr('   minpoly deg %d, leading coeff %s, palindromic %s, algebraic integer(<=>unit) %s'%(f.degree(),lc,pal,abs(lc)==1))
    if pal:
        pr('   c=xi*+1/xi* = 2+4/sinh^2 r* has minpoly deg %d: %s'%(Q.degree(),Q.as_expr()))
    pr('   minpoly(xi*) =',f.as_expr())
    pr('   disc(minpoly xi*) = %s'%fp.discriminant())
    pr('   energy-momentum (Hessian PD) threshold: mode %s, xi_EM = %s  (diff from xi* = %.3e)'%(em_mode,emb[0].str(30,radius=False),float((emb[0]-z).mid())))
open('certify.out','w').write('\n'.join(out)+'\n')
