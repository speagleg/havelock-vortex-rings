# Rigorous certification for the tanh model, N = 2..16  (v2, 2026-10-10).
# Changes from certify.py (second-reader must-fix list, LOG "second read of the C5 tanh result"):
#   (1) no floating-point root selection anywhere: Arb real roots are selected by the EXACT test
#       imag == 0 (arb_fmpz_poly_complex_roots returns real roots with exactly zero imaginary part)
#       plus exact rational containment of the whole ball in the rational isolating interval;
#   (2) containment, not midpoint: [mid-rad, mid+rad] (converted exactly to rationals) must lie inside
#       the rational isolating interval [a,b] of width < 1e-45;
#   (3) explicit exact sign change of A_crit at xi*: A_crit(a) < 0 < A_crit(b) at the rational endpoints,
#       with xi* the only root of A_crit in (0,1); hence a_crit < 0 (unstable, since b_crit > 0) on (0,xi*)
#       and every a_m > 0 (all a_m b_m > 0, spectrally stable) on (xi*,1);
#   (4) energy-momentum threshold: interval disjointness and the EM vs spectral ordering are asserted exactly.
# Everything else is as in certify.py:
#  - exact integer polynomials A_m, B_m, D_m = A_m B_m - (1-x)^2 G_m^2 (sympy, from exact.py)
#  - root counting / isolation on (0,1) with exact rational arithmetic (sympy real-root isolation)
#  - Arb ball enclosure of xi*(N) (python-flint acb roots of the irreducible factor, 400 bits)
# Sign facts used:  a_m = A_m / (2(1-x)^2(1-x^N)^2),  b_m = B_m / (2(1-x^N)^2),
#   a_m b_m - g_m^2 = D_m / (4 (1-x)^2 (1-x^N)^4)  -- all denominators > 0 on (0,1).
# Spectral stability <=> a_m b_m > 0 for all m=1..N-1 ; energy-momentum (sufficient) <=> a_m>0, a_m b_m-g_m^2>0.
import sympy as sp, flint, sys
from fractions import Fraction as Fr
from exact import polys, x
flint.ctx.prec = 400
def strip(P):  # remove factors x and (x-1)
    P=sp.Poly(P,x)
    while P.eval(0)==0: P=sp.quo(P,sp.Poly(x,x))
    while P.eval(1)==0: P=sp.quo(P,sp.Poly(x-1,x))
    return P
def Q(q): return Fr(int(sp.numer(q)),int(sp.denom(q)))
def S(q): return sp.Rational(q.numerator,q.denominator)
def roots01(P):
    """exact isolating intervals (rational) of the distinct real roots of P in (0,1)"""
    P=strip(P); Sq=sp.Poly(sp.quo(P,sp.gcd(P,P.diff(x))),x)
    return [(Q(a),Q(b)) for (a,b),k in Sq.intervals(inf=0,sup=1,eps=sp.Rational(1,10**15)) if b>0 and a<1], P
def sign_at(P,q): return sp.sign(sp.Poly(P,x).eval(S(q)))
def arb_exact(z):
    """exact rational endpoints [mid-rad, mid+rad] of a real arb ball"""
    m,e=z.mid().man_exp(); r,f=z.rad().man_exp()
    mid=Fr(int(m))*(Fr(2)**int(e)); rad=Fr(int(r))*(Fr(2)**int(f))
    return mid-rad, mid+rad
def real_root_in(fz, lo, hi):
    """the unique Arb real root of integer poly fz whose WHOLE ball lies in [lo,hi] (exact tests only)"""
    fp=flint.fmpz_poly([int(c) for c in reversed(sp.Poly(fz,x).all_coeffs())])
    sel=[]
    for r,mult in fp.complex_roots():
        if not r.imag.is_zero(): continue          # exact: real roots are returned with imag exactly 0
        L,U=arb_exact(r.real)
        if lo<=L and U<=hi: sel.append(r.real)
        else: assert U<lo or L>hi, 'Arb ball straddles an isolating-interval endpoint'
    assert len(sel)==1, ('expected exactly one Arb real root in interval', len(sel))
    return sel[0], fp
def refine(f, lo, hi, eps=sp.Rational(1,10**45)):
    I=[iv for iv,k in sp.Poly(f,x).intervals(eps=eps,inf=S(lo),sup=S(hi))]
    assert len(I)==1; return Q(I[0][0]),Q(I[0][1])
def palin_reduce(p):  # p palindromic degree 2k -> Q(c) with p(x)=x^k Q(x+1/x)
    c=sp.symbols('c'); k=p.degree()//2
    e=sp.expand(p.as_expr()/x**k); Qc=0
    for j in range(k,-1,-1):
        co=sp.Poly(sp.expand(e*x**k),x).coeff_monomial(x**(k+j))
        Qc+=co*c**j; e=sp.expand(e-co*(x+1/x)**j)
    assert sp.simplify(e)==0
    return sp.Poly(Qc,c)
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
        data[m]=(rA,As,rD,Ds,A)
    mc=N//2
    if all(data[m][0]==[] for m in data):
        assert all(sign_at(data[m][1],Fr(1,2))>0 and data[m][2]==[] and sign_at(data[m][3],Fr(1,2))>0 for m in data)
        pr('N=%d: b_m>0, a_m>0, det>0 on all of (0,1) for every m  => spectrally AND energy-momentum stable for all r (certified, exact Sturm)'%N)
        continue
    tops={m:(data[m][0][-1] if data[m][0] else None) for m in data}
    lo,hi=tops[mc]
    for m in data:
        if m!=mc and tops[m] is not None: assert tops[m][1]<lo, ('other mode crosses later',N,m)
        t=tops[m][1] if tops[m] else Fr(0)
        assert sign_at(data[m][1],(t+1)/2)>0
    assert len(data[mc][0])==1, 'critical mode has several roots'
    unstable_modes=[m for m in data if tops[m] is not None]
    # irreducible factor carrying the root (exact selection by isolating intervals)
    Afull=data[mc][4]; facs=Afull.factor_list()[1]
    fac=[f for f,e in facs if roots01(f)[0]]
    assert len(fac)==1; f=fac[0]
    assert sp.Poly(f,x).is_irreducible
    # exact rational isolating interval of width < 1e-45 inside (lo,hi)
    a,b=refine(f,lo,hi)
    assert lo<=a<b<=hi and b-a<Fr(1,10**45)
    # (3) exact sign change of A_crit (whole numerator; denominators are > 0 on (0,1))
    sa,sb=sign_at(Afull,a),sign_at(Afull,b)
    assert sa<0 and sb>0, ('A_crit sign change', N, sa, sb)
    # (1)+(2) Arb ball, selected and checked by exact containment in [a,b]
    z,fp=real_root_in(f,a,b)
    L,U=arb_exact(z); rad=float(z.rad())
    assert a<=L and U<=b
    pal = f.all_coeffs()==f.all_coeffs()[::-1]
    lc=f.LC(); Qc=palin_reduce(f) if pal else None
    # energy-momentum threshold: largest root over m of D_m, intervals refined until disjoint from the rest
    cands=[(data[m][2][-1],m) for m in data if data[m][2]]
    em,emm=max(cands,key=lambda t:t[0][0])
    em_mode=[m for iv,m in cands if iv==em]
    for iv,m in cands:
        if iv!=em: assert iv[1]<=em[1]  # exact
    Dm=data[em_mode[0]][3]
    Dfac=[g for g,e in Dm.factor_list()[1] if any(em[0]<=iv[0] and iv[1]<=em[1] or (iv[0]<=em[1] and em[0]<=iv[1]) for iv in roots01(g)[0])]
    assert len(Dfac)==1
    ea,eb=refine(Dfac[0],em[0],em[1])
    for iv,m in cands:  # every other D_m's top root is certainly below the EM root
        if m not in em_mode: assert iv[1]<ea, ('EM root not separated from mode',N,m)
    emz,_=real_root_in(Dfac[0],ea,eb)
    if N%2==0:
        assert sp.Poly(Dfac[0],x)==sp.Poly(f,x) or sp.Poly(Dfac[0],x)==-sp.Poly(f,x)
        rel='equal to xi* (same irreducible factor; g_crit = 0 for even N)'
    else:
        assert ea>b, 'EM threshold not certainly above xi*'
        rel='strictly above xi*: rational intervals disjoint, [%.3e, %.3e] above'%(float(ea-b),float(eb-a))
    pr('N=%d critical mode m=%d; modes with a_m<0 somewhere in (0,1): %s; b_m>0 all m; other modes cross strictly earlier: certified'%(N,mc,unstable_modes))
    pr('   xi*  = %s   (Arb radius %.1e; Arb ball contained in exact rational isolating interval of width %.1e: certified)'%(z.str(55,radius=False),rad,float(b-a)))
    pr('   A_crit sign change at xi*: A_crit(a)<0<A_crit(b) exactly, unique root in (0,1) => unstable on (0,xi*), stable on (xi*,1): certified')
    pr('   r_E* = %s   r* = 2 artanh sqrt(xi*) = %s'%(z.sqrt().str(30,radius=False),(2*z.sqrt().atanh()).str(30,radius=False)))
    pr('   minpoly deg %d, leading coeff %s, palindromic %s, algebraic integer(<=>unit) %s'%(f.degree(),lc,pal,abs(lc)==1))
    if pal:
        pr('   c=xi*+1/xi* = 2+4/sinh^2 r* has minpoly deg %d: %s'%(Qc.degree(),Qc.as_expr()))
    pr('   minpoly(xi*) =',f.as_expr())
    pr('   disc(minpoly xi*) = %s'%fp.discriminant())
    pr('   energy-momentum (Hessian PD) threshold: mode %s, xi_EM = %s  (diff from xi* = %.3e); %s'%(em_mode,emz.str(30,radius=False),float((emz-z).mid()),rel))
open('certify_v2.out','w').write('\n'.join(out)+'\n')
