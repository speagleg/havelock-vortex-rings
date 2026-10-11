import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))   # tanh/exact.py
import sympy as sp
from fractions import Fraction as Fr
from math import comb, lcm
from exact import polys, x

def coeffs(P):  # const first, ints
    P=sp.Poly(P,x); n=P.degree()
    return [int(P.coeff_monomial(x**i)) for i in range(n+1)]
def strip(P):
    P=sp.Poly(P,x); a=b=0
    while P.eval(0)==0: P=sp.quo(P,sp.Poly(x,x)); a+=1
    while P.eval(1)==0: P=sp.quo(P,sp.Poly(x-1,x)); b+=1
    return P,a,b
def ev(c,q): 
    r=Fr(0)
    for a in reversed(c): r=r*q+a
    return r
def bern(c,u,v):
    """b_k with P(x)=sum_k b_k (x-u)^k (v-x)^(n-k)"""
    n=len(c)-1; u=Fr(u); v=Fr(v); h=v-u
    # P(u + h t) coefficients in t
    t=[Fr(0)]*(n+1)
    # expand via Taylor at u
    cc=[Fr(a) for a in c]
    for k in range(n+1):
        # k-th Taylor coeff = sum_{i>=k} c_i C(i,k) u^(i-k)
        t[k]=sum(cc[i]*comb(i,k)*u**(i-k) for i in range(k,n+1))
    # P = sum_k t_k (x-u)^k ; write in t-variable s=(x-u)/h: sum t_k h^k s^k; Bernstein on [0,1] in s:
    a=[t[k]*h**k for k in range(n+1)]
    bb=[sum(Fr(comb(k,i),comb(n,i))*a[i] for i in range(k+1)) for k in range(n+1)]
    # s^k (1-s)^(n-k) * C(n,k) = C(n,k) (x-u)^k (v-x)^(n-k) / h^n
    return [bb[k]*comb(n,k)/h**n for k in range(n+1)]
def nice_mid(u,v):
    # simple rational strictly between u and v
    u=Fr(u); v=Fr(v)
    d=1
    while True:
        lo=u*d; k=int(lo)+1
        if Fr(k,d)<v and Fr(k,d)>u:
            # choose closest to midpoint among denominators d
            mid=(u+v)/2; k=round(mid*d)
            if u<Fr(k,d)<v: return Fr(k,d)
        d*=2 if d<16 else 10
def cover(c,u,v,depth=0):
    """list of (u,v,bcoeffs) pieces + split points; all bernstein coeffs >=0 with some >0"""
    b=bern(c,u,v)
    if all(z>=0 for z in b) and any(z>0 for z in b): return [('piece',Fr(u),Fr(v),b)]
    assert depth<40, 'no cover'
    w=nice_mid(u,v)
    assert ev(c,w)>0, ('nonpositive at split',w)
    return cover(c,u,w,depth+1)+[('pt',w)]+cover(c,w,v,depth+1)
def lit(q):
    q=Fr(q)
    if q.denominator==1: return '(%d : ℝ)'%q.numerator
    return '((%d : ℝ) / %d)'%(q.numerator,q.denominator)
def L(c): return '[' + ', '.join(str(a) for a in c) + ']'

out=[]; cnt=[0]
def emit(s): out.append(s)

def pos_lemma(name,c,u,v):
    """emit theorem name (x) (hu : u < x) (hv : x < v) : 0 < peval c x"""
    if len(c)==1:
        assert c[0]>0
        emit(f"""theorem {name} (x : ℝ) (hu : {lit(u)} < x) (hv : x < {lit(v)}) :
    0 < peval {L(c)} x := by norm_num [peval]
""")
        return 1
    pieces=cover(c,u,v)
    subs=[]
    for i,p in enumerate(pieces):
        if p[0]!='piece': continue
        _,pu,pv,b=p
        K=lcm(*[z.denominator for z in b])
        ints=[int(z*K) for z in b]; n=len(b)-1
        terms=[]
        for k,ck in enumerate(ints):
            if ck==0: continue
            t=str(ck)
            if k>0: t+=' * p ^ %d'%k
            if n-k>0: t+=' * q ^ %d'%(n-k)
            terms.append(t)
        expr=' + '.join(terms)
        sub='%s_%d'%(name,len(subs))
        subs.append(sub)
        emit(f"""theorem {sub} (x : ℝ) (hu : {lit(pu)} < x) (hv : x < {lit(pv)}) :
    0 < peval {L(c)} x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < {expr} := by
    intro p q hp hq; positivity
  have key : ({K} : ℝ) * peval {L(c)} x =
      (fun p q : ℝ => {expr}) (x - {lit(pu)}) ({lit(pv)} - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < ({K} : ℝ) * peval {L(c)} x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)
""")
    pts=[p[1] for p in pieces if p[0]=='pt']
    if not pts:
        emit(f"""theorem {name} (x : ℝ) (hu : {lit(u)} < x) (hv : x < {lit(v)}) :
    0 < peval {L(c)} x := {subs[0]} x hu hv
""")
        return 1
    s=f"""theorem {name} (x : ℝ) (hu : {lit(u)} < x) (hv : x < {lit(v)}) :
    0 < peval {L(c)} x := by
"""
    ind='  '
    for j,w in enumerate(pts):
        low='hu' if j==0 else 'h%d'%(j-1)
        s+=f"{ind}rcases lt_trichotomy x {lit(w)} with h{j} | h{j} | h{j}\n"
        s+=f"{ind}· exact {subs[j]} x {low} h{j}\n"
        s+=f"{ind}· rw [h{j}]; norm_num [peval]\n"
        if j==len(pts)-1:
            s+=f"{ind}· exact {subs[j+1]} x h{j} hv\n"
        else:
            s+=f"{ind}· "
            # next rcases goes on same line
            ind2=ind+'  '
            # hack: we put the next line content after "· "
            s+='__NEXT__'
            ind=ind2
    # fix "· __NEXT__<ind>rcases" -> "· rcases"
    import re
    s=re.sub(r'__NEXT__ +', '', s)
    emit(s)
    return len(subs)


def numer(N,m,kind):
    A,B,G=polys(N,m)
    if kind=='A': return A
    if kind=='B': return B
    return A*B-sp.Poly((1-x)**2,x)*G*G
def link_lemma(name,N,m,kind,Pfull):
    """emit: block * den = x^a * (x-1)^b * peval S x   (S has the sign folded in). returns (a,b,S)"""
    S,a,b=strip(Pfull)
    # Pfull = x^a (x-1)^b S
    Sc=coeffs(S.as_expr())
    nat_rw=[]
    for (l,r) in [(f'{N} - {m}',N-m),(f'2 * {N} - {m}',2*N-m),(f'{N} + {m}',N+m),(f'{N} - {m} + {N}',2*N-m),(f'{m} + {N}',N+m),(f'{N} + {N}',2*N)]:
        nat_rw.append(f'show ({l} : ℕ) = {r} by norm_num')
    if kind=='A':
        lhs=f'aBlock {N} {m} x * (2 * (1 - x) ^ 2 * (1 - x ^ {N}) ^ 2)'; unf='aBlock cR'
    elif kind=='B':
        lhs=f'bBlock {N} {m} x * (2 * (1 - x ^ {N}) ^ 2)'; unf='bBlock W'
    else:
        lhs=f'(aBlock {N} {m} x * bBlock {N} {m} x - gBlock {N} {m} x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ {N}) ^ 4)'; unf='aBlock bBlock gBlock W cR'
    rhs=(f'x ^ {a} * ' if a else '')+f'(x - 1) ^ {b} * peval {L(Sc)} x'
    emit(f"""theorem {name} (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    {lhs} = {rhs} := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ {N} ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold {unf}
  simp only [{', '.join(nat_rw)}]
  push_cast
  field_simp
  simp only [peval]
  ring
""")
    return a,b,Sc

import mpmath as mp
mp.mp.dps=60
def link_lemma2(name,N,m,kind,Pfull):
    """like link_lemma but with (1-x)^b; requires b even"""
    S,a,b=strip(Pfull); assert b%2==0
    Sc=coeffs(S.as_expr())
    nat_rw=[]
    for (l,r) in [(f'{N} - {m}',N-m),(f'2 * {N} - {m}',2*N-m),(f'{N} + {m}',N+m),(f'{N} - {m} + {N}',2*N-m),(f'{m} + {N}',N+m),(f'{N} + {N}',2*N)]:
        nat_rw.append(f'show ({l} : ℕ) = {r} by norm_num')
    if kind=='A':
        lhs=f'aBlock {N} {m} x * (2 * (1 - x) ^ 2 * (1 - x ^ {N}) ^ 2)'; unf='aBlock cR'
    elif kind=='B':
        lhs=f'bBlock {N} {m} x * (2 * (1 - x ^ {N}) ^ 2)'; unf='bBlock W'
    else:
        lhs=f'(aBlock {N} {m} x * bBlock {N} {m} x - gBlock {N} {m} x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ {N}) ^ 4)'; unf='aBlock bBlock gBlock W cR'
    rhs=(f'x ^ {a} * ' if a else '')+f'(1 - x) ^ {b} * peval {L(Sc)} x'
    emit(f"""theorem {name} (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    {lhs} = {rhs} := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ {N} ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold {unf}
  simp only [{', '.join(nat_rw)}]
  push_cast
  field_simp
  simp only [peval]
  ring
""")
    return a,b,Sc
def den(N,kind):
    if kind=='A': return f'(mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))'
    if kind=='B': return f'(mul_pos two_pos (pow_pos hX 2))'
    return f'(mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))'
def blockexpr(N,m,kind):
    if kind=='A': return f'aBlock {N} {m} x'
    if kind=='B': return f'bBlock {N} {m} x'
    return f'aBlock {N} {m} x * bBlock {N} {m} x - gBlock {N} {m} x ^ 2'
def block_pos(N,m,kind,lo='0',lohyp='h0',posname=None,linkname=None,a=0,b=4,name=None):
    """0 < block on (lo,1) given link lemma and pos lemma for S on (lo,1)"""
    numer=(f'mul_pos (mul_pos (pow_pos h0 {a}) (pow_pos h1x {b})) ({posname} x {lohyp} h1)') if a else (f'mul_pos (pow_pos h1x {b}) ({posname} x {lohyp} h1)')
    emit(f"""theorem {name} (x : ℝ) (h0 : 0 < x) {'' if lo=='0' else f'(hlo : {lo} < x) '}(h1 : x < 1) :
    0 < {blockexpr(N,m,kind)} := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ {N} := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := {numer}
  rw [← {linkname} x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right {den(N,kind)}).mp hn
""")

def sym_rw(N,m,mp_,kinds):
    r=[]
    if 'A' in kinds: r.append(f'show aBlock {N} {m} x = aBlock {N} {mp_} x from aBlock_symm {N} {mp_} (by norm_num) x')
    if 'B' in kinds: r.append(f'show bBlock {N} {m} x = bBlock {N} {mp_} x from bBlock_symm {N} {mp_} (by norm_num) x')
    if 'D' in kinds: r.append(f'show gBlock {N} {m} x = -gBlock {N} {mp_} x from gBlock_symm {N} {mp_} (by norm_num) x, neg_sq')
    return 'rw ['+', '.join(r)+']'

def gen_small(N):
    P=f'N{N}_'
    for m in range(1,N//2+1):
        for k in 'ABD':
            a,b,S=link_lemma2(P+f'link{k}{m}',N,m,k,numer(N,m,k))
            pos_lemma(P+f'pos{k}{m}',S,0,1)
            block_pos(N,m,k,posname=P+f'pos{k}{m}',linkname=P+f'link{k}{m}',a=a,b=b,name=P+f'{k}{m}_pos')
    s=f"""/-- Note Thm 7(a) for `N = {N}`: on all of `(0,1)`, every Fourier block has `a_m > 0`,
`b_m > 0` and `a_m b_m − g_m² > 0` (spectrally stable, positive-definite non-rotational Hessian). -/
theorem stable_N{N} (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (m : ℕ) (hm1 : 1 ≤ m) (hm2 : m ≤ {N-1}) :
    0 < aBlock {N} m x ∧ 0 < bBlock {N} m x ∧
      0 < aBlock {N} m x * bBlock {N} m x - gBlock {N} m x ^ 2 := by
  interval_cases m
"""
    for m in range(1,N):
        mm=m if m<=N//2 else N-m
        t=f'⟨{P}A{mm}_pos x h0 h1, {P}B{mm}_pos x h0 h1, {P}D{mm}_pos x h0 h1⟩'
        if m!=mm: s+=f'  · {sym_rw(N,m,mm,"ABD")}\n    exact {t}\n'
        else: s+=f'  · exact {t}\n'
    emit(s)

def froot(f,guess):
    fx=lambda z: sum(c*z**i for i,c in enumerate(f))
    return mp.findroot(fx, mp.mpf(guess))

XISTAR={8:'0.0673562230040812622684374765662',9:'0.1128710632366027404396940190565',10:'0.1661142415262391180817586997849',
 11:'0.2046011768826692906355358331548',12:'0.2474177768040677091739650887199',13:'0.2798613886439066870299965646857',
 14:'0.3149389631155066233935054633811',15:'0.3424585835121702566010226386713',16:'0.3716818139060454566367287813219'}

def gen_big(N):
    P=f'N{N}_'; mc=N//2
    A=numer(N,mc,'A'); S,a0,b0=strip(A); assert a0==0 and b0==4
    facs=A.factor_list()[1]
    f=[sp.Poly(g,x) for g,e in facs if g.degree()>=10][0]
    fc=coeffs(f.as_expr())
    if fc[0]<0: fc=[-c for c in fc]
    fpoly=sp.Poly(sum(c*x**i for i,c in enumerate(fc)),x)
    Hn=sp.Poly(-sp.quo(S,fpoly),x); assert sp.rem(S,fpoly).is_zero
    Hc=coeffs(Hn.as_expr()); assert all(c>=0 for c in Hc)
    xs=froot(fc,XISTAR[N]); assert abs(xs-mp.mpf(XISTAR[N]))<1e-25
    ta=Fr(int(mp.floor(xs*10**20)),10**20); tb=ta+Fr(1,10**20)
    assert ev(fc,ta)>0 and ev(fc,tb)<0
    df=[i*fc[i] for i in range(1,len(fc))]; ndf=[-c for c in df]
    dr=[z for z in sp.Poly(sum(c*x**i for i,c in enumerate(df)),x).real_roots() if 0<z<1]
    Lm=Fr(0) if not dr else Fr(1,16)
    assert all(z<Lm for z in dr)
    # alpha0: coarse dyadic lower bound for x > xi* used for the non-crit modes
    mx=Fr(0)
    for m in range(1,mc):
        Sm,_,_=strip(numer(N,m,'A'))
        for (lo,hi),k in sp.Poly(Sm,x).intervals(inf=0,sup=1,eps=sp.Rational(1,10**12)):
            mx=max(mx,Fr(int(sp.numer(hi)),int(sp.denom(hi))))
    k=4
    while True:
        al0=Fr(int(mp.floor(xs*2**k)),2**k)
        if al0>mx and al0<ta: break
        k+=1
    nonroot={}
    for m in range(1,mc+1):
        a,b,SB=link_lemma2(P+f'linkB{m}',N,m,'B',numer(N,m,'B'))
        pos_lemma(P+f'posB{m}',SB,0,1)
        block_pos(N,m,'B',posname=P+f'posB{m}',linkname=P+f'linkB{m}',a=a,b=b,name=P+f'B{m}_pos')
        if m==mc: continue
        a,b,SA=link_lemma2(P+f'linkA{m}',N,m,'A',numer(N,m,'A'))
        rts=sp.Poly(sum(c*x**i for i,c in enumerate(SA)),x).intervals(inf=0,sup=1)
        if not rts:
            nonroot[m]=True
            pos_lemma(P+f'posA{m}',SA,0,1)
            block_pos(N,m,'A',posname=P+f'posA{m}',linkname=P+f'linkA{m}',a=a,b=b,name=P+f'A{m}_pos')
        else:
            nonroot[m]=False
            pos_lemma(P+f'posA{m}',SA,al0,1)
            block_pos(N,m,'A',lo=lit(al0),lohyp='hlo',posname=P+f'posA{m}',linkname=P+f'linkA{m}',a=a,b=b,name=P+f'A{m}_pos')
    a,b,SA=link_lemma2(P+f'linkA{mc}',N,mc,'A',A)
    emit(f"""/-- The minimal polynomial `f_{N}` of `ξ*({N})` (constant term first). -/
def f{N} : List ℤ := {L(fc)}

theorem {P}crit_factor (x : ℝ) : peval {L(SA)} x = -(peval {L(Hc)} x * peval f{N} x) := by
  simp only [f{N}, peval]; ring
""")
    pos_lemma(P+'posH',Hc,0,1)
    pos_lemma(P+'dneg',ndf,Lm,1)
    if Lm>0:
        pos_lemma(P+'fleft',fc,0,Lm)
        hl=f"""(fun x hx0 hx => by
      rcases hx.lt_or_eq with h | h
      · exact {P}fleft x hx0 h
      · rw [h]; norm_num [f{N}, peval])"""
    else:
        hl="(fun x hx0 hx => absurd hx (not_le.mpr hx0))"
    emit(f"""theorem {P}dlink (x : ℝ) : dpeval f{N} x = -peval {L(ndf)} x := by
  simp only [f{N}, dpeval, peval]; ring

/-- `f_{N}` changes sign exactly once in `(0,1)`, at a point of `({ta}, {tb})`. -/
theorem {P}root : ∃ ξs, {lit(ta)} < ξs ∧ ξs < {lit(tb)} ∧ peval f{N} ξs = 0 ∧
    ∀ x, 0 < x → x < 1 → (x < ξs → 0 < peval f{N} x) ∧ (ξs < x → peval f{N} x < 0) := by
  apply sign_change (peval f{N}) {lit(Lm)} _ _ (by norm_num) (by norm_num) (by norm_num)
    (continuous_peval _).continuousOn
  · apply peval_strictAntiOn
    intro x hu hv
    rw [{P}dlink]
    linarith [{P}dneg x hu hv]
  · exact {hl}
  · norm_num [f{N}, peval]
  · norm_num [f{N}, peval]

theorem {P}A{mc}_neg (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : 0 < peval f{N} x) :
    aBlock {N} {mc} x < 0 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ {N} := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : (1 - x) ^ 4 * peval {L(SA)} x < 0 := by
    rw [{P}crit_factor]
    have := mul_pos ({P}posH x h0 h1) hf
    nlinarith [pow_pos h1x 4]
  rw [← {P}linkA{mc} x h0 h1] at hn
  have hd := {den(N,'A')}
  by_contra hc
  push Not at hc
  nlinarith [mul_nonneg hc hd.le]

theorem {P}A{mc}_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : peval f{N} x < 0) :
    0 < aBlock {N} {mc} x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ {N} := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : 0 < (1 - x) ^ 4 * peval {L(SA)} x := by
    rw [{P}crit_factor]
    have := mul_pos ({P}posH x h0 h1) (neg_pos.mpr hf)
    have h4 := pow_pos h1x 4
    nlinarith
  rw [← {P}linkA{mc} x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right {den(N,'A')}).mp hn
""")
    s=f"""/-- Note Thm 7(b) for `N = {N}`: there is a unique `ξ* ∈ (0,1)` (the unique root of `f_{N}` in
`(0,1)`, enclosed in `({ta}, {tb})`) such that, for every `ξ ∈ (0,1)`: all `b_m > 0`;
`a_{mc} < 0` for `ξ < ξ*` (spectrally unstable, since `b_{mc} > 0`); and all `a_m > 0` for `ξ > ξ*`. -/
theorem cert_N{N} : ∃ ξs : ℝ, {lit(ta)} < ξs ∧ ξs < {lit(tb)} ∧ peval f{N} ξs = 0 ∧
    (∀ y, 0 < y → y < 1 → peval f{N} y = 0 → y = ξs) ∧
    ∀ x, 0 < x → x < 1 →
      (∀ m, 1 ≤ m → m ≤ {N-1} → 0 < bBlock {N} m x) ∧
      (x < ξs → aBlock {N} {mc} x < 0) ∧
      (ξs < x → ∀ m, 1 ≤ m → m ≤ {N-1} → 0 < aBlock {N} m x) := by
  obtain ⟨ξs, hta, htb, hf0, hs⟩ := {P}root
  refine ⟨ξs, hta, htb, hf0, ?_, ?_⟩
  · intro y hy0 hy1 hy
    rcases lt_trichotomy y ξs with h | h | h
    · linarith [(hs y hy0 hy1).1 h]
    · exact h
    · linarith [(hs y hy0 hy1).2 h]
  intro x h0 h1
  refine ⟨?_, fun hlt => {P}A{mc}_neg x h0 h1 ((hs x h0 h1).1 hlt), ?_⟩
  · intro m hm1 hm2
    interval_cases m
"""
    for m in range(1,N):
        mm=m if m<=mc else N-m
        if m!=mm: s+=f'    · {sym_rw(N,m,mm,"B")}\n      exact {P}B{mm}_pos x h0 h1\n'
        else: s+=f'    · exact {P}B{mm}_pos x h0 h1\n'
    s+=f"""  · intro hgt m hm1 hm2
    have hlo : {lit(al0)} < x := lt_trans (lt_trans (by norm_num) hta) hgt
    have hfn := (hs x h0 h1).2 hgt
    interval_cases m
"""
    for m in range(1,N):
        mm=m if m<=mc else N-m
        if mm==mc: t=f'{P}A{mc}_pos x h0 h1 hfn'
        elif nonroot[mm]: t=f'{P}A{mm}_pos x h0 h1'
        else: t=f'{P}A{mm}_pos x h0 hlo h1'
        if m!=mm: s+=f'    · {sym_rw(N,m,mm,"A")}\n      exact {t}\n'
        else: s+=f'    · exact {t}\n'
    emit(s)
    return dict(f=fc,al0=al0,Lm=Lm,ta=ta,tb=tb,xs=xs,nonroot=nonroot,Hc=Hc)

def build(Ns, path):
    out.clear()
    info={}
    for N in Ns:
        emit(f'/-! ### N = {N} -/\n')
        if N<=7: gen_small(N)
        else: info[N]=gen_big(N)
    g=open(os.path.join(HERE,'General.lean')).read().replace('end Havelock.Tanh','\n'.join(out)+'\nend Havelock.Tanh')
    g=g.replace('noncomputable section','set_option linter.unusedSimpArgs false\n\nnoncomputable section',1)
    open(path,'w').write(g)
    return info

def gen_extras(info):
    # factorisation lemmas for each N>=8
    for N,d in info.items():
        mc=N//2
        A=numer(N,mc,'A'); S,_,_=strip(A)
        fp=sp.Poly(sum(c*x**i for i,c in enumerate(d['f'])),x)
        H=sp.Poly(-sp.quo(S,fp),x)
        if N%2==0:
            n=N//2; e=1 if N%4==2 else 0
            base=sp.Poly((x**n+1)*(x+1)**e,x)
            c=sp.quo(H,base); assert sp.rem(H,base).is_zero and c.degree()==0
            c=int(c.as_expr())
            hexpr=f'{c} * (x ^ {n} + 1)'+(f' * (x + 1)' if e else '')
            # A = (1-x)^4 * S = -(1-x)^4 H f ; (1-x)^4=(x-1)^4
            emit(f"""/-- Factorisation of the critical numerator for `N = {N}` (note §6, "Structure of the numerator"):
`A_{mc} = 2(1−ξ)²(1−ξ^{N})² a_{mc} = −{c}·(ξ−1)⁴ (ξ^{n}+1){'(ξ+1)' if e else ''} f_{N}(ξ)`. -/
theorem A_crit_factor_N{N} (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock {N} {mc} x * (2 * (1 - x) ^ 2 * (1 - x ^ {N}) ^ 2)
      = -({hexpr}) * (x - 1) ^ 4 * peval f{N} x := by
  rw [N{N}_linkA{mc} x h0 h1, N{N}_crit_factor]
  simp only [peval]; ring
""")
        else:
            assert H.degree()==0; c=int(H.as_expr())
            emit(f"""/-- Factorisation of the critical numerator for `N = {N}`:
`A_{mc} = −{c}·(ξ−1)⁴ f_{N}(ξ)`. -/
theorem A_crit_factor_N{N} (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock {N} {mc} x * (2 * (1 - x) ^ 2 * (1 - x ^ {N}) ^ 2) = -{c} * (x - 1) ^ 4 * peval f{N} x := by
  rw [N{N}_linkA{mc} x h0 h1, N{N}_crit_factor]
  simp only [peval]; ring
""")
    # palindromy / degree / leading coefficient
    for N,d in info.items():
        f=d['f']
        emit(f"""/-- `f_{N}` is palindromic, of degree {len(f)-1}, with leading coefficient {f[-1]} (note Table 1). -/
theorem f{N}_data : f{N}.reverse = f{N} ∧ f{N}.length = {len(f)} ∧ f{N}.getLast? = some {f[-1]} := by
  refine ⟨?_, ?_, ?_⟩ <;> rfl
""")
    # Krein window, odd N
    for N,d in info.items():
        if N%2==0: continue
        mc=N//2
        xs=d['xs']; w=Fr(int(mp.ceil(xs*10**10)),10**10)
        assert w>d['tb'] and w<1
        A,B,G=polys(N,mc); D=A*B-sp.Poly((1-x)**2,x)*G*G
        assert D.eval(sp.Rational(w.numerator,w.denominator))<0
        W=lit(w)
        emit(f"""theorem krein_D_N{N} : aBlock {N} {mc} {W} * bBlock {N} {mc} {W} - gBlock {N} {mc} {W} ^ 2 < 0 := by
  norm_num [aBlock, bBlock, gBlock, W, cR]

/-- Odd-`N` Krein window for `N = {N}` (note Prop 10, qualitative part): just above the spectral
threshold `ξ*` there is a point where the ring is spectrally stable (all `a_m, b_m > 0`) but the
critical block is indefinite (`a b − g² < 0`); hence `ξ_EM > ξ*`. -/
theorem krein_N{N} : ∃ ξs : ℝ, peval f{N} ξs = 0 ∧ {lit(d['ta'])} < ξs ∧ ξs < {lit(d['tb'])} ∧
    (∀ y, 0 < y → y < 1 → peval f{N} y = 0 → y = ξs) ∧
    ∃ x : ℝ, ξs < x ∧ x < 1 ∧
      (∀ m, 1 ≤ m → m ≤ {N-1} → 0 < aBlock {N} m x ∧ 0 < bBlock {N} m x) ∧
      aBlock {N} {mc} x * bBlock {N} {mc} x - gBlock {N} {mc} x ^ 2 < 0 := by
  obtain ⟨ξs, hta, htb, hf0, huniq, hs⟩ := cert_N{N}
  have hw : ξs < {W} := lt_trans htb (by norm_num)
  refine ⟨ξs, hf0, hta, htb, huniq, {W}, hw, by norm_num, fun m hm1 hm2 => ?_, krein_D_N{N}⟩
  have h := hs {W} (by norm_num) (by norm_num)
  exact ⟨h.2.2 hw m hm1 hm2, h.1 m hm1 hm2⟩
""")

N8_EXTRA = r'''
/-! ### `ξ*(8)`: an algebraic unit (note Cor 8) -/

open Polynomial in
/-- The polynomial of note Cor 8, as an element of `ℤ[X]`. -/
def P8 : ℤ[X] :=
  X ^ 10 - 12 * X ^ 9 - 37 * X ^ 8 - 72 * X ^ 7 - 84 * X ^ 6 - 88 * X ^ 5 - 84 * X ^ 4
    - 72 * X ^ 3 - 37 * X ^ 2 - 12 * X + 1

theorem P8_monic : P8.Monic := by unfold P8; monicity!

theorem P8_natDegree : P8.natDegree = 10 := by unfold P8; compute_degree!

theorem aeval_P8 (ξ : ℝ) : Polynomial.aeval ξ P8 = peval f8 ξ := by
  simp [P8, f8, peval, map_ofNat] <;> ring

theorem peval_f8_inv (ξ : ℝ) (h0 : ξ ≠ 0) : peval f8 ξ⁻¹ = peval f8 ξ / ξ ^ 10 := by
  simp only [f8, peval]; field_simp; push_cast; ring

/-- Note Cor 8: any real root of `P_8` (in particular `ξ*(8)`) is a nonzero algebraic integer whose
inverse is also an algebraic integer, i.e. an algebraic unit. -/
theorem f8_root_unit (ξ : ℝ) (h : peval f8 ξ = 0) :
    ξ ≠ 0 ∧ IsIntegral ℤ ξ ∧ IsIntegral ℤ ξ⁻¹ ∧ ξ * ξ⁻¹ = 1 := by
  have h0 : ξ ≠ 0 := by
    intro h0; rw [h0] at h; norm_num [f8, peval] at h
  refine ⟨h0, ⟨P8, P8_monic, ?_⟩, ⟨P8, P8_monic, ?_⟩, mul_inv_cancel₀ h0⟩
  · rw [← Polynomial.aeval_def, aeval_P8, h]
  · rw [← Polynomial.aeval_def, aeval_P8, peval_f8_inv ξ h0, h, zero_div]

/-- Note Cor 8: `C = ξ + 1/ξ` satisfies `C⁵ − 12C⁴ − 42C³ − 24C² + 32C + 32 = 0`. -/
theorem f8_C_poly (ξ : ℝ) (h : peval f8 ξ = 0) :
    (ξ + ξ⁻¹) ^ 5 - 12 * (ξ + ξ⁻¹) ^ 4 - 42 * (ξ + ξ⁻¹) ^ 3 - 24 * (ξ + ξ⁻¹) ^ 2
      + 32 * (ξ + ξ⁻¹) + 32 = 0 := by
  have h0 : ξ ≠ 0 := (f8_root_unit ξ h).1
  have key : (ξ + ξ⁻¹) ^ 5 - 12 * (ξ + ξ⁻¹) ^ 4 - 42 * (ξ + ξ⁻¹) ^ 3 - 24 * (ξ + ξ⁻¹) ^ 2
      + 32 * (ξ + ξ⁻¹) + 32 = peval f8 ξ / ξ ^ 5 := by
    simp only [f8, peval]; field_simp; push_cast; ring
  rw [key, h, zero_div]

/-- Note Cor 8 / Thm 7 for `N = 8`, packaged: the spectral threshold `ξ*(8)` is an algebraic unit,
root of the monic palindromic degree-10 `P_8`, and `ξ* + 1/ξ*` is a root of the quintic. -/
theorem xi8_unit : ∃ ξs : ℝ, ((3367811150204063113 : ℝ) / 50000000000000000000) < ξs ∧
    ξs < ((6735622300408126227 : ℝ) / 100000000000000000000) ∧
    Polynomial.aeval ξs P8 = 0 ∧ IsIntegral ℤ ξs ∧ IsIntegral ℤ ξs⁻¹ ∧
    (ξs + ξs⁻¹) ^ 5 - 12 * (ξs + ξs⁻¹) ^ 4 - 42 * (ξs + ξs⁻¹) ^ 3 - 24 * (ξs + ξs⁻¹) ^ 2
      + 32 * (ξs + ξs⁻¹) + 32 = 0 ∧
    ∀ x, 0 < x → x < 1 →
      (x < ξs → aBlock 8 4 x < 0) ∧ (ξs < x → ∀ m, 1 ≤ m → m ≤ 7 → 0 < aBlock 8 m x) := by
  obtain ⟨ξs, hta, htb, hf0, _, hs⟩ := cert_N8
  obtain ⟨_, hi, hi', _⟩ := f8_root_unit ξs hf0
  exact ⟨ξs, hta, htb, by rw [aeval_P8, hf0], hi, hi', f8_C_poly ξs hf0,
    fun x h0 h1 => ⟨(hs x h0 h1).2.1, (hs x h0 h1).2.2⟩⟩

/-! ### Relation to the hyperbolic radius -/

/-- With `ξ = tanh²(r/2)` (`r > 0`): `c = (1+ξ)/(1−ξ) = cosh r` and `ξ + 1/ξ = 2 + 4/sinh² r`
(note §3.1 and §5, "The natural invariant"). -/
theorem tanh_relations (r : ℝ) (hr : 0 < r) :
    cR (Real.tanh (r / 2) ^ 2) = Real.cosh r ∧
      Real.tanh (r / 2) ^ 2 + (Real.tanh (r / 2) ^ 2)⁻¹ = 2 + 4 / Real.sinh r ^ 2 := by
  have hs : 0 < Real.sinh (r / 2) := Real.sinh_pos_iff.mpr (by linarith)
  have hc : 0 < Real.cosh (r / 2) := Real.cosh_pos _
  have hid : Real.cosh (r / 2) ^ 2 - Real.sinh (r / 2) ^ 2 = 1 := Real.cosh_sq_sub_sinh_sq _
  have e1 : Real.cosh r = Real.cosh (r / 2) ^ 2 + Real.sinh (r / 2) ^ 2 := by
    rw [← Real.cosh_two_mul]; ring_nf
  have e2 : Real.sinh r = 2 * Real.sinh (r / 2) * Real.cosh (r / 2) := by
    rw [← Real.sinh_two_mul]; ring_nf
  rw [Real.tanh_eq_sinh_div_cosh, e1, e2]
  have hsc : Real.cosh (r / 2) ^ 2 - Real.sinh (r / 2) ^ 2 ≠ 0 := by rw [hid]; norm_num
  constructor
  · unfold cR
    rw [div_pow]
    field_simp
    nlinarith [hid]
  · rw [div_pow]
    field_simp
    nlinarith [hid]
'''

def build2(Ns, path, extras=True):
    out.clear(); info={}
    for N in Ns:
        emit(f'/-! ### N = {N} -/\n')
        if N<=7: gen_small(N)
        else: info[N]=gen_big(N)
    if extras:
        emit('/-! ### Factorisations, palindromy, Krein window -/\n')
        gen_extras(info)
        if 8 in info: emit(N8_EXTRA)
    g=open(os.path.join(HERE,'General.lean')).read().replace('end Havelock.Tanh','\n'.join(out)+'\nend Havelock.Tanh')
    g=g.replace('noncomputable section','set_option linter.unusedSimpArgs false\n\nnoncomputable section',1)
    open(path,'w').write(g)
    return info


if __name__=='__main__':
    # Regenerate the single-file Lean certificate:  python gen.py [OUT]
    # OUT defaults to ./TanhAlgebra.generated.lean; nothing under lean/ is written. OUT is header.txt, then
    # General.lean with the N = 2..16 certificates and the extras spliced in, then footer.txt (the
    # #print axioms audit). It is the pre-split TanhAlgebra.lean (commit 73c86a7) byte for byte;
    # `python lean/scripts/split_tanh.py OUT <leandir>` maps it onto the committed module layout
    # (lean/Havelock/TanhAlgebra.lean + lean/Havelock/Tanh/*.lean). No other inputs are needed (the old
    # info.pkl was only a by-product of build2 and is recomputed in memory).
    out_path = sys.argv[1] if len(sys.argv) > 1 else 'TanhAlgebra.generated.lean'
    tmp = out_path + '.body'
    build2(list(range(2, 17)), tmp)
    with open(out_path, 'w') as f:
        f.write(open(os.path.join(HERE, 'header.txt')).read() + open(tmp).read()
                + open(os.path.join(HERE, 'footer.txt')).read())
    os.remove(tmp)
    print('wrote', out_path)
