import os
import re, sys, sympy as sp, mpmath as mp
mp.mp.dps=40
x=sp.symbols('x')
claimed={8:'0.06735622300408126226',9:'0.11287106323660274043',10:'0.16611424152623911808',11:'0.20460117688266929063',12:'0.24741777680406770917',13:'0.27986138864390668702',14:'0.31493896311550662339',15:'0.34245858351217025660',16:'0.37168181390604545663'}
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tanh', 'irreducibility', 'polys.gp')
for line in open(path):
    m=re.match(r'f\[(\d+)\] = Pol\(\[(.*)\]\);',line)
    if not m: continue
    N=int(m.group(1)); c=[int(t) for t in m.group(2).split(',')]
    f=sp.Poly(c,x)
    pal = c==c[::-1]
    irr = len(sp.factor_list(f.as_expr())[1])==1 and sp.factor_list(f.as_expr())[1][0][1]==1
    rts=[r for r in sp.Poly(f).nroots(n=30) if abs(sp.im(r))<1e-25 and 0<sp.re(r)<1]
    rts=[sp.re(r) for r in rts]
    match = any(abs(r-sp.Float(claimed[N],30))<1e-19 for r in rts)
    # Salem test on 1/xi*: all other roots on unit circle?
    allr=sp.Poly(f).nroots(n=30)
    off=[r for r in allr if abs(abs(r)-1)>1e-15]
    print(f"N={N:2d} deg={f.degree():2d} lead={c[0]:3d} palindromic={pal} irreducible(Q)={irr} unit={abs(c[0])==1} roots in (0,1)={len(rts)} match_table={match}  #roots off unit circle={len(off)}")
