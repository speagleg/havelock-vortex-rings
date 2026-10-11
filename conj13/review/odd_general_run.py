"""Driver for odd_general.py: term-by-term offset sums.
E_k = sum over Taylor monomials x^i y^j z^k (coef[(i,j,k)](c,sn,lam)) times the lattice average of the
substituted Fourier monomials; each piece (lambda-power) x (p,s-monomial) is summed over offsets separately
by the residue formula in odd_general.py, then everything is recombined (sums are linear)."""
import pickle, time
exec(open('odd_general.py').read().split("if __name__=='__main__':")[0])
def split_lambda(R):
    R=sp.cancel(sp.together(R)); num,den=sp.fraction(R)
    cst,fl=sp.factor_list(den)
    dl=sp.Integer(1); dt=sp.Integer(1)*cst
    for f,e in fl:
        if f.has(lam):
            assert not f.has(t) and not f.has(sig), ('mixed denominator factor',f)
            dl*=f**e
        else: dt*=f**e
    Pn=sp.Poly(sp.expand(num),lam)
    return [(lam**k[0]/dl, co/dt) for k,co in Pn.terms()]
def full_offset_sum(R):
    out=0; md=0
    for lf,tr in split_lambda(R):
        v,m=offset_sum(tr); out+=lf*v; md=max(md,m)
    return out,md
xvars=(xa_i,xa_j,xb_i,xb_j,zz)
res={}; MD=0
for K in (2,4):
    tot={}
    # Fpoly coefficient of A^K, organised by Taylor monomial
    for (i,j,k),co in coef.items():
        piece=sp.expand(((A*xa_i+A**2*xb_i)**i*(A*xa_j+A**2*xb_j)**j*(A**2*zz)**k)).coeff(A,K)
        if piece==0: continue
        av=0
        for term in sp.Add.make_args(piece):
            cf,mono=term.as_coeff_Mul()
            av+=cf*avg_mono(mono)
        av=sp.expand(av)
        if av==0: continue
        Pl=sp.Poly(av,p,s)
        for mon,ac in Pl.terms():
            R=to_t(co*ac)
            val,md=full_offset_sum(R); MD=max(MD,md)
            tot[mon]=tot.get(mon,0)+val
        log('  A^%d term %s done'%(K,(i,j,k)))
    for mon,v in tot.items():
        res[(K,)+mon]=sp.factor(sp.cancel(n/2*v))
        log('A^%d p^%d s^%d :'%((K,)+mon), res[(K,)+mon])
log('max Laurent degree at 0/inf over all pieces:',MD,' (residue formula valid for n >',MD,')')
pickle.dump(res,open('odd_general.pkl','wb'))
log('done')
