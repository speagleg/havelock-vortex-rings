"""Fig 5b (N=6): min and max eigenvalue x 4pi at a few K, for the tanh law (her eq 17) and the sinh(D) law (Fig 5a fit)."""
import sympy as sp, mpmath as mp
import boatto_hessian as bh
mp.mp.dps = 20
bh.FORMS[('sinhD',-1)] = bh.c*(sp.log(bh.d-1)+sp.log(bh.d+1)); bh.FORMS[('sinhD',1)] = bh.FORMS[('boatto',1)]
for K in ['-0.75','-0.5','-0.25','0.000001','0.25','0.5','0.75']:
    out=[]
    for model in ['boatto','sinhD']:
        mn,mx,_,_,_ = bh.min_nontrivial(6, mp.mpf(K), model)
        out.append('%s: min %7.3f max %7.3f' % (model, 4*mp.pi*mn, 4*mp.pi*mx))
    print('K=%s  ' % K + ' | '.join(out))
