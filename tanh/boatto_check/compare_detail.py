"""Pointwise comparison of digitized Fig 5a with candidate models (values x 4 pi; negative values x2 per the
observed kink artifact).  K<0: c[ln(d-1)+beta ln(d+1)];  K>0: c[ln(1-d)+gamma ln(1+d)]."""
import sys, sympy as sp, mpmath as mp, numpy as np
import boatto_hessian as bh
mp.mp.dps = 25
dig = np.load('fig5a_digitized.npy', allow_pickle=True).item()
def setmodel(beta, gamma):
    bh.FORMS[('x',-1)] = bh.c*(sp.log(bh.d-1) + sp.nsimplify(beta)*sp.log(bh.d+1))
    bh.FORMS[('x',+1)] = bh.c*(sp.log(1-bh.d) + sp.nsimplify(gamma)*sp.log(1+bh.d))
    bh._cache.pop(('x',-1),None); bh._cache.pop(('x',+1),None)
def val(N,K):
    return float(bh.min_nontrivial(N, mp.mpf(K), 'x')[0])*4*np.pi
def crossing(N, lo, hi):
    lo, hi = mp.mpf(lo), mp.mpf(hi); flo = val(N, lo)
    for _ in range(45):
        mid=(lo+hi)/2
        if val(N,mid)*flo>0: lo=mid
        else: hi=mid
    return float((lo+hi)/2)
beta, gamma = sys.argv[1], sys.argv[2]
setmodel(beta, gamma)
print('model: K<0 beta=%s, K>0 gamma=%s' % (beta, gamma))
for N in [5,6,7,8,9]:
    rows=[]
    for Kt in [-0.55,-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4,0.5,0.55]:
        near=[v for K,v in dig[N] if abs(K-Kt)<0.003]
        if not near: continue
        m = val(N,Kt); m2 = 2*m if m<0 else m
        rows.append('K=%+.2f fig=%6.3f model=%6.3f model(neg x2)=%6.3f' % (Kt, np.median(near), m, m2))
    print('N=%d'%N); print('   '+'\n   '.join(rows))
for N,(lo,hi) in {5:(0.3,0.6),6:(0.1,0.3),8:(-0.6,-0.1),9:(-0.6,-0.1)}.items():
    try: print('N=%d model crossing %.6f' % (N, crossing(N,lo,hi)))
    except Exception as e: print('N=%d no crossing in [%s,%s]'%(N,lo,hi))
