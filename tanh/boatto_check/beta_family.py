"""Family H = c sum [ln(d-1) + beta ln(d+1)] (K<0).  beta=-1: tanh (Boatto eq 17); beta=0: sinh(D/2) model.
Find the N=8,9 zero crossings vs beta, and compare N=6..9 values (x4pi, negatives doubled) with digitized Fig 5a."""
import sympy as sp, mpmath as mp, numpy as np
import boatto_hessian as bh
mp.mp.dps = 25
dig = np.load('fig5a_digitized.npy', allow_pickle=True).item()
def setbeta(beta):
    f = bh.c*(sp.log(bh.d-1) + sp.Rational(beta).limit_denominator(1000)*sp.log(bh.d+1))
    bh.FORMS[('beta',-1)] = f; bh._cache.pop(('beta',-1), None)
def minev(N,K): return bh.min_nontrivial(N, K, 'beta')[0]
def crossing(N, lo, hi):
    flo = minev(N, lo)
    for _ in range(40):
        mid=(lo+hi)/2
        if minev(N,mid)*flo>0: lo=mid
        else: hi=mid
    return (lo+hi)/2
import sys
for beta in (sys.argv[1:] or ['-1','-0.75','-0.5','-0.25','0','0.5','0.75','1','1.25']):
    setbeta(beta)
    c8 = crossing(8, mp.mpf('-0.6'), mp.mpf('-0.01')); c9 = crossing(9, mp.mpf('-0.8'), mp.mpf('-0.01'))
    res=[]
    for N in [6,7,8,9]:
        for K,v in dig[N]:
            if -0.58 < K < -0.02 and int(round(K*1000))%50==0:
                t = float(minev(N, mp.mpf(K)))*4*np.pi
                if t<0: t*=2
                res.append(v-t)
    print('beta=%s  N=8 cross %s  N=9 cross %s   rms(fig - model, K<0, N=6..9) = %.3f   median|resid| = %.4f' % (beta, mp.nstr(c8,6), mp.nstr(c9,6), np.sqrt(np.mean(np.square(res))), np.median(np.abs(res))))
