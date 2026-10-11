# Numerical scan of the tanh-model Fourier blocks: sign structure, critical mode, spectral vs energy-momentum thresholds.
import mpmath as mp
from blocks import blocks
mp.mp.dps=30
def status(N,xi):
    spec=True; pd=True; bad=[]; badpd=[]
    for m in range(1,N):
        a,b,g,_=blocks(N,xi,m)
        if a*b<=0: spec=False; bad.append(m)
        if not (a>0 and a*b-g*g>0): pd=False; badpd.append(m)
    return spec,pd,bad,badpd
for N in range(2,21):
    grid=[mp.mpf(k)/2000 for k in range(1,2000)]
    sp_tr=[];pd_tr=[];minb=mp.inf;modes=set();pdmodes=set()
    prev=None
    for xi in grid:
        s,p,bad,badpd=status(N,xi); modes|=set(bad); pdmodes|=set(badpd)
        for m in range(1,N): minb=min(minb,blocks(N,xi,m)[1])
        if prev is not None:
            if s!=prev[0]: sp_tr.append(float(xi))
            if p!=prev[1]: pd_tr.append(float(xi))
        prev=(s,p)
    print(N,'spectral transitions',sp_tr,'PD transitions',pd_tr,'modes ab<=0',sorted(modes),'modes not PD',sorted(pdmodes),'min b_m on grid %.4g'%float(minb))
