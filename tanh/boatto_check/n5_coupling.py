"""Side issue (K>0, N=5): does over-weighting the Q-P coupling G by a factor k reproduce her N=5 curve near K~0.45?
Plotted model = min over modes of 1/2[(Q+P) - sqrt((Q-P)^2 + 4 k^2 |G|^2)] (x 4pi), negatives doubled."""
import mpmath as mp, numpy as np
from map_and_verify import modeblock
mp.mp.dps = 20
fig = {0.30:0.916, 0.40:0.368, 0.44:0.070, 0.52:-0.780, 0.56:-1.309, 0.58:-1.581}
blocks = {K: [modeblock(5, mp.mpf(K), 'boatto', m) for m in (1,2)] for K in fig}
for k in [1, 1.5, 2, 2.5, 3]:
    row=[]; 
    for K,v in fig.items():
        vals=[]
        for Q,P,G in blocks[K]:
            q,p,g = float(mp.re(Q))*4*np.pi, float(mp.re(P))*4*np.pi, float(abs(G))*4*np.pi
            vals.append(0.5*((q+p)-np.sqrt((q-p)**2+4*k*k*g*g)))
        m=min(vals); m = 2*m if m<0 else m
        row.append('K=%.2f fig %.3f model %.3f' % (K, v, m))
    print('k=%.1f: ' % k + ' | '.join(row))
for K in [0.3,0.4,0.45]:
    print('K=%.2f modes (4piQ,4piP,4pi|G|):' % K, [(round(float(mp.re(Q))*4*np.pi,3), round(float(mp.re(P))*4*np.pi,3), round(float(abs(G))*4*np.pi,3)) for Q,P,G in [modeblock(5, mp.mpf(K),'boatto',m) for m in (1,2)]])
