"""Digitize Fig 5b (N=6 min/max eigenvalue) and compare with tanh law and sinh(D) law (K<0), tan law (K>0).
Plotted-value model: min -> min(lam, 2 lam), max -> 2 lam  (the 'min over {lam_i} U {2 lam_i}' artefact)."""
import numpy as np, sympy as sp, mpmath as mp
from PIL import Image
import boatto_hessian as bh
mp.mp.dps = 20
bh.FORMS[('sinhD',-1)] = bh.c*(sp.log(bh.d-1)+sp.log(bh.d+1)); bh.FORMS[('sinhD',1)] = bh.FORMS[('boatto',1)]
im = np.asarray(Image.open('fig/fig5b_crop.png').convert('RGB')).astype(int)
R,G,B = im[...,0],im[...,1],im[...,2]
dark=(R<90)&(G<90)&(B<90); red=(R>180)&(G<90)&(B<90)
H,W=R.shape
cols=dark.sum(0); rows=dark.sum(1)
vl=[x for x in range(W) if cols[x]>0.5*H*0.7]; hl=[y for y in range(H) if rows[y]>0.5*W*0.7]
x0,x1,y0,y1=min(vl),max(vl),min(hl),max(hl); print('frame',x0,x1,y0,y1)
Kof=lambda x:-0.8+1.6*(x-x0)/(x1-x0); vof=lambda y:25-35*(y-y0)/(y1-y0); xK=lambda K:int(round(x0+(K+0.8)/1.6*(x1-x0)))
for K in [-0.75,-0.5,-0.25,-0.1,0.1,0.25,0.4,0.5,0.6,0.7]:
    x=xK(K)
    ysr=[y for y in np.nonzero(red[y0+5:y1-5,x])[0]+y0+5 if not (K>-0.72 and K<-0.05 and vof(y)>4.3)]  # skip label text
    ysd=[y for y in np.nonzero(dark[y0+5:y1-5,x])[0]+y0+5 if vof(y)>15 and vof(y)<(20.9 if K<0.1 else 23.5)]  # skip label text
    fmin=vof(np.median(ysr)) if len(ysr) else np.nan; fmax=vof(np.median(ysd)) if ysd else np.nan
    out=[]
    for model in ['boatto','sinhD']:
        mn,mx,_,_,_=bh.min_nontrivial(6, mp.mpf(K), model); mn=float(mn)*4*np.pi; mx=float(mx)*4*np.pi
        out.append('%s min %6.2f max %6.2f' % (model, min(mn,2*mn), 2*mx))
    print('K=%+.2f fig min %6.2f max %6.2f | %s' % (K,fmin,fmax,' | '.join(out)))
