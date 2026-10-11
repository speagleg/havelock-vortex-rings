"""N=5 (black curve) re-digitization avoiding the frame, dotted lines and the 'N=5' label; compare with model."""
import numpy as np, mpmath as mp, sympy as sp
from PIL import Image
import boatto_hessian as bh
im = np.asarray(Image.open('fig/fig5a_crop.png').convert('RGB')).astype(int)
dark = (im[...,0]<90)&(im[...,1]<90)&(im[...,2]<90)
x0,x1,y0,y1 = 186,1179,195,1188
Kof=lambda x:-0.6+1.2*(x-x0)/(x1-x0); vof=lambda y:5-10*(y-y0)/(y1-y0)
xK=lambda K: x0+(K+0.6)/1.2*(x1-x0); yv=lambda v: y0+(5-v)/10*(y1-y0)
bh.FORMS[('x',-1)] = bh.c*(sp.log(bh.d-1)+sp.log(bh.d+1)); bh.FORMS[('x',1)] = bh.c*(sp.log(1-bh.d)-sp.log(1+bh.d))
for K in [-0.55,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.44,0.46,0.48,0.52,0.56,0.58]:
    m = float(bh.min_nontrivial(5, mp.mpf(K), 'x')[0])*4*np.pi; m2 = 2*m if m<0 else m
    x=int(round(xK(K)))
    # search dark pixels within +-60 px of predicted model y, excluding the dotted zero line row band
    yp = yv(m2); ys=[y for y in range(max(y0+5,int(yp)-60), min(y1-5,int(yp)+60)) if dark[y,x] and abs(y-yv(0))>3]
    fig = vof(np.median(ys)) if ys else float('nan')
    print('K=%+.2f  fig=%6.3f  model=%6.3f  model(neg x2)=%6.3f' % (K, fig, m, m2))
