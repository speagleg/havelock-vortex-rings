import sys, mpmath as mp
sys.path.insert(0,'..'); from blocks import blocks
mp.mp.dps=40
f=lambda k: k-k**2/2+2*k*mp.e**(-2/k)/(1-mp.e**(-2/k))-mp.mpf(1)/8-mp.e**(-1/k)/(2*(1-mp.e**(-1/k))**2)
k=mp.findroot(f,0.134); print('kappa_inf =',mp.nstr(k,25))
for N in [16,64,256,1024]:
    xs=mp.findroot(lambda x: blocks(N,x,N//2)[0], 1-2/(k*N))
    c=(1+xs)/(1-xs); print(N,'cosh r*/N =',mp.nstr(c/N,15),' (cosh r* - kappa N)*N =',mp.nstr((c-k*N)*N,6))
