# Large-N asymptotics of the tanh threshold (NUMERICAL, 40 digits), and model comparison table.
import mpmath as mp
from blocks import blocks
mp.mp.dps=40
# limit equation in kappa = cosh(r*)/N, u = exp(-1/kappa)
def lim(k):
    u=mp.e**(-1/k); return k-k**2/2+2*k*u**2/(1-u**2)-mp.mpf(1)/8-u/(2*(1-u)**2)
k0=1-mp.sqrt(3)/2
kinf=mp.findroot(lim,k0)
print('kappa_inf (u-corrected) =',kinf,'  algebraic part 1-sqrt3/2 =',k0,'  u=e^{-1/kappa}=',mp.e**(-1/kinf))
def xistar(N):
    m=N//2; f=lambda x: blocks(N,x,m)[0]
    c=kinf*N; g=(c-1)/(c+1) if N>30 else mp.mpf('0.3')
    return mp.findroot(f,g)
def sinh_xi(N):
    q=mp.mpf(N-3)/4 if N%2 else mp.mpf((N-2)**2)/(4*(N-1)); return (mp.sqrt(q)-1)/(mp.sqrt(q)+1)
print('\nN   xi*_tanh            cosh r*_tanh      cosh r*/N     cosh r*-kinf*N   | xi*_sinh   cosh r*_sinh  | min_m a_m/N^2 at 1.001*xi* (all modes)')
prev=None
for N in [8,9,10,12,16,20,30,40,60,80,100,150,200,300,400,600,800]:
    xs=xistar(N); c=(1+xs)/(1-xs); xsh=sinh_xi(N); csh=(1+xsh)/(1-xsh)
    # all-mode check just above threshold (numerical)
    x1=xs*(1+mp.mpf('1e-3')) if xs<0.99 else 1-(1-xs)*(1-mp.mpf('1e-3'))
    mn=min(blocks(N,x1,m)[0] for m in range(1,N))
    mb=min(blocks(N,x1,m)[1] for m in range(1,N))
    print('%-4d %s %s %s %s | %s %s | %s  min b %s'%(N,mp.nstr(xs,18),mp.nstr(c,14),mp.nstr(c/N,10),mp.nstr(c-kinf*N,10),mp.nstr(xsh,10),mp.nstr(csh,8),mp.nstr(mn/N**2,4),mp.nstr(mb,4)))
