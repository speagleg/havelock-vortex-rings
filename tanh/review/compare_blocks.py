"""Compare candidate Fourier blocks against Fourier projection of MY Cartesian Hessian (indep_lin.py)."""
import sys, mpmath as mp
sys.path.insert(0,'..')
from indep_lin import linop
from blocks import blocks
mp.mp.dps=60
worst=0
for N in [3,5,7,8,9,10,12,13]:
    for xi in ['0.03','0.11','0.4','0.8']:
        xi=mp.mpf(xi); R=mp.sqrt(xi)
        L,HF,Om,res=linop(N,xi,True)
        for m in range(1,N):
            V=mp.zeros(2*N,2)
            for j in range(N):
                th=2*mp.pi*j/N; e=mp.expj(m*th)
                V[2*j,0]=e*R*mp.cos(th); V[2*j+1,0]=e*R*mp.sin(th)      # d/ds
                V[2*j,1]=-e*R*mp.sin(th); V[2*j+1,1]=e*R*mp.cos(th)     # d/dtheta
            B=(V.H*HF*V)/N
            a,b,g,E=blocks(N,xi,m)
            d=max(abs(B[0,0]-a),abs(B[1,1]-b),abs(B[0,1]-1j*g),abs(B[1,0]+1j*g))
            d2=max(abs(B[0,0]-a),abs(B[1,1]-b),abs(B[0,1]+1j*g),abs(B[1,0]-1j*g))
            worst=max(worst,min(d,d2))
        kap=4*xi/(1-xi)**2
        worst=max(worst,abs(Om*kap-E))
print('max |my Fourier block - candidate block| (up to sign convention of g), and Omega check:',mp.nstr(worst,5))
