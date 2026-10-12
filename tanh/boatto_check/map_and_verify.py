"""Precise zero crossings of Boatto's S~ (reimplemented from her equations, her variables), mapped to
xi = tanh^2(sqrt|K| r_o / 2) with r_o = 1, and compared with ../blocks.py (tanh).
Mode-resolved: for mode m we use det of the Hermitian 2x2 block (Q_m, G_m; G_m*, P_m), whose sign change
(with Q_m > 0) is exactly the definiteness boundary."""
import sys, sympy as sp, mpmath as mp
sys.path.insert(0, '..')
import boatto_hessian as bh
from blocks import blocks
mp.mp.dps = 40
def modeblock(N, K, model, m):
    S,_ = bh.hessian(N, K, model)
    v = [mp.exp(2j*mp.pi*m*a/N)/mp.sqrt(N) for a in range(N)]
    blk = lambda r0,c0: sum(mp.conj(v[a])*S[r0+a,c0+b]*v[b] for a in range(N) for b in range(N))
    return blk(0,0), blk(N,N), blk(0,N)
def det(N,K,model,m):
    Q,P,G = modeblock(N,K,model,m); return mp.re(Q*P) - abs(G)**2
def Kcross(N, model, guess):
    m = N//2
    return mp.findroot(lambda K: det(N,K,model,m), (mp.mpf(guess)*mp.mpf('1.25'), mp.mpf(guess)*mp.mpf('0.7')), solver='anderson', tol=mp.mpf(10)**-30)
xi_of = lambda K: mp.tanh(mp.sqrt(-K)/2)**2
K_of = lambda xi: -(2*mp.atanh(mp.sqrt(xi)))**2
print('Map (r_o = 1 geodesic, her Remark 2 p.361-362):  xi = tanh^2(sqrt(-K)/2),  K = -(2 artanh sqrt(xi))^2')
for N, g in [(8,-0.28),(9,-0.49),(10,-0.7)]:
    m = N//2
    xs = mp.findroot(lambda x: blocks(N,x,m)[0], {8:0.067,9:0.11,10:0.166}[N])
    Kt = Kcross(N,'boatto',g)
    print('N=%d' % N)
    print('   ours (blocks.py) xi* = %s  ->  K = %s' % (mp.nstr(xs,25), mp.nstr(K_of(xs),20)))
    print('   reimpl. tanh (her eq.17)  K* = %s  ->  xi = %s   |diff xi| = %s' % (mp.nstr(Kt,20), mp.nstr(xi_of(Kt),25), mp.nstr(abs(xi_of(Kt)-xs),3)))
# per-mode check: sign/ratio of reimplemented P_m vs our a_m at K=-0.25, N=8 (G_4 = 0 so P_4 is the eigenvalue)
K = mp.mpf('-0.25'); xi = xi_of(K)
print('K=-0.25 -> xi = %s' % mp.nstr(xi,15))
for N in [8,9]:
    for m in range(1, N//2+1):
        Q,P,G = modeblock(N,K,'boatto',m); a,b,g,_ = blocks(N,xi,m)
        print('   N=%d m=%d  ours a=%+.6f b=%+.6f g=%+.6f a*b-g^2=%+.6f | Boatto-vars 4piP=%+.6f 4piQ=%+.6f 4pi|G|=%.6f det-sign %+d  ratio P/a=%.6f Q/b=%.6f' % (
            N,m,a,b,g,a*b-g*g, 4*mp.pi*mp.re(P), 4*mp.pi*mp.re(Q), 4*mp.pi*abs(G), int(mp.sign(mp.re(Q*P)-abs(G)**2)), mp.re(P)/a, mp.re(Q)/b))
