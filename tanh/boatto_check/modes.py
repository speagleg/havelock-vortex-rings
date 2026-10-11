"""Per-Fourier-mode eigenvalues of S~ (Boatto variables), x 4 pi, to see which mode is the minimum."""
import sys, mpmath as mp
from boatto_hessian import hessian
mp.mp.dps = 30
def mode_eigs(N, K, model='boatto'):
    S, _ = hessian(N, K, model)
    out = []
    for m in range(0, N//2+1):
        # complex Fourier vector e_m on q and on p
        v = [mp.exp(2j*mp.pi*m*a/N)/mp.sqrt(N) for a in range(N)]
        def blk(r0, c0):
            return sum(mp.conj(v[a])*S[r0+a, c0+b]*v[b] for a in range(N) for b in range(N))
        Qm, Pm, Gm = blk(0,0), blk(N,N), blk(0,N)
        M = mp.matrix([[Qm, Gm],[mp.conj(Gm), Pm]])
        E = mp.eighe(M, eigvals_only=True) if hasattr(mp,'eighe') else None
        tr = mp.re(Qm+Pm); dt = mp.re(Qm*Pm) - abs(Gm)**2
        lo = (tr - mp.sqrt(tr**2-4*dt))/2; hi=(tr + mp.sqrt(tr**2-4*dt))/2
        out.append((m, 4*mp.pi*mp.re(Qm), 4*mp.pi*mp.re(Pm), 4*mp.pi*abs(Gm), 4*mp.pi*lo, 4*mp.pi*hi))
    return out
if __name__ == '__main__':
    model = sys.argv[1] if len(sys.argv)>1 else 'boatto'
    for N in [6,8,9]:
        for K in ['-0.5','-0.3','-0.25','-0.1','0.1','0.3']:
            print('N=%d K=%s' % (N, K))
            for m,q,p,g,lo,hi in mode_eigs(N, mp.mpf(K), model):
                print('   m=%d  4piQ=%9.4f 4piP=%9.4f 4pi|G|=%8.4f  eig: %9.4f %9.4f' % (m,q,p,g,lo,hi))
