"""Independent linearisation of the tanh-model regular N-gon in the Poincare disk.
Cartesian coordinates, exact (sympy-generated) Hessian of the pair term, mpmath eigenvalues.
F = H - Omega*J, L = Lambda^{-1} S Hess(F), S = [[0,1],[-1,0]] per vortex, Lambda = 4/(1-|z|^2)^2."""
import sys, sympy as sp, mpmath as mp
mp.mp.dps = 60
x1,y1,x2,y2 = sp.symbols('x1 y1 x2 y2', real=True)
# pair term: -ln|z-w| + ln|1 - z conj(w)|
g = -sp.log((x1-x2)**2+(y1-y2)**2)/2 + sp.log((1-(x1*x2+y1*y2))**2 + (x2*y1-x1*y2)**2)/2
v = [x1,y1,x2,y2]
gH = sp.lambdify(v, sp.hessian(g, v), 'mpmath')
gG = sp.lambdify(v, [sp.diff(g,a) for a in v], 'mpmath')
x,y = sp.symbols('x y', real=True)
Jf = (1+x**2+y**2)/(1-x**2-y**2)
jH = sp.lambdify([x,y], sp.hessian(Jf,[x,y]), 'mpmath')
jG = sp.lambdify([x,y], [sp.diff(Jf,x), sp.diff(Jf,y)], 'mpmath')

def ngon(N, xi):
    R = mp.sqrt(xi)
    return [(R*mp.cos(2*mp.pi*k/N), R*mp.sin(2*mp.pi*k/N)) for k in range(N)]

def linop(N, xi, return_parts=False):
    P = ngon(N, xi); n = 2*N
    gradH = [mp.mpf(0)]*n; HH = mp.zeros(n,n); HJ = mp.zeros(n,n); gradJ=[mp.mpf(0)]*n
    for j in range(N):
        for k in range(j+1,N):
            a = list(P[j])+list(P[k]); idx=[2*j,2*j+1,2*k,2*k+1]
            Hs = gH(*a); gs = gG(*a)
            for p in range(4):
                gradH[idx[p]] += gs[p]
                for q in range(4): HH[idx[p],idx[q]] += Hs[p,q]
        h = jH(*P[j]); gj = jG(*P[j])
        for p in range(2):
            gradJ[2*j+p] = gj[p]
            for q in range(2): HJ[2*j+p,2*j+q] = h[p,q]
    Om = gradH[0]/gradJ[0]
    resid = max(abs(gradH[i]-Om*gradJ[i]) for i in range(n))
    HF = HH - Om*HJ
    L = mp.zeros(n,n)
    for j in range(N):
        lam = 4/(1-xi)**2
        # rows: xdot = (1/lam) dF/dy, ydot = -(1/lam) dF/dx
        for q in range(n):
            L[2*j,q] = HF[2*j+1,q]/lam
            L[2*j+1,q] = -HF[2*j,q]/lam
    if return_parts: return L, HF, Om, resid
    return L, Om, resid

def maxre(N, xi):
    L,_,_ = linop(N, xi)
    ev = mp.eig(L, left=False, right=False)
    return max(abs(mp.re(e)) for e in ev)

def threshold(N, lo, hi, tol=mp.mpf(10)**-22):
    # stable for xi > xi*: maxre(lo)>tol, maxre(hi)<tol
    for _ in range(140):
        mid = (lo+hi)/2
        if maxre(N, mid) > tol: lo = mid
        else: hi = mid
    return (lo+hi)/2

if __name__ == '__main__':
    # coarse scan to locate sign change then bisect
    for N in [int(a) for a in sys.argv[1:]] or [8,9,10,12]:
        grid = [mp.mpf(i)/200 for i in range(1,200)]
        prev=None; br=None
        for xi in grid:
            s = maxre(N, xi) > mp.mpf(10)**-22
            if prev is not None and prev and not s: br=(xi-mp.mpf(1)/200, xi)
            if prev is not None and (not prev) and s: print(N,'re-destabilises near',xi)
            prev=s
        L,Om,res = linop(N, br[0])
        print('N=%d bracket %s resid %s' % (N, mp.nstr(br[0],5), mp.nstr(res,3)))
        xs = threshold(N, *br)
        print('N=%d xi* = %s' % (N, mp.nstr(xs, 45)))
        sys.stdout.flush()
