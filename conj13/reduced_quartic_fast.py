"""Conj13, odd n, MT family: fast reduced quartic.  Same as reduced_quartic.py, but the minimisation is
restricted to the Fourier modes that can couple to a^2 at cubic order (D_n-invariance, MT Prop. 12(iii)
plus the V_0 radial direction): modes 0 and 1 (radial and tangential) plus the second kernel
vector.  The other modes enter f_eff only at O(eps^6), so c4 is unchanged (checked against the full
minimisation of reduced_quartic.py for n = 5..11).  Kernel = mode (n-1)/2 radial (MT: the Hessian is
diagonal in (r, theta) for this model, eq. (2.2)/(2.6)); verified numerically: gradient and Hessian
action on the kernel vectors are printed.  NUMERICAL (mpmath 40 digits).
usage: python reduced_quartic_fast.py n [n ...]
"""
import sys, mpmath as mp
mp.mp.dps = 40
from engine import Ring, fit_even_poly
from reduced_quartic import b_threshold

def modevecs(n, l):
    out = []
    for trig in (mp.cos, mp.sin):
        for kind in ('r', 't'):
            v = [mp.mpf(0)]*(2*n)
            for k in range(n):
                ph = 2*mp.pi*k/n
                c = trig(l*ph)
                if kind == 'r':
                    v[k] += c*mp.cos(ph); v[n+k] += c*mp.sin(ph)
                else:
                    v[k] += -c*mp.sin(ph); v[n+k] += c*mp.cos(ph)
            out.append(v)
    return out

def run(n):
    b = b_threshold(n)
    lam, r0 = (b, 1) if b != 0 else (0, 1)
    R = Ring(n, lam, r0, 'MT')
    lc = (n-1)//2
    kc, ks, _, _ = [mp.matrix(v) for v in modevecs(n, lc)][0], *[mp.matrix(v) for v in modevecs(n, lc)][1:]
    vc = modevecs(n, lc)
    e1 = mp.matrix(vc[0]); e2 = mp.matrix(vc[2])        # cos-radial, sin-radial
    e1 /= mp.norm(e1); e2 /= mp.norm(e2)
    # sanity: kernel of the Hessian of Hhat (Hessian-vector product by central differences of the gradient)
    h = mp.mpf(10)**-15
    def hv(v):
        qp = [R.q0[i] + h*v[i] for i in range(2*n)]; qm = [R.q0[i] - h*v[i] for i in range(2*n)]
        gp, gm = R.gradH(qp), R.gradH(qm)
        return mp.norm(mp.matrix([(gp[i]-gm[i])/(2*h) for i in range(2*n)]))
    print(f"== n={n} b*={mp.nstr(b, 20)}: |Hess.e1|={mp.nstr(hv(e1), 3)} |Hess.e2|={mp.nstr(hv(e2), 3)} |grad Hhat(q0)|={mp.nstr(mp.norm(mp.matrix(R.gradH(R.q0))), 3)}")
    G = R.gradJ(R.q0); t = R.rot_gen(R.q0)
    Nb = Ring._gs([t] + G)
    normals = Nb[1:]
    S = modevecs(n, 0) + modevecs(n, 1) + [list(e2)]
    F = Ring._gs(S, basis=Nb + [e1])
    for name, e, extra in (('e1', e1, []),):
        eps = [mp.mpf(k)/2000 for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5)]
        vals = [R.f_eff(ep, e, F, normals)[0] for ep in eps]
        c = fit_even_poly(eps, vals, list(range(2, 12)))
        print(f"   constrained J=mu : c2={mp.nstr(c[0], 5)} c3={mp.nstr(c[1], 5)} c4={mp.nstr(c[2], 18)}  "
              f"identify: {mp.identify(c[2])}")
        sys.stdout.flush()

if __name__ == '__main__':
    for a in sys.argv[1:]:
        run(int(a))
