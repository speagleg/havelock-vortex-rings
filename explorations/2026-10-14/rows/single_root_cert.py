"""Computer-assisted proof (Arb ball arithmetic, python-flint) that for the tanh model the threshold function F(pi) of a
hypercycle row has exactly ONE root in l' for each kappa in (1/sqrt2, 1): the stability boundary is a single curve.
Criterion: at every zero of F with alpha = arccos(kappa) in (0, pi/4), dF/dl' (alpha fixed) > 0. Then F(alpha, .) can only
cross zero upwards, so it has at most one root; it has one, since F < 0 for small l' (C1 below) and F > 0 for
cosh l' >= 1/(2 kappa^2 - 1) (stability_bound.py). For alpha >= pi/4 there are no zeros (instability_proof_check.py).
Notation: r = e^{-l'}, k = cos(alpha), s = sin(alpha); G = F/s^2 (same zeros, same sign of dF/dl' at zeros) has the Lambert form
  G = sum_n a_n r^n/(1-r^n) - b_n r^n/(1+r^n),  a_n = 2k U_{2n-1}(k) - 2n sin^2(n alpha),  b_n = 2n cos^2(n alpha),
with |a_n| <= 6n, |b_n| <= 2n (U Chebyshev of the 2nd kind; U_{2n-1}(cos a) = 2n sinc(2na)/sinc(a) is used,
since the three-term recurrence inflates ball radii). Regions covering every zero:
  A  l' >= 6: analytic with explicit remainder: G = P1 r + P2 r^2 + E, P1 = a1-b1, P2 = a1+b1+a2-b2; at a zero
     dG/dl' = -P2 r^2 + E - r E_r, and -P2 > (|E| + |r E_r|)/r^2 is checked in Arb (zeros have 2k^2-1 < 1/cosh l').
  B  0.25 <= l' <= 6, 0 <= alpha <= pi/4: adaptive boxes; on each, G excludes 0 or dG/dl' > 0 (with tail bounds).
  C  l' < 0.25, eta = pi alpha/l', using the dual series F = F0 - R1 - R2 of dual_series.py (R1, R2 >= 0):
     C1 eta >= 7: F <= F0 < 0 analytically.   C2 eta <= 2.5: F >= h(eta) - eta^2 [(2t^2/3)(1/2+eta) + 2t/pi] - R > 0,
     t = 0.25/pi, h = horocycle function (h/eta^2 bounded below by mean-value enclosures; near 0 by the Bernoulli series).   C3 2.5 <= eta <= 7 (alpha <= 7*0.25/pi): adaptive boxes in (alpha, eta);
     on each, F0 - Rmax > 0, or F0 < 0, or dF0/deta + max|dR/deta| < 0 (dF/dl' = -(eta/l') dF/deta at fixed alpha).
Usage: python single_root_cert.py   (about a minute)"""
import math, time
from flint import arb, ctx
ctx.prec = 96
PI = arb.pi()

def iv(a, b):
    a, b = float(a), float(b); m = (a + b)/2; r = (b - a)/2 + 4e-16*max(abs(a), abs(b), 1e-300)
    return arb(m, r)
def lo(x): return float(x.mid()) - float(x.rad())
def hi(x): return float(x.mid()) + float(x.rad())
def pos(x): return x > 0
def neg(x): return x < 0

def tail_pow(M, k, r):
    """upper bound for sum_{n>=M} n^k r^n (ratio bound), r a float upper bound."""
    q = r*(1 + 1/M)**k
    assert q < 1
    return M**k * r**M/(1 - q)

# ---------------- Lambert form of G and dG/dl' on a box ----------------
def G_and_dG(al, lp, N):
    k, s = al.cos(), al.sin()
    r = (-lp).exp(); rmax = math.exp(-lo(lp))
    sa = al.sinc()                          # U_{2n-1}(cos a) = sin(2na)/sin(a) = 2n sinc(2na)/sinc(a): no recurrence blow-up
    G = arb(0); dG = arb(0)
    for n in range(1, N+1):
        sn, cn = (n*al).sin(), (n*al).cos()
        a = 2*k*(2*n)*(2*n*al).sinc()/sa - 2*n*sn*sn; b = 2*n*cn*cn
        rn = r**n
        G += a*rn/(1 - rn) - b*rn/(1 + rn)
        dG += -a*n*rn/(1 - rn)**2 + b*n*rn/(1 + rn)**2
    tG = 8*tail_pow(N+1, 1, rmax)/(1 - rmax); tD = 8*tail_pow(N+1, 2, rmax)/(1 - rmax)**2
    return G + arb(0, tG), dG + arb(0, tD)

def region_B(l0=0.25, l1=6.0, a1=math.pi/4, maxdepth=40):
    stack = [(0.0, a1, l0, l1, 0)]; nbox = 0; nderiv = 0
    while stack:
        A0, A1, L0, L1, d = stack.pop(); nbox += 1
        assert d < maxdepth, ('depth exceeded', A0, A1, L0, L1)
        N = int(40/L0) + 8
        G, dG = G_and_dG(iv(A0, A1), iv(L0, L1), N)
        if pos(G) or neg(G): continue
        if pos(dG): nderiv += 1; continue
        # split the dimension with larger scaled width
        if (A1 - A0)*8 > (L1 - L0):
            m = (A0 + A1)/2; stack += [(A0, m, L0, L1, d+1), (m, A1, L0, L1, d+1)]
        else:
            m = (L0 + L1)/2; stack += [(A0, A1, L0, m, d+1), (A0, A1, m, L1, d+1)]
    return nbox, nderiv

# ---------------- region A ----------------
def region_A(L=6.0):
    r = math.exp(-L)
    Kmax = 0.5 + 0.5/math.cosh(L) + 1e-9          # widened against float rounding
    al = iv(math.acos(math.sqrt(Kmax)), math.pi/4)
    k = al.cos(); U1 = 2*k; U3 = 8*k**3 - 4*k
    a1 = 2*k*U1 - 2*al.sin()**2; b1 = 2*al.cos()**2
    a2 = 2*k*U3 - 4*(2*al).sin()**2; b2 = 4*(2*al).cos()**2
    P2 = a1 + b1 + a2 - b2
    E = r**3*(6/(1-r) + 2) + 8*tail_pow(3, 1, r) + (6/(1-r**2) + 2)*tail_pow(2, 1, r*r)
    rE = 3*r**3*(6/(1-r) + 2) + r**4*(6/(1-r)**2 + 2) + 8*tail_pow(3, 2, r) + 24/(1-r*r)**2*tail_pow(2, 2, r*r)
    need = (E + rE)/r**2
    return P2, need

# ---------------- region C ----------------
def S_enc(x0, x1):
    """enclosure of sinc on [x0, x1] subset [0, pi] (decreasing)."""
    v1 = arb(x1).sinc() if x1 > 0 else arb(1); v0 = arb(x0).sinc() if x0 > 0 else arb(1)
    return arb((lo(v1) + hi(v0))/2, (hi(v0) - lo(v1))/2 + 1e-30)
def F0_dF0(A0, A1, E0, E1):
    al, eta = iv(A0, A1), iv(E0, E1)
    S2, S1 = S_enc(2*A0, 2*A1), S_enc(A0, A1)
    c2, c1 = (2*eta).coth(), eta.coth()
    F0 = eta*S2*c2 - eta**2/4*S1**2*c1**2 - eta/PI*(2*al).sin() - (2*al).cos()/4
    dF0 = S2*(c2 - 2*eta/(2*eta).sinh()**2) - S1**2/4*(2*eta*c1**2 - 2*eta**2*c1/eta.sinh()**2) - (2*al).sin()/PI
    return F0, dF0
def R_bounds(lmax, amax):
    """bounds on R1 + R2 and on |dR/deta| (any alpha <= amax, any l' <= lmax), from the dual series."""
    rho1 = math.exp(-4*math.pi*(math.pi - amax)/lmax); rho2 = math.exp(-2*math.pi*(math.pi - amax)/lmax)
    R1 = (2*math.pi/lmax)*(0.5 + 4*math.pi/lmax)*tail_pow(1, 1, rho1)/(1 - rho1)
    R2 = (4*math.pi**2/lmax**2)*rho2/(1 - rho2)**3
    fac = 3/lmax + 8*math.pi**2*2/lmax**2          # generous |d/dl'| factor (n-weighted sums absorbed by n^2 tails)
    R1d = (2*math.pi/lmax)*(0.5 + 4*math.pi/lmax)*tail_pow(1, 2, rho1)/(1 - rho1)*fac
    R2d = (4*math.pi**2/lmax**2)*tail_pow(1, 2, rho2)/(1 - rho2)*fac
    # dR/deta = -(l'/eta) dR/dl', and l'/eta = pi alpha/eta^2 <= pi amax/2.5^2
    return R1 + R2, (R1d + R2d)*math.pi*amax/2.5**2
def region_C3(lc=0.25, e0=2.5, e1=7.0, maxdepth=40):
    amax = 7*lc/math.pi; lmax = math.pi*amax/e0
    Rmax, dRmax = R_bounds(lmax, amax)
    stack = [(0.0, amax, e0, e1, 0)]; nbox = 0; nderiv = 0
    while stack:
        A0, A1, E0, E1, d = stack.pop(); nbox += 1
        assert d < maxdepth, ('depth exceeded', A0, A1, E0, E1)
        F0, dF0 = F0_dF0(A0, A1, E0, E1)
        if pos(F0 - Rmax) or neg(F0): continue
        if neg(dF0 + dRmax): nderiv += 1; continue
        if (A1 - A0)*8 > (E1 - E0):
            m = (A0 + A1)/2; stack += [(A0, m, E0, E1, d+1), (m, A1, E0, E1, d+1)]
        else:
            m = (E0 + E1)/2; stack += [(A0, A1, E0, m, d+1), (A0, A1, m, E1, d+1)]
    return nbox, nderiv, Rmax, dRmax
def _g(e):        # h(eta)/eta^2, h = -eta^2 coth^2(eta)/4 + eta coth(2 eta) - 1/4
    return -e.coth()**2/4 + (2*e).coth()/e - 1/(4*e**2)
def _gp(e):       # derivative of _g
    return e.coth()/(2*e.sinh()**2) - 2/((2*e).sinh()**2*e) - (2*e).coth()/e**2 + 1/(2*e**3)
def region_C2(lc=0.25, e1=2.5, d0=0.05):
    """h(eta)/eta^2 > (2t^2/3)(1/2+eta) + 2t/pi + R/eta^2 on (0, e1], t = lc/pi. Returns the smallest margin found
    (before subtracting R) and the R/eta^2 bound."""
    t = lc/math.pi
    # (0, d0]: h/eta^2 >= 1/2 - (1/36 + 8/45) eta^2 from the alternating Bernoulli series of x coth x (|x| < pi)
    worst = 0.5 - (1/36 + 8/45)*d0**2 - ((2*t*t/3)*(0.5 + d0) + 2*t/math.pi)
    stack = [(d0, e1)]
    while stack:
        x0, x1 = stack.pop()
        e = iv(x0, x1); m = arb((x0 + x1)/2)
        gl = _g(m) + _gp(e)*(e - m)                  # mean-value enclosure
        margin = lo(gl) - ((2*t*t/3)*(0.5 + x1) + 2*t/math.pi)
        if margin > 0.02 or x1 - x0 < 1e-4:
            worst = min(worst, margin); continue
        stack += [(x0, (x0 + x1)/2), ((x0 + x1)/2, x1)]
    # R1/eta^2 <= 16 sum n e^{-4n(pi/t - eta)}/(1 - Q1), R2/eta^2 <= 4 rho/(1-rho)^3, rho = e^{-(2pi/t - 2 eta)}
    R1c = 2*16*tail_pow(1, 1, math.exp(-4*(math.pi/t - e1))); rho = math.exp(-(2*math.pi/t - 2*e1)); R2c = 4*rho/(1 - rho)**3
    return worst, R1c + R2c
if __name__ == '__main__':
    t0 = time.time()
    P2, need = region_A()
    print(f'A  (l\' >= 6):  -P2 in [{lo(-P2):.6f}, {hi(-P2):.6f}],  needed > {need:.4f}:  {lo(-P2) > need}')
    nb, nd = region_B()
    print(f'B  (0.25 <= l\' <= 6): {nb} boxes, {nd} decided by dG/dl\' > 0, all others by the sign of G  [{time.time()-t0:.0f} s]')
    nb, nd, Rm, dRm = region_C3()
    print(f'C3 (2.5 <= eta <= 7, l\' < 0.25): {nb} boxes, {nd} by dF0/deta < 0; R <= {Rm:.2e}, |dR/deta| <= {dRm:.2e}')
    w, Rc = region_C2()
    print(f'C2 (eta <= 2.5): min over eta of h/eta^2 - bound >= {w:.4f} > R/eta^2 <= {Rc:.2e}: {w > Rc}')
    print('C1 (eta >= 7): F <= F0 <= (pi s/l\')[coth 2eta - (eta/2pi) coth^2 eta] < 0, since 7/(2pi) > coth 14:',
          7/(2*math.pi) > 1/math.tanh(14))
    print(f'total {time.time()-t0:.0f} s')
