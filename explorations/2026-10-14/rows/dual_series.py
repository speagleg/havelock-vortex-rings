"""Hypercycle rows, tanh model: the threshold function F(pi) in theta functions, its modular (dual) series, and an
elementary closed-form threshold.
(1) Theta form (r = e^{-l'}):  using R^+(r) = R^-(r) - 2 R^-(r^2) the Lambert series of qpi_proof.py collapses to
      F(pi) = sin2a S(a; r) - 2 s^2 [E(r) - E(r^2)] + 2 s^2 C(a; r^2),
    S = sum sin(2na) r^n/(1-r^n) = [(ln th1)'(a | e^{-l'/2}) - cot a]/4,  C(a; r^2) = [(ln th1)''(a | e^{-l'}) + csc^2 a]/8,
    E(r) = sum n r^n/(1-r^n)  (th1(z|q) Jacobi theta with nome q).
(2) Modular transformation (tau -> -1/tau) gives the exact dual series, with eta = pi a/l', Q1 = e^{-4pi^2/l'}, Q2 = e^{-2pi^2/l'}:
      F(pi) = F0 - (2pi/l') sum_n Q1^n/(1-Q1^n) [sin2a sinh(4 pi n a/l') - 4 pi s^2 n/l']
                 - (2 pi^2 s^2/l'^2) sum_n n Q2^n/(1-Q2^n) [1 + cosh(2 pi n a/l')],
      F0 = (pi sin2a/(2 l')) coth(2 eta) - (pi^2 s^2/(4 l'^2)) coth^2(eta) - a sin2a/l' - cos2a/4.
    The corrections are O(e^{-2 pi (pi - a)/l'}). For a <= pi/4 (kappa >= 1/sqrt2) both sums are positive, because
    sinh x >= x and a sin2a - s^2 = s (2 a cos a - sin a) > 0; hence F(pi) < F0 there.
(3) Consequences: the elementary equation F0 = 0 gives l*(kappa) to 1e-4 at kappa = 0.8 and to 1e-11 for kappa >= 0.9;
    since F < F0, every root of F0 is a rigorous lower bound on the true threshold wherever F0 < 0 below it;
    as a -> 0 at fixed arc spacing, F0 -> -eta^2 coth^2(eta)/4 + eta coth(2 eta) - 1/4 (the horocycle criterion, Cor. 20).
Usage: python dual_series.py"""
import mpmath as mp
mp.mp.dps = 40

def Fpi_series(kappa, lp, nmax=None):
    a = mp.acos(kappa); s2 = mp.sin(a)**2; lp = mp.mpf(lp)
    nmax = nmax or int(120/float(lp)) + 10; tot = 0
    for n in range(1, nmax+1):
        r = mp.e**(-n*lp)
        tot += (mp.sin(2*a)*mp.sin(2*n*a) - 2*n*s2*mp.sin(n*a)**2)*r/(1-r) - 2*n*s2*mp.cos(n*a)**2*r/(1+r)
    return tot
def _dlog1(z, q): return mp.jtheta(1, z, q, 1)/mp.jtheta(1, z, q)
def _dlog2(z, q):
    t0, t1, t2 = mp.jtheta(1, z, q), mp.jtheta(1, z, q, 1), mp.jtheta(1, z, q, 2)
    return t2/t0 - (t1/t0)**2
def _E(r): return mp.nsum(lambda n: n*r**n/(1-r**n), [1, mp.inf])
def Fpi_theta(kappa, lp):
    a = mp.acos(kappa); s2 = mp.sin(a)**2; lp = mp.mpf(lp); r = mp.e**(-lp)
    return (mp.sin(2*a)/4*(_dlog1(a, mp.e**(-lp/2)) - mp.cot(a)) + s2/4*(_dlog2(a, r) + 1/mp.sin(a)**2)
            - 2*s2*(_E(r) - _E(r**2)))
def F0(kappa, lp):
    a = mp.acos(kappa); lp = mp.mpf(lp); eta = mp.pi*a/lp
    return (mp.pi*mp.sin(2*a)/(2*lp))*mp.coth(2*eta) - (mp.pi**2*mp.sin(a)**2/(4*lp**2))*mp.coth(eta)**2 \
           - a*mp.sin(2*a)/lp - mp.cos(2*a)/4
def dual_corrections(kappa, lp, nmax=60):
    a = mp.acos(kappa); s2 = mp.sin(a)**2; lp = mp.mpf(lp)
    Q1, Q2 = mp.e**(-4*mp.pi**2/lp), mp.e**(-2*mp.pi**2/lp)
    S1 = mp.fsum(Q1**n/(1-Q1**n)*(mp.sin(2*a)*mp.sinh(4*mp.pi*n*a/lp) - 4*mp.pi*s2*n/lp) for n in range(1, nmax))
    S2 = mp.fsum(n*Q2**n/(1-Q2**n)*(1 + mp.cosh(2*mp.pi*n*a/lp)) for n in range(1, nmax))
    return 2*mp.pi/lp*S1, 2*mp.pi**2*s2/lp**2*S2
def Fpi_dual(kappa, lp):
    c1, c2 = dual_corrections(kappa, lp); return F0(kappa, lp) - c1 - c2

if __name__ == '__main__':
    print('1. three forms of F(pi): Lambert series / theta functions / dual series')
    for kappa, lp in [(0.3, 5.0), (0.72, 3.0), (0.72, 8.0), (0.8, 1.5), (0.95, 0.3)]:
        x, y, z = Fpi_series(kappa, lp), Fpi_theta(kappa, lp), Fpi_dual(kappa, lp)
        print(f'   kappa {kappa:4}, l\' {lp:3}: {mp.nstr(x, 16):>22}  rel diff theta {mp.nstr(abs(y-x)/abs(x), 2):>8}  dual {mp.nstr(abs(z-x)/abs(x), 2):>8}')
    print('2. both correction sums positive for kappa >= 1/sqrt 2 (grid): ', end='')
    ok = all(min(dual_corrections(k, l)) > 0 for k in [0.7072, 0.75, 0.8, 0.9, 0.99, 0.9999] for l in [0.05, 0.5, 2, 6, 20])
    print(ok)
    print('3. arc-length threshold: root of the elementary F0 vs root of the exact series')
    for kappa in [0.72, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999]:
        s = mp.sqrt(1-mp.mpf(kappa)**2)
        guess = {0.72: 4.76, 0.75: 3.19, 0.8: 2.26, 0.85: 1.77, 0.9: 1.44, 0.95: 1.17, 0.99: 0.954, 0.999: 0.874}[kappa]*s
        ex = mp.findroot(lambda l: Fpi_dual(kappa, l), (guess*0.95, guess*1.05), solver='anderson')
        r0 = mp.findroot(lambda l: F0(kappa, l), (ex*0.7, ex*1.01), solver='anderson')
        print(f'   kappa {kappa:6}: l*_F0 = {mp.nstr(r0/s, 12):>15}   l* = {mp.nstr(ex/s, 15):>18}   relative gap {mp.nstr((ex-r0)/ex, 2)}')
