"""Table check for the paper: runs the third reader's Cartesian order-by-order KKT engine
(../review3/kkt_exact.py, no Fourier ansatz, no symmetry reduction) exactly in GF(p) (3 primes x both
conjugates of lam) and in 40-digit floats, for n = 4..25, and compares c4 with the closed forms
  even n: (n-2)(m^4+40m^2+48m+16)/(192 n (n-1)), m = n-2;   odd n: (n-1)(n-3)(n^2-9n+38)/(128 n).
Even-n threshold: lam = (n^2 +- 4(n-2) sqrt(n-1)) / (8(n-1) - n^2)."""
import sys, os, time, mpmath as mp
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'review3'))
from sympy import isprime, sqrt_mod, primitive_root
from kkt_exact import run, Fp
mp.mp.dps = 40

def formula(n):
    if n % 2: return Fraction((n-1)*(n-3)*(n*n-9*n+38), 128*n)
    m = n-2; return Fraction((n-2)*(m**4+40*m*m+48*m+16), 192*n*(n-1))
def disc(n): return n-3 if n % 2 else n-1
def lam_roots(n, sq):
    if n == 7: return [0*sq]
    if n % 2: return [(-(n+1)+4*sq)/(n-7), (-(n+1)-4*sq)/(n-7)]
    A = 8*(n-1)-n*n
    return [(n*n+4*(n-2)*sq)/A, (n*n-4*(n-2)*sq)/A]
def float_run(n):
    l = n//2; ph = [2*mp.pi*k/n for k in range(n)]; out = []
    for lam in lam_roots(n, mp.sqrt(disc(n))):
        if -1 < lam <= 1:
            out.append((lam, run(n, lam, [mp.cos(x) for x in ph], [mp.sin(x) for x in ph],
                                 [mp.cos(l*x) for x in ph], [mp.sin(l*x) for x in ph], mp.mpf(1))))
    return out
def modp_run(n, p, sign):
    Fp.p = p; g = primitive_root(p); z4 = Fp(pow(g, (p-1)//(4*n), p))
    zn = z4*z4*z4*z4; I = Fp(pow(z4.v, n, p)); assert (I*I).v == p-1
    l = n//2
    def pw(a, k): return Fp(pow(a.v, k % (p-1), p))
    cs = [(pw(zn, k)+pw(zn, -k))/2 for k in range(n)]; sn = [(pw(zn, k)-pw(zn, -k))/(2*I) for k in range(n)]
    cl = [(pw(zn, l*k)+pw(zn, -l*k))/2 for k in range(n)]; sl = [(pw(zn, l*k)-pw(zn, -l*k))/(2*I) for k in range(n)]
    sq = Fp(sqrt_mod(disc(n), p)); assert (sq*sq).v == disc(n) % p
    return run(n, lam_roots(n, sign*sq)[0], cs, sn, cl, sl, Fp(1))
def primes_for(n, count, start=2**61):
    out = []; k = start//(4*n)
    while len(out) < count:
        k += 1; p = 4*n*k+1
        if isprime(p) and pow(disc(n), (p-1)//2, p) == 1: out.append(p)
    return out
if __name__ == '__main__':
    ns = [int(a) for a in sys.argv[1:]] or list(range(4, 26))
    for n in ns:
        T = time.time(); F = formula(n)
        for lam, r in float_run(n):
            print('n=%d b_n=%s c4(float)=%s formula=%s=%s diff=%.1e nu1=%.1e' % (n, mp.nstr(lam, 12),
                  mp.nstr(r['c4'], 25), F, mp.nstr(mp.mpf(F.numerator)/F.denominator, 25),
                  float(r['c4']-mp.mpf(F.numerator)/F.denominator), float(abs(r['nue1']))))
        ok = True
        for p in primes_for(n, 3):
            for sign in ((1,) if n == 7 else (1, -1)):
                r = modp_run(n, p, sign); Ffp = Fp(F.numerator)/F.denominator
                ok &= r['nue1'].iszero() and r['nue2'].iszero() and (r['c4']-Ffp).iszero()
        print('n=%d EXACT GF(p) (3 primes x conjugates): threshold, no cubic, c4 == %s : %s [%.0fs]' % (n, F, ok, time.time()-T), flush=True)
