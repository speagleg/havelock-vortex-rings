"""Driver: route B (kkt_exact.py) for given n, in mpmath (root b_n in (-1,1)) and exactly in GF(p)
for several primes and both conjugates of lam.  Compares with
  odd n:  (n-1)(n-3)(n^2-9n+38)/(128 n)    even n = 8: 95/56."""
import sys, time, mpmath as mp
from fractions import Fraction
from sympy import isprime, sqrt_mod, primitive_root
from kkt_exact import run, Fp
mp.mp.dps = 40

def formula(n):
    if n % 2: return Fraction((n - 1) * (n - 3) * (n * n - 9 * n + 38), 128 * n)
    return {8: Fraction(95, 56)}[n]

def lam_roots(n, sq):     # sq = sqrt(n-3) (odd) or sqrt(7) (n=8)
    if n == 7: return [0 * sq]
    if n == 8: return [-8 + 3 * sq, -8 - 3 * sq]
    return [(-(n + 1) + 4 * sq) / (n - 7), (-(n + 1) - 4 * sq) / (n - 7)]

def float_run(n):
    l = n // 2; ph = [2 * mp.pi * k / n for k in range(n)]
    sq = mp.sqrt(n - 3 if n % 2 else 7)
    out = []
    for lam in lam_roots(n, sq):
        if not (-1 < lam <= 1): continue
        r = run(n, lam, [mp.cos(x) for x in ph], [mp.sin(x) for x in ph],
                [mp.cos(l * x) for x in ph], [mp.sin(l * x) for x in ph], mp.mpf(1))
        out.append((lam, r))
    return out

def modp_run(n, p, sign):
    Fp.p = p
    g = primitive_root(p); z4 = Fp(pow(g, (p - 1) // (4 * n), p))
    zn = z4 * z4 * z4 * z4; I = Fp(pow(z4.v, n, p)); assert (I * I).v == p - 1
    l = n // 2
    def pw(a, k): return Fp(pow(a.v, k % (p - 1), p))
    cs = [(pw(zn, k) + pw(zn, -k)) / 2 for k in range(n)]
    sn = [(pw(zn, k) - pw(zn, -k)) / (2 * I) for k in range(n)]
    cl = [(pw(zn, l * k) + pw(zn, -l * k)) / 2 for k in range(n)]
    sl = [(pw(zn, l * k) - pw(zn, -l * k)) / (2 * I) for k in range(n)]
    d = n - 3 if n % 2 else 7
    sq = Fp(sqrt_mod(d, p)); assert (sq * sq).v == d % p
    lam = lam_roots(n, sign * sq)[0]
    return run(n, lam, cs, sn, cl, sl, Fp(1))

def primes_for(n, count, start=2**61):
    out = []; k = start // (4 * n)
    while len(out) < count:
        k += 1; p = 4 * n * k + 1
        if isprime(p) and pow(n - 3 if n % 2 else 7, (p - 1) // 2, p) == 1: out.append(p)
    return out

if __name__ == '__main__':
    ns = [int(a) for a in sys.argv[1:]] or [5, 7, 9, 11, 13, 15, 17, 19]
    for n in ns:
        T = time.time(); F = formula(n)
        for lam, r in float_run(n):
            print('n=%d float: b_n=%s nu_e1=%.2e nu_e2=%.2e c4=%s formula=%s diff=%.2e' % (
                n, mp.nstr(lam, 15), float(abs(r['nue1'])), float(abs(r['nue2'])), mp.nstr(r['c4'], 30),
                mp.nstr(mp.mpf(F.numerator) / F.denominator, 30), float(r['c4'] - mp.mpf(F.numerator) / F.denominator)))
        ok = True; vals = []
        for p in primes_for(n, 3):
            for sign in ((1,) if n == 7 else (1, -1)):
                r = modp_run(n, p, sign)
                Ffp = Fp(F.numerator) / F.denominator
                good = r['nue1'].iszero() and r['nue2'].iszero() and (r['c4'] - Ffp).iszero()
                ok &= good; vals.append(good)
        print('n=%d EXACT GF(p): 3 primes x %d conjugate(s) of lam: nu_e1=nu_e2=0 and c4 == %s : %s  [%.0fs]' % (
            n, 1 if n == 7 else 2, F, ok, time.time() - T), flush=True)
