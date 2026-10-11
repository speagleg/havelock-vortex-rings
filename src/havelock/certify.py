"""Exact certification with python-flint: factorisation over Q, Arb root isolation, exact Sturm counts, and
mode-sign checks. No floating-point comparison enters any pass/fail decision: Arb balls are converted to exact
dyadic rationals and every sign is evaluated exactly at a rational point.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

import flint
import mpmath as mp

from .blocks import numerators
from .errors import CertificationFailed
from .limits import PRECISION_LADDER, WITNESS_MAX_PRIME


class _Unresolved(Exception):
    """The current Arb precision could not separate the root balls; try the next rung of the ladder."""


def fmpq_of(q: Fraction) -> flint.fmpq:
    return flint.fmpq(q.numerator, q.denominator)


def sign_at(p: flint.fmpz_poly, q: Fraction) -> int:
    v = p(fmpq_of(q))
    return (v > 0) - (v < 0)


def squarefree(p: flint.fmpz_poly) -> flint.fmpz_poly:
    g = p.gcd(p.derivative())
    return p // g if g.degree() > 0 else p


def strip_roots(p: flint.fmpz_poly, roots) -> flint.fmpz_poly:
    """Divide out every factor (q x - p) for the rational roots p/q listed, as often as it divides."""
    for r in roots:
        r = Fraction(r)
        lin = flint.fmpz_poly([-r.numerator, r.denominator])
        while not p.is_zero() and p.degree() > 0 and p(fmpq_of(r)) == 0:
            p = p // lin
    return p


def sturm_count(p: flint.fmpz_poly, lo: Fraction, hi: Fraction) -> int:
    """Number of distinct real roots of p in the half-open interval (lo, hi] (exact Sturm sequence)."""
    if p.degree() <= 0:
        return 0
    s0 = flint.fmpq_poly(squarefree(p).coeffs())
    seq = [s0, s0.derivative()]
    while seq[-1].degree() > 0:
        r = -(seq[-2] % seq[-1])
        if r.is_zero():
            break
        seq.append(r / abs(r.coeffs()[-1]))

    def changes(q: Fraction) -> int:
        qq = fmpq_of(q)
        signs = [v > 0 for v in (s(qq) for s in seq) if v != 0]
        return sum(1 for a, b in zip(signs, signs[1:]) if a != b)

    return changes(lo) - changes(hi)


def _ball(z) -> tuple[Fraction, Fraction]:
    m, e = z.mid().man_exp()
    r, f = z.rad().man_exp()
    mid = Fraction(int(m)) * Fraction(2) ** int(e)
    rad = Fraction(int(r)) * Fraction(2) ** int(f)
    return mid - rad, mid + rad


def real_root_balls(p: flint.fmpz_poly, lo: Fraction, hi: Fraction, bits: int):
    """Exact rational enclosures [L, U] of the real roots of p inside (lo, hi), from Arb (selected by the exact
    test imag == 0). Raises _Unresolved if a real ball meets lo or hi."""
    if p.degree() <= 0:
        return []
    old = flint.ctx.prec
    flint.ctx.prec = bits
    try:
        roots = squarefree(p).complex_roots()
    finally:
        flint.ctx.prec = old
    out = []
    for z, _ in roots:
        if not z.imag.is_zero():
            continue
        L, U = _ball(z.real)
        if lo < L and U < hi:
            out.append((L, U))
        elif not (U < lo or L > hi):
            raise _Unresolved
    return sorted(out)


def canonical_factors(p: flint.fmpz_poly) -> list[flint.fmpz_poly]:
    """Distinct irreducible factors of p over Z: primitive, positive leading coefficient, sorted by
    (degree, coefficients). This order defines `factor_index` in certificates."""
    _, facs = p.factor()
    out = []
    for f, _ in facs:
        if f.coeffs()[-1] < 0:
            f = -f
        out.append(f)
    return sorted(out, key=lambda f: (f.degree(), [int(c) for c in f.coeffs()]))


def _is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def mod_p_witness(f: flint.fmpz_poly, max_prime: int = WITNESS_MAX_PRIME):
    """Smallest prime p not dividing lead(f) with f irreducible mod p (this proves f irreducible over Q), as
    {"p": p, "degrees": [deg f]}; None if no such prime <= max_prime."""
    lead = int(f.coeffs()[-1])
    n = f.degree()
    for p in range(2, max_prime + 1):
        if not _is_prime(p) or lead % p == 0:
            continue
        _, fs = flint.nmod_poly([int(c) % p for c in f.coeffs()], p).factor()
        degs = sorted(g.degree() for g, e in fs for _ in range(e))
        if degs == [n]:
            return {"p": p, "degrees": degs}
    return None


@dataclass(frozen=True)
class IsolatedRoot:
    """A real algebraic number: the root of the integer polynomial `poly` (constant term first) in [lo, hi]."""

    poly: tuple
    lo: Fraction
    hi: Fraction

    def __str__(self) -> str:
        return (f"the root in [{float(self.lo):.17g}, {float(self.hi):.17g}] (width {float(self.hi - self.lo):.1e}) "
                f"of the degree-{len(self.poly) - 1} polynomial {list(self.poly)} (constant term first)")

    def approx(self, dps: int = 50):
        with mp.workdps(dps + 5):
            return (mp.mpf(self.lo.numerator) / self.lo.denominator + mp.mpf(self.hi.numerator) / self.hi.denominator) / 2


@dataclass(frozen=True)
class Certification:
    """Everything a threshold certificate records (see havelock.certificates for the JSON form)."""

    family: str
    N: int
    mode: int
    polynomial: flint.fmpz_poly
    factor: flint.fmpz_poly
    factor_index: int
    interval: tuple
    endpoint_signs: tuple
    sturm: int
    points: tuple
    mode_signs: dict
    stable_side: dict
    bits: int
    irreducibility_witness: dict | None = field(default=None)


def _domain_points(dom):
    lo, hi, _, _ = dom
    return Fraction(lo), Fraction(hi)


def coupled(family: str, N: int) -> bool:
    """True if some g_m (1 <= m <= N//2) is not identically zero; then certificates also record a_m b_m - g_m^2."""
    return any(not numerators(family, N, m).G.is_zero() for m in range(1, N // 2 + 1))


def _signs(family: str, N: int, m: int, q: Fraction, with_d: bool) -> dict:
    nm = numerators(family, N, m)
    s = {"A": sign_at(nm.A, q), "B": sign_at(nm.B, q)}
    if with_d:
        s["D"] = sign_at(nm.D, q)
    return s


def certify_threshold(family: str, N: int, domain, stable_side: str, start_bits: int | None = None,
                      witness: bool = False) -> Certification | None:
    """Certify the stability threshold of the N-ring, or return None if no a_m changes sign in the domain.

    Steps: (1) Arb-isolate every real root of every A_m (m = 1..N//2) inside the domain; the threshold is the
    largest root (stable_side "above") or the smallest (stable_side "below"), and its ball must be disjoint from
    all others; (2) select the irreducible factor of A_crit carrying it; (3) exact checks: A_crit changes sign
    across the interval and has exactly one root in it (Sturm); (4) signs of A_m, B_m (and D_m when g != 0) at
    the rational points 2 lo - hi and 2 hi - lo; (5) on the stable side of the interval no A_m or B_m has a root
    (Sturm), and all of them are positive there. Raises CertificationFailed if a check fails."""
    dlo, dhi = _domain_points(domain)
    modes = range(1, N // 2 + 1)
    ladder = [b for b in PRECISION_LADDER if start_bits is None or b >= start_bits] or [start_bits]
    for bits in ladder:
        try:
            balls = []
            for m in modes:
                A = strip_roots(numerators(family, N, m).A, (dlo, dhi))
                balls += [(L, U, m) for (L, U) in real_root_balls(A, dlo, dhi, bits)]
            if not balls:
                return None
            pick = max(balls, key=lambda t: t[1]) if stable_side == "above" else min(balls, key=lambda t: t[0])
            L, U, mode = pick
            if L == U:  # Arb returns an exact rational root as a point ball; widen it to an interval
                L, U = L - Fraction(1, 2**bits), U + Fraction(1, 2**bits)
            for (L2, U2, m2) in balls:
                if (L2, U2, m2) != pick and not (U2 < L or L2 > U):
                    raise _Unresolved
            break
        except _Unresolved:
            continue
    else:
        raise CertificationFailed(f"Arb root isolation did not separate the roots for N={N} "
                                  f"(ladder {list(PRECISION_LADDER)} bits)")
    P = numerators(family, N, mode).A
    factors = canonical_factors(P)
    carrying = [i for i, f in enumerate(factors) if sign_at(f, L) * sign_at(f, U) < 0]
    if len(carrying) != 1:
        raise CertificationFailed("could not identify the irreducible factor carrying the threshold")
    idx = carrying[0]
    es = (sign_at(P, L), sign_at(P, U))
    if es[0] * es[1] >= 0:
        raise CertificationFailed("A_crit does not change sign across the isolating interval")
    st = sturm_count(P, L, U)
    if st != 1:
        raise CertificationFailed(f"Sturm count of A_crit in the interval is {st}, expected 1")
    w = U - L
    below, above = L - w, U + w
    if not (dlo < below and above < dhi):
        raise CertificationFailed("sample points fall outside the domain")
    cd = coupled(family, N)
    mode_signs = {m: {"below": _signs(family, N, m, below, cd), "above": _signs(family, N, m, above, cd)}
                  for m in modes}
    side_lo, side_hi = (U, dhi) if stable_side == "above" else (dlo, L)
    counts = {}
    for m in modes:
        nm = numerators(family, N, m)
        counts[m] = {"A": sturm_count(strip_roots(nm.A, (dlo, dhi)), side_lo, side_hi),
                     "B": sturm_count(strip_roots(nm.B, (dlo, dhi)), side_lo, side_hi)}
    stable_pt, unstable_pt = ("above", "below") if stable_side == "above" else ("below", "above")
    ok_stable = all(mode_signs[m][stable_pt]["A"] > 0 and mode_signs[m][stable_pt]["B"] > 0 for m in modes)
    ok_unstable = mode_signs[mode][unstable_pt]["A"] * mode_signs[mode][unstable_pt]["B"] < 0
    ok_counts = all(c["A"] == 0 and c["B"] == 0 for c in counts.values())
    if not (ok_stable and ok_unstable and ok_counts):
        raise CertificationFailed("the sign pattern on the two sides of the threshold is not the expected one "
                                  "(stable side: every a_m, b_m > 0; unstable side: a_crit b_crit < 0)")
    f = factors[idx]
    return Certification(family=family, N=N, mode=mode, polynomial=P, factor=f, factor_index=idx,
                         interval=(L, U), endpoint_signs=es, sturm=st, points=(below, above),
                         mode_signs=mode_signs,
                         stable_side={"side": stable_side, "interval": (side_lo, side_hi), "root_counts": counts},
                         bits=bits, irreducibility_witness=mod_p_witness(f) if witness else None)


def certify_no_threshold(family: str, N: int, domain) -> bool:
    """True iff every A_m, B_m (and D_m when g != 0), m = 1..N//2, has no root inside the domain and is positive
    at its midpoint (exact Sturm counts): the ring is stable, with positive-definite blocks, at every radius."""
    dlo, dhi = _domain_points(domain)
    mid = (dlo + dhi) / 2
    for m in range(1, N // 2 + 1):
        nm = numerators(family, N, m)
        polys = [nm.A, nm.B] + ([nm.D] if not nm.G.is_zero() else [])
        for p in polys:
            q = strip_roots(p, (dlo, dhi))
            if sturm_count(q, dlo, dhi) != 0 or sign_at(p, mid) <= 0:
                return False
    return True
