"""Independent exact verifier for threshold certificates (`havelock verify`).

Deliberately self-contained: it rebuilds the block numerators from the papers' closed forms with plain Python
integer lists and re-implements the Sturm count, instead of importing havelock.blocks or havelock.certify. Only
exact integer/rational arithmetic is used (python-flint integer factorisation, fmpq_poly Sturm sequences, Fraction
evaluation); no floating point enters any check.

Every field is type-checked before use: a wrong-typed or malformed certificate yields an INVALID report, never a
traceback.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from fractions import Fraction

import flint

SCHEMA = "havelock-threshold-certificate/1"
_FAMILY = {"mt": "mt", "sphere": "mt", "plane": "mt", "sinh": "mt", "tanh": "tanh"}
_DOMAIN = {"mt": (Fraction(-1), Fraction(1)), "tanh": (Fraction(0), Fraction(1))}
_CLOSED = {"mt": (False, True), "tanh": (False, False)}          # (lo_closed, hi_closed) of the family domain
_STABLE = {"mt": "below", "tanh": "above"}
# Geometry view of each model name: (lo, hi, lo_closed, hi_closed) of lam_r2 (MT family) or xi (tanh). The certified
# root must lie in the view, so a certificate relabelled to a view that does not contain the threshold is INVALID.
_VIEW = {"mt": (Fraction(-1), Fraction(1), False, True), "sphere": (Fraction(0), Fraction(1), False, True),
         "plane": (Fraction(0), Fraction(0), True, True), "sinh": (Fraction(-1), Fraction(0), False, False),
         "tanh": (Fraction(0), Fraction(1), False, False)}
# The next three literals mirror havelock.limits.VERIFY_MAX_N, WITNESS_MAX_PRIME and TANH_EXACT_MAX_N; they are not
# imported, to keep this module independent of the rest of the package (tests enforce the mirror).
_MAX_N = 100000                          # refuse absurd sizes instead of hanging
_MAX_PRIME = 20000                       # largest irreducibility-witness prime accepted
_TANH_MAX_N = 64                         # tanh certificates exist only for N <= 64; larger N would stall
# Each integer in a rational string has at most this many digits: below Python's 4300-digit int/str conversion limit,
# and above the ~2470 digits of an endpoint produced at the top of the precision ladder (4096 bits). Mirrors
# havelock.limits.CERT_MAX_DIGITS (tests enforce the mirror).
_MAX_DIGITS = 2600
_RATIONAL = re.compile(r"-?[0-9]+(/[0-9]+)?")
# The fields of schema havelock-threshold-certificate/1; any other key is rejected.
_REQUIRED = ("schema", "model", "N", "mode", "polynomial", "factor_index", "interval", "endpoint_signs", "sturm",
             "points", "mode_signs", "stable_side", "interval_clear")
_OPTIONAL = ("irreducibility_witness", "provenance")
_STABLE_SIDE_KEYS = {"side", "interval", "root_counts"}
_WITNESS_KEYS = {"p", "degrees"}


# ---- integer polynomial lists (constant term first) ----
def _add(*ps):
    n = max(len(p) for p in ps)
    return [sum(p[i] for p in ps if i < len(p)) for i in range(n)]


def _mul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                out[i + j] += a * b
    return out


def _sc(c, p):
    return [c * a for a in p]


def _mono(k):
    return [0] * k + [1]


def _trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def _mt_numerators(N, m):
    """(A, B, D, G) of the MT family (paper Proposition 3), constant term first."""
    T = m * (N - m)
    A = _add(_sc(2 * (N - 1), [1, 0, 1]), _sc(-T, [1, 2, 1]))
    return _trim(A), [T], _trim(_sc(T, A)), [0]


def _tanh_numerators(N, m):
    """(A, B, D, G) of the tanh model (tanh/note_tanh.tex Thm 2 cleared of denominators), constant term first."""
    X = _mono(N)
    omx, omX, opx = [1, -1], _add([1], _sc(-1, X)), [1, 1]
    omX2, omx2 = _mul(omX, omX), _mul(omx, omx)
    A = _add(_sc(2 * (N - 1), _mul(_mul(opx, omx), omX2)), _sc(-4, _mul(_mono(2), omX2)),
             _sc(-m * (N - m), _mul(omx2, omX2)), _sc(4 * N, _mul(_mul(_mul(opx, omx), X), omX)),
             _sc(-2 * N * N, _mul(X, omx2)),
             _sc(-N, _mul(omx2, _add(_sc(m, _mono(m)), _sc(N - m, _mono(N - m)), _sc(N - m, _mono(N + m)),
                                     _sc(m, _mono(2 * N - m))))))
    Wm = _add(_sc(m, _mul(_mono(m), omX)), _sc(N, _mono(m + N)))
    WNm = _add(_sc(N - m, _mul(_mono(N - m), omX)), _sc(N, _mono(2 * N - m)))
    B = _add(_sc(m * (N - m), omX2), _sc(2 * N * N, X), _sc(-N, _add(Wm, WNm)))
    G = _sc(N, _add(Wm, _sc(-1, WNm)))
    D = _add(_mul(A, B), _sc(-1, _mul(omx2, _mul(G, G))))
    return _trim(A), _trim(B), _trim(D), _trim(G)


_REBUILD = {"mt": _mt_numerators, "tanh": _tanh_numerators}


def _sign(p, q: Fraction) -> int:
    """Exact sign of p(q) for an integer list p (constant term first): for q = a/b with b > 0 this is the sign of
    the homogenised integer sum c_k a^k b^(d-k), computed by Horner's rule in integers (no gcd at each step)."""
    a, b = q.numerator, q.denominator
    r, bp = 0, 1
    for c in reversed(p):
        r = r * a + c * bp
        bp *= b
    return (r > 0) - (r < 0)


def _fz(p):
    return flint.fmpz_poly([int(a) for a in p])


def _strip(p, roots):
    P = _fz(p)
    for r in roots:
        lin = flint.fmpz_poly([-r.numerator, r.denominator])
        while P.degree() > 0 and P(flint.fmpq(r.numerator, r.denominator)) == 0:
            P = P // lin
    return P


def _sturm_sequence(P: flint.fmpz_poly) -> list:
    """Sturm sequence of the square-free part of P (empty for a constant)."""
    if P.degree() <= 0:
        return []
    g = P.gcd(P.derivative())
    if g.degree() > 0:
        P = P // g
    s0 = flint.fmpq_poly(P.coeffs())
    seq = [s0, s0.derivative()]
    while seq[-1].degree() > 0:
        r = -(seq[-2] % seq[-1])
        if r == 0:
            break
        seq.append(r / abs(r.coeffs()[-1]))
    return seq


def _variations(seq, q: Fraction) -> int:
    qq = flint.fmpq(q.numerator, q.denominator)
    s = [v > 0 for v in (t(qq) for t in seq) if v != 0]
    return sum(1 for a, b in zip(s, s[1:]) if a != b)


def _sturm_exact(seq, lo: Fraction, hi: Fraction) -> int:
    return _variations(seq, lo) - _variations(seq, hi) if seq else 0


def _dyadic(q: Fraction, k: int, up: bool) -> Fraction:
    """q rounded down (up=False) or up (up=True) to a multiple of 2^-k."""
    s = q * (1 << k)
    return Fraction(math.ceil(s) if up else math.floor(s), 1 << k)


def _bits(*qs: Fraction) -> int:
    return max(max(q.numerator.bit_length(), q.denominator.bit_length()) for q in qs)


def _sturm_count(seq, lo: Fraction, hi: Fraction) -> int:
    """Number of distinct real roots in (lo, hi] of the polynomial whose Sturm sequence is `seq`, exactly.

    Evaluating a long Sturm sequence at endpoints with thousands of digits is slow, so the count is first
    sandwiched between the counts on dyadic intervals of 2^-k resolution rounded inwards and outwards
    (inner <= exact <= outer, since (lo, hi] contains the inner and lies in the outer half-open interval); when the
    two agree they give the exact count. k doubles from 64; at the endpoints' own size the exact count is used."""
    if not seq:
        return 0
    k, top = 64, _bits(lo, hi)
    while k < top:
        olo, ohi = _dyadic(lo, k, False), _dyadic(hi, k, True)
        ilo, ihi = _dyadic(lo, k, True), _dyadic(hi, k, False)
        outer = _sturm_exact(seq, olo, ohi)
        inner = _sturm_exact(seq, ilo, ihi) if ilo < ihi else 0
        if inner == outer:
            return outer
        k *= 2
    return _sturm_exact(seq, lo, hi)


def _sturm(P: flint.fmpz_poly, lo: Fraction, hi: Fraction, cache: dict | None = None) -> int:
    """distinct real roots of P in (lo, hi]; `cache` maps coefficient tuples to Sturm sequences"""
    key = tuple(int(c) for c in P.coeffs())
    if cache is None:
        seq = _sturm_sequence(P)
    else:
        seq = cache.get(key)
        if seq is None:
            seq = cache[key] = _sturm_sequence(P)
    return _sturm_count(seq, lo, hi)


def _canonical_factors(P: flint.fmpz_poly):
    _, facs = P.factor()
    out = [(-f if f.coeffs()[-1] < 0 else f) for f, _ in facs]
    return sorted(out, key=lambda f: (f.degree(), [int(c) for c in f.coeffs()]))


def _is_prime(n):
    return isinstance(n, int) and not isinstance(n, bool) and 2 <= n <= _MAX_PRIME \
        and all(n % d for d in range(2, int(n**0.5) + 1))


# ---- explicit field validation (a failure raises ValueError and is reported as "malformed certificate") ----
def _int(x, what):
    if not isinstance(x, int) or isinstance(x, bool):
        raise ValueError(f"{what} must be an integer, got {x!r}")
    return x


def _int_list(x, what, length=None):
    if not isinstance(x, list) or (length is not None and len(x) != length) or not x and length is None:
        raise ValueError(f"{what} must be a list of integers" + (f" of length {length}" if length else ""))
    return [_int(a, f"{what} entry") for a in x]


def _strict_eq(a, b):
    """== that also demands equal types throughout (so True, 1.0, [2.0] never match integers)."""
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_strict_eq(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_strict_eq(x, y) for x, y in zip(a, b))
    return type(a) is type(b) and a == b


def _dict(x, what):
    if not isinstance(x, dict):
        raise ValueError(f"{what} must be an object, got {type(x).__name__}")
    return x


def _frac(s, what):
    if not isinstance(s, str) or not _RATIONAL.fullmatch(s):
        raise ValueError(f"{what} must be a rational string 'p' or 'p/q', got {s[:80] if isinstance(s, str) else s!r}")
    if any(len(part.lstrip("-")) > _MAX_DIGITS for part in s.split("/")):
        raise ValueError(f"{what}: a numerator or denominator has more than {_MAX_DIGITS} digits (the verifier's cap "
                         f"on certificate rationals; genuine endpoints have at most about 2470)")
    return Fraction(s)


def _frac_pair(x, what):
    if not isinstance(x, list) or len(x) != 2:
        raise ValueError(f"{what} must be a list of two rational strings")
    return _frac(x[0], what), _frac(x[1], what)


@dataclass
class VerifyReport:
    ok: bool
    failures: list = field(default_factory=list)
    checks: int = 0

    def __bool__(self):
        return self.ok


def verify_dict(d: dict) -> VerifyReport:
    """Re-check every mathematical field of a certificate dict (and that the model's view contains the certified
    root); provenance is not checked. Returns a report listing every failed check."""
    fails, n = [], 0

    def check(cond, msg):
        nonlocal n
        n += 1
        if not cond:
            fails.append(msg)
        return cond

    try:
        if not check(isinstance(d, dict), "certificate must be a JSON object"):
            return VerifyReport(False, fails, n)
        if not check(d.get("schema") == SCHEMA, f"schema must be {SCHEMA!r}"):
            return VerifyReport(False, fails, n)
        unknown = sorted(str(k) for k in d if k not in _REQUIRED and k not in _OPTIONAL)
        if not check(not unknown, f"unknown field(s) {unknown[:5]} (not in schema {SCHEMA})"):
            return VerifyReport(False, fails, n)
        model = d["model"]
        if not check(isinstance(model, str) and model in _FAMILY, f"unknown model {model!r}"):
            return VerifyReport(False, fails, n)
        N, mode = d["N"], d["mode"]
        if not check(isinstance(N, int) and not isinstance(N, bool) and isinstance(mode, int)
                     and not isinstance(mode, bool) and 2 <= N <= _MAX_N and 1 <= mode <= N // 2,
                     f"N and mode must be integers with 2 <= N <= {_MAX_N} and 1 <= mode <= N//2"):
            return VerifyReport(False, fails, n)
        fam = _FAMILY[model]
        if not check(fam != "tanh" or N <= _TANH_MAX_N,
                     f"tanh certificates are limited to N <= {_TANH_MAX_N} (the exact-certification size cap)"):
            return VerifyReport(False, fails, n)
        dlo, dhi = _DOMAIN[fam]
        lo_closed, hi_closed = _CLOSED[fam]
        cache = {}                       # Sturm sequences, shared by the interval and stable-side counts
        nums = {m: _REBUILD[fam](N, m) for m in range(1, N // 2 + 1)}
        poly = _int_list(d["polynomial"], "polynomial")
        if not check(_trim(poly) == nums[mode][0], "polynomial is not the numerator A_mode of the critical block"):
            return VerifyReport(False, fails, n)
        lo, hi = _frac_pair(d["interval"], "interval")
        check(dlo < lo < hi < dhi, "interval must satisfy domain_lo < lo < hi < domain_hi")
        # the certified root (the unique root of `polynomial` in (lo, hi), checked below) must lie in the view
        vlo, vhi, vlc, vhc = _VIEW[model]
        if vlo == vhi:     # a one-point view (the plane, lam_r2 = 0): the root must be exactly that point
            check(lo < vlo < hi and _sign(poly, vlo) == 0,
                  f"the threshold is not in the {model!r} view {{{vlo}}}: need lo < {vlo} < hi and "
                  f"polynomial({vlo}) = 0")
        else:              # vlo <= lo < root < hi <= vhi: the root is strictly inside the view
            check(vlo <= lo and hi <= vhi,
                  f"the threshold interval ({lo}, {hi}) is not inside the {model!r} view "
                  f"{'[' if vlc else '('}{vlo}, {vhi}{']' if vhc else ')'}")
        factors = _canonical_factors(_fz(poly))
        k = _int(d["factor_index"], "factor_index")
        if check(0 <= k < len(factors), "factor_index out of range"):
            f = factors[k]
            fl, fh = f(flint.fmpq(lo.numerator, lo.denominator)), f(flint.fmpq(hi.numerator, hi.denominator))
            check(fl * fh < 0, "the selected irreducible factor does not change sign on the interval")
        es = [_sign(poly, lo), _sign(poly, hi)]
        check(_int_list(d["endpoint_signs"], "endpoint_signs", 2) == es and es[0] * es[1] < 0,
              "endpoint_signs do not match, or the polynomial does not change sign")
        check(_int(d["sturm"], "sturm") == 1 and _sturm(_fz(poly), lo, hi, cache) == 1, "Sturm count on the interval is not 1")

        # no other block vanishes in (lo, hi]: B_m for every m, A_m for m != mode
        expected_clear = {}
        for m, (A, B, _, _) in nums.items():
            expected_clear[f"B_{m}"] = _sturm(_strip(B, (dlo, dhi)), lo, hi, cache)
            if m != mode:
                expected_clear[f"A_{m}"] = _sturm(_strip(A, (dlo, dhi)), lo, hi, cache)
        got_clear = d.get("interval_clear")
        check(isinstance(got_clear, dict) and all(isinstance(v, int) and not isinstance(v, bool)
                                                  for v in got_clear.values()) and got_clear == expected_clear
              and all(v == 0 for v in expected_clear.values()),
              f"interval_clear: recorded {got_clear!r}, recomputed {expected_clear!r}; every count must be 0 "
              "(no other block may vanish in the threshold interval)")

        below, above = _frac_pair(d["points"], "points")
        check(dlo < below < lo and hi < above < dhi, "sample points must lie just outside the interval")
        ms = _dict(d["mode_signs"], "mode_signs")
        check(sorted(ms) == sorted(str(m) for m in nums), "mode_signs must list every mode 1..N//2 (keys \"1\", \"2\", ...)")
        stable_pt = "above" if _STABLE[fam] == "above" else "below"
        unstable_pt = "below" if stable_pt == "above" else "above"
        coupled = any(any(c != 0 for c in G) for (_, _, _, G) in nums.values())
        for m, (A, B, D, _) in nums.items():
            rec = _dict(ms.get(str(m), {}), f"mode_signs[{m}]")
            check(set(rec) <= {"below", "above"}, f"mode_signs[{m}] has unknown fields {sorted(set(rec) - {'below', 'above'})}")
            for where, q in (("below", below), ("above", above)):
                exp = {"A": _sign(A, q), "B": _sign(B, q)}
                if coupled:
                    exp["D"] = _sign(D, q)
                check(_strict_eq(rec.get(where), exp), f"mode_signs[{m}][{where}] = {rec.get(where)} != recomputed {exp}")
                if where == stable_pt:
                    check(exp["A"] > 0 and exp["B"] > 0, f"mode {m} is not stable at the {where} point")
        uA, uB = _sign(nums[mode][0], {"below": below, "above": above}[unstable_pt]), \
            _sign(nums[mode][1], {"below": below, "above": above}[unstable_pt])
        check(uA * uB < 0, "the critical mode is not unstable at the unstable-side point")
        side = _dict(d["stable_side"], "stable_side")
        check(set(side) == _STABLE_SIDE_KEYS, f"stable_side must have exactly the fields {sorted(_STABLE_SIDE_KEYS)}")
        exp_iv = [str(hi), str(dhi)] if _STABLE[fam] == "above" else [str(dlo), str(lo)]
        check(side.get("side") == _STABLE[fam] and [str(x) for x in _frac_pair(side["interval"], "stable_side.interval")]
              == exp_iv, "stable_side must be the open interval between the threshold interval and the domain end")
        s_lo, s_hi = (Fraction(s) for s in exp_iv)
        rc = _dict(side["root_counts"], "stable_side.root_counts")
        check(sorted(rc) == sorted(str(m) for m in nums), "stable_side.root_counts must list exactly the modes 1..N//2")
        for m, (A, B, _, _) in nums.items():
            got = {"A": _sturm(_strip(A, (dlo, dhi)), s_lo, s_hi, cache),
                   "B": _sturm(_strip(B, (dlo, dhi)), s_lo, s_hi, cache)}
            check(_strict_eq(rc.get(str(m)), got) and got == {"A": 0, "B": 0},
                  f"stable-side root counts for mode {m}: recorded {rc.get(str(m))}, got {got}")

        # a closed domain endpoint on the stable side must itself be stable (A_m, B_m > 0, exact signs)
        end = dhi if _STABLE[fam] == "above" and hi_closed else dlo if _STABLE[fam] == "below" and lo_closed else None
        if end is not None:
            bad = [m for m, (A, B, _, _) in nums.items() if _sign(A, end) <= 0 or _sign(B, end) <= 0]
            check(not bad, f"at the closed domain endpoint {end}, a block is not positive (modes {bad[:5]})")

        w = d.get("irreducibility_witness")
        if w is not None and check(0 <= k < len(factors), "witness needs a valid factor"):
            w = _dict(w, "irreducibility_witness")
            check(set(w) == _WITNESS_KEYS, f"irreducibility_witness must have exactly the fields {sorted(_WITNESS_KEYS)}")
            f = factors[k]
            p = w.get("p")
            ok_p = check(_is_prime(p) and int(f.coeffs()[-1]) % p != 0, "witness prime must be a prime "
                         f"<= {_MAX_PRIME} not dividing lead(f)")
            if ok_p:
                _, fs = flint.nmod_poly([int(c) % p for c in f.coeffs()], p).factor()
                degs = sorted(g.degree() for g, e in fs for _ in range(e))
                check(_strict_eq(w.get("degrees"), degs) and degs == [f.degree()], "f is not irreducible modulo the witness prime")
    except (KeyError, TypeError, ValueError, ZeroDivisionError, AttributeError, IndexError, OverflowError) as e:
        fails.append(f"malformed certificate: {e!r}")
    return VerifyReport(not fails, fails, n)
