"""Independent exact verifier for threshold certificates (`havelock verify`).

Deliberately self-contained: it rebuilds the block numerators from the papers' closed forms with plain Python
integer lists and re-implements the Sturm count, instead of importing havelock.blocks or havelock.certify. Only
exact integer/rational arithmetic is used (python-flint integer factorisation, fmpq_poly Sturm sequences, Fraction
evaluation); no floating point enters any check.

Every field is type-checked before use: a wrong-typed or malformed certificate yields an INVALID report, never a
traceback.
"""
from __future__ import annotations

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
_RATIONAL = re.compile(r"-?[0-9]{1,5000}(/[0-9]{1,5000})?")


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


def _ev(p, q: Fraction) -> Fraction:
    r = Fraction(0)
    for a in reversed(p):
        r = r * q + a
    return r


def _sign(p, q):
    v = _ev(p, q)
    return (v > 0) - (v < 0)


def _fz(p):
    return flint.fmpz_poly([int(a) for a in p])


def _strip(p, roots):
    P = _fz(p)
    for r in roots:
        lin = flint.fmpz_poly([-r.numerator, r.denominator])
        while P.degree() > 0 and P(flint.fmpq(r.numerator, r.denominator)) == 0:
            P = P // lin
    return P


def _sturm(P: flint.fmpz_poly, lo: Fraction, hi: Fraction) -> int:
    """distinct real roots of P in (lo, hi]"""
    if P.degree() <= 0:
        return 0
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

    def V(q):
        qq = flint.fmpq(q.numerator, q.denominator)
        s = [v > 0 for v in (t(qq) for t in seq) if v != 0]
        return sum(1 for a, b in zip(s, s[1:]) if a != b)

    return V(lo) - V(hi)


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
        raise ValueError(f"{what} must be a rational string 'p' or 'p/q', got {s!r}")
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
        nums = {m: _REBUILD[fam](N, m) for m in range(1, N // 2 + 1)}
        poly = _int_list(d["polynomial"], "polynomial")
        if not check(_trim(poly) == nums[mode][0], "polynomial is not the numerator A_mode of the critical block"):
            return VerifyReport(False, fails, n)
        lo, hi = _frac_pair(d["interval"], "interval")
        check(dlo < lo < hi < dhi, "interval must satisfy domain_lo < lo < hi < domain_hi")
        # F1: the certified root (the unique root of `polynomial` in (lo, hi), checked below) must lie in the view
        vlo, vhi, vlc, vhc = _VIEW[model]
        if vlo == vhi:     # a one-point view (the plane, lam_r2 = 0): the root must be exactly that point
            check(lo < vlo < hi and _ev(poly, vlo) == 0,
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
        check(_int(d["sturm"], "sturm") == 1 and _sturm(_fz(poly), lo, hi) == 1, "Sturm count on the interval is not 1")

        # T4-I1: no other block vanishes in (lo, hi]: B_m for every m, A_m for m != mode
        expected_clear = {}
        for m, (A, B, _, _) in nums.items():
            expected_clear[f"B_{m}"] = _sturm(_strip(B, (dlo, dhi)), lo, hi)
            if m != mode:
                expected_clear[f"A_{m}"] = _sturm(_strip(A, (dlo, dhi)), lo, hi)
        got_clear = d.get("interval_clear")
        check(isinstance(got_clear, dict) and all(isinstance(v, int) and not isinstance(v, bool)
                                                  for v in got_clear.values()) and got_clear == expected_clear
              and all(v == 0 for v in expected_clear.values()),
              f"interval_clear: recorded {got_clear!r}, recomputed {expected_clear!r}; every count must be 0 "
              "(no other block may vanish in the threshold interval)")

        below, above = _frac_pair(d["points"], "points")
        check(dlo < below < lo and hi < above < dhi, "sample points must lie just outside the interval")
        ms = _dict(d["mode_signs"], "mode_signs")
        check(sorted(int(m) for m in ms) == list(nums), "mode_signs must list every mode 1..N//2")
        stable_pt = "above" if _STABLE[fam] == "above" else "below"
        unstable_pt = "below" if stable_pt == "above" else "above"
        coupled = any(any(c != 0 for c in G) for (_, _, _, G) in nums.values())
        for m, (A, B, D, _) in nums.items():
            rec = _dict(ms.get(str(m), {}), f"mode_signs[{m}]")
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
        exp_iv = [str(hi), str(dhi)] if _STABLE[fam] == "above" else [str(dlo), str(lo)]
        check(side.get("side") == _STABLE[fam] and [str(x) for x in _frac_pair(side["interval"], "stable_side.interval")]
              == exp_iv, "stable_side must be the open interval between the threshold interval and the domain end")
        s_lo, s_hi = (Fraction(s) for s in exp_iv)
        rc = _dict(side["root_counts"], "stable_side.root_counts")
        for m, (A, B, _, _) in nums.items():
            got = {"A": _sturm(_strip(A, (dlo, dhi)), s_lo, s_hi), "B": _sturm(_strip(B, (dlo, dhi)), s_lo, s_hi)}
            check(_strict_eq(rc.get(str(m)), got) and got == {"A": 0, "B": 0},
                  f"stable-side root counts for mode {m}: recorded {rc.get(str(m))}, got {got}")

        # T4-I2: a closed domain endpoint on the stable side must itself be stable (A_m, B_m > 0, exact signs)
        end = dhi if _STABLE[fam] == "above" and hi_closed else dlo if _STABLE[fam] == "below" and lo_closed else None
        if end is not None:
            bad = [m for m, (A, B, _, _) in nums.items() if _sign(A, end) <= 0 or _sign(B, end) <= 0]
            check(not bad, f"at the closed domain endpoint {end}, a block is not positive (modes {bad[:5]})")

        w = d.get("irreducibility_witness")
        if w is not None and check(0 <= k < len(factors), "witness needs a valid factor"):
            w = _dict(w, "irreducibility_witness")
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
