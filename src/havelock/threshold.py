"""threshold(model, N): where the regular N-ring loses stability, which mode goes first, and the basis."""
from __future__ import annotations

import time

import mpmath as mp
import sympy as sp

from .blocks import blocks_mp
from .certify import IsolatedRoot, certify_no_threshold, certify_threshold
from .errors import CertificationFailed, InvalidInput, OutOfValidityRange
from .limits import TANH_EXACT_MAX_N, WORK_DPS
from .models import check_N, get_model
from .models.mt import bn_exact
from .status import NoThreshold, Status, Threshold, make_provenance

_MT_BASIS = ("paper eq. (2) and Proposition 4 (conj13/paper/main.pdf; Montaldi–Tokieda Thm 6, "
             "MT (3.6)): a_m = C1 - m(n-m) with C1 = 2(n-1)(1+b^2)/(1+b)^2 decreasing on (-1, 1], so the ring is "
             "stable for lam_r2 < b_n and the mode floor(n/2) goes first; Lean Havelock.Conj13.bn_spec")
_SINH_BASIS = ("cosh^2 r* = q(N), q = (N-3)/4 (N odd), (N-2)^2/(4(N-1)) (N even) (sinh/note.tex, Section 3; "
               "MT (2.11)); Lean Havelock.Sinh")


def _views_values_relative_to_bn(view: str, bn) -> str:
    """'stable' / 'unstable' if every lam_r2 in the view lies below / above b_n, else 'contains'."""
    if view == "plane":
        return "contains" if bn == 0 else ("stable" if bn > 0 else "unstable")
    if view == "sinh":
        return "contains" if bn < 0 else "stable"
    if view == "sphere":
        return "contains" if bn > 0 else "unstable"
    return "contains"


def _mt_threshold(mdl, N, t0, inputs, strict, precision):
    if N <= 3:
        if mdl.name in ("mt", "sphere"):
            raise OutOfValidityRange(
                f"N = {N}: paper equation (2) has its only root at lam_r2 = 1 (the equator), where the momentum is "
                "fixed by all of SO(3) and the threshold analysis (|b| < 1) does not apply; Montaldi–Tokieda treat this "
                "case separately (MT p. 15)")
        return NoThreshold(status=Status.PROVED, basis=f"C1(b) > floor(N^2/4) for every b < 1 when N = {N} "
                           "(paper Proposition 3)", provenance=make_provenance("threshold", inputs, t0),
                           model=mdl.name, N=N, verdict="stable at every radius",
                           geometry={"plane": "plane", "sinh": "H2"}[mdl.name])
    bn = bn_exact(N)
    where = _views_values_relative_to_bn(mdl.name, bn)
    geometry = {"mt": mdl.geometry(bn), "sphere": "sphere", "plane": "plane", "sinh": "H2"}[mdl.name]
    if where != "contains":
        return NoThreshold(status=Status.PROVED,
                           basis=(f"b_{N} = {str(bn).replace('**', '^')} lies outside the {mdl.name} range "
                                  f"{mdl.domain_str()}; " + _MT_BASIS),
                           provenance=make_provenance("threshold", inputs, t0), model=mdl.name, N=N,
                           verdict=f"{where} at every radius", geometry=geometry)
    with mp.workdps(WORK_DPS):
        value = mp.mpf(str(sp.N(bn, WORK_DPS + 10)))
        r, xi, cosh_r = mdl.radius(value), mdl.xi(value), mdl.cosh_r(value)
    basis = _SINH_BASIS if mdl.name == "sinh" else _MT_BASIS
    side = "below"
    if mdl.name == "plane":
        # a one-point view: there is no stable side inside it
        side = None
        basis = ("degenerate view: the plane is the single value lam_r2 = 0 (no length scale), and every ring with "
                 "N = 7 sits exactly at the bifurcation point b_7 = 0 (spectrally stable, a_3 = 0); the stable side "
                 "lam_r2 < 0 of the mt family lies in the sinh geometry, the unstable side lam_r2 > 0 on the sphere; "
                 + basis)
    return Threshold(status=Status.PROVED, basis=basis, provenance=make_provenance("threshold", inputs, t0),
                     model=mdl.name, N=N, value=value, exact=bn, mode=N // 2, geometry=geometry,
                     stable_side=side, xi=xi, r=r, cosh_r=cosh_r)


def _tanh_numerical(N: int):
    """NUMERICAL fallback: the largest sign change of a_floor(N/2) on (0, 1), located on a grid in
    u = artanh(sqrt(xi)) and refined by bisection in mpmath (no claim about the other modes)."""
    m = N // 2
    with mp.workdps(WORK_DPS):
        f = lambda xi: blocks_mp("tanh", N, m, xi)[0]  # noqa: E731
        us = [mp.mpf(k) / 400 for k in range(1, 4001)]       # geodesic radius r = 2u up to 20
        xs = [mp.tanh(u) ** 2 for u in us]
        vals = [f(x) for x in xs]
        idx = max((i for i in range(len(xs) - 1) if vals[i] < 0 <= vals[i + 1]), default=None)
        if idx is None:
            raise CertificationFailed(f"no sign change of a_{m} found on the grid for N={N}")
        lo, hi = xs[idx], xs[idx + 1]
        for _ in range(200):
            mid = (lo + hi) / 2
            if f(mid) < 0:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2, m


def _tanh_theory(N: int) -> str:
    """What tanh/note_tanh.tex proves about the threshold of the N-ring (N >= 8), with its statement numbers."""
    if N % 2 == 0:
        unique = "note_tanh Proposition 15 (every even N >= 8)"
    elif N <= 15:
        unique = "note_tanh Theorem 7(b) (certified for 8 <= N <= 16)"
    else:
        unique = "note_tanh Proposition 21 (every odd N >= 17)"
    text = (f"the only sign change of a_floor(N/2) in (0, 1), {unique}; critical mode floor(N/2) for every N, "
            "Corollary 14")
    if N % 2 == 0:
        return (text + "; for every even N the energy-momentum threshold equals xi* (Theorem 18(a)), so the ring is "
                "EM-stable for xi > xi*")
    return (text + "; for every odd N >= 9 the Krein window is non-empty, xi_EM > xi* (Theorem 18(b)): just above xi* "
            "the ring is spectrally stable but the Hessian is indefinite, and nonlinear stability there is undecided")


def _tanh_threshold(mdl, N, t0, inputs, strict, precision):
    if N > TANH_EXACT_MAX_N:
        if strict:
            raise CertificationFailed(f"N = {N} exceeds the exact-certification size cap N <= {TANH_EXACT_MAX_N} "
                                      "(hv.limits()); use strict=False for a NUMERICAL value")
        return _tanh_fallback(mdl, N, t0, inputs,
                              f"N = {N} exceeds the exact-certification size cap N <= {TANH_EXACT_MAX_N}")
    try:
        if N <= 7 and certify_no_threshold("tanh", N, mdl.family_domain):
            return NoThreshold(status=Status.CERTIFIED,
                               basis=("exact Sturm counts: a_m, b_m and a_m b_m - g_m^2 have no root in (0, 1) and are "
                                      "positive (note_tanh Theorem 7(a)); Lean-checked: "
                                      f"Havelock.Tanh.stable_N{N}"),
                               provenance=make_provenance("threshold", inputs, t0), model="tanh", N=N,
                               verdict="stable at every radius", geometry="H2")
        c = certify_threshold("tanh", N, mdl.family_domain, "above", start_bits=precision, witness=True)
        if c is None:   # fail closed (any N): no a_m sign change proves no EM claim, and N <= 7 was not certified
            raise CertificationFailed(f"N = {N}: no a_m changes sign in (0, 1), but 'stable at every radius' was "
                                      "not certified (certify_no_threshold); no certified verdict")
    except CertificationFailed as e:
        if strict:
            raise
        return _tanh_fallback(mdl, N, t0, inputs, str(e))
    lo, hi = c.interval
    exact = IsolatedRoot(poly=tuple(int(a) for a in c.factor.coeffs()), lo=lo, hi=hi)
    with mp.workdps(WORK_DPS):
        value = exact.approx(WORK_DPS)
        r, cosh_r = mdl.radius(value), mdl.cosh_r(value)
    lean = f"; Lean-checked: Havelock.Tanh.cert_N{N} (note_tanh Theorem 7(b))" if 8 <= N <= 16 else ""
    # odd N: report the run-time evidence for the Krein window when the certificate has it, i.e. the exact sign of
    # D_crit = a b - g^2 is negative at the rational point just above the interval
    above_signs = c.mode_signs.get(c.mode, {}).get("above", {})
    seen = (f" (seen here at run time: exact D_{c.mode} < 0 at the rational point {float(c.points[1]):.6g} just above "
            "the interval)" if N % 2 and above_signs.get("D") == -1 else "")
    basis = (f"spectral threshold xi* ({_tanh_theory(N)}){seen}. CERTIFIED at run time: Arb root isolation "
             f"({c.bits} bits) inside an exact rational interval of width "
             f"< 2^-{c.bits - 8}, exact Sturm count 1, no B_m and no A_m (m != {c.mode}) vanishing in the interval "
             f"(exact Sturm counts), exact signs of every block on both sides, and no root of any "
             f"a_m, b_m on ({float(hi):.6g}, 1); blocks of tanh/note_tanh.tex Thm 2 (PROVED); critical mode {c.mode}"
             + lean)
    return Threshold(status=Status.CERTIFIED, basis=basis, provenance=make_provenance("threshold", inputs, t0),
                     model="tanh", N=N, value=value, exact=exact, mode=c.mode, geometry="H2", stable_side="above",
                     xi=value, r=r, cosh_r=cosh_r, _certification=c)


def _tanh_fallback(mdl, N, t0, inputs, reason):
    value, m = _tanh_numerical(N)
    with mp.workdps(WORK_DPS):
        r, cosh_r = mdl.radius(value), mdl.cosh_r(value)
    return Threshold(status=Status.NUMERICAL, fallback_from=Status.CERTIFIED,
                     basis=(f"NUMERICAL (strict=False): {reason}. The value is the largest sign change of a_{m} on a "
                            f"grid, refined by bisection; only a_{m} is evaluated at run time. Spectral threshold xi* "
                            f"is {_tanh_theory(N)}" if N >= 8 else
                            f"NUMERICAL (strict=False): {reason}. Largest sign change of a_{m} on a grid, refined by "
                            f"bisection; only a_{m} is evaluated"),
                     provenance=make_provenance("threshold", inputs, t0), model="tanh", N=N, value=value, exact=None,
                     mode=m, geometry="H2", stable_side="above", xi=value, r=r, cosh_r=cosh_r)


_HANDLERS = {"mt": _mt_threshold, "tanh": _tanh_threshold}


def threshold(model: str, N: int, *, strict: bool = True, precision: int | None = None):
    """The stability threshold of the regular N-ring in `model`, as a Threshold, or a NoThreshold when the ring
    is stable (or unstable) at every radius. `precision` is the starting Arb precision in bits for certified
    models; `strict=False` allows a NUMERICAL fallback (recorded in `fallback_from`) instead of raising."""
    t0 = time.perf_counter()
    mdl = get_model(model)
    N = check_N(N)
    if precision is not None and (isinstance(precision, bool) or not isinstance(precision, int) or precision < 64):
        raise InvalidInput(f"precision must be an int >= 64 (bits), got {precision!r}")
    if not isinstance(strict, bool):
        raise InvalidInput("strict must be a bool")
    inputs = {"model": mdl.name, "N": N, "strict": strict, "precision": precision}
    return _HANDLERS[mdl.family](mdl, N, t0, inputs, strict, precision)
