"""threshold(model, N): where the regular N-ring loses stability, which mode goes first, and the basis."""
from __future__ import annotations

import time

import mpmath as mp
import sympy as sp

from .errors import InvalidInput, OutOfValidityRange
from .limits import WORK_DPS
from .models import check_N, get_model
from .models.mt import bn_exact
from .status import NoThreshold, Status, Threshold, make_provenance

_MT_BASIS = ("paper eq. (bn) and Proposition prop:spectral (conj13/paper/main.tex; Montaldi–Tokieda Thm 6, "
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
                f"N = {N}: equation (bn) has its only root at lam_r2 = 1 (the equator), where the momentum is fixed "
                "by all of SO(3) and the threshold analysis (|b| < 1) does not apply; Montaldi–Tokieda treat this "
                "case separately (MT p. 15)")
        return NoThreshold(status=Status.PROVED, basis=f"C1(b) > floor(N^2/4) for every b < 1 when N = {N} "
                           "(paper Proposition prop:blocks)", provenance=make_provenance("threshold", inputs, t0),
                           model=mdl.name, N=N, verdict="stable at every radius",
                           geometry={"plane": "plane", "sinh": "H2"}[mdl.name])
    bn = bn_exact(N)
    where = _views_values_relative_to_bn(mdl.name, bn)
    geometry = {"mt": mdl.geometry(bn), "sphere": "sphere", "plane": "plane", "sinh": "H2"}[mdl.name]
    if where != "contains":
        return NoThreshold(status=Status.PROVED,
                           basis=f"b_{N} = {bn} lies outside the {mdl.name} range {mdl.domain_str()}; " + _MT_BASIS,
                           provenance=make_provenance("threshold", inputs, t0), model=mdl.name, N=N,
                           verdict=f"{where} at every radius", geometry=geometry)
    with mp.workdps(WORK_DPS):
        value = mp.mpf(str(sp.N(bn, WORK_DPS + 10)))
        r, xi, cosh_r = mdl.radius(value), mdl.xi(value), mdl.cosh_r(value)
    basis = _SINH_BASIS if mdl.name == "sinh" else _MT_BASIS
    if mdl.name == "plane":
        basis = "the plane has no length scale: every ring with N = 7 sits at the bifurcation point b_7 = 0; " + basis
    return Threshold(status=Status.PROVED, basis=basis, provenance=make_provenance("threshold", inputs, t0),
                     model=mdl.name, N=N, value=value, exact=bn, mode=N // 2, geometry=geometry,
                     stable_side="below", xi=xi, r=r, cosh_r=cosh_r)


_HANDLERS = {"mt": _mt_threshold}


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
