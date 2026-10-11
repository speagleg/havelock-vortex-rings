"""stability(model, N, value): the verdict at one parameter value, from exact signs of the Fourier blocks."""
from __future__ import annotations

import time

import mpmath as mp

from .blocks import blocks_exact
from .limits import WORK_DPS
from .models import check_N, get_model, to_fraction
from .status import ModeBlock, Stability, Status, make_provenance

_BASIS = {
    "mt": ("exact rational evaluation of a_m = C1 - m(n-m), b_m = m(n-m), g_m = 0 (paper Proposition prop:blocks, "
           "MT (2.6)); EM-stable = a_m > 0 and b_m > 0 for m = 1..floor(N/2), which makes the reduced Hessian "
           "positive definite (paper Proposition prop:spectral, MT Thm 2)"),
    "tanh": ("exact rational evaluation of the blocks of tanh/note_tanh.tex Thm 2 (Fourier blocks, PROVED); "
             "criteria of Prop. lin: spectrally stable iff every a_m b_m >= 0, EM-stable iff a_m > 0 and "
             "a_m b_m - g_m^2 > 0 for every m (Lyapunov conclusions only modulo SO(2,1), note Remark em)"),
}


def _annotate(basis: str, blocks, verdict: str) -> str:
    """Append a phrase for every reason a "spectrally-stable-only" verdict is not EM-stable (unstable and
    EM-stable verdicts need none). The reasons are exhaustive: a_m = 0; b_m = 0 with a_m != 0; a_m b_m > 0 with
    a_m b_m - g_m^2 <= 0 (Krein window); a_m < 0, b_m < 0 with a_m b_m - g_m^2 > 0 (negative definite)."""
    if verdict != "spectrally-stable-only":
        return basis
    zero = [k.m for k in blocks if k.a == 0]
    if zero:
        basis += f"; a_m = 0 exactly for m = {zero}: a degenerate (bifurcation) point, see havelock.nonlinear"
    bzero = [k.m for k in blocks if k.b == 0 and k.a != 0]
    if bzero:
        basis += (f"; b_m = 0 exactly for m = {bzero} (a_m != 0): the Hessian is singular on that mode, "
                  "so it is not positive definite and nonlinear stability is undecided")
    window = [k.m for k in blocks if k.a * k.b > 0 and k.a * k.b - k.g * k.g <= 0]
    if window:
        basis += f"; a_m b_m > 0 but a_m b_m - g_m^2 <= 0 for m = {window}: Krein window, nonlinear stability undecided"
    neg = [k.m for k in blocks if k.a < 0 and k.b < 0 and k.a * k.b - k.g * k.g > 0]
    if neg:
        basis += (f"; the block is negative definite for m = {neg} (a_m < 0, b_m < 0, a_m b_m - g_m^2 > 0): "
                  "spectrally stable but not an energy minimum, so no energy-momentum conclusion")
    return basis


def _eigs(a, b, g, kappa0):
    with mp.workdps(WORK_DPS):
        root = mp.sqrt(-mp.mpf(a.numerator) / a.denominator * mp.mpf(b.numerator) / b.denominator + 0j)
        ig = 1j * mp.mpf(g.numerator) / g.denominator
        return (complex((-ig + root) / kappa0), complex((-ig - root) / kappa0))


def stability(model: str, N: int, value) -> Stability:
    """Verdict at lam_r2 (MT family) or xi (tanh) = value. "EM-stable": a_m > 0 and a_m b_m - g_m^2 > 0 for
    every mode m >= 1 (positive-definite energy–momentum Hessian); "spectrally-stable-only": every a_m b_m >= 0
    but the Hessian is not positive definite (a_m = 0 at a bifurcation point; b_m = 0 with a_m != 0; a Krein
    window a_m b_m > 0 >= a_m b_m - g_m^2; or a negative-definite block a_m, b_m < 0, a_m b_m - g_m^2 > 0);
    `basis` names each case that occurs;
    "unstable": some a_m b_m < 0 (a real pair of eigenvalues). The value is used exactly (floats by their
    binary value), so every sign is decided exactly."""
    t0 = time.perf_counter()
    mdl = get_model(model)
    N = check_N(N)
    v = to_fraction(value)
    mdl.check_value(v)
    with mp.workdps(WORK_DPS):
        vm = mp.mpf(v.numerator) / v.denominator
        kappa0, radius = mdl.kappa0(vm), mdl.radius(vm)
    blocks = []
    for m in range(1, N // 2 + 1):
        a, b, g = blocks_exact(mdl.family, N, m, v)
        blocks.append(ModeBlock(m=m, a=a, b=b, g=g, eigenvalues=_eigs(a, b, g, kappa0)))
    em = all(k.a > 0 and k.a * k.b - k.g * k.g > 0 for k in blocks)
    spectral = all(k.a * k.b >= 0 for k in blocks)
    verdict = "EM-stable" if em else ("spectrally-stable-only" if spectral else "unstable")
    basis = _BASIS[mdl.family]
    basis = _annotate(basis, blocks, verdict)
    status = Status.PROVED if mdl.family == "mt" else Status.CERTIFIED
    return Stability(status=status, basis=basis,
                     provenance=make_provenance("stability", {"model": mdl.name, "N": N, "value": str(v)}, t0),
                     model=mdl.name, N=N, value=v, verdict=verdict, geometry=mdl.geometry(v), blocks=tuple(blocks),
                     r=radius)
