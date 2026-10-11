"""nonlinear(model, N): the quartic coefficient c4 of the reduced Hamiltonian along the kernel at the bifurcation
point, and the verdict it gives. MT family: closed forms (PROVED). The level-set engine below (a port of
conj13/review/c4_levelset.py) evaluates the same quantity numerically for any model; it is used for the tanh model
and to cross-check the closed forms."""
from __future__ import annotations

import functools
import time

import mpmath as mp
import sympy as sp

from .errors import OutOfValidityRange
from .limits import LEVELSET
from .models import check_N, get_model
from .models.mt import c4_closed
from .status import Nonlinear, NoThreshold, Status, make_provenance
from .threshold import threshold

_S = sp.symbols("x1 y1 x2 y2 p n")


@functools.lru_cache(maxsize=None)
def _compiled(name: str):
    x1, y1, x2, y2, p, n = _S
    pair, one, Js = get_model(name).sympy_kernel(x1, y1, x2, y2, p, n)
    v4, v2 = [x1, y1, x2, y2], [x1, y1]

    def lam(e, vs, extra):
        return (sp.lambdify(vs + extra, e, "mpmath"),
                sp.lambdify(vs + extra, [sp.diff(e, a) for a in vs], "mpmath"),
                sp.lambdify(vs + extra, [[sp.diff(e, a, b) for b in vs] for a in vs], "mpmath"))

    return lam(pair, v4, [p]), lam(one, v2, [p, n]), [lam(J, v2, [p, n]) for J in Js]


def _H(name, q, N, p):
    pair, one, _ = _compiled(name)
    dim = 2 * N
    val, g, Hs = mp.mpf(0), [mp.mpf(0)] * dim, mp.zeros(dim, dim)
    for i in range(N):
        for j in range(i + 1, N):
            a, ids = (q[2 * i], q[2 * i + 1], q[2 * j], q[2 * j + 1]), (2 * i, 2 * i + 1, 2 * j, 2 * j + 1)
            val += pair[0](*a, p)
            d, h = pair[1](*a, p), pair[2](*a, p)
            for r in range(4):
                g[ids[r]] += d[r]
                for c in range(4):
                    Hs[ids[r], ids[c]] += h[r][c]
    for i in range(N):
        a = (q[2 * i], q[2 * i + 1])
        val += one[0](*a, p, N)
        d, h = one[1](*a, p, N), one[2](*a, p, N)
        for r in range(2):
            g[2 * i + r] += d[r]
            for c in range(2):
                Hs[2 * i + r, 2 * i + c] += h[r][c]
    return val, g, Hs


def _J(name, q, N, p):
    _, _, Js = _compiled(name)
    dim = 2 * N
    vals, grads, hess = [], [], []
    for J in Js:
        v, g, Hs = mp.mpf(0), [mp.mpf(0)] * dim, mp.zeros(dim, dim)
        for i in range(N):
            a = (q[2 * i], q[2 * i + 1])
            v += J[0](*a, p, N)
            d, h = J[1](*a, p, N), J[2](*a, p, N)
            for r in range(2):
                g[2 * i + r] += d[r]
                for c in range(2):
                    Hs[2 * i + r, 2 * i + c] += h[r][c]
        vals.append(v), grads.append(g), hess.append(Hs)
    return vals, grads, hess


def _vectors(N):
    l, dim = N // 2, 2 * N
    e, e2, t = [mp.mpf(0)] * dim, [mp.mpf(0)] * dim, [mp.mpf(0)] * dim
    for k in range(N):
        ph = 2 * mp.pi * k / N
        e[2 * k], e[2 * k + 1] = mp.cos(l * ph) * mp.cos(ph), mp.cos(l * ph) * mp.sin(ph)
        e2[2 * k], e2[2 * k + 1] = mp.sin(l * ph) * mp.cos(ph), mp.sin(l * ph) * mp.sin(ph)
        t[2 * k], t[2 * k + 1] = -mp.sin(ph), mp.cos(ph)
    unit = lambda v: [a / mp.sqrt(mp.fsum(b * b for b in v)) for a in v]  # noqa: E731
    return unit(e), unit(e2), unit(t)


def _level_set_min(name, N, p, eps, max_iter=60):
    """min { H(q) - H(q0) : J(q) = J(q0), <q - q0, e> = eps, <q - q0, t> = 0, <q - q0, e'> = 0 (N odd) } by
    Newton's method on the KKT system (Cartesian chart coordinates, Euclidean-unit kernel vector e)."""
    mdl = get_model(name)
    dim = 2 * N
    rho = mdl.chart_radius(p)
    q0 = []
    for k in range(N):
        ph = 2 * mp.pi * k / N
        q0 += [rho * mp.cos(ph), rho * mp.sin(ph)]
    e, e2, t = _vectors(N)
    J0, G0, _ = _J(name, q0, N, p)
    _, gH0, _ = _H(name, q0, N, p)
    omega = mp.fsum(gH0[i] * G0[2][i] for i in range(dim)) / mp.fsum(G0[2][i] ** 2 for i in range(dim))
    lin = [e, t] + ([e2] if N % 2 else [])
    rhs = [eps, 0] + ([0] if N % 2 else [])
    nc = len(lin) + 3
    mult = [mp.mpf(0)] * nc
    mult[len(lin) + 2] = omega
    q = [q0[i] + eps * e[i] for i in range(dim)]
    tol = mp.mpf(10) ** (-mp.mp.dps + 8)
    for _ in range(max_iter):
        _, g, Hs = _H(name, q, N, p)
        Jv, Jg, Jh = _J(name, q, N, p)
        cons = [mp.fsum(lin[a][i] * (q[i] - q0[i]) for i in range(dim)) - rhs[a] for a in range(len(lin))]
        cons += [Jv[k] - J0[k] for k in range(3)]
        C = lin + Jg
        L = Hs.copy()
        for k in range(3):
            L -= mult[len(lin) + k] * Jh[k]
        K, r = mp.zeros(dim + nc, dim + nc), mp.zeros(dim + nc, 1)
        for i in range(dim):
            for j in range(dim):
                K[i, j] = L[i, j]
            r[i] = g[i] - mp.fsum(mult[a] * C[a][i] for a in range(nc))
        for a in range(nc):
            for i in range(dim):
                K[i, dim + a], K[dim + a, i] = -C[a][i], C[a][i]
            r[dim + a] = cons[a]
        d = mp.lu_solve(K, -r)
        q = [q[i] + d[i] for i in range(dim)]
        mult = [mult[a] + d[dim + a] for a in range(nc)]
        if mp.norm(d) < tol:
            break
    else:
        raise OutOfValidityRange(f"level-set Newton iteration did not converge (N={N}, eps={eps})")
    return _H(name, q, N, p)[0] - _H(name, q0, N, p)[0]


def c4_level_set(name: str, N: int, p, *, dps: int = LEVELSET["dps"], steps: int = LEVELSET["steps"],
                 eps_denominator: int = LEVELSET["eps_denominator"]):
    """Level-set effective quartic (paper Definition def:c4) by Newton on the KKT system at eps = k/eps_denominator,
    k = 1..steps, symmetrised in +-eps, then Neville extrapolation of f(eps)/eps^4 to eps = 0. NUMERICAL."""
    with mp.workdps(dps):
        p = mp.mpf(p)
        epss = [mp.mpf(k) / eps_denominator for k in range(1, steps + 1)]
        fs = [(_level_set_min(name, N, p, s) + _level_set_min(name, N, p, -s)) / 2 / s**4 for s in epss]
        xs = [s * s for s in epss]
        P = list(fs)
        for k in range(1, len(xs)):
            for i in range(len(xs) - k):
                P[i] = (xs[i + k] * P[i] - xs[i] * P[i + 1]) / (xs[i + k] - xs[i])
        return +P[0]


def _mt_nonlinear(mdl, N, t, t0, inputs):
    c4 = c4_closed(N)
    basis = ("paper Theorem thm:even (N even) / thm:odd (N odd): c4 = "
             + ("(n-2)(m^4+40m^2+48m+16)/(192n(n-1)), m = n-2" if N % 2 == 0
                else "(n-1)(n-3)(n^2-9n+38)/(128n)")
             + " > 0, so the ring is Lyapunov stable modulo rotations at b_n (paper Theorem thm:main); "
             "Lean Havelock.Conj13.c4_at_bn (the algebra downstream of the series coefficients)")
    return Nonlinear(status=Status.PROVED, basis=basis, provenance=make_provenance("nonlinear", inputs, t0),
                     model=mdl.name, N=N, value=t.value, mode=t.mode, c4=c4, verdict="stable",
                     method="closed form", r=t.r)


def _tanh_nonlinear(mdl, N, t, t0, inputs):
    if N % 2:
        raise OutOfValidityRange(
            f"N = {N} is odd: at the spectral threshold xi* the critical block has a_m = 0, so its determinant is "
            "a b - g^2 = -g^2 <= 0 and the block is not positive definite (tanh/note_tanh.tex; for N = 9, 11, 13, 15 "
            "the note proves a Krein window above xi*, where the ring is linearly stable but the Hessian is "
            "indefinite). The energy-momentum argument behind the quartic (c4) test needs a definite block, so c4 "
            "is not defined here; nonlinear stability at xi* is undecided")
    c4 = c4_level_set("tanh", N, t.value)
    verdict = "stable" if c4 > mp.mpf("1e-8") else ("unstable" if c4 < -mp.mpf("1e-8") else "undecided")
    basis = ("NUMERICAL: level-set engine (Newton on the KKT system on J^{-1}(mu), eps = k/200, k = 1..6, Neville "
             "extrapolation; port of conj13/review/c4_levelset.py), Euclidean-unit kernel vector in the Poincare-disk "
             "chart, H as in tanh/note_tanh.tex eq. (H). c4 > 0 indicates a strict (quartic) minimum of the reduced "
             "Hamiltonian; the energy–momentum argument of conj13/paper Section nonlinear has not been checked for "
             "this model, and any stability conclusion is modulo SO(2,1). No closed form for tanh c4 is known")
    return Nonlinear(status=Status.NUMERICAL, basis=basis, provenance=make_provenance("nonlinear", inputs, t0),
                     model=mdl.name, N=N, value=t.value, mode=t.mode, c4=c4, verdict=verdict,
                     method="level-set engine", r=t.r)


_HANDLERS = {"mt": _mt_nonlinear, "tanh": _tanh_nonlinear}


def nonlinear(model: str, N: int) -> Nonlinear:
    """c4 at the bifurcation point of the N-ring and the nonlinear verdict ("stable", "unstable", "undecided")."""
    t0 = time.perf_counter()
    mdl = get_model(model)
    N = check_N(N)
    inputs = {"model": mdl.name, "N": N}
    t = threshold(mdl.name, N)
    if isinstance(t, NoThreshold):
        raise OutOfValidityRange(f"model {mdl.name!r}, N = {N} has no bifurcation point ({t.verdict})")
    return _HANDLERS[mdl.family](mdl, N, t, t0, inputs)
