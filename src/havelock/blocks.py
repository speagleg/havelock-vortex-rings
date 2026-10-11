"""Exact Fourier blocks of the energy–momentum Hessian at the regular N-ring, and an mpmath fast path.

In log-polar coordinates z_k = rho e^{2 pi i k/N} e^{s_k + i theta_k}, with the Euclidean structure
sum (ds_k^2 + dtheta_k^2), the Hessian of F = H - Omega J at the ring is block-diagonal over the Fourier modes
m = 0..N-1. On mode m it is the Hermitian matrix [[a_m, i g_m], [-i g_m, b_m]]. Mode 0 carries the rotation
(b_0 = g_0 = 0). Each family supplies integer numerators:
    a_m = A / den_a,   b_m = B / den_b,   g_m = G / den_b,   a_m b_m - g_m^2 = D / (den_a den_b),
with den_a, den_b > 0 on the family's domain and D = A B - (den_a / den_b) G^2.

MT family (paper Proposition prop:blocks, MT (2.6)): a_m = C1 - m(n-m), b_m = m(n-m), g_m = 0,
C1 = 2(n-1)(1+b^2)/(1+b)^2.
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from fractions import Fraction

import flint
import mpmath as mp
import sympy as sp

from .errors import InvalidInput

X = sp.Symbol("x")
_x = flint.fmpz_poly([0, 1])
_one = flint.fmpz_poly([1])


@dataclass(frozen=True)
class Numerators:
    A: flint.fmpz_poly
    B: flint.fmpz_poly
    G: flint.fmpz_poly
    D: flint.fmpz_poly
    den_a: flint.fmpz_poly
    den_b: flint.fmpz_poly


def _mt_numerators(N: int, m: int):
    T = m * (N - m)
    A = 2 * (N - 1) * (_one + _x**2) - T * (_one + _x) ** 2
    return A, flint.fmpz_poly([T]), flint.fmpz_poly([]), (_one + _x) ** 2, _one


def _mt_exprs(N: int, m: int):
    C1 = 2 * (N - 1) * (1 + X**2) / (1 + X) ** 2
    return C1 - m * (N - m), sp.Integer(m * (N - m)), sp.Integer(0)


def _mono(k: int) -> flint.fmpz_poly:
    return flint.fmpz_poly([0] * k + [1])


def _tanh_numerators(N: int, m: int):
    """note Thm 2 cleared of denominators: A = 2(1-x)^2(1-X)^2 a_m, B = 2(1-X)^2 b_m, G = 2(1-X)^2 g_m, X = x^N."""
    X = _mono(N)
    omx, omX, opx = _one - _x, _one - X, _one + _x
    A = (2 * (N - 1) * opx * omx * omX**2 - 4 * _x**2 * omX**2 - m * (N - m) * omx**2 * omX**2
         + 4 * N * opx * omx * X * omX - 2 * N * N * X * omx**2
         - N * omx**2 * (m * _mono(m) + (N - m) * _mono(N - m) + (N - m) * _mono(N + m) + m * _mono(2 * N - m)))
    Wm = m * _mono(m) * omX + N * _mono(m + N)            # (1-X)^2 W_m
    WNm = (N - m) * _mono(N - m) * omX + N * _mono(2 * N - m)  # (1-X)^2 W_{N-m}
    B = m * (N - m) * omX**2 + 2 * N * N * X - N * (Wm + WNm)
    G = N * (Wm - WNm)
    return A, B, G, 2 * omx**2 * omX**2, 2 * omX**2


def _tanh_exprs(N: int, m: int):
    """note Thm 2, transcribed literally."""
    c = (1 + X) / (1 - X)
    XN = X**N

    def W(k):
        return k * X**k / (1 - XN) + N * X ** (k + N) / (1 - XN) ** 2

    a = ((N - 1) * c - (c - 1) ** 2 / 2 - sp.Rational(m * (N - m), 2) + 2 * N * c * XN / (1 - XN)
         - N**2 * XN / (1 - XN) ** 2
         - sp.Rational(N, 2) * (m * X**m + (N - m) * X ** (N - m) + (N - m) * X ** (N + m) + m * X ** (2 * N - m))
         / (1 - XN) ** 2)
    b = sp.Rational(m * (N - m), 2) + N * W(N) - sp.Rational(N, 2) * (W(m) + W(N - m))
    g = sp.Rational(N, 2) * (W(m) - W(N - m))
    return a, b, g


_NUMERATORS = {"mt": _mt_numerators, "tanh": _tanh_numerators}
_EXPRS = {"mt": _mt_exprs, "tanh": _tanh_exprs}


def _check(family: str, N: int, m: int) -> None:
    if family not in _NUMERATORS:
        raise InvalidInput(f"no Fourier blocks for family {family!r}")
    if not 0 <= m <= N - 1:
        raise InvalidInput(f"mode m must satisfy 0 <= m <= N-1, got m={m}, N={N}")


@functools.lru_cache(maxsize=4096)
def numerators(family: str, N: int, m: int) -> Numerators:
    """Integer numerators and positive denominators of a_m, b_m, g_m and a_m b_m - g_m^2."""
    _check(family, N, m)
    A, B, G, den_a, den_b = _NUMERATORS[family](N, m)
    q, r = divmod(den_a, den_b)
    assert r == 0
    return Numerators(A=A, B=B, G=G, D=A * B - q * G * G, den_a=den_a, den_b=den_b)


def block_exprs(family: str, N: int, m: int):
    """(a_m, b_m, g_m) as sympy rational functions of X (the native parameter), transcribed from the papers."""
    _check(family, N, m)
    return _EXPRS[family](N, m)


def _fq(v: Fraction) -> flint.fmpq:
    return flint.fmpq(v.numerator, v.denominator)


def _frac(q) -> Fraction:
    return Fraction(int(q.p), int(q.q))


def blocks_exact(family: str, N: int, m: int, v: Fraction):
    """(a_m, b_m, g_m) at the rational point v, as exact Fractions."""
    nm = numerators(family, N, m)
    q = _fq(Fraction(v))
    da, db = nm.den_a(q), nm.den_b(q)
    if da == 0 or db == 0:
        raise InvalidInput(f"the blocks are singular at {v}")
    return _frac(nm.A(q) / da), _frac(nm.B(q) / db), _frac(nm.G(q) / db)


def _horner(p: flint.fmpz_poly, v):
    acc = mp.mpf(0)
    for c in reversed(p.coeffs()):
        acc = acc * v + int(c)
    return acc


def blocks_mp(family: str, N: int, m: int, v):
    """(a_m, b_m, g_m) at v in the current mpmath precision (fast path; m = 0 gives the rotation block)."""
    nm = numerators(family, N, m)
    v = mp.mpf(v)
    da, db = _horner(nm.den_a, v), _horner(nm.den_b, v)
    return _horner(nm.A, v) / da, _horner(nm.B, v) / db, _horner(nm.G, v) / db
