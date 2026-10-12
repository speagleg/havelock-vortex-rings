"""Montaldi–Tokieda family [MT (1.12)] and its geometry-restricted views.

Chart z in C with sigma = 1 + lam|z|^2, symplectic form sum 2/sigma^2 dx^dy, ring of chart radius r0 = 1 and
native parameter b = lam r0^2 = lam_r2 (paper conj13/paper/main.tex, eq. (1)):
    H = -sum_{i<j} log|z_i - z_j|^2 + (n-1) sum_k log sigma_k,   J = sum_k (x_k, y_k, |z_k|^2) / sigma_k.
b > 0 is the sphere (curvature 4 lam), b = 0 the plane, b < 0 the hyperbolic plane with the sinh interaction.
The domain is (-1, 1]: b and 1/b describe the same ring (reflection through the equator), so b > 1 is excluded.
"""
from __future__ import annotations

from fractions import Fraction

import mpmath as mp
import numpy as np
import sympy as sp

from ..errors import InvalidInput

# name -> (lo, hi, lo_closed, hi_closed)
_DOMAINS = {
    "mt": (Fraction(-1), Fraction(1), False, True),
    "sphere": (Fraction(0), Fraction(1), False, True),
    "plane": (Fraction(0), Fraction(0), True, True),
    "sinh": (Fraction(-1), Fraction(0), False, False),
}


def bn_exact(n: int) -> sp.Expr:
    """b_n of paper eq. (2) in closed form (Lean Havelock.Conj13.bn): even n, t = sqrt(n-1):
    (1 + 2t - t^2)/(t^2 + 2t - 1); odd n, s = sqrt(n-3): (2 - s)/(2 + s). b_2 = b_3 = 1 (the equator)."""
    if n % 2 == 0:
        t = sp.sqrt(n - 1)
        return sp.radsimp((1 + 2 * t - t**2) / (t**2 + 2 * t - 1))
    s = sp.sqrt(n - 3)
    return sp.radsimp((2 - s) / (2 + s))


def sinh_q(N: int) -> Fraction:
    """q(N) = cosh^2 r* of the sinh model (sinh/note.tex, Section 3)."""
    if N % 2:
        return Fraction(N - 3, 4)
    return Fraction((N - 2) ** 2, 4 * (N - 1))


def c4_closed(n: int) -> Fraction:
    """Level-set quartic coefficient at b_n (paper Theorems even/odd; Lean Havelock.Conj13.c4closed)."""
    if n % 2 == 0:
        m = n - 2
        return Fraction((n - 2) * (m**4 + 40 * m**2 + 48 * m + 16), 192 * n * (n - 1))
    return Fraction((n - 1) * (n - 3) * (n * n - 9 * n + 38), 128 * n)


class MTView:
    family = "mt"
    param = "lam_r2"
    stable_side = "below"
    family_domain = _DOMAINS["mt"]

    def __init__(self, name: str):
        self.name = name
        self.domain = _DOMAINS[name]

    def domain_str(self) -> str:
        lo, hi, lc, hc = self.domain
        if lo == hi:
            return f"{{{lo}}}"
        return f"{'[' if lc else '('}{lo}, {hi}{']' if hc else ')'}"

    def check_value(self, v: Fraction) -> None:
        lo, hi, lc, hc = self.domain
        if not ((lo < v or (lc and v == lo)) and (v < hi or (hc and v == hi))):
            raise InvalidInput(f"{self.param} = {v} is outside the domain {self.domain_str()} of model {self.name!r}")

    @staticmethod
    def geometry(v) -> str:
        return "sphere" if v > 0 else ("plane" if v == 0 else "H2")

    @staticmethod
    def radius(v):
        """Geodesic ring radius with curvature normalised to +-1 (sphere: colatitude; H2: hyperbolic radius);
        for the plane the chart radius r0 = 1."""
        v = mp.mpf(v)
        if v > 0:
            return 2 * mp.atan(mp.sqrt(v))
        if v == 0:
            return mp.mpf(1)
        return 2 * mp.atanh(mp.sqrt(-v))

    @staticmethod
    def xi(v):
        """tanh^2(r/2) on H2 (equal to -lam_r2); None on the sphere and the plane."""
        v = mp.mpf(v)
        return -v if v < 0 else None

    @staticmethod
    def cosh_r(v):
        v = mp.mpf(v)
        return (1 - v) / (1 + v) if v < 0 else None

    @staticmethod
    def kappa0(v):
        """Symplectic density in log-polar coordinates at the ring: omega = kappa0 ds^dtheta."""
        return 2 / (1 + mp.mpf(v)) ** 2

    @staticmethod
    def chart_radius(v):
        return mp.mpf(1)

    @staticmethod
    def H_mp(z, v):
        n = len(z)
        v = mp.mpf(v)
        h = mp.mpf(0)
        for i in range(n):
            for j in range(i + 1, n):
                h -= mp.log(abs(z[i] - z[j]) ** 2)
        for w in z:
            h += (n - 1) * mp.log(1 + v * abs(w) ** 2)
        return h

    @staticmethod
    def J_mp(z, v):
        v = mp.mpf(v)
        s = [1 + v * abs(w) ** 2 for w in z]
        return (mp.fsum(w.real / si for w, si in zip(z, s)), mp.fsum(w.imag / si for w, si in zip(z, s)),
                mp.fsum(abs(w) ** 2 / si for w, si in zip(z, s)))

    @staticmethod
    def sympy_kernel(x1, y1, x2, y2, p, n):
        """(pair term, one-body term, (J1, J2, J3)) as sympy expressions for the level-set engine."""
        u = x1**2 + y1**2
        s = 1 + p * u
        return (-sp.log((x1 - x2) ** 2 + (y1 - y2) ** 2), (n - 1) * sp.log(s), (x1 / s, y1 / s, u / s))

    @staticmethod
    def numpy_kernel(N: int, v: float):
        """Vectorised dH/d(conj z), inverse symplectic density, H, J and the pair invariant |zi-zj|^2/(si sj)."""
        iu = np.triu_indices(N, 1)

        def dH_dzbar(z):
            dz = z[:, None] - z[None, :]
            np.fill_diagonal(dz, 1.0)
            inv = 1.0 / np.conj(dz)
            np.fill_diagonal(inv, 0.0)
            return -inv.sum(axis=1) + (N - 1) * v * z / (1 + v * np.abs(z) ** 2)

        def winv(z):
            return (1 + v * np.abs(z) ** 2) ** 2 / 2

        def H(z):
            d2 = np.abs(z[:, None] - z[None, :])[iu] ** 2
            return -np.log(d2).sum() + (N - 1) * np.log(1 + v * np.abs(z) ** 2).sum()

        def J(z):
            s = 1 + v * np.abs(z) ** 2
            return np.array([(z.real / s).sum(), (z.imag / s).sum(), (np.abs(z) ** 2 / s).sum()])

        def pairq(z):
            s = 1 + v * np.abs(z) ** 2
            return (np.abs(z[:, None] - z[None, :]) ** 2 / (s[:, None] * s[None, :]))[iu]

        return dH_dzbar, winv, H, J, pairq
