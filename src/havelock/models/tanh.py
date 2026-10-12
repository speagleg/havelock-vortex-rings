"""The tanh (Kimura) model on H^2 (tanh/note_tanh.tex, eq. (1)).

Poincare disk, metric 4|dz|^2/(1-|z|^2)^2 (curvature -1), symplectic form sum 4/(1-|z_k|^2)^2 dx^dy,
    H = -sum_{j<k} ln tanh(d_jk/2) = sum_{j<k} ( -ln|z_j - z_k| + ln|1 - z_j conj z_k| ),
    J = sum_k ( 2x_k, 2y_k, 1 + |z_k|^2 ) / (1 - |z_k|^2)      (J3 = sum cosh r_k generates rotations).
Native parameter xi = tanh^2(r/2) in (0, 1): the ring has chart radius sqrt(xi) and geodesic radius r.
"""
from __future__ import annotations

from fractions import Fraction

import mpmath as mp
import numpy as np
import sympy as sp

from ..errors import InvalidInput


class TanhModel:
    name = "tanh"
    family = "tanh"
    param = "xi"
    stable_side = "above"
    family_domain = (Fraction(0), Fraction(1), False, False)
    domain = family_domain

    @staticmethod
    def domain_str() -> str:
        return "(0, 1)"

    def check_value(self, v: Fraction) -> None:
        if not 0 < v < 1:
            raise InvalidInput(f"xi = {v} is outside the domain (0, 1) of model 'tanh'")

    @staticmethod
    def geometry(v) -> str:
        return "H2"

    @staticmethod
    def radius(v):
        return 2 * mp.atanh(mp.sqrt(mp.mpf(v)))

    @staticmethod
    def xi(v):
        return mp.mpf(v)

    @staticmethod
    def cosh_r(v):
        v = mp.mpf(v)
        return (1 + v) / (1 - v)

    @staticmethod
    def kappa0(v):
        v = mp.mpf(v)
        return 4 * v / (1 - v) ** 2

    @staticmethod
    def chart_radius(v):
        return mp.sqrt(mp.mpf(v))

    @staticmethod
    def H_mp(z, v):
        n = len(z)
        h = mp.mpf(0)
        for i in range(n):
            for j in range(i + 1, n):
                h += -mp.log(abs(z[i] - z[j])) + mp.log(abs(1 - z[i] * mp.conj(z[j])))
        return h

    @staticmethod
    def J_mp(z, v):
        u = [abs(w) ** 2 for w in z]
        return (mp.fsum(2 * w.real / (1 - ui) for w, ui in zip(z, u)),
                mp.fsum(2 * w.imag / (1 - ui) for w, ui in zip(z, u)),
                mp.fsum((1 + ui) / (1 - ui) for ui in u))

    @staticmethod
    def sympy_kernel(x1, y1, x2, y2, p, n):
        u = x1**2 + y1**2
        re, im = 1 - (x1 * x2 + y1 * y2), y1 * x2 - x1 * y2
        pair = -sp.log((x1 - x2) ** 2 + (y1 - y2) ** 2) / 2 + sp.log(re**2 + im**2) / 2
        return pair, sp.Integer(0), (2 * x1 / (1 - u), 2 * y1 / (1 - u), (1 + u) / (1 - u))

    @staticmethod
    def numpy_kernel(N: int, v: float):
        iu = np.triu_indices(N, 1)

        def dH_dzbar(z):
            dz = z[:, None] - z[None, :]
            np.fill_diagonal(dz, 1.0)
            inv = 1.0 / np.conj(dz)
            np.fill_diagonal(inv, 0.0)
            T = -0.5 * z[None, :] / (1 - np.conj(z)[:, None] * z[None, :])
            np.fill_diagonal(T, 0.0)
            return -0.5 * inv.sum(axis=1) + T.sum(axis=1)

        def winv(z):
            return (1 - np.abs(z) ** 2) ** 2 / 4

        def H(z):
            d = np.abs(z[:, None] - z[None, :])[iu]
            return -np.log(d).sum() + np.log(np.abs(1 - z[:, None] * np.conj(z)[None, :])[iu]).sum()

        def J(z):
            u = np.abs(z) ** 2
            return np.array([(2 * z.real / (1 - u)).sum(), (2 * z.imag / (1 - u)).sum(), ((1 + u) / (1 - u)).sum()])

        def pairq(z):
            return (np.abs(z[:, None] - z[None, :]) ** 2 / np.abs(1 - z[:, None] * np.conj(z)[None, :]) ** 2)[iu]

        return dH_dzbar, winv, H, J, pairq
