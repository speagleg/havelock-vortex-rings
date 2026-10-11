"""Model registry, input validation and exact parameter parsing.

A model object exposes: name, family ("mt" or "tanh"), param (native parameter name), stable_side
("below": stable for values below the threshold; "above": stable above it), family_domain (lo, hi, lo_closed,
hi_closed) used by certification, check_value(v), geometry(v), radius(v), xi(v), cosh_r(v), kappa0(v),
chart_radius(v), H_mp(z, v), J_mp(z, v), numpy_kernel(N, v) and sympy_kernel(x1, y1, x2, y2, p, n).
"""
from __future__ import annotations

import math
import numbers
from fractions import Fraction

import mpmath as mp
import sympy as sp

from ..errors import InvalidInput
from .mt import MTView
from .tanh import TanhModel

MODELS = {name: MTView(name) for name in ("mt", "sphere", "plane", "sinh")}
MODELS["tanh"] = TanhModel()


def get_model(name):
    """The model object registered under `name` (case-insensitive)."""
    if not isinstance(name, str) or name.lower() not in MODELS:
        raise InvalidInput(f"unknown model {name!r}; known models: {', '.join(MODELS)}")
    return MODELS[name.lower()]


def check_N(N) -> int:
    """N as an int >= 2 (bools and non-integral values are rejected)."""
    if isinstance(N, bool) or not isinstance(N, numbers.Integral):
        raise InvalidInput(f"N must be an integer, got {N!r}")
    if N < 2:
        raise InvalidInput(f"N must be >= 2, got {N}")
    return int(N)


def to_fraction(value) -> Fraction:
    """Exact rational value of a parameter. Floats are converted exactly (their binary value); strings may be
    decimals or p/q; sympy Rationals and mpmath mpf are accepted. Irrational or non-finite values are rejected."""
    if isinstance(value, bool):
        raise InvalidInput("a parameter value cannot be a bool")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, numbers.Integral):
        return Fraction(int(value))
    if isinstance(value, float):
        if not math.isfinite(value):
            raise InvalidInput(f"parameter value must be finite, got {value!r}")
        return Fraction(value)
    if isinstance(value, str):
        try:
            return Fraction(value.strip())
        except (ValueError, ZeroDivisionError):
            raise InvalidInput(f"cannot parse {value!r} as a rational number") from None
    if isinstance(value, sp.Rational):
        return Fraction(int(value.p), int(value.q))
    if isinstance(value, mp.mpf):
        if not mp.isfinite(value):
            raise InvalidInput(f"parameter value must be finite, got {value!r}")
        man, exp = value.man_exp
        return Fraction(int(man)) * Fraction(2) ** int(exp)
    raise InvalidInput(f"unsupported parameter type {type(value).__name__}; use int, float, str, Fraction, "
                       "sympy.Rational or mpmath.mpf")
