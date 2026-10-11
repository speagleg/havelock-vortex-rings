"""Status labels, provenance, and the result classes returned by every public function."""
from __future__ import annotations

import enum
import json
import math
import platform
import time
from dataclasses import dataclass, field, fields, is_dataclass
from fractions import Fraction
from typing import Any

import mpmath as mp
import numpy as np
import sympy as sp


class Status(str, enum.Enum):
    """How an answer is known. PROVED: a paper theorem covers this exact case. CERTIFIED: exact or ball
    arithmetic computed at run time. NUMERICAL: floating point with an error estimate."""

    PROVED = "PROVED"
    CERTIFIED = "CERTIFIED"
    NUMERICAL = "NUMERICAL"

    def __str__(self) -> str:
        return self.value


VERDICTS = ("EM-stable", "spectrally-stable-only", "unstable")


@dataclass(frozen=True)
class Provenance:
    function: str
    inputs: dict
    version: str
    elapsed_s: float
    python: str
    libraries: dict


def make_provenance(function: str, inputs: dict, t0: float) -> Provenance:
    """Provenance record for a call that started at time.perf_counter() == t0."""
    import flint
    import scipy

    from . import __version__

    return Provenance(function=function, inputs=dict(inputs), version=__version__,
                      elapsed_s=time.perf_counter() - t0, python=platform.python_version(),
                      libraries={"mpmath": mp.__version__, "sympy": sp.__version__, "numpy": np.__version__,
                                 "scipy": scipy.__version__, "python-flint": flint.__version__})


def jsonable(v: Any) -> Any:
    """Convert a result (or any value inside one) to plain JSON types. Fractions become "p/q" strings and
    mpmath numbers become decimal strings with all their digits, so no precision is lost."""
    if isinstance(v, enum.Enum):
        return v.value
    if v is None or isinstance(v, (bool, int, str)):
        return v
    if isinstance(v, float):
        return v if math.isfinite(v) else str(v)
    if isinstance(v, Fraction):
        return str(v)
    if isinstance(v, mp.mpf):
        if not mp.isfinite(v):
            return str(v)
        man, _ = v.man_exp
        return mp.nstr(v, max(15, int(abs(int(man)).bit_length() * 0.30103) + 1))
    if isinstance(v, (complex, mp.mpc)):
        return [jsonable(v.real), jsonable(v.imag)]
    if isinstance(v, np.ndarray):
        return [jsonable(x) for x in v.tolist()]
    if isinstance(v, np.generic):
        return jsonable(v.item())
    if is_dataclass(v) and not isinstance(v, type):
        return {f.name: jsonable(getattr(v, f.name)) for f in fields(v) if not f.name.startswith("_")}
    if isinstance(v, dict):
        return {str(k): jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [jsonable(x) for x in v]
    if isinstance(v, sp.Basic):
        return str(v)
    return str(v)


@dataclass(kw_only=True)
class Result:
    """Base class of every public result. `status` is the label, `basis` says why, `fallback_from` records a
    downgrade (for example CERTIFIED -> NUMERICAL when strict=False and a size cap was hit)."""

    status: Status
    basis: str
    provenance: Provenance
    fallback_from: Status | None = None

    def to_dict(self) -> dict:
        d = jsonable(self)
        d["type"] = type(self).__name__
        return d

    def to_json(self, **kwargs) -> str:
        return json.dumps(self.to_dict(), **kwargs)


@dataclass(kw_only=True)
class Threshold(Result):
    """A stability threshold in the model's native parameter. `stable_side` is "below" (MT family: stable for
    lam_r2 < value) or "above" (tanh: stable for xi > value). `exact` is a sympy expression or an
    IsolatedRoot (minimal polynomial with an isolating interval)."""

    model: str
    N: int
    value: Any
    exact: Any
    mode: int
    geometry: str
    stable_side: str
    xi: Any = None
    r: Any = None
    cosh_r: Any = None
    _certification: Any = field(default=None, repr=False, compare=False)


@dataclass(kw_only=True)
class NoThreshold(Result):
    """The ring has no threshold in this model/geometry: `verdict` is "stable at every radius" or
    "unstable at every radius"."""

    model: str
    N: int
    verdict: str
    geometry: str


@dataclass(frozen=True)
class ModeBlock:
    """Fourier block of mode m at a point: the Hermitian matrix [[a, i g], [-i g, b]] and the two eigenvalues
    sigma = (-i g +- sqrt(-a b)) / kappa0 of the linearised flow on that mode."""

    m: int
    a: Fraction
    b: Fraction
    g: Fraction
    eigenvalues: tuple


@dataclass(kw_only=True)
class Stability(Result):
    """Verdict at one parameter value. `r` is the geodesic ring radius there (curvature normalised to +-1; the
    chart radius 1 for the plane)."""

    model: str
    N: int
    value: Fraction
    verdict: str
    geometry: str
    blocks: tuple
    r: Any


@dataclass(kw_only=True)
class Nonlinear(Result):
    """Quartic coefficient c4 of the reduced Hamiltonian along the kernel at the bifurcation point. `r` is the
    geodesic ring radius at the bifurcation point (the threshold's r)."""

    model: str
    N: int
    value: Any
    mode: int
    c4: Any
    verdict: str
    method: str
    r: Any


@dataclass(kw_only=True)
class Simulation(Result):
    """Direct integration of the 2N-dimensional system (always NUMERICAL). `trajectory[i, k]` is vortex k at
    `times[i]` in the model's chart. Deviation is measured by the model's pair invariant (a function of the
    geodesic distance), so it is blind to rigid motions. `r` is the geodesic ring radius at `value` (curvature
    normalised to +-1; the chart radius 1 for the plane)."""

    model: str
    N: int
    value: Fraction
    times: Any
    trajectory: Any
    initial_deviation: float
    max_deviation: float
    growth: float
    energy_drift: float
    momentum_drift: float
    drift_tolerance: float
    r: Any
