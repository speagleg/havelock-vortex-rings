"""simulate(model, N, value): direct integration of the 2N-dimensional point-vortex system (always NUMERICAL).
Port of conj13/integrate.py to every model. Hamilton's equations x' = H_y / w, y' = -H_x / w with w the
symplectic density, i.e. z' = -2i (dH/d conj z) / w."""
from __future__ import annotations

import math
import time

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

from .errors import InvalidInput
from .limits import DRIFT_FACTOR, SIM_SAMPLES, WORK_DPS
from .models import check_N, get_model, to_fraction
from .status import Simulation, Status, make_provenance


def _positive(name, x, upper=math.inf):
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or not 0 < x <= upper:
        raise InvalidInput(f"{name} must be a finite number in (0, {upper}], got {x!r}")
    return float(x)


def simulate(model: str, N: int, value, *, perturb: float = 1e-4, t_end: float = 100.0, rtol: float = 1e-12,
             seed: int | None = None) -> Simulation:
    """Integrate the ring at `value`, perturbed by relative amplitude `perturb`: along the critical radial mode
    floor(N/2) when seed is None, otherwise by a Gaussian perturbation of every coordinate (numpy seed). DOP853
    with atol = rtol/100. The deviation D(t) is max over pairs of |q_ij(t) - q_ij(ring)| for the model's pair
    invariant q; `growth` = max D / D(0) (NaN when perturb = 0). Energy and momentum drifts are relative to max(1, |initial|) and are
    expected below drift_tolerance = 1e3 rtol."""
    t0 = time.perf_counter()
    mdl = get_model(model)
    N = check_N(N)
    v = to_fraction(value)
    mdl.check_value(v)
    if isinstance(perturb, bool) or not isinstance(perturb, (int, float)) or not math.isfinite(perturb) \
            or not 0 <= perturb < 0.5:
        raise InvalidInput(f"perturb must be in [0, 0.5), got {perturb!r}")
    t_end = _positive("t_end", t_end)
    rtol = _positive("rtol", rtol, 1e-3)
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
        raise InvalidInput(f"seed must be an int or None, got {seed!r}")
    vf = float(v)
    with mp.workdps(WORK_DPS):
        radius = mdl.radius(mp.mpf(v.numerator) / v.denominator)
    dH, winv, H, J, pairq = mdl.numpy_kernel(N, vf)
    rho = float(mdl.chart_radius(vf))
    phi = 2 * np.pi * np.arange(N) / N
    ring = rho * np.exp(1j * phi)
    if seed is None:
        z0 = ring * (1 + perturb * np.cos((N // 2) * phi))
    else:
        rng = np.random.default_rng(seed)
        z0 = ring + perturb * rho * (rng.standard_normal(N) + 1j * rng.standard_normal(N))

    def rhs(t, y):
        z = y[:N] + 1j * y[N:]
        zd = -2j * dH(z) * winv(z)
        return np.concatenate([zd.real, zd.imag])

    ts = np.linspace(0.0, t_end, SIM_SAMPLES)
    sol = solve_ivp(rhs, (0.0, t_end), np.concatenate([z0.real, z0.imag]), method="DOP853", rtol=rtol,
                    atol=rtol * 1e-2, t_eval=ts)
    if not sol.success:
        raise InvalidInput(f"integration failed: {sol.message}")
    Z = sol.y[:N].T + 1j * sol.y[N:].T
    q0 = pairq(ring)
    D = np.array([np.max(np.abs(pairq(z) - q0)) for z in Z])
    D[0] = np.max(np.abs(pairq(z0) - q0))  # from the initial condition itself, not the solver's first sample
    Hs = np.array([H(z) for z in Z])
    Js = np.array([J(z) for z in Z])
    e_drift = float(np.max(np.abs(Hs - Hs[0])) / max(1.0, abs(Hs[0])))
    j_drift = float(np.max(np.linalg.norm(Js - Js[0], axis=1)) / max(1.0, float(np.linalg.norm(Js[0]))))
    growth = float(D.max() / D[0]) if D[0] > 0 else float("nan")
    return Simulation(status=Status.NUMERICAL,
                      basis=("scipy DOP853, rtol = %g, atol = rtol/100; Hamilton's equations of the model in its "
                             "chart (conj13/integrate.py)" % rtol),
                      provenance=make_provenance("simulate", {"model": mdl.name, "N": N, "value": str(v),
                                                              "perturb": perturb, "t_end": t_end, "rtol": rtol,
                                                              "seed": seed}, t0),
                      model=mdl.name, N=N, value=v, times=sol.t, trajectory=Z, initial_deviation=float(D[0]),
                      max_deviation=float(D.max()), growth=growth, energy_drift=e_drift, momentum_drift=j_drift,
                      drift_tolerance=DRIFT_FACTOR * rtol, r=radius)
