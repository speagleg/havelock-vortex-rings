import math

import numpy as np
import pytest

import havelock as hv
from havelock.simulate import simulate

# Recorded while planning (t_end = 100, rtol = 1e-12, perturb = 1e-4 along mode floor(N/2)):
#   sinh N=9 at -0.3:  growth 1.00, energy drift 6.8e-12, momentum drift 6.6e-12
#   sinh N=9 at -0.05: growth 1.0e4
#   sphere N=4 at 0.5: growth 1.0e4;  sphere N=4 at 0.1: growth 1.00


def test_em_stable_ring_stays_close():
    s = simulate("sinh", 9, "-0.3")
    assert s.status == hv.Status.NUMERICAL
    assert s.growth < 1.5 and s.max_deviation < 1.5 * s.initial_deviation
    assert s.energy_drift < s.drift_tolerance and s.momentum_drift < s.drift_tolerance
    assert s.drift_tolerance == pytest.approx(1e-9)
    assert s.trajectory.shape == (2001, 9) and s.times[-1] == pytest.approx(100.0)


def test_ring_below_threshold_grows():
    s = simulate("sinh", 9, "-0.05")
    assert s.growth > 100
    assert s.energy_drift < s.drift_tolerance and s.momentum_drift < s.drift_tolerance


def test_sphere_both_sides():
    assert simulate("sphere", 4, "0.1").growth < 1.5
    assert simulate("sphere", 4, "0.5").growth > 100


def test_unperturbed_ring_is_a_relative_equilibrium():
    s = simulate("sinh", 9, "-0.3", perturb=0.0)
    assert s.initial_deviation == 0.0 and s.max_deviation < 1e-8 and math.isnan(s.growth)  # recorded: 6.8e-11


def test_seeded_perturbation_is_reproducible():
    a = simulate("plane", 6, 0, seed=7, t_end=10.0)
    b = simulate("plane", 6, 0, seed=7, t_end=10.0)
    assert np.array_equal(a.trajectory, b.trajectory)


@pytest.mark.parametrize("kw", [{"perturb": -1e-3}, {"perturb": 0.7}, {"t_end": 0}, {"rtol": 1e-2},
                                {"seed": 1.5}, {"t_end": float("nan")}])
def test_invalid_arguments(kw):
    with pytest.raises(hv.InvalidInput):
        simulate("sinh", 9, "-0.3", **kw)


@pytest.mark.slow
def test_long_run_stays_bounded():
    s = simulate("sinh", 9, "-0.3", t_end=2000.0)
    assert s.growth < 1.5 and s.energy_drift < s.drift_tolerance


def test_simulation_reports_r():
    import mpmath as mp
    s = simulate("sinh", 9, "-0.3", t_end=1.0)
    assert "r" in s.to_dict()
    with mp.workdps(50):
        assert abs(s.r - 2 * mp.atanh(mp.sqrt(mp.mpf(3) / 10))) < mp.mpf("1e-40")


# ---- the perturbed state and the trajectory must stay in the phase space ----
@pytest.mark.parametrize("args,kw", [
    (("tanh", 8, "0.99"), {"perturb": 0.1, "t_end": 1.0}),             # max |z0| = 1.094: outside the Poincare disk
    (("tanh", 8, "0.99"), {"perturb": 0.3, "t_end": 1.0, "seed": 1}),
    (("sinh", 8, "-0.9"), {"perturb": 0.1, "t_end": 1.0}),             # |z0| = 1.1 > 1/sqrt(0.9) = 1.054
])
def test_initial_state_outside_the_phase_space_is_rejected(args, kw):
    with pytest.raises(hv.InvalidInput, match="phase space"):
        simulate(*args, **kw)


def test_initial_state_inside_the_phase_space_is_accepted():
    s = simulate("tanh", 8, "0.5", perturb=0.1, t_end=1.0)
    assert np.all(np.abs(s.trajectory) < 1)


def _fake_solver(bad):
    import scipy.integrate

    def fake(fun, t_span, y0, t_eval, **kw):
        sol = scipy.integrate.solve_ivp(fun, t_span, y0, t_eval=t_eval, **kw)
        n = len(y0) // 2
        if bad == "outside":
            sol.y[0, -1] = 5.0                       # vortex 0 leaves the disk at the last sample
        else:
            sol.y[:n, -1] = sol.y[0, -1]             # every vortex at the same point: H is not finite
            sol.y[n:, -1] = sol.y[n, -1]
        return sol
    return fake


@pytest.mark.parametrize("bad", ["outside", "collision"])
def test_trajectory_leaving_the_phase_space_or_non_finite_drift_is_an_error(monkeypatch, bad):
    import sys
    S = sys.modules["havelock.simulate"]
    monkeypatch.setattr(S, "solve_ivp", _fake_solver(bad))
    with pytest.raises(hv.OutOfValidityRange):
        simulate("tanh", 8, "0.5", t_end=1.0)
