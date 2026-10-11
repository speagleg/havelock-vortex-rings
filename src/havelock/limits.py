"""Documented limits (hv.limits())."""
from __future__ import annotations

import copy

#: decimal digits carried by every mpmath value the package returns
WORK_DPS = 50
#: Arb working precisions in bits, tried in order until the certification succeeds
PRECISION_LADDER = (256, 512, 1024, 2048, 4096)
#: simulate(): the reported drift tolerance is DRIFT_FACTOR * rtol
DRIFT_FACTOR = 1e3
#: simulate(): number of output samples
SIM_SAMPLES = 2001
#: nonlinear(): level-set engine settings (decimal digits, number of eps values, eps = k/EPS_DENOM)
LEVELSET = {"dps": 40, "steps": 6, "eps_denominator": 200}
#: largest prime tried when searching for a mod-p irreducibility witness
WITNESS_MAX_PRIME = 20000
#: largest N for which certificate() issues a certificate (the independent verifier refuses larger N)
VERIFY_MAX_N = 100000

_LIMITS = {
    "work_dps": WORK_DPS,
    "precision_ladder_bits": list(PRECISION_LADDER),
    "simulation": {"drift_tolerance": "1e3 * rtol", "default_rtol": 1e-12, "samples": SIM_SAMPLES},
    "nonlinear_levelset": dict(LEVELSET),
    "irreducibility_witness_max_prime": WITNESS_MAX_PRIME,
    "certificate_max_N": VERIFY_MAX_N,
}


def limits() -> dict:
    """A copy of the documented limits: size caps, precision ladder and tolerances."""
    return copy.deepcopy(_LIMITS)
