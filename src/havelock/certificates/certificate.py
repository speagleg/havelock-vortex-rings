"""Certificate dataclass: JSON export, Lean export (tanh), and verification."""
from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

from ..certify import certify_threshold
from ..errors import CertificationFailed, InvalidInput
from ..limits import VERIFY_MAX_N
from ..models import get_model
from ..status import Status, Threshold, jsonable, make_provenance
from .verify import SCHEMA, VerifyReport, verify_dict


def _to_data(t: Threshold, c, prov) -> dict:
    return {
        "schema": SCHEMA, "model": t.model, "N": t.N, "mode": c.mode,
        "polynomial": [int(a) for a in c.polynomial.coeffs()],
        "factor_index": c.factor_index,
        "interval": [str(c.interval[0]), str(c.interval[1])],
        "endpoint_signs": list(c.endpoint_signs),
        "sturm": c.sturm,
        "points": [str(p) for p in c.points],
        "mode_signs": {str(m): v for m, v in c.mode_signs.items()},
        "stable_side": {"side": c.stable_side["side"], "interval": [str(a) for a in c.stable_side["interval"]],
                        "root_counts": {str(m): v for m, v in c.stable_side["root_counts"].items()}},
        "irreducibility_witness": c.irreducibility_witness,
        "interval_clear": dict(c.interval_clear),
        "provenance": jsonable(prov),
    }


@dataclass(frozen=True)
class Certificate:
    """A machine-checkable threshold certificate (format: docs/certificates.md)."""

    data: dict

    def to_json(self, indent: int | None = 1) -> str:
        return json.dumps(self.data, indent=indent)

    def to_lean(self, path) -> Path:
        """Write a Lean 4 module proving the threshold (tanh model only; it imports Havelock.Tanh.General)."""
        raise InvalidInput("Lean export is available for the tanh model only; MT-family thresholds are proved in "
                           "the Lean library itself (Havelock.Conj13.bn_spec, Havelock.Sinh)")

    @classmethod
    def verify(cls, data) -> VerifyReport:
        """Re-check a certificate (a Certificate, a dict, a JSON string, or a path to a JSON file) exactly.
        A string starting with "{" or "[" is JSON text; any other string is a path."""
        if isinstance(data, Certificate):
            data = data.data
        elif isinstance(data, (str, Path)):
            is_text = isinstance(data, str) and data.lstrip()[:1] in ("{", "[")
            text = data if is_text else Path(data).read_text()
            try:
                data = json.loads(text)
            except (ValueError, RecursionError) as e:   # JSONDecodeError, int-digit limit, deep nesting
                return VerifyReport(False, [f"not valid JSON: {e}"], 0)
        if not isinstance(data, dict):
            return VerifyReport(False, ["certificate must be a JSON object"], 0)
        return verify_dict(data)


def certificate(t) -> Certificate:
    """Certificate for a Threshold returned by havelock.threshold (any model)."""
    if not isinstance(t, Threshold):
        raise InvalidInput("certificate() needs a Threshold (a NoThreshold has nothing to certify)")
    if t.status == Status.NUMERICAL or t.fallback_from is not None:
        raise CertificationFailed("no certificate for a NUMERICAL result; use strict=True within the size cap")
    if t.N > VERIFY_MAX_N:
        raise CertificationFailed(f"certificates are limited to N <= {VERIFY_MAX_N} so they can be verified "
                                  f"(got N={t.N})")
    t0 = time.perf_counter()
    mdl = get_model(t.model)
    c = t._certification
    if c is None:
        c = certify_threshold(mdl.family, t.N, mdl.family_domain, mdl.stable_side, witness=True)
    prov = make_provenance("certificate", {"model": t.model, "N": t.N}, t0)
    return Certificate(_to_data(t, c, prov))
