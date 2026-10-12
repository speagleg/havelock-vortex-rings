"""Certificate dataclass: JSON export, Lean export (tanh), and verification."""
from __future__ import annotations

import json
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from ..certify import certify_threshold
from ..errors import CertificationFailed, InvalidInput
from ..limits import CERT_MAX_DIGITS, VERIFY_MAX_N
from ..models import get_model
from ..status import Status, Threshold, jsonable, make_provenance
from .verify import SCHEMA, VerifyReport, verify_dict

#: Lean export range (tanh): N = 8..16 is emitted and tested; N <= 7 has no threshold, and from N = 17 on the
#: emitter fails (f' has a root in (1/16, 1)).
LEAN_EXPORT_MIN_N, LEAN_EXPORT_MAX_N = 8, 16


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


def _no_duplicate_keys(pairs):
    """json object_pairs_hook: a repeated key would let a later value silently replace an earlier one."""
    out = {}
    for k, v in pairs:
        if k in out:
            raise ValueError(f"duplicate key {k!r} in a JSON object")
        out[k] = v
    return out


@dataclass(frozen=True)
class Certificate:
    """A machine-checkable threshold certificate (format: docs/certificates.md)."""

    data: dict

    def to_json(self, indent: int | None = 1) -> str:
        return json.dumps(self.data, indent=indent)

    def to_lean(self, path) -> Path:
        """Write a Lean 4 module proving the threshold (tanh model, N = 8..16 only; it imports
        Havelock.Tanh.General). Other tanh N raise CertificationFailed."""
        if self.data["model"] != "tanh":
            raise InvalidInput("Lean export is available for the tanh model only; MT-family thresholds are proved "
                               "in the Lean library itself (Havelock.Conj13.bn_spec, Havelock.Sinh)")
        N = self.data["N"]
        if isinstance(N, bool) or not isinstance(N, int) or not LEAN_EXPORT_MIN_N <= N <= LEAN_EXPORT_MAX_N:
            raise CertificationFailed(f"Lean export is available for tanh N = {LEAN_EXPORT_MIN_N}..{LEAN_EXPORT_MAX_N} "
                                      f"only (got N = {N!r}); other N raise CertificationFailed (N <= 7 has no "
                                      "threshold, and from N = 17 the emitter's proof template, which needs f' to "
                                      "have no root in (1/16, 1), does not apply)")
        import flint
        import mpmath as mp

        from ..blocks import numerators
        from ..certify import canonical_factors
        from .lean import emit_tanh_module

        f = canonical_factors(flint.fmpz_poly(self.data["polynomial"]))[self.data["factor_index"]]
        lo, hi = (Fraction(s) for s in self.data["interval"])
        with mp.workdps(60):
            xs = (mp.mpf(lo.numerator) / lo.denominator + mp.mpf(hi.numerator) / hi.denominator) / 2

        def numer(N_, m, kind):
            nm = numerators("tanh", N_, m)
            return [int(c) for c in {"A": nm.A, "B": nm.B}[kind].coeffs()]

        try:
            text = emit_tanh_module(N, numer, [int(c) for c in f.coeffs()], xs)
        except ValueError as e:
            raise CertificationFailed(f"Lean emitter: {e}") from None
        path = Path(path)
        path.write_text(text, encoding="utf-8")
        return path

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
                data = json.loads(text, object_pairs_hook=_no_duplicate_keys)
            except (ValueError, RecursionError) as e:   # JSONDecodeError, duplicate key, int-digit limit, deep nesting
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
    cap = 10 ** CERT_MAX_DIGITS
    if any(abs(q.numerator) >= cap or q.denominator >= cap for q in (*c.interval, *c.points)):
        raise CertificationFailed(f"an endpoint has more than {CERT_MAX_DIGITS} digits, above the verifier's cap "
                                  "(hv.limits()['certificate_max_digits']); use a smaller starting precision")
    return Certificate(_to_data(t, c, prov))
