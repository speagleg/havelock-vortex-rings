"""Command-line interface: havelock threshold|stability|nonlinear|simulate|cert|verify|limits."""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import mpmath as mp

from . import __version__
from .certificates import Certificate, certificate
from .errors import HavelockError
from .models import get_model
from .limits import limits
from .nonlinear import nonlinear
from .simulate import simulate
from .stability import stability
from .status import NoThreshold, jsonable
from .threshold import threshold


def _num(x, digits=25):
    return mp.nstr(x, digits) if isinstance(x, mp.mpf) else str(x)


def _print_threshold(t):
    if isinstance(t, NoThreshold):
        print(f"{t.model}, N = {t.N}: {t.verdict} ({t.geometry})")
    else:
        p = get_model(t.model).param
        if t.stable_side is None:
            side = "degenerate one-point view: the ring sits at the bifurcation point"
        elif t.stable_side == "below":
            side = f"stable for {p} < threshold"
        elif t.N % 2:
            side = f"spectrally stable for {p} > threshold; Krein window just above it"
        else:
            side = f"stable for {p} > threshold"
        print(f"{t.model}, N = {t.N}: threshold {p} = {_num(t.value)} ({side}), critical mode "
              f"{t.mode}, geometry {t.geometry}")
        extra = [f"{k} = {_num(getattr(t, k))}" for k in ("xi", "r", "cosh_r")
                 if getattr(t, k) is not None and k != p]
        if extra:
            print("  " + ", ".join(extra))
        print(f"  exact: {t.exact}")
    print(f"  status: {t.status}" + (f" (fallback from {t.fallback_from})" if t.fallback_from else ""))
    print(f"  basis: {t.basis}")


def _emit(args, result, human):
    if args.json:
        print(json.dumps(result if isinstance(result, dict) else result.to_dict(), indent=1))
    else:
        human()


def _cmd_threshold(a):
    t = threshold(a.model, a.N, strict=not a.no_strict, precision=a.precision)
    _emit(a, t, lambda: _print_threshold(t))


def _cmd_stability(a):
    s = stability(a.model, a.N, a.value)

    def human():
        print(f"{s.model}, N = {s.N}, value {s.value}: {s.verdict}  [{s.status}]")
        for k in s.blocks:
            print(f"  m = {k.m}: a = {float(k.a):.6g}, b = {float(k.b):.6g}, g = {float(k.g):.6g}")
        print(f"  basis: {s.basis}")

    _emit(a, s, human)


def _cmd_nonlinear(a):
    r = nonlinear(a.model, a.N)
    _emit(a, r, lambda: print(f"{r.model}, N = {r.N}: c4 = {_num(r.c4)} at {_num(r.value)}, verdict "
                              f"{r.verdict}  [{r.status}]\n  basis: {r.basis}"))


def _cmd_simulate(a):
    s = simulate(a.model, a.N, a.value, perturb=a.perturb, t_end=a.t_end, rtol=a.rtol, seed=a.seed)
    summary = {k: v for k, v in s.to_dict().items() if k not in ("times", "trajectory")}
    _emit(a, summary, lambda: print(
        f"{s.model}, N = {s.N}, value {s.value}: growth {s.growth:.4g} (max/initial deviation; max "
        f"{s.max_deviation:.3e}), "
        f"energy drift {s.energy_drift:.2e}, momentum drift {s.momentum_drift:.2e} "
        f"(tolerance {s.drift_tolerance:.1e})  [NUMERICAL]"))


def _cmd_cert(a):
    t = threshold(a.model, a.N)
    c = certificate(t)
    if a.json == "-":
        print(c.to_json())
    elif a.json:
        Path(a.json).write_text(c.to_json())
        print(f"wrote {a.json}")
    if a.lean:
        p = c.to_lean(a.lean)
        print(f"wrote {p}")
        if a.check:
            proj = Path(a.lean_project).resolve()
            r = subprocess.run(["lake", "env", "lean", str(Path(p).resolve())], cwd=proj)
            if r.returncode != 0:
                print("Lean check FAILED", file=sys.stderr)
                return 1
            print("Lean check passed")
    if not a.json and not a.lean:
        print(c.to_json())
    return 0


def _cmd_verify(a):
    try:
        rep = Certificate.verify(Path(a.cert))
    except (OSError, ValueError) as e:     # missing/unreadable file, non-UTF-8 bytes (UnicodeDecodeError)
        print(f"havelock: {type(e).__name__}: {e}", file=sys.stderr)
        return 2
    if a.json:
        print(json.dumps(jsonable(rep), indent=1))
    else:
        print(("VALID" if rep.ok else "INVALID") + f" ({rep.checks} checks)")
        for f in rep.failures:
            print("  - " + f)
    return 0 if rep.ok else 1


def _cmd_limits(a):
    lim = limits()
    _emit(a, lim, lambda: print(json.dumps(lim, indent=1)))


def _add_value(s):
    """VALUE positionally, or --value VALUE; exactly one. Negative values (-0.3, -1/7, -1e-3) work either way."""
    s.add_argument("value", nargs="?", default=None, help="parameter value, e.g. -0.3, -1/7 or -1e-3")
    s.add_argument("--value", dest="value_opt", default=None, metavar="VALUE",
                   help="parameter value (alternative to the positional VALUE)")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="havelock", description="Stability of regular N-vortex rings")
    p.add_argument("--version", action="version", version=f"havelock {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    def add(name, fn, help_):
        s = sub.add_parser(name, help=help_)
        s.set_defaults(fn=fn, _sub=s)
        return s

    s = add("threshold", _cmd_threshold, "stability threshold of the N-ring")
    s.add_argument("model"), s.add_argument("N", type=int)
    s.add_argument("--no-strict", action="store_true", help="allow a NUMERICAL fallback instead of raising")
    s.add_argument("--precision", type=int, default=None, help="starting Arb precision in bits")
    s.add_argument("--json", action="store_true")
    s = add("stability", _cmd_stability, "verdict at one parameter value")
    s.add_argument("model"), s.add_argument("N", type=int)
    _add_value(s)
    s.add_argument("--json", action="store_true")
    s = add("nonlinear", _cmd_nonlinear, "c4 and the nonlinear verdict at the bifurcation point")
    s.add_argument("model"), s.add_argument("N", type=int)
    s.add_argument("--json", action="store_true")
    s = add("simulate", _cmd_simulate, "integrate the perturbed ring (NUMERICAL)")
    s.add_argument("model"), s.add_argument("N", type=int)
    _add_value(s)
    s.add_argument("--perturb", type=float, default=1e-4), s.add_argument("--t-end", type=float, default=100.0)
    s.add_argument("--rtol", type=float, default=1e-12), s.add_argument("--seed", type=int, default=None)
    s.add_argument("--json", action="store_true")
    s = add("cert", _cmd_cert, "threshold certificate (JSON, optionally Lean)")
    s.add_argument("model"), s.add_argument("N", type=int)
    s.add_argument("--json", nargs="?", const="-", default=None, metavar="PATH",
                   help="write the certificate to PATH (stdout if PATH is omitted)")
    s.add_argument("--lean", default=None, metavar="PATH",
                   help="write a Lean 4 module proving the threshold (tanh N = 8..16 only; other N raise "
                        "CertificationFailed and other models InvalidInput, both exit 2)")
    s.add_argument("--check", action="store_true", help="compile the Lean module with `lake env lean`")
    s.add_argument("--lean-project", default="lean", help="Lean project directory (default: ./lean)")
    s = add("verify", _cmd_verify, "re-check a certificate exactly")
    s.add_argument("cert"), s.add_argument("--json", action="store_true")
    s = add("limits", _cmd_limits, "documented limits")
    s.add_argument("--json", action="store_true")
    return p


# a negative number: integer, decimal, scientific notation (-1e-3, -.3, -2.5E+4) or rational (-1/7)
_NEG_VALUE = re.compile(r"^-(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$|^-\d+/\d+$")
# options of `stability` / `simulate` that take a value (argparse would read a following '-1e-3' as an option)
_VALUE_OPTIONS = ("--value", "--perturb", "--t-end", "--rtol", "--seed")


def _negative_values_to_option(argv):
    """argparse reads '-1/7' or '-1e-3' as an unknown option. A negative number right after an option that takes a
    value is attached to it ('--t-end -1' -> '--t-end=-1'); a negative number in the VALUE position is passed on as
    --value=... (unless --value is also given, which argparse then reports)."""
    if len(argv) < 2 or argv[0] not in ("stability", "simulate"):
        return argv
    has_value_opt = any(x == "--value" or x.startswith("--value=") for x in argv)
    out = []
    for x in argv:
        if _NEG_VALUE.match(x) and out and out[-1] in _VALUE_OPTIONS:
            out[-1] = f"{out[-1]}={x}"
        elif _NEG_VALUE.match(x) and not has_value_opt:
            out.append(f"--value={x}")
        else:
            out.append(x)
    return out


def main(argv=None) -> int:
    argv = _negative_values_to_option(list(sys.argv[1:] if argv is None else argv))
    args = build_parser().parse_args(argv)
    if hasattr(args, "value_opt"):
        if (args.value is None) == (args.value_opt is None):
            args._sub.error("give exactly one of the positional VALUE or --value VALUE")
        if args.value is None:
            args.value = args.value_opt
    try:
        rc = args.fn(args)
    except (HavelockError, OSError) as e:
        print(f"havelock: {type(e).__name__}: {e}", file=sys.stderr)
        return 2
    return rc or 0


if __name__ == "__main__":
    sys.exit(main())
