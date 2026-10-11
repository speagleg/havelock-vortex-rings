# havelock quickstart

Install (Python >= 3.10):

    pip install git+<repo-url>          # or, from a checkout: pip install -e .

Every answer carries a status: PROVED (a paper theorem covers this exact case), CERTIFIED (exact or ball arithmetic
at run time) or NUMERICAL (floating point), plus a `basis` string that says why.

## Command line

```console
$ havelock threshold sinh 12
$ havelock threshold sphere 9
$ havelock stability sinh 9 -0.3
$ havelock stability sinh 9 --value=-1/7
$ havelock nonlinear mt 13
$ havelock simulate sinh 9 -0.05 --t-end 50
$ havelock cert sinh 8 --json cert.json
$ havelock verify cert.json
$ havelock limits
```

Add `--json` to any command for machine-readable output.

`stability` and `simulate` take the parameter value either positionally or as `--value VALUE` (exactly one). A negative
decimal such as `-0.3` works positionally, but argparse cannot read a negative rational such as `-1/7` there: write
`--value=-1/7`.

Exit codes: 0 success; 1 an invalid certificate (`verify`) or a failed Lean check; 2 an error (a `HavelockError`, a
missing, unreadable or non-UTF-8 file) reported on stderr without a traceback.

## Python

```python
import havelock as hv

t = hv.threshold("sinh", 12)      # Threshold: t.value, t.exact, t.mode, t.r, t.status, t.basis
s = hv.stability("sinh", 9, "-0.3")  # verdict "EM-stable"; values are used exactly ("-0.3" is -3/10)
c = hv.nonlinear("mt", 13)        # c.c4 == Fraction(675, 104), c.verdict == "stable"
sim = hv.simulate("sinh", 9, "-0.05")  # sim.growth, sim.energy_drift, sim.drift_tolerance (NUMERICAL)
cert = hv.certificate(t)          # cert.to_json(); hv.Certificate.verify(cert.to_json()).ok
```

Floats are converted to their exact binary value; pass strings such as `"-0.3"` or `"1/7"` for decimal or rational
values.
