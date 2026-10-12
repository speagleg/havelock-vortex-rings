# havelock quickstart

Install (Python >= 3.10):

    pip install git+https://github.com/speagleg/havelock-vortex-rings

or, from a checkout, into a virtual environment:

    python3 -m venv .venv && .venv/bin/pip install -e ".[dev]"
    source .venv/bin/activate          # so that `havelock` is on the PATH

Every answer carries a status: PROVED (a paper theorem covers this exact case), CERTIFIED (exact or ball arithmetic
at run time) or NUMERICAL (floating point), plus a `basis` string that says why.

## Command line

```console
$ havelock threshold sinh 12
sinh, N = 12: threshold lam_r2 = -0.2024108640318572506321909 (stable for lam_r2 < threshold), critical mode 6, geometry H2
  xi = 0.2024108640318572506321909, r = 0.9691521676443146127214524, cosh_r = 1.50755672288881811323406
  exact: (-18 + 5*sqrt(11))/7
  status: PROVED
  basis: cosh^2 r* = q(N), q = (N-3)/4 (N odd), (N-2)^2/(4(N-1)) (N even) (sinh/note.tex, Section 3; MT (2.11)); Lean Havelock.Sinh
```

Each model has its own native parameter: `lam_r2` = λr₀² for the Montaldi–Tokieda family (`mt`, `sphere`, `plane`,
`sinh`; the ring is stable *below* the threshold), and `xi` = tanh²(r/2) for the tanh model (stable *above* it). The
geodesic radius `r` and `cosh_r` are printed for every model so the models can be compared.

More examples:

```console
$ havelock threshold sphere 9                  # "unstable at every radius": b_9 lies outside the sphere's range
$ havelock stability sinh 9 -0.3               # EM-stable [PROVED]
$ havelock stability sinh 9 -1/7               # rationals are used exactly
$ havelock nonlinear mt 13                     # c4 = 675/104 > 0: stable at the bifurcation point
$ havelock simulate sinh 9 -0.3 --t-end 50     # stable side (threshold -0.101): growth ~1
$ havelock simulate sinh 9 -0.05 --t-end 50    # unstable side: growth ~10^4
$ havelock cert sinh 8 --json cert.json
$ havelock verify cert.json                    # VALID (… checks)
$ havelock threshold tanh 12                   # CERTIFIED at run time (seconds)
$ havelock stability tanh 9 0.11287108         # spectrally-stable-only: inside the N = 9 Krein window
$ havelock cert tanh 8 --json cert8.json --lean CertN8.lean
$ havelock limits
```

Notes on the examples:

- `simulate` reports `growth` = max deviation / initial deviation of the inter-vortex distances. Values near 1 mean
  the perturbation stayed small; values of 10² or more mean it grew. The simulation is NUMERICAL evidence only; the
  verdicts come from `stability` and `nonlinear`.
- `stability tanh 9 0.11287108` evaluates at a point about 2·10⁻⁸ above the spectral threshold ξ*(9), so the critical
  block has a tiny a₄ ≈ 3·10⁻⁷. All a_m b_m are positive (linearly stable), but a₄b₄ − g₄² < 0, so the ring sits in
  the Krein window, where nonlinear stability is undecided. The value is echoed as the exact rational it denotes.
- Certified tanh thresholds take from a few seconds (small N) to about 20 s (N = 64, the size cap).

Add `--json` to any command for machine-readable output.

`stability` and `simulate` take the parameter value either positionally or as `--value VALUE` (exactly one).
Negative decimals (`-0.3`) and negative rationals (`-1/7`) both work positionally.

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
x = hv.threshold("tanh", 12)      # CERTIFIED; x.exact is the minimal polynomial with an isolating interval
x8 = hv.certificate(hv.threshold("tanh", 8))
x8.to_lean("CertN8.lean")         # compile with: cd lean && lake env lean ../CertN8.lean
```

Floats are converted to their exact binary value; pass strings such as `"-0.3"` or `"1/7"` for decimal or rational
values.
