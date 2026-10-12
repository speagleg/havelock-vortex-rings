# `havelock` v0.1: design spec

Date: 2026-10-11. Status: implemented as `havelock` 0.1.0 (2026-10-12). Kept unedited as the design record; the user
documentation is [docs/quickstart.md](../quickstart.md), [docs/models.md](../models.md) and
[docs/certificates.md](../certificates.md), and where the built package differs from this spec, those documents
and the code are authoritative.

## 1. Purpose

The package is a Python library and command-line tool for researchers in point-vortex dynamics. Given a vortex model and a
number N of vortices, it answers these questions about the regular N-vortex ring:

- Where does the ring lose stability?
- Which Fourier mode goes first?
- Is the ring nonlinearly stable at the bifurcation point?
- What does the motion look like?

Every answer carries an honest label (PROVED, CERTIFIED or NUMERICAL) and the reason for that label. Certified answers can
be exported as machine-checkable certificates, including a Lean file.

The package turns the two research results in this repo into a tool others can query for any N, not only the tabulated
ones:

- the Conjecture 13 paper (`conj13/paper`);
- the tanh-model note (`tanh/note_tanh.tex`).

**Success criteria:**

1. A researcher can run `pip install git+<repo>` and then `havelock threshold tanh 12`, and get a certified threshold
   together with its basis.
2. Every number the package prints that also appears in the papers matches the paper. Tests enforce this.
3. No result is labelled PROVED unless a paper theorem covers that exact case.

## 2. Global constraints

- Python ≥ 3.10. Runtime dependencies: `mpmath`, `sympy`, `numpy`, `scipy`, `python-flint` (the minimum versions are
  fixed in the plan from the versions tested). No other runtime dependencies. Lean is optional and never required at
  runtime.
- The package lives in `src/havelock/` in this repo, built with `pyproject.toml` (src layout) and exposing the
  console script `havelock`. The version is 0.1.0.
- The package does not depend on `bolzakit`. It borrows bolzakit's conventions (status labels, `fallback_from`,
  exception style) but imports nothing from it.
- The research scripts (`conj13/`, `tanh/`, `sinh/`) stay where they are, as the record behind the papers. The
  package ports and cleans up their logic; it does not import them.
- Nothing is published (no PyPI, no public repo) without the user's explicit decision.
- **Gate:** the tanh model stays out of any release until the Boatto (2008) discrepancy is resolved. The investigation
  is in `tanh/boatto_check/`. If our tanh blocks are wrong, the tanh model is fixed or withdrawn first.

## 3. Models and parameters

| Model name | Hamiltonian | Native parameter | Geometry |
|---|---|---|---|
| `"mt"` | Montaldi–Tokieda family [MT (1.12)] | `lam_r2` = λr₀² (stable for λr₀² < b_n) | sphere (>0), plane (=0), H² (<0) |
| `"sphere"`, `"plane"`, `"sinh"` | the same family, restricted to one geometry | `lam_r2` within that geometry's range | as named |
| `"tanh"` | H = −Σ ln tanh(d/2) on H² (Kimura's Hamiltonian) | `xi` = tanh²(r/2) ∈ (0,1) | H² |

- Every result also reports the geodesic ring radius (with curvature normalised to ±1, or the chart radius for the
  plane), so the models can be compared.
- A geometry-restricted view returns the family's threshold if b_n lies in that geometry. Otherwise it returns a
  `NoThreshold` result stating "stable at every radius" or "unstable at every radius", with its status. For example,
  `threshold("sphere", 9)` returns "unstable at every radius".

## 4. Public API

All public functions return a subclass of `Result`. No public function returns a bare float.

```python
import havelock as hv

t = hv.threshold(model, N, *, strict=True, precision=None)   # -> Threshold | NoThreshold
#   .value (native parameter), .xi / .r / .cosh_r where meaningful, .exact (closed form: a sympy expression, or a
#   minimal polynomial with an isolating interval), .mode (critical Fourier mode), .status, .basis, .provenance,
#   .fallback_from

hv.stability(model, N, value)        # -> Stability: verdict in {"EM-stable", "spectrally-stable-only", "unstable"},
                                     #    the per-mode blocks (a_m, b_m, g_m) and eigenvalues, status, basis
hv.nonlinear(model, N)               # -> Nonlinear: c4 at the bifurcation point, verdict, status, basis
hv.simulate(model, N, value, *, perturb=1e-4, t_end=100.0, rtol=1e-12, seed=None)
                                     # -> Simulation: times, trajectory, max deviation, relative energy and momentum
                                     #    drift; always NUMERICAL
hv.certificate(t)                    # -> Certificate, with .to_json(), .to_lean(path), and classmethod .verify(json)
hv.limits()                          # -> dict of documented limits (size cap, tolerances)
```

**CLI.** `havelock threshold MODEL N`, `havelock stability MODEL N VALUE`, `havelock nonlinear MODEL N`,
`havelock simulate MODEL N VALUE [...]`, `havelock cert MODEL N [--json PATH] [--lean PATH [--check]]`,
`havelock verify CERT.json`. Every command accepts `--json`, which prints machine-readable output.

## 5. Status rule

| Case | Status | Basis string cites |
|---|---|---|
| MT family thresholds b_n, critical mode, n ≥ 4 | PROVED | paper (eq:bn); Lean `Havelock.Conj13` |
| MT family c₄, n ≥ 4 | PROVED | paper Thms even/odd; Lean (the algebra downstream of the series coefficients) |
| MT family, n = 3 | — | raises `OutOfValidityRange` (c₄ = 0; this case rests on MT's own argument) |
| sinh-view threshold for any N | PROVED | cosh² r* = q(N), `sinh/note.tex`; Lean `Havelock.Sinh` |
| tanh block formulas, any N | PROVED (paper) | note Thm 2 |
| tanh threshold, any N | CERTIFIED | computed at runtime: Arb root isolation plus exact Sturm count; the basis adds "Lean-checked" for N ≤ 16 |
| tanh c₄ | NUMERICAL | level-set engine; no closed form is known |
| `stability` at a point | PROVED/CERTIFIED if the verdict follows from exact sign checks, otherwise NUMERICAL | |
| `simulate` | NUMERICAL | always |

If precision or size limits force a fallback from CERTIFIED to NUMERICAL, the result records it in `fallback_from`.
This happens only with `strict=False`; with `strict=True` the call raises.

## 6. Components

```
src/havelock/
  __init__.py      public API re-exports
  status.py        Result base class; Status enum; provenance (version, inputs, timing)
  errors.py        HavelockError > InvalidInput, OutOfValidityRange, CertificationFailed
  models/          one module per family: mt.py, tanh.py. Each holds the kernel, ring geometry and parameter conversions
  blocks.py        exact Fourier blocks a_m, b_m, g_m (sympy rational functions of the parameter); mpmath fast path
  threshold.py     builds the threshold equation, solves it, and calls certify
  certify.py       python-flint: square-free and irreducible factorisation over Q; Arb root isolation; exact Sturm
                   sign checks; mode-sign checks
  nonlinear.py     c4: MT closed forms; numerical level-set engine for tanh (ported from conj13/engine.py)
  simulate.py      integrator ported from conj13/integrate.py, with energy and momentum drift diagnostics
  certificates/    Certificate dataclass, JSON schema (versioned), exact verifier, Lean emitter
                   (ported from tanh/lean_gen)
  cli.py
```

**Data flow:** model → exact blocks → threshold polynomial (cleared numerator of the critical block) → factor over Q →
isolate the relevant root of the irreducible factor → certify the sign pattern of every mode on each side → Result.

**Lean library additions:**

- A mod-p irreducibility certificate for f_N, N ≤ 16. PARI has found small primes for each N
  (see `tanh/irreducibility/`). If this proves impractical in Mathlib, it is dropped from v0.1 and the degrees stay
  labelled "CERTIFIED by computer algebra"; it does not block the release.
- The Lean emitter targets `Havelock.Tanh.General`, whatever module structure the CI split settles on.

## 7. Certificate format

The JSON has these fields:

- `schema` (a version);
- `model`, `N`;
- `polynomial`, as integer coefficients with the constant term first;
- `factor_index`: which irreducible factor holds the threshold;
- `interval`, as rational endpoints [p/q, p'/q'];
- `endpoint_signs`;
- `sturm`: the number of real roots in the interval, which must be 1;
- `mode_signs`: for each Fourier mode, the sign of a_m and b_m (and of a_m b_m − g_m² when it is relevant) at two
  rational points, one just below the interval and one just above;
- `irreducibility_witness` (optional): a prime p and the factorisation degrees mod p;
- `provenance`.

`havelock verify` re-checks every field using exact integer and rational arithmetic only. It uses no floating point,
and it never needs the code that produced the certificate.

## 8. Errors and limits

- **Invalid input:** N < 2, a parameter outside its domain, or an unknown model raises `InvalidInput` before any
  computation starts.
- **Certification failure:** precision is raised geometrically up to a cap. If the result still cannot be certified,
  the call raises `CertificationFailed`, or returns NUMERICAL with `fallback_from` set when `strict=False`.
- **Out of validity:** where an answer would be wrong rather than merely slow, the call raises `OutOfValidityRange`
  with an explanation. This covers MT n = 3, and any case outside a closed form's proven range.
- **Size cap:** exact tanh certification has a cap on N, set from benchmarks in the plan (expected to be somewhere
  around N ≈ 40–80). Above the cap the call follows the strict/non-strict rule above.
- **Visible limits:** `hv.limits()` reports the cap, the simulation tolerances and the precision ladder.

## 9. Testing

1. **Paper values.**
   - the c₄ table, n = 4–25;
   - the MT b_n closed forms;
   - the sinh q(N) closed form;
   - the tanh table, N = 8–16: ξ* to 25 digits, degrees, leading coefficients and critical modes.
2. **Independent cross-check.** For several N and parameter values in every model, the exact block eigenvalues must equal
   those of a brute-force finite-difference Hessian of the full 2N-dimensional system, with momentum handled
   correctly. This is the guard against a Boatto-type modelling error.
3. **Certificate round-trip.** Generating a certificate and then verifying it must pass. A corrupted certificate (a
   shifted endpoint, a changed coefficient, or a flipped mode sign) must fail.
4. **Lean export.** In CI, generate the N = 8 certificate and compile it against the Lean library.
5. **Simulation.**
   - An EM-stable ring (sinh N = 9 well past its threshold) stays within a bounded deviation.
   - A ring below its threshold shows growth.
   - Energy and momentum drift stay below the reported tolerance.
6. **CLI.** Smoke tests for every subcommand; `--json` output must parse.
7. **Research reproduction.** The ported code reproduces the recorded outputs that the reproducibility pass
   identifies (`REPRODUCE.md`).
8. **Speed tiers.** A slow tier (large-N certification, long simulations) is marked and runs only on manual dispatch.

CI runs Python 3.10–3.13 on Linux. The Lean workflow runs alongside it.

## 10. Documentation

- `README.md` gets a short "Using the package" section.
- `docs/quickstart.md`.
- `docs/models.md`: one section per model, saying what is proved, certified or numerical, with citations to the papers.
- `docs/certificates.md`: the format and how to verify a certificate.

## 11. Out of scope for v0.1

- higher-genus surfaces (bolzakit covers those);
- unequal vortex strengths;
- rings not centred on a symmetry point;
- physical units or BEC parameters;
- a closed form for tanh c₄;
- publishing.
