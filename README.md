# Vortex-ring stability beyond the plane

[![lean](https://github.com/speagleg/havelock-vortex-rings/actions/workflows/lean.yml/badge.svg)](https://github.com/speagleg/havelock-vortex-rings/actions/workflows/lean.yml)
[![python](https://github.com/speagleg/havelock-vortex-rings/actions/workflows/python.yml/badge.svg)](https://github.com/speagleg/havelock-vortex-rings/actions/workflows/python.yml)
[![licence: MIT (code), CC BY 4.0 (papers)](https://img.shields.io/badge/licence-MIT%20%2F%20CC%20BY%204.0-blue)](#licence)

Gordon Speagle, Independent Researcher, Boise, ID, USA.

The regular ring of N identical point vortices is a classical object: Thomson posed the question in 1883, Havelock
answered it in 1931: in the plane the ring is linearly stable exactly when N ≤ 7 (N = 7 only marginally; its
nonlinear stability was settled later). This repository studies the same question
on curved surfaces: the sphere, the plane and the hyperbolic plane. It contains two research results, the proofs and
computations behind them, a Lean 4 check of their exact algebra, and `havelock`, a Python package and command-line tool
that answers the same questions for any N.

    pip install git+https://github.com/speagleg/havelock-vortex-rings
    havelock threshold tanh 12

> **Status.** These are drafts. They have been checked by exact and interval computer algebra, by an independent
> from-scratch numerical implementation, and by a Lean build in CI. **No human specialist has reviewed them yet.**
> Each claim is labelled PROVED, CERTIFIED or NUMERICAL in [docs/results.md](docs/results.md). The work was developed
> with AI assistance (Claude); each paper's acknowledgement says how.

## The results

**1. A proof of Montaldi–Tokieda's Conjecture 13 for n ≥ 4.**

Montaldi and Tokieda (2013) conjectured that the ring stays stable *at* the curvature where it bifurcates. They treat
even n by a closed-form symbolic computation (Maple), odd n = 5, 9, 11 numerically, n = 7 with Maple, and n = 3 by a
separate remark.

- We give a proof for every n ≥ 4, uniform in n, in which the decisive quartic coefficient is computed on the
  momentum level set; in particular it covers odd n ≥ 13.
- It rests on exact closed forms for that coefficient:
  - c₄ = (n−2)(m⁴+40m²+48m+16)/(192n(n−1)) with m = n−2, for even n;
  - c₄ = (n−1)(n−3)(n²−9n+38)/(128n), for odd n.
- It also includes a self-contained nonlinear-stability argument.

Paper: [conj13/paper/main.pdf](conj13/paper/main.pdf).

**2. Exact, certified thresholds for the tanh (Kimura) model on H².**

Boatto (2008) studied this model's ring stability numerically. Here the stability blocks are obtained in closed form.

- The blocks are explicit rational functions of ξ = tanh²(r/2), r the ring's hyperbolic radius.
- For N = 8–16 each threshold is the unique root in (0, 1) of an irreducible palindromic integer polynomial, and is
  certified: 25 digits are printed in [docs/results.md](docs/results.md), `tanh/certify_v2.out` gives 55 digits from
  Arb balls of radius < 3·10⁻¹⁵³ inside exact isolating intervals, and Lean checks intervals of width 10⁻²⁰.
- ξ*(8) is an algebraic unit of degree 10.
- For every N: the tangential stiffness b_m is positive and the critical mode is ⌊N/2⌋. Every N ≥ 8 has exactly one
  threshold (proved for even N and for odd N ≥ 17, certified for odd N = 9–15). For even N it is also the
  energy–momentum threshold; every odd N ≥ 9 has a non-empty Krein window, of width O(N⁻³).
- Large rings: cosh r* = κ∞N + O(1/N) with κ∞ = 0.13431… defined by (coth η − 1/η)² = 2 tanh(η)/η, η = 1/(2κ∞).
  Physically, a large ring stabilises once neighbouring vortices are more than ℓ* ≈ 0.844 curvature radii apart,
  measured along the ring. That is exactly the critical spacing of an infinite vortex row on a horocycle, the
  ring's N → ∞ limit.
- An independent implementation of Boatto's published equations agrees with these thresholds
  ([tanh/boatto_check/](tanh/boatto_check/README.md)).

Note: [tanh/note_tanh.pdf](tanh/note_tanh.pdf).

The sinh model, whose thresholds are Montaldi and Tokieda's, is written up in closed form (an expository note) in
[sinh/note.pdf](sinh/note.pdf).

## Checking the work yourself

| Check | How | Time |
|---|---|---|
| Lean formalisation of the exact algebra | `(cd lean && lake exe cache get)`, then `bash scripts/build_lean.sh` (every module one at a time, as CI does), then `(cd lean && lake env lean AxiomAudit.lean)`. See [docs/lean.md](docs/lean.md). | about 40 min, plus about 3 min for the Mathlib cache |
| Package test suite | `pip install -e ".[dev]" && pytest -q`. Expected: `399 passed, 6 skipped`; the six skipped tests are marked slow (run them with `pytest --runslow`; one of them also needs `lake` and a built `lean/`). | about 5 min |
| Re-run the computations behind the papers | `pip install -r requirements.txt`, then `python scripts/reproduce.py --tier fast`. See [REPRODUCE.md](REPRODUCE.md). All 98 recorded outputs reproduce byte for byte. | fast tier (59 outputs + the Lean-generator check) ~40 min; fast + medium (all 98) ~2.5 h |
| Irreducibility certificates (tanh) | `cd tanh/irreducibility && gp -q irreducibility.gp` | seconds |
| Cross-check against Boatto's published equations | [tanh/boatto_check/README.md](tanh/boatto_check/README.md) | under 1 min |

Prerequisites: Python ≥ 3.10 (the recorded outputs were reproduced with 3.12.3 and the versions pinned in
`requirements.txt`); [PARI/GP](https://pari.math.u-bordeaux.fr/) (`gp` on the PATH, version 2.15 used) for the
irreducibility certificates; [elan](https://github.com/leanprover/elan), which installs the Lean toolchain pinned in
`lean/lean-toolchain`, for the Lean check; and optionally `pdflatex` to rebuild the PDFs. For the reproduction, install
the pinned versions first and the package second:
`python3 -m venv .venv && .venv/bin/pip install -r requirements.txt && .venv/bin/pip install -e ".[dev]"`. (Installing
the package first also works, but pulls the newest numpy, which the pinned install then replaces.)

CI builds the Lean library on every push. The build fails on any `sorry`, `admit`, `native_decide` or extra axiom,
and on any main theorem that depends on an axiom other than `propext`, `Classical.choice` and `Quot.sound`.

## Repository layout

| Path | Contents |
|---|---|
| `conj13/` | Conjecture 13: computation engines, the checks in `review/`, `review3/` and `gaps/`, their outputs, and the paper in `paper/` |
| `tanh/` | The tanh model: block and certification scripts, checks in `review/`, `irreducibility/`, `boatto_check/`, the Lean certificate generator in `lean_gen/`, and the note |
| `sinh/` | The sinh-model note and its checks |
| `lean/` | Lean 4 + Mathlib library `Havelock` ([docs/lean.md](docs/lean.md)) |
| `docs/` | [results.md](docs/results.md): every claim with its status. [lean.md](docs/lean.md). Package docs: [quickstart.md](docs/quickstart.md), [models.md](docs/models.md), [certificates.md](docs/certificates.md). |
| `src/` | The `havelock` Python package (`src/havelock`) |
| `tests/` | The package's pytest suite |
| `independent_check/` | An independent from-scratch implementation (no shared code) that checks both results numerically, with its outputs; see its README |
| `explorations/` | Dated exploratory work (NUMERICAL unless marked): independent review checks, screened and shell-BEC rings, Rossby/strip roll-up, Saturn jets, acoustic ergoregions, and rows on curves of constant geodesic curvature (now in the tanh note as Prop. 22). Each dated folder has its own README |
| `literature.md` | The works cited, each with one sentence on what it contains |
| `scripts/`, `REPRODUCE.md` | The reproduction driver and the table of every recorded output |

## Using the `havelock` Python package

The `havelock` package (`src/havelock`, version 0.1.0) answers four questions about the regular N-vortex ring in a
given model: where it loses stability, which Fourier mode goes first, whether it is nonlinearly stable at the
bifurcation point, and what the motion looks like. Every answer carries a status (PROVED, CERTIFIED or NUMERICAL)
and a `basis` string that says why. Certified thresholds can be exported as machine-checkable certificates.

    python3 -m venv .venv && .venv/bin/pip install -e ".[dev]"
    .venv/bin/havelock threshold sinh 12
    .venv/bin/python -m pytest -q

What 0.1.0 covers:

- **MT family** (`mt`, `sphere`, `plane`, `sinh`): thresholds, critical mode and the quartic coefficient c4 are PROVED
  for n >= 4 (closed forms: even n >= 4, odd n >= 5). n = 2, 3 raise `OutOfValidityRange` on the sphere and the full
  family, and are stable at every radius in the plane and sinh views. Stability verdicts are PROVED for
  −1 < lam_r2 < 1 and CERTIFIED at the equator lam_r2 = 1.
- **tanh model** (`tanh`): block formulas PROVED; thresholds CERTIFIED at run time for N <= 64 (exact Sturm counts and
  Arb root isolation; "Lean-checked" for 8 <= N <= 16); stability verdicts CERTIFIED from exact signs; c4 NUMERICAL and
  for even N only. The certified tanh threshold is the spectral one. For every even N it is also the
  energy–momentum threshold (note_tanh Theorem 18(a)). For every odd N >= 9 the energy–momentum threshold lies
  strictly above it (Theorem 18(b)), in a Krein window where the ring is linearly stable and nonlinear stability is
  undecided; the result's basis says so, and also reports the run-time evidence when the certificate has it (exact
  sign D_crit < 0 just above the threshold interval).
- **Simulation** (`simulate`): NUMERICAL, all models.
- **Certificates**: JSON threshold certificates for every model, re-checked by an independent exact verifier
  (`havelock verify`); Lean 4 export for tanh N = 8..16 only (`havelock cert tanh 8 --lean CertN8.lean`); other N
  raise `CertificationFailed`.

See [docs/quickstart.md](docs/quickstart.md), [docs/models.md](docs/models.md) and
[docs/certificates.md](docs/certificates.md). The package does not import the research scripts in `conj13/`, `tanh/`
and `sinh/`; it ports their logic, and its tests compare against their recorded outputs.

## Citing and contact

Please cite the papers in `conj13/paper` and `tanh/` as drafts by Gordon Speagle (2026); `CITATION.cff` has the
details (GitHub's "Cite this repository" button reads it). Contact: gspeagle@gmail.com. Corrections and comments are
welcome, by email or as GitHub issues.

## Licence

Code (the `havelock` package, the research scripts and the Lean library) is under the [MIT licence](LICENSE). The papers,
notes and documentation are under [CC BY 4.0](LICENSE-docs.md).
