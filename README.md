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

> **Status.** These are drafts. They have been checked by independent AI readers, by exact and interval computer
> algebra, and by a Lean build in CI. **No human specialist has reviewed them yet.** Each claim is labelled PROVED,
> CERTIFIED or NUMERICAL in [docs/results.md](docs/results.md). The work was developed with AI assistance (Claude);
> each paper's acknowledgement says how.

## The results

**1. Montaldi–Tokieda's Conjecture 13, completed.**

Montaldi and Tokieda (2013) conjectured that the ring stays stable *at* the curvature where it bifurcates. They proved
it for even n and for n = 3, 5, 7, 9, 11, using Maple and numerics.

- We give the first proof for odd n ≥ 13. The argument is uniform in n, so it also re-proves 4 ≤ n ≤ 11.
- It rests on exact closed forms for the decisive quartic coefficient:
  - c₄ = (n−2)(m⁴+40m²+48m+16)/(192n(n−1)) with m = n−2, for even n;
  - c₄ = (n−1)(n−3)(n²−9n+38)/(128n), for odd n.
- It also includes a self-contained nonlinear-stability argument.

Paper: [conj13/paper/main.pdf](conj13/paper/main.pdf).

**2. Exact, certified thresholds for the tanh (Kimura) model on H².**

Boatto (2008) computed this model's ring stability numerically. Montaldi and Tokieda later called it "no longer
tractable analytically".

- The stability blocks are explicit rational functions of ξ = tanh²(r/2).
- Each threshold is a root of an irreducible palindromic integer polynomial. The thresholds are certified to 25 digits
  for N ≤ 16.
- ξ*(8) is an algebraic unit of degree 10.
- For every N, not just N ≤ 16: the tangential stiffness b_m is positive and the critical mode is ⌊N/2⌋; every even
  N ≥ 8 has exactly one threshold, which is also its energy–momentum threshold; every odd N ≥ 9 has a non-empty Krein
  window, of width O(N⁻³).
- Large rings: cosh r* = κ∞N + O(1/N) with κ∞ = 0.13431… defined by (coth η − 1/η)² = 2 tanh(η)/η, η = 1/(2κ∞).
  Physically, a large ring stabilises once neighbouring vortices are more than ℓ* ≈ 0.844 curvature radii apart,
  measured along the ring. That is exactly the critical spacing of an infinite vortex row on a horocycle, the
  ring's N → ∞ limit.
- Implemented as written, Boatto's equations reproduce our thresholds exactly.

Note: [tanh/note_tanh.pdf](tanh/note_tanh.pdf).

The sinh model, whose thresholds Montaldi and Tokieda already know, is written up in closed form in
[sinh/note.pdf](sinh/note.pdf).

## Checking the work yourself

| Check | How | Time |
|---|---|---|
| Lean formalisation of the exact algebra | `cd lean && lake exe cache get && lake build`, then `lake env lean AxiomAudit.lean`. With less than 16 GB of RAM, build module by module as in [docs/lean.md](docs/lean.md). | about 25 min |
| Package test suite | `pip install -e ".[dev]" && pytest -q` | ~3 min |
| Re-run the computations behind the papers | `pip install -r requirements.txt`, then `python scripts/reproduce.py --tier fast`. See [REPRODUCE.md](REPRODUCE.md). All 70 recorded outputs reproduce byte for byte. | fast tier ~25 min; medium ~2 h |
| Irreducibility certificates (tanh) | `cd tanh/irreducibility && gp -q irreducibility.gp` | seconds |
| Reimplementation of Boatto's equations | [tanh/boatto_check/README.md](tanh/boatto_check/README.md) | — |

CI builds the Lean library on every push. The build fails on any `sorry`, `admit`, `native_decide` or extra axiom,
and on any main theorem that depends on an axiom other than `propext`, `Classical.choice` and `Quot.sound`.

## Repository layout

| Path | Contents |
|---|---|
| `conj13/` | Conjecture 13: derivations (`derivation.txt`, `gaps/PROOF.txt`), computation engines, review checks, and the paper in `paper/` |
| `tanh/` | The tanh model: derivation, block and certification scripts, `irreducibility/`, `boatto_check/`, the Lean certificate generator in `lean_gen/`, and the note |
| `sinh/` | The sinh-model note and its checks |
| `lean/` | Lean 4 + Mathlib library `Havelock` ([docs/lean.md](docs/lean.md)) |
| `docs/` | [results.md](docs/results.md): every claim with its status. [lean.md](docs/lean.md). Package docs: [quickstart.md](docs/quickstart.md), [models.md](docs/models.md), [certificates.md](docs/certificates.md). |
| `src/` | The `havelock` Python package (`src/havelock`) |
| `tests/` | The package's pytest suite |
| `independent_check/` | A from-scratch second-reader check of both results (no shared code), with its outputs; see its README |
| `literature.md` | What we read, and how much of it (full text or abstract only) |
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
  family, and are stable at every radius in the plane and sinh views.
- **tanh model** (`tanh`): block formulas PROVED; thresholds CERTIFIED at run time for N <= 64 (exact Sturm counts and
  Arb root isolation; "Lean-checked" for 8 <= N <= 16); stability verdicts CERTIFIED from exact signs; c4 NUMERICAL and
  for even N only. The certified tanh threshold is the spectral one. For even N <= 16 it is also the
  energy–momentum threshold (note_tanh); for N > 16 only the spectral threshold is certified. For odd N the
  energy–momentum threshold can lie strictly above it, in a Krein window where the ring is linearly stable and
  nonlinear stability is undecided (the note proves the window non-empty for N = 9, 11, 13, 15); the result's basis
  states this gap only when the certificate proves it (exact sign D_crit < 0 just above the threshold interval).
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
