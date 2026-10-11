# Vortex-ring stability beyond the plane

[![lean](https://github.com/speagleg/havelock-vortex-rings/actions/workflows/lean.yml/badge.svg)](https://github.com/speagleg/havelock-vortex-rings/actions/workflows/lean.yml)

Gordon Speagle, Independent Researcher, Boise, ID, USA.

The regular ring of N identical point vortices is a classical object: Thomson posed the question in 1883, Havelock
answered it in 1931: in the plane the ring is linearly stable exactly when N ≤ 7 (N = 7 only marginally; its
nonlinear stability was settled later). This repository studies the same question
on curved surfaces: the sphere, the plane and the hyperbolic plane. It contains two research results, the proofs and
computations behind them, and a Lean 4 check of their exact algebra.

> **Status.** These are drafts. They have been checked by independent AI readers, by exact and interval computer
> algebra, and by a Lean build in CI. **No human specialist has reviewed them yet.** Each claim is labelled PROVED,
> CERTIFIED or NUMERICAL in [docs/results.md](docs/results.md).

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
- Boatto's equations reproduce our values exactly, but her published hyperbolic-plane curves do not. The note explains
  why.

Note: [tanh/note_tanh.pdf](tanh/note_tanh.pdf).

The sinh model, whose thresholds Montaldi and Tokieda already know, is written up in closed form in
[sinh/note.pdf](sinh/note.pdf).

## Checking the work yourself

| Check | How | Time |
|---|---|---|
| Lean formalisation of the exact algebra | `cd lean && lake exe cache get && lake build`, then `lake env lean AxiomAudit.lean`. See [docs/lean.md](docs/lean.md). | about 25 min |
| Re-run the computations behind the papers | `pip install -r requirements.txt`, then `python scripts/reproduce.py --tier fast`. See [REPRODUCE.md](REPRODUCE.md). All 72 recorded outputs reproduce byte for byte. | fast tier ~25 min; medium ~2 h |
| Irreducibility certificates (tanh) | `cd tanh/irreducibility && gp -q irreducibility.gp` | seconds |
| Boatto comparison | [tanh/boatto_check/RESULT.md](tanh/boatto_check/RESULT.md) | — |

CI builds the Lean library on every push. The build fails on any `sorry`, `admit`, `native_decide` or extra axiom,
and on any main theorem that depends on an axiom other than `propext`, `Classical.choice` and `Quot.sound`.

## Repository layout

| Path | Contents |
|---|---|
| `conj13/` | Conjecture 13: derivations (`derivation.txt`, `gaps/PROOF.txt`), computation engines, review checks, and the paper in `paper/` |
| `tanh/` | The tanh model: derivation, block and certification scripts, `irreducibility/`, `boatto_check/`, the Lean certificate generator in `lean_gen/`, and the note |
| `sinh/` | The sinh-model note and its checks |
| `lean/` | Lean 4 + Mathlib library `Havelock` ([docs/lean.md](docs/lean.md)) |
| `docs/` | [results.md](docs/results.md): every claim with its status. [lean.md](docs/lean.md). `design/`: the design spec for the planned Python package |
| `literature.md` | What we read, and how much of it (full text or abstract only) |
| `directions-2026-10-10.md` | Working notes on further directions, dated. Some items have since been resolved; for example, Conjecture 13 is item 7 |

## Coming next: the `havelock` Python package

A library and CLI is planned; its design is in [docs/design/](docs/design/). It will answer, for any N and any of the
four models, where the ring loses stability, which mode goes first, and whether the ring is nonlinearly stable at that
point. Every answer will carry its PROVED, CERTIFIED or NUMERICAL label and can be exported as a machine-checkable
certificate. It is not built yet.

## Citing and contact

Please cite the papers in `conj13/paper` and `tanh/` as drafts by Gordon Speagle (2026). Contact: gspeagle@gmail.com.
The work was developed with AI assistance (Claude), as each paper's acknowledgement states.
