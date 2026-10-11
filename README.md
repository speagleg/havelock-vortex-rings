# Vortex-ring stability: Havelock thresholds beyond the plane

Gordon Speagle, Independent Researcher, Boise, ID, USA.

Research on the stability of the regular N-vortex ring (the Havelock problem). It covers the sphere, the plane and
the hyperbolic plane H², and two vortex models on H². Everything here is a **draft**. It has been reviewed by
independent AI readers and checked by computer algebra; **no human specialist has reviewed it yet**.

## Contents

| Directory | What | Status |
|---|---|---|
| `conj13/` | Proof of Montaldi–Tokieda's Conjecture 13 (arXiv:1210.5662, p. 23) for n ≥ 4: nonlinear stability of the ring on the reduced level set J⁻¹(μ)/SO(2) past the bifurcation. Repairs three errors in the published argument (their conclusions stand). Paper draft in `conj13/paper/`. | Draft; AI-reviewed. n = 3 rests on Montaldi–Tokieda's own statement. |
| `tanh/` | Thresholds for the tanh (Kimura) model H = −Σ ln tanh(d/2) on H². The Fourier blocks are rational in ξ = tanh²(r/2), the threshold polynomials are palindromic, and the table is certified for N ≤ 16. ξ*(8) is a degree-10 unit. Large-N asymptotics cosh r* ≈ κ∞N. Note in `tanh/note_tanh.pdf`. | Draft; AI-reviewed. Prior-work check incomplete (Kimura 1999 and Hwang–Kim 2013 not yet read). |
| `sinh/` | The sinh model (Montaldi–Tokieda (1.12)), written up in closed form as cosh² r* = q(N), with the classification of unit thresholds. | Expository. The thresholds are known (Montaldi–Tokieda 2013). |
| `lean/` | Lean 4 + Mathlib formalisation of the exact-algebra layer. The dynamics are not formalised; each file's header says exactly what is. | See below. |
| `literature.md` | Literature notes and the prior-work check. | |
| `directions-2026-10-10.md` | Open directions. | |

## Lean

```bash
export PATH=$HOME/.elan/bin:$PATH
cd lean && lake build            # Mathlib v4.34.1 (lake exe cache get on a fresh machine)
```

- `Havelock/SinhArithmetic.lean`: the sinh-model arithmetic, with no `sorry`. It covers the integrality classification c(N) ∈ ℤ ⇔ N ∈ {8, 9, 11, 15, 23}, the unit property, ξ*(8) = 8 − 3√7 and ξ*(10) = 1/7.
- The Conjecture 13 and tanh algebra layers are in progress.

## Reproducing the computations

The Python scripts are plain numpy/sympy/mpmath. Each `*.out` file is the recorded output of the script with the same name.
