# Lean formalisation

The Lean 4 + Mathlib library in `lean/` checks the exact-algebra layer of both results: closed forms, thresholds,
certified root intervals and sign patterns. It does **not** formalise the dynamics. Each file's header comment lists
exactly what it proves and what it takes as input; the summary is in [results.md](results.md).

## Building and checking

[![lean](https://github.com/speagleg/havelock-vortex-rings/actions/workflows/lean.yml/badge.svg)](https://github.com/speagleg/havelock-vortex-rings/actions/workflows/lean.yml)

```bash
export PATH=$HOME/.elan/bin:$PATH
cd lean && lake exe cache get && lake build     # Lean and Mathlib v4.34.1
lake env lean AxiomAudit.lean                   # #print axioms for the main theorems of all three areas
```

CI (`.github/workflows/lean.yml`) runs on every push and pull request. It fails if any `sorry`, `admit`,
`native_decide` or `axiom` appears in the Lean code (comments excluded), if the build fails, or if a theorem listed in
`AxiomAudit.lean` depends on any axiom other than `propext`, `Classical.choice` and `Quot.sound`. CI builds the modules
one at a time.

The generated tanh certificates are split into modules of at most about 90 kB each, so no module needs much memory.
Each module peaks at **under 4.5 GB resident**, and only up to about 2 GB of that is anonymous memory; the rest is
mapped Mathlib `.olean` files. The slowest module takes about 4 minutes, and a full sequential check takes about
25 minutes. A parallel `lake build` checks several modules at once, so on a machine with less than about 16 GB of RAM,
build the modules one at a time, as CI does.

## Module layout (`lean/`)

- `Havelock.lean`: imports everything.
- `AxiomAudit.lean`: the `#print axioms` list that CI checks. It is not part of the library.
- `Havelock/SinhArithmetic.lean`, `Havelock/Conj13Algebra.lean`: the sinh and Conjecture 13 layers.
- `Havelock/TanhAlgebra.lean`: a thin module that imports `Havelock/Tanh/*` and prints the tanh axiom audit. All
  names are `Havelock.Tanh.*`, as before the split.
  - `Tanh/General.lean`: the definitions, the all-N lemmas, `peval` and the sign-change lemma.
  - `Tanh/StableN2to7.lean`: the certificates for N = 2–7 (`stable_N2` … `stable_N7`).
  - `Tanh/CertN8.lean` … `Tanh/CertN11.lean`, plus `Tanh/CertN12a,b`, `N13a,b`, `N14a,b`, `N15a,b` and
    `N16a,b,c`: the certificates for N ≥ 8 (`cert_N8` … `cert_N16`). The larger N are split into parts, and each part
    imports the one before it.
  - `Tanh/Extras.lean`: the factorisations, palindromy, the Krein window, ξ*(8) being a unit, and `tanh_relations`.
- `Havelock/TanhLargeN.lean`: imports `Havelock/Tanh/LargeN/*` (the tanh note's large-N section) in dependency order;
  its header lists each module and the note results it covers.

## What each file proves

- `Havelock/SinhArithmetic.lean`: the sinh-model arithmetic, with no `sorry`. It covers the integrality classification c(N) ∈ ℤ ⇔ N ∈ {8, 9, 11, 15, 23}, the unit property, ξ*(8) = 8 − 3√7 and ξ*(10) = 1/7.
- `Havelock/Conj13Algebra.lean`: the exact-algebra layer of the Conjecture 13 paper, with no `sorry`. Only the axioms `propext`, `Classical.choice` and `Quot.sound` are used.
  - **Thresholds b_n for all n ≥ 4.** It proves the closed form, that b_n is the unique root in (−1, 1], and the geometry split: sphere for n ≤ 6, plane at n = 7, H² for n ≥ 8.
  - **Spectral algebra.** It proves the Fejér sum identity Σ(1−cos 2πml/n)/(1−cos 2πl/n) = m(n−m) and finds the kernel mode and spectral gap, including 3/2 at n = 4 and 8/3 at n = 5.
  - **The c₄ closed forms, symbolically in n.** The even and odd closed forms are derived from the printed level-set coefficients, and c₄ > 0 at the true threshold for every n ≥ 4. Table 1 (n = 4–25) is exact over ℚ.
  - **The comparison with Montaldi–Tokieda (paper §7).** The P1 comparisons are exact: 9/14 < 6/7 at n = 7, and (199 − 27√7)/56 > 95/56 at n = 8. The P2 completing-the-square identity with Lemma 8's weights is also proved.
  - **Not formalised:** the Taylor-coefficient inputs (h₂, h₄, f₂, E_ij) are taken from the paper as definitions; they were derived by computer algebra and checked mod p. The energy–momentum and dynamics argument (§ nonlinear) is not formalised either.
- `Havelock/TanhAlgebra.lean` (with `Havelock/Tanh/*`): the exact-algebra layer of the tanh note, with no `sorry`. Only the standard axioms are used.
  - **For all N and m:** the W_k closed form (`HasSum`), the compact a_m, the reciprocity a_m(1/ξ) = a_m(ξ), the even-N threshold equation, and the linearisation's characteristic polynomial.
  - **For every N = 2–16:** the certification theorem. N ≤ 7 is stable for all r. For N ≥ 8, f_N has a unique root ξ* in (0, 1), located in a rational interval of width 10⁻²⁰, and the critical block changes sign there. The proofs use Bernstein positivity certificates.
  - **ξ*(8):** it is an algebraic integer and a unit.
  - **Odd N = 9, 11, 13, 15:** the Krein window is non-empty.
  - **Not formalised:** the derivation of the blocks (they are taken as definitions), the dynamics, and irreducibility (certified separately with PARI/GP, see `tanh/irreducibility/`; not in Lean). The asymptotics are in `TanhLargeN` (below).
  - The certificates are generated by `tanh/lean_gen/gen.py`. It writes one file.
    `python lean/scripts/split_tanh.py OUT lean` splits that file into `Havelock/Tanh/*` and the thin
    `TanhAlgebra.lean`, reproducing the committed files exactly.

- `Havelock/TanhLargeN.lean` (with `Havelock/Tanh/LargeN/*`): the large-N section of the tanh note (Theorems 11–18,
  Proposition 19, Corollary 20), with no `sorry`, stated about the same block definitions `aBlock`, `bBlock`, `gBlock`.
  - **For all N, m and radii:** the exact symbol form of the blocks (Thm 11), b_m > 0 (Lemma 12), the critical mode
    ⌊N/2⌋ (Lemma 13), spectral stability iff a_⌊N/2⌋ ≥ 0 (Cor 14), and the EM determinant identity with G ≥ 0 (Lemma 17).
  - **Even N:** the closed form and its monotonicity, one threshold for every even N ≥ 8, stability at every radius
    for N ≤ 6, and EM threshold = spectral threshold (Prop 15, Thm 18(a)).
  - **The limit:** Φ in Langevin form, strictly decreasing, its root τ∞ with 7.4453 < τ∞ < 7.4454; the thresholds
    `xiStar`, `xiEM` defined as infima; |cosh r* − κ∞N| ≤ C/N and the same for r_EM (Thms 16, 18(c)); the arc spacing
    ℓ_N = ℓ* + O(N⁻²); the odd-N Krein window non-empty for every odd N ≥ 9, with width O(N⁻³) (Thm 18(b)).
  - **The horocyclic row:** the four lattice sums (the Poisson-type sum via Fourier coefficients on the circle), the
    Bloch symbol identity (Prop 19), and the critical spacing for stability and definiteness (Cor 20).
  - **Not formalised:** as above, the blocks themselves are the starting point; for the row, the per-pair Hessian
    entries are taken as definitions (the note's proof states them; their derivation from the pair term is not in
    Lean); the dynamics (stability is stated through the block and Bloch-matrix criteria of note Prop 3); and the
    second-order constants C₁ and W, which the note itself labels as derived.

The axiom audit (`#print axioms`) shows only `propext`, `Classical.choice` and `Quot.sound`, or a subset of them.
