# Threshold certificates

`hv.certificate(threshold)` (CLI: `havelock cert MODEL N --json PATH`) produces a JSON object. `havelock verify PATH`
re-checks every field with exact integer and rational arithmetic only (python-flint factorisation and Sturm
sequences, `fractions.Fraction` evaluation). The verifier rebuilds the block numerators from the closed forms itself;
it does not use the code that produced the certificate.

## Limit on N

Certificates are issued only for N <= 100000 (`hv.limits()["certificate_max_N"]`), so that the independent verifier
can always re-check them. `hv.certificate` raises `CertificationFailed` for a larger N. `verify` reports a malformed
file (for example invalid JSON) as an invalid certificate rather than crashing; on the command line a missing,
unreadable or non-UTF-8 file exits with code 2 and a message.

## Fields (schema `havelock-threshold-certificate/1`)

| field | meaning | how `verify` checks it |
|---|---|---|
| `schema` | format version | equality |
| `model`, `N` | model name and number of vortices | known model, N >= 2 |
| `mode` | critical Fourier mode | 1 <= mode <= N//2 |
| `polynomial` | integer numerator A_mode of a_mode (constant term first) | equals the rebuilt numerator |
| `factor_index` | which irreducible factor carries the threshold (factors primitive, leading coefficient > 0, sorted by degree then coefficients) | the factor changes sign on `interval` |
| `interval` | rational endpoints `"p/q"` | inside the domain; `polynomial` has opposite signs at the ends |
| `endpoint_signs` | signs of `polynomial` at the endpoints | recomputed |
| `sturm` | number of distinct real roots of `polynomial` in the interval | must be 1; recomputed (Sturm) |
| `interval_clear` | for every mode m, the number of distinct real roots inside `interval` of B_m, and of A_m for m other than the critical mode | recomputed (Sturm); every count must be 0, so no other block vanishes in the threshold interval |
| `points` | rational points just below and just above the interval | outside the interval, inside the domain |
| `mode_signs` | for every mode m: signs of A_m, B_m (and D_m = numerator of a_m b_m − g_m², when some g_m ≠ 0) at both points | recomputed; all positive on the stable side; a_mode b_mode < 0 on the other |
| `stable_side` | the open interval between the threshold interval and the domain end on the stable side, with the number of roots of every A_m and B_m in it | recomputed (Sturm); all zero. If the domain endpoint on the stable side is closed, the blocks must also be positive there |
| `irreducibility_witness` | optional prime p with the factor irreducible mod p | p prime, p ∤ lead, factorisation degrees mod p equal [deg] |
| `provenance` | package version, inputs, timing | not checked |

Together these prove: the ring is stable on the whole stable side (every a_m, b_m > 0 there), the critical mode is
unstable just across the interval, and the threshold lies in the interval as the unique root there of the
irreducible factor.

`points`, `stable_side` and `interval_clear` extend the field list of the design spec (Section 7); without them the
verifier could check the signs only at two points, not on the whole stable side.
