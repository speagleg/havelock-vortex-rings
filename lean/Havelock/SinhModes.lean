/-
Havelock — sinh model: derivation of the mode decomposition and the threshold (sinh/note.tex,
Proposition prop:modes and §3), reusing the tanh machinery of `TanhBlocks` (the sinh pair terms are the
tanh ones without the image term, plus a one-body background term).

* `Sinh.Hessian`: `sinh_ring_critical(_iff)` (the ring is critical exactly for `sinhOmega`),
  `sinh_blocks_eigen`, `sinh_blocks_hermitian` (on mode m: a_m = (N−1)(1+c²)/2 − m(N−m)/2, b_m = m(N−m)/2,
  g_m = 0).
* `Sinh.Threshold`: `threshold` (spectral stability ⟺ a_{⌊N/2⌋} ≥ 0 ⟺ q(N) ≤ cosh² r), `threshold_strict`
  (the strict form, definiteness), `threshold_cosh`, `stable_small` (N ≤ 7), `threshold_xiStar` (N ≥ 8),
  `xiStar_add_inv` (link to `SinhArithmetic.c`).
-/
import Havelock.Sinh.Hessian
import Havelock.Sinh.Threshold
