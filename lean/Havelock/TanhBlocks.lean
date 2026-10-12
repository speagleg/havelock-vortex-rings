/-
Havelock — tanh model: derivation of the Fourier blocks (tanh/note_tanh.tex, Lemma 1 and Theorem 2) from the
Hamiltonian. `ringF N ξ Ω s θ` is H − ΩJ for N unit vortices on the Poincaré disk at
z_j = √ξ w^j exp(s_j + iθ_j); derivatives are taken along lines t ↦ F(ts, tθ), and the Hessian entries by
polarisation (see the module docstrings).

* `Tanh.BlocksSums`, `Tanh.BlocksCalc`: roots-of-unity filter, the complex Fejér identity `fejerC`, and the
  calculus of the pair, image and one-body terms.
* `Tanh.Blocks`: `ring_critical` / `ring_critical_iff` (note Lemma 1: the ring is critical exactly for the
  rotation rate Ω with Ωκ₀ = `omegaKappa`), `grad_formula`.
* `Tanh.BlocksHess`: `hessQ_eq`, `hessSS_eq`, `hessTT_eq`, `hessST_eq`, `hessTS_eq` (the Hessian as explicit
  circulant matrices).
* `Tanh.BlocksFourier`: `blocks_eigen` (note Thm 2: on the Fourier vector (w^{mb})_b the Hessian acts by the
  blocks `aBlock`, `bBlock`, `gBlock` of `General.lean`), `blocks_eigen_zero`, `blocks_hermitian` (no coupling
  between different modes).
-/
import Havelock.Tanh.BlocksSums
import Havelock.Tanh.BlocksCalc
import Havelock.Tanh.Blocks
import Havelock.Tanh.BlocksHess
import Havelock.Tanh.BlocksFourier
