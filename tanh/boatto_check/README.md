# Reimplementation of Boatto (2008)

S. Boatto, "Curvature perturbations and stability of a ring of vortices", DCDS-B 10 (2008), no. 2&3, 349–375,
doi:10.3934/dcdsb.2008.10.349. Boatto formulated this stability problem for the image ("tanh") Green's function on
the hyperbolic surface, gave the circulant block structure of the Hessian with a coupling block, used the
Hessian-definiteness criterion, and found numerically that rings with N ≥ 7 are stabilised by sufficiently negative
curvature. Credit for the question, the block structure and the qualitative picture is hers.

The scripts here implement her equations as written, independently of `../blocks.py`, and compare them with the
thresholds of the tanh note.

- `boatto_hessian.py`: her canonical variables (eqs. (11), (21)), the pair quantity d (eq. (22)), the Hamiltonian
  (eqs. (16)–(18)) and the full 2N×2N Hessian of H̃ = H + νM at the ring (eqs. (29), (31)), at 40 digits. Her ring has
  fixed geodesic radius r_o = 1 (Remark 2) and the curvature K varies, so ξ = tanh²(√(−K)/2), K = −(2 artanh √ξ)².
- `map_and_verify.py` → `map_and_verify.out`: the definiteness threshold K* from her equations for N = 8, 9, 10,
  mapped to ξ, compared with `../blocks.py`; and a mode-by-mode comparison of her Fourier blocks (Q_m, P_m, G_m)
  with ours (b_m, a_m, g_m) at K = −0.25.
- `modes.py` → `modes_boatto.out`: per-mode Q, P, G and block eigenvalues for N = 6, 8, 9.

Result: her equations reproduce our thresholds exactly, K*(8) = −0.28218454857650 (ξ*(8) to 10⁻⁴¹) and
K*(10) = −0.74890006954582 (ξ*(10) to 10⁻⁴²). For odd N her criterion (definiteness) is our energy–momentum
threshold: K_EM(9) = −0.48869837, i.e. ξ_EM(9), including the certified gap ξ_EM − ξ* = 3.22·10⁻⁸. Her blocks equal
ours up to m-independent factors (Q_m/b_m = 1/2π, P_m/a_m = 0.134906 for every m at K = −0.25).

The paper itself is not included (third-party copyright); fetch it from the publisher.
