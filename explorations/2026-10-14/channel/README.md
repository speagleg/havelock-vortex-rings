# channel/ — the same vortex row in a flat channel: walls versus curvature

In strip coordinates z = u + i phi (0 < phi < pi), the hyperbolic plane's tanh kernel is exactly the Green's function of
a flat channel with impermeable walls at phi = 0 and pi. This is the conformal invariance of the Laplace Green's
function under z -> e^z (strip -> upper half-plane):

    W = -ln tanh(d/2) = ln |(e^z - e^{conj w}) / (e^z - e^w)|.

So a row of vortices at phi = alpha has the same pair interaction in both settings. Only the dynamics differ:

| | H^2 (hypercycle row) | flat channel (row parallel to the walls) |
|---|---|---|
| equations | u' = sin^2(phi) dH/dphi, phi' = -sin^2(phi) dH/du | u' = dH/dphi, phi' = -dH/du |
| Hamiltonian | sum_{i<j} W_ij | sum_{i<j} W_ij + (1/2) sum_i ln(2 sin phi_i) (Robin / image self-term) |

The channel's Robin term (1/2) ln sin(phi) is minus half the log conformal factor of H^2, sigma = -ln sin(phi). On H^2
the geometric self-energy of a vortex cancels the walls' image self-energy, and what remains is the mobility
sin^2(phi) = e^{-2 sigma}.

| Script | What it does |
|---|---|
| `channel_vs_h2.py` | Direct supercell linearisation of both systems on a grid of (alpha, spacing) -> `channel_vs_h2.out` |
| `channel_symbol.py` | Exact Bloch symbols of both, the identity between them, a 25,600-point scan -> `channel_symbol.out` |

**Results.**

1. *Exact identity* (checked to 3e-14). With V_ch the channel row's translation speed,
   B_H2(q) = sin^2(a) B_ch(q) + sin(2a) V_ch - cos(2a)/2  and  C_H2(q) = sin^2(a) C_ch(q).
   The cos 2a of the 45 degree rule appears here as the explicit difference between the two systems.
2. *The channel row is never stable.* It is unstable at all 25,600 points with alpha in (0, pi/2] and spacing in
   [0.02, 30] (exact symbol), and at every point of the direct supercell grid. The along-row stiffness satisfies
   min_q B_ch <= -1/2 everywhere. At large spacing B_ch(pi) -> -csc^2(a)/2: the wall-image self-term alone
   destabilises the zigzag mode. On H^2 the same quantity B_H2(pi) tends to 0 with the sign of cos 2a (the 45 degree rule).
3. *Reading.* Same interaction, opposite outcome: walls destabilise a vortex row, whereas the constant negative
   curvature of H^2 (which is the channel without the image term, but with mobility sin^2 phi) stabilises it once
   kappa > 1/sqrt 2 and the spacing exceeds l*(kappa).

Not proved: that min_q B_ch <= -1/2 for every alpha and spacing (only the grid). Literature: double vortex rows in
channels go back to Rosenhead (1929); the single-row channel computation here is a check made for this comparison,
not a claim of novelty.
