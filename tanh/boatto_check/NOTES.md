# Boatto (2008): reading notes

> The paper PDFs (`dl/`), their text layers (`boatto2008.txt`, `mparc06-311.txt`) and the figure renders (`fig/`) named
> below are not committed (third-party copyright; regenerable from the publisher's PDF). See `RESULT.md` for the
> conclusions and for which scripts need these files.

Source: S. Boatto, "Curvature perturbations and stability of a ring of vortices", DCDS-B 10 (2008), no. 2&3, 349–375.
- Publisher PDF: `dl/boatto2008.pdf`, from `https://www.aimsciences.org/data/article/export-pdf?id=c627c07f-95aa-4501-9c88-aaad38ec9bdb`. Text is in `boatto2008.txt`.
- Preprint: mp_arc 06-311 (Nov 2006), `dl/mparc06-311.pdf`. Its text layer is font-garbled. The preprint's Fig. 5 (p. 17, rendered in `fig/preprint-17.png`) is the same figure as the published one.

Page numbers below are the journal's.

## Coordinates and metric
- p. 353, eq. (10): metric "dz dz̄ / (1 + K|z|²)²". As printed this has curvature 4K, because the factor 4 is missing. This does not matter: everything she computes uses geodesic distance r and √|K| r.
- p. 353, eq. (11): Darboux coordinates "q_α = φ_α, p_α^K = (Γ_α/K)(cosh(√|K| r_α) − 1), K < 0". For the sphere, "φ_α and r_α are respectively the longitude and the geodesic distance".
- p. 361, Remark 2: "r_α being the geodesic distance of the αth vortex from the center of the ring".

## How the ring radius is held fixed
- p. 362, Remark 2: "Since we are interested in studying the stability of the relative equilibrium of a ring, q_α = 2π(α−1)/N, p_α = p_o = (1/K)(cosh(√|K| r_o) − 1) (K<0) …, and how the stability properties change with the curvature of the surface K, there is no loss of generality in fixing r_o = 1."
- So the **geodesic radius is fixed, r_o = 1**, and K varies. The dimensionless radius in curvature −1 units is r = √|K|·r_o = √|K|.
- Caveat: in eqs. (39)–(40) (p. 361) the symbol r_o is reused for the embedded radius, "z_o = cos(√K), r_o = √(1 − z_o²)", i.e. sin √K. The geodesic radius is still 1 there.

## Green's function and Hamiltonian
- p. 354, eq. (14), hyperbolic surface: "Ψ_K = −(1/2π) ln[tanh((√|K|/2) r)]". No Robin or background term is added on the hyperbolic surface; a background term appears only for the whole sphere, eq. (15).
- p. 354, eq. (17): "H_{K<0} = c Σ_{α<β} Γ_αΓ_β ln((d_αβ − 1)/(d_αβ + 1)) (hyperbolic surface)". Eq. (16) is the punctured sphere, ln((1−d)/(1+d)). Eq. (18) is the whole sphere S², ln(1 − d).
- p. 355: "with c = −1/(4π)".
- p. 355, eq. (19), and p. 355, eq. (22) with Γ = 1: d_αβ = (1+Kp_α)(1+Kp_β) + K√(p_α p_β (2+Kp_α)(2+Kp_β)) cos(q_α − q_β).
  For K<0 this is cosh(√|K| D_αβ) by the hyperbolic law of cosines. So (d−1)/(d+1) = tanh²(√|K| D/2), and **(17) is exactly our tanh model**: H = (1/2π)·[−Σ ln tanh(D/2)].
- p. 358, eq. (28): "H = H1 + H2, H1 = c Σ log(1 − d_αβ), H2 = −c Σ log(1 + d_αβ)".

## Conserved quantity and rotating frame
- p. 357: M_z = (1/√|K|) Σ (1 + K p_α), which is LINEAR in p. The x- and y-components contain √(Kp(Kp+2)).
- p. 358, eq. (29): "In the corotating reference system, the Hamiltonian takes the form H̃ = H + νM".
- p. 356, eq. (27) gives ν_{K<0}.
- Because M_z is linear in the canonical p, the Hessian of H̃ in (q, p) equals the Hessian of H in (q, p). The ring is a critical point of H̃. We checked this numerically: ∂H/∂q = 0 and ∂H/∂p_α is the same for all α, to 1e-41 (`boatto_hessian.py`).

## Which matrix she monitors
- p. 358, (31): "S̃ is the Hessian of H̃ evaluated at the equilibrium X*". This is the full 2N×2N Hessian in the canonical (q, p), with no reduction.
- p. 360, ii): "It is preferable to work in the full system (as opposed to the reduced system) … S̃ … has a zero eigenvalue, due to the invariance under rotation … if all the remaining eigenvalues have the same sign we shall consider the Hessian as definite … stability is guaranteed when the minimum eigenvalue is positive."
- p. 359, b): "S̃ = [[Q, G],[Gᵀ, P]] where … Q and P are symmetric circulant matrices, while G is an asymmetric circulant matrix".
- p. 360–361, eqs. (38)–(40): explicit Q, P, G for the punctured sphere only. Two flaws there: "(D2)_1j" in t2 is never defined, and (38) is written with "H̃1" in every entry.
  p. 360: "The calculations for the hyperbolic case are done is a similar way." The hyperbolic blocks are **never written out**.
- p. 363, eqs. (44)–(45): λ_j = ½[(λ_Qj + λ_Pj) ∓ √((λ_Qj − λ_Pj)² + 4|λ²_Gj|)].
- p. 363: "we (numerically) monitor the maximum and the minimum eigenvalues of S̃ … stability is assured when the minimum eigenvalue of S is positive … for N ≥ 7, we can only assure stability for sufficiently negative K."
- Fig. 5a caption (p. 363): "minimum eigenvalue of the Hessian S versus the curvature K".
- Her criterion counts every non-rotational eigenvalue, including the breathing mode m = 0. That mode is always positive here (4πP_0 ≈ 9–16). So her criterion is the same as our energy–momentum (Hessian-definiteness) test.
- Small typographical point: the published remark reads "λ_G1 = 0 & λ_P1 = 0 → λ_1 = 0"; the rotation zero is in Q, and the preprint, p. 17, prints λ_Q1.

## Map K ↔ ξ (task 2)
r_o = 1 geodesic, so the curvature −1 radius is r = √(−K) and

    ξ = tanh²(√(−K)/2),        K = −(2 artanh √ξ)² = −r².

This is the map used in the comparison; it was re-derived independently here.

## Figure digitization
- Fig. 5a: rendered at 600 dpi (`fig/hi-15.png`, crop `fig/fig5a_crop.png`) and digitized by curve colour (`digitize_fig5a.py`). The axis frame is K ∈ [−0.6, 0.6], values in [−5, 5].
- Measured zero crossings: N=6 at +0.2024, N=7 at −0.0011, N=8 at −0.2465, N=9 at −0.3964. N=5 is near +0.449; label overlap makes it less reliable.
- Fig. 5b: the same procedure is in `fig5b_digitize.py`.
- Accuracy is about ±0.001–0.002 in K, judged from the N=7 crossing, which must be at 0 (planar heptagon).
