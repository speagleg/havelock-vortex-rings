# Item 01 — literature check (Rule L), 2026-10-07

Scope: core claims C1–C4 (N-ring thresholds on H²). For each source: whether I obtained and read it,
and what it says about N-ring stability thresholds on H². Local copies were kept in the session scratchpad
and are not committed.

## Verdict (short)

- **C1 (mode decomposition, C₁(ξ) = (N−1)(1+ξ²)/(1−ξ)²) is KNOWN.** It is Montaldi & Tokieda (2013),
  §2.2 eq. (2.6), together with the stability criterion (2.11) / Theorem 2, for the
  "background-vorticity" Green's function (1.12) continued to λ < 0. On the hyperbolic plane that model is
  exactly the sinh model of this item (pair term ln[|z−w|²/((1−|z|²)(1−|w|²))] = 2 ln sinh(d/2) + const).
- **C2 (marginal threshold) is KNOWN.** Montaldi–Tokieda (2.11) reads
  (1+λ²r₀⁴)/(1+λr₀²)² > ⌊n²/4⌋/(2(n−1)). With ξ = −λr₀² this is (N−1)(1+ξ²)/(1−ξ)² > ⌊N²/4⌋/2 = T,
  i.e. C1/C2. Their table of bifurcation values λr₀² for n = 8…13 is −0.0627, −0.101, −0.143, −0.172,
  −0.202, −0.225, i.e. −ξ*(N) (compare work/verify_core.out: 0.06275, 0.10102, 0.14286, 0.17157, 0.20241,
  0.22515). Their worked example "n = 15 stable for r₀² > (2−√3)/|λ|" is ξ*(15) = 2−√3.
  The **cosh² r* = q(N) form** is a one-line rewriting of (2.11) through (1+ξ²)/(1−ξ)² = (1+cosh²r)/2.
  It does not appear in that form in what I read. It is not a new result, and should be presented as a
  reformulation.
- **C3 (unit classification N ∈ {8,9,11,15,23}) was not found in any source I read.** It is an elementary
  arithmetic remark about the known thresholds: POSSIBLY-NEW but trivial. It should not be sold as a result
  about vortex dynamics.
- **Model-attribution problem (new finding).** Two independent Montaldi-group sources say that Kimura (1999)
  used the **tanh** (image) Green's function on H², not the sinh one:
  - Montaldi & Nava-Gaxiola (2014), eq. (12): H = −(1/4π) Σ ΓᵢΓⱼ ln[(⟨Xᵢ,Xⱼ⟩+1)/(⟨Xᵢ,Xⱼ⟩−1)] "as constructed
    by Kimura [10]". With ⟨X,Y⟩ = x₁x₂+y₁y₂−z₁z₂ = −cosh d (Nava-Gaxiola thesis eq. 2.5), the ratio is
    tanh²(d/2). They also note that the pair energy tends to 0 as d → ∞, which holds for tanh and not for sinh.
  - Montaldi & Tokieda (2013) §4.2: for G = log |z−w|²/|1+λzw̄|², "On the hyperbolic plane λ < 0 this reduces
    to the Hamiltonian adopted by Kimura [10]". They add: "The calculations based on this option … get so
    cumbersome that the stability problem seems no longer tractable analytically".
  So the source corpus's name "Kimura's sinh model" is very probably a misattribution. **I could not obtain
  Kimura 1999 itself** (royalsocietypublishing.org returns 403/Cloudflare; DOI 10.1098/rspa.1999.0311 confirmed
  via Crossref), so under Rule C this stays a strong secondary-source finding until a human checks the
  original. In this item's note the model is therefore called the **sinh model**
  (= Montaldi–Tokieda's family (1.12) at λ < 0), never "Kimura's model".
  Consequence: the threshold problem for the model that the literature attributes to Kimura is the **tanh**
  problem (C5/B3). Its closed form is open and Montaldi–Tokieda explicitly left it untreated.

## Per-source record

| Source | Obtained / read? | What it says about N-ring thresholds on H² |
|---|---|---|
| **Montaldi, J. & Tokieda, T.** "Deformation of geometry and bifurcations of vortex rings", in *Recent Trends in Dynamical Systems*, Springer Proc. Math. Stat. 35 (2013) 335–370; arXiv:1210.5662 | **Yes, full text (arXiv PDF)** | Smooth family M_λ (curvature 4λ), symplectic form 2dx∧dy/σ², σ = 1+λ\|z\|². Green's function (1.12) log[\|z−w\|²/(σ_zσ_w)] (sinh model at λ<0). §2: ring angular velocity (2.1); Hessian in (r,θ) block-diagonal with no r–θ coupling (2.2); Fourier eigenvalues (2.6) ε_r^(ℓ) = κ²/(4πr₀²)[2(n−1)(1+λ²r₀⁴)/σ² − ℓ(n−ℓ)], ε_θ^(ℓ) = κ²ℓ(n−ℓ)/4π; symplectic slice (2.9); **Theorem 2**: Lyapunov-stable (modulo G_µ) if (2.11) holds; Remark 4: linearly unstable if (2.11) is reversed. "Hyperbolic planes: when n ≤ 7 the ring is always stable. When n > 7, it is stable (for a given λ) for a sufficiently large." Table of bifurcation values for n = 4…13. §3: bifurcations at the threshold and nonlinear stability of the heptagon (Cor. 7); Conjecture 13 (n-gon stable at the bifurcation point). §4.1: the criterion for G = log\|z−w\|² (1 + (5/6)λ²r₀⁴ in place of 1+λ²r₀⁴). §4.2: the tanh/Kimura Green's function is "no longer tractable analytically". Uses identity (2.4) Σ cos(2πℓj/n)/(1−cos(2πj/n)) = (n²−1)/6 − ℓ(n−ℓ), citing Hansen's table. **⇒ C1, C2 KNOWN.** |
| **Montaldi, J. & Nava-Gaxiola, C.** "Point vortices on the hyperbolic plane", J. Math. Phys. 55, 102702 (2014); arXiv:1403.2138 | **Yes, full text** | Only N = 2, 3 are analysed (relative equilibria, SL(2,ℝ) momentum types, stability). The intro summarises Boatto and Montaldi–Tokieda: "for a given radius of ring, the ring becomes Lyapunov stable as the curvature decreases and for N > 7 stability only occurs for negative curvature". Rings on H² are "mentioned briefly in [3]" (Aref et al. 2003); "a more in-depth study is provided by Hwang and Kim [8]", which presents conditions for relative equilibria of rings. Eq. (12): Kimura's Hamiltonian in tanh form (see above). No explicit H² ring threshold formula. |
| **Nava-Gaxiola, C.** PhD thesis "Point vortices on the hyperboloid", Manchester 2013 (eprints.maths.manchester.ac.uk/2022) | Yes (intro, ch. 2, ch. 4 read; searched for "ring") | Same Hamiltonian (4.18) attributed to Kimura (tanh form). States that Hwang–Kim 2009 used a *different* Hamiltonian from Hwang–Kim 2013, and that the 2013 one coincides with Kimura's (tanh). No ring thresholds. |
| **Kimura, Y.** "Vortex motion on surfaces with constant curvature", Proc. R. Soc. Lond. A 455 (1999) 245–259, doi:10.1098/rspa.1999.0311 | **No** (403 / Cloudflare; no open copy found) | Per the two Montaldi-group sources: unified S²/H² formulation, tanh-type Hamiltonian on H², dipole motion (geodesic conjecture). Nothing on ring thresholds is reported by the secondary sources. **Unverified (Rule C): do not quote anything from Kimura directly.** |
| **Hwang, S. & Kim, S.-C.** "Point vortices on hyperbolic sphere", J. Geom. Phys. 59 (2009) 475–488, doi:10.1016/j.geomphys.2009.01.003 | **No** (ScienceDirect 403) | Per Montaldi–Nava-Gaxiola and the thesis: fixed and relative equilibria of three vortices, collapse. Different Hamiltonian from the 2013 paper. |
| **Hwang, S. & Kim, S.-C.** "Relative equilibria of point vortices on the hyperbolic sphere", J. Math. Phys. 54, 063101 (2013), doi:10.1063/1.4811454 | **Abstract only** (Crossref; AIP 403) | Abstract: "…In particular, the configuration of rings of vortices is investigated in detail." Per secondary sources: conditions for ring relative equilibria, any 2-vortex configuration is a relative equilibrium, Kimura's (tanh) Hamiltonian. **Unknown whether it contains a stability threshold.** If it does, it is for the tanh model and would bear on C5, not C1. Needs a human with journal access. |
| **Boatto, S.** "Curvature perturbations and stability of a ring of vortices", DCDS-B 10 (2008) 349–375 (preprint mp_arc 06-311) | Abstract read; preprint PDF downloaded but its text is unextractable (Type-3 fonts); TeX source link returns HTTP 500 | Abstract: stability of N-gons as a function of N and curvature K, Dirichlet criterion; "Stability ranges are the K-intervals for which the Hessian … is positive or negative definite." Montaldi–Tokieda and Montaldi–Nava-Gaxiola credit it, jointly with Montaldi–Tokieda, for "N > 7 stable only for negative curvature". Its explicit formula and Green's function are **unverified**. It may well contain C1/C2 earlier than 2013. Treat as a probable prior source. |
| **Polvani, L. M. & Dritschel, D. G.** J. Fluid Mech. 255 (1993) 35–64 | Not read (paywalled) | Cited by Montaldi–Tokieda as the first linear-stability study of rings on **S²**. Nothing on H². |
| **Boatto, S. & Cabral, H. E.** SIAM J. Appl. Math. 64 (2003) 216–230 | Not read | Cited by Montaldi–Tokieda as the nonlinear stability of a latitudinal ring on **S²** ("for spheres is due to [2]"). Nothing on H². |
| **Laurent-Polz, Montaldi, Roberts** J. Geom. Mech. 3 (2011) 439–486 (arXiv math/0402430) | Located, not read in detail | S² symmetric relative equilibria. Not H². |
| **Boatto, S. & Koiller, J.** "Vortices on closed surfaces", arXiv:0802.4313 | Yes (searched) | General Green's-function/Robin formulation on surfaces. No H² ring thresholds. |
| Searches | — | Web: "point vortex ring hyperbolic plane stability threshold N>7 relative equilibria" (extended); "Montaldi Tokieda Deformation of geometry…"; "Boatto Curvature perturbations…"; "Hwang Kim Relative equilibria … ring stability"; "polygonal/regular polygon point vortices constant negative curvature Lyapunov stability ring radius threshold explicit" (extended); "vortex ring hyperbolic plane stability threshold 8−3√7 / algebraic unit / cosh²". Crossref for DOIs. None of the hits states the cosh² form or the unit classification. |

## Consequences for the item

1. C1 and C2: relabel **KNOWN** (Montaldi–Tokieda 2013, eq. (2.6), (2.11), Thm 2 and table; probably also
   Boatto 2008, unverified). The note is **expository** for these. Its value is a short self-contained
   derivation in log-polar coordinates (Hessian diagonal per Fourier mode), a high-precision independent
   dynamics check, and the correct model name.
2. C3: elementary; novelty POSSIBLY-NEW but immaterial. Present it as a remark.
3. The README's "Kimura sinh model" wording must change: the sinh model is Montaldi–Tokieda's (1.12) family.
   Kimura's model is (per two secondary sources) the tanh model. This should be recorded in README §3 / FIXES
   by the coordinator, since this session's write scope is README Status + §2 only.
