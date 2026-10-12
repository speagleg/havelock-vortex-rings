# Literature notes

Two literature checks, recorded as they were made: one for the tanh-model note (`tanh/note_tanh.tex`,
2026-10-11) and an earlier one for the sinh-model note (`sinh/note.tex`, 2026-10-07). Local copies of the papers
are not committed (copyright).

## Tanh model (`tanh/note_tanh.tex`), 2026-10-11

Question: does any accessible source already contain ring-stability thresholds for the **tanh model**
(pair interaction −ln tanh(d/2) on H², the Hamiltonian attributed to Kimura 1999)? "Full text" means the whole
PDF was obtained and the relevant sections were read; "abstract only" means nothing beyond the abstract was seen.

### Verdict (short)

- **OVERLAP FOUND: Boatto (2008) treats the tanh model numerically.** Boatto, "Curvature perturbations and
  stability of a ring of vortices", DCDS-B 10 (2008) 349–375 (open full text from the publisher). For the
  hyperbolic surface she uses the Green's function "ΨK = −(1/2π) ln tanh(√|K| r/2)" (her eq. after (14)) and the
  Hamiltonian "HK<0 = c Σ ΓαΓβ ln((dαβ − 1)/(dαβ + 1))" with c = −1/4π (eq. (17)), where dαβ = cosh of the
  geodesic distance in her Darboux coordinates (11), (19). That is exactly our model. She writes the Hessian of
  the rotating-frame Hamiltonian as S̃ = [[Q, G], [Gᵀ, P]] with "Q and P … symmetric circulant matrices, while G
  is an asymmetric circulant matrix" (§5, item b), i.e. the same Fourier-block structure with a
  radial–tangential coupling that our Theorem 1 has. Her block eigenvalues are (44)–(45). She then says:
  "we (numerically) monitor the maximum and the minimum eigenvalues of S̃ … For a given N, Figure 5 shows that
  stability is guaranteed for sufficiently low K. In particular for N ≥ 7, we can only assure stability for
  sufficiently negative K." Figure 5a plots the minimum Hessian eigenvalue against K ∈ [−0.6, 0.6] (ring of
  geodesic radius r₀ = 1) for N = 5,…,9. Her explicit matrix entries (39)–(40) are written out only for the
  spherical case ("The calculations for the hyperbolic case are done is a similar way").
  **What overlaps with our note:** the model; the stability question as a function of curvature (equivalently,
  of radius); the circulant/Fourier block structure with a non-zero skew coupling G; and the qualitative
  conclusion that, on H², rings with N ≥ 8 are Hessian-definite only for sufficiently negative curvature
  (sufficiently large radius). Her criterion is Hessian definiteness (Dirichlet), which is our energy–momentum
  threshold ξ_EM, not the spectral threshold.
  **What we did not find in it:** closed forms of the block entries (hers are finite trigonometric sums
  evaluated numerically); any threshold value in closed form or to more than plotting accuracy; algebraicity,
  degrees, minimal polynomials, palindromy or the ξ*(8) unit; the separation of spectral and
  energy–momentum thresholds (the odd-N Krein window); the statement "N ≤ 7 stable at every radius" (her Fig.
  5a is consistent with it on the range plotted); large-N asymptotics.
  **Fig. 5a crossings (resolved).** Read by eye, the zero crossings in her Fig. 5a are at K ≈ −0.25 (N = 8) and
  K ≈ −0.40 (N = 9), while our thresholds predict K = −r*² = −0.282 and −0.489. A reimplementation of her
  stated equations reproduces our ξ*(8) to 1e-41 and our ξ_EM(9) exactly, so her equations and our thresholds
  agree. Her plotted K < 0 curves appear to be matched instead by the pair law −ln sinh D, i.e. her eq. (17) with the
  sign of one term flipped; her qualitative conclusions are unaffected. Details: `tanh/boatto_check/RESULT.md`.
- **Montaldi & Nava-Gaxiola (2014)** and **Nava-Gaxiola's thesis (2013)**: full text; same (tanh/Kimura)
  Hamiltonian; only N = 2, 3 analysed; no ring thresholds. Both send the reader to Hwang–Kim (2013) for rings.
- **Hwang & Kim (2013)**: **abstract only**. Its abstract says "the configuration of rings of vortices is
  investigated in detail"; the thesis says it gives "relative equilibria conditions for rings of vortices" with
  the same Hamiltonian as ours. Whether it contains stability thresholds is **unknown**. This is the main
  unread item.
- **Kimura (1999)**: **abstract only**. The abstract mentions the Hamiltonian formalism on S² and H² and
  dipole motion; nothing about rings.
- **Boatto & Simó (2008)** "Thomson's heptagon: a case of bifurcation at infinity", Physica D 237: **not
  obtained** (no open copy found). From Boatto (2008, Theorem 5.3, "see Boatto and Simó [5]") it concerns the
  sphere and plane, not H².
- **Boatto & Koiller**, "Vortices on closed surfaces", arXiv:0802.4313: full text searched; closed surfaces
  only; H² appears only in a historical sentence; no ring thresholds.
- **Montaldi & Tokieda (2013)**, arXiv:1210.5662v1, §4.2, p. 25 (full text): for G = log|z−w|²/|1+λzw̄|²,
  "On the hyperbolic plane λ < 0 this reduces to the Hamiltonian adopted by Kimura [10]" and "The calculations
  based on this option of Green's function get so cumbersome that the stability problem seems no longer
  tractable analytically. It seems likely that the results are similar to those in section 3, though the details
  of where the bifurcations occur will differ." They cite Boatto (2008) as [1] in the introduction ("part of it
  has leaked into the literature [1]") but do not mention her hyperbolic computation in §4.2.

**Consequence for the note.** "Likely new" can be claimed only for the exact/algebraic content (closed-form
rational blocks, certified thresholds, minimal polynomials, palindromy, the unit, the Krein window, the
asymptotic constant), and only subject to Hwang–Kim (2013), which remains unread. The stability question
itself, the block structure and the qualitative N ≥ 8 picture are in Boatto (2008) and must be credited.

### Per-source record (tanh model)

| Source | Accessible? | Tanh-model ring thresholds? |
|---|---|---|
| Boatto, DCDS-B 10 (2008) 349–375, doi:10.3934/dcdsb.2008.10.349 | **Full text** (publisher PDF, bronze OA via aimsciences.org) | **Yes, numerically** (Fig. 5a, N = 5..9, Hessian definiteness). Same model; Q/P/G circulant block structure. No closed forms. See verdict. |
| Montaldi & Nava-Gaxiola, J. Math. Phys. 55 (2014) 102702, arXiv:1403.2138 | **Full text** (arXiv) | No. N = 2, 3 only. Eq. (12) is Kimura's Hamiltonian: (⟨X,Y⟩+1)/(⟨X,Y⟩−1) with ⟨X,Y⟩ = −cosh d (eq. (2)), i.e. tanh²(d/2). Intro: "a more in-depth study is provided by Hwang and Kim [8]. In that paper, they present conditions for relative equilibria for rings of vortices on the hyperbolic plane". Credits Boatto [5] and MT [18] for "for N > 7 stability only occurs for negative curvature". |
| Nava-Gaxiola, PhD thesis "Point vortices on the hyperboloid", Manchester 2013, eprints.maths.manchester.ac.uk/2022 | **Full text** (searched for ring/polygon/Hwang; intro read) | No. N = 2, 3 only. "In contrast with [11], the Hamiltonian used in [12] [Hwang–Kim 2013] coincides with the our choice of Hamiltonian." |
| Hwang & Kim, J. Math. Phys. 54 (2013) 063101, doi:10.1063/1.4811454 | **Abstract only** (OpenAlex/Semantic Scholar: closed; AIP behind Cloudflare; not on arXiv; author page has no PDF) | **Unknown.** Rings "investigated in detail"; secondary sources mention only relative-equilibrium conditions. |
| Hwang & Kim, J. Geom. Phys. 59 (2009) 475–488 | Not obtained (ScienceDirect returns an HTML challenge) | Per MNG and the thesis: three vortices, collapse; different Hamiltonian. Not a ring paper. |
| Kimura, Proc. R. Soc. Lond. A 455 (1999) 245–259 | **Abstract only** | Abstract: Hamiltonian formalism on S²/H², first integrals, dipole geodesics. Nothing on rings in the abstract. |
| Boatto & Simó, Physica D 237 (2008) | Not obtained | Sphere/plane (per Boatto 2008). |
| Boatto & Koiller, arXiv:0802.4313 | Full text (searched) | No. |

Searches (2026-10-11): arXiv API (title/author queries for Hwang–Kim, "hyperbolic sphere" vortices); OpenAlex
and Semantic Scholar by DOI for open-access locations; web search for Hwang–Kim 2013 full text; S.-C. Kim's
publication page (cau.ac.kr/~kimsc). Not searched: MathSciNet/zbMATH reviews of Hwang–Kim 2013, and papers
citing Boatto (2008) or Hwang–Kim (2013) (e.g. Newton, "Point vortex dynamics in the post-Aref era", 2014).

---

## Sinh model (`sinh/note.tex`), 2026-10-07

Scope: the sinh note's core claims C1–C4 (N-ring thresholds on H²). For each source: whether it was obtained and
read, and what it says about N-ring stability thresholds on H².

### Verdict (short)

- **C1 (mode decomposition, C₁(ξ) = (N−1)(1+ξ²)/(1−ξ)²) is KNOWN.** It is Montaldi & Tokieda (2013),
  §2.2 eq. (2.6), together with the stability criterion (2.11) / Theorem 2, for the
  "background-vorticity" Green's function (1.12) continued to λ < 0. On the hyperbolic plane that model is
  exactly the sinh model of the note (pair term ln[|z−w|²/((1−|z|²)(1−|w|²))] = 2 ln sinh(d/2) + const).
- **C2 (marginal threshold) is KNOWN.** Montaldi–Tokieda (2.11) reads
  (1+λ²r₀⁴)/(1+λr₀²)² > ⌊n²/4⌋/(2(n−1)). With ξ = −λr₀² this is (N−1)(1+ξ²)/(1−ξ)² > ⌊N²/4⌋/2 = T,
  i.e. C1/C2. Their table of bifurcation values λr₀² for n = 8…13 is −0.0627, −0.101, −0.143, −0.172,
  −0.202, −0.225, i.e. −ξ*(N) (compare `sinh/verify_core.out`: 0.06275, 0.10102, 0.14286, 0.17157, 0.20241,
  0.22515). Their worked example "n = 15 stable for r₀² > (2−√3)/|λ|" is ξ*(15) = 2−√3.
  The **cosh² r* = q(N) form** is a one-line rewriting of (2.11) through (1+ξ²)/(1−ξ)² = (1+cosh²r)/2.
  It does not appear in that form in the sources read. It is not a new result, and should be presented as a
  reformulation.
- **C3 (unit classification N ∈ {8,9,11,15,23}) was not found in any source read.** It is an elementary
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
  So the name "Kimura's sinh model", used for the sinh Hamiltonian in this project's earlier drafts, is very
  probably a misattribution. **Kimura 1999 itself could not be obtained** (royalsocietypublishing.org returns
  403/Cloudflare; DOI 10.1098/rspa.1999.0311 confirmed via Crossref), so this stays a strong secondary-source
  finding until someone checks the original. In the sinh note the model is therefore called the **sinh model**
  (= Montaldi–Tokieda's family (1.12) at λ < 0), never "Kimura's model".
  Consequence: the threshold problem for the model that the literature attributes to Kimura is the **tanh**
  problem, treated in `tanh/note_tanh.tex`. At the time of this check its closed form was open, and
  Montaldi–Tokieda explicitly left it untreated.

### Per-source record

| Source | Obtained / read? | What it says about N-ring thresholds on H² |
|---|---|---|
| **Montaldi, J. & Tokieda, T.** "Deformation of geometry and bifurcations of vortex rings", in *Recent Trends in Dynamical Systems*, Springer Proc. Math. Stat. 35 (2013) 335–370; arXiv:1210.5662 | **Yes, full text (arXiv PDF)** | Smooth family M_λ (curvature 4λ), symplectic form 2dx∧dy/σ², σ = 1+λ\|z\|². Green's function (1.12) log[\|z−w\|²/(σ_zσ_w)] (sinh model at λ<0). §2: ring angular velocity (2.1); Hessian in (r,θ) block-diagonal with no r–θ coupling (2.2); Fourier eigenvalues (2.6) ε_r^(ℓ) = κ²/(4πr₀²)[2(n−1)(1+λ²r₀⁴)/σ² − ℓ(n−ℓ)], ε_θ^(ℓ) = κ²ℓ(n−ℓ)/4π; symplectic slice (2.9); **Theorem 2**: Lyapunov-stable (modulo G_µ) if (2.11) holds; Remark 4: linearly unstable if (2.11) is reversed. "Hyperbolic planes: when n ≤ 7 the ring is always stable. When n > 7, it is stable (for a given λ) for a sufficiently large." Table of bifurcation values for n = 4…13. §3: bifurcations at the threshold and nonlinear stability of the heptagon (Cor. 7); Conjecture 13 (n-gon stable at the bifurcation point). §4.1: the criterion for G = log\|z−w\|² (1 + (5/6)λ²r₀⁴ in place of 1+λ²r₀⁴). §4.2: the tanh/Kimura Green's function is "no longer tractable analytically". Uses identity (2.4) Σ cos(2πℓj/n)/(1−cos(2πj/n)) = (n²−1)/6 − ℓ(n−ℓ), citing Hansen's table. **⇒ C1, C2 KNOWN.** |
| **Montaldi, J. & Nava-Gaxiola, C.** "Point vortices on the hyperbolic plane", J. Math. Phys. 55, 102702 (2014); arXiv:1403.2138 | **Yes, full text** | Only N = 2, 3 are analysed (relative equilibria, SL(2,ℝ) momentum types, stability). The intro summarises Boatto and Montaldi–Tokieda: "for a given radius of ring, the ring becomes Lyapunov stable as the curvature decreases and for N > 7 stability only occurs for negative curvature". Rings on H² are "mentioned briefly in [3]" (Aref et al. 2003); "a more in-depth study is provided by Hwang and Kim [8]", which presents conditions for relative equilibria of rings. Eq. (12): Kimura's Hamiltonian in tanh form (see above). No explicit H² ring threshold formula. |
| **Nava-Gaxiola, C.** PhD thesis "Point vortices on the hyperboloid", Manchester 2013 (eprints.maths.manchester.ac.uk/2022) | Yes (intro, ch. 2, ch. 4 read; searched for "ring") | Same Hamiltonian (4.18) attributed to Kimura (tanh form). States that Hwang–Kim 2009 used a *different* Hamiltonian from Hwang–Kim 2013, and that the 2013 one coincides with Kimura's (tanh). No ring thresholds. |
| **Kimura, Y.** "Vortex motion on surfaces with constant curvature", Proc. R. Soc. Lond. A 455 (1999) 245–259, doi:10.1098/rspa.1999.0311 | **No** (403 / Cloudflare; no open copy found) | Per the two Montaldi-group sources: unified S²/H² formulation, tanh-type Hamiltonian on H², dipole motion (geodesic conjecture). Nothing on ring thresholds is reported by the secondary sources. **Unverified (original not read): nothing is quoted from Kimura directly.** |
| **Hwang, S. & Kim, S.-C.** "Point vortices on hyperbolic sphere", J. Geom. Phys. 59 (2009) 475–488, doi:10.1016/j.geomphys.2009.01.003 | **No** (ScienceDirect 403) | Per Montaldi–Nava-Gaxiola and the thesis: fixed and relative equilibria of three vortices, collapse. Different Hamiltonian from the 2013 paper. |
| **Hwang, S. & Kim, S.-C.** "Relative equilibria of point vortices on the hyperbolic sphere", J. Math. Phys. 54, 063101 (2013), doi:10.1063/1.4811454 | **Abstract only** (Crossref; AIP 403) | Abstract: "…In particular, the configuration of rings of vortices is investigated in detail." Per secondary sources: conditions for ring relative equilibria, any 2-vortex configuration is a relative equilibrium, Kimura's (tanh) Hamiltonian. **Unknown whether it contains a stability threshold.** If it does, it is for the tanh model and would bear on the tanh note, not C1. Needs a human with journal access. |
| **Boatto, S.** "Curvature perturbations and stability of a ring of vortices", DCDS-B 10 (2008) 349–375 (preprint mp_arc 06-311) | Abstract read; preprint PDF downloaded but its text is unextractable (Type-3 fonts); TeX source link returns HTTP 500 | Abstract: stability of N-gons as a function of N and curvature K, Dirichlet criterion; "Stability ranges are the K-intervals for which the Hessian … is positive or negative definite." Montaldi–Tokieda and Montaldi–Nava-Gaxiola credit it, jointly with Montaldi–Tokieda, for "N > 7 stable only for negative curvature". Its explicit formula and Green's function are **unverified**. It may well contain C1/C2 earlier than 2013. Treat as a probable prior source. **Update 2026-10-11:** the published version is open (see the tanh section at the top). On H² it uses the *tanh* Green's function, not the sinh one, and its results are numerical (Fig. 5a); so it is prior work for the tanh note, not for the sinh closed form. |
| **Polvani, L. M. & Dritschel, D. G.** J. Fluid Mech. 255 (1993) 35–64 | Not read (paywalled) | Cited by Montaldi–Tokieda as the first linear-stability study of rings on **S²**. Nothing on H². |
| **Boatto, S. & Cabral, H. E.** SIAM J. Appl. Math. 64 (2003) 216–230 | Not read | Cited by Montaldi–Tokieda as the nonlinear stability of a latitudinal ring on **S²** ("for spheres is due to [2]"). Nothing on H². |
| **Laurent-Polz, Montaldi, Roberts** J. Geom. Mech. 3 (2011) 439–486 (arXiv math/0402430) | Located, not read in detail | S² symmetric relative equilibria. Not H². |
| **Boatto, S. & Koiller, J.** "Vortices on closed surfaces", arXiv:0802.4313 | Yes (searched) | General Green's-function/Robin formulation on surfaces. No H² ring thresholds. |
| Searches | — | Web: "point vortex ring hyperbolic plane stability threshold N>7 relative equilibria" (extended); "Montaldi Tokieda Deformation of geometry…"; "Boatto Curvature perturbations…"; "Hwang Kim Relative equilibria … ring stability"; "polygonal/regular polygon point vortices constant negative curvature Lyapunov stability ring radius threshold explicit" (extended); "vortex ring hyperbolic plane stability threshold 8−3√7 / algebraic unit / cosh²". Crossref for DOIs. None of the hits states the cosh² form or the unit classification. |

### Consequences for the sinh note

1. C1 and C2: relabel **KNOWN** (Montaldi–Tokieda 2013, eq. (2.6), (2.11), Thm 2 and table; probably also
   Boatto 2008, unverified). The note is **expository** for these. Its value is a short self-contained
   derivation in log-polar coordinates (Hessian diagonal per Fourier mode), a high-precision independent
   dynamics check, and the correct model name.
2. C3: elementary; novelty POSSIBLY-NEW but immaterial. Present it as a remark.
3. Naming: the sinh model is Montaldi–Tokieda's (1.12) family, not "Kimura's sinh model". Kimura's model is (per
   two secondary sources) the tanh model. The repository uses "sinh model" and "tanh model" accordingly.
