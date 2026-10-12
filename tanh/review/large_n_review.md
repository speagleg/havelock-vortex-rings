> Independent AI reader's report (2026-10-12) on the large-N section of `tanh/note_tanh.tex`. All GAP and WORDING items
> below were applied before the section was committed; line references are to the draft it reviewed.

# Referee report: note_tanh.tex, Section "Large N" (sec:asym)

**Verdict: ACCEPT WITH MINOR REVISIONS.** No ERRORs found. Every PROVED statement was re-derived by hand and is
correct. Two small GAPs (one-line fixes) and several WORDING points follow.

Scripts: the reader's own (none import large_n.py); kept locally, not committed.
- `check1.py`: blocks from the statement of thm:blocks vs blocks.py; the symbol identity; the bpos identity; crit's v-formula;
  lambda'' bound; the evenN formula, monotonicity and limits; the odd-N P12+P3 formula and its planar limit; eq:phi; eta_inf; C1; W.
- `check2.py`: evenN monotonicity in detail (derivative sign).
- `check3.py`: thresholds for N = 8..801 from the statement-form blocks; scans of stability below/above tau_*; a_m >= a_crit.
- `horo.py`: independent horocyclic-row Bloch symbol. Lattice sums use closed-form Clausen sums for the 1/n^2 and 1/n^3 parts.
- `em.py`: limit-symbol definiteness near x = 0.

## Numerical summary (own code)
| check | result |
|---|---|
| statement-form blocks vs blocks.py (300 random N<=200, incl. xi near 0 and 1) | rel. err 2e-30 (dps 50) |
| thm:symbol identity, same sample | rel. err 3e-27 |
| lem:bpos identity, 2000 random (x,tau), tau in [1e-4, 40] | rel. err 5e-30; beta>0 throughout |
| lem:crit v-formula; max v/sinh h over 3000 samples | err 4e-50; max = -1.73 < 0 |
| sup lambda'' | 1/3 (attained as t->0) |
| prop:evenN formula / limits / dF/deta < 0 | err 1e-46; limits exact; derivative strictly negative (N=4..12) |
| odd formula 8a/N^2 = P12+P3-3/N^2, planar limit -(N-1)(N-7)/N^2 | err 4e-48; limits exact |
| eq:phi | err 2e-47 |
| eta_inf, kappa_inf, ell_* | 3.72266962247450934361..., 0.13431221427262813906..., 0.84390853129... (all digits agree) |
| Phi'(tau_inf) | -0.0153242426 (nonzero) |
| C1 even / odd | 0.5412188400 / 0.4010382403 (agree) |
| (cosh r_* - kappa N) N at N=800 / 801 | 0.541220596 / 0.401039312 (agree) |
| stable for tau in (0,tau_*), unstable just above, N in {8,...,801} | yes for all |
| (ell_N - ell_*) N^2 | about -20 (O(N^-2) confirmed; ell_N < ell_*) |
| min (a_m - a_crit)/N^2, noncritical m, N<=60 | 9.5e-5 > 0 |
| horocycle row symbol vs (alpha,beta,|gamma|), 6 spacings x 4 wavenumbers | max discrepancy 2.7e-16 |
| W | 0.0037089 (agrees) |

## Statement-by-statement

**Setup.** tau = N artanh(1/c) and c = coth(tau/N) are correct, and xi^N = e^{-2tau}.

**thm:symbol (PROVED): correct.** Re-derived W_m/N, W_{N-m}/N and NW_N = N^2/(4 sinh^2 tau). The identity
theta e^{x tau}+(1-theta)e^{-x tau} = cosh x tau - x sinh x tau holds, and sinh tau + e^{-tau} = cosh tau gives Q.
V'' = Nc coth tau - c^2/2 - 1/2 is checked. N=2 and all m are covered.

**lem:bpos (PROVED): correct.** b cosh a + a cosh b = 2tau(cosh tau cosh x tau - x sinh tau sinh x tau) = 2 tau Q. The
sinh^2 rewriting is correct, f(2tau) = sinh^2 tau/(2 tau), and the identity follows. The superadditivity argument is valid: sinh(u)/u is
strictly increasing, so f(t)/t is too, and a, b > 0 for |x| < 1. The endpoints x = ±1 are m = 0, N, which are excluded.

**lem:crit (PROVED): correct.** I checked each step:
- 4 sinh^2 tau psi_x = sinh tau v(h): verified by hand (x sinh tau = h S(tau)).
- T - h coth h = D + I, since (t coth t)' = lambda' + t lambda''. The reduction to I < 1 + e^D - D is equivalent.
- The RHS is increasing in D and > 2 when h < tau, because lambda' > 0 means D > 0. At h = tau, I = 0 < 2.
- lambda'' <= 1/3: sinh^2 t - t^2 = sum_{j>=2} 2^{2j-1}t^{2j}/(2j)!, and t^2 sinh^2 t has t^{2j} coefficient
  2^{2j-3}/(2j-2)!. The ratio (2j)(2j-1)/12 is >= 1 for j >= 2 (equality at j=2), and all coefficients are >= 0, so the bound holds.
- The integral bound gives tau^2/6 <= 1/2 (tau <= sqrt3), else 1/2 + ln(tau/sqrt3), and this is <= 2 iff tau <= sqrt3 e^{3/2} = 7.7622.
  Since RHS > 2 strictly, the strict inequality follows.
- tau > 7.76, h >= tau/2: I <= ln 2. h < tau/2: e^D >= S(tau)/S(tau/2) = cosh(tau/2). The bound
  1 + E - ln E >= 1 + e^{tau/2}/2 - tau/2 is valid because each piece is bounded separately. That quantity exceeds 1/2 + ln(tau/sqrt3) at 7.76
  (21.3 vs 2) and the gap is increasing.
- The consequence for a_m (evenness, then minimal |x| at m = floor(N/2), N - floor(N/2)) is correct, including the equality case.

**cor:alln (PROVED): correct.** Prop. lin plus b_m > 0 gives: stable iff all a_m >= 0 iff a_crit >= 0. If a_crit < 0, then
sqrt(-ab) is real and nonzero, so there is exponential growth. At a_crit = 0 the ring is spectrally (not necessarily linearly) stable. That is consistent with the
definitions.

**prop:evenN (PROVED): correct, with one small gap.** Formula: at x=0, psi + 1/(4 sinh^2 tau) = 1/8 + 1/(8 sinh^2 eta)
(using 1 + cosh tau = 2cosh^2 eta). The rest follows as stated. Monotonicity: d/dv ln(tanh(kv)/tanh v) < 0 for k = N/2 > 1. Also
mu <= coth eta (v coth v increasing, v = 2eta/N <= eta), and (coth eta - mu)' >= 0. The first term's derivative is *strictly*
negative, so the derivative of the whole is strictly negative and the sign change really is simple. Limits are checked:
eta->0 gives 2 - 4/N^2; eta->inf gives 4/N - (1-2/N)^2 - 4/N^2 = -(N^2-8N+8)/N^2.

**eq:phi and monotonicity of Phi: correct.** 8Phi = 2(coth eta + tanh eta)/eta - 1/eta^2 - coth^2 eta = 2 tanh eta/eta - L^2.
The limits are 1/4 and -1/8, and the zero is unique.

**thm:largeN (PROVED): correct, with one small gap.**
- (d1) Uniformity as tau -> 0: delta = (u coth u - 1)/tau is in (0, tau/(3N^2)). The formula
  Delta_N = delta(coth tau - 1/tau - delta/2) - 1/(2N^2) is checked, and |Delta_N| <= C_T/N^2 on (0,T] holds uniformly. The stable direction uses only
  epsilon_N >= 0 (Lemma crit), so no uniformity of the psi_xx bound is needed near tau = 0. Good.
- (d2) Threshold definition: with tau_* = sup{tau0 : stable on (0,tau0)}, stability on (0, tau_inf - d] gives
  tau_* >= tau_inf - d, and instability at the single point tau_inf + d gives tau_* <= tau_inf + d. This is correct, and it needs no
  uniqueness of the threshold (relevant for odd N).
- (d3) Stability needs only a_m >= 0 and b_m > 0: correct for *spectral* stability via Prop. lin (g_m only shifts frequencies).
- (d4) Conversions are correct: cosh r_* = coth(tau_*/N) = N/tau_* + O(1/N); sinh r_* = 1/sinh(tau_*/N); ell_N = 2pi/tau_* + O(N^-2).
  Numerically, (ell_N - ell_*)N^2 is about -20.
- Remarks: r_* = ln(2 kappa N) + o(1) is correct (in fact O(N^-2)), and 1 - xi_* ~ 14.8907/N is correct.

**Second order (DERIVED): derivation correct.** delta tau = D/(N^2 |Phi'|); cosh r_* = N/tau_* + tau_*/(3N) gives C1 as stated.
Labelled correctly (no remainder bounds).

**Physical interpretation:** wording is mostly at the right strength. My independent check confirms the horocycle identification to 3e-16.

## Issues

1. **GAP (thm:largeN proof).** "let c_0 = 1/2|Phi'(tau_inf)| > 0". Strict monotonicity of Phi does not imply
   Phi'(tau_inf) != 0. Fix: add "Phi' < 0 everywhere, since (tanh eta/eta)' = (2eta - sinh 2eta)/(2eta^2 cosh^2 eta) < 0 and
   (L^2)' = 2 L L' > 0". (Numerically Phi'(tau_inf) = -0.015324.) Add the same word to the eq:phi paragraph: "strictly decreasing, with Phi' < 0".
2. **GAP (prop:evenN).** The statement says "every even N <= 6 is stable at every radius", but the monotonicity argument needs k = N/2 > 1, so it excludes
   N = 2. Fix: add "For N = 2, mu = coth eta and the right side is identically 1."
3. **WORDING (prop:evenN proof).** "for k = N/2 > 1" reuses k, which thm:symbol defines as c/N. Rename it (e.g. p = N/2).
4. **WORDING (thm:largeN, "In words").** "neighbouring vortices on the threshold ring are ell_N -> ell_* curvature radii apart."
   ell_N is the *arc length along the ring* per vortex. The geodesic distance between neighbours tends to 2 arsinh(ell_*/2) = 0.8207, not 0.844.
   Fix: "the arc length between neighbouring vortices along the threshold ring is ell_N -> ell_* curvature radii".
5. **WORDING (Physical interpretation, "A critical spacing").** "becomes spectrally stable once neighbouring vortices are more than ell_* R apart"
   holds only asymptotically. Since ell_N ~ ell_* - 20/N^2 < ell_*, finite rings are already stable slightly below ell_*.
   Fix: add "asymptotically (the finite-N threshold spacing ell_N approaches ell_* from below, ell_N = ell_* + O(N^-2))".
6. **WORDING (odd-N remark after prop:evenN).** "below 0.61% of |P12'| for N >= 9 (NUMERICAL)". large_n.out samples only
   N in {9,11,15,21,51,101,1001} on a finite v-grid. Fix: "for the sampled N in {9,...,1001}". Also note in the text that P12 is strictly decreasing
   for odd N by the same proof as prop:evenN (k = N/2 > 1 still holds). That is proved; only the P3 comparison is numerical.
7. **WORDING (after eq:kappa).** "eta_inf ~ 2+sqrt3 and kappa_inf ~ 1-sqrt3/2, within 3.4e-4". Only kappa is within 3.4e-4.
   eta_inf - (2+sqrt3) = -9.4e-3. Fix: "kappa_inf ~ 1 - sqrt3/2 within 3.4e-4 (eta_inf within 9.4e-3)".
8. **WORDING (cor:alln).** "This removes the restriction N <= 16 from the statements ... of Thm cert(b) and Conjecture conj:deg." The current
   conj:deg no longer contains those statements. Fix: "... and from an earlier version of Conjecture conj:deg".
9. **WORDING (horocycle paragraph).** "The critical mode theta = 1/2 displaces neighbours in opposite normal directions." The
   theta = 1/2 Fourier mode alternates both components. The negative stiffness is in the normal entry, but the growing eigenvector also has a
   tangential part, and that part makes the vortices bunch in pairs. Fix: "the negative stiffness of the theta = 1/2 mode is in the alternating
   normal displacement; the growing eigenvector couples it to alternating tangential displacement, which pairs the vortices".
10. **WORDING ("Why curvature stabilises").** The screening explanation is heuristic. Mark it as such ("Heuristically, ...").
11. **Suggestion (not an error).** The horocycle identification looks provable. In horo.py the slowly decaying parts are exact
    Clausen sums (sum cos(n x)/n^2, sum sin(n x)/n^3), and the remainder matches the limit symbol to 3e-16. A Poisson-summation or
    partial-fraction derivation (sum 1/(n^2 l^2 + 4) via coth) would let the label be upgraded from NUMERICAL. The current label is honest.

## (e) PROVED vs numerical
No statement labelled PROVED depends on numerical evidence. thm:symbol, lem:bpos, lem:crit, cor:alln, prop:evenN and
thm:largeN are fully analytic given thm:blocks, once issue 1 is fixed. The constants eta_inf, kappa_inf and ell_* are numerical values of an analytically defined root. The
NUMERICAL labels on odd-N uniqueness, the horocycle identification, the energy-momentum/Krein-window asymptotics and the
second-order remainders are correctly placed. The energy-momentum claim checks out: the x^2 coefficient of alpha beta - gamma^2 at
tau_inf is +0.01485, and the grid minimum is 1.19e-5 > 0.
