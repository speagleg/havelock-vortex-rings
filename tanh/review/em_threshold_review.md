> Independent AI reader's report (2026-10-12) on Lemma 17 and Theorem 18 of `tanh/note_tanh.tex`. All GAP and WORDING items
> below were applied before commit; line references are to the draft it reviewed. The reader's scripts are kept locally.

# Referee report: "The energy--momentum threshold", Lemma lem:em, Theorem thm:em (tanh/note_tanh.tex, Sec. Large N)

**Verdict: ACCEPT WITH MINOR REVISIONS.**

Lemma lem:em and parts (a) and (c) of Theorem thm:em are correct. The core of part (b) is also correct: the window is
non-empty for every odd N >= 9, and its width is O(N^-3). One clause of the statement of (b) is not proved: "on
(xi_*, xi_EM) ... the critical block is indefinite". The proof covers only a neighbourhood of xi_*. My numerics support
the clause, but it must be weakened or labelled NUMERICAL. The proof also has two small gaps and there are a few wording
points. I found no false mathematical statement.

The reader's scripts (`ids.py`, `window.py`, `window2.py`, `phizeros.py`, `wconst.py`). They use
mpmath at 50-120 digits and `tanh/blocks.py` for the exact blocks.

## 1. Hand re-derivations (all confirmed)

* **gamma = R/(4 sinh^2 tau).** With theta = (1-x)/2, theta e^{x tau} - (1-theta) e^{-x tau} = sinh(x tau) - x cosh(x tau). Then
  4S gamma = sinh tau (sinh x tau - x cosh x tau) + e^{-tau} sinh x tau = cosh tau sinh x tau - x sinh tau cosh x tau = R,
  using sinh tau + e^{-tau} = cosh tau. Numerical relative error: 1.7e-55.
* **A, B, gamma^2/beta.** psi(0) - psi(x) = x^2/8 + (c - Q)/(4S), so 8S(psi(0) - psi(x)) = A. Also 8S beta = B and
  gamma^2/beta = R^2/(2SB), so 8SB G = AB - 4R^2.
* **Q^2 - R^2 = c^2 - x^2 S.** The cross terms cancel, leaving (c^2 - x^2 S)(cosh^2 x tau - sinh^2 x tau). Confirmed.
* **AB - 4R^2 = 2(c+1)^2 K.** Expanding gives
  AB - 4R^2 = x^2(1-x^2)S^2 - 2x^2 S + 2c(1-x^2)S + 4c(1+c) - 2(c+1)^2 Q.
  With S = (c-1)(c+1), divide by 2(c+1)^2. The middle terms collapse to c - x^2(c-1), which leaves
  c - Q - x^2(c-1) + x^2(1-x^2)(c-1)^2/2. This is exactly K. Hence beta G = (AB - 4R^2)/(64S^2) = (c+1)^2 K/(32 S^2).
  Numerical relative error: 4e-44 for beta G and 7e-58 for AB - 4R^2.
* **dK/dtau.** dQ/dtau = (1-x^2) sinh tau cosh x tau, and
  dK/dtau = sinh tau (1-x^2)[x^2(cosh tau - 1) - (cosh x tau - 1)]. This matches the note (finite-difference check 1e-47).
  The sign argument is right: the k = 1 terms are equal, and for k >= 2, x^{2k} < x^2 with tau^4 > 0. With K(x,0) = 0 this
  gives K > 0. A grid of tau = 10^{-4..2.5} by 201 x-values, at 120 digits, gives K > 0 and dK > 0 everywhere.
* **Final identity.** a_m b_m - g_m^2 = N^4 beta (Phi_N + G). Checked against `blocks.py` at 200 random (N <= 200, m, xi):
  relative error 1.7e-57.
* **R(1/N, tau) > 0.** R/(cosh tau cosh x tau) = tanh(x tau) - x tanh(tau). Strict concavity of tanh on [0, inf) with
  tanh 0 = 0 makes this positive for 0 < x < 1. Correct, and confirmed on the grid.

## 2. Theorem thm:em

**(a) Correct.** a_m >= a_{N/2} = N^2 Phi_N, and Lemma lem:em gives the determinants. The "iff" holds for every even N,
including N = 2, where m = N/2 is the only mode.

* "The thresholds coincide" needs more than "PD iff a > 0, stable iff a >= 0". With the inf-definitions, the two
  thresholds could differ if a_{N/2} touched 0 without changing sign. Prop. prop:evenN excludes this: a_{N/2} is strictly
  monotone with a simple zero (N >= 8) or no zero (N <= 6). See GAP-1.
* For N = 2, 4, 6 both thresholds are 0, because a_{N/2} > 0 at every radius. Prop. evenN gives limits 1, 1/2 and 1/9 times
  1/8. The statement is true but vacuous for these N. See WORDING-1.

**(b) Non-emptiness: correct for every odd N >= 9.**

* As r -> inf, 8 a_crit/N^2 -> 2 - 2/N^2 > 0. I derived this by small-tau expansion and checked it numerically.
* The planar limit -(N-1)(N-7)/N^2 is negative for N >= 9.
* So xi_* = sup{a_crit < 0} lies in (0,1), with a_crit(xi_*) = 0. This needs only the inf/sup definition, not uniqueness
  of the zero.
* At xi_*, the determinant equals -g^2 < 0, by R > 0. By continuity, xi_EM > xi_*.

**(b) "On (xi_*, xi_EM) ... the critical block is indefinite": not proved on the whole interval (GAP-2).** The proof shows
only two things. The critical determinant is negative on some right-neighbourhood of xi_*. And, by the definition of the
infimum, non-PD points accumulate at xi_EM. Two things are missing:

* Some other block m could fail definiteness inside the window. On the window a_m > 0, but Phi_N may be negative there, of
  size up to eps_N. So one needs G(x_m, tau) > eps_N(tau) for |x_m| >= 3/N. Since G ~ g_2 x^2 and eps_N ~ p_2/N^2, this is a
  quantitative inequality, roughly 9 g_2 > p_2, and it does not follow from anything proved.
* The non-PD set could be disconnected.

The note itself flags this assumption for W: "it assumes that the critical mode is the last to become definite". It is
inconsistent to use the same assumption silently in a PROVED statement.

My numerics support the claim. For N in {9, 11, 13, 15, 17, 21, 31, 61, 121, 241}:

* a_crit has exactly one sign change on tau in (0, 3 tau_inf].
* The critical determinant has exactly one sign change on [tau_Phi, tau_*].
* At 200 points inside the window, the critical block is indefinite and every other block is PD.
* At 400 points on [tau_Phi, tau_EM), all blocks are PD.

The underlying cause is that Prop. prop:kreincert (N <= 15) is certified, but the general-N clause is not.

**(b) O(N^-3) width: correct modulo GAP-3.** I checked each step:

* Phi_N(0+) = 1/4 - 1/(2N^2) > 0.
* Phi_N >= c_0 d - C_T/N^2 on (0, tau_inf - d] (Phi decreasing), so tau_Phi is within K/N^2 of tau_inf.
* tau_Phi <= tau_EM <= tau_*: Lemma lem:em gives PD on (0, tau_Phi), and xi_EM >= xi_*.
* Phi_N(tau_*) = -eps_N(tau_*), and |Phi_N'| >= c_0 near tau_inf. The mean value theorem then gives
  tau_* - tau_Phi <= eps_N(tau_*)/c_0.
* |e^{-2 tau_EM/N} - e^{-2 tau_*/N}| <= (2/N)(tau_* - tau_EM).

The only unproved ingredient is |Delta_N'| = O(N^-2). The proof of thm:largeN bounds only |Delta_N|.

The objects borrowed from the proof of thm:largeN are used consistently:

* Phi_N = Phi + Delta_N, with Delta_N = E_N - E_inf.
* a_crit/N^2 = Phi_N + eps_N, with eps_N = psi(0) - psi(1/N) >= 0.
* c_0 = |Phi'(tau_inf)|/2.

One caveat: "first zero tau_Phi" is the right notion. For small N, Phi_N's only zero lies far from tau_inf (tau = 9.57 at
N = 9 and 8.62 at N = 11, against tau_inf = 7.445). So the asymptotic argument says nothing for small N. That is fine for
an O(.) claim, but the remark's numbers at small N rely on locating the correct zero; `large_n.py` happens to find it.

**(c) Correct.** tau_EM lies in [tau_Phi, tau_*], both are tau_inf + O(N^-2), and cosh r = coth(tau/N) = N/tau + O(1/N). For
even N, (c) is immediate from (a) and thm:largeN.

## 3. Numbers in the remark (checked against large_n.out and my own code)

| Quantity | Note | large_n.out / mine | Status |
|---|---|---|---|
| bound at N = 241 | <= 14.7 | N^3(xi_Phi - xi_*) = 14.6228 (both) | true, but 14.62 rounds to 14.6 (WORDING-3) |
| actual N^3(xi_EM - xi_*), N = 241 | ~0.0035 | 0.0034787 (both) | OK |
| W | ~0.0037 | 0.00370893 (both, from tau_inf, R_x = tau cosh tau - sinh tau, Phi' = -0.0153242) | OK |
| N = 61 | 0.0028 | 0.00280024 | OK |
| N = 121 | 0.0032 | 0.00324945 | OK |

The "actual" xi_EM in `large_n.py` is computed as a root of the **critical** determinant only. My full all-blocks scan
confirms that it equals the true infimum for the sampled N.

The quoted 14.7 is the intermediate bound xi_EM <= xi_Phi. The final displayed estimates of the proof are looser. At
N = 241:

* (2/N)(tau_* - tau_Phi) N^3 = 15.56.
* (2/N) eps_N/c_0 N^3 = 31.1.

## 4. Issues and fixes

* **GAP-1 (thm:em(a), proof).** "Spectral stability is a_{N/2} >= 0, so the two thresholds coincide" needs a_{N/2} not to
  touch zero without a sign change.
  *Fix:* add "and by Prop. prop:evenN, a_{N/2} has no zero (N <= 6) or a single simple sign change (N >= 8), so both
  infima equal that zero (or 0)".
* **GAP-2 (thm:em(b), statement).** "On (xi_*, xi_EM) the ring is spectrally stable while the critical block is
  indefinite" is proved only on a right-neighbourhood of xi_*.
  *Fix (preferred):* restate as "the ring is spectrally stable on (xi_*, xi_EM), and the Hessian is not positive definite
  on a right-neighbourhood of xi_* (there the critical block is indefinite)". Optionally add: "that the critical block is
  indefinite on the whole window, and every other block PD, is CERTIFIED for N = 9..15 (Prop. kreincert) and NUMERICAL for
  N in {17, 21, 31, 61, 121, 241}".
  Make the same change in the status table (line ~110) if it implies the whole-window statement. The Krein-signature
  sentence in Sec. krein should likewise be read as N <= 15 or near xi_*.
* **GAP-3 (thm:em(b), proof).** |Delta_N'| = O(N^-2) is asserted, but thm:largeN proves only |Delta_N| = O(N^-2).
  *Fix (one line):* delta = L(tau/N)/N, with L the Langevin function. So delta' = L'(tau/N)/N^2 lies in (0, 1/(3N^2)], and
  Delta_N' = delta'(L(tau) - delta) + delta L'(tau). Hence |Delta_N'| <= (1 + tau/3)/(3N^2).
  Also say explicitly that |Phi'| >= (3/2)c_0 on a fixed neighbourhood of tau_inf, so |Phi_N'| >= c_0 there for large N.
* **WORDING-1 (thm:em(a)).** For N = 2, 4, 6 there is no threshold: xi_* = xi_EM = 0.
  *Fix:* "for N <= 6 the Hessian is PD at every radius; for even N >= 8, xi_EM = xi_*".
* **WORDING-2 (lem:em proof).** "Write c = cosh tau" reuses c, which is defined globally as cosh r = coth(tau/N).
  *Fix:* use C or ch for cosh tau.
  In the same proof, also cite Lemma lem:bpos (beta > 0) for "G >= 0 with equality only at x = 0" and for dividing by beta.
* **WORDING-3 (remark after thm:em).** "the proof gives N^3(xi_EM - xi_*) <= 14.7 at N = 241". The value is 14.62, and it
  comes from the intermediate bound xi_EM <= xi_Phi.
  *Fix:* "the bound xi_EM <= xi_Phi gives N^3(xi_EM - xi_*) <= 14.62 at N = 241 (the cruder estimates at the end of the
  proof give about 16 and 31)".
* **WORDING-4 (minor).** In (b), "its positive limit as r -> inf" can be stated explicitly: 8 a_crit/N^2 -> 2 - 2/N^2.
  The comparison table (line 769) "width ~1e-7" refers to N <= 15. For large N the width is about 0.0037/N^3, so qualify
  it with "(N <= 15)".
