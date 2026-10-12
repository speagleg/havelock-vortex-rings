/-
Havelock — Montaldi–Tokieda Conjecture 13 (vortex ring stable AT its bifurcation point):
the exact-algebra layer of the proof in conj13/paper/main.tex.

Conventions as in the paper: `H = −Σ log(|z_i−z_j|²/(σ_iσ_j))`, `σ = 1 + λ|z|²`, ring radius
`r₀ = 1`, so `b = λ`. "Paper (k)" refers to the numbered displays/theorems of main.tex
(by LaTeX label: Prop. blocks = `prop:blocks`, Prop. spectral = `prop:spectral`, Thm even = `thm:even`,
Thm odd = `thm:odd`, Def. c4 = `def:c4`, eq. (fejer), eq. (bracket), eq. (bn)).

FORMALISED (no `sorry`; axioms: propext, Classical.choice, Quot.sound only)

§1 The bifurcation point (paper eq. (bn) = MT (3.6) and the explicit `b_n` below it)
* `bn n`              closed form: even `(1+2t−t²)/(t²+2t−1)`, `t = √(n−1)`; odd `(2−s)/(2+s)`, `s = √(n−3)`.
* `bn_spec`           ∀ n ≥ 4: `−1 < b_n < 1` and `(1+b_n²)/(1+b_n)² = ⌊n²/4⌋/(2(n−1))`.
* `bn_unique`         ∀ n ≥ 4: `b_n` is the ONLY root of (3.6) in `(−1, 1]` (via `rfun_injOn`).
* `bn_even_paper`, `bn_odd_paper`: equal to the paper's `(4(n−2)√(n−1)−n²)/(n²−8n+8)` and
                      `(4√(n−3)−(n+1))/(n−7)`; `bn_seven : b₇ = 0`.
* `bn_four/five/six/eight/ten`: `2−√3`, `3−2√2`, `9−4√5`, `−(8−3√7)`, `−1/7`.
* `bn_sign`           ∀ n ≥ 4: `b_n > 0` for `n ≤ 6` (sphere), `b_n < 0` for `n ≥ 8` (H²).
§2 Fourier blocks / reduced space (paper §3: Prop. blocks, item (c) eq. (bracket), Prop. spectral and the remark after it)
* `phi_crit_and_C1`   the one-body function `φ(s)` of §3.1 (real `exp`/`log`, genuine derivatives):
                      with `ω₀ = −(n−1)(1−b²)/2`, `φ'(0) = 0` and `φ''(0) = C₁ = 2(n−1)(1+b²)/(1+b)²`.
* `C1_eq_iff`         `C₁ = ⌊n²/4⌋ ⇔ (3.6)`.
* `am_factor`, `am_nonneg_and_kernel`, `am_next`: ∀ n, m ≤ ⌊n/2⌋:
                      `⌊n²/4⌋ − m(n−m) = (⌊n/2⌋−m)(⌊n/2⌋−m + n mod 2) ≥ 0`, `= 0 ⇔ m = ⌊n/2⌋`,
                      next value `1` (even n) / `2` (odd n).
* `V1'_form`, `V1'_kappa_rho`: the mode-1 Hessian restricted to `V₁' = V₁ ∩ ker dJ₊` (given in
                      coordinates by `Ax − Bv = 0, Ay + Bu = 0`) is `κ·Id`, `κ = (a₁+b₁ρ²)/(1+ρ²)`.
* `V1'_n4`, `V1'_n5`  `κ = 3/2` at n = 4 and `κ = 8/3` at n = 5 (paper, after Prop. spectral).
* `bracket_value`     `{J₁,J₂}(x₀) = n/2 − λJ₃(x₀) = n(1−b)/(2(1+b)) ≠ 0` (algebra only).
§3 Even n = 2M, Thm even
* `h2_eq_zero_iff`    the stated `h₂` vanishes ⇔ (3.6) for even n.
* `c4even_diff`       exact identity (no threshold assumed):
                      `h₄(2/(1+b)²)⁴/n² − c4even(n) = (M²(1+b)² − 2(2M−1)(1+b²))·S/(8M(2M−1)(1+b)⁴)`.
* `c4even_reduction`  on (3.6): `c₄ = (n−2)(m⁴+40m²+48m+16)/(192n(n−1))`, m = n−2 (all M, both roots).
§4 Odd n, Thm odd
* `f2_eq_zero_iff`, `odd_rel_iff`: the stated `f₂` vanishes ⇔ `(n−7)(1+λ²)+2(n+1)λ = 0` ⇔ (3.6) for odd n.
* `Q_values`          `Q = √(w(1−λw))` at `w₀ = 1/(1+λ)`: `Q = 1/(1+λ)`, `Q' = (1−λ)/2`,
                      `Q'' = −(1+λ)³/4` (genuine derivatives of `Real.sqrt`).
* `J1con_iff`, `E_vertex`, `a2_pos`: eliminating `s` by the `J₁` constraint leaves
                      `a₂(p−p*)² + const` with `a₂ = n(n−1)(1−λ)²(1+λ)²/8 > 0`.
* `c4odd_reduction`   on the odd relation: `min{64/(n²(1+λ)⁸)·E(p,s) : J₁ = 0}` exists (IsLeast) and
                      equals `(n−1)(n−3)(n²−9n+38)/(128n)` — for all real n > 1 incl. n = 7, λ = 0.
§5 At b_n, positivity, Table 1
* `c4_even_at_bn`, `c4_odd_at_bn`, `c4_at_bn`: ∀ n ≥ 4, the paper's defining expression evaluated
                      at the true irrational `b_n` equals the closed form, and it is > 0.
* `c4even_pos`, `c4odd_pos`, `even_quartic_pos`, `odd_quadratic_disc`: sign facts.
* `table_values`      all 22 rationals of Table 1 (n = 4..25), exactly over ℚ;
                      `table_n7` (`9/14` at b₇ = 0), `table_n8` (`95/56` at b₈ = 3√7−8).
§6 Root-of-unity sums (all n)
* `fejer`             paper eq. (fejer): `Σ_{l=1}^{n−1}(1−cos(2πml/n))/(1−cos(2πl/n)) = m(n−m)`, 0 ≤ m ≤ n.
* `sum_cos_eq_zero`   `Σ_{k<n} cos(2πjk/n) = 0` for `n ∤ j`.
* `sum_cos_sq_mul_cos` `Σ_k cos²(lφ_k)cos φ_k = n/4`, `l = (n−1)/2`, odd n ≥ 3 (used in the `J₁` constraint).

In this module the Taylor coefficients `h₂, h₄` (even n) and `f₂, E₀₀…E₀₂` (odd n) and the normalisation
factors `ε ↔ τ`, `ε ↔ A` (in `c4evenOf`, `c4oddSet`) are DEFINITIONS copied from main.tex; `Conj13Series.lean`
derives all of them from the Hamiltonian, including the absence of an A³ term along the odd ansatz.
`Conj13Stability.lean` has sec:nonlinear (Lemmas B1–B4, Theorem 17) in an abstract setting.

NOT FORMALISED anywhere in the library (on paper / sympy scripts only)
* Prop. blocks / Prop. spectral beyond the integer and 2×2 algebra above: the Hessian of F at the ring
  (pair-log second derivatives and block-diagonalisation), `a_m = C₁ − m(n−m)` and the eigenvalues `b_m`
  as Hessian eigenvalues, the kernel being the radial mode V_{n/2}, and positivity on E. Here `a₁, b₁`
  are inputs to `V1'_n4`, `V1'_n5`; only the one-body part (`phi_crit_and_C1`) is derived.
* The Fourier sums `Σ e^{iφ}cos φ = n/2` etc. used to write `dJ₊` on `V₁`, and the derivative computation
  of `{J₁, J₂}` (`bracket_value` is the algebra only).
* Symmetric criticality (even n), completeness of the odd ansatz, and the identification of the slaved
  minimum with the reduced Hamiltonian's quartic (Def. c4).
* For sec:nonlinear: Lemma A (the slice), the manifold chart on S/K, existence of the slaving map β and the
  slaving inequality (hypotheses `hβ`, `hslave` of `lemma_B1`), and the hand-off from B3 to B4, which goes
  through the chart (a hypothesis); the verification of (A1)–(A4) for the ring.
* The dynamics (existence and global existence of the flow, conservation of H and J), the reduction
  `J⁻¹(μ)/K` as a manifold, and the numerical integrations.
-/
import Mathlib.Algebra.Order.Ring.Star
import Mathlib.Algebra.Order.Star.Real
import Mathlib.Algebra.Ring.IsFormallyReal
import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Analysis.LocallyConvex.AbsConvexOpen
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Data.Int.Star
import Mathlib.RingTheory.RootsOfUnity.Complex
import Mathlib.Tactic.NormNum.RealSqrt
import Mathlib.Topology.Algebra.Module.ModuleTopology
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.IntervalCases
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

namespace Havelock.Conj13

open Real

/-! ## §1. The bifurcation point `b_n` (MT (3.6) = paper eq. (2)) -/

/-- MT's threshold function `r(b) = (1+b²)/(1+b)²`. -/
noncomputable def rfun (b : ℝ) : ℝ := (1 + b ^ 2) / (1 + b) ^ 2

/-- The key factorisation behind uniqueness of `b_n`. -/
theorem rfun_cross (x y : ℝ) :
    (1 + x ^ 2) * (1 + y) ^ 2 - (1 + y ^ 2) * (1 + x) ^ 2 = 2 * (y - x) * (1 - x * y) := by ring

/-- `r` is injective on `(-1, 1]`: MT (3.6) has at most one root there. -/
theorem rfun_injOn (x y : ℝ) (hx : -1 < x) (hx1 : x ≤ 1) (hy : -1 < y) (hy1 : y ≤ 1)
    (h : rfun x = rfun y) : x = y := by
  unfold rfun at h
  have hx' : (1 + x) ^ 2 ≠ 0 := pow_ne_zero 2 (by linarith)
  have hy' : (1 + y) ^ 2 ≠ 0 := pow_ne_zero 2 (by linarith)
  rw [div_eq_div_iff hx' hy'] at h
  have h2 : 2 * (y - x) * (1 - x * y) = 0 := by rw [← rfun_cross]; linarith
  rcases mul_eq_zero.mp h2 with h3 | h3
  · linarith
  · -- x*y = 1 with x, y ∈ (-1,1]: (1-x)(1+y) ≥ 0 and (1-y)(1+x) ≥ 0 give y ≥ x ≥ y
    nlinarith [mul_nonneg (sub_nonneg.2 hx1) (by linarith : (0:ℝ) ≤ 1 + y),
      mul_nonneg (sub_nonneg.2 hy1) (by linarith : (0:ℝ) ≤ 1 + x)]

/-- `r(p/q) = (q²+p²)/(q+p)²`. -/
theorem rfun_div (p q : ℝ) (hq : q ≠ 0) (hpq : q + p ≠ 0) :
    rfun (p / q) = (q ^ 2 + p ^ 2) / (q + p) ^ 2 := by
  unfold rfun
  have h1 : 1 + p / q = (q + p) / q := by field_simp
  have h2 : 1 + (p / q) ^ 2 = (q ^ 2 + p ^ 2) / q ^ 2 := by field_simp
  rw [h1, h2, div_pow, div_div_div_cancel_right₀ (pow_ne_zero 2 hq)]

/-- `b_n` in closed form. Even `n`: with `t = √(n−1)`, `b_n = (1+2t−t²)/(t²+2t−1)`
(equal to the paper's `(4(n−2)√(n−1) − n²)/(n² − 8n + 8)`, see `bn_even_paper`).
Odd `n`: with `s = √(n−3)`, `b_n = (2−s)/(2+s)` (equal to the paper's
`(4√(n−3) − (n+1))/(n−7)` for `n ≠ 7`, and to `0` at `n = 7`, see `bn_odd_paper`). -/
noncomputable def bn (n : ℕ) : ℝ :=
  if n % 2 = 0 then
    (1 + 2 * √((n : ℝ) - 1) - (√((n : ℝ) - 1)) ^ 2) / ((√((n : ℝ) - 1)) ^ 2 + 2 * √((n : ℝ) - 1) - 1)
  else (2 - √((n : ℝ) - 3)) / (2 + √((n : ℝ) - 3))

/-- `⌊n²/4⌋ / (2(n−1))`, the right-hand side of MT (3.6). -/
noncomputable def rhs36 (n : ℕ) : ℝ := ((n ^ 2 / 4 : ℕ) : ℝ) / (2 * ((n : ℝ) - 1))

theorem floor_sq_div_four_even (M : ℕ) : (2 * M) ^ 2 / 4 = M ^ 2 := by
  have : (2 * M) ^ 2 = 4 * M ^ 2 := by ring
  rw [this]; omega

theorem floor_sq_div_four_odd (M : ℕ) : (2 * M + 1) ^ 2 / 4 = M ^ 2 + M := by
  have : (2 * M + 1) ^ 2 = 4 * (M ^ 2 + M) + 1 := by ring
  rw [this]; omega

/-- Even-`n` form of `b_n`, as an explicit function of `t = √(n−1)`. -/
theorem bn_even_facts (n : ℕ) (hn : 4 ≤ n) (he : n % 2 = 0) :
    let t := √((n : ℝ) - 1)
    bn n = (1 + 2 * t - t ^ 2) / (t ^ 2 + 2 * t - 1) ∧ t ^ 2 = (n : ℝ) - 1 ∧ 3 ≤ t ^ 2 ∧
      0 < t ^ 2 + 2 * t - 1 := by
  intro t
  have hn' : (4 : ℝ) ≤ n := by exact_mod_cast hn
  have ht2 : t ^ 2 = (n : ℝ) - 1 := Real.sq_sqrt (by linarith)
  have ht0 : 0 ≤ t := Real.sqrt_nonneg _
  refine ⟨by simp [bn, he, t], ht2, by linarith, by nlinarith⟩

/-- Odd-`n` form of `b_n`, as an explicit function of `s = √(n−3)`. -/
theorem bn_odd_facts (n : ℕ) (hn : 4 ≤ n) (ho : n % 2 = 1) :
    let s := √((n : ℝ) - 3)
    bn n = (2 - s) / (2 + s) ∧ s ^ 2 = (n : ℝ) - 3 ∧ 0 < s := by
  intro s
  have hn' : (4 : ℝ) ≤ n := by exact_mod_cast hn
  refine ⟨by simp [bn, ho, s], Real.sq_sqrt (by linarith), Real.sqrt_pos.mpr (by linarith)⟩

/-- **MT (3.6), existence + location**: for every `n ≥ 4`, `b_n ∈ (−1, 1)` and
`(1+b_n²)/(1+b_n)² = ⌊n²/4⌋/(2(n−1))`. -/
theorem bn_spec (n : ℕ) (hn : 4 ≤ n) :
    -1 < bn n ∧ bn n < 1 ∧ rfun (bn n) = rhs36 n := by
  have hn' : (4 : ℝ) ≤ n := by exact_mod_cast hn
  rcases Nat.mod_two_eq_zero_or_one n with he | ho
  · obtain ⟨hb, ht2, ht3, hden⟩ := bn_even_facts n hn he
    set t := √((n : ℝ) - 1)
    have ht0 : 0 ≤ t := Real.sqrt_nonneg _
    have ht1 : 1 < t := by nlinarith
    obtain ⟨M, rfl⟩ : ∃ M, n = 2 * M := ⟨n / 2, by omega⟩
    have hM : ((2 * M : ℕ) : ℝ) = 2 * (M : ℝ) := by push_cast; ring
    rw [hb]
    refine ⟨?_, ?_, ?_⟩
    · rw [lt_div_iff₀ hden]; nlinarith
    · rw [div_lt_iff₀ hden]; nlinarith
    · rw [rfun_div _ _ hden.ne' (by nlinarith)]
      unfold rhs36
      rw [floor_sq_div_four_even]
      have hMt : (2 * (M : ℝ)) = t ^ 2 + 1 := by rw [← hM, ht2]; ring
      have hM' : (M : ℝ) = (t ^ 2 + 1) / 2 := by linarith
      rw [div_eq_div_iff (by nlinarith) (by rw [hM]; nlinarith)]
      push_cast
      rw [hM']
      ring
  · obtain ⟨hb, hs2, hs0⟩ := bn_odd_facts n hn ho
    set s := √((n : ℝ) - 3)
    obtain ⟨M, rfl⟩ : ∃ M, n = 2 * M + 1 := ⟨n / 2, by omega⟩
    rw [hb]
    have hden : 0 < 2 + s := by linarith
    refine ⟨?_, ?_, ?_⟩
    · rw [lt_div_iff₀ hden]; nlinarith
    · rw [div_lt_iff₀ hden]; nlinarith
    · rw [rfun_div _ _ hden.ne' (by nlinarith)]
      unfold rhs36
      rw [floor_sq_div_four_odd]
      push_cast at hs2 ⊢
      have hM' : (M : ℝ) = (s ^ 2) / 2 + 1 := by linarith
      rw [div_eq_div_iff (by nlinarith) (by rw [hM']; nlinarith)]
      rw [hM']
      ring

/-- **MT (3.6), uniqueness**: `b_n` is the unique root of (3.6) in `(−1, 1]`. -/
theorem bn_unique (n : ℕ) (hn : 4 ≤ n) (x : ℝ) (hx : -1 < x) (hx1 : x ≤ 1)
    (h : rfun x = rhs36 n) : x = bn n := by
  obtain ⟨h1, h2, h3⟩ := bn_spec n hn
  exact rfun_injOn x (bn n) hx hx1 h1 h2.le (h.trans h3.symm)

/-- The paper's even closed form `b_n = (4(n−2)√(n−1) − n²)/(n² − 8n + 8)`. -/
theorem bn_even_paper (n : ℕ) (hn : 4 ≤ n) (he : n % 2 = 0) :
    bn n = (4 * ((n : ℝ) - 2) * √((n : ℝ) - 1) - (n : ℝ) ^ 2) / ((n : ℝ) ^ 2 - 8 * n + 8) := by
  obtain ⟨hb, ht2, ht3, hden⟩ := bn_even_facts n hn he
  set t := √((n : ℝ) - 1)
  have ht0 : 0 ≤ t := Real.sqrt_nonneg _
  have hN : (n : ℝ) = t ^ 2 + 1 := by linarith
  -- n² − 8n + 8 has no integer root
  have hq : ((n : ℤ) ^ 2 - 8 * n + 8) ≠ 0 := by
    intro h
    have : (n : ℤ) ≤ 8 := by nlinarith
    have : n ≤ 8 := by omega
    interval_cases n <;> simp_all
  have hq' : ((n : ℝ) ^ 2 - 8 * n + 8) ≠ 0 := by exact_mod_cast hq
  -- (t²−2t−1)(t²+2t−1) = n² − 8n + 8
  have hfac : (t ^ 2 - 2 * t - 1) * (t ^ 2 + 2 * t - 1) = (n : ℝ) ^ 2 - 8 * n + 8 := by
    rw [hN]; ring
  have hm : t ^ 2 - 2 * t - 1 ≠ 0 := by
    intro h0; apply hq'; rw [← hfac, h0, zero_mul]
  rw [hb, div_eq_div_iff hden.ne' hq', ← hfac, hN]
  ring

/-- The paper's odd closed form `b_n = (4√(n−3) − (n+1))/(n−7)`, `n ≠ 7`. -/
theorem bn_odd_paper (n : ℕ) (hn : 4 ≤ n) (ho : n % 2 = 1) (h7 : n ≠ 7) :
    bn n = (4 * √((n : ℝ) - 3) - ((n : ℝ) + 1)) / ((n : ℝ) - 7) := by
  obtain ⟨hb, hs2, hs0⟩ := bn_odd_facts n hn ho
  set s := √((n : ℝ) - 3)
  have hN : (n : ℝ) = s ^ 2 + 3 := by linarith
  have h7' : (n : ℝ) - 7 ≠ 0 := by
    intro h; apply h7; have : (n : ℝ) = 7 := by linarith
    exact_mod_cast this
  rw [hb, div_eq_div_iff (by linarith) h7', hN]
  ring

theorem bn_seven : bn 7 = 0 := by
  obtain ⟨hb, hs2, hs0⟩ := bn_odd_facts 7 (by norm_num) (by norm_num)
  rw [hb]
  have : √((7 : ℕ) - 3 : ℝ) = 2 := by
    rw [show ((7 : ℕ) : ℝ) - 3 = 2 ^ 2 by norm_num]; exact Real.sqrt_sq (by norm_num)
  rw [this]; norm_num

theorem bn_four : bn 4 = 2 - √3 := by
  rw [bn_even_paper 4 le_rfl rfl]
  norm_num
  have h3 : √3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  rw [div_eq_iff (by norm_num)]; ring

theorem bn_five : bn 5 = 3 - 2 * √2 := by
  obtain ⟨hb, hs2, hs0⟩ := bn_odd_facts 5 (by norm_num) (by norm_num)
  rw [hb]
  have : √((5 : ℕ) - 3 : ℝ) = √2 := by norm_num
  rw [this]
  have h2 : √2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have h0 : 0 < √2 := by positivity
  rw [div_eq_iff (by linarith)]; nlinarith

theorem bn_six : bn 6 = 9 - 4 * √5 := by
  rw [bn_even_paper 6 (by norm_num) rfl]
  norm_num
  rw [div_eq_iff (by norm_num)]; ring

theorem bn_eight : bn 8 = -(8 - 3 * √7) := by
  rw [bn_even_paper 8 (by norm_num) rfl]
  norm_num
  rw [div_eq_iff (by norm_num)]; ring

theorem bn_ten : bn 10 = -1 / 7 := by
  rw [bn_even_paper 10 (by norm_num) rfl]
  have : √((10 : ℕ) - 1 : ℝ) = 3 := by
    rw [show ((10 : ℕ) : ℝ) - 1 = 3 ^ 2 by norm_num]; exact Real.sqrt_sq (by norm_num)
  rw [this]; norm_num

/-- Geometry of the bifurcation point: sphere (`b_n > 0`) for `n = 4, 5, 6`,
plane (`b_7 = 0`), hyperbolic plane (`b_n < 0`) for every `n ≥ 8`. -/
theorem bn_sign (n : ℕ) (hn : 4 ≤ n) :
    (n ≤ 6 → 0 < bn n) ∧ (8 ≤ n → bn n < 0) := by
  have hn' : (4 : ℝ) ≤ n := by exact_mod_cast hn
  rcases Nat.mod_two_eq_zero_or_one n with he | ho
  · obtain ⟨hb, ht2, ht3, hden⟩ := bn_even_facts n hn he
    set t := √((n : ℝ) - 1)
    have ht0 : 0 ≤ t := Real.sqrt_nonneg _
    rw [hb]
    refine ⟨fun h6 => ?_, fun h8 => ?_⟩
    · have : (n : ℝ) ≤ 6 := by exact_mod_cast h6
      apply div_pos _ hden; nlinarith
    · have : (8 : ℝ) ≤ n := by exact_mod_cast h8
      apply div_neg_of_neg_of_pos _ hden; nlinarith
  · obtain ⟨hb, hs2, hs0⟩ := bn_odd_facts n hn ho
    set s := √((n : ℝ) - 3)
    rw [hb]
    refine ⟨fun h6 => ?_, fun h8 => ?_⟩
    · have : (n : ℝ) ≤ 6 := by exact_mod_cast h6
      apply div_pos _ (by linarith); nlinarith
    · have : (8 : ℝ) ≤ n := by exact_mod_cast h8
      apply div_neg_of_neg_of_pos _ (by linarith); nlinarith


/-! ## §2. Fourier blocks, reduced-space algebra (paper §3, Prop. blocks / Prop. spectral) -/

/-- Paper §3.1: the one-body critical point and second derivative.
With `φ'(0) = −(n−1) + 2(n−1)b/(1+b) − 2ω₀/(1+b)²` and
`C₁ = φ''(0) = 4(n−1)b/(1+b)² − ω₀·4(1−b)/(1+b)³`, the choice `ω₀ = −(n−1)(1−b²)/2` makes
`φ'(0) = 0` and gives `C₁ = 2(n−1)(1+b²)/(1+b)²`. (The derivative formulas themselves are
`phi_hasDerivAt`, `phi'_hasDerivAt` below.) -/
theorem omega0_C1 (n b : ℝ) (hb : 1 + b ≠ 0) :
    let ω₀ := -(n - 1) * (1 - b ^ 2) / 2
    (-(n - 1) + 2 * (n - 1) * b / (1 + b) - ω₀ * 2 / (1 + b) ^ 2 = 0) ∧
    (4 * (n - 1) * b / (1 + b) ^ 2 - ω₀ * 4 * (1 - b) / (1 + b) ^ 3
        = 2 * (n - 1) * (1 + b ^ 2) / (1 + b) ^ 2) := by
  intro ω₀
  constructor <;> (simp only [ω₀]; field_simp; ring)

/-- The one-body part `φ(s) = −(n−1)s + (n−1) log(1 + b e^{2s}) − ω e^{2s}/(1 + b e^{2s})`. -/
noncomputable def phi (n b ω s : ℝ) : ℝ :=
  -(n - 1) * s + (n - 1) * Real.log (1 + b * Real.exp (2 * s)) -
    ω * (Real.exp (2 * s) / (1 + b * Real.exp (2 * s)))

/-- Its derivative, valid wherever `1 + b e^{2s} > 0`. -/
noncomputable def phi' (n b ω s : ℝ) : ℝ :=
  -(n - 1) + (n - 1) * (2 * b * Real.exp (2 * s)) / (1 + b * Real.exp (2 * s)) -
    ω * (2 * Real.exp (2 * s) / (1 + b * Real.exp (2 * s)) ^ 2)

theorem phi_hasDerivAt (n b ω s : ℝ) (hpos : 0 < 1 + b * Real.exp (2 * s)) :
    HasDerivAt (phi n b ω) (phi' n b ω s) s := by
  have he : HasDerivAt (fun s => Real.exp (2 * s)) (Real.exp (2 * s) * 2) s := by
    simpa using ((hasDerivAt_id s).const_mul 2).exp
  have hu : HasDerivAt (fun s => 1 + b * Real.exp (2 * s)) (b * (Real.exp (2 * s) * 2)) s :=
    (he.const_mul b).const_add 1
  have hlog := hu.log hpos.ne'
  have hdiv := he.div hu hpos.ne'
  have hlin : HasDerivAt (fun s => -(n - 1) * s) (-(n - 1)) s := by
    simpa using (hasDerivAt_id s).const_mul (-(n - 1))
  have := (hlin.add (hlog.const_mul (n - 1))).sub (hdiv.const_mul ω)
  convert this using 1
  · funext x; simp only [phi, Pi.sub_apply, Pi.add_apply, Pi.div_apply]
  · unfold phi'
    field_simp
    ring

theorem phi'_hasDerivAt_zero (n b ω : ℝ) (hpos : 0 < 1 + b) :
    HasDerivAt (phi' n b ω)
      (4 * (n - 1) * b / (1 + b) ^ 2 - ω * 4 * (1 - b) / (1 + b) ^ 3) 0 := by
  have he : HasDerivAt (fun s => Real.exp (2 * s)) (Real.exp (2 * 0) * 2) 0 := by
    simpa using ((hasDerivAt_id (0:ℝ)).const_mul 2).exp
  have hu : HasDerivAt (fun s => 1 + b * Real.exp (2 * s)) (b * (Real.exp (2 * 0) * 2)) 0 :=
    (he.const_mul b).const_add 1
  have h0 : (1 + b * Real.exp (2 * 0)) ≠ 0 := by simp; linarith
  have h1 := ((he.const_mul (2 * b)).div hu h0).const_mul (n - 1)
  have h2 := ((he.const_mul 2).div (hu.pow 2) (pow_ne_zero 2 h0)).const_mul ω
  have := (h1.const_add (-(n - 1))).sub h2
  convert this using 1
  · funext x; simp only [phi', Pi.sub_apply, Pi.div_apply, Pi.pow_apply]; ring
  · simp
    have : 1 + b ≠ 0 := by linarith
    field_simp
    ring

/-- **Paper §3.1 (one-body part of `F = H − ω₀J₃`)**: with `ω₀ = −(n−1)(1−b²)/2`, the
one-body function `φ` has a critical point at `s = 0` and `φ''(0) = C₁ = 2(n−1)(1+b²)/(1+b)²`. -/
theorem phi_crit_and_C1 (n b : ℝ) (hb : -1 < b) :
    let ω₀ := -(n - 1) * (1 - b ^ 2) / 2
    HasDerivAt (phi n b ω₀) 0 0 ∧
    HasDerivAt (phi' n b ω₀) (2 * (n - 1) * (1 + b ^ 2) / (1 + b) ^ 2) 0 ∧
    (∀ᶠ s in nhds (0:ℝ), HasDerivAt (phi n b ω₀) (phi' n b ω₀ s) s) := by
  intro ω₀
  have hpos : 0 < 1 + b := by linarith
  obtain ⟨hc, hC⟩ := omega0_C1 n b hpos.ne'
  refine ⟨?_, ?_, ?_⟩
  · have := phi_hasDerivAt n b ω₀ 0 (by simp; linarith)
    convert this using 1
    unfold phi'; simp; rw [← hc]; ring
  · rw [← hC]; exact phi'_hasDerivAt_zero n b ω₀ hpos
  · have hcont : Continuous (fun s : ℝ => 1 + b * Real.exp (2 * s)) := by fun_prop
    have h0 : (0:ℝ) < 1 + b * Real.exp (2 * 0) := by simp; linarith
    have : ∀ᶠ s in nhds (0:ℝ), 0 < 1 + b * Real.exp (2 * s) :=
      continuousAt_const.eventually_lt hcont.continuousAt h0
    filter_upwards [this] with s hs
    exact phi_hasDerivAt n b ω₀ s hs

/-- **Threshold ⇔ MT (3.6)** (Prop. blocks): `C₁ = ⌊n²/4⌋` iff `(1+b²)/(1+b)² = ⌊n²/4⌋/(2(n−1))`. -/
theorem C1_eq_iff (n : ℕ) (hn : 2 ≤ n) (b : ℝ) (hb : 1 + b ≠ 0) :
    2 * ((n : ℝ) - 1) * (1 + b ^ 2) / (1 + b) ^ 2 = ((n ^ 2 / 4 : ℕ) : ℝ) ↔ rfun b = rhs36 n := by
  have hn1 : (n : ℝ) - 1 ≠ 0 := by
    have : (2 : ℝ) ≤ n := by exact_mod_cast hn
    linarith
  unfold rfun rhs36
  constructor
  · intro h; rw [← h]; field_simp
  · intro h; rw [mul_div_assoc, h]; field_simp

/-- **Radial Fourier eigenvalues** `a_m = ⌊n²/4⌋ − m(n−m)` at `b_n` (Prop. blocks / Prop. spectral):
for `m ≤ ⌊n/2⌋`, `a_m = (⌊n/2⌋ − m)(⌊n/2⌋ − m + (n mod 2))`. Hence `a_m ≥ 0`, `a_m = 0` iff
`m = ⌊n/2⌋`, and the next value `a_{⌊n/2⌋−1}` is `1` (even `n`) or `2` (odd `n`). -/
theorem am_factor (n m : ℕ) (hm : m ≤ n / 2) :
    ((n ^ 2 / 4 : ℕ) : ℤ) - (m : ℤ) * ((n : ℤ) - m)
      = ((n / 2 : ℕ) - (m : ℤ)) * ((n / 2 : ℕ) - (m : ℤ) + (n % 2 : ℕ)) := by
  rcases Nat.mod_two_eq_zero_or_one n with he | ho
  · obtain ⟨M, rfl⟩ : ∃ M, n = 2 * M := ⟨n / 2, by omega⟩
    rw [floor_sq_div_four_even, show 2 * M / 2 = M by omega, show 2 * M % 2 = 0 by omega]
    push_cast; ring
  · obtain ⟨M, rfl⟩ : ∃ M, n = 2 * M + 1 := ⟨n / 2, by omega⟩
    rw [floor_sq_div_four_odd, show (2 * M + 1) / 2 = M by omega,
      show (2 * M + 1) % 2 = 1 by omega]
    push_cast; ring

theorem am_nonneg_and_kernel (n m : ℕ) (hm : m ≤ n / 2) :
    0 ≤ ((n ^ 2 / 4 : ℕ) : ℤ) - (m : ℤ) * ((n : ℤ) - m) ∧
    (((n ^ 2 / 4 : ℕ) : ℤ) - (m : ℤ) * ((n : ℤ) - m) = 0 ↔ m = n / 2) := by
  rw [am_factor n m hm]
  have h1 : (0 : ℤ) ≤ ((n / 2 : ℕ) : ℤ) - m := by omega
  have h2 : (0 : ℤ) ≤ ((n % 2 : ℕ) : ℤ) := by positivity
  refine ⟨by nlinarith, ?_⟩
  constructor
  · intro h
    rcases mul_eq_zero.mp h with h | h <;> omega
  · intro h; subst h; simp

theorem am_next (n : ℕ) (hn : 4 ≤ n) :
    ((n ^ 2 / 4 : ℕ) : ℤ) - ((n / 2 - 1 : ℕ) : ℤ) * ((n : ℤ) - ((n / 2 - 1 : ℕ) : ℤ))
      = 1 + (n % 2 : ℕ) := by
  rw [am_factor n (n / 2 - 1) (by omega)]
  have : ((n / 2 - 1 : ℕ) : ℤ) = ((n / 2 : ℕ) : ℤ) - 1 := by omega
  rw [this]; ring

/-- **Mode-1 block restricted to `V₁' = V₁ ∩ ker dJ₊`** (Prop. spectral, "mode 1"; paper remark after
Prop. spectral). In the coordinates `s = x cos φ + y sin φ`, `θ = u cos φ + v sin φ` on `V₁`,
`dJ₊ = 0` reads `A x − B v = 0`, `A y + B u = 0` (this coordinate reduction of `dJ₊` uses
`Σ_k e^{iφ_k} cos φ_k = n/2`, `Σ_k e^{iφ_k} sin φ_k = i n/2`, see `sum_cexp_mul_cos`),
and on it the Hessian `a₁(x²+y²) + b₁(u²+v²)` equals `κ (x²+y²+u²+v²)` with
`κ = (a₁ B² + b₁ A²)/(A² + B²)`. -/
theorem V1'_form (a₁ b₁ A B x y u v : ℝ) (hAB : A ^ 2 + B ^ 2 ≠ 0)
    (h1 : A * x - B * v = 0) (h2 : A * y + B * u = 0) :
    a₁ * (x ^ 2 + y ^ 2) + b₁ * (u ^ 2 + v ^ 2)
      = (a₁ * B ^ 2 + b₁ * A ^ 2) / (A ^ 2 + B ^ 2) * (x ^ 2 + y ^ 2 + u ^ 2 + v ^ 2) := by
  rw [div_mul_eq_mul_div, eq_div_iff hAB]
  linear_combination (a₁ - b₁) * (A * x + B * v) * h1 + (a₁ - b₁) * (A * y - B * u) * h2

/-- With `A = (1−b)/(1+b)²`, `B = 1/(1+b)` (paper §3.2(b)), `κ = (a₁ + b₁ρ²)/(1+ρ²)`,
`ρ = (1−b)/(1+b)`. -/
theorem V1'_kappa_rho (a₁ b₁ b : ℝ) (hb : 1 + b ≠ 0) :
    let A := (1 - b) / (1 + b) ^ 2
    let B := 1 / (1 + b)
    let ρ := (1 - b) / (1 + b)
    (a₁ * B ^ 2 + b₁ * A ^ 2) / (A ^ 2 + B ^ 2) = (a₁ + b₁ * ρ ^ 2) / (1 + ρ ^ 2) := by
  intro A B ρ
  simp only [A, B, ρ]
  have h2 : (1 + b) ^ 2 + (1 - b) ^ 2 ≠ 0 := by
    intro h; have : b ^ 2 = -1 := by nlinarith
    nlinarith [sq_nonneg b]
  field_simp
  ring

/-- **Paper, after Prop. spectral**: for `n = 4` (`a₁ = 1, b₁ = 3`, `b₄ = 2 − √3`) the `V₁'`
eigenvalue is `3/2`. -/
theorem V1'_n4 :
    let b := bn 4
    let ρ := (1 - b) / (1 + b)
    ((((4:ℕ) ^ 2 / 4 : ℕ) : ℝ) - (4 - 1) + (4 - 1) * ρ ^ 2) / (1 + ρ ^ 2) = 3 / 2 := by
  intro b ρ
  have h3 : √3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  have h3p : 1 < √3 := by
    rw [show (1:ℝ) = √1 by simp]; exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
  have hρ : ρ ^ 2 = 1 / 3 := by
    simp only [ρ, b, bn_four]
    rw [div_pow, div_eq_div_iff (by nlinarith) (by norm_num)]
    nlinarith
  rw [hρ]; norm_num

/-- **Paper, after Prop. spectral**: for `n = 5` (`a₁ = 2, b₁ = 4`, `b₅ = 3 − 2√2`) the `V₁'`
eigenvalue is `8/3`. -/
theorem V1'_n5 :
    let b := bn 5
    let ρ := (1 - b) / (1 + b)
    ((((5:ℕ) ^ 2 / 4 : ℕ) : ℝ) - (5 - 1) + (5 - 1) * ρ ^ 2) / (1 + ρ ^ 2) = 8 / 3 := by
  intro b ρ
  have h2 : √2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have h2p : √2 < 2 := by
    rw [show (2:ℝ) = √4 by rw [show (4:ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]]
    exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
  have hρ : ρ ^ 2 = 1 / 2 := by
    simp only [ρ, b, bn_five]
    rw [div_pow, div_eq_div_iff (by nlinarith) (by norm_num)]
    nlinarith
  rw [hρ]; norm_num

/-- **Paper eq. (bracket)**: `{J₁,J₂}(x₀) = n/2 − λJ₃(x₀)` with `J₃(x₀) = n/(1+b)` equals
`n(1−b)/(2(1+b))`, which is nonzero for `b ∈ (−1, 1)`. -/
theorem bracket_value (n b : ℝ) (hn : n ≠ 0) (hb : -1 < b) (hb1 : b < 1) :
    n / 2 - b * (n / (1 + b)) = n * (1 - b) / (2 * (1 + b)) ∧ n * (1 - b) / (2 * (1 + b)) ≠ 0 := by
  have h : 1 + b ≠ 0 := by linarith
  refine ⟨by field_simp; ring, ?_⟩
  apply div_ne_zero (mul_ne_zero hn (by linarith)) (by positivity)

/-! ## §3. Even `n = 2M`: the reduced function on `Fix(D_M) ∩ J⁻¹(μ)` (Thm even) -/

/-- Paper Thm even proof: the `τ²` coefficient `h₂`. -/
noncomputable def h2 (M b : ℝ) : ℝ :=
  -(M / 4) * (1 + b) ^ 2 * (M ^ 2 * (1 + b) ^ 2 - (4 * M - 2) * (1 + b ^ 2))

/-- Paper Thm even proof: the `τ⁴` coefficient `h₄`. -/
noncomputable def h4 (M b : ℝ) : ℝ :=
  M * (1 + b) ^ 4 / 96 *
    (M ^ 4 * (1 + b) ^ 4 - 16 * M ^ 2 * (1 + b) * (1 + b ^ 3) + 24 * (2 * M - 1) * (1 + b ^ 4))

/-- Paper Thm even proof: `c₄ = h₄ (2/(1+λ)²)⁴ / n²`, `n = 2M` (unit kernel vector normalisation). -/
noncomputable def c4evenOf (M b : ℝ) : ℝ := h4 M b * (2 / (1 + b) ^ 2) ^ 4 / (2 * M) ^ 2

/-- The even closed form `(n−2)(m⁴+40m²+48m+16)/(192 n (n−1))`, `m = n − 2`. -/
def c4even {K : Type*} [Field K] (n : K) : K :=
  (n - 2) * ((n - 2) ^ 4 + 40 * (n - 2) ^ 2 + 48 * (n - 2) + 16) / (192 * n * (n - 1))

/-- The even closed form of paper Thm even (`conj13/paper/main.tex`), expanded as a polynomial in `n`. -/
theorem c4even_alt {K : Type*} [Field K] (n : K) :
    c4even n = (n - 2) * (n ^ 4 - 8 * n ^ 3 + 64 * n ^ 2 - 144 * n + 96) / (192 * n * (n - 1)) := by
  unfold c4even; ring

/-- `h₂ = 0` iff MT (3.6) for even `n = 2M`, i.e. `2(2M−1)(1+b²) = M²(1+b)²`. -/
theorem h2_eq_zero_iff (M b : ℝ) (hM : M ≠ 0) (hb : 1 + b ≠ 0) :
    h2 M b = 0 ↔ 2 * (2 * M - 1) * (1 + b ^ 2) = M ^ 2 * (1 + b) ^ 2 := by
  unfold h2
  constructor
  · intro h
    have : M ^ 2 * (1 + b) ^ 2 - (4 * M - 2) * (1 + b ^ 2) = 0 := by
      rcases mul_eq_zero.mp h with h | h
      · rcases mul_eq_zero.mp h with h | h
        · exfalso; apply hM; linarith
        · exact absurd (pow_eq_zero_iff (n := 2) (by norm_num) |>.mp h) hb
      · exact h
    linarith
  · intro h
    have : M ^ 2 * (1 + b) ^ 2 - (4 * M - 2) * (1 + b ^ 2) = 0 := by linarith
    rw [this, mul_zero]

/-- The exact difference identity behind Thm even (no threshold assumed):
`c4evenOf − c4even(2M) = −R·S / (8M(2M−1)(1+b)⁴)` with `R` the threshold polynomial. -/
theorem c4even_diff (M b : ℝ) (hM : M ≠ 0) (hM1 : 2 * M - 1 ≠ 0) (hb : 1 + b ≠ 0) :
    c4evenOf M b - c4even (2 * M)
      = (M ^ 2 * (1 + b) ^ 2 - 2 * (2 * M - 1) * (1 + b ^ 2)) *
          (3 * M ^ 2 * (1 + b) ^ 2 - 12 * M * b ^ 2 - 16 * M * b - 12 * M + 6 * b ^ 2 + 8 * b + 6) /
          (8 * M * (2 * M - 1) * (1 + b) ^ 4) := by
  unfold c4evenOf h4 c4even
  have h2M : 2 * M ≠ 0 := mul_ne_zero two_ne_zero hM
  have hM1' : M * 2 - 1 ≠ 0 := by rwa [mul_comm]
  field_simp
  ring

/-- **Thm even (even `n = 2M`), reduction step**: on the threshold
`2(2M−1)(1+b²) = M²(1+b)²` (MT (3.6) for even `n`), `h₂ = 0` and
`c₄ = h₄(2/(1+b)²)⁴/n² = (n−2)(m⁴+40m²+48m+16)/(192 n(n−1))`, `n = 2M`, `m = n−2`. -/
theorem c4even_reduction (M b : ℝ) (hM : M ≠ 0) (hM1 : 2 * M - 1 ≠ 0) (hb : 1 + b ≠ 0)
    (hrel : 2 * (2 * M - 1) * (1 + b ^ 2) = M ^ 2 * (1 + b) ^ 2) :
    h2 M b = 0 ∧ c4evenOf M b = c4even (2 * M) := by
  refine ⟨(h2_eq_zero_iff M b hM hb).mpr hrel, ?_⟩
  have := c4even_diff M b hM hM1 hb
  rw [show M ^ 2 * (1 + b) ^ 2 - 2 * (2 * M - 1) * (1 + b ^ 2) = 0 by linarith] at this
  simpa [sub_eq_zero] using this


/-! ## §4. Odd `n`: the `A⁴` coefficient `E(p,s)` and the `J₁` constraint (Thm odd) -/

/-- Paper Thm odd proof: the kernel `A²` coefficient `f₂`. -/
noncomputable def f2 (n b : ℝ) : ℝ :=
  -(n * (n - 1) * (1 + b) ^ 2 / 64) * ((n - 7) * (1 + b ^ 2) + 2 * (n + 1) * b)

/-- Paper Thm odd proof: the coefficients `E₀₀, E₁₀, E₀₁, E₂₀, E₀₂` (as stated in the paper). -/
noncomputable def E00 (n b : ℝ) : ℝ :=
  n * (n - 1) * (1 + b) ^ 4 / 8192 *
    ((n ^ 3 + n ^ 2 - 69 * n + 291) * (1 + b ^ 4) + 4 * (n ^ 3 + n ^ 2 - 21 * n - 13) * (b + b ^ 3)
      + (6 * n ^ 3 + 6 * n ^ 2 - 30 * n + 82) * b ^ 2)
noncomputable def E10 (n b : ℝ) : ℝ :=
  n * (n - 1) * (1 - b) * (1 + b) ^ 3 / 64 * ((n - 5) * (1 + b ^ 2) + 2 * (n - 1) * b)
noncomputable def E01 (n b : ℝ) : ℝ := n * (n - 1) ^ 2 * (1 + b) ^ 4 / 64
noncomputable def E20 (n b : ℝ) : ℝ := n * (n - 1) * (1 - b) ^ 2 * (1 + b) ^ 2 / 16
noncomputable def E02 (n : ℝ) : ℝ := n * (n - 1) / 4

/-- `E(p,s) = E₀₀ + E₁₀p + E₀₁s + E₂₀p² + E₀₂s²` (no `ps` term). -/
noncomputable def Efun (n b p s : ℝ) : ℝ :=
  E00 n b + E10 n b * p + E01 n b * s + E20 n b * p ^ 2 + E02 n * s ^ 2

/-- The `O(A²)` part of the constraint `J₁ = 0`, divided by `(n/2)A²`:
`Q'(w₀)p − Q(w₀)s + Q''(w₀)/4` with `Q(w₀) = 1/(1+b)`, `Q'(w₀) = (1−b)/2`,
`Q''(w₀) = −(1+b)³/4` (these three values are `Q_values`). -/
noncomputable def J1con (b p s : ℝ) : ℝ := (1 - b) / 2 * p - 1 / (1 + b) * s + (-(1 + b) ^ 3 / 4) / 4

/-- The slaved `s` as a function of `p` (solution of `J1con = 0`). -/
noncomputable def sOf (b p : ℝ) : ℝ := (1 + b) * ((1 - b) * p / 2 - (1 + b) ^ 3 / 16)

/-- `a₂`, the coefficient of `p²` after eliminating `s` (paper: `n(n−1)(1−λ)²(1+λ)²/8`). -/
noncomputable def a2 (n b : ℝ) : ℝ := n * (n - 1) * (1 - b) ^ 2 * (1 + b) ^ 2 / 8

/-- The minimising `p`. -/
noncomputable def pstar (n b : ℝ) : ℝ :=
  (1 + b) * ((3 * n - 13) * (1 + b ^ 2) + (6 * n - 10) * b) / (32 * (b - 1))

/-- The odd closed form `(n−1)(n−3)(n²−9n+38)/(128 n)`. -/
def c4odd {K : Type*} [Field K] (n : K) : K := (n - 1) * (n - 3) * (n ^ 2 - 9 * n + 38) / (128 * n)

/-- The set of values of `64/(n²(1+b)⁸)·E(p,s)` over the constraint `J₁ = 0`; its least element
is the level-set quartic `c₄` (paper: `c₄ = 64/(n²(1+λ)⁸) min_p E(p, s(p))`). -/
noncomputable def c4oddSet (n b : ℝ) : Set ℝ :=
  {c | ∃ p s, J1con b p s = 0 ∧ c = 64 / (n ^ 2 * (1 + b) ^ 8) * Efun n b p s}

/-- Odd threshold relation ⇔ MT (3.6) for odd `n`: `(1+b²)/(1+b)² = (n+1)/8`. -/
theorem odd_rel_iff (n b : ℝ) (hb : 1 + b ≠ 0) :
    (n - 7) * (1 + b ^ 2) + 2 * (n + 1) * b = 0 ↔ rfun b = (n + 1) / 8 := by
  unfold rfun
  rw [div_eq_div_iff (pow_ne_zero 2 hb) (by norm_num)]
  constructor <;> intro h <;> linarith

/-- `f₂ = 0` iff the odd threshold relation. -/
theorem f2_eq_zero_iff (n b : ℝ) (hn : n ≠ 0) (hn1 : n ≠ 1) (hb : 1 + b ≠ 0) :
    f2 n b = 0 ↔ (n - 7) * (1 + b ^ 2) + 2 * (n + 1) * b = 0 := by
  unfold f2
  have hc : n * (n - 1) * (1 + b) ^ 2 / 64 ≠ 0 :=
    div_ne_zero (mul_ne_zero (mul_ne_zero hn (sub_ne_zero.mpr hn1)) (pow_ne_zero 2 hb)) (by norm_num)
  rw [neg_mul, neg_eq_zero, mul_eq_zero]
  simp [hc]

theorem J1con_iff (b p s : ℝ) (hb : 1 + b ≠ 0) : J1con b p s = 0 ↔ s = sOf b p := by
  unfold J1con sOf
  constructor
  · intro h; field_simp at h ⊢; linarith
  · intro h; rw [h]; field_simp; ring

/-- Eliminating `s` leaves a quadratic in `p` with leading coefficient `a₂` and vertex `p*`. -/
theorem E_vertex (n b p : ℝ) (hb1 : b - 1 ≠ 0) :
    Efun n b p (sOf b p) = Efun n b (pstar n b) (sOf b (pstar n b)) + a2 n b * (p - pstar n b) ^ 2 := by
  unfold Efun E00 E10 E01 E20 E02 sOf a2 pstar
  field_simp
  ring

theorem a2_pos (n b : ℝ) (hn : 1 < n) (hb : -1 < b) (hb1 : b < 1) : 0 < a2 n b := by
  unfold a2
  have h1 : 0 < (1 - b) ^ 2 * (1 + b) ^ 2 := by
    apply mul_pos <;> apply pow_pos <;> linarith
  have h2 : 0 < n * (n - 1) := mul_pos (by linarith) (by linarith)
  have : n * (n - 1) * (1 - b) ^ 2 * (1 + b) ^ 2 = (n * (n - 1)) * ((1 - b) ^ 2 * (1 + b) ^ 2) := by ring
  rw [this]
  exact div_pos (mul_pos h2 h1) (by norm_num)

/-- The minimum value, reduced modulo the odd threshold relation. -/
theorem odd_min_value (n b : ℝ) (hb : -1 < b) (hb1 : b < 1)
    (hrel : (n - 7) * (1 + b ^ 2) + 2 * (n + 1) * b = 0) :
    64 / (n ^ 2 * (1 + b) ^ 8) * Efun n b (pstar n b) (sOf b (pstar n b)) = c4odd n := by
  have h1 : 1 + b ≠ 0 := by linarith
  have h2 : b - 1 ≠ 0 := by linarith
  have hn : n = (7 - 2 * b + 7 * b ^ 2) / (1 + b) ^ 2 := by
    rw [eq_div_iff (pow_ne_zero 2 h1)]; linarith
  have h7 : 7 - 2 * b + 7 * b ^ 2 ≠ 0 := by nlinarith [sq_nonneg b]
  subst hn
  unfold Efun E00 E10 E01 E20 E02 sOf pstar c4odd
  field_simp
  ring

/-- **Thm odd (odd `n`), reduction step**: on the odd threshold relation
`(n−7)(1+b²) + 2(n+1)b = 0` with `b ∈ (−1,1)`, `n > 1`: `f₂ = 0`, and
`min {64/(n²(1+b)⁸)·E(p,s) : J₁-constraint} = (n−1)(n−3)(n²−9n+38)/(128n)`,
attained at `p = p*`. -/
theorem c4odd_reduction (n b : ℝ) (hn : 1 < n) (hb : -1 < b) (hb1 : b < 1)
    (hrel : (n - 7) * (1 + b ^ 2) + 2 * (n + 1) * b = 0) :
    f2 n b = 0 ∧ IsLeast (c4oddSet n b) (c4odd n) := by
  have h1 : 1 + b ≠ 0 := by linarith
  have h2 : b - 1 ≠ 0 := by linarith
  refine ⟨(f2_eq_zero_iff n b (by linarith) (by linarith) h1).mpr hrel, ?_, ?_⟩
  · refine ⟨pstar n b, sOf b (pstar n b), (J1con_iff _ _ _ h1).mpr rfl, ?_⟩
    rw [odd_min_value n b hb hb1 hrel]
  · rintro c ⟨p, s, hc, rfl⟩
    rw [(J1con_iff _ _ _ h1).mp hc, E_vertex n b p h2, ← odd_min_value n b hb hb1 hrel]
    have hk : 0 ≤ 64 / (n ^ 2 * (1 + b) ^ 8) := by positivity
    have := a2_pos n b hn hb hb1
    nlinarith [mul_nonneg hk (mul_nonneg this.le (sq_nonneg (p - pstar n b)))]


/-- Paper Thm odd proof: `Q = √P`, `P(w) = w(1−λw)`, at `w₀ = 1/(1+λ)` (`λ ∈ (−1,1)`):
`Q(w₀) = 1/(1+λ)`, `Q'(w₀) = (1−λ)/2`, `Q''(w₀) = −(1+λ)³/4`. Stated as: `Q` has derivative
`Q₁(w) = (1−2λw)/(2√P(w))` wherever `P(w) > 0`, and `Q₁` has derivative `−(1+λ)³/4` at `w₀`. -/
theorem Q_values (b : ℝ) (hb : -1 < b) :
    let P : ℝ → ℝ := fun w => w * (1 - b * w)
    let Q : ℝ → ℝ := fun w => √(P w)
    let Q₁ : ℝ → ℝ := fun w => (1 - 2 * b * w) / (2 * √(P w))
    let w₀ := 1 / (1 + b)
    Q w₀ = 1 / (1 + b) ∧ Q₁ w₀ = (1 - b) / 2 ∧
    (∀ w, 0 < P w → HasDerivAt Q (Q₁ w) w) ∧ HasDerivAt Q₁ (-(1 + b) ^ 3 / 4) w₀ := by
  intro P Q Q₁ w₀
  have h1 : 0 < 1 + b := by linarith
  have hP0 : P w₀ = (1 / (1 + b)) ^ 2 := by
    simp only [P, w₀]; field_simp; ring
  have hsq : √(P w₀) = 1 / (1 + b) := by rw [hP0]; exact Real.sqrt_sq (by positivity)
  have hPd : ∀ w, HasDerivAt P (1 - 2 * b * w) w := by
    intro w
    have := (hasDerivAt_id w).mul ((hasDerivAt_id w).const_mul b |>.const_sub 1)
    convert this using 1
    · funext x; simp [P]
    · simp; ring
  have hQd : ∀ w, 0 < P w → HasDerivAt Q (Q₁ w) w := by
    intro w hw
    have := (hPd w).sqrt hw.ne'
    convert this using 1
  refine ⟨hsq, ?_, hQd, ?_⟩
  · simp only [Q₁]; rw [hsq]; simp only [w₀]; field_simp; ring
  · have hpos : 0 < P w₀ := by rw [hP0]; positivity
    have hnum : HasDerivAt (fun w => 1 - 2 * b * w) (-(2 * b)) w₀ := by
      simpa using ((hasDerivAt_id w₀).const_mul (2 * b)).const_sub 1
    have hden : HasDerivAt (fun w => 2 * √(P w)) (2 * ((1 - 2 * b * w₀) / (2 * √(P w₀)))) w₀ :=
      ((hPd w₀).sqrt hpos.ne').const_mul 2
    have hne : 2 * √(P w₀) ≠ 0 := by rw [hsq]; positivity
    have := hnum.div hden hne
    convert this using 1
    rw [hsq]; simp only [w₀]; field_simp; ring

/-! ## §5. At the actual threshold `b_n`: closed forms, positivity, the table -/

/-- For even `n = 2M ≥ 4`, `b_n` satisfies `2(2M−1)(1+b²) = M²(1+b)²`. -/
theorem even_rel_bn (M : ℕ) (hM : 2 ≤ M) :
    2 * (2 * (M : ℝ) - 1) * (1 + bn (2 * M) ^ 2) = (M : ℝ) ^ 2 * (1 + bn (2 * M)) ^ 2 := by
  obtain ⟨h1, h2, h3⟩ := bn_spec (2 * M) (by omega)
  unfold rfun rhs36 at h3
  rw [floor_sq_div_four_even] at h3
  have hb : (1 + bn (2 * M)) ^ 2 ≠ 0 := pow_ne_zero 2 (by linarith)
  have hM' : (2 : ℝ) ≤ M := by exact_mod_cast hM
  push_cast at h3
  rw [div_eq_div_iff hb (by nlinarith)] at h3
  linarith

/-- For odd `n ≥ 5`, `b_n` satisfies `(n−7)(1+b²) + 2(n+1)b = 0`. -/
theorem odd_rel_bn (M : ℕ) (hM : 2 ≤ M) :
    (((2 * M + 1 : ℕ) : ℝ) - 7) * (1 + bn (2 * M + 1) ^ 2) + 2 * (((2 * M + 1 : ℕ) : ℝ) + 1) * bn (2 * M + 1) = 0 := by
  obtain ⟨h1, h2, h3⟩ := bn_spec (2 * M + 1) (by omega)
  have hb : 1 + bn (2 * M + 1) ≠ 0 := by linarith
  rw [odd_rel_iff _ _ hb, h3]
  unfold rhs36
  rw [floor_sq_div_four_odd]
  have hM' : (2 : ℝ) ≤ M := by exact_mod_cast hM
  push_cast
  rw [div_eq_div_iff (by nlinarith) (by norm_num)]
  ring

/-- **Thm even at `b_n`** (all even `n ≥ 4`): `h₂(b_n) = 0` and
`c₄ = h₄(2/(1+b_n)²)⁴/n² = (n−2)(m⁴+40m²+48m+16)/(192 n(n−1))`. -/
theorem c4_even_at_bn (M : ℕ) (hM : 2 ≤ M) :
    h2 M (bn (2 * M)) = 0 ∧ c4evenOf M (bn (2 * M)) = c4even ((2 * M : ℕ) : ℝ) := by
  obtain ⟨h1, h2, -⟩ := bn_spec (2 * M) (by omega)
  have hM' : (2 : ℝ) ≤ M := by exact_mod_cast hM
  have := c4even_reduction M (bn (2 * M)) (by positivity) (by linarith) (by linarith) (even_rel_bn M hM)
  push_cast
  exact this

/-- **Thm odd at `b_n`** (all odd `n ≥ 5`): `f₂(b_n) = 0` and the constrained minimum of
`64/(n²(1+b_n)⁸)E` is `(n−1)(n−3)(n²−9n+38)/(128n)`. -/
theorem c4_odd_at_bn (M : ℕ) (hM : 2 ≤ M) :
    f2 ((2 * M + 1 : ℕ) : ℝ) (bn (2 * M + 1)) = 0 ∧
    IsLeast (c4oddSet ((2 * M + 1 : ℕ) : ℝ) (bn (2 * M + 1))) (c4odd ((2 * M + 1 : ℕ) : ℝ)) := by
  obtain ⟨h1, h2, -⟩ := bn_spec (2 * M + 1) (by omega)
  have hM' : (2 : ℝ) ≤ M := by exact_mod_cast hM
  exact c4odd_reduction _ _ (by push_cast; linarith) h1 h2 (odd_rel_bn M hM)

/-- **Sign facts (Thm even, Thm odd: "in particular `c₄ > 0`")**. -/
theorem c4even_pos (n : ℝ) (hn : 2 < n) : 0 < c4even n := by
  unfold c4even
  have hm : 0 < n - 2 := by linarith
  apply div_pos
  · apply mul_pos hm; positivity
  · apply mul_pos (mul_pos (by norm_num) (by linarith)) (by linarith)

theorem c4odd_pos (n : ℝ) (hn : 3 < n) : 0 < c4odd n := by
  unfold c4odd
  have hq : 0 < n ^ 2 - 9 * n + 38 := by nlinarith [sq_nonneg (n - 9 / 2)]
  apply div_pos
  · apply mul_pos (mul_pos (by linarith) (by linarith)) hq
  · linarith

/-- `n² − 9n + 38` has discriminant `−71` (paper Thm odd). -/
theorem odd_quadratic_disc : (9 : ℤ) ^ 2 - 4 * 38 = -71 := by norm_num

/-- The even-`n` quartic factor in `m = n − 2` is positive for every `m ≥ 0`. -/
theorem even_quartic_pos (m : ℝ) (hm : 0 ≤ m) : 0 < m ^ 4 + 40 * m ^ 2 + 48 * m + 16 := by positivity

/-- The closed form as a function of `n : ℕ`. -/
def c4closed (n : ℕ) : ℚ := if n % 2 = 0 then c4even (n : ℚ) else c4odd (n : ℚ)

theorem c4even_cast (n : ℕ) : ((c4even (n : ℚ) : ℚ) : ℝ) = c4even (n : ℝ) := by
  unfold c4even; push_cast; ring
theorem c4odd_cast (n : ℕ) : ((c4odd (n : ℚ) : ℚ) : ℝ) = c4odd (n : ℝ) := by
  unfold c4odd; push_cast; ring

/-- "The paper's defining expression for `c₄` at `b_n` equals `c`": for even `n` this is
`h₄(2/(1+b_n)²)⁴/n² = c`, for odd `n` it is "`c` is the least value of `64/(n²(1+b_n)⁸)E(p,s)`
on the `J₁` constraint". -/
def C4AtThreshold (n : ℕ) (c : ℝ) : Prop :=
  if n % 2 = 0 then c4evenOf ((n / 2 : ℕ) : ℝ) (bn n) = c
  else IsLeast (c4oddSet (n : ℝ) (bn n)) c

/-- **Thm even + Thm odd, for every `n ≥ 4`**: the defining expression at `b_n` equals the closed form,
and it is positive. -/
theorem c4_at_bn (n : ℕ) (hn : 4 ≤ n) :
    C4AtThreshold n (c4closed n : ℝ) ∧ 0 < c4closed n := by
  unfold C4AtThreshold c4closed
  rcases Nat.mod_two_eq_zero_or_one n with he | ho
  · obtain ⟨M, rfl⟩ : ∃ M, n = 2 * M := ⟨n / 2, by omega⟩
    simp only [he, ↓reduceIte, show 2 * M / 2 = M by omega]
    refine ⟨?_, ?_⟩
    · rw [c4even_cast]; exact (c4_even_at_bn M (by omega)).2
    · have : (0 : ℝ) < c4even ((2 * M : ℕ) : ℝ) := c4even_pos _ (by
        have : (4 : ℝ) ≤ ((2 * M : ℕ) : ℝ) := by exact_mod_cast hn
        linarith)
      rw [← c4even_cast] at this; exact_mod_cast this
  · obtain ⟨M, rfl⟩ : ∃ M, n = 2 * M + 1 := ⟨n / 2, by omega⟩
    simp only [ho, one_ne_zero, ↓reduceIte]
    refine ⟨?_, ?_⟩
    · rw [c4odd_cast]; exact (c4_odd_at_bn M (by omega)).2
    · have : (0 : ℝ) < c4odd ((2 * M + 1 : ℕ) : ℝ) := c4odd_pos _ (by
        have : (5 : ℝ) ≤ ((2 * M + 1 : ℕ) : ℝ) := by exact_mod_cast (by omega : 5 ≤ 2 * M + 1)
        linarith)
      rw [← c4odd_cast] at this; exact_mod_cast this

/-- **Table 1** (paper), exact over `ℚ`: the closed form at `n = 4, …, 25`. -/
theorem table_values :
    c4closed 4 = 1/4 ∧ c4closed 5 = 9/40 ∧ c4closed 6 = 23/30 ∧ c4closed 7 = 9/14 ∧
    c4closed 8 = 95/56 ∧ c4closed 9 = 19/12 ∧ c4closed 10 = 49/15 ∧ c4closed 11 = 75/22 ∧
    c4closed 12 = 755/132 ∧ c4closed 13 = 675/104 ∧ c4closed 14 = 1693/182 ∧ c4closed 15 = 56/5 ∧
    c4closed 16 = 1141/80 ∧ c4closed 17 = 609/34 ∧ c4closed 18 = 3190/153 ∧ c4closed 19 = 27 ∧
    c4closed 20 = 11139/380 ∧ c4closed 21 = 2175/56 ∧ c4closed 22 = 6145/154 ∧
    c4closed 23 = 2475/46 ∧ c4closed 24 = 29183/552 ∧ c4closed 25 = 7227/100 := by
  simp only [c4closed, c4even, c4odd]
  norm_num

/-- **Table 1 at the threshold**: e.g. at `n = 8` (hyperbolic plane, `b₈ = −(8 − 3√7)`) the
defining expression `h₄(2/(1+b₈)²)⁴/64` is exactly `95/56`, and at `n = 7` (plane, `b₇ = 0`) the
constrained minimum is exactly `9/14`. -/
theorem table_n8 : C4AtThreshold 8 (95 / 56) := by
  have := (c4_at_bn 8 (by norm_num)).1
  rw [table_values.2.2.2.2.1] at this
  convert this using 1; norm_num

theorem table_n7 : C4AtThreshold 7 (9 / 14) := by
  have := (c4_at_bn 7 (by norm_num)).1
  rw [table_values.2.2.2.1] at this
  convert this using 1; norm_num

/-! ## §6. Root-of-unity sums used in §3 and in the proof of Thm odd (all `n`) -/

section RootsOfUnity

open Complex Finset

theorem normSq_one_sub_exp (θ : ℝ) :
    Complex.normSq (1 - Complex.exp (θ * I)) = 2 - 2 * Real.cos θ := by
  rw [Complex.exp_mul_I, Complex.normSq_apply]
  simp [Complex.cos_ofReal_re, Complex.sin_ofReal_re, Complex.cos_ofReal_im, Complex.sin_ofReal_im]
  nlinarith [Real.sin_sq_add_cos_sq θ]

theorem sum_zpow_geom {ζ : ℂ} {n : ℕ} (hζ : IsPrimitiveRoot ζ n) (k : ℤ) (hk : ¬ (n : ℤ) ∣ k) :
    ∑ l ∈ range n, (ζ ^ k) ^ l = 0 := by
  have h1 : ζ ^ k ≠ 1 := by rwa [Ne, hζ.zpow_eq_one_iff_dvd]
  have h2 : (ζ ^ k) ^ n = 1 := by
    rw [← zpow_natCast, ← zpow_mul, mul_comm, zpow_mul, zpow_natCast, hζ.pow_eq_one, one_zpow]
  rw [geom_sum_eq h1, h2, sub_self, zero_div]

theorem sum_normSq_geom (n m : ℕ) (hn : 0 < n) (hm : m ≤ n) :
    ∑ l ∈ range n, Complex.normSq (∑ a ∈ range m, (Complex.exp (2 * π * I / n) ^ l) ^ a)
      = (n : ℝ) * m := by
  set ζ := Complex.exp (2 * π * I / n) with hζdef
  have hζ : IsPrimitiveRoot ζ n := Complex.isPrimitiveRoot_exp n (by omega)
  have hζ0 : ζ ≠ 0 := Complex.exp_ne_zero _
  have hconj : (starRingEnd ℂ) ζ = ζ⁻¹ := by
    rw [hζdef, ← Complex.exp_conj, ← Complex.exp_neg]
    congr 1
    simp [map_div₀, map_ofNat, neg_div]
  apply Complex.ofReal_injective
  push_cast
  have key : ∀ l, ((Complex.normSq (∑ a ∈ range m, (ζ ^ l) ^ a) : ℝ) : ℂ)
      = ∑ a ∈ range m, ∑ b ∈ range m, (ζ ^ ((a : ℤ) - b)) ^ l := by
    intro l
    rw [← Complex.mul_conj, map_sum, Finset.sum_mul_sum]
    apply Finset.sum_congr rfl; intro a _
    apply Finset.sum_congr rfl; intro b _
    rw [map_pow, map_pow, hconj, zpow_sub₀ hζ0, zpow_natCast, zpow_natCast, div_eq_mul_inv,
      mul_pow, ← pow_mul, ← pow_mul, inv_pow, inv_pow, ← pow_mul, ← pow_mul, mul_comm l a, mul_comm l b]
  simp_rw [key]
  rw [Finset.sum_comm]
  simp_rw [Finset.sum_comm (s := range n)]
  have inner : ∀ a ∈ range m, ∑ b ∈ range m, ∑ l ∈ range n, (ζ ^ ((a : ℤ) - b)) ^ l = (n : ℂ) := by
    intro a ha
    rw [Finset.sum_eq_single a]
    · simp
    · intro b hb hba
      apply sum_zpow_geom hζ
      intro hdvd
      have ha' := Finset.mem_range.mp ha
      have hb' := Finset.mem_range.mp hb
      have : (a : ℤ) - b = 0 := by
        apply Int.eq_zero_of_abs_lt_dvd hdvd
        rw [abs_lt]; constructor <;> omega
      omega
    · intro h; exact absurd ha h
  rw [Finset.sum_congr rfl inner]
  simp [mul_comm]

/-- **Fejér-type identity** (paper eq. (fejer), equivalent to MT (2.4)), for every `n ≥ 1` and
`0 ≤ m ≤ n`: `Σ_{l=1}^{n−1} (1 − cos(2πml/n))/(1 − cos(2πl/n)) = m(n − m)`. -/
theorem fejer (n m : ℕ) (hn : 0 < n) (hm : m ≤ n) :
    ∑ l ∈ Ico 1 n, (1 - Real.cos (2 * π * m * l / n)) / (1 - Real.cos (2 * π * l / n))
      = (m : ℝ) * ((n : ℝ) - m) := by
  set ζ := Complex.exp (2 * π * I / n) with hζdef
  have hζ : IsPrimitiveRoot ζ n := Complex.isPrimitiveRoot_exp n (by omega)
  have hpow : ∀ k : ℕ, ∀ l : ℕ, ζ ^ (k * l) = Complex.exp (((2 * π * k * l / n : ℝ)) * I) := by
    intro k l
    rw [hζdef, ← Complex.exp_nat_mul]
    congr 1; push_cast; ring
  have term : ∀ l ∈ Ico 1 n, (1 - Real.cos (2 * π * m * l / n)) / (1 - Real.cos (2 * π * l / n))
      = Complex.normSq (∑ a ∈ range m, (ζ ^ l) ^ a) := by
    intro l hl
    have hl' := Finset.mem_Ico.mp hl
    have hw : ζ ^ l ≠ 1 := by
      intro h
      have := hζ.dvd_of_pow_eq_one l h
      have := Nat.le_of_dvd (by omega) this
      omega
    have e1 : ζ ^ l = Complex.exp (((2 * π * l / n : ℝ)) * I) := by
      have := hpow 1 l; simp only [one_mul, Nat.cast_one, mul_one] at this; exact this
    have e2 : (ζ ^ l) ^ m = Complex.exp (((2 * π * m * l / n : ℝ)) * I) := by
      rw [← pow_mul, mul_comm, hpow]
    have hsum : ∑ a ∈ range m, (ζ ^ l) ^ a = (1 - (ζ ^ l) ^ m) / (1 - ζ ^ l) := by
      rw [eq_div_iff (sub_ne_zero.mpr (Ne.symm hw)), geom_sum_mul_neg]
    rw [hsum, Complex.normSq_div, e2, e1, normSq_one_sub_exp, normSq_one_sub_exp]
    rw [show (2:ℝ) - 2 * Real.cos (2 * π * m * l / n) = 2 * (1 - Real.cos (2 * π * m * l / n)) by ring,
      show (2:ℝ) - 2 * Real.cos (2 * π * l / n) = 2 * (1 - Real.cos (2 * π * l / n)) by ring,
      mul_div_mul_left _ _ two_ne_zero]
  rw [Finset.sum_congr rfl term]
  have htot := sum_normSq_geom n m hn hm
  rw [Finset.range_eq_Ico, Finset.sum_eq_sum_Ico_succ_bot hn] at htot
  simp only [pow_zero, one_pow, Finset.sum_const, Finset.card_range, nsmul_eq_mul, mul_one] at htot
  rw [← hζdef] at htot
  have h0 : Complex.normSq (m : ℂ) = (m : ℝ) ^ 2 := by
    rw [Complex.normSq_natCast]; ring
  rw [h0] at htot
  linear_combination htot


theorem sum_cos_eq_zero (n j : ℕ) (hn : 0 < n) (hj : ¬ n ∣ j) :
    ∑ k ∈ range n, Real.cos (2 * π * j * k / n) = 0 := by
  have hζ : IsPrimitiveRoot (Complex.exp (2 * π * I / n)) n := Complex.isPrimitiveRoot_exp n (by omega)
  have h := sum_zpow_geom hζ (j : ℤ) (by exact_mod_cast hj)
  have h2 : ∀ k : ℕ, (Complex.exp (2 * π * I / n) ^ (j : ℤ)) ^ k
      = Complex.exp (((2 * π * j * k / n : ℝ)) * I) := by
    intro k
    rw [zpow_natCast, ← pow_mul, ← Complex.exp_nat_mul]
    congr 1; push_cast; ring
  simp_rw [h2] at h
  have := congrArg Complex.re h
  rw [Complex.re_sum] at this
  simp only [Complex.exp_ofReal_mul_I_re] at this
  exact this

/-- Paper Thm odd proof (`J₁` constraint): `Σ_k cos²(lφ_k) cos φ_k = n/4`, `l = (n−1)/2`, odd `n ≥ 3`. -/
theorem sum_cos_sq_mul_cos (n : ℕ) (ho : n % 2 = 1) (h3 : 3 ≤ n) :
    ∑ k ∈ range n, Real.cos (2 * π * ((n - 1) / 2 : ℕ) * k / n) ^ 2 * Real.cos (2 * π * k / n)
      = n / 4 := by
  have hn : (0 : ℝ) < n := by positivity
  have hl : (((n - 1) / 2 : ℕ) : ℝ) = ((n : ℝ) - 1) / 2 := by
    obtain ⟨M, rfl⟩ : ∃ M, n = 2 * M + 1 := ⟨n / 2, by omega⟩
    rw [show (2 * M + 1 - 1) / 2 = M by omega]; push_cast; ring
  have term : ∀ k ∈ range n,
      Real.cos (2 * π * ((n - 1) / 2 : ℕ) * k / n) ^ 2 * Real.cos (2 * π * k / n)
        = Real.cos (2 * π * (1 : ℕ) * k / n) / 2 + 1 / 4
          + Real.cos (2 * π * ((n - 2 : ℕ)) * k / n) / 4 := by
    intro k _
    set x := 2 * π * ((n - 1) / 2 : ℕ) * k / n
    set y := 2 * π * k / n
    have hxy1 : 2 * x + y = (k : ℕ) * (2 * π) := by
      simp only [x, y]; rw [hl]; field_simp; ring
    have hxy2 : 2 * x - y = 2 * π * ((n - 2 : ℕ) : ℝ) * k / n := by
      simp only [x, y]; rw [hl, Nat.cast_sub (by omega)]; field_simp; push_cast; ring
    have hc : Real.cos (2 * x) * Real.cos y = (Real.cos (2 * x + y) + Real.cos (2 * x - y)) / 2 := by
      rw [Real.cos_add, Real.cos_sub]; ring
    rw [Real.cos_sq x, add_mul, mul_comm (Real.cos (2 * x) / 2), mul_div_assoc', mul_comm _ (Real.cos (2 * x)), hc,
      hxy1, hxy2, Real.cos_nat_mul_two_pi]
    simp only [y, Nat.cast_one, mul_one]
    ring
  have e1 := sum_cos_eq_zero n 1 (by omega) (by intro h; have := Nat.le_of_dvd one_pos h; omega)
  have e2 := sum_cos_eq_zero n (n - 2) (by omega)
    (by intro h; have := Nat.le_of_dvd (by omega) h; omega)
  rw [Finset.sum_congr rfl term, Finset.sum_add_distrib, Finset.sum_add_distrib, ← Finset.sum_div,
    e1, ← Finset.sum_div (range n) (fun k : ℕ => Real.cos (2 * π * ((n - 2 : ℕ) : ℝ) * k / n)) 4, e2]
  simp
  ring

end RootsOfUnity

end Havelock.Conj13
