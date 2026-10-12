import Havelock.Sinh.Hessian
import Havelock.Tanh.LargeN.Dynamics
import Havelock.SinhArithmetic

/-! # Sinh model: the stability threshold `cosh² r* = q(N)` (note §2 after Prop. `prop:modes`, §3)

On Fourier mode `m` the Hessian of `F = H − ΩJ` at the ring is `diag(a_m, b_m)` with
`a_m = sinhA N m ξ = (N−1)(1+c²)/2 − m(N−m)/2`, `b_m = sinhB N m = m(N−m)/2`, `c = cR ξ = (1+ξ)/(1−ξ)`
(`Sinh/Hessian.lean`, `sinh_blocks_eigen`). The linearised flow on mode `m` is `κ₀⁻¹ J₀ diag(a_m, b_m)`
(note Prop. `prop:modes`); we use the tanh development's `flowBlock a b g = J₀ [[a, i g], [−i g, b]]` with
`g = 0` (the positive factor `1/κ₀` is dropped, as in `Havelock/Tanh/LargeN/Dynamics.lean`).

* `SpecStable N ξ`: on every mode `1 ≤ m ≤ N−1` all eigenvalues of the linearised flow have zero real part
  (mode `0`, the rotation direction, has `b_0 = 0` and only the eigenvalue `0`).
  `StrictStable N ξ`: `a_m > 0` for all `1 ≤ m ≤ N−1` (Hessian positive definite transverse to rotations).
* `threshold` (note §2–3): for `N ≥ 2`, `0 < ξ < 1`:
  `SpecStable N ξ ↔ ∀ m ∈ [1,N−1], a_m b_m ≥ 0 ↔ ∀ m ∈ [1,N−1], a_m ≥ 0 ↔ a_{⌊N/2⌋} ≥ 0
   ↔ ⌊N/2⌋⌈N/2⌉ ≤ (N−1)(1+c²) ↔ q(N) ≤ c² ↔ ⌊N/2⌋⌈N/2⌉/2 ≤ (N−1)(1+ξ²)/(1−ξ)²`
  (`⌊N/2⌋ = N/2`, `⌈N/2⌉ = (N+1)/2` in `ℕ`; the last form is the note's criterion, with `≥`).
  `threshold_strict`: the same with strict inequalities for `StrictStable` (the note's `>` criterion).
* `cR_tanh_sq`, `threshold_cosh`: with `ξ = tanh²(r/2)`, `c = cosh r`, so `SpecStable ↔ q(N) ≤ cosh² r`.
* `q_le_one`, `one_lt_q`: `q(N) ≤ 1` for `2 ≤ N ≤ 7` (every ring stable, `stable_small`), `q(N) > 1`
  for `N ≥ 8`; `threshold_xiStar`: for `N ≥ 8`, `SpecStable N ξ ↔ ξ*(N) ≤ ξ` with
  `ξ* = (√q − 1)/(√q + 1)` (`xiStar`, the number of note §4 / `SinhArithmetic`), and
  `xiStar_add_inv`: `ξ* + 1/ξ* = c(N)` (`SinhArithmetic.c`).
-/

noncomputable section

namespace Havelock.Sinh

open Havelock.Tanh

/-! ## Stability notions -/

/-- Spectral stability of the ring: on every Fourier mode `1 ≤ m ≤ N−1` every eigenvalue of the linearised
flow `J₀ diag(a_m, b_m)` (up to the factor `1/κ₀ > 0`) has zero real part. -/
def SpecStable (N : ℕ) (ξ : ℝ) : Prop :=
  ∀ m : ℕ, 1 ≤ m → m < N → ∀ σ : ℂ,
    Module.End.HasEigenvalue (Matrix.toLin' (LargeN.flowBlock (sinhA N m ξ) (sinhB N m) 0)) σ →
      σ.re = 0

/-- Strict (energy–momentum) stability condition: `a_m > 0` on every mode `1 ≤ m ≤ N−1`
(with `b_m > 0`, the Hessian of `F` is positive definite transverse to the rotation orbit). -/
def StrictStable (N : ℕ) (ξ : ℝ) : Prop := ∀ m : ℕ, 1 ≤ m → m < N → 0 < sinhA N m ξ

/-! ## Elementary facts -/

theorem sinhB_pos {N m : ℕ} (hm1 : 1 ≤ m) (hmN : m < N) : 0 < sinhB N m := by
  unfold sinhB
  have h1 : (1 : ℝ) ≤ m := by exact_mod_cast hm1
  have h2 : (m : ℝ) + 1 ≤ N := by exact_mod_cast hmN
  have : 0 < (m : ℝ) * ((N : ℝ) - m) := mul_pos (by linarith) (by linarith)
  linarith

/-- `m(N−m) ≤ ⌊N/2⌋(N − ⌊N/2⌋)` for every natural `m` (the most unstable mode is `⌊N/2⌋`). -/
theorem mode_le_half (N m : ℕ) :
    (m : ℝ) * ((N : ℝ) - m) ≤ ((N / 2 : ℕ) : ℝ) * ((N : ℝ) - ((N / 2 : ℕ) : ℝ)) := by
  rcases Nat.even_or_odd' N with ⟨k, rfl | rfl⟩
  · have hk : 2 * k / 2 = k := by omega
    rw [hk]; push_cast
    nlinarith [sq_nonneg ((k : ℝ) - m)]
  · have hk : (2 * k + 1) / 2 = k := by omega
    rw [hk]; push_cast
    rcases le_or_gt m k with h | h
    · have h' : (m : ℝ) ≤ k := by exact_mod_cast h
      nlinarith
    · have h' : (k : ℝ) + 1 ≤ m := by exact_mod_cast h
      nlinarith

theorem sinhA_half_le (N m : ℕ) (ξ : ℝ) : sinhA N (N / 2) ξ ≤ sinhA N m ξ := by
  unfold sinhA
  have := mode_le_half N m
  linarith

/-- `⌈N/2⌉ = N − ⌊N/2⌋`. -/
theorem ceil_half_cast (N : ℕ) : (((N + 1) / 2 : ℕ) : ℝ) = (N : ℝ) - ((N / 2 : ℕ) : ℝ) := by
  have : (N + 1) / 2 = N - N / 2 := by omega
  rw [this, Nat.cast_sub (Nat.div_le_self N 2)]

/-- `⌊N/2⌋⌈N/2⌉ = (N−1)(q(N)+1)` (i.e. `q(N) = ⌊N/2⌋⌈N/2⌉/(N−1) − 1`, note §3). -/
theorem floor_ceil_eq_q (N : ℕ) :
    ((N / 2 : ℕ) : ℝ) * (((N + 1) / 2 : ℕ) : ℝ) = ((N : ℝ) - 1) * ((q N : ℝ) + 1) := by
  unfold q
  rcases Nat.even_or_odd' N with ⟨k, rfl | rfl⟩
  · have hk : 2 * k / 2 = k := by omega
    have hk' : (2 * k + 1) / 2 = k := by omega
    have hodd : ¬ (2 * k % 2 = 1) := by omega
    rw [hk, hk']
    simp only [hodd, ↓reduceIte]
    have hne : (2 * (k : ℝ) - 1) ≠ 0 := by
      intro h
      have h' : (2 * (k : ℝ)) = 1 := by linarith
      norm_cast at h'
      omega
    have h4 : (4 * (2 * (k : ℝ) - 1)) ≠ 0 := mul_ne_zero four_ne_zero hne
    push_cast
    rw [div_add_one h4, mul_div_assoc', eq_div_iff h4]
    ring
  · have hk : (2 * k + 1) / 2 = k := by omega
    have hk' : (2 * k + 1 + 1) / 2 = k + 1 := by omega
    have hodd : (2 * k + 1) % 2 = 1 := by omega
    rw [hk, hk']
    simp only [hodd, ↓reduceIte]
    push_cast
    ring

/-! ## The threshold -/

/-- **Note §2–3 (threshold, `≥` form).** For `N ≥ 2` and `0 < ξ < 1`, with `c = (1+ξ)/(1−ξ)`,
`a_m = (N−1)(1+c²)/2 − m(N−m)/2`, `b_m = m(N−m)/2`, `⌊N/2⌋ = N/2`, `⌈N/2⌉ = (N+1)/2`:
spectral stability ⟺ `a_m b_m ≥ 0` for all `1 ≤ m ≤ N−1` ⟺ `a_m ≥ 0` for all `1 ≤ m ≤ N−1`
⟺ `a_{⌊N/2⌋} ≥ 0` ⟺ `⌊N/2⌋⌈N/2⌉ ≤ (N−1)(1+c²)` ⟺ `q(N) ≤ c²`
⟺ `⌊N/2⌋⌈N/2⌉/2 ≤ (N−1)(1+ξ²)/(1−ξ)²` (the note's criterion `C₁(ξ) ≥ T`). -/
theorem threshold {N : ℕ} (hN : 2 ≤ N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) :
    (SpecStable N ξ ↔ ∀ m : ℕ, 1 ≤ m → m < N → 0 ≤ sinhA N m ξ * sinhB N m)
    ∧ ((∀ m : ℕ, 1 ≤ m → m < N → 0 ≤ sinhA N m ξ * sinhB N m)
        ↔ ∀ m : ℕ, 1 ≤ m → m < N → 0 ≤ sinhA N m ξ)
    ∧ ((∀ m : ℕ, 1 ≤ m → m < N → 0 ≤ sinhA N m ξ) ↔ 0 ≤ sinhA N (N / 2) ξ)
    ∧ (0 ≤ sinhA N (N / 2) ξ
        ↔ ((N / 2 : ℕ) : ℝ) * (((N + 1) / 2 : ℕ) : ℝ) ≤ ((N : ℝ) - 1) * (1 + cR ξ ^ 2))
    ∧ (((N / 2 : ℕ) : ℝ) * (((N + 1) / 2 : ℕ) : ℝ) ≤ ((N : ℝ) - 1) * (1 + cR ξ ^ 2)
        ↔ (q N : ℝ) ≤ cR ξ ^ 2)
    ∧ ((q N : ℝ) ≤ cR ξ ^ 2
        ↔ ((N / 2 : ℕ) : ℝ) * (((N + 1) / 2 : ℕ) : ℝ) / 2 ≤ sinhC1 N ξ) := by
  have hN1 : (1 : ℝ) < N := by exact_mod_cast hN
  have hhalf1 : 1 ≤ N / 2 := by omega
  have hhalfN : N / 2 < N := by omega
  have hC1 : sinhC1 N ξ = ((N : ℝ) - 1) * (1 + cR ξ ^ 2) / 2 := by
    have hne : (1 : ℝ) - ξ ≠ 0 := by linarith
    unfold sinhC1 cR; field_simp; ring
  have hA : sinhA N (N / 2) ξ
      = (((N : ℝ) - 1) * (1 + cR ξ ^ 2) - ((N / 2 : ℕ) : ℝ) * (((N + 1) / 2 : ℕ) : ℝ)) / 2 := by
    unfold sinhA; rw [ceil_half_cast]; ring
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · unfold SpecStable
    exact forall_congr' fun m => forall_congr' fun hm1 => forall_congr' fun hmN =>
      LargeN.flowBlock_re_eigenvalues_iff _ _ 0
  · exact forall_congr' fun m => forall_congr' fun hm1 => forall_congr' fun hmN =>
      ⟨fun h => nonneg_of_mul_nonneg_left h (sinhB_pos hm1 hmN),
       fun h => mul_nonneg h (sinhB_pos hm1 hmN).le⟩
  · exact ⟨fun h => h _ hhalf1 hhalfN, fun h m _ _ => h.trans (sinhA_half_le N m ξ)⟩
  · rw [hA]; constructor <;> intro h <;> linarith
  · rw [floor_ceil_eq_q]
    constructor <;> intro h <;> nlinarith
  · rw [hC1, floor_ceil_eq_q]
    constructor <;> intro h <;> nlinarith

/-- **Note §2–3 (threshold, strict form).** For `N ≥ 2`, `0 < ξ < 1`:
`a_m > 0` for all `1 ≤ m ≤ N−1` ⟺ `a_{⌊N/2⌋} > 0` ⟺ `⌊N/2⌋⌈N/2⌉ < (N−1)(1+c²)` ⟺ `q(N) < c²`
⟺ `(N−1)(1+ξ²)/(1−ξ)² > ⌊N/2⌋⌈N/2⌉/2` (the note's abstract criterion). -/
theorem threshold_strict {N : ℕ} (hN : 2 ≤ N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) :
    (StrictStable N ξ ↔ 0 < sinhA N (N / 2) ξ)
    ∧ (0 < sinhA N (N / 2) ξ
        ↔ ((N / 2 : ℕ) : ℝ) * (((N + 1) / 2 : ℕ) : ℝ) < ((N : ℝ) - 1) * (1 + cR ξ ^ 2))
    ∧ (((N / 2 : ℕ) : ℝ) * (((N + 1) / 2 : ℕ) : ℝ) < ((N : ℝ) - 1) * (1 + cR ξ ^ 2)
        ↔ (q N : ℝ) < cR ξ ^ 2)
    ∧ ((q N : ℝ) < cR ξ ^ 2
        ↔ ((N / 2 : ℕ) : ℝ) * (((N + 1) / 2 : ℕ) : ℝ) / 2 < sinhC1 N ξ) := by
  have hN1 : (1 : ℝ) < N := by exact_mod_cast hN
  have hhalf1 : 1 ≤ N / 2 := by omega
  have hhalfN : N / 2 < N := by omega
  have hC1 : sinhC1 N ξ = ((N : ℝ) - 1) * (1 + cR ξ ^ 2) / 2 := by
    have hne : (1 : ℝ) - ξ ≠ 0 := by linarith
    unfold sinhC1 cR; field_simp; ring
  have hA : sinhA N (N / 2) ξ
      = (((N : ℝ) - 1) * (1 + cR ξ ^ 2) - ((N / 2 : ℕ) : ℝ) * (((N + 1) / 2 : ℕ) : ℝ)) / 2 := by
    unfold sinhA; rw [ceil_half_cast]; ring
  refine ⟨?_, ?_, ?_, ?_⟩
  · exact ⟨fun h => h _ hhalf1 hhalfN, fun h m _ _ => h.trans_le (sinhA_half_le N m ξ)⟩
  · rw [hA]; constructor <;> intro h <;> linarith
  · rw [floor_ceil_eq_q]
    constructor <;> intro h <;> nlinarith
  · rw [hC1, floor_ceil_eq_q]
    constructor <;> intro h <;> nlinarith

/-- Spectral stability ⟺ `q(N) ≤ c²` (the chain of `threshold`, collapsed). -/
theorem specStable_iff {N : ℕ} (hN : 2 ≤ N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) :
    SpecStable N ξ ↔ (q N : ℝ) ≤ cR ξ ^ 2 := by
  obtain ⟨e1, e2, e3, e4, e5, -⟩ := threshold hN h0 h1
  exact e1.trans (e2.trans (e3.trans (e4.trans e5)))

/-- Strict stability ⟺ `q(N) < c²`. -/
theorem strictStable_iff {N : ℕ} (hN : 2 ≤ N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) :
    StrictStable N ξ ↔ (q N : ℝ) < cR ξ ^ 2 := by
  obtain ⟨e1, e2, e3, -⟩ := threshold_strict hN h0 h1
  exact e1.trans (e2.trans e3)

/-! ## The `cosh²` form -/

/-- With `ξ = tanh²(r/2)` (`r` the hyperbolic ring radius), `c = (1+ξ)/(1−ξ) = cosh r` (note §3). -/
theorem cR_tanh_sq (r : ℝ) : cR (Real.tanh (r / 2) ^ 2) = Real.cosh r := by
  have hc : 0 < Real.cosh (r / 2) := Real.cosh_pos _
  have hcs : Real.cosh (r / 2) ^ 2 = Real.sinh (r / 2) ^ 2 + 1 := Real.cosh_sq _
  have h2 : Real.cosh r = Real.cosh (r / 2) ^ 2 + Real.sinh (r / 2) ^ 2 := by
    rw [← Real.cosh_two_mul]; ring_nf
  unfold cR
  rw [Real.tanh_eq_sinh_div_cosh, div_pow, h2]
  have hne : Real.cosh (r / 2) ^ 2 ≠ 0 := by positivity
  have hne2 : 1 - Real.sinh (r / 2) ^ 2 / Real.cosh (r / 2) ^ 2 ≠ 0 := by
    rw [one_sub_div hne, hcs, add_sub_cancel_left]
    exact div_ne_zero one_ne_zero (by positivity)
  rw [div_eq_iff hne2]
  field_simp
  rw [hcs]; ring

theorem tanh_sq_mem (r : ℝ) (hr : 0 < r) :
    0 < Real.tanh (r / 2) ^ 2 ∧ Real.tanh (r / 2) ^ 2 < 1 := by
  have hp : 0 < Real.tanh (r / 2) := by
    rw [Real.tanh_eq_sinh_div_cosh]
    exact div_pos (Real.sinh_pos_iff.mpr (by linarith)) (Real.cosh_pos _)
  have hl : Real.tanh (r / 2) < 1 := Real.tanh_lt_one _
  exact ⟨by positivity, by nlinarith⟩

/-- **Note §3 (`cosh²` form).** For `N ≥ 2` and hyperbolic radius `r > 0` (`ξ = tanh²(r/2)`):
spectral stability ⟺ `cosh² r ≥ q(N)`, and strict stability ⟺ `cosh² r > q(N)`. -/
theorem threshold_cosh {N : ℕ} (hN : 2 ≤ N) {r : ℝ} (hr : 0 < r) :
    (SpecStable N (Real.tanh (r / 2) ^ 2) ↔ (q N : ℝ) ≤ Real.cosh r ^ 2)
    ∧ (StrictStable N (Real.tanh (r / 2) ^ 2) ↔ (q N : ℝ) < Real.cosh r ^ 2) := by
  obtain ⟨h0, h1⟩ := tanh_sq_mem r hr
  rw [specStable_iff hN h0 h1, strictStable_iff hN h0 h1, cR_tanh_sq]
  exact ⟨Iff.rfl, Iff.rfl⟩

/-! ## Connection with `SinhArithmetic`: `q ≤ 1` for `N ≤ 7`, `ξ*` for `N ≥ 8` -/

theorem one_lt_cR {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) : 1 < cR ξ := by
  unfold cR
  rw [one_lt_div (by linarith)]
  linarith

/-- `q(N) ≤ 1` for `2 ≤ N ≤ 7` (note §3). -/
theorem q_le_one {N : ℕ} (h2 : 2 ≤ N) (h7 : N ≤ 7) : q N ≤ 1 := by
  interval_cases N <;> norm_num [q]

/-- `q(N) > 1` for `N ≥ 8` (note §3). -/
theorem one_lt_q {N : ℕ} (h8 : 8 ≤ N) : 1 < q N := by
  unfold q
  split_ifs with h
  · have : (9 : ℚ) ≤ N := by
      have : 9 ≤ N := by omega
      exact_mod_cast this
    rw [lt_div_iff₀ (by norm_num)]
    linarith
  · have hN : (8 : ℚ) ≤ N := by exact_mod_cast h8
    rw [lt_div_iff₀ (by nlinarith)]
    nlinarith

/-- **Note §3.** For `2 ≤ N ≤ 7` every ring `0 < ξ < 1` is spectrally stable. -/
theorem stable_small {N : ℕ} (h2 : 2 ≤ N) (h7 : N ≤ 7) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) :
    SpecStable N ξ := by
  rw [specStable_iff h2 h0 h1]
  have hq : (q N : ℝ) ≤ 1 := by exact_mod_cast q_le_one h2 h7
  have hc := one_lt_cR h0 h1
  nlinarith

/-- `ξ*(N) = (√q − 1)/(√q + 1)` (note §3; the number studied in note §4 and `SinhArithmetic`). -/
def xiStar (N : ℕ) : ℝ := (Real.sqrt (q N) - 1) / (Real.sqrt (q N) + 1)

/-- `ξ* + 1/ξ* = c(N) = 2 + 4/(q(N) − 1)` for `N ≥ 8` (`SinhArithmetic.xi_add_inv`, `SinhArithmetic.c`). -/
theorem xiStar_add_inv {N : ℕ} (h8 : 8 ≤ N) : xiStar N + (xiStar N)⁻¹ = (c N : ℝ) := by
  have hq : (1 : ℝ) < q N := by exact_mod_cast one_lt_q h8
  have := xi_add_inv (q N : ℝ) hq
  simp only at this
  unfold xiStar c
  rw [this]
  push_cast
  ring

theorem cR_le_iff {ξ s : ℝ} (h1 : ξ < 1) (hs : 0 ≤ s) :
    s ≤ cR ξ ↔ (s - 1) / (s + 1) ≤ ξ := by
  unfold cR
  rw [le_div_iff₀ (by linarith), div_le_iff₀ (by linarith)]
  constructor <;> intro h <;> linarith

theorem cR_lt_iff {ξ s : ℝ} (h1 : ξ < 1) (hs : 0 ≤ s) :
    s < cR ξ ↔ (s - 1) / (s + 1) < ξ := by
  unfold cR
  rw [lt_div_iff₀ (by linarith), div_lt_iff₀ (by linarith)]
  constructor <;> intro h <;> linarith

/-- **Note §3 (`ξ*` form).** For `N ≥ 8` and `0 < ξ < 1`: the ring is spectrally stable iff
`ξ ≥ ξ*(N)`, and strictly stable iff `ξ > ξ*(N)`, `ξ* = (√q(N) − 1)/(√q(N) + 1)`. -/
theorem threshold_xiStar {N : ℕ} (h8 : 8 ≤ N) {ξ : ℝ} (h0 : 0 < ξ) (h1 : ξ < 1) :
    (SpecStable N ξ ↔ xiStar N ≤ ξ) ∧ (StrictStable N ξ ↔ xiStar N < ξ) := by
  have h2 : 2 ≤ N := by omega
  have hc : 0 ≤ cR ξ := (zero_lt_one.trans (one_lt_cR h0 h1)).le
  have hq0 : (0 : ℝ) ≤ q N := by
    have : (1 : ℝ) < q N := by exact_mod_cast one_lt_q h8
    linarith
  have hs := Real.sqrt_nonneg (q N : ℝ)
  have hsq := Real.sq_sqrt hq0
  rw [specStable_iff h2 h0 h1, strictStable_iff h2 h0 h1, xiStar,
    ← cR_le_iff h1 hs, ← cR_lt_iff h1 hs]
  constructor
  · constructor
    · intro h; nlinarith
    · intro h; nlinarith
  · constructor
    · intro h; nlinarith
    · intro h; nlinarith

end Havelock.Sinh
