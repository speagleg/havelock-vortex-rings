/-
Havelock — Vortex-ring stability thresholds on H² (sinh model): the arithmetic core.

Setting (proved on paper in sinh/note.tex, NOT formalised here):
for the sinh model the N-ring is marginally stable at cosh² r* = q(N) with
  q(N) = (N−3)/4 (N odd),  q(N) = (N−2)²/(4(N−1)) (N even),
and ξ* = tanh²(r*/2) = (√q − 1)/(√q + 1) satisfies ξ* + 1/ξ* = c(N) := 2 + 4/(q(N) − 1).

Formalised here (no `sorry`):
* `xi_add_inv`            : for real q > 1, ξ = (√q−1)/(√q+1) satisfies ξ + ξ⁻¹ = 2 + 4/(q−1).
* `four_div_odd`          : odd N > 7:  4/(q−1) = 16/(N−7)   (in ℚ).
* `four_div_even`         : even N ≥ 8: 4/(q−1) = 16(N−1)/(N²−8N+8)   (in ℚ).
* `rat_div_int_iff_dvd`   : for d > 0, n/d ∈ ℤ ↔ d ∣ n.
* `odd_dvd_iff`           : odd N > 7:  (N−7) ∣ 16 ↔ N ∈ {9, 11, 15, 23}.
* `even_dvd_iff`          : even N ≥ 8: (N²−8N+8) ∣ 16(N−1) ↔ N = 8.
* `c_int_iff`             : for every N > 7, c(N) ∈ ℤ ↔ N ∈ {8, 9, 11, 15, 23}
                            (with c(8),…,c(23) = 16, 10, 6, 4, 3: `c_values`).
* `isIntegral_of_c_int`   : if c ∈ ℤ and ξ + ξ⁻¹ = c (ξ ≠ 0), then ξ is an algebraic integer
                            (root of the monic X² − cX + 1) and a unit (ξ · (c − ξ) = 1).
* `xi8_eq`, `xi8_root`    : ξ*(8) = (√(9/7) − 1)/(√(9/7) + 1) = 8 − 3√7, and it solves ξ² − 16ξ + 1 = 0.
* `xi10_eq`               : ξ*(10) = 1/7 (q(10) = 16/9) — a norm-1 threshold that is not an algebraic
                            integer (rationality ⇒ not integral is the rational root theorem; stated
                            on paper, `not_isIntegral_one_seventh` formalises it for 1/7).

NOT formalised: the dynamics (Hessian, Fourier block-diagonalisation, C₁ formula), and the converse
"ξ* algebraic integer ⇒ c ∈ ℤ" for irrational ξ* (minimal-polynomial argument, on paper).
-/
import Mathlib

namespace Havelock.Sinh

/-- q(N) as a rational number. -/
def q (N : ℕ) : ℚ :=
  if N % 2 = 1 then ((N : ℚ) - 3) / 4 else ((N : ℚ) - 2) ^ 2 / (4 * ((N : ℚ) - 1))

/-- c(N) = ξ* + 1/ξ* = 2 + 4/(q(N) − 1). -/
def c (N : ℕ) : ℚ := 2 + 4 / (q N - 1)

theorem xi_add_inv (q : ℝ) (hq : 1 < q) :
    let ξ := (Real.sqrt q - 1) / (Real.sqrt q + 1)
    ξ + ξ⁻¹ = 2 + 4 / (q - 1) := by
  intro ξ
  have hs : 1 < Real.sqrt q := by
    rw [show (1:ℝ) = Real.sqrt 1 by simp]
    exact Real.sqrt_lt_sqrt (by norm_num) hq
  have hsq : Real.sqrt q ^ 2 = q := Real.sq_sqrt (by linarith)
  have h1 : Real.sqrt q - 1 ≠ 0 := by linarith
  have h2 : Real.sqrt q + 1 ≠ 0 := by linarith
  have h3 : q - 1 ≠ 0 := by linarith
  simp only [ξ, inv_div]
  rw [div_add_div _ _ h2 h1]
  have h4 : (Real.sqrt q + 1) * (Real.sqrt q - 1) = q - 1 := by nlinarith
  have h5 : (Real.sqrt q - 1) * (Real.sqrt q - 1) + (Real.sqrt q + 1) * (Real.sqrt q + 1)
      = 2 * (q - 1) + 4 := by nlinarith
  rw [h4, h5, add_div, mul_div_assoc, div_self h3, mul_one]

theorem four_div_odd (N : ℕ) (h7 : 7 < N) :
    4 / (((N : ℚ) - 3) / 4 - 1) = 16 / ((N : ℚ) - 7) := by
  have h : (N : ℚ) - 7 ≠ 0 := by
    have : (7 : ℚ) < N := by exact_mod_cast h7
    linarith
  have h' : ((N : ℚ) - 3) / 4 - 1 = ((N : ℚ) - 7) / 4 := by ring
  rw [h', div_div_eq_mul_div]
  ring

theorem four_div_even (N : ℕ) (h8 : 8 ≤ N) :
    4 / (((N : ℚ) - 2) ^ 2 / (4 * ((N : ℚ) - 1)) - 1)
      = 16 * ((N : ℚ) - 1) / ((N : ℚ) ^ 2 - 8 * N + 8) := by
  have hN : (8 : ℚ) ≤ N := by exact_mod_cast h8
  have h1 : (N : ℚ) - 1 ≠ 0 := by linarith
  have h2 : (N : ℚ) ^ 2 - 8 * N + 8 ≠ 0 := by nlinarith
  have h' : ((N : ℚ) - 2) ^ 2 / (4 * ((N : ℚ) - 1)) - 1
      = ((N : ℚ) ^ 2 - 8 * N + 8) / (4 * ((N : ℚ) - 1)) := by
    field_simp
    ring
  rw [h', div_div_eq_mul_div]
  ring

theorem rat_div_int_iff_dvd (n d : ℤ) (hd : 0 < d) :
    (∃ k : ℤ, (n : ℚ) / d = k) ↔ d ∣ n := by
  have hd' : (d : ℚ) ≠ 0 := by exact_mod_cast hd.ne'
  constructor
  · rintro ⟨k, hk⟩
    refine ⟨k, ?_⟩
    have : (n : ℚ) = d * k := by
      rw [div_eq_iff hd'] at hk
      linarith
    exact_mod_cast this
  · rintro ⟨k, rfl⟩
    exact ⟨k, by push_cast; field_simp⟩

theorem odd_dvd_iff (N : ℕ) (h7 : 7 < N) (hodd : N % 2 = 1) :
    ((N : ℤ) - 7) ∣ 16 ↔ N = 9 ∨ N = 11 ∨ N = 15 ∨ N = 23 := by
  constructor
  · intro h
    have hle : (N : ℤ) - 7 ≤ 16 := Int.le_of_dvd (by norm_num) h
    have hN : N ≤ 23 := by omega
    interval_cases N <;> simp_all
  · rintro (rfl | rfl | rfl | rfl) <;> decide

theorem even_dvd_iff (N : ℕ) (h8 : 8 ≤ N) (heven : N % 2 = 0) :
    ((N : ℤ) ^ 2 - 8 * N + 8) ∣ 16 * ((N : ℤ) - 1) ↔ N = 8 := by
  constructor
  · intro h
    by_contra hne
    by_cases hbig : 23 ≤ N
    · have hz : (23 : ℤ) ≤ N := by exact_mod_cast hbig
      have hle := Int.le_of_dvd (by nlinarith) h
      nlinarith
    · have hN : N ≤ 22 := by omega
      interval_cases N <;> simp_all
  · rintro rfl
    decide

/-- Main arithmetic classification: for N > 7, c(N) = ξ* + 1/ξ* is an integer
iff N ∈ {8, 9, 11, 15, 23}. -/
theorem c_int_iff (N : ℕ) (h7 : 7 < N) :
    (∃ k : ℤ, c N = k) ↔ N = 8 ∨ N = 9 ∨ N = 11 ∨ N = 15 ∨ N = 23 := by
  -- c N = 2 + t, so c N ∈ ℤ ↔ t ∈ ℤ
  have shift : ∀ t : ℚ, (∃ k : ℤ, 2 + t = k) ↔ (∃ k : ℤ, t = k) := by
    intro t
    constructor
    · rintro ⟨k, hk⟩; exact ⟨k - 2, by push_cast; linarith⟩
    · rintro ⟨k, hk⟩; exact ⟨k + 2, by push_cast; linarith⟩
  unfold c q
  rw [shift]
  rcases Nat.mod_two_eq_zero_or_one N with h | h
  · -- even N
    have h8 : 8 ≤ N := by omega
    simp only [h, Nat.zero_ne_one, ↓reduceIte]
    rw [four_div_even N h8]
    have hpos : (0 : ℤ) < (N : ℤ) ^ 2 - 8 * N + 8 := by
      have : (8 : ℤ) ≤ N := by exact_mod_cast h8
      nlinarith
    have e := rat_div_int_iff_dvd (16 * ((N : ℤ) - 1)) ((N : ℤ) ^ 2 - 8 * N + 8) hpos
    push_cast at e
    rw [e, even_dvd_iff N h8 h]
    omega
  · -- odd N
    simp only [h, ↓reduceIte]
    rw [four_div_odd N h7]
    have hpos : (0 : ℤ) < (N : ℤ) - 7 := by
      have : (7 : ℤ) < N := by exact_mod_cast h7
      linarith
    have e := rat_div_int_iff_dvd 16 ((N : ℤ) - 7) hpos
    push_cast at e
    rw [e, odd_dvd_iff N h7 h]
    omega

theorem c_values : c 8 = 16 ∧ c 9 = 10 ∧ c 11 = 6 ∧ c 15 = 4 ∧ c 23 = 3 := by
  unfold c q
  norm_num

/-- ξ* is an algebraic integer, and a unit, as soon as c = ξ + ξ⁻¹ is an integer. -/
theorem isIntegral_of_c_int (ξ : ℝ) (hξ : ξ ≠ 0) (k : ℤ) (hk : ξ + ξ⁻¹ = k) :
    IsIntegral ℤ ξ ∧ ξ * ((k : ℝ) - ξ) = 1 := by
  have hroot : ξ ^ 2 - (k : ℝ) * ξ + 1 = 0 := by
    rw [← hk]; field_simp; ring
  refine ⟨⟨Polynomial.X ^ 2 - Polynomial.C k * Polynomial.X + 1, ?_, ?_⟩, ?_⟩
  · -- monic
    have : (Polynomial.X ^ 2 - Polynomial.C k * Polynomial.X + 1 : Polynomial ℤ)
        = Polynomial.X ^ 2 + (Polynomial.C (-k) * Polynomial.X + Polynomial.C 1) := by
      simp [sub_eq_add_neg]; ring
    rw [this]
    exact Polynomial.monic_X_pow_add (by
      calc (Polynomial.C (-k) * Polynomial.X + Polynomial.C 1 : Polynomial ℤ).degree
          ≤ 1 := Polynomial.degree_linear_le
        _ < 2 := by norm_num)
  · simp only [Polynomial.eval₂_add, Polynomial.eval₂_sub, Polynomial.eval₂_mul,
      Polynomial.eval₂_C, Polynomial.eval₂_X, Polynomial.eval₂_pow, Polynomial.eval₂_one]
    simp only [eq_intCast]
    linear_combination hroot
  · linear_combination -hroot

theorem sqrt_nine_sevenths : Real.sqrt (9 / 7) = 3 * Real.sqrt 7 / 7 := by
  have h7 : (0 : ℝ) < Real.sqrt 7 := Real.sqrt_pos.mpr (by norm_num)
  have hs : Real.sqrt 7 ^ 2 = 7 := Real.sq_sqrt (by norm_num)
  rw [Real.sqrt_eq_iff_mul_self_eq (by norm_num) (by positivity)]
  nlinarith

theorem xi8_eq :
    (Real.sqrt (9 / 7) - 1) / (Real.sqrt (9 / 7) + 1) = 8 - 3 * Real.sqrt 7 := by
  rw [sqrt_nine_sevenths]
  have h7 : (0 : ℝ) < Real.sqrt 7 := Real.sqrt_pos.mpr (by norm_num)
  have hs : Real.sqrt 7 ^ 2 = 7 := Real.sq_sqrt (by norm_num)
  rw [div_eq_iff (by positivity)]
  nlinarith

theorem q8 : q 8 = 9 / 7 := by unfold q; norm_num

theorem xi8_root : (8 - 3 * Real.sqrt 7) ^ 2 - 16 * (8 - 3 * Real.sqrt 7) + 1 = 0 := by
  have hs : Real.sqrt 7 ^ 2 = 7 := Real.sq_sqrt (by norm_num)
  nlinarith

theorem q10 : q 10 = 16 / 9 := by unfold q; norm_num

theorem xi10_eq : (Real.sqrt (16 / 9) - 1) / (Real.sqrt (16 / 9) + 1) = 1 / 7 := by
  have : Real.sqrt (16 / 9) = 4 / 3 := by
    rw [Real.sqrt_eq_iff_mul_self_eq (by norm_num) (by norm_num)]; norm_num
  rw [this]; norm_num

/-- 1/7 is not an algebraic integer (ℤ is integrally closed). -/
theorem not_isIntegral_one_seventh : ¬ IsIntegral ℤ (1 / 7 : ℚ) := by
  intro h
  obtain ⟨k, hk⟩ := IsIntegrallyClosed.isIntegral_iff.mp h
  have : (7 : ℚ) * (algebraMap ℤ ℚ k) = 1 := by rw [hk]; norm_num
  simp at this
  have : (7 : ℤ) * k = 1 := by exact_mod_cast this
  omega

end Havelock.Sinh
