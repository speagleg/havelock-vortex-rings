import Mathlib.NumberTheory.ZetaValues
import Mathlib.Analysis.Complex.RealDeriv
import Havelock.Tanh.LargeN.Defs

/-! # Tanh model, large `N`: lattice sums for the horocyclic row (note Prop. `prop:row`)

The four classical sums used in the proof of note Prop. `prop:row` (subsection "Physical
interpretation"), all proved here without contour integrals:

1. `∑_{k≠0} 1/k² = π²/3` (`hasSum_inv_sq_int`);
2. `∑_{k≥1} cos(qk)/k² = π²/6 − πq/2 + q²/4` for `0 ≤ q ≤ 2π` (`hasSum_cos_div_sq`), and its
   `ℤ`-indexed exponential form `hasSum_inv_sq_fourier` (both from Mathlib's Bernoulli Fourier series);
3. `∑_{k∈ℤ} 1/(k² + a²) = (π/a) coth(πa)` (`hasSum_inv_sq_add_sq'`, equivalently
   `∑_n 1/(π²n² + τ²) = coth τ/τ`, `hasSum_inv_sq_add_sq`);
4. `∑_{k∈ℤ} e^{2πikθ}/(kℓ − 2i)² = −τ² w(θ)` for `0 ≤ θ ≤ 1`, `τ = 2π/ℓ`, with
   `w(θ) = (θ/(1−X) + X/(1−X)²) e^{−2τθ}`, `X = e^{−2τ}` (`hasSum_rowP`). With `b = 2/ℓ`, `q = 2πθ` this
   is the note's `∑_k e^{iqk}/(k − ib)² = −2π ∑_{n≥0}(q + 2πn) e^{−b(q+2πn)}` divided by `ℓ²` (the
   geometric series in `n` summed in closed form).

Sums 3 and 4 are obtained from one generic statement, `hasSum_of_fourierCoeff`: if `F` is continuous
with `F 0 = F 1`, and its Fourier coefficients on `[0,1]` are an absolutely summable `c`, then
`∑ c_n e^{2πinθ} = F θ` on `[0,1]` (Mathlib's `has_pointwise_sum_fourier_series_of_summable` on the unit
circle). The coefficients are elementary integrals of (affine) × exponential (`integral_exp_affine`).
-/

noncomputable section
namespace Havelock.Tanh.LargeN

open Real
open Complex (I)

/-- Pointwise Fourier series on `[0,1]` from absolutely summable coefficients: if `F` is continuous,
`F 0 = F 1`, and `∫₀¹ e^{−2πinx} F x dx = c n` with `c` summable, then `∑_n c n e^{2πinθ} = F θ`
for every `θ ∈ [0,1]`. -/
theorem hasSum_of_fourierCoeff {F : ℝ → ℂ} (hF : Continuous F) (h01 : F 0 = F 1) {c : ℤ → ℂ}
    (hcoef : ∀ n : ℤ, ∫ x in (0 : ℝ)..1, Complex.exp (-(2 * π * I * n * x)) * F x = c n)
    (hs : Summable c) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (fun n : ℤ => c n * Complex.exp (2 * π * I * n * θ)) (F θ) := by
  suffices key : ∀ {y : ℝ}, y ∈ Set.Ico (0 : ℝ) 1 →
      HasSum (fun n : ℤ => c n * Complex.exp (2 * π * I * n * y)) (F y) by
    rw [← Set.Ico_insert_right (zero_le_one' ℝ), Set.mem_insert_iff, or_comm] at hθ
    rcases hθ with hθ | rfl
    · exact key hθ
    · have h0 := key (Set.left_mem_Ico.mpr zero_lt_one)
      rw [← h01]
      convert h0 using 2 with n
      have := Complex.exp_int_mul_two_pi_mul_I n
      push_cast
      rw [mul_one, mul_zero, Complex.exp_zero, ← this]
      ring_nf
  intro y hy
  let f : C(UnitAddCircle, ℂ) :=
    ⟨AddCircle.liftIco 1 0 F, AddCircle.liftIco_zero_continuous h01 hF.continuousOn⟩
  have hc : ∀ n : ℤ, fourierCoeff f n = c n := by
    intro n
    rw [← hcoef n]
    change fourierCoeff (AddCircle.liftIco 1 0 F) n = _
    rw [fourierCoeff_liftIco_eq, fourierCoeffOn_eq_integral]
    simp only [fourier_coe_apply, smul_eq_mul, zero_add, sub_zero, div_one, one_smul]
    congr 1
    ext x
    congr 2
    push_cast
    ring
  have h := has_pointwise_sum_fourier_series_of_summable (f := f)
    (hs.congr fun n => (hc n).symm) (y : UnitAddCircle)
  simp_rw [hc] at h
  convert h using 1
  · funext n
    rw [smul_eq_mul, fourier_coe_apply, Complex.ofReal_one, div_one]
  · change F y = AddCircle.liftIco 1 0 F (y : UnitAddCircle)
    rw [AddCircle.liftIco_zero_coe_apply hy]

/-- `∫₀¹ e^{−cx}(Ax + B) dx` in closed form (`c ≠ 0`). -/
theorem integral_exp_affine {c A B : ℂ} (hc : c ≠ 0) :
    ∫ x in (0 : ℝ)..1, Complex.exp (-(c * x)) * (A * x + B) =
      A / c ^ 2 + B / c - Complex.exp (-c) * (A / c + A / c ^ 2 + B / c) := by
  have hd : ∀ x : ℝ, HasDerivAt
      (fun x : ℝ => -Complex.exp (-(c * x)) * (A * x / c + A / c ^ 2 + B / c))
      (Complex.exp (-(c * x)) * (A * x + B)) x := by
    intro x
    have h1 : HasDerivAt (fun x : ℝ => (x : ℂ)) 1 x := by
      simpa using (hasDerivAt_id x).ofReal_comp
    have h2 : HasDerivAt (fun x : ℝ => Complex.exp (-(c * x)))
        (Complex.exp (-(c * x)) * (-(c * 1))) x :=
      ((h1.const_mul c).neg).cexp
    have h3 : HasDerivAt (fun x : ℝ => A * x / c + A / c ^ 2 + B / c) (A * 1 / c) x :=
      (((h1.const_mul A).div_const c).add_const _).add_const _
    convert (h2.neg).mul h3 using 1
    simp only [Pi.neg_apply]
    field_simp
    ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ => hd x)]
  · simp only [Complex.ofReal_one, mul_one, Complex.ofReal_zero, mul_zero, neg_zero,
      Complex.exp_zero, zero_div, zero_add]
    ring
  · apply Continuous.intervalIntegrable
    fun_prop

/-- Remove the `k = 0` term of an `ℤ`-indexed sum. -/
theorem hasSum_ite_zero {E : Type*} [AddCommGroup E] [TopologicalSpace E] [IsTopologicalAddGroup E]
    {f : ℤ → E} {a : E} (hf : HasSum f a) :
    HasSum (fun k : ℤ => if k = 0 then 0 else f k) (a - f 0) := by
  convert hf.update 0 0 using 1
  · funext k
    by_cases hk : k = 0
    · subst hk; simp
    · simp [hk]
  · abel

/-- Summability of `1/(a n² + b)` over `ℤ`, `a, b > 0`. -/
theorem summable_one_div_sq_add {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    Summable (fun n : ℤ => 1 / (a * (n : ℝ) ^ 2 + b)) := by
  have h1 : Summable (fun n : ℤ => 1 / (n : ℝ) ^ 2) := Real.summable_one_div_int_pow.mpr one_lt_two
  have h2 : Summable (fun n : ℤ => if n = 0 then 1 / b else 0) := (hasSum_ite_eq 0 _).summable
  refine Summable.of_nonneg_of_le (fun n => by positivity) (fun n => ?_) ((h1.mul_left (1 / a)).add h2)
  by_cases hn : n = 0
  · subst hn; simp
  · have hn2 : (0 : ℝ) < (n : ℝ) ^ 2 := by positivity
    simp only [hn, ite_false, add_zero]
    rw [div_mul_div_comm, one_mul]
    apply one_div_le_one_div_of_le (by positivity)
    linarith

/-- The closed form `w(θ) = (θ/(1−X) + X/(1−X)²) e^{−2τθ}`, `X = e^{−2τ}`, of the bracket in the note's
`P(θ)` (proof of note Prop. `prop:row`). -/
def rowW (τ θ : ℝ) : ℝ :=
  (θ / (1 - Real.exp (-2 * τ)) + Real.exp (-2 * τ) / (1 - Real.exp (-2 * τ)) ^ 2) * Real.exp (-2 * τ * θ)

/-- `1 − e^{−2τ} ≠ 0` for `τ > 0`. -/
theorem one_sub_exp_ne {τ : ℝ} (hτ : 0 < τ) : 1 - Real.exp (-2 * τ) ≠ 0 := by
  have : Real.exp (-2 * τ) < 1 := by simpa using Real.exp_lt_exp.mpr (show -2 * τ < 0 by linarith)
  linarith

/-- `w(0) = w(1)`: the periodic extension of `w` is continuous. -/
theorem rowW_zero_eq_one {τ : ℝ} (hτ : 0 < τ) : rowW τ 0 = rowW τ 1 := by
  have h := one_sub_exp_ne hτ
  simp only [rowW, mul_zero, Real.exp_zero, mul_one, zero_div, zero_add]
  field_simp
  ring

/-- Lattice sum (4) of the note (proof of Prop. `prop:row`), via the Fourier series of the continuous
`1`-periodic extension of `−τ² w(θ)`: for `ℓ > 0`, `τ = 2π/ℓ`, `0 ≤ θ ≤ 1`,
`∑_{k∈ℤ} e^{2πikθ}/(kℓ − 2i)² = −τ² w(θ)`. (With `b = 2/ℓ`, `q = 2πθ` this is
`∑_k e^{iqk}/(k − ib)² = −2π ∑_{n≥0}(q+2πn) e^{−b(q+2πn)}` divided by `ℓ²`.) -/
theorem hasSum_rowP {ℓ : ℝ} (hℓ : 0 < ℓ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (fun k : ℤ => 1 / ((k : ℂ) * ℓ - 2 * I) ^ 2 * Complex.exp (2 * π * I * k * θ))
      ((-(2 * π / ℓ) ^ 2 * rowW (2 * π / ℓ) θ : ℝ) : ℂ) := by
  set τ := 2 * π / ℓ with hτdef
  have hτ : 0 < τ := by positivity
  have h1X := one_sub_exp_ne hτ
  set X := Real.exp (-2 * τ) with hX
  refine hasSum_of_fourierCoeff (F := fun x => ((-τ ^ 2 * rowW τ x : ℝ) : ℂ)) ?_ ?_ ?_ ?_ hθ
  · apply Complex.continuous_ofReal.comp
    unfold rowW; fun_prop
  · simp only [rowW_zero_eq_one hτ]
  · intro n
    set c : ℂ := 2 * τ + 2 * π * I * n with hc
    have hc0 : c ≠ 0 := by
      intro h0
      have := congrArg Complex.re h0
      simp [hc] at this
      linarith
    have hint : ∀ x : ℝ, Complex.exp (-(2 * π * I * n * x)) * ((-τ ^ 2 * rowW τ x : ℝ) : ℂ) =
        Complex.exp (-(c * x)) * ((-τ ^ 2 / (1 - X) : ℝ) * x + (-τ ^ 2 * X / (1 - X) ^ 2 : ℝ)) := by
      intro x
      simp only [rowW, ← hX]
      push_cast
      rw [Complex.ofReal_exp]
      have : -(c * x) = -(2 * π * I * n * x) + ((-2 * τ * x : ℝ) : ℂ) := by
        rw [hc]; push_cast; ring
      rw [this, Complex.exp_add]
      push_cast
      field_simp
      ring
    simp_rw [hint]
    rw [integral_exp_affine hc0]
    have hec : Complex.exp (-c) = X := by
      have : -c = ((-2 * τ : ℝ) : ℂ) + (-n : ℤ) * (2 * π * I) := by rw [hc]; push_cast; ring
      rw [this, Complex.exp_add, Complex.exp_int_mul_two_pi_mul_I, mul_one, ← Complex.ofReal_exp]
    rw [hec]
    have h1X' : (1 : ℂ) - X ≠ 0 := by exact_mod_cast h1X
    have step : ((-τ ^ 2 / (1 - X) : ℝ) : ℂ) / c ^ 2 + ((-τ ^ 2 * X / (1 - X) ^ 2 : ℝ) : ℂ) / c
        - (X : ℂ) * (((-τ ^ 2 / (1 - X) : ℝ) : ℂ) / c + ((-τ ^ 2 / (1 - X) : ℝ) : ℂ) / c ^ 2
          + ((-τ ^ 2 * X / (1 - X) ^ 2 : ℝ) : ℂ) / c) = -(τ : ℂ) ^ 2 / c ^ 2 := by
      push_cast
      field_simp
      ring
    rw [step]
    have hℓ' : (ℓ : ℂ) ≠ 0 := by exact_mod_cast hℓ.ne'
    have hcc : c = (τ : ℂ) * I * ((n : ℂ) * ℓ - 2 * I) := by
      rw [hc, hτdef]; push_cast; field_simp; ring_nf; rw [Complex.I_sq]; ring
    have hτ' : (τ : ℂ) ≠ 0 := by exact_mod_cast hτ.ne'
    have hm : (n : ℂ) * ℓ - 2 * I ≠ 0 := by
      intro h0
      have := congrArg Complex.im h0
      simp at this
    rw [hcc, mul_pow, mul_pow, Complex.I_sq]
    field_simp
  · refine (summable_one_div_sq_add (a := ℓ ^ 2) (b := 4) (by positivity) (by norm_num)).of_norm_bounded
      (fun n => le_of_eq ?_)
    rw [norm_div, norm_one, norm_pow, ← Complex.normSq_eq_norm_sq, Complex.normSq_apply]
    simp
    ring

/-- Lattice sum (3) of the note in the form used for the velocity: for `τ > 0`,
`∑_{n∈ℤ} τ sinh τ/(π²n² + τ²) e^{2πinθ} = cosh(τ(1−2θ))` on `[0,1]`. -/
theorem hasSum_cosh_fourier {τ : ℝ} (hτ : 0 < τ) {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (fun n : ℤ => ((τ * Real.sinh τ / (π ^ 2 * n ^ 2 + τ ^ 2) : ℝ) : ℂ) *
      Complex.exp (2 * π * I * n * θ)) ((Real.cosh (τ * (1 - 2 * θ)) : ℝ) : ℂ) := by
  refine hasSum_of_fourierCoeff (F := fun x => ((Real.cosh (τ * (1 - 2 * x)) : ℝ) : ℂ)) ?_ ?_ ?_ ?_ hθ
  · fun_prop
  · simp only [mul_zero, sub_zero, mul_one]
    rw [show τ * (1 - 2) = -τ by ring, Real.cosh_neg]
  · intro n
    set z : ℂ := π * n * I with hz
    have hz2 : z ^ 2 = -((π : ℂ) * n) ^ 2 := by rw [hz, mul_pow, Complex.I_sq]; ring
    set E : ℂ := Complex.exp τ with hE
    have hE0 : E ≠ 0 := Complex.exp_ne_zero _
    have hc1 : (2 * τ + 2 * z : ℂ) ≠ 0 := by
      intro h0; have := congrArg Complex.re h0; simp [hz] at this; linarith
    have hc2 : (-2 * τ + 2 * z : ℂ) ≠ 0 := by
      intro h0; have := congrArg Complex.re h0; simp [hz] at this; linarith
    have hint : ∀ x : ℝ, Complex.exp (-(2 * π * I * n * x)) * ((Real.cosh (τ * (1 - 2 * x)) : ℝ) : ℂ) =
        Complex.exp (-((2 * τ + 2 * z) * x)) * (0 * x + E / 2) +
          Complex.exp (-((-2 * τ + 2 * z) * x)) * (0 * x + E⁻¹ / 2) := by
      intro x
      have e1 : Complex.exp (-((2 * τ + 2 * z) * x)) * E =
          Complex.exp (-(2 * π * I * n * x)) * Complex.exp (τ * (1 - 2 * x)) := by
        rw [hE, ← Complex.exp_add, ← Complex.exp_add]; congr 1; rw [hz]; ring
      have e2 : Complex.exp (-((-2 * τ + 2 * z) * x)) * E⁻¹ =
          Complex.exp (-(2 * π * I * n * x)) * Complex.exp (-(τ * (1 - 2 * x))) := by
        rw [hE, ← Complex.exp_neg, ← Complex.exp_add, ← Complex.exp_add]; congr 1; rw [hz]; ring
      rw [Real.cosh_eq]
      push_cast
      linear_combination (-1 / 2 : ℂ) * e1 + (-1 / 2 : ℂ) * e2
    simp_rw [hint]
    rw [intervalIntegral.integral_add (Continuous.intervalIntegrable (by fun_prop) _ _)
      (Continuous.intervalIntegrable (by fun_prop) _ _), integral_exp_affine hc1,
      integral_exp_affine hc2]
    have he1 : Complex.exp (-(2 * τ + 2 * z)) = E⁻¹ ^ 2 := by
      have : -(2 * τ + 2 * z) = 2 * (-(τ : ℂ)) + (-n : ℤ) * (2 * π * I) := by rw [hz]; push_cast; ring
      rw [this, Complex.exp_add, Complex.exp_int_mul_two_pi_mul_I, mul_one, hE, ← Complex.exp_neg,
        ← Complex.exp_nat_mul]
      push_cast; ring_nf
    have he2 : Complex.exp (-(-2 * τ + 2 * z)) = E ^ 2 := by
      have : -(-2 * τ + 2 * z) = 2 * (τ : ℂ) + (-n : ℤ) * (2 * π * I) := by rw [hz]; push_cast; ring
      rw [this, Complex.exp_add, Complex.exp_int_mul_two_pi_mul_I, mul_one, hE,
        ← Complex.exp_nat_mul]
      push_cast; ring_nf
    rw [he1, he2]
    have hden : ((π ^ 2 * n ^ 2 + τ ^ 2 : ℝ) : ℂ) = τ ^ 2 - z ^ 2 := by
      rw [hz2]; push_cast; ring
    have hden0 : (τ : ℂ) ^ 2 - z ^ 2 ≠ 0 := by
      rw [← hden]; exact_mod_cast (by positivity : (0 : ℝ) < π ^ 2 * n ^ 2 + τ ^ 2).ne'
    have hsh : ((Real.sinh τ : ℝ) : ℂ) = (E - E⁻¹) / 2 := by
      rw [Real.sinh_eq]; push_cast; rw [hE, Complex.exp_neg]
    rw [Complex.ofReal_div, Complex.ofReal_mul, hden, hsh]
    have h2 : (2 * τ + 2 * z : ℂ) = 2 * (τ + z) := by ring
    have h3 : (-2 * τ + 2 * z : ℂ) = -(2 * (τ - z)) := by ring
    have h4 : (τ : ℂ) ^ 2 - z ^ 2 = (τ + z) * (τ - z) := by ring
    have hp : (τ : ℂ) + z ≠ 0 := by intro h0; apply hc1; rw [h2, h0, mul_zero]
    have hm : (τ : ℂ) - z ≠ 0 := by intro h0; apply hc2; rw [h3, h0, mul_zero, neg_zero]
    rw [h2, h3, h4]
    field_simp
    ring
  · have hs := (summable_one_div_sq_add (a := π ^ 2) (b := τ ^ 2) (by positivity)
      (by positivity)).mul_left (τ * Real.sinh τ)
    refine (Complex.summable_ofReal.mpr hs).congr (fun n => ?_)
    push_cast; ring

/-- Lattice sum (3) of the note: `∑_{n∈ℤ} 1/(π²n² + τ²) = coth τ / τ` for `τ > 0`. -/
theorem hasSum_inv_sq_add_sq {τ : ℝ} (hτ : 0 < τ) :
    HasSum (fun n : ℤ => 1 / (π ^ 2 * (n : ℝ) ^ 2 + τ ^ 2)) (coth τ / τ) := by
  have h := hasSum_cosh_fourier hτ (θ := 0) (by simp)
  simp only [Complex.ofReal_zero, mul_zero, Complex.exp_zero, mul_one, sub_zero] at h
  have h' := (Complex.hasSum_ofReal.mp h).div_const (τ * Real.sinh τ)
  have hs : Real.sinh τ ≠ 0 := (Real.sinh_pos_iff.mpr hτ).ne'
  convert h' using 1
  · funext n; field_simp
  · rw [coth]; field_simp

/-- Lattice sum (3) in the classical form: `∑_{k∈ℤ} 1/(k² + a²) = (π/a) coth(πa)` for `a > 0`. -/
theorem hasSum_inv_sq_add_sq' {a : ℝ} (ha : 0 < a) :
    HasSum (fun k : ℤ => 1 / ((k : ℝ) ^ 2 + a ^ 2)) (π / a * coth (π * a)) := by
  have h := (hasSum_inv_sq_add_sq (mul_pos pi_pos ha)).mul_left (π ^ 2)
  convert h using 1
  · funext k; field_simp
  · field_simp

/-- Lattice sums (1)–(2) of the note over `ℤ` (Mathlib's Fourier series of `B₂`): for `0 ≤ θ ≤ 1`,
`∑_{k∈ℤ, k≠0} e^{2πikθ}/k² = 2π²(θ² − θ + 1/6)` (the `k = 0` term is `1/0 = 0` in Lean). -/
theorem hasSum_inv_sq_fourier {θ : ℝ} (hθ : θ ∈ Set.Icc (0 : ℝ) 1) :
    HasSum (fun k : ℤ => 1 / (k : ℂ) ^ 2 * Complex.exp (2 * π * I * k * θ))
      ((2 * π ^ 2 * (θ ^ 2 - θ + 1 / 6) : ℝ) : ℂ) := by
  have h := hasSum_one_div_pow_mul_fourier_mul_bernoulliFun (k := 2) le_rfl hθ
  convert h using 1
  · funext k; rw [fourier_coe_apply, Complex.ofReal_one, div_one]
  · rw [bernoulliFun_two, mul_pow, mul_pow, Complex.I_sq]; push_cast; ring

/-- Lattice sum (1) of the note: `∑_{k≠0} 1/k² = π²/3` (sum over `ℤ`; the `k = 0` term is `0`). -/
theorem hasSum_inv_sq_int : HasSum (fun k : ℤ => 1 / (k : ℝ) ^ 2) (π ^ 2 / 3) := by
  have h := hasSum_inv_sq_fourier (θ := 0) (by simp)
  simp only [Complex.ofReal_zero, mul_zero, Complex.exp_zero, mul_one] at h
  have e : (fun k : ℤ => 1 / (k : ℂ) ^ 2) = fun k : ℤ => ((1 / (k : ℝ) ^ 2 : ℝ) : ℂ) := by
    funext k; push_cast; rfl
  rw [e] at h
  have h' := Complex.hasSum_ofReal.mp h
  convert h' using 1; ring

/-- Lattice sum (2) of the note: `∑_{k≥1} cos(qk)/k² = π²/6 − πq/2 + q²/4` for `0 ≤ q ≤ 2π`
(the `k = 0` term is `1/0 = 0` in Lean). -/
theorem hasSum_cos_div_sq {q : ℝ} (hq0 : 0 ≤ q) (hq : q ≤ 2 * π) :
    HasSum (fun k : ℕ => Real.cos (q * k) / (k : ℝ) ^ 2) (π ^ 2 / 6 - π * q / 2 + q ^ 2 / 4) := by
  have hθ : q / (2 * π) ∈ Set.Icc (0 : ℝ) 1 :=
    ⟨by positivity, (div_le_one (by positivity)).mpr hq⟩
  have h := hasSum_one_div_nat_pow_mul_cos (k := 1) one_ne_zero hθ
  have h2 : (Polynomial.map (algebraMap ℚ ℝ) (Polynomial.bernoulli 2)).eval (q / (2 * π))
      = bernoulliFun 2 (q / (2 * π)) := rfl
  simp only [mul_one, h2, bernoulliFun_two] at h
  convert h using 1
  · funext k
    rw [show 2 * π * (k : ℝ) * (q / (2 * π)) = q * k by field_simp]
    ring
  · field_simp; norm_num; ring

end Havelock.Tanh.LargeN
