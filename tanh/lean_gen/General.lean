import Mathlib

noncomputable section

namespace Havelock.Tanh

/-- `c = cosh r = (1+ξ)/(1-ξ)`, with `ξ = tanh²(r/2)`. -/
def cR (ξ : ℝ) : ℝ := (1 + ξ) / (1 - ξ)

/-- `W_k = Σ_{p ≥ 1, p ≡ k (N)} p ξ^p` in closed form (note §3.1). -/
def W (N k : ℕ) (ξ : ℝ) : ℝ :=
  k * ξ ^ k / (1 - ξ ^ N) + N * ξ ^ (k + N) / (1 - ξ ^ N) ^ 2

/-- `Ω κ₀` at the relative equilibrium (note Lemma 1). -/
def omegaKappa (N : ℕ) (ξ : ℝ) : ℝ :=
  -((N : ℝ) - 1) / 2 + ξ / (1 - ξ) - N * ξ ^ N / (1 - ξ ^ N)

/-- `V''(0) = 2ξ/(1-ξ)² − 2 c Ω κ₀` (proof of note Thm 2, one-body terms). -/
def Vpp (N : ℕ) (ξ : ℝ) : ℝ := 2 * ξ / (1 - ξ) ^ 2 - 2 * cR ξ * omegaKappa N ξ

/-- `a_m` as the sum of the three Hessian contributions (proof of note Thm 2). -/
def aSum (N m : ℕ) (ξ : ℝ) : ℝ :=
  Vpp N ξ - m * ((N : ℝ) - m) / 2 - N * W N N ξ - N / 2 * (W N m ξ + W N (N - m) ξ)

/-- `a_m`, compact form of note Thm 2. -/
def aBlock (N m : ℕ) (ξ : ℝ) : ℝ :=
  ((N : ℝ) - 1) * cR ξ - (cR ξ - 1) ^ 2 / 2 - m * ((N : ℝ) - m) / 2
    + 2 * N * cR ξ * ξ ^ N / (1 - ξ ^ N) - (N : ℝ) ^ 2 * ξ ^ N / (1 - ξ ^ N) ^ 2
    - N / 2 * ((m * ξ ^ m + ((N : ℝ) - m) * ξ ^ (N - m) + ((N : ℝ) - m) * ξ ^ (N + m)
        + m * ξ ^ (2 * N - m)) / (1 - ξ ^ N) ^ 2)

/-- `b_m` (note Thm 2). -/
def bBlock (N m : ℕ) (ξ : ℝ) : ℝ :=
  m * ((N : ℝ) - m) / 2 + N * W N N ξ - N / 2 * (W N m ξ + W N (N - m) ξ)

/-- `g_m` (note Thm 2). -/
def gBlock (N m : ℕ) (ξ : ℝ) : ℝ := N / 2 * (W N m ξ - W N (N - m) ξ)

theorem W_hasSum (N k : ℕ) {ξ : ℝ} (h0 : 0 ≤ ξ) (h1 : ξ < 1) (hN : 0 < N) :
    HasSum (fun q : ℕ => ((k + q * N : ℕ) : ℝ) * ξ ^ (k + q * N)) (W N k ξ) := by
  have hX0 : 0 ≤ ξ ^ N := pow_nonneg h0 N
  have hX1 : ξ ^ N < 1 := pow_lt_one₀ h0 h1 hN.ne'
  have hg := hasSum_geometric_of_lt_one hX0 hX1
  have hc := hasSum_coe_mul_geometric_of_norm_lt_one
    (show ‖ξ ^ N‖ < 1 by rw [Real.norm_of_nonneg hX0]; exact hX1)
  have := ((hg.mul_left (k : ℝ)).add (hc.mul_left (N : ℝ))).mul_left (ξ ^ k)
  have hne : (1 : ℝ) - ξ ^ N ≠ 0 := by linarith
  convert this using 1
  · funext q
    push_cast
    rw [pow_add, pow_mul']
    ring
  · unfold W
    rw [pow_add]
    field_simp

theorem W_self (N : ℕ) (ξ : ℝ) (hX : ξ ^ N ≠ 1) :
    W N N ξ = N * ξ ^ N / (1 - ξ ^ N) ^ 2 := by
  have hne : (1 : ℝ) - ξ ^ N ≠ 0 := sub_ne_zero.mpr (Ne.symm hX)
  unfold W
  rw [pow_add]
  field_simp
  ring

theorem Vpp_eq (N : ℕ) (ξ : ℝ) (h1 : ξ ≠ 1) (hX : ξ ^ N ≠ 1) :
    Vpp N ξ = ((N : ℝ) - 1) * cR ξ - (cR ξ - 1) ^ 2 / 2 + 2 * N * cR ξ * ξ ^ N / (1 - ξ ^ N) := by
  have hne : (1 : ℝ) - ξ ^ N ≠ 0 := sub_ne_zero.mpr (Ne.symm hX)
  have hne1 : (1 : ℝ) - ξ ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  unfold Vpp omegaKappa cR
  field_simp
  ring

theorem aSum_eq_aBlock (N m : ℕ) (hm : m ≤ N) (ξ : ℝ) (h1 : ξ ≠ 1) (hX : ξ ^ N ≠ 1) :
    aSum N m ξ = aBlock N m ξ := by
  obtain ⟨k, rfl⟩ : ∃ k, N = m + k := ⟨N - m, by omega⟩
  have hne : (1 : ℝ) - ξ ^ (m + k) ≠ 0 := sub_ne_zero.mpr (Ne.symm hX)
  have hne1 : (1 : ℝ) - ξ ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  have e1 : m + k - m = k := by omega
  have e2 : 2 * (m + k) - m = m + k + k := by omega
  unfold aSum aBlock Vpp omegaKappa W cR
  rw [e1, e2]
  push_cast
  simp only [pow_add]
  field_simp
  ring


/-- Reciprocity (note Prop 5): `a_m(1/ξ) = a_m(ξ)` as rational functions. -/
theorem aBlock_recip (N m : ℕ) (hm : m ≤ N) (ξ : ℝ) (h0 : ξ ≠ 0) (h1 : ξ ≠ 1)
    (hX : ξ ^ N ≠ 1) : aBlock N m ξ⁻¹ = aBlock N m ξ := by
  obtain ⟨k, rfl⟩ : ∃ k, N = m + k := ⟨N - m, by omega⟩
  have hne : (1 : ℝ) - ξ ^ m * ξ ^ k ≠ 0 := by rw [← pow_add]; exact sub_ne_zero.mpr (Ne.symm hX)
  have hne' : ξ ^ m * ξ ^ k - 1 ≠ 0 := by rw [← pow_add]; exact sub_ne_zero.mpr hX
  have hne1 : (1 : ℝ) - ξ ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  have hne1' : ξ - 1 ≠ 0 := sub_ne_zero.mpr h1
  have e1 : m + k - m = k := by omega
  have e2 : 2 * (m + k) - m = m + k + k := by omega
  unfold aBlock cR
  rw [e1, e2]
  push_cast
  simp only [inv_pow, pow_add]
  field_simp
  ring

/-- `a_{N-m} = a_m` (note Thm 2). -/
theorem aBlock_symm (N m : ℕ) (hm : m ≤ N) (ξ : ℝ) : aBlock N (N - m) ξ = aBlock N m ξ := by
  obtain ⟨k, rfl⟩ : ∃ k, N = m + k := ⟨N - m, by omega⟩
  have e1 : m + k - m = k := by omega
  have e2 : m + k - k = m := by omega
  have e3 : 2 * (m + k) - m = m + k + k := by omega
  have e4 : 2 * (m + k) - k = m + k + m := by omega
  unfold aBlock
  rw [e1, e2, e3, e4]
  push_cast
  ring_nf

/-- `b_{N-m} = b_m` (note Thm 2). -/
theorem bBlock_symm (N m : ℕ) (hm : m ≤ N) (ξ : ℝ) : bBlock N (N - m) ξ = bBlock N m ξ := by
  have e : N - (N - m) = m := by omega
  unfold bBlock
  rw [e, Nat.cast_sub hm]
  ring

/-- `g_{N-m} = -g_m` (note Thm 2). -/
theorem gBlock_symm (N m : ℕ) (hm : m ≤ N) (ξ : ℝ) : gBlock N (N - m) ξ = -gBlock N m ξ := by
  have e : N - (N - m) = m := by omega
  unfold gBlock
  rw [e]
  ring

/-- `g_{N/2} = 0` for even `N` (note Thm 2). -/
theorem gBlock_half (n : ℕ) (ξ : ℝ) : gBlock (2 * n) n ξ = 0 := by
  have e : 2 * n - n = n := by omega
  unfold gBlock
  rw [e]
  ring

/-- Even-`N` threshold equation (note eq. (even)): with `y = ξ^{N/2}`,
`a_{N/2} = (N−1)c − (c−1)²/2 + 2Ncy²/(1−y²) − N²/8 − N²y/(2(1−y)²)`. -/
theorem aBlock_even (n : ℕ) (ξ : ℝ) (hy : ξ ^ n ≠ 1) (hy' : ξ ^ n ≠ -1) :
    aBlock (2 * n) n ξ =
      ((2 * n : ℕ) - 1 : ℝ) * cR ξ - (cR ξ - 1) ^ 2 / 2
        + 2 * (2 * n : ℕ) * cR ξ * (ξ ^ n) ^ 2 / (1 - (ξ ^ n) ^ 2)
        - ((2 * n : ℕ) : ℝ) ^ 2 / 8 - ((2 * n : ℕ) : ℝ) ^ 2 * ξ ^ n / (2 * (1 - ξ ^ n) ^ 2) := by
  have e1 : 2 * n - n = n := by omega
  have e2 : 2 * (2 * n) - n = n + n + n := by omega
  have e3 : 2 * n + n = n + n + n := by omega
  have h1 : (1 : ℝ) - ξ ^ n ≠ 0 := sub_ne_zero.mpr (Ne.symm hy)
  have h2 : (1 : ℝ) + ξ ^ n ≠ 0 := by
    intro h; apply hy'; linarith
  have hX : (1 : ℝ) - ξ ^ n * ξ ^ n ≠ 0 := by
    have : (1 : ℝ) - ξ ^ n * ξ ^ n = (1 - ξ ^ n) * (1 + ξ ^ n) := by ring
    rw [this]; exact mul_ne_zero h1 h2
  have hX2 : (1 : ℝ) - (ξ ^ n) ^ 2 ≠ 0 := by rw [sq]; exact hX
  unfold aBlock
  rw [e1, e2, e3, two_mul n, pow_add, pow_add]
  push_cast
  field_simp
  ring

/-- Linearisation (note Prop 3): the characteristic polynomial of `J M_m`,
`J = [[0,1],[-1,0]]`, `M_m = [[a, i g],[-i g, b]]`, is `(λ + i g)² + a b`. -/
theorem charpoly_block (a b g l : ℂ) :
    Matrix.det (!![0, 1; -1, 0] * !![a, Complex.I * g; -Complex.I * g, b] - l • (1 : Matrix (Fin 2) (Fin 2) ℂ))
      = (l + Complex.I * g) ^ 2 + a * b := by
  simp [Matrix.det_fin_two]
  ring_nf

/-- Self-reciprocity of the numerator `A_m = 2(1−ξ)²(1−ξ^N)² a_m` (note §5):
`A_m(1/ξ) = ξ^{−(2N+2)} A_m(ξ)`. -/
theorem numerA_recip (N m : ℕ) (hm : m ≤ N) (ξ : ℝ) (h0 : ξ ≠ 0) (h1 : ξ ≠ 1)
    (hX : ξ ^ N ≠ 1) :
    2 * (1 - ξ⁻¹) ^ 2 * (1 - ξ⁻¹ ^ N) ^ 2 * aBlock N m ξ⁻¹
      = (ξ ^ (2 * N + 2))⁻¹ * (2 * (1 - ξ) ^ 2 * (1 - ξ ^ N) ^ 2 * aBlock N m ξ) := by
  rw [aBlock_recip N m hm ξ h0 h1 hX]
  have hXn : ξ ^ N ≠ 0 := pow_ne_zero _ h0
  rw [inv_pow, pow_add, pow_mul]
  field_simp
  ring


/-! ### Exact polynomial certificates -/

/-- Horner evaluation of an integer coefficient list (constant term first). -/
def peval : List ℤ → ℝ → ℝ
  | [], _ => 0
  | c :: l, x => c + x * peval l x

/-- Formal derivative of `peval l`, evaluated at `x`. -/
def dpeval : List ℤ → ℝ → ℝ
  | [], _ => 0
  | _ :: l, x => peval l x + x * dpeval l x

theorem hasDerivAt_peval : ∀ (l : List ℤ) (x : ℝ), HasDerivAt (peval l) (dpeval l x) x
  | [], x => by
    have e : peval [] = fun _ => (0 : ℝ) := by funext y; rfl
    rw [e]; simpa [dpeval] using hasDerivAt_const x (0 : ℝ)
  | c :: l, x => by
    have h := ((hasDerivAt_id x).mul (hasDerivAt_peval l x)).const_add (c : ℝ)
    have e : peval (c :: l) = fun y => (c : ℝ) + id y * peval l y := by funext y; rfl
    rw [e]; simpa [dpeval] using h

theorem continuous_peval (l : List ℤ) : Continuous (peval l) :=
  continuous_iff_continuousAt.mpr fun x => (hasDerivAt_peval l x).continuousAt

theorem peval_strictAntiOn (l : List ℤ) (u v : ℝ)
    (h : ∀ x, u < x → x < v → dpeval l x < 0) : StrictAntiOn (peval l) (Set.Icc u v) := by
  apply strictAntiOn_of_deriv_neg (convex_Icc u v) (continuous_peval l).continuousOn
  intro x hx
  rw [interior_Icc] at hx
  rw [(hasDerivAt_peval l x).deriv]
  exact h x hx.1 hx.2

/-- Generic sign-change lemma: if `f` is positive on `(0,L]`, continuous and strictly decreasing on
`[L,1]`, and `f α > 0 > f β` with `L ≤ α < β < 1`, then `f` has a unique zero `ξs` in `(0,1)`; it lies
in `(α,β)`, `f > 0` on `(0,ξs)` and `f < 0` on `(ξs,1)`. -/
theorem sign_change (f : ℝ → ℝ) (L α β : ℝ) (hLα : L ≤ α) (hαβ : α < β) (hβ1 : β < 1)
    (hc : ContinuousOn f (Set.Icc L 1)) (hm : StrictAntiOn f (Set.Icc L 1))
    (hl : ∀ x, 0 < x → x ≤ L → 0 < f x) (hfa : 0 < f α) (hfb : f β < 0) :
    ∃ ξs, α < ξs ∧ ξs < β ∧ f ξs = 0 ∧
      ∀ x, 0 < x → x < 1 → (x < ξs → 0 < f x) ∧ (ξs < x → f x < 0) := by
  have hsub : Set.Icc α β ⊆ Set.Icc L 1 := Set.Icc_subset_Icc hLα hβ1.le
  obtain ⟨ξs, hξ, hf⟩ : ∃ ξs ∈ Set.Icc α β, f ξs = 0 :=
    intermediate_value_Icc' hαβ.le (hc.mono hsub) ⟨hfb.le, hfa.le⟩
  have hξL : ξs ∈ Set.Icc L 1 := hsub hξ
  have hne1 : ξs ≠ α := by rintro rfl; linarith
  have hne2 : ξs ≠ β := by rintro rfl; linarith
  refine ⟨ξs, lt_of_le_of_ne hξ.1 (Ne.symm hne1), lt_of_le_of_ne hξ.2 hne2, hf,
    fun x hx0 hx1 => ⟨fun hlt => ?_, fun hgt => ?_⟩⟩
  · rcases le_or_gt x L with h | h
    · exact hl x hx0 h
    · have := hm ⟨h.le, hx1.le⟩ hξL hlt
      linarith
  · have := hm hξL ⟨by linarith [hξL.1], hx1.le⟩ hgt
    linarith

/-- Positivity test used by the certificates. -/
example (p q : ℝ) (hp : 0 < p) (hq : 0 < q) : 0 < 3 * p ^ 2 * q + 5 * q ^ 3 + p := by positivity

end Havelock.Tanh
