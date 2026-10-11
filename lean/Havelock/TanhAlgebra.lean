/-
Havelock — Vortex rings on H² in the tanh (Kimura) model, H = −Σ ln tanh(d/2): the exact-algebra layer.

Source: tanh/note_tanh.tex (theorem numbers below are the LaTeX
counter values: Lemma 1 = lem:re, Thm 2 = thm:blocks, Prop 3 = prop:lin, Prop 5 = prop:recip,
Thm 6 = thm:cert, Cor 7 = "N = 8: an algebraic unit", Prop 9 = "Odd N: the Krein window"),
derivation.txt, exact.py, identities.py, certify_v2.py.

Conventions. ξ = tanh²(r/2) ∈ (0,1), X = ξ^N, c = (1+ξ)/(1−ξ) (`cR`). The Fourier-block entries
are DEFINED here by the closed forms of note Thm 2 (`aBlock`, `bBlock`, `gBlock`, with `W` the closed
form of W_k). Integer polynomials are coefficient lists, constant term first, evaluated by `peval`.

FORMALISED (no `sorry`):
§1 Rational-function identities (all N, all m ≤ N, unless stated)
* `W_hasSum`        : Σ_{q≥0} (k+qN) ξ^{k+qN} = k ξ^k/(1−X) + N ξ^{k+N}/(1−X)²  (W_k closed form, §3.1);
                      `W_self`: W_N = N X/(1−X)².
* `Vpp_eq`          : 2ξ/(1−ξ)² − 2c·Ωκ₀ = (N−1)c − (c−1)²/2 + 2NcX/(1−X), with Ωκ₀ from Lemma 1
                      (proof of Thm 2, one-body terms).
* `aSum_eq_aBlock`  : V''(0) − m(N−m)/2 − N W_N − (N/2)(W_m+W_{N−m}) equals the compact a_m of Thm 2
                      (the "expanding over the common denominator" step; identities.py).
* `aBlock_symm`, `bBlock_symm`, `gBlock_symm`, `gBlock_half` : a_{N−m}=a_m, b_{N−m}=b_m, g_{N−m}=−g_m,
                      g_{N/2}=0 (Thm 2, last sentence).
* `aBlock_even`     : even-N threshold form, note eq. (even), with y = ξ^{N/2}.
* `aBlock_recip`    : reciprocity a_m(1/ξ) = a_m(ξ) (Prop 5);
  `numerA_recip`    : A_m(1/ξ) = ξ^{−(2N+2)} A_m(ξ) for A_m = 2(1−ξ)²(1−X)² a_m (§5, after Prop 5).
* `charpoly_block`  : det(J M_m − λ) = (λ + i g_m)² + a_m b_m (Prop 3, eigenvalue formula).
* `tanh_relations`  : for ξ = tanh²(r/2): c = cosh r and ξ + 1/ξ = 2 + 4/sinh² r (§5, natural invariant).
§2 Certified statements for N = 2..16 (exact Bernstein-basis positivity certificates, checked by
   `ring` + `positivity`; IVT + strict monotonicity via the derivative for the threshold root)
* `stable_N2` … `stable_N7` : Thm 6(a): for every ξ ∈ (0,1) and 1 ≤ m ≤ N−1:
                      a_m > 0, b_m > 0, a_m b_m − g_m² > 0.
* `cert_N8` … `cert_N16`    : Thm 6(b): ∃! ξ* ∈ (0,1) root of f_N, enclosed in a rational interval of
                      width 10⁻²⁰ that matches Table 1; for all ξ ∈ (0,1): every b_m > 0; a_{⌊N/2⌋} < 0
                      for ξ < ξ* (unstable mode ⌊N/2⌋); every a_m > 0 for ξ > ξ* (spectrally stable).
                      Uniqueness = f_N has no other root in (0,1).
* `A_crit_factor_N8` … `_N16` : exact factorisation of the critical numerator
                      2(1−ξ)²(1−X)² a_{⌊N/2⌋} = −k (ξ−1)⁴ · [ (ξ^{N/2}+1) (ξ+1)^{[N≡2 mod 4]} for even N ] · f_N(ξ)
                      (k = 2 for N ≡ 0 mod 4, 1 otherwise). NOTE the single factor (ξ+1), see DISCREPANCY.
* `fN_data`         : f_N is palindromic, has the degree and leading coefficient of Table 1.
* `krein_N9, _N11, _N13, _N15` : Prop 9 (qualitative): just above ξ* there is a rational point where
                      all a_m, b_m > 0 (spectrally stable) but a b − g² < 0 on the critical block,
                      hence ξ_EM > ξ*.
§3 N = 8 (Cor 7)
* `P8_monic`, `P8_natDegree`, `aeval_P8` : P_8 = x¹⁰−12x⁹−…−12x+1 ∈ ℤ[X] is monic of degree 10 and is f_8.
* `f8_root_unit`, `xi8_unit` : ξ*(8) is an algebraic integer whose inverse is an algebraic integer
                      (a unit); C = ξ*+1/ξ* satisfies C⁵−12C⁴−42C³−24C²+32C+32 = 0 (`f8_C_poly`).

NOT FORMALISED (on paper / numerical only):
* The derivation of the Hessian and its Fourier block-diagonalisation (Thm 2 itself: the power-sum
  decomposition (eq:decomp), the series for ln|1−w|, the root-of-unity filter Σ_k w^{(j−k)p} w^{mk}
  = N w^{mj}[p ≡ m], the pair-term Hessian from note.tex). Here the blocks are taken as definitions;
  only the algebra connecting the three contributions and the closed forms is checked.
* Lemma 1 (relative equilibrium) and the dynamics; Prop 3's stability criteria beyond the
  characteristic polynomial; Remark 4 (energy–momentum modulo SO(2,1)); nonlinear stability.
* Irreducibility of f_N over ℚ (hence "minimal polynomial", "degree of ξ*", the degree pattern
  (eq:deg), and "ξ*(N) not an algebraic integer for N = 9..16" are NOT formalised; only that ξ*(N)
  is a root of the stated palindromic f_N with the stated degree/leading coefficient).
* The Krein-window gap sizes and the location of ξ_EM; even-N coincidence of the two thresholds.
* Large-N asymptotics (κ∞) and all numerical fits (§7); the Arb balls and 50-digit values (we give
  20-digit rational enclosures).
-/
import Mathlib

set_option linter.unusedSimpArgs false

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

/-! ### N = 2 -/

theorem N2_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 2 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 2) ^ 2) = (1 - x) ^ 4 * peval [1, 2, 1] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 2 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (2 - 1 : ℕ) = 1 by norm_num, show (2 * 2 - 1 : ℕ) = 3 by norm_num, show (2 + 1 : ℕ) = 3 by norm_num, show (2 - 1 + 2 : ℕ) = 3 by norm_num, show (1 + 2 : ℕ) = 3 by norm_num, show (2 + 2 : ℕ) = 4 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N2_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1, 2, 1] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 1 * q ^ 2 + 4 * p ^ 1 * q ^ 1 + 4 * p ^ 2 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [1, 2, 1] x =
      (fun p q : ℝ => 1 * q ^ 2 + 4 * p ^ 1 * q ^ 1 + 4 * p ^ 2) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [1, 2, 1] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N2_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1, 2, 1] x := N2_posA1_0 x hu hv

theorem N2_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 2 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 2 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N2_posA1 x h0 h1)
  rw [← N2_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N2_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 2 1 x * (2 * (1 - x ^ 2) ^ 2) = (1 - x) ^ 4 * peval [1] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 2 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (2 - 1 : ℕ) = 1 by norm_num, show (2 * 2 - 1 : ℕ) = 3 by norm_num, show (2 + 1 : ℕ) = 3 by norm_num, show (2 - 1 + 2 : ℕ) = 3 by norm_num, show (1 + 2 : ℕ) = 3 by norm_num, show (2 + 2 : ℕ) = 4 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N2_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1] x := by norm_num [peval]

theorem N2_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 2 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 2 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N2_posB1 x h0 h1)
  rw [← N2_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N2_linkD1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    (aBlock 2 1 x * bBlock 2 1 x - gBlock 2 1 x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ 2) ^ 4) = (1 - x) ^ 8 * peval [1, 2, 1] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 2 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock bBlock gBlock W cR
  simp only [show (2 - 1 : ℕ) = 1 by norm_num, show (2 * 2 - 1 : ℕ) = 3 by norm_num, show (2 + 1 : ℕ) = 3 by norm_num, show (2 - 1 + 2 : ℕ) = 3 by norm_num, show (1 + 2 : ℕ) = 3 by norm_num, show (2 + 2 : ℕ) = 4 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N2_posD1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1, 2, 1] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 1 * q ^ 2 + 4 * p ^ 1 * q ^ 1 + 4 * p ^ 2 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [1, 2, 1] x =
      (fun p q : ℝ => 1 * q ^ 2 + 4 * p ^ 1 * q ^ 1 + 4 * p ^ 2) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [1, 2, 1] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N2_posD1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1, 2, 1] x := N2_posD1_0 x hu hv

theorem N2_D1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 2 1 x * bBlock 2 1 x - gBlock 2 1 x ^ 2 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 2 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 8) (N2_posD1 x h0 h1)
  rw [← N2_linkD1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))).mp hn

/-- Note Thm 6(a) for `N = 2`: on all of `(0,1)`, every Fourier block has `a_m > 0`,
`b_m > 0` and `a_m b_m − g_m² > 0` (spectrally stable, positive-definite non-rotational Hessian). -/
theorem stable_N2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (m : ℕ) (hm1 : 1 ≤ m) (hm2 : m ≤ 1) :
    0 < aBlock 2 m x ∧ 0 < bBlock 2 m x ∧
      0 < aBlock 2 m x * bBlock 2 m x - gBlock 2 m x ^ 2 := by
  interval_cases m
  · exact ⟨N2_A1_pos x h0 h1, N2_B1_pos x h0 h1, N2_D1_pos x h0 h1⟩

/-! ### N = 3 -/

theorem N3_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 3 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 3) ^ 2) = (1 - x) ^ 4 * peval [2, 9, 14, 9, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 3 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (3 - 1 : ℕ) = 2 by norm_num, show (2 * 3 - 1 : ℕ) = 5 by norm_num, show (3 + 1 : ℕ) = 4 by norm_num, show (3 - 1 + 3 : ℕ) = 5 by norm_num, show (1 + 3 : ℕ) = 4 by norm_num, show (3 + 3 : ℕ) = 6 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N3_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 9, 14, 9, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 4 + 17 * p ^ 1 * q ^ 3 + 53 * p ^ 2 * q ^ 2 + 72 * p ^ 3 * q ^ 1 + 36 * p ^ 4 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 9, 14, 9, 2] x =
      (fun p q : ℝ => 2 * q ^ 4 + 17 * p ^ 1 * q ^ 3 + 53 * p ^ 2 * q ^ 2 + 72 * p ^ 3 * q ^ 1 + 36 * p ^ 4) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 9, 14, 9, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N3_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 9, 14, 9, 2] x := N3_posA1_0 x hu hv

theorem N3_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 3 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 3 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N3_posA1 x h0 h1)
  rw [← N3_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N3_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 3 1 x * (2 * (1 - x ^ 3) ^ 2) = (1 - x) ^ 4 * peval [2, 5, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 3 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (3 - 1 : ℕ) = 2 by norm_num, show (2 * 3 - 1 : ℕ) = 5 by norm_num, show (3 + 1 : ℕ) = 4 by norm_num, show (3 - 1 + 3 : ℕ) = 5 by norm_num, show (1 + 3 : ℕ) = 4 by norm_num, show (3 + 3 : ℕ) = 6 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N3_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 5, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 2 + 9 * p ^ 1 * q ^ 1 + 9 * p ^ 2 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 5, 2] x =
      (fun p q : ℝ => 2 * q ^ 2 + 9 * p ^ 1 * q ^ 1 + 9 * p ^ 2) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 5, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N3_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 5, 2] x := N3_posB1_0 x hu hv

theorem N3_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 3 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 3 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N3_posB1 x h0 h1)
  rw [← N3_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N3_linkD1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    (aBlock 3 1 x * bBlock 3 1 x - gBlock 3 1 x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ 3) ^ 4) = (1 - x) ^ 8 * peval [4, 28, 68, 88, 68, 28, 4] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 3 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock bBlock gBlock W cR
  simp only [show (3 - 1 : ℕ) = 2 by norm_num, show (2 * 3 - 1 : ℕ) = 5 by norm_num, show (3 + 1 : ℕ) = 4 by norm_num, show (3 - 1 + 3 : ℕ) = 5 by norm_num, show (1 + 3 : ℕ) = 4 by norm_num, show (3 + 3 : ℕ) = 6 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N3_posD1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [4, 28, 68, 88, 68, 28, 4] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 4 * q ^ 6 + 52 * p ^ 1 * q ^ 5 + 268 * p ^ 2 * q ^ 4 + 720 * p ^ 3 * q ^ 3 + 1080 * p ^ 4 * q ^ 2 + 864 * p ^ 5 * q ^ 1 + 288 * p ^ 6 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [4, 28, 68, 88, 68, 28, 4] x =
      (fun p q : ℝ => 4 * q ^ 6 + 52 * p ^ 1 * q ^ 5 + 268 * p ^ 2 * q ^ 4 + 720 * p ^ 3 * q ^ 3 + 1080 * p ^ 4 * q ^ 2 + 864 * p ^ 5 * q ^ 1 + 288 * p ^ 6) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [4, 28, 68, 88, 68, 28, 4] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N3_posD1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [4, 28, 68, 88, 68, 28, 4] x := N3_posD1_0 x hu hv

theorem N3_D1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 3 1 x * bBlock 3 1 x - gBlock 3 1 x ^ 2 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 3 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 8) (N3_posD1 x h0 h1)
  rw [← N3_linkD1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))).mp hn

/-- Note Thm 6(a) for `N = 3`: on all of `(0,1)`, every Fourier block has `a_m > 0`,
`b_m > 0` and `a_m b_m − g_m² > 0` (spectrally stable, positive-definite non-rotational Hessian). -/
theorem stable_N3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (m : ℕ) (hm1 : 1 ≤ m) (hm2 : m ≤ 2) :
    0 < aBlock 3 m x ∧ 0 < bBlock 3 m x ∧
      0 < aBlock 3 m x * bBlock 3 m x - gBlock 3 m x ^ 2 := by
  interval_cases m
  · exact ⟨N3_A1_pos x h0 h1, N3_B1_pos x h0 h1, N3_D1_pos x h0 h1⟩
  · rw [show aBlock 3 2 x = aBlock 3 1 x from aBlock_symm 3 1 (by norm_num) x, show bBlock 3 2 x = bBlock 3 1 x from bBlock_symm 3 1 (by norm_num) x, show gBlock 3 2 x = -gBlock 3 1 x from gBlock_symm 3 1 (by norm_num) x, neg_sq]
    exact ⟨N3_A1_pos x h0 h1, N3_B1_pos x h0 h1, N3_D1_pos x h0 h1⟩

/-! ### N = 4 -/

theorem N4_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 4 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 4) ^ 2) = (1 - x) ^ 4 * peval [3, 14, 33, 44, 33, 14, 3] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 4 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (4 - 1 : ℕ) = 3 by norm_num, show (2 * 4 - 1 : ℕ) = 7 by norm_num, show (4 + 1 : ℕ) = 5 by norm_num, show (4 - 1 + 4 : ℕ) = 7 by norm_num, show (1 + 4 : ℕ) = 5 by norm_num, show (4 + 4 : ℕ) = 8 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N4_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [3, 14, 33, 44, 33, 14, 3] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 3 * q ^ 6 + 32 * p ^ 1 * q ^ 5 + 148 * p ^ 2 * q ^ 4 + 376 * p ^ 3 * q ^ 3 + 548 * p ^ 4 * q ^ 2 + 432 * p ^ 5 * q ^ 1 + 144 * p ^ 6 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [3, 14, 33, 44, 33, 14, 3] x =
      (fun p q : ℝ => 3 * q ^ 6 + 32 * p ^ 1 * q ^ 5 + 148 * p ^ 2 * q ^ 4 + 376 * p ^ 3 * q ^ 3 + 548 * p ^ 4 * q ^ 2 + 432 * p ^ 5 * q ^ 1 + 144 * p ^ 6) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [3, 14, 33, 44, 33, 14, 3] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N4_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [3, 14, 33, 44, 33, 14, 3] x := N4_posA1_0 x hu hv

theorem N4_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 4 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 4 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N4_posA1 x h0 h1)
  rw [← N4_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N4_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 4 1 x * (2 * (1 - x ^ 4) ^ 2) = (1 - x) ^ 4 * peval [3, 8, 14, 8, 3] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 4 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (4 - 1 : ℕ) = 3 by norm_num, show (2 * 4 - 1 : ℕ) = 7 by norm_num, show (4 + 1 : ℕ) = 5 by norm_num, show (4 - 1 + 4 : ℕ) = 7 by norm_num, show (1 + 4 : ℕ) = 5 by norm_num, show (4 + 4 : ℕ) = 8 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N4_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [3, 8, 14, 8, 3] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 3 * q ^ 4 + 20 * p ^ 1 * q ^ 3 + 56 * p ^ 2 * q ^ 2 + 72 * p ^ 3 * q ^ 1 + 36 * p ^ 4 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [3, 8, 14, 8, 3] x =
      (fun p q : ℝ => 3 * q ^ 4 + 20 * p ^ 1 * q ^ 3 + 56 * p ^ 2 * q ^ 2 + 72 * p ^ 3 * q ^ 1 + 36 * p ^ 4) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [3, 8, 14, 8, 3] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N4_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [3, 8, 14, 8, 3] x := N4_posB1_0 x hu hv

theorem N4_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 4 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 4 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N4_posB1 x h0 h1)
  rw [← N4_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N4_linkD1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    (aBlock 4 1 x * bBlock 4 1 x - gBlock 4 1 x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ 4) ^ 4) = (1 - x) ^ 8 * peval [9, 66, 237, 520, 794, 908, 794, 520, 237, 66, 9] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 4 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock bBlock gBlock W cR
  simp only [show (4 - 1 : ℕ) = 3 by norm_num, show (2 * 4 - 1 : ℕ) = 7 by norm_num, show (4 + 1 : ℕ) = 5 by norm_num, show (4 - 1 + 4 : ℕ) = 7 by norm_num, show (1 + 4 : ℕ) = 5 by norm_num, show (4 + 4 : ℕ) = 8 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N4_posD1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [9, 66, 237, 520, 794, 908, 794, 520, 237, 66, 9] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 9 * q ^ 10 + 156 * p ^ 1 * q ^ 9 + 1236 * p ^ 2 * q ^ 8 + 5872 * p ^ 3 * q ^ 7 + 18504 * p ^ 4 * q ^ 6 + 40448 * p ^ 5 * q ^ 5 + 62240 * p ^ 6 * q ^ 4 + 66752 * p ^ 7 * q ^ 3 + 47888 * p ^ 8 * q ^ 2 + 20800 * p ^ 9 * q ^ 1 + 4160 * p ^ 10 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [9, 66, 237, 520, 794, 908, 794, 520, 237, 66, 9] x =
      (fun p q : ℝ => 9 * q ^ 10 + 156 * p ^ 1 * q ^ 9 + 1236 * p ^ 2 * q ^ 8 + 5872 * p ^ 3 * q ^ 7 + 18504 * p ^ 4 * q ^ 6 + 40448 * p ^ 5 * q ^ 5 + 62240 * p ^ 6 * q ^ 4 + 66752 * p ^ 7 * q ^ 3 + 47888 * p ^ 8 * q ^ 2 + 20800 * p ^ 9 * q ^ 1 + 4160 * p ^ 10) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [9, 66, 237, 520, 794, 908, 794, 520, 237, 66, 9] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N4_posD1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [9, 66, 237, 520, 794, 908, 794, 520, 237, 66, 9] x := N4_posD1_0 x hu hv

theorem N4_D1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 4 1 x * bBlock 4 1 x - gBlock 4 1 x ^ 2 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 4 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 8) (N4_posD1 x h0 h1)
  rw [← N4_linkD1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))).mp hn

theorem N4_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 4 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 4) ^ 2) = (1 - x) ^ 4 * peval [2, 16, 22, 32, 22, 16, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 4 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (4 - 2 : ℕ) = 2 by norm_num, show (2 * 4 - 2 : ℕ) = 6 by norm_num, show (4 + 2 : ℕ) = 6 by norm_num, show (4 - 2 + 4 : ℕ) = 6 by norm_num, show (2 + 4 : ℕ) = 6 by norm_num, show (4 + 4 : ℕ) = 8 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N4_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 16, 22, 32, 22, 16, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 6 + 28 * p ^ 1 * q ^ 5 + 132 * p ^ 2 * q ^ 4 + 320 * p ^ 3 * q ^ 3 + 440 * p ^ 4 * q ^ 2 + 336 * p ^ 5 * q ^ 1 + 112 * p ^ 6 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 16, 22, 32, 22, 16, 2] x =
      (fun p q : ℝ => 2 * q ^ 6 + 28 * p ^ 1 * q ^ 5 + 132 * p ^ 2 * q ^ 4 + 320 * p ^ 3 * q ^ 3 + 440 * p ^ 4 * q ^ 2 + 336 * p ^ 5 * q ^ 1 + 112 * p ^ 6) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 16, 22, 32, 22, 16, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N4_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 16, 22, 32, 22, 16, 2] x := N4_posA2_0 x hu hv

theorem N4_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 4 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 4 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N4_posA2 x h0 h1)
  rw [← N4_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N4_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 4 2 x * (2 * (1 - x ^ 4) ^ 2) = (1 - x) ^ 4 * peval [4, 16, 24, 16, 4] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 4 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (4 - 2 : ℕ) = 2 by norm_num, show (2 * 4 - 2 : ℕ) = 6 by norm_num, show (4 + 2 : ℕ) = 6 by norm_num, show (4 - 2 + 4 : ℕ) = 6 by norm_num, show (2 + 4 : ℕ) = 6 by norm_num, show (4 + 4 : ℕ) = 8 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N4_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [4, 16, 24, 16, 4] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 4 * q ^ 4 + 32 * p ^ 1 * q ^ 3 + 96 * p ^ 2 * q ^ 2 + 128 * p ^ 3 * q ^ 1 + 64 * p ^ 4 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [4, 16, 24, 16, 4] x =
      (fun p q : ℝ => 4 * q ^ 4 + 32 * p ^ 1 * q ^ 3 + 96 * p ^ 2 * q ^ 2 + 128 * p ^ 3 * q ^ 1 + 64 * p ^ 4) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [4, 16, 24, 16, 4] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N4_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [4, 16, 24, 16, 4] x := N4_posB2_0 x hu hv

theorem N4_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 4 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 4 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N4_posB2 x h0 h1)
  rw [← N4_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N4_linkD2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    (aBlock 4 2 x * bBlock 4 2 x - gBlock 4 2 x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ 4) ^ 4) = (1 - x) ^ 8 * peval [8, 96, 392, 896, 1392, 1600, 1392, 896, 392, 96, 8] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 4 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock bBlock gBlock W cR
  simp only [show (4 - 2 : ℕ) = 2 by norm_num, show (2 * 4 - 2 : ℕ) = 6 by norm_num, show (4 + 2 : ℕ) = 6 by norm_num, show (4 - 2 + 4 : ℕ) = 6 by norm_num, show (2 + 4 : ℕ) = 6 by norm_num, show (4 + 4 : ℕ) = 8 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N4_posD2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [8, 96, 392, 896, 1392, 1600, 1392, 896, 392, 96, 8] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 8 * q ^ 10 + 176 * p ^ 1 * q ^ 9 + 1616 * p ^ 2 * q ^ 8 + 8448 * p ^ 3 * q ^ 7 + 28384 * p ^ 4 * q ^ 6 + 64832 * p ^ 5 * q ^ 5 + 102848 * p ^ 6 * q ^ 4 + 112640 * p ^ 7 * q ^ 3 + 81920 * p ^ 8 * q ^ 2 + 35840 * p ^ 9 * q ^ 1 + 7168 * p ^ 10 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [8, 96, 392, 896, 1392, 1600, 1392, 896, 392, 96, 8] x =
      (fun p q : ℝ => 8 * q ^ 10 + 176 * p ^ 1 * q ^ 9 + 1616 * p ^ 2 * q ^ 8 + 8448 * p ^ 3 * q ^ 7 + 28384 * p ^ 4 * q ^ 6 + 64832 * p ^ 5 * q ^ 5 + 102848 * p ^ 6 * q ^ 4 + 112640 * p ^ 7 * q ^ 3 + 81920 * p ^ 8 * q ^ 2 + 35840 * p ^ 9 * q ^ 1 + 7168 * p ^ 10) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [8, 96, 392, 896, 1392, 1600, 1392, 896, 392, 96, 8] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N4_posD2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [8, 96, 392, 896, 1392, 1600, 1392, 896, 392, 96, 8] x := N4_posD2_0 x hu hv

theorem N4_D2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 4 2 x * bBlock 4 2 x - gBlock 4 2 x ^ 2 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 4 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 8) (N4_posD2 x h0 h1)
  rw [← N4_linkD2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))).mp hn

/-- Note Thm 6(a) for `N = 4`: on all of `(0,1)`, every Fourier block has `a_m > 0`,
`b_m > 0` and `a_m b_m − g_m² > 0` (spectrally stable, positive-definite non-rotational Hessian). -/
theorem stable_N4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (m : ℕ) (hm1 : 1 ≤ m) (hm2 : m ≤ 3) :
    0 < aBlock 4 m x ∧ 0 < bBlock 4 m x ∧
      0 < aBlock 4 m x * bBlock 4 m x - gBlock 4 m x ^ 2 := by
  interval_cases m
  · exact ⟨N4_A1_pos x h0 h1, N4_B1_pos x h0 h1, N4_D1_pos x h0 h1⟩
  · exact ⟨N4_A2_pos x h0 h1, N4_B2_pos x h0 h1, N4_D2_pos x h0 h1⟩
  · rw [show aBlock 4 3 x = aBlock 4 1 x from aBlock_symm 4 1 (by norm_num) x, show bBlock 4 3 x = bBlock 4 1 x from bBlock_symm 4 1 (by norm_num) x, show gBlock 4 3 x = -gBlock 4 1 x from gBlock_symm 4 1 (by norm_num) x, neg_sq]
    exact ⟨N4_A1_pos x h0 h1, N4_B1_pos x h0 h1, N4_D1_pos x h0 h1⟩

/-! ### N = 5 -/

theorem N5_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 5 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 5) ^ 2) = (1 - x) ^ 4 * peval [4, 19, 46, 81, 100, 81, 46, 19, 4] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 5 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (5 - 1 : ℕ) = 4 by norm_num, show (2 * 5 - 1 : ℕ) = 9 by norm_num, show (5 + 1 : ℕ) = 6 by norm_num, show (5 - 1 + 5 : ℕ) = 9 by norm_num, show (1 + 5 : ℕ) = 6 by norm_num, show (5 + 5 : ℕ) = 10 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N5_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [4, 19, 46, 81, 100, 81, 46, 19, 4] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 4 * q ^ 8 + 51 * p ^ 1 * q ^ 7 + 291 * p ^ 2 * q ^ 6 + 980 * p ^ 3 * q ^ 5 + 2140 * p ^ 4 * q ^ 4 + 3100 * p ^ 5 * q ^ 3 + 2900 * p ^ 6 * q ^ 2 + 1600 * p ^ 7 * q ^ 1 + 400 * p ^ 8 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [4, 19, 46, 81, 100, 81, 46, 19, 4] x =
      (fun p q : ℝ => 4 * q ^ 8 + 51 * p ^ 1 * q ^ 7 + 291 * p ^ 2 * q ^ 6 + 980 * p ^ 3 * q ^ 5 + 2140 * p ^ 4 * q ^ 4 + 3100 * p ^ 5 * q ^ 3 + 2900 * p ^ 6 * q ^ 2 + 1600 * p ^ 7 * q ^ 1 + 400 * p ^ 8) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [4, 19, 46, 81, 100, 81, 46, 19, 4] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N5_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [4, 19, 46, 81, 100, 81, 46, 19, 4] x := N5_posA1_0 x hu hv

theorem N5_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 5 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 5 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N5_posA1 x h0 h1)
  rw [← N5_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N5_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 5 1 x * (2 * (1 - x ^ 5) ^ 2) = (1 - x) ^ 4 * peval [4, 11, 20, 30, 20, 11, 4] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 5 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (5 - 1 : ℕ) = 4 by norm_num, show (2 * 5 - 1 : ℕ) = 9 by norm_num, show (5 + 1 : ℕ) = 6 by norm_num, show (5 - 1 + 5 : ℕ) = 9 by norm_num, show (1 + 5 : ℕ) = 6 by norm_num, show (5 + 5 : ℕ) = 10 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N5_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [4, 11, 20, 30, 20, 11, 4] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 4 * q ^ 6 + 35 * p ^ 1 * q ^ 5 + 135 * p ^ 2 * q ^ 4 + 300 * p ^ 3 * q ^ 3 + 400 * p ^ 4 * q ^ 2 + 300 * p ^ 5 * q ^ 1 + 100 * p ^ 6 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [4, 11, 20, 30, 20, 11, 4] x =
      (fun p q : ℝ => 4 * q ^ 6 + 35 * p ^ 1 * q ^ 5 + 135 * p ^ 2 * q ^ 4 + 300 * p ^ 3 * q ^ 3 + 400 * p ^ 4 * q ^ 2 + 300 * p ^ 5 * q ^ 1 + 100 * p ^ 6) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [4, 11, 20, 30, 20, 11, 4] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N5_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [4, 11, 20, 30, 20, 11, 4] x := N5_posB1_0 x hu hv

theorem N5_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 5 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 5 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N5_posB1 x h0 h1)
  rw [← N5_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N5_linkD1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    (aBlock 5 1 x * bBlock 5 1 x - gBlock 5 1 x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ 5) ^ 4) = (1 - x) ^ 8 * peval [16, 120, 448, 1180, 2336, 3648, 4700, 5104, 4700, 3648, 2336, 1180, 448, 120, 16] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 5 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock bBlock gBlock W cR
  simp only [show (5 - 1 : ℕ) = 4 by norm_num, show (2 * 5 - 1 : ℕ) = 9 by norm_num, show (5 + 1 : ℕ) = 6 by norm_num, show (5 - 1 + 5 : ℕ) = 9 by norm_num, show (1 + 5 : ℕ) = 6 by norm_num, show (5 + 5 : ℕ) = 10 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N5_posD1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [16, 120, 448, 1180, 2336, 3648, 4700, 5104, 4700, 3648, 2336, 1180, 448, 120, 16] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 16 * q ^ 14 + 344 * p ^ 1 * q ^ 13 + 3464 * p ^ 2 * q ^ 12 + 21740 * p ^ 3 * q ^ 11 + 95220 * p ^ 4 * q ^ 10 + 308300 * p ^ 5 * q ^ 9 + 761600 * p ^ 6 * q ^ 8 + 1459400 * p ^ 7 * q ^ 7 + 2182100 * p ^ 8 * q ^ 6 + 2537000 * p ^ 9 * q ^ 5 + 2261900 * p ^ 10 * q ^ 4 + 1503000 * p ^ 11 * q ^ 3 + 705500 * p ^ 12 * q ^ 2 + 210000 * p ^ 13 * q ^ 1 + 30000 * p ^ 14 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [16, 120, 448, 1180, 2336, 3648, 4700, 5104, 4700, 3648, 2336, 1180, 448, 120, 16] x =
      (fun p q : ℝ => 16 * q ^ 14 + 344 * p ^ 1 * q ^ 13 + 3464 * p ^ 2 * q ^ 12 + 21740 * p ^ 3 * q ^ 11 + 95220 * p ^ 4 * q ^ 10 + 308300 * p ^ 5 * q ^ 9 + 761600 * p ^ 6 * q ^ 8 + 1459400 * p ^ 7 * q ^ 7 + 2182100 * p ^ 8 * q ^ 6 + 2537000 * p ^ 9 * q ^ 5 + 2261900 * p ^ 10 * q ^ 4 + 1503000 * p ^ 11 * q ^ 3 + 705500 * p ^ 12 * q ^ 2 + 210000 * p ^ 13 * q ^ 1 + 30000 * p ^ 14) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [16, 120, 448, 1180, 2336, 3648, 4700, 5104, 4700, 3648, 2336, 1180, 448, 120, 16] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N5_posD1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [16, 120, 448, 1180, 2336, 3648, 4700, 5104, 4700, 3648, 2336, 1180, 448, 120, 16] x := N5_posD1_0 x hu hv

theorem N5_D1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 5 1 x * bBlock 5 1 x - gBlock 5 1 x ^ 2 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 5 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 8) (N5_posD1 x h0 h1)
  rw [← N5_linkD1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))).mp hn

theorem N5_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 5 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 5) ^ 2) = (1 - x) ^ 4 * peval [2, 20, 40, 53, 70, 53, 40, 20, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 5 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (5 - 2 : ℕ) = 3 by norm_num, show (2 * 5 - 2 : ℕ) = 8 by norm_num, show (5 + 2 : ℕ) = 7 by norm_num, show (5 - 2 + 5 : ℕ) = 8 by norm_num, show (2 + 5 : ℕ) = 7 by norm_num, show (5 + 5 : ℕ) = 10 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N5_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 20, 40, 53, 70, 53, 40, 20, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 8 + 36 * p ^ 1 * q ^ 7 + 236 * p ^ 2 * q ^ 6 + 825 * p ^ 3 * q ^ 5 + 1775 * p ^ 4 * q ^ 4 + 2475 * p ^ 5 * q ^ 3 + 2225 * p ^ 6 * q ^ 2 + 1200 * p ^ 7 * q ^ 1 + 300 * p ^ 8 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 20, 40, 53, 70, 53, 40, 20, 2] x =
      (fun p q : ℝ => 2 * q ^ 8 + 36 * p ^ 1 * q ^ 7 + 236 * p ^ 2 * q ^ 6 + 825 * p ^ 3 * q ^ 5 + 1775 * p ^ 4 * q ^ 4 + 2475 * p ^ 5 * q ^ 3 + 2225 * p ^ 6 * q ^ 2 + 1200 * p ^ 7 * q ^ 1 + 300 * p ^ 8) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 20, 40, 53, 70, 53, 40, 20, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N5_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 20, 40, 53, 70, 53, 40, 20, 2] x := N5_posA2_0 x hu hv

theorem N5_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 5 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 5 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N5_posA2 x h0 h1)
  rw [← N5_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N5_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 5 2 x * (2 * (1 - x ^ 5) ^ 2) = (1 - x) ^ 4 * peval [6, 24, 50, 65, 50, 24, 6] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 5 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (5 - 2 : ℕ) = 3 by norm_num, show (2 * 5 - 2 : ℕ) = 8 by norm_num, show (5 + 2 : ℕ) = 7 by norm_num, show (5 - 2 + 5 : ℕ) = 8 by norm_num, show (2 + 5 : ℕ) = 7 by norm_num, show (5 + 5 : ℕ) = 10 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N5_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [6, 24, 50, 65, 50, 24, 6] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 6 * q ^ 6 + 60 * p ^ 1 * q ^ 5 + 260 * p ^ 2 * q ^ 4 + 625 * p ^ 3 * q ^ 3 + 875 * p ^ 4 * q ^ 2 + 675 * p ^ 5 * q ^ 1 + 225 * p ^ 6 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [6, 24, 50, 65, 50, 24, 6] x =
      (fun p q : ℝ => 6 * q ^ 6 + 60 * p ^ 1 * q ^ 5 + 260 * p ^ 2 * q ^ 4 + 625 * p ^ 3 * q ^ 3 + 875 * p ^ 4 * q ^ 2 + 675 * p ^ 5 * q ^ 1 + 225 * p ^ 6) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [6, 24, 50, 65, 50, 24, 6] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N5_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [6, 24, 50, 65, 50, 24, 6] x := N5_posB2_0 x hu hv

theorem N5_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 5 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 5 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N5_posB2 x h0 h1)
  rw [← N5_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N5_linkD2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    (aBlock 5 2 x * bBlock 5 2 x - gBlock 5 2 x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ 5) ^ 4) = (1 - x) ^ 8 * peval [12, 168, 820, 2408, 4992, 7996, 10424, 11360, 10424, 7996, 4992, 2408, 820, 168, 12] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 5 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock bBlock gBlock W cR
  simp only [show (5 - 2 : ℕ) = 3 by norm_num, show (2 * 5 - 2 : ℕ) = 8 by norm_num, show (5 + 2 : ℕ) = 7 by norm_num, show (5 - 2 + 5 : ℕ) = 8 by norm_num, show (2 + 5 : ℕ) = 7 by norm_num, show (5 + 5 : ℕ) = 10 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N5_posD2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 168, 820, 2408, 4992, 7996, 10424, 11360, 10424, 7996, 4992, 2408, 820, 168, 12] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 12 * q ^ 14 + 336 * p ^ 1 * q ^ 13 + 4096 * p ^ 2 * q ^ 12 + 29720 * p ^ 3 * q ^ 11 + 145660 * p ^ 4 * q ^ 10 + 514900 * p ^ 5 * q ^ 9 + 1362500 * p ^ 6 * q ^ 8 + 2755200 * p ^ 7 * q ^ 7 + 4296300 * p ^ 8 * q ^ 6 + 5160500 * p ^ 9 * q ^ 5 + 4717100 * p ^ 10 * q ^ 4 + 3193000 * p ^ 11 * q ^ 3 + 1518000 * p ^ 12 * q ^ 2 + 455000 * p ^ 13 * q ^ 1 + 65000 * p ^ 14 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [12, 168, 820, 2408, 4992, 7996, 10424, 11360, 10424, 7996, 4992, 2408, 820, 168, 12] x =
      (fun p q : ℝ => 12 * q ^ 14 + 336 * p ^ 1 * q ^ 13 + 4096 * p ^ 2 * q ^ 12 + 29720 * p ^ 3 * q ^ 11 + 145660 * p ^ 4 * q ^ 10 + 514900 * p ^ 5 * q ^ 9 + 1362500 * p ^ 6 * q ^ 8 + 2755200 * p ^ 7 * q ^ 7 + 4296300 * p ^ 8 * q ^ 6 + 5160500 * p ^ 9 * q ^ 5 + 4717100 * p ^ 10 * q ^ 4 + 3193000 * p ^ 11 * q ^ 3 + 1518000 * p ^ 12 * q ^ 2 + 455000 * p ^ 13 * q ^ 1 + 65000 * p ^ 14) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [12, 168, 820, 2408, 4992, 7996, 10424, 11360, 10424, 7996, 4992, 2408, 820, 168, 12] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N5_posD2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 168, 820, 2408, 4992, 7996, 10424, 11360, 10424, 7996, 4992, 2408, 820, 168, 12] x := N5_posD2_0 x hu hv

theorem N5_D2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 5 2 x * bBlock 5 2 x - gBlock 5 2 x ^ 2 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 5 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 8) (N5_posD2 x h0 h1)
  rw [← N5_linkD2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))).mp hn

/-- Note Thm 6(a) for `N = 5`: on all of `(0,1)`, every Fourier block has `a_m > 0`,
`b_m > 0` and `a_m b_m − g_m² > 0` (spectrally stable, positive-definite non-rotational Hessian). -/
theorem stable_N5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (m : ℕ) (hm1 : 1 ≤ m) (hm2 : m ≤ 4) :
    0 < aBlock 5 m x ∧ 0 < bBlock 5 m x ∧
      0 < aBlock 5 m x * bBlock 5 m x - gBlock 5 m x ^ 2 := by
  interval_cases m
  · exact ⟨N5_A1_pos x h0 h1, N5_B1_pos x h0 h1, N5_D1_pos x h0 h1⟩
  · exact ⟨N5_A2_pos x h0 h1, N5_B2_pos x h0 h1, N5_D2_pos x h0 h1⟩
  · rw [show aBlock 5 3 x = aBlock 5 2 x from aBlock_symm 5 2 (by norm_num) x, show bBlock 5 3 x = bBlock 5 2 x from bBlock_symm 5 2 (by norm_num) x, show gBlock 5 3 x = -gBlock 5 2 x from gBlock_symm 5 2 (by norm_num) x, neg_sq]
    exact ⟨N5_A2_pos x h0 h1, N5_B2_pos x h0 h1, N5_D2_pos x h0 h1⟩
  · rw [show aBlock 5 4 x = aBlock 5 1 x from aBlock_symm 5 1 (by norm_num) x, show bBlock 5 4 x = bBlock 5 1 x from bBlock_symm 5 1 (by norm_num) x, show gBlock 5 4 x = -gBlock 5 1 x from gBlock_symm 5 1 (by norm_num) x, neg_sq]
    exact ⟨N5_A1_pos x h0 h1, N5_B1_pos x h0 h1, N5_D1_pos x h0 h1⟩

/-! ### N = 6 -/

theorem N6_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 6 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 6) ^ 2) = (1 - x) ^ 4 * peval [5, 24, 59, 106, 161, 190, 161, 106, 59, 24, 5] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 6 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (6 - 1 : ℕ) = 5 by norm_num, show (2 * 6 - 1 : ℕ) = 11 by norm_num, show (6 + 1 : ℕ) = 7 by norm_num, show (6 - 1 + 6 : ℕ) = 11 by norm_num, show (1 + 6 : ℕ) = 7 by norm_num, show (6 + 6 : ℕ) = 12 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N6_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [5, 24, 59, 106, 161, 190, 161, 106, 59, 24, 5] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 5 * q ^ 10 + 74 * p ^ 1 * q ^ 9 + 500 * p ^ 2 * q ^ 8 + 2042 * p ^ 3 * q ^ 7 + 5621 * p ^ 4 * q ^ 6 + 10970 * p ^ 5 * q ^ 5 + 15440 * p ^ 6 * q ^ 4 + 15500 * p ^ 7 * q ^ 3 + 10625 * p ^ 8 * q ^ 2 + 4500 * p ^ 9 * q ^ 1 + 900 * p ^ 10 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [5, 24, 59, 106, 161, 190, 161, 106, 59, 24, 5] x =
      (fun p q : ℝ => 5 * q ^ 10 + 74 * p ^ 1 * q ^ 9 + 500 * p ^ 2 * q ^ 8 + 2042 * p ^ 3 * q ^ 7 + 5621 * p ^ 4 * q ^ 6 + 10970 * p ^ 5 * q ^ 5 + 15440 * p ^ 6 * q ^ 4 + 15500 * p ^ 7 * q ^ 3 + 10625 * p ^ 8 * q ^ 2 + 4500 * p ^ 9 * q ^ 1 + 900 * p ^ 10) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [5, 24, 59, 106, 161, 190, 161, 106, 59, 24, 5] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N6_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [5, 24, 59, 106, 161, 190, 161, 106, 59, 24, 5] x := N6_posA1_0 x hu hv

theorem N6_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 6 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 6 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N6_posA1 x h0 h1)
  rw [← N6_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N6_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 6 1 x * (2 * (1 - x ^ 6) ^ 2) = (1 - x) ^ 4 * peval [5, 14, 26, 40, 55, 40, 26, 14, 5] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 6 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (6 - 1 : ℕ) = 5 by norm_num, show (2 * 6 - 1 : ℕ) = 11 by norm_num, show (6 + 1 : ℕ) = 7 by norm_num, show (6 - 1 + 6 : ℕ) = 11 by norm_num, show (1 + 6 : ℕ) = 7 by norm_num, show (6 + 6 : ℕ) = 12 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N6_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [5, 14, 26, 40, 55, 40, 26, 14, 5] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 5 * q ^ 8 + 54 * p ^ 1 * q ^ 7 + 264 * p ^ 2 * q ^ 6 + 770 * p ^ 3 * q ^ 5 + 1485 * p ^ 4 * q ^ 4 + 1950 * p ^ 5 * q ^ 3 + 1700 * p ^ 6 * q ^ 2 + 900 * p ^ 7 * q ^ 1 + 225 * p ^ 8 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [5, 14, 26, 40, 55, 40, 26, 14, 5] x =
      (fun p q : ℝ => 5 * q ^ 8 + 54 * p ^ 1 * q ^ 7 + 264 * p ^ 2 * q ^ 6 + 770 * p ^ 3 * q ^ 5 + 1485 * p ^ 4 * q ^ 4 + 1950 * p ^ 5 * q ^ 3 + 1700 * p ^ 6 * q ^ 2 + 900 * p ^ 7 * q ^ 1 + 225 * p ^ 8) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [5, 14, 26, 40, 55, 40, 26, 14, 5] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N6_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [5, 14, 26, 40, 55, 40, 26, 14, 5] x := N6_posB1_0 x hu hv

theorem N6_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 6 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 6 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N6_posB1 x h0 h1)
  rw [← N6_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N6_linkD1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    (aBlock 6 1 x * bBlock 6 1 x - gBlock 6 1 x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ 6) ^ 4) = (1 - x) ^ 8 * peval [25, 190, 725, 1964, 4302, 7824, 12050, 16136, 19123, 20222, 19123, 16136, 12050, 7824, 4302, 1964, 725, 190, 25] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 6 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock bBlock gBlock W cR
  simp only [show (6 - 1 : ℕ) = 5 by norm_num, show (2 * 6 - 1 : ℕ) = 11 by norm_num, show (6 + 1 : ℕ) = 7 by norm_num, show (6 - 1 + 6 : ℕ) = 11 by norm_num, show (1 + 6 : ℕ) = 7 by norm_num, show (6 + 6 : ℕ) = 12 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N6_posD1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [25, 190, 725, 1964, 4302, 7824, 12050, 16136, 19123, 20222, 19123, 16136, 12050, 7824, 4302, 1964, 725, 190, 25] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 25 * q ^ 18 + 640 * p ^ 1 * q ^ 17 + 7780 * p ^ 2 * q ^ 16 + 59804 * p ^ 3 * q ^ 15 + 326462 * p ^ 4 * q ^ 14 + 1346672 * p ^ 5 * q ^ 13 + 4358184 * p ^ 6 * q ^ 12 + 11331636 * p ^ 7 * q ^ 11 + 24028647 * p ^ 8 * q ^ 10 + 41914916 * p ^ 9 * q ^ 9 + 60359180 * p ^ 10 * q ^ 8 + 71673860 * p ^ 11 * q ^ 7 + 69760006 * p ^ 12 * q ^ 6 + 54991972 * p ^ 13 * q ^ 5 + 34410496 * p ^ 14 * q ^ 4 + 16533720 * p ^ 15 * q ^ 3 + 5761665 * p ^ 16 * q ^ 2 + 1304100 * p ^ 17 * q ^ 1 + 144900 * p ^ 18 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [25, 190, 725, 1964, 4302, 7824, 12050, 16136, 19123, 20222, 19123, 16136, 12050, 7824, 4302, 1964, 725, 190, 25] x =
      (fun p q : ℝ => 25 * q ^ 18 + 640 * p ^ 1 * q ^ 17 + 7780 * p ^ 2 * q ^ 16 + 59804 * p ^ 3 * q ^ 15 + 326462 * p ^ 4 * q ^ 14 + 1346672 * p ^ 5 * q ^ 13 + 4358184 * p ^ 6 * q ^ 12 + 11331636 * p ^ 7 * q ^ 11 + 24028647 * p ^ 8 * q ^ 10 + 41914916 * p ^ 9 * q ^ 9 + 60359180 * p ^ 10 * q ^ 8 + 71673860 * p ^ 11 * q ^ 7 + 69760006 * p ^ 12 * q ^ 6 + 54991972 * p ^ 13 * q ^ 5 + 34410496 * p ^ 14 * q ^ 4 + 16533720 * p ^ 15 * q ^ 3 + 5761665 * p ^ 16 * q ^ 2 + 1304100 * p ^ 17 * q ^ 1 + 144900 * p ^ 18) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [25, 190, 725, 1964, 4302, 7824, 12050, 16136, 19123, 20222, 19123, 16136, 12050, 7824, 4302, 1964, 725, 190, 25] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N6_posD1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [25, 190, 725, 1964, 4302, 7824, 12050, 16136, 19123, 20222, 19123, 16136, 12050, 7824, 4302, 1964, 725, 190, 25] x := N6_posD1_0 x hu hv

theorem N6_D1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 6 1 x * bBlock 6 1 x - gBlock 6 1 x ^ 2 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 6 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 8) (N6_posD1 x h0 h1)
  rw [← N6_linkD1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))).mp hn

theorem N6_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 6 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 6) ^ 2) = (1 - x) ^ 4 * peval [2, 24, 50, 88, 110, 136, 110, 88, 50, 24, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 6 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (6 - 2 : ℕ) = 4 by norm_num, show (2 * 6 - 2 : ℕ) = 10 by norm_num, show (6 + 2 : ℕ) = 8 by norm_num, show (6 - 2 + 6 : ℕ) = 10 by norm_num, show (2 + 6 : ℕ) = 8 by norm_num, show (6 + 6 : ℕ) = 12 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N6_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 24, 50, 88, 110, 136, 110, 88, 50, 24, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 10 + 44 * p ^ 1 * q ^ 9 + 356 * p ^ 2 * q ^ 8 + 1592 * p ^ 3 * q ^ 7 + 4562 * p ^ 4 * q ^ 6 + 8972 * p ^ 5 * q ^ 5 + 12464 * p ^ 6 * q ^ 4 + 12224 * p ^ 7 * q ^ 3 + 8186 * p ^ 8 * q ^ 2 + 3420 * p ^ 9 * q ^ 1 + 684 * p ^ 10 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 24, 50, 88, 110, 136, 110, 88, 50, 24, 2] x =
      (fun p q : ℝ => 2 * q ^ 10 + 44 * p ^ 1 * q ^ 9 + 356 * p ^ 2 * q ^ 8 + 1592 * p ^ 3 * q ^ 7 + 4562 * p ^ 4 * q ^ 6 + 8972 * p ^ 5 * q ^ 5 + 12464 * p ^ 6 * q ^ 4 + 12224 * p ^ 7 * q ^ 3 + 8186 * p ^ 8 * q ^ 2 + 3420 * p ^ 9 * q ^ 1 + 684 * p ^ 10) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 24, 50, 88, 110, 136, 110, 88, 50, 24, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N6_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 24, 50, 88, 110, 136, 110, 88, 50, 24, 2] x := N6_posA2_0 x hu hv

theorem N6_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 6 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 6 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N6_posA2 x h0 h1)
  rw [← N6_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N6_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 6 2 x * (2 * (1 - x ^ 6) ^ 2) = (1 - x) ^ 4 * peval [8, 32, 68, 112, 136, 112, 68, 32, 8] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 6 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (6 - 2 : ℕ) = 4 by norm_num, show (2 * 6 - 2 : ℕ) = 10 by norm_num, show (6 + 2 : ℕ) = 8 by norm_num, show (6 - 2 + 6 : ℕ) = 10 by norm_num, show (2 + 6 : ℕ) = 8 by norm_num, show (6 + 6 : ℕ) = 12 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N6_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [8, 32, 68, 112, 136, 112, 68, 32, 8] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 8 * q ^ 8 + 96 * p ^ 1 * q ^ 7 + 516 * p ^ 2 * q ^ 6 + 1640 * p ^ 3 * q ^ 5 + 3396 * p ^ 4 * q ^ 4 + 4704 * p ^ 5 * q ^ 3 + 4256 * p ^ 6 * q ^ 2 + 2304 * p ^ 7 * q ^ 1 + 576 * p ^ 8 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [8, 32, 68, 112, 136, 112, 68, 32, 8] x =
      (fun p q : ℝ => 8 * q ^ 8 + 96 * p ^ 1 * q ^ 7 + 516 * p ^ 2 * q ^ 6 + 1640 * p ^ 3 * q ^ 5 + 3396 * p ^ 4 * q ^ 4 + 4704 * p ^ 5 * q ^ 3 + 4256 * p ^ 6 * q ^ 2 + 2304 * p ^ 7 * q ^ 1 + 576 * p ^ 8) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [8, 32, 68, 112, 136, 112, 68, 32, 8] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N6_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [8, 32, 68, 112, 136, 112, 68, 32, 8] x := N6_posB2_0 x hu hv

theorem N6_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 6 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 6 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N6_posB2 x h0 h1)
  rw [← N6_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N6_linkD2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    (aBlock 6 2 x * bBlock 6 2 x - gBlock 6 2 x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ 6) ^ 4) = (1 - x) ^ 8 * peval [16, 256, 1304, 4160, 9912, 18816, 29744, 40448, 48304, 51200, 48304, 40448, 29744, 18816, 9912, 4160, 1304, 256, 16] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 6 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock bBlock gBlock W cR
  simp only [show (6 - 2 : ℕ) = 4 by norm_num, show (2 * 6 - 2 : ℕ) = 10 by norm_num, show (6 + 2 : ℕ) = 8 by norm_num, show (6 - 2 + 6 : ℕ) = 10 by norm_num, show (2 + 6 : ℕ) = 8 by norm_num, show (6 + 6 : ℕ) = 12 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N6_posD2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [16, 256, 1304, 4160, 9912, 18816, 29744, 40448, 48304, 51200, 48304, 40448, 29744, 18816, 9912, 4160, 1304, 256, 16] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 16 * q ^ 18 + 544 * p ^ 1 * q ^ 17 + 8104 * p ^ 2 * q ^ 16 + 72896 * p ^ 3 * q ^ 15 + 451832 * p ^ 4 * q ^ 14 + 2070992 * p ^ 5 * q ^ 13 + 7323576 * p ^ 6 * q ^ 12 + 20524704 * p ^ 7 * q ^ 11 + 46373352 * p ^ 8 * q ^ 10 + 85339664 * p ^ 9 * q ^ 9 + 128537480 * p ^ 10 * q ^ 8 + 158453120 * p ^ 11 * q ^ 7 + 159062968 * p ^ 12 * q ^ 6 + 128589328 * p ^ 13 * q ^ 5 + 82098544 * p ^ 14 * q ^ 4 + 40061952 * p ^ 15 * q ^ 3 + 14114304 * p ^ 16 * q ^ 2 + 3214080 * p ^ 17 * q ^ 1 + 357120 * p ^ 18 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [16, 256, 1304, 4160, 9912, 18816, 29744, 40448, 48304, 51200, 48304, 40448, 29744, 18816, 9912, 4160, 1304, 256, 16] x =
      (fun p q : ℝ => 16 * q ^ 18 + 544 * p ^ 1 * q ^ 17 + 8104 * p ^ 2 * q ^ 16 + 72896 * p ^ 3 * q ^ 15 + 451832 * p ^ 4 * q ^ 14 + 2070992 * p ^ 5 * q ^ 13 + 7323576 * p ^ 6 * q ^ 12 + 20524704 * p ^ 7 * q ^ 11 + 46373352 * p ^ 8 * q ^ 10 + 85339664 * p ^ 9 * q ^ 9 + 128537480 * p ^ 10 * q ^ 8 + 158453120 * p ^ 11 * q ^ 7 + 159062968 * p ^ 12 * q ^ 6 + 128589328 * p ^ 13 * q ^ 5 + 82098544 * p ^ 14 * q ^ 4 + 40061952 * p ^ 15 * q ^ 3 + 14114304 * p ^ 16 * q ^ 2 + 3214080 * p ^ 17 * q ^ 1 + 357120 * p ^ 18) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [16, 256, 1304, 4160, 9912, 18816, 29744, 40448, 48304, 51200, 48304, 40448, 29744, 18816, 9912, 4160, 1304, 256, 16] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N6_posD2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [16, 256, 1304, 4160, 9912, 18816, 29744, 40448, 48304, 51200, 48304, 40448, 29744, 18816, 9912, 4160, 1304, 256, 16] x := N6_posD2_0 x hu hv

theorem N6_D2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 6 2 x * bBlock 6 2 x - gBlock 6 2 x ^ 2 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 6 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 8) (N6_posD2 x h0 h1)
  rw [← N6_linkD2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))).mp hn

theorem N6_linkA3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 6 3 x * (2 * (1 - x) ^ 2 * (1 - x ^ 6) ^ 2) = (1 - x) ^ 4 * peval [1, 22, 59, 72, 93, 118, 93, 72, 59, 22, 1] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 6 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (6 - 3 : ℕ) = 3 by norm_num, show (2 * 6 - 3 : ℕ) = 9 by norm_num, show (6 + 3 : ℕ) = 9 by norm_num, show (6 - 3 + 6 : ℕ) = 9 by norm_num, show (3 + 6 : ℕ) = 9 by norm_num, show (6 + 6 : ℕ) = 12 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N6_posA3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1, 22, 59, 72, 93, 118, 93, 72, 59, 22, 1] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 1 * q ^ 10 + 32 * p ^ 1 * q ^ 9 + 302 * p ^ 2 * q ^ 8 + 1456 * p ^ 3 * q ^ 7 + 4307 * p ^ 4 * q ^ 6 + 8516 * p ^ 5 * q ^ 5 + 11710 * p ^ 6 * q ^ 4 + 11276 * p ^ 7 * q ^ 3 + 7409 * p ^ 8 * q ^ 2 + 3060 * p ^ 9 * q ^ 1 + 612 * p ^ 10 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [1, 22, 59, 72, 93, 118, 93, 72, 59, 22, 1] x =
      (fun p q : ℝ => 1 * q ^ 10 + 32 * p ^ 1 * q ^ 9 + 302 * p ^ 2 * q ^ 8 + 1456 * p ^ 3 * q ^ 7 + 4307 * p ^ 4 * q ^ 6 + 8516 * p ^ 5 * q ^ 5 + 11710 * p ^ 6 * q ^ 4 + 11276 * p ^ 7 * q ^ 3 + 7409 * p ^ 8 * q ^ 2 + 3060 * p ^ 9 * q ^ 1 + 612 * p ^ 10) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [1, 22, 59, 72, 93, 118, 93, 72, 59, 22, 1] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N6_posA3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1, 22, 59, 72, 93, 118, 93, 72, 59, 22, 1] x := N6_posA3_0 x hu hv

theorem N6_A3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 6 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 6 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N6_posA3 x h0 h1)
  rw [← N6_linkA3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N6_linkB3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 6 3 x * (2 * (1 - x ^ 6) ^ 2) = (1 - x) ^ 4 * peval [9, 36, 90, 144, 171, 144, 90, 36, 9] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 6 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (6 - 3 : ℕ) = 3 by norm_num, show (2 * 6 - 3 : ℕ) = 9 by norm_num, show (6 + 3 : ℕ) = 9 by norm_num, show (6 - 3 + 6 : ℕ) = 9 by norm_num, show (3 + 6 : ℕ) = 9 by norm_num, show (6 + 6 : ℕ) = 12 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N6_posB3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [9, 36, 90, 144, 171, 144, 90, 36, 9] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 9 * q ^ 8 + 108 * p ^ 1 * q ^ 7 + 594 * p ^ 2 * q ^ 6 + 1944 * p ^ 3 * q ^ 5 + 4131 * p ^ 4 * q ^ 4 + 5832 * p ^ 5 * q ^ 3 + 5346 * p ^ 6 * q ^ 2 + 2916 * p ^ 7 * q ^ 1 + 729 * p ^ 8 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [9, 36, 90, 144, 171, 144, 90, 36, 9] x =
      (fun p q : ℝ => 9 * q ^ 8 + 108 * p ^ 1 * q ^ 7 + 594 * p ^ 2 * q ^ 6 + 1944 * p ^ 3 * q ^ 5 + 4131 * p ^ 4 * q ^ 4 + 5832 * p ^ 5 * q ^ 3 + 5346 * p ^ 6 * q ^ 2 + 2916 * p ^ 7 * q ^ 1 + 729 * p ^ 8) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [9, 36, 90, 144, 171, 144, 90, 36, 9] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N6_posB3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [9, 36, 90, 144, 171, 144, 90, 36, 9] x := N6_posB3_0 x hu hv

theorem N6_B3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 6 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 6 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N6_posB3 x h0 h1)
  rw [← N6_linkB3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N6_linkD3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    (aBlock 6 3 x * bBlock 6 3 x - gBlock 6 3 x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ 6) ^ 4) = (1 - x) ^ 8 * peval [9, 234, 1413, 4896, 12078, 23292, 37170, 50832, 60867, 64566, 60867, 50832, 37170, 23292, 12078, 4896, 1413, 234, 9] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 6 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock bBlock gBlock W cR
  simp only [show (6 - 3 : ℕ) = 3 by norm_num, show (2 * 6 - 3 : ℕ) = 9 by norm_num, show (6 + 3 : ℕ) = 9 by norm_num, show (6 - 3 + 6 : ℕ) = 9 by norm_num, show (3 + 6 : ℕ) = 9 by norm_num, show (6 + 6 : ℕ) = 12 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N6_posD3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [9, 234, 1413, 4896, 12078, 23292, 37170, 50832, 60867, 64566, 60867, 50832, 37170, 23292, 12078, 4896, 1413, 234, 9] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 9 * q ^ 18 + 396 * p ^ 1 * q ^ 17 + 6768 * p ^ 2 * q ^ 16 + 66672 * p ^ 3 * q ^ 15 + 441738 * p ^ 4 * q ^ 14 + 2131776 * p ^ 5 * q ^ 13 + 7853472 * p ^ 6 * q ^ 12 + 22747464 * p ^ 7 * q ^ 11 + 52787475 * p ^ 8 * q ^ 10 + 99275652 * p ^ 9 * q ^ 9 + 152187876 * p ^ 10 * q ^ 8 + 190307880 * p ^ 11 * q ^ 7 + 193254498 * p ^ 12 * q ^ 6 + 157675896 * p ^ 13 * q ^ 5 + 101400012 * p ^ 14 * q ^ 4 + 49752792 * p ^ 15 * q ^ 3 + 17595873 * p ^ 16 * q ^ 2 + 4015332 * p ^ 17 * q ^ 1 + 446148 * p ^ 18 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [9, 234, 1413, 4896, 12078, 23292, 37170, 50832, 60867, 64566, 60867, 50832, 37170, 23292, 12078, 4896, 1413, 234, 9] x =
      (fun p q : ℝ => 9 * q ^ 18 + 396 * p ^ 1 * q ^ 17 + 6768 * p ^ 2 * q ^ 16 + 66672 * p ^ 3 * q ^ 15 + 441738 * p ^ 4 * q ^ 14 + 2131776 * p ^ 5 * q ^ 13 + 7853472 * p ^ 6 * q ^ 12 + 22747464 * p ^ 7 * q ^ 11 + 52787475 * p ^ 8 * q ^ 10 + 99275652 * p ^ 9 * q ^ 9 + 152187876 * p ^ 10 * q ^ 8 + 190307880 * p ^ 11 * q ^ 7 + 193254498 * p ^ 12 * q ^ 6 + 157675896 * p ^ 13 * q ^ 5 + 101400012 * p ^ 14 * q ^ 4 + 49752792 * p ^ 15 * q ^ 3 + 17595873 * p ^ 16 * q ^ 2 + 4015332 * p ^ 17 * q ^ 1 + 446148 * p ^ 18) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [9, 234, 1413, 4896, 12078, 23292, 37170, 50832, 60867, 64566, 60867, 50832, 37170, 23292, 12078, 4896, 1413, 234, 9] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N6_posD3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [9, 234, 1413, 4896, 12078, 23292, 37170, 50832, 60867, 64566, 60867, 50832, 37170, 23292, 12078, 4896, 1413, 234, 9] x := N6_posD3_0 x hu hv

theorem N6_D3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 6 3 x * bBlock 6 3 x - gBlock 6 3 x ^ 2 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 6 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 8) (N6_posD3 x h0 h1)
  rw [← N6_linkD3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))).mp hn

/-- Note Thm 6(a) for `N = 6`: on all of `(0,1)`, every Fourier block has `a_m > 0`,
`b_m > 0` and `a_m b_m − g_m² > 0` (spectrally stable, positive-definite non-rotational Hessian). -/
theorem stable_N6 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (m : ℕ) (hm1 : 1 ≤ m) (hm2 : m ≤ 5) :
    0 < aBlock 6 m x ∧ 0 < bBlock 6 m x ∧
      0 < aBlock 6 m x * bBlock 6 m x - gBlock 6 m x ^ 2 := by
  interval_cases m
  · exact ⟨N6_A1_pos x h0 h1, N6_B1_pos x h0 h1, N6_D1_pos x h0 h1⟩
  · exact ⟨N6_A2_pos x h0 h1, N6_B2_pos x h0 h1, N6_D2_pos x h0 h1⟩
  · exact ⟨N6_A3_pos x h0 h1, N6_B3_pos x h0 h1, N6_D3_pos x h0 h1⟩
  · rw [show aBlock 6 4 x = aBlock 6 2 x from aBlock_symm 6 2 (by norm_num) x, show bBlock 6 4 x = bBlock 6 2 x from bBlock_symm 6 2 (by norm_num) x, show gBlock 6 4 x = -gBlock 6 2 x from gBlock_symm 6 2 (by norm_num) x, neg_sq]
    exact ⟨N6_A2_pos x h0 h1, N6_B2_pos x h0 h1, N6_D2_pos x h0 h1⟩
  · rw [show aBlock 6 5 x = aBlock 6 1 x from aBlock_symm 6 1 (by norm_num) x, show bBlock 6 5 x = bBlock 6 1 x from bBlock_symm 6 1 (by norm_num) x, show gBlock 6 5 x = -gBlock 6 1 x from gBlock_symm 6 1 (by norm_num) x, neg_sq]
    exact ⟨N6_A1_pos x h0 h1, N6_B1_pos x h0 h1, N6_D1_pos x h0 h1⟩

/-! ### N = 7 -/

theorem N7_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 7 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 7) ^ 2) = (1 - x) ^ 4 * peval [6, 29, 72, 131, 202, 281, 322, 281, 202, 131, 72, 29, 6] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 7 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (7 - 1 : ℕ) = 6 by norm_num, show (2 * 7 - 1 : ℕ) = 13 by norm_num, show (7 + 1 : ℕ) = 8 by norm_num, show (7 - 1 + 7 : ℕ) = 13 by norm_num, show (1 + 7 : ℕ) = 8 by norm_num, show (7 + 7 : ℕ) = 14 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N7_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [6, 29, 72, 131, 202, 281, 322, 281, 202, 131, 72, 29, 6] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 6 * q ^ 12 + 101 * p ^ 1 * q ^ 11 + 787 * p ^ 2 * q ^ 10 + 3766 * p ^ 3 * q ^ 9 + 12376 * p ^ 4 * q ^ 8 + 29575 * p ^ 5 * q ^ 7 + 53011 * p ^ 6 * q ^ 6 + 72226 * p ^ 7 * q ^ 5 + 74578 * p ^ 8 * q ^ 4 + 57085 * p ^ 9 * q ^ 3 + 30821 * p ^ 10 * q ^ 2 + 10584 * p ^ 11 * q ^ 1 + 1764 * p ^ 12 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [6, 29, 72, 131, 202, 281, 322, 281, 202, 131, 72, 29, 6] x =
      (fun p q : ℝ => 6 * q ^ 12 + 101 * p ^ 1 * q ^ 11 + 787 * p ^ 2 * q ^ 10 + 3766 * p ^ 3 * q ^ 9 + 12376 * p ^ 4 * q ^ 8 + 29575 * p ^ 5 * q ^ 7 + 53011 * p ^ 6 * q ^ 6 + 72226 * p ^ 7 * q ^ 5 + 74578 * p ^ 8 * q ^ 4 + 57085 * p ^ 9 * q ^ 3 + 30821 * p ^ 10 * q ^ 2 + 10584 * p ^ 11 * q ^ 1 + 1764 * p ^ 12) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [6, 29, 72, 131, 202, 281, 322, 281, 202, 131, 72, 29, 6] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N7_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [6, 29, 72, 131, 202, 281, 322, 281, 202, 131, 72, 29, 6] x := N7_posA1_0 x hu hv

theorem N7_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 7 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 7 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N7_posA1 x h0 h1)
  rw [← N7_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N7_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 7 1 x * (2 * (1 - x ^ 7) ^ 2) = (1 - x) ^ 4 * peval [6, 17, 32, 50, 70, 91, 70, 50, 32, 17, 6] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 7 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (7 - 1 : ℕ) = 6 by norm_num, show (2 * 7 - 1 : ℕ) = 13 by norm_num, show (7 + 1 : ℕ) = 8 by norm_num, show (7 - 1 + 7 : ℕ) = 13 by norm_num, show (1 + 7 : ℕ) = 8 by norm_num, show (7 + 7 : ℕ) = 14 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N7_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [6, 17, 32, 50, 70, 91, 70, 50, 32, 17, 6] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 6 * q ^ 10 + 77 * p ^ 1 * q ^ 9 + 455 * p ^ 2 * q ^ 8 + 1638 * p ^ 3 * q ^ 7 + 4004 * p ^ 4 * q ^ 6 + 7007 * p ^ 5 * q ^ 5 + 8967 * p ^ 6 * q ^ 4 + 8330 * p ^ 7 * q ^ 3 + 5390 * p ^ 8 * q ^ 2 + 2205 * p ^ 9 * q ^ 1 + 441 * p ^ 10 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [6, 17, 32, 50, 70, 91, 70, 50, 32, 17, 6] x =
      (fun p q : ℝ => 6 * q ^ 10 + 77 * p ^ 1 * q ^ 9 + 455 * p ^ 2 * q ^ 8 + 1638 * p ^ 3 * q ^ 7 + 4004 * p ^ 4 * q ^ 6 + 7007 * p ^ 5 * q ^ 5 + 8967 * p ^ 6 * q ^ 4 + 8330 * p ^ 7 * q ^ 3 + 5390 * p ^ 8 * q ^ 2 + 2205 * p ^ 9 * q ^ 1 + 441 * p ^ 10) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [6, 17, 32, 50, 70, 91, 70, 50, 32, 17, 6] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N7_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [6, 17, 32, 50, 70, 91, 70, 50, 32, 17, 6] x := N7_posB1_0 x hu hv

theorem N7_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 7 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 7 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N7_posB1 x h0 h1)
  rw [← N7_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N7_linkD1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    (aBlock 7 1 x * bBlock 7 1 x - gBlock 7 1 x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ 7) ^ 4) = (1 - x) ^ 8 * peval [36, 276, 1068, 2944, 6584, 12744, 21648, 32544, 43976, 54116, 61144, 63664, 61144, 54116, 43976, 32544, 21648, 12744, 6584, 2944, 1068, 276, 36] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 7 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock bBlock gBlock W cR
  simp only [show (7 - 1 : ℕ) = 6 by norm_num, show (2 * 7 - 1 : ℕ) = 13 by norm_num, show (7 + 1 : ℕ) = 8 by norm_num, show (7 - 1 + 7 : ℕ) = 13 by norm_num, show (1 + 7 : ℕ) = 8 by norm_num, show (7 + 7 : ℕ) = 14 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N7_posD1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [36, 276, 1068, 2944, 6584, 12744, 21648, 32544, 43976, 54116, 61144, 63664, 61144, 54116, 43976, 32544, 21648, 12744, 6584, 2944, 1068, 276, 36] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 36 * q ^ 22 + 1068 * p ^ 1 * q ^ 21 + 15180 * p ^ 2 * q ^ 20 + 137704 * p ^ 3 * q ^ 19 + 895860 * p ^ 4 * q ^ 18 + 4452084 * p ^ 5 * q ^ 17 + 17575236 * p ^ 6 * q ^ 16 + 56570304 * p ^ 7 * q ^ 15 + 151176368 * p ^ 8 * q ^ 14 + 339690540 * p ^ 9 * q ^ 13 + 647340764 * p ^ 10 * q ^ 12 + 1051895544 * p ^ 11 * q ^ 11 + 1461238016 * p ^ 12 * q ^ 10 + 1735182120 * p ^ 13 * q ^ 9 + 1756424600 * p ^ 14 * q ^ 8 + 1506708448 * p ^ 15 * q ^ 7 + 1084819820 * p ^ 16 * q ^ 6 + 645956808 * p ^ 17 * q ^ 5 + 311059840 * p ^ 18 * q ^ 4 + 116976720 * p ^ 19 * q ^ 3 + 32403896 * p ^ 20 * q ^ 2 + 5916064 * p ^ 21 * q ^ 1 + 537824 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [36, 276, 1068, 2944, 6584, 12744, 21648, 32544, 43976, 54116, 61144, 63664, 61144, 54116, 43976, 32544, 21648, 12744, 6584, 2944, 1068, 276, 36] x =
      (fun p q : ℝ => 36 * q ^ 22 + 1068 * p ^ 1 * q ^ 21 + 15180 * p ^ 2 * q ^ 20 + 137704 * p ^ 3 * q ^ 19 + 895860 * p ^ 4 * q ^ 18 + 4452084 * p ^ 5 * q ^ 17 + 17575236 * p ^ 6 * q ^ 16 + 56570304 * p ^ 7 * q ^ 15 + 151176368 * p ^ 8 * q ^ 14 + 339690540 * p ^ 9 * q ^ 13 + 647340764 * p ^ 10 * q ^ 12 + 1051895544 * p ^ 11 * q ^ 11 + 1461238016 * p ^ 12 * q ^ 10 + 1735182120 * p ^ 13 * q ^ 9 + 1756424600 * p ^ 14 * q ^ 8 + 1506708448 * p ^ 15 * q ^ 7 + 1084819820 * p ^ 16 * q ^ 6 + 645956808 * p ^ 17 * q ^ 5 + 311059840 * p ^ 18 * q ^ 4 + 116976720 * p ^ 19 * q ^ 3 + 32403896 * p ^ 20 * q ^ 2 + 5916064 * p ^ 21 * q ^ 1 + 537824 * p ^ 22) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [36, 276, 1068, 2944, 6584, 12744, 21648, 32544, 43976, 54116, 61144, 63664, 61144, 54116, 43976, 32544, 21648, 12744, 6584, 2944, 1068, 276, 36] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N7_posD1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [36, 276, 1068, 2944, 6584, 12744, 21648, 32544, 43976, 54116, 61144, 63664, 61144, 54116, 43976, 32544, 21648, 12744, 6584, 2944, 1068, 276, 36] x := N7_posD1_0 x hu hv

theorem N7_D1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 7 1 x * bBlock 7 1 x - gBlock 7 1 x ^ 2 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 7 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 8) (N7_posD1 x h0 h1)
  rw [← N7_linkD1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))).mp hn

theorem N7_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 7 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 7) ^ 2) = (1 - x) ^ 4 * peval [2, 28, 60, 108, 168, 201, 238, 201, 168, 108, 60, 28, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 7 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (7 - 2 : ℕ) = 5 by norm_num, show (2 * 7 - 2 : ℕ) = 12 by norm_num, show (7 + 2 : ℕ) = 9 by norm_num, show (7 - 2 + 7 : ℕ) = 12 by norm_num, show (2 + 7 : ℕ) = 9 by norm_num, show (7 + 7 : ℕ) = 14 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N7_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 28, 60, 108, 168, 201, 238, 201, 168, 108, 60, 28, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 12 + 52 * p ^ 1 * q ^ 11 + 500 * p ^ 2 * q ^ 10 + 2688 * p ^ 3 * q ^ 9 + 9450 * p ^ 4 * q ^ 8 + 23457 * p ^ 5 * q ^ 7 + 42805 * p ^ 6 * q ^ 6 + 58506 * p ^ 7 * q ^ 5 + 59976 * p ^ 8 * q ^ 4 + 45325 * p ^ 9 * q ^ 3 + 24157 * p ^ 10 * q ^ 2 + 8232 * p ^ 11 * q ^ 1 + 1372 * p ^ 12 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 28, 60, 108, 168, 201, 238, 201, 168, 108, 60, 28, 2] x =
      (fun p q : ℝ => 2 * q ^ 12 + 52 * p ^ 1 * q ^ 11 + 500 * p ^ 2 * q ^ 10 + 2688 * p ^ 3 * q ^ 9 + 9450 * p ^ 4 * q ^ 8 + 23457 * p ^ 5 * q ^ 7 + 42805 * p ^ 6 * q ^ 6 + 58506 * p ^ 7 * q ^ 5 + 59976 * p ^ 8 * q ^ 4 + 45325 * p ^ 9 * q ^ 3 + 24157 * p ^ 10 * q ^ 2 + 8232 * p ^ 11 * q ^ 1 + 1372 * p ^ 12) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 28, 60, 108, 168, 201, 238, 201, 168, 108, 60, 28, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N7_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 28, 60, 108, 168, 201, 238, 201, 168, 108, 60, 28, 2] x := N7_posA2_0 x hu hv

theorem N7_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 7 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 7 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N7_posA2 x h0 h1)
  rw [← N7_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N7_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 7 2 x * (2 * (1 - x ^ 7) ^ 2) = (1 - x) ^ 4 * peval [10, 40, 86, 144, 210, 245, 210, 144, 86, 40, 10] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 7 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (7 - 2 : ℕ) = 5 by norm_num, show (2 * 7 - 2 : ℕ) = 12 by norm_num, show (7 + 2 : ℕ) = 9 by norm_num, show (7 - 2 + 7 : ℕ) = 12 by norm_num, show (2 + 7 : ℕ) = 9 by norm_num, show (7 + 7 : ℕ) = 14 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N7_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [10, 40, 86, 144, 210, 245, 210, 144, 86, 40, 10] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 10 * q ^ 10 + 140 * p ^ 1 * q ^ 9 + 896 * p ^ 2 * q ^ 8 + 3472 * p ^ 3 * q ^ 7 + 9086 * p ^ 4 * q ^ 6 + 16905 * p ^ 5 * q ^ 5 + 22785 * p ^ 6 * q ^ 4 + 22050 * p ^ 7 * q ^ 3 + 14700 * p ^ 8 * q ^ 2 + 6125 * p ^ 9 * q ^ 1 + 1225 * p ^ 10 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [10, 40, 86, 144, 210, 245, 210, 144, 86, 40, 10] x =
      (fun p q : ℝ => 10 * q ^ 10 + 140 * p ^ 1 * q ^ 9 + 896 * p ^ 2 * q ^ 8 + 3472 * p ^ 3 * q ^ 7 + 9086 * p ^ 4 * q ^ 6 + 16905 * p ^ 5 * q ^ 5 + 22785 * p ^ 6 * q ^ 4 + 22050 * p ^ 7 * q ^ 3 + 14700 * p ^ 8 * q ^ 2 + 6125 * p ^ 9 * q ^ 1 + 1225 * p ^ 10) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [10, 40, 86, 144, 210, 245, 210, 144, 86, 40, 10] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N7_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [10, 40, 86, 144, 210, 245, 210, 144, 86, 40, 10] x := N7_posB2_0 x hu hv

theorem N7_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 7 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 7 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N7_posB2 x h0 h1)
  rw [← N7_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N7_linkD2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    (aBlock 7 2 x * bBlock 7 2 x - gBlock 7 2 x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ 7) ^ 4) = (1 - x) ^ 8 * peval [20, 360, 1892, 6176, 15416, 31852, 56184, 86560, 118860, 147724, 167772, 174968, 167772, 147724, 118860, 86560, 56184, 31852, 15416, 6176, 1892, 360, 20] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 7 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock bBlock gBlock W cR
  simp only [show (7 - 2 : ℕ) = 5 by norm_num, show (2 * 7 - 2 : ℕ) = 12 by norm_num, show (7 + 2 : ℕ) = 9 by norm_num, show (7 - 2 + 7 : ℕ) = 12 by norm_num, show (2 + 7 : ℕ) = 9 by norm_num, show (7 + 7 : ℕ) = 14 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N7_posD2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [20, 360, 1892, 6176, 15416, 31852, 56184, 86560, 118860, 147724, 167772, 174968, 167772, 147724, 118860, 86560, 56184, 31852, 15416, 6176, 1892, 360, 20] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 20 * q ^ 22 + 800 * p ^ 1 * q ^ 21 + 14072 * p ^ 2 * q ^ 20 + 150416 * p ^ 3 * q ^ 19 + 1117340 * p ^ 4 * q ^ 18 + 6203596 * p ^ 5 * q ^ 17 + 26925500 * p ^ 6 * q ^ 16 + 94114496 * p ^ 7 * q ^ 15 + 270396308 * p ^ 8 * q ^ 14 + 647695524 * p ^ 9 * q ^ 13 + 1306106172 * p ^ 10 * q ^ 12 + 2231002144 * p ^ 11 * q ^ 11 + 3238281032 * p ^ 12 * q ^ 10 + 3995834752 * p ^ 13 * q ^ 9 + 4181750944 * p ^ 14 * q ^ 8 + 3691475760 * p ^ 15 * q ^ 7 + 2723384328 * p ^ 16 * q ^ 6 + 1655060064 * p ^ 17 * q ^ 5 + 810396496 * p ^ 18 * q ^ 4 + 308768600 * p ^ 19 * q ^ 3 + 86339960 * p ^ 20 * q ^ 2 + 15846600 * p ^ 21 * q ^ 1 + 1440600 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [20, 360, 1892, 6176, 15416, 31852, 56184, 86560, 118860, 147724, 167772, 174968, 167772, 147724, 118860, 86560, 56184, 31852, 15416, 6176, 1892, 360, 20] x =
      (fun p q : ℝ => 20 * q ^ 22 + 800 * p ^ 1 * q ^ 21 + 14072 * p ^ 2 * q ^ 20 + 150416 * p ^ 3 * q ^ 19 + 1117340 * p ^ 4 * q ^ 18 + 6203596 * p ^ 5 * q ^ 17 + 26925500 * p ^ 6 * q ^ 16 + 94114496 * p ^ 7 * q ^ 15 + 270396308 * p ^ 8 * q ^ 14 + 647695524 * p ^ 9 * q ^ 13 + 1306106172 * p ^ 10 * q ^ 12 + 2231002144 * p ^ 11 * q ^ 11 + 3238281032 * p ^ 12 * q ^ 10 + 3995834752 * p ^ 13 * q ^ 9 + 4181750944 * p ^ 14 * q ^ 8 + 3691475760 * p ^ 15 * q ^ 7 + 2723384328 * p ^ 16 * q ^ 6 + 1655060064 * p ^ 17 * q ^ 5 + 810396496 * p ^ 18 * q ^ 4 + 308768600 * p ^ 19 * q ^ 3 + 86339960 * p ^ 20 * q ^ 2 + 15846600 * p ^ 21 * q ^ 1 + 1440600 * p ^ 22) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [20, 360, 1892, 6176, 15416, 31852, 56184, 86560, 118860, 147724, 167772, 174968, 167772, 147724, 118860, 86560, 56184, 31852, 15416, 6176, 1892, 360, 20] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N7_posD2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [20, 360, 1892, 6176, 15416, 31852, 56184, 86560, 118860, 147724, 167772, 174968, 167772, 147724, 118860, 86560, 56184, 31852, 15416, 6176, 1892, 360, 20] x := N7_posD2_0 x hu hv

theorem N7_D2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 7 2 x * bBlock 7 2 x - gBlock 7 2 x ^ 2 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 7 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 8) (N7_posD2 x h0 h1)
  rw [← N7_linkD2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))).mp hn

theorem N7_linkA3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 7 3 x * (2 * (1 - x) ^ 2 * (1 - x ^ 7) ^ 2) = x ^ 1 * (1 - x) ^ 4 * peval [24, 68, 107, 130, 161, 196, 161, 130, 107, 68, 24] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 7 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (7 - 3 : ℕ) = 4 by norm_num, show (2 * 7 - 3 : ℕ) = 11 by norm_num, show (7 + 3 : ℕ) = 10 by norm_num, show (7 - 3 + 7 : ℕ) = 11 by norm_num, show (3 + 7 : ℕ) = 10 by norm_num, show (7 + 7 : ℕ) = 14 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N7_posA3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [24, 68, 107, 130, 161, 196, 161, 130, 107, 68, 24] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 24 * q ^ 10 + 308 * p ^ 1 * q ^ 9 + 1799 * p ^ 2 * q ^ 8 + 6314 * p ^ 3 * q ^ 7 + 14819 * p ^ 4 * q ^ 6 + 24500 * p ^ 5 * q ^ 5 + 29204 * p ^ 6 * q ^ 4 + 25088 * p ^ 7 * q ^ 3 + 15092 * p ^ 8 * q ^ 2 + 5880 * p ^ 9 * q ^ 1 + 1176 * p ^ 10 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [24, 68, 107, 130, 161, 196, 161, 130, 107, 68, 24] x =
      (fun p q : ℝ => 24 * q ^ 10 + 308 * p ^ 1 * q ^ 9 + 1799 * p ^ 2 * q ^ 8 + 6314 * p ^ 3 * q ^ 7 + 14819 * p ^ 4 * q ^ 6 + 24500 * p ^ 5 * q ^ 5 + 29204 * p ^ 6 * q ^ 4 + 25088 * p ^ 7 * q ^ 3 + 15092 * p ^ 8 * q ^ 2 + 5880 * p ^ 9 * q ^ 1 + 1176 * p ^ 10) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [24, 68, 107, 130, 161, 196, 161, 130, 107, 68, 24] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N7_posA3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [24, 68, 107, 130, 161, 196, 161, 130, 107, 68, 24] x := N7_posA3_0 x hu hv

theorem N7_A3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 7 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 7 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (mul_pos (pow_pos h0 1) (pow_pos h1x 4)) (N7_posA3 x h0 h1)
  rw [← N7_linkA3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N7_linkB3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 7 3 x * (2 * (1 - x ^ 7) ^ 2) = (1 - x) ^ 4 * peval [12, 48, 120, 219, 308, 350, 308, 219, 120, 48, 12] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 7 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (7 - 3 : ℕ) = 4 by norm_num, show (2 * 7 - 3 : ℕ) = 11 by norm_num, show (7 + 3 : ℕ) = 10 by norm_num, show (7 - 3 + 7 : ℕ) = 11 by norm_num, show (3 + 7 : ℕ) = 10 by norm_num, show (7 + 7 : ℕ) = 14 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N7_posB3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 48, 120, 219, 308, 350, 308, 219, 120, 48, 12] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 12 * q ^ 10 + 168 * p ^ 1 * q ^ 9 + 1092 * p ^ 2 * q ^ 8 + 4347 * p ^ 3 * q ^ 7 + 11753 * p ^ 4 * q ^ 6 + 22589 * p ^ 5 * q ^ 5 + 31311 * p ^ 6 * q ^ 4 + 30968 * p ^ 7 * q ^ 3 + 20972 * p ^ 8 * q ^ 2 + 8820 * p ^ 9 * q ^ 1 + 1764 * p ^ 10 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [12, 48, 120, 219, 308, 350, 308, 219, 120, 48, 12] x =
      (fun p q : ℝ => 12 * q ^ 10 + 168 * p ^ 1 * q ^ 9 + 1092 * p ^ 2 * q ^ 8 + 4347 * p ^ 3 * q ^ 7 + 11753 * p ^ 4 * q ^ 6 + 22589 * p ^ 5 * q ^ 5 + 31311 * p ^ 6 * q ^ 4 + 30968 * p ^ 7 * q ^ 3 + 20972 * p ^ 8 * q ^ 2 + 8820 * p ^ 9 * q ^ 1 + 1764 * p ^ 10) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [12, 48, 120, 219, 308, 350, 308, 219, 120, 48, 12] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N7_posB3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 48, 120, 219, 308, 350, 308, 219, 120, 48, 12] x := N7_posB3_0 x hu hv

theorem N7_B3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 7 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 7 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N7_posB3 x h0 h1)
  rw [← N7_linkB3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N7_linkD3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    (aBlock 7 3 x * bBlock 7 3 x - gBlock 7 3 x ^ 2) * (4 * (1 - x) ^ 2 * (1 - x ^ 7) ^ 4) = x ^ 1 * (1 - x) ^ 8 * peval [288, 1968, 7428, 20112, 43296, 78016, 121808, 168768, 210880, 240140, 250640, 240140, 210880, 168768, 121808, 78016, 43296, 20112, 7428, 1968, 288] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 7 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock bBlock gBlock W cR
  simp only [show (7 - 3 : ℕ) = 4 by norm_num, show (2 * 7 - 3 : ℕ) = 11 by norm_num, show (7 + 3 : ℕ) = 10 by norm_num, show (7 - 3 + 7 : ℕ) = 11 by norm_num, show (3 + 7 : ℕ) = 10 by norm_num, show (7 + 7 : ℕ) = 14 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N7_posD3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [288, 1968, 7428, 20112, 43296, 78016, 121808, 168768, 210880, 240140, 250640, 240140, 210880, 168768, 121808, 78016, 43296, 20112, 7428, 1968, 288] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 288 * q ^ 20 + 7728 * p ^ 1 * q ^ 19 + 99540 * p ^ 2 * q ^ 18 + 818664 * p ^ 3 * q ^ 17 + 4824036 * p ^ 4 * q ^ 16 + 21660352 * p ^ 5 * q ^ 15 + 76940192 * p ^ 6 * q ^ 14 + 221542720 * p ^ 7 * q ^ 13 + 525574784 * p ^ 8 * q ^ 12 + 1038294124 * p ^ 9 * q ^ 11 + 1719186364 * p ^ 10 * q ^ 10 + 2392897360 * p ^ 11 * q ^ 9 + 2798887840 * p ^ 12 * q ^ 8 + 2740847144 * p ^ 13 * q ^ 7 + 2229749704 * p ^ 14 * q ^ 6 + 1487451056 * p ^ 15 * q ^ 5 + 796961872 * p ^ 16 * q ^ 4 + 331712556 * p ^ 17 * q ^ 3 + 101331804 * p ^ 18 * q ^ 2 + 20360480 * p ^ 19 * q ^ 1 + 2036048 * p ^ 20 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [288, 1968, 7428, 20112, 43296, 78016, 121808, 168768, 210880, 240140, 250640, 240140, 210880, 168768, 121808, 78016, 43296, 20112, 7428, 1968, 288] x =
      (fun p q : ℝ => 288 * q ^ 20 + 7728 * p ^ 1 * q ^ 19 + 99540 * p ^ 2 * q ^ 18 + 818664 * p ^ 3 * q ^ 17 + 4824036 * p ^ 4 * q ^ 16 + 21660352 * p ^ 5 * q ^ 15 + 76940192 * p ^ 6 * q ^ 14 + 221542720 * p ^ 7 * q ^ 13 + 525574784 * p ^ 8 * q ^ 12 + 1038294124 * p ^ 9 * q ^ 11 + 1719186364 * p ^ 10 * q ^ 10 + 2392897360 * p ^ 11 * q ^ 9 + 2798887840 * p ^ 12 * q ^ 8 + 2740847144 * p ^ 13 * q ^ 7 + 2229749704 * p ^ 14 * q ^ 6 + 1487451056 * p ^ 15 * q ^ 5 + 796961872 * p ^ 16 * q ^ 4 + 331712556 * p ^ 17 * q ^ 3 + 101331804 * p ^ 18 * q ^ 2 + 20360480 * p ^ 19 * q ^ 1 + 2036048 * p ^ 20) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [288, 1968, 7428, 20112, 43296, 78016, 121808, 168768, 210880, 240140, 250640, 240140, 210880, 168768, 121808, 78016, 43296, 20112, 7428, 1968, 288] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N7_posD3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [288, 1968, 7428, 20112, 43296, 78016, 121808, 168768, 210880, 240140, 250640, 240140, 210880, 168768, 121808, 78016, 43296, 20112, 7428, 1968, 288] x := N7_posD3_0 x hu hv

theorem N7_D3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 7 3 x * bBlock 7 3 x - gBlock 7 3 x ^ 2 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 7 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (mul_pos (pow_pos h0 1) (pow_pos h1x 8)) (N7_posD3 x h0 h1)
  rw [← N7_linkD3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos four_pos (pow_pos h1x 2)) (pow_pos hX 4))).mp hn

/-- Note Thm 6(a) for `N = 7`: on all of `(0,1)`, every Fourier block has `a_m > 0`,
`b_m > 0` and `a_m b_m − g_m² > 0` (spectrally stable, positive-definite non-rotational Hessian). -/
theorem stable_N7 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (m : ℕ) (hm1 : 1 ≤ m) (hm2 : m ≤ 6) :
    0 < aBlock 7 m x ∧ 0 < bBlock 7 m x ∧
      0 < aBlock 7 m x * bBlock 7 m x - gBlock 7 m x ^ 2 := by
  interval_cases m
  · exact ⟨N7_A1_pos x h0 h1, N7_B1_pos x h0 h1, N7_D1_pos x h0 h1⟩
  · exact ⟨N7_A2_pos x h0 h1, N7_B2_pos x h0 h1, N7_D2_pos x h0 h1⟩
  · exact ⟨N7_A3_pos x h0 h1, N7_B3_pos x h0 h1, N7_D3_pos x h0 h1⟩
  · rw [show aBlock 7 4 x = aBlock 7 3 x from aBlock_symm 7 3 (by norm_num) x, show bBlock 7 4 x = bBlock 7 3 x from bBlock_symm 7 3 (by norm_num) x, show gBlock 7 4 x = -gBlock 7 3 x from gBlock_symm 7 3 (by norm_num) x, neg_sq]
    exact ⟨N7_A3_pos x h0 h1, N7_B3_pos x h0 h1, N7_D3_pos x h0 h1⟩
  · rw [show aBlock 7 5 x = aBlock 7 2 x from aBlock_symm 7 2 (by norm_num) x, show bBlock 7 5 x = bBlock 7 2 x from bBlock_symm 7 2 (by norm_num) x, show gBlock 7 5 x = -gBlock 7 2 x from gBlock_symm 7 2 (by norm_num) x, neg_sq]
    exact ⟨N7_A2_pos x h0 h1, N7_B2_pos x h0 h1, N7_D2_pos x h0 h1⟩
  · rw [show aBlock 7 6 x = aBlock 7 1 x from aBlock_symm 7 1 (by norm_num) x, show bBlock 7 6 x = bBlock 7 1 x from bBlock_symm 7 1 (by norm_num) x, show gBlock 7 6 x = -gBlock 7 1 x from gBlock_symm 7 1 (by norm_num) x, neg_sq]
    exact ⟨N7_A1_pos x h0 h1, N7_B1_pos x h0 h1, N7_D1_pos x h0 h1⟩

/-! ### N = 8 -/

theorem N8_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 8 1 x * (2 * (1 - x ^ 8) ^ 2) = (1 - x) ^ 4 * peval [7, 20, 38, 60, 85, 112, 140, 112, 85, 60, 38, 20, 7] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 8 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (8 - 1 : ℕ) = 7 by norm_num, show (2 * 8 - 1 : ℕ) = 15 by norm_num, show (8 + 1 : ℕ) = 9 by norm_num, show (8 - 1 + 8 : ℕ) = 15 by norm_num, show (1 + 8 : ℕ) = 9 by norm_num, show (8 + 8 : ℕ) = 16 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N8_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [7, 20, 38, 60, 85, 112, 140, 112, 85, 60, 38, 20, 7] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 7 * q ^ 12 + 104 * p ^ 1 * q ^ 11 + 720 * p ^ 2 * q ^ 10 + 3080 * p ^ 3 * q ^ 9 + 9100 * p ^ 4 * q ^ 8 + 19656 * p ^ 5 * q ^ 7 + 32032 * p ^ 6 * q ^ 6 + 39984 * p ^ 7 * q ^ 5 + 38220 * p ^ 8 * q ^ 4 + 27440 * p ^ 9 * q ^ 3 + 14112 * p ^ 10 * q ^ 2 + 4704 * p ^ 11 * q ^ 1 + 784 * p ^ 12 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [7, 20, 38, 60, 85, 112, 140, 112, 85, 60, 38, 20, 7] x =
      (fun p q : ℝ => 7 * q ^ 12 + 104 * p ^ 1 * q ^ 11 + 720 * p ^ 2 * q ^ 10 + 3080 * p ^ 3 * q ^ 9 + 9100 * p ^ 4 * q ^ 8 + 19656 * p ^ 5 * q ^ 7 + 32032 * p ^ 6 * q ^ 6 + 39984 * p ^ 7 * q ^ 5 + 38220 * p ^ 8 * q ^ 4 + 27440 * p ^ 9 * q ^ 3 + 14112 * p ^ 10 * q ^ 2 + 4704 * p ^ 11 * q ^ 1 + 784 * p ^ 12) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [7, 20, 38, 60, 85, 112, 140, 112, 85, 60, 38, 20, 7] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N8_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [7, 20, 38, 60, 85, 112, 140, 112, 85, 60, 38, 20, 7] x := N8_posB1_0 x hu hv

theorem N8_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 8 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 8 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N8_posB1 x h0 h1)
  rw [← N8_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N8_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 8 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 8) ^ 2) = (1 - x) ^ 4 * peval [7, 34, 85, 156, 243, 342, 449, 504, 449, 342, 243, 156, 85, 34, 7] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 8 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (8 - 1 : ℕ) = 7 by norm_num, show (2 * 8 - 1 : ℕ) = 15 by norm_num, show (8 + 1 : ℕ) = 9 by norm_num, show (8 - 1 + 8 : ℕ) = 15 by norm_num, show (1 + 8 : ℕ) = 9 by norm_num, show (8 + 8 : ℕ) = 16 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N8_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [7, 34, 85, 156, 243, 342, 449, 504, 449, 342, 243, 156, 85, 34, 7] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 7 * q ^ 14 + 132 * p ^ 1 * q ^ 13 + 1164 * p ^ 2 * q ^ 12 + 6376 * p ^ 3 * q ^ 11 + 24300 * p ^ 4 * q ^ 10 + 68376 * p ^ 5 * q ^ 9 + 147056 * p ^ 6 * q ^ 8 + 246736 * p ^ 7 * q ^ 7 + 326284 * p ^ 8 * q ^ 6 + 340256 * p ^ 9 * q ^ 5 + 276752 * p ^ 10 * q ^ 4 + 170912 * p ^ 11 * q ^ 3 + 76048 * p ^ 12 * q ^ 2 + 21952 * p ^ 13 * q ^ 1 + 3136 * p ^ 14 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [7, 34, 85, 156, 243, 342, 449, 504, 449, 342, 243, 156, 85, 34, 7] x =
      (fun p q : ℝ => 7 * q ^ 14 + 132 * p ^ 1 * q ^ 13 + 1164 * p ^ 2 * q ^ 12 + 6376 * p ^ 3 * q ^ 11 + 24300 * p ^ 4 * q ^ 10 + 68376 * p ^ 5 * q ^ 9 + 147056 * p ^ 6 * q ^ 8 + 246736 * p ^ 7 * q ^ 7 + 326284 * p ^ 8 * q ^ 6 + 340256 * p ^ 9 * q ^ 5 + 276752 * p ^ 10 * q ^ 4 + 170912 * p ^ 11 * q ^ 3 + 76048 * p ^ 12 * q ^ 2 + 21952 * p ^ 13 * q ^ 1 + 3136 * p ^ 14) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [7, 34, 85, 156, 243, 342, 449, 504, 449, 342, 243, 156, 85, 34, 7] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N8_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [7, 34, 85, 156, 243, 342, 449, 504, 449, 342, 243, 156, 85, 34, 7] x := N8_posA1_0 x hu hv

theorem N8_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 8 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 8 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N8_posA1 x h0 h1)
  rw [← N8_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N8_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 8 2 x * (2 * (1 - x ^ 8) ^ 2) = (1 - x) ^ 4 * peval [12, 48, 104, 176, 260, 352, 400, 352, 260, 176, 104, 48, 12] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 8 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (8 - 2 : ℕ) = 6 by norm_num, show (2 * 8 - 2 : ℕ) = 14 by norm_num, show (8 + 2 : ℕ) = 10 by norm_num, show (8 - 2 + 8 : ℕ) = 14 by norm_num, show (2 + 8 : ℕ) = 10 by norm_num, show (8 + 8 : ℕ) = 16 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N8_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 48, 104, 176, 260, 352, 400, 352, 260, 176, 104, 48, 12] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 12 * q ^ 12 + 192 * p ^ 1 * q ^ 11 + 1424 * p ^ 2 * q ^ 10 + 6496 * p ^ 3 * q ^ 9 + 20384 * p ^ 4 * q ^ 8 + 46592 * p ^ 5 * q ^ 7 + 80032 * p ^ 6 * q ^ 6 + 104768 * p ^ 7 * q ^ 5 + 104336 * p ^ 8 * q ^ 4 + 77440 * p ^ 9 * q ^ 3 + 40832 * p ^ 10 * q ^ 2 + 13824 * p ^ 11 * q ^ 1 + 2304 * p ^ 12 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [12, 48, 104, 176, 260, 352, 400, 352, 260, 176, 104, 48, 12] x =
      (fun p q : ℝ => 12 * q ^ 12 + 192 * p ^ 1 * q ^ 11 + 1424 * p ^ 2 * q ^ 10 + 6496 * p ^ 3 * q ^ 9 + 20384 * p ^ 4 * q ^ 8 + 46592 * p ^ 5 * q ^ 7 + 80032 * p ^ 6 * q ^ 6 + 104768 * p ^ 7 * q ^ 5 + 104336 * p ^ 8 * q ^ 4 + 77440 * p ^ 9 * q ^ 3 + 40832 * p ^ 10 * q ^ 2 + 13824 * p ^ 11 * q ^ 1 + 2304 * p ^ 12) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [12, 48, 104, 176, 260, 352, 400, 352, 260, 176, 104, 48, 12] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N8_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 48, 104, 176, 260, 352, 400, 352, 260, 176, 104, 48, 12] x := N8_posB2_0 x hu hv

theorem N8_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 8 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 8 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N8_posB2 x h0 h1)
  rw [← N8_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N8_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 8 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 8) ^ 2) = (1 - x) ^ 4 * peval [2, 32, 70, 128, 202, 288, 334, 384, 334, 288, 202, 128, 70, 32, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 8 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (8 - 2 : ℕ) = 6 by norm_num, show (2 * 8 - 2 : ℕ) = 14 by norm_num, show (8 + 2 : ℕ) = 10 by norm_num, show (8 - 2 + 8 : ℕ) = 14 by norm_num, show (2 + 8 : ℕ) = 10 by norm_num, show (8 + 8 : ℕ) = 16 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N8_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 32, 70, 128, 202, 288, 334, 384, 334, 288, 202, 128, 70, 32, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 14 + 60 * p ^ 1 * q ^ 13 + 668 * p ^ 2 * q ^ 12 + 4192 * p ^ 3 * q ^ 11 + 17384 * p ^ 4 * q ^ 10 + 51632 * p ^ 5 * q ^ 9 + 114976 * p ^ 6 * q ^ 8 + 197120 * p ^ 7 * q ^ 7 + 263720 * p ^ 8 * q ^ 6 + 276016 * p ^ 9 * q ^ 5 + 223952 * p ^ 10 * q ^ 4 + 137472 * p ^ 11 * q ^ 3 + 60768 * p ^ 12 * q ^ 2 + 17472 * p ^ 13 * q ^ 1 + 2496 * p ^ 14 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 32, 70, 128, 202, 288, 334, 384, 334, 288, 202, 128, 70, 32, 2] x =
      (fun p q : ℝ => 2 * q ^ 14 + 60 * p ^ 1 * q ^ 13 + 668 * p ^ 2 * q ^ 12 + 4192 * p ^ 3 * q ^ 11 + 17384 * p ^ 4 * q ^ 10 + 51632 * p ^ 5 * q ^ 9 + 114976 * p ^ 6 * q ^ 8 + 197120 * p ^ 7 * q ^ 7 + 263720 * p ^ 8 * q ^ 6 + 276016 * p ^ 9 * q ^ 5 + 223952 * p ^ 10 * q ^ 4 + 137472 * p ^ 11 * q ^ 3 + 60768 * p ^ 12 * q ^ 2 + 17472 * p ^ 13 * q ^ 1 + 2496 * p ^ 14) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 32, 70, 128, 202, 288, 334, 384, 334, 288, 202, 128, 70, 32, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N8_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 32, 70, 128, 202, 288, 334, 384, 334, 288, 202, 128, 70, 32, 2] x := N8_posA2_0 x hu hv

theorem N8_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 8 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 8 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N8_posA2 x h0 h1)
  rw [← N8_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N8_linkB3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 8 3 x * (2 * (1 - x ^ 8) ^ 2) = (1 - x) ^ 4 * peval [15, 60, 150, 276, 429, 560, 620, 560, 429, 276, 150, 60, 15] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 8 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (8 - 3 : ℕ) = 5 by norm_num, show (2 * 8 - 3 : ℕ) = 13 by norm_num, show (8 + 3 : ℕ) = 11 by norm_num, show (8 - 3 + 8 : ℕ) = 13 by norm_num, show (3 + 8 : ℕ) = 11 by norm_num, show (8 + 8 : ℕ) = 16 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N8_posB3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [15, 60, 150, 276, 429, 560, 620, 560, 429, 276, 150, 60, 15] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 15 * q ^ 12 + 240 * p ^ 1 * q ^ 11 + 1800 * p ^ 2 * q ^ 10 + 8376 * p ^ 3 * q ^ 9 + 26988 * p ^ 4 * q ^ 8 + 63608 * p ^ 5 * q ^ 7 + 112816 * p ^ 6 * q ^ 6 + 152240 * p ^ 7 * q ^ 5 + 155660 * p ^ 8 * q ^ 4 + 118000 * p ^ 9 * q ^ 3 + 63200 * p ^ 10 * q ^ 2 + 21600 * p ^ 11 * q ^ 1 + 3600 * p ^ 12 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [15, 60, 150, 276, 429, 560, 620, 560, 429, 276, 150, 60, 15] x =
      (fun p q : ℝ => 15 * q ^ 12 + 240 * p ^ 1 * q ^ 11 + 1800 * p ^ 2 * q ^ 10 + 8376 * p ^ 3 * q ^ 9 + 26988 * p ^ 4 * q ^ 8 + 63608 * p ^ 5 * q ^ 7 + 112816 * p ^ 6 * q ^ 6 + 152240 * p ^ 7 * q ^ 5 + 155660 * p ^ 8 * q ^ 4 + 118000 * p ^ 9 * q ^ 3 + 63200 * p ^ 10 * q ^ 2 + 21600 * p ^ 11 * q ^ 1 + 3600 * p ^ 12) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [15, 60, 150, 276, 429, 560, 620, 560, 429, 276, 150, 60, 15] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N8_posB3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [15, 60, 150, 276, 429, 560, 620, 560, 429, 276, 150, 60, 15] x := N8_posB3_0 x hu hv

theorem N8_B3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 8 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 8 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N8_posB3 x h0 h1)
  rw [← N8_linkB3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N8_linkA3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 8 3 x * (2 * (1 - x) ^ 2 * (1 - x ^ 8) ^ 2) = (1 - x) ^ 4 * peval [-1, 26, 77, 124, 187, 222, 265, 312, 265, 222, 187, 124, 77, 26, -1] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 8 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (8 - 3 : ℕ) = 5 by norm_num, show (2 * 8 - 3 : ℕ) = 13 by norm_num, show (8 + 3 : ℕ) = 11 by norm_num, show (8 - 3 + 8 : ℕ) = 13 by norm_num, show (3 + 8 : ℕ) = 11 by norm_num, show (8 + 8 : ℕ) = 16 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N8_posA3_0 (x : ℝ) (hu : ((1 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-1, 26, 77, 124, 187, 222, 265, 312, 265, 222, 187, 124, 77, 26, -1] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 23037697185962549 * q ^ 14 + 1161983742116078656 * p ^ 1 * q ^ 13 + 15231086235985439744 * p ^ 2 * q ^ 12 + 104104210923103748096 * p ^ 3 * q ^ 11 + 454261449494102278144 * p ^ 4 * q ^ 10 + 1389131640027006107648 * p ^ 5 * q ^ 9 + 3131622602020761894912 * p ^ 6 * q ^ 8 + 5356714817818919436288 * p ^ 7 * q ^ 7 + 7057267384849627348992 * p ^ 8 * q ^ 6 + 7189251487945969369088 * p ^ 9 * q ^ 5 + 5620827272588478644224 * p ^ 10 * q ^ 4 + 3298246652951848288256 * p ^ 11 * q ^ 3 + 1385933245727369068544 * p ^ 12 * q ^ 2 + 377293562382590672896 * p ^ 13 * q ^ 1 + 50728546202701266944 * p ^ 14 := by
    intro p q hp hq; positivity
  have key : (9730975341796875 : ℝ) * peval [-1, 26, 77, 124, 187, 222, 265, 312, 265, 222, 187, 124, 77, 26, -1] x =
      (fun p q : ℝ => 23037697185962549 * q ^ 14 + 1161983742116078656 * p ^ 1 * q ^ 13 + 15231086235985439744 * p ^ 2 * q ^ 12 + 104104210923103748096 * p ^ 3 * q ^ 11 + 454261449494102278144 * p ^ 4 * q ^ 10 + 1389131640027006107648 * p ^ 5 * q ^ 9 + 3131622602020761894912 * p ^ 6 * q ^ 8 + 5356714817818919436288 * p ^ 7 * q ^ 7 + 7057267384849627348992 * p ^ 8 * q ^ 6 + 7189251487945969369088 * p ^ 9 * q ^ 5 + 5620827272588478644224 * p ^ 10 * q ^ 4 + 3298246652951848288256 * p ^ 11 * q ^ 3 + 1385933245727369068544 * p ^ 12 * q ^ 2 + 377293562382590672896 * p ^ 13 * q ^ 1 + 50728546202701266944 * p ^ 14) (x - ((1 : ℝ) / 16)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (9730975341796875 : ℝ) * peval [-1, 26, 77, 124, 187, 222, 265, 312, 265, 222, 187, 124, 77, 26, -1] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N8_posA3 (x : ℝ) (hu : ((1 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-1, 26, 77, 124, 187, 222, 265, 312, 265, 222, 187, 124, 77, 26, -1] x := N8_posA3_0 x hu hv

theorem N8_A3_pos (x : ℝ) (h0 : 0 < x) (hlo : ((1 : ℝ) / 16) < x) (h1 : x < 1) :
    0 < aBlock 8 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 8 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N8_posA3 x hlo h1)
  rw [← N8_linkA3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N8_linkB4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 8 4 x * (2 * (1 - x ^ 8) ^ 2) = (1 - x) ^ 4 * peval [16, 64, 160, 320, 496, 640, 704, 640, 496, 320, 160, 64, 16] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 8 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (8 - 4 : ℕ) = 4 by norm_num, show (2 * 8 - 4 : ℕ) = 12 by norm_num, show (8 + 4 : ℕ) = 12 by norm_num, show (8 - 4 + 8 : ℕ) = 12 by norm_num, show (4 + 8 : ℕ) = 12 by norm_num, show (8 + 8 : ℕ) = 16 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N8_posB4_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [16, 64, 160, 320, 496, 640, 704, 640, 496, 320, 160, 64, 16] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 16 * q ^ 12 + 256 * p ^ 1 * q ^ 11 + 1920 * p ^ 2 * q ^ 10 + 8960 * p ^ 3 * q ^ 9 + 29056 * p ^ 4 * q ^ 8 + 69120 * p ^ 5 * q ^ 7 + 123904 * p ^ 6 * q ^ 6 + 168960 * p ^ 7 * q ^ 5 + 174336 * p ^ 8 * q ^ 4 + 133120 * p ^ 9 * q ^ 3 + 71680 * p ^ 10 * q ^ 2 + 24576 * p ^ 11 * q ^ 1 + 4096 * p ^ 12 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [16, 64, 160, 320, 496, 640, 704, 640, 496, 320, 160, 64, 16] x =
      (fun p q : ℝ => 16 * q ^ 12 + 256 * p ^ 1 * q ^ 11 + 1920 * p ^ 2 * q ^ 10 + 8960 * p ^ 3 * q ^ 9 + 29056 * p ^ 4 * q ^ 8 + 69120 * p ^ 5 * q ^ 7 + 123904 * p ^ 6 * q ^ 6 + 168960 * p ^ 7 * q ^ 5 + 174336 * p ^ 8 * q ^ 4 + 133120 * p ^ 9 * q ^ 3 + 71680 * p ^ 10 * q ^ 2 + 24576 * p ^ 11 * q ^ 1 + 4096 * p ^ 12) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [16, 64, 160, 320, 496, 640, 704, 640, 496, 320, 160, 64, 16] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N8_posB4 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [16, 64, 160, 320, 496, 640, 704, 640, 496, 320, 160, 64, 16] x := N8_posB4_0 x hu hv

theorem N8_B4_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 8 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 8 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N8_posB4 x h0 h1)
  rw [← N8_linkB4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N8_linkA4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 8 4 x * (2 * (1 - x) ^ 2 * (1 - x ^ 8) ^ 2) = (1 - x) ^ 4 * peval [-2, 24, 74, 144, 166, 200, 242, 288, 242, 200, 166, 144, 74, 24, -2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 8 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (8 - 4 : ℕ) = 4 by norm_num, show (2 * 8 - 4 : ℕ) = 12 by norm_num, show (8 + 4 : ℕ) = 12 by norm_num, show (8 - 4 + 8 : ℕ) = 12 by norm_num, show (4 + 8 : ℕ) = 12 by norm_num, show (8 + 8 : ℕ) = 16 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

/-- The minimal polynomial `f_8` of `ξ*(8)` (constant term first). -/
def f8 : List ℤ := [1, -12, -37, -72, -84, -88, -84, -72, -37, -12, 1]

theorem N8_crit_factor (x : ℝ) : peval [-2, 24, 74, 144, 166, 200, 242, 288, 242, 200, 166, 144, 74, 24, -2] x = -(peval [2, 0, 0, 0, 2] x * peval f8 x) := by
  simp only [f8, peval]; ring

theorem N8_posH_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 0, 0, 0, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 4 + 8 * p ^ 1 * q ^ 3 + 12 * p ^ 2 * q ^ 2 + 8 * p ^ 3 * q ^ 1 + 4 * p ^ 4 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 0, 0, 0, 2] x =
      (fun p q : ℝ => 2 * q ^ 4 + 8 * p ^ 1 * q ^ 3 + 12 * p ^ 2 * q ^ 2 + 8 * p ^ 3 * q ^ 1 + 4 * p ^ 4) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 0, 0, 0, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N8_posH (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 0, 0, 0, 2] x := N8_posH_0 x hu hv

theorem N8_dneg_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 74, 216, 336, 440, 504, 504, 296, 108, -10] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 12 * q ^ 9 + 182 * p ^ 1 * q ^ 8 + 1240 * p ^ 2 * q ^ 7 + 4928 * p ^ 3 * q ^ 6 + 12648 * p ^ 4 * q ^ 5 + 21996 * p ^ 5 * q ^ 4 + 26352 * p ^ 6 * q ^ 3 + 21312 * p ^ 7 * q ^ 2 + 10656 * p ^ 8 * q ^ 1 + 2480 * p ^ 9 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [12, 74, 216, 336, 440, 504, 504, 296, 108, -10] x =
      (fun p q : ℝ => 12 * q ^ 9 + 182 * p ^ 1 * q ^ 8 + 1240 * p ^ 2 * q ^ 7 + 4928 * p ^ 3 * q ^ 6 + 12648 * p ^ 4 * q ^ 5 + 21996 * p ^ 5 * q ^ 4 + 26352 * p ^ 6 * q ^ 3 + 21312 * p ^ 7 * q ^ 2 + 10656 * p ^ 8 * q ^ 1 + 2480 * p ^ 9) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [12, 74, 216, 336, 440, 504, 504, 296, 108, -10] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N8_dneg (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 74, 216, 336, 440, 504, 504, 296, 108, -10] x := N8_dneg_0 x hu hv

theorem N8_dlink (x : ℝ) : dpeval f8 x = -peval [12, 74, 216, 336, 440, 504, 504, 296, 108, -10] x := by
  simp only [f8, dpeval, peval]; ring

/-- `f_8` changes sign exactly once in `(0,1)`, at a point of `(3367811150204063113/50000000000000000000, 6735622300408126227/100000000000000000000)`. -/
theorem N8_root : ∃ ξs, ((3367811150204063113 : ℝ) / 50000000000000000000) < ξs ∧ ξs < ((6735622300408126227 : ℝ) / 100000000000000000000) ∧ peval f8 ξs = 0 ∧
    ∀ x, 0 < x → x < 1 → (x < ξs → 0 < peval f8 x) ∧ (ξs < x → peval f8 x < 0) := by
  apply sign_change (peval f8) (0 : ℝ) _ _ (by norm_num) (by norm_num) (by norm_num)
    (continuous_peval _).continuousOn
  · apply peval_strictAntiOn
    intro x hu hv
    rw [N8_dlink]
    linarith [N8_dneg x hu hv]
  · exact (fun x hx0 hx => absurd hx (not_le.mpr hx0))
  · norm_num [f8, peval]
  · norm_num [f8, peval]

theorem N8_A4_neg (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : 0 < peval f8 x) :
    aBlock 8 4 x < 0 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 8 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : (1 - x) ^ 4 * peval [-2, 24, 74, 144, 166, 200, 242, 288, 242, 200, 166, 144, 74, 24, -2] x < 0 := by
    rw [N8_crit_factor]
    have := mul_pos (N8_posH x h0 h1) hf
    nlinarith [pow_pos h1x 4]
  rw [← N8_linkA4 x h0 h1] at hn
  have hd := (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))
  by_contra hc
  push Not at hc
  nlinarith [mul_nonneg hc hd.le]

theorem N8_A4_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : peval f8 x < 0) :
    0 < aBlock 8 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 8 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : 0 < (1 - x) ^ 4 * peval [-2, 24, 74, 144, 166, 200, 242, 288, 242, 200, 166, 144, 74, 24, -2] x := by
    rw [N8_crit_factor]
    have := mul_pos (N8_posH x h0 h1) (neg_pos.mpr hf)
    have h4 := pow_pos h1x 4
    nlinarith
  rw [← N8_linkA4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

/-- Note Thm 6(b) for `N = 8`: there is a unique `ξ* ∈ (0,1)` (the unique root of `f_8` in
`(0,1)`, enclosed in `(3367811150204063113/50000000000000000000, 6735622300408126227/100000000000000000000)`) such that, for every `ξ ∈ (0,1)`: all `b_m > 0`;
`a_4 < 0` for `ξ < ξ*` (spectrally unstable, since `b_4 > 0`); and all `a_m > 0` for `ξ > ξ*`. -/
theorem cert_N8 : ∃ ξs : ℝ, ((3367811150204063113 : ℝ) / 50000000000000000000) < ξs ∧ ξs < ((6735622300408126227 : ℝ) / 100000000000000000000) ∧ peval f8 ξs = 0 ∧
    (∀ y, 0 < y → y < 1 → peval f8 y = 0 → y = ξs) ∧
    ∀ x, 0 < x → x < 1 →
      (∀ m, 1 ≤ m → m ≤ 7 → 0 < bBlock 8 m x) ∧
      (x < ξs → aBlock 8 4 x < 0) ∧
      (ξs < x → ∀ m, 1 ≤ m → m ≤ 7 → 0 < aBlock 8 m x) := by
  obtain ⟨ξs, hta, htb, hf0, hs⟩ := N8_root
  refine ⟨ξs, hta, htb, hf0, ?_, ?_⟩
  · intro y hy0 hy1 hy
    rcases lt_trichotomy y ξs with h | h | h
    · linarith [(hs y hy0 hy1).1 h]
    · exact h
    · linarith [(hs y hy0 hy1).2 h]
  intro x h0 h1
  refine ⟨?_, fun hlt => N8_A4_neg x h0 h1 ((hs x h0 h1).1 hlt), ?_⟩
  · intro m hm1 hm2
    interval_cases m
    · exact N8_B1_pos x h0 h1
    · exact N8_B2_pos x h0 h1
    · exact N8_B3_pos x h0 h1
    · exact N8_B4_pos x h0 h1
    · rw [show bBlock 8 5 x = bBlock 8 3 x from bBlock_symm 8 3 (by norm_num) x]
      exact N8_B3_pos x h0 h1
    · rw [show bBlock 8 6 x = bBlock 8 2 x from bBlock_symm 8 2 (by norm_num) x]
      exact N8_B2_pos x h0 h1
    · rw [show bBlock 8 7 x = bBlock 8 1 x from bBlock_symm 8 1 (by norm_num) x]
      exact N8_B1_pos x h0 h1
  · intro hgt m hm1 hm2
    have hlo : ((1 : ℝ) / 16) < x := lt_trans (lt_trans (by norm_num) hta) hgt
    have hfn := (hs x h0 h1).2 hgt
    interval_cases m
    · exact N8_A1_pos x h0 h1
    · exact N8_A2_pos x h0 h1
    · exact N8_A3_pos x h0 hlo h1
    · exact N8_A4_pos x h0 h1 hfn
    · rw [show aBlock 8 5 x = aBlock 8 3 x from aBlock_symm 8 3 (by norm_num) x]
      exact N8_A3_pos x h0 hlo h1
    · rw [show aBlock 8 6 x = aBlock 8 2 x from aBlock_symm 8 2 (by norm_num) x]
      exact N8_A2_pos x h0 h1
    · rw [show aBlock 8 7 x = aBlock 8 1 x from aBlock_symm 8 1 (by norm_num) x]
      exact N8_A1_pos x h0 h1

/-! ### N = 9 -/

theorem N9_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 9 1 x * (2 * (1 - x ^ 9) ^ 2) = (1 - x) ^ 4 * peval [8, 23, 44, 70, 100, 133, 168, 204, 168, 133, 100, 70, 44, 23, 8] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 9 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (9 - 1 : ℕ) = 8 by norm_num, show (2 * 9 - 1 : ℕ) = 17 by norm_num, show (9 + 1 : ℕ) = 10 by norm_num, show (9 - 1 + 9 : ℕ) = 17 by norm_num, show (1 + 9 : ℕ) = 10 by norm_num, show (9 + 9 : ℕ) = 18 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N9_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [8, 23, 44, 70, 100, 133, 168, 204, 168, 133, 100, 70, 44, 23, 8] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 8 * q ^ 14 + 135 * p ^ 1 * q ^ 13 + 1071 * p ^ 2 * q ^ 12 + 5304 * p ^ 3 * q ^ 11 + 18360 * p ^ 4 * q ^ 10 + 47124 * p ^ 5 * q ^ 9 + 92820 * p ^ 6 * q ^ 8 + 143208 * p ^ 7 * q ^ 7 + 174960 * p ^ 8 * q ^ 6 + 169596 * p ^ 9 * q ^ 5 + 129276 * p ^ 10 * q ^ 4 + 75600 * p ^ 11 * q ^ 3 + 32256 * p ^ 12 * q ^ 2 + 9072 * p ^ 13 * q ^ 1 + 1296 * p ^ 14 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [8, 23, 44, 70, 100, 133, 168, 204, 168, 133, 100, 70, 44, 23, 8] x =
      (fun p q : ℝ => 8 * q ^ 14 + 135 * p ^ 1 * q ^ 13 + 1071 * p ^ 2 * q ^ 12 + 5304 * p ^ 3 * q ^ 11 + 18360 * p ^ 4 * q ^ 10 + 47124 * p ^ 5 * q ^ 9 + 92820 * p ^ 6 * q ^ 8 + 143208 * p ^ 7 * q ^ 7 + 174960 * p ^ 8 * q ^ 6 + 169596 * p ^ 9 * q ^ 5 + 129276 * p ^ 10 * q ^ 4 + 75600 * p ^ 11 * q ^ 3 + 32256 * p ^ 12 * q ^ 2 + 9072 * p ^ 13 * q ^ 1 + 1296 * p ^ 14) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [8, 23, 44, 70, 100, 133, 168, 204, 168, 133, 100, 70, 44, 23, 8] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N9_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [8, 23, 44, 70, 100, 133, 168, 204, 168, 133, 100, 70, 44, 23, 8] x := N9_posB1_0 x hu hv

theorem N9_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 9 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 9 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N9_posB1 x h0 h1)
  rw [← N9_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N9_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 9 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 9) ^ 2) = (1 - x) ^ 4 * peval [8, 39, 98, 181, 284, 403, 534, 673, 744, 673, 534, 403, 284, 181, 98, 39, 8] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 9 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (9 - 1 : ℕ) = 8 by norm_num, show (2 * 9 - 1 : ℕ) = 17 by norm_num, show (9 + 1 : ℕ) = 10 by norm_num, show (9 - 1 + 9 : ℕ) = 17 by norm_num, show (1 + 9 : ℕ) = 10 by norm_num, show (9 + 9 : ℕ) = 18 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N9_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [8, 39, 98, 181, 284, 403, 534, 673, 744, 673, 534, 403, 284, 181, 98, 39, 8] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 8 * q ^ 16 + 167 * p ^ 1 * q ^ 15 + 1643 * p ^ 2 * q ^ 14 + 10128 * p ^ 3 * q ^ 13 + 43860 * p ^ 4 * q ^ 12 + 141780 * p ^ 5 * q ^ 11 + 354756 * p ^ 6 * q ^ 10 + 702984 * p ^ 7 * q ^ 9 + 1119072 * p ^ 8 * q ^ 8 + 1442268 * p ^ 9 * q ^ 7 + 1507500 * p ^ 10 * q ^ 6 + 1271088 * p ^ 11 * q ^ 5 + 851760 * p ^ 12 * q ^ 4 + 440496 * p ^ 13 * q ^ 3 + 166608 * p ^ 14 * q ^ 2 + 41472 * p ^ 15 * q ^ 1 + 5184 * p ^ 16 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [8, 39, 98, 181, 284, 403, 534, 673, 744, 673, 534, 403, 284, 181, 98, 39, 8] x =
      (fun p q : ℝ => 8 * q ^ 16 + 167 * p ^ 1 * q ^ 15 + 1643 * p ^ 2 * q ^ 14 + 10128 * p ^ 3 * q ^ 13 + 43860 * p ^ 4 * q ^ 12 + 141780 * p ^ 5 * q ^ 11 + 354756 * p ^ 6 * q ^ 10 + 702984 * p ^ 7 * q ^ 9 + 1119072 * p ^ 8 * q ^ 8 + 1442268 * p ^ 9 * q ^ 7 + 1507500 * p ^ 10 * q ^ 6 + 1271088 * p ^ 11 * q ^ 5 + 851760 * p ^ 12 * q ^ 4 + 440496 * p ^ 13 * q ^ 3 + 166608 * p ^ 14 * q ^ 2 + 41472 * p ^ 15 * q ^ 1 + 5184 * p ^ 16) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [8, 39, 98, 181, 284, 403, 534, 673, 744, 673, 534, 403, 284, 181, 98, 39, 8] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N9_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [8, 39, 98, 181, 284, 403, 534, 673, 744, 673, 534, 403, 284, 181, 98, 39, 8] x := N9_posA1_0 x hu hv

theorem N9_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 9 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 9 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N9_posA1 x h0 h1)
  rw [← N9_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N9_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 9 2 x * (2 * (1 - x ^ 9) ^ 2) = (1 - x) ^ 4 * peval [14, 56, 122, 208, 310, 424, 546, 609, 546, 424, 310, 208, 122, 56, 14] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 9 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (9 - 2 : ℕ) = 7 by norm_num, show (2 * 9 - 2 : ℕ) = 16 by norm_num, show (9 + 2 : ℕ) = 11 by norm_num, show (9 - 2 + 9 : ℕ) = 16 by norm_num, show (2 + 9 : ℕ) = 11 by norm_num, show (9 + 9 : ℕ) = 18 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N9_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [14, 56, 122, 208, 310, 424, 546, 609, 546, 424, 310, 208, 122, 56, 14] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 14 * q ^ 14 + 252 * p ^ 1 * q ^ 13 + 2124 * p ^ 2 * q ^ 12 + 11136 * p ^ 3 * q ^ 11 + 40680 * p ^ 4 * q ^ 10 + 109872 * p ^ 5 * q ^ 9 + 227136 * p ^ 6 * q ^ 8 + 366849 * p ^ 7 * q ^ 7 + 467775 * p ^ 8 * q ^ 6 + 471429 * p ^ 9 * q ^ 5 + 371763 * p ^ 10 * q ^ 4 + 223587 * p ^ 11 * q ^ 3 + 97461 * p ^ 12 * q ^ 2 + 27783 * p ^ 13 * q ^ 1 + 3969 * p ^ 14 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [14, 56, 122, 208, 310, 424, 546, 609, 546, 424, 310, 208, 122, 56, 14] x =
      (fun p q : ℝ => 14 * q ^ 14 + 252 * p ^ 1 * q ^ 13 + 2124 * p ^ 2 * q ^ 12 + 11136 * p ^ 3 * q ^ 11 + 40680 * p ^ 4 * q ^ 10 + 109872 * p ^ 5 * q ^ 9 + 227136 * p ^ 6 * q ^ 8 + 366849 * p ^ 7 * q ^ 7 + 467775 * p ^ 8 * q ^ 6 + 471429 * p ^ 9 * q ^ 5 + 371763 * p ^ 10 * q ^ 4 + 223587 * p ^ 11 * q ^ 3 + 97461 * p ^ 12 * q ^ 2 + 27783 * p ^ 13 * q ^ 1 + 3969 * p ^ 14) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [14, 56, 122, 208, 310, 424, 546, 609, 546, 424, 310, 208, 122, 56, 14] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N9_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [14, 56, 122, 208, 310, 424, 546, 609, 546, 424, 310, 208, 122, 56, 14] x := N9_posB2_0 x hu hv

theorem N9_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 9 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 9 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N9_posB2 x h0 h1)
  rw [← N9_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N9_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 9 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 9) ^ 2) = (1 - x) ^ 4 * peval [2, 36, 80, 148, 236, 340, 456, 517, 582, 517, 456, 340, 236, 148, 80, 36, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 9 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (9 - 2 : ℕ) = 7 by norm_num, show (2 * 9 - 2 : ℕ) = 16 by norm_num, show (9 + 2 : ℕ) = 11 by norm_num, show (9 - 2 + 9 : ℕ) = 16 by norm_num, show (2 + 9 : ℕ) = 11 by norm_num, show (9 + 9 : ℕ) = 18 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N9_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 36, 80, 148, 236, 340, 456, 517, 582, 517, 456, 340, 236, 148, 80, 36, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 16 + 68 * p ^ 1 * q ^ 15 + 860 * p ^ 2 * q ^ 14 + 6168 * p ^ 3 * q ^ 13 + 29460 * p ^ 4 * q ^ 12 + 101712 * p ^ 5 * q ^ 11 + 266304 * p ^ 6 * q ^ 10 + 544737 * p ^ 7 * q ^ 9 + 886791 * p ^ 8 * q ^ 8 + 1160685 * p ^ 9 * q ^ 7 + 1225107 * p ^ 10 * q ^ 6 + 1038051 * p ^ 11 * q ^ 5 + 696213 * p ^ 12 * q ^ 4 + 359415 * p ^ 13 * q ^ 3 + 135585 * p ^ 14 * q ^ 2 + 33696 * p ^ 15 * q ^ 1 + 4212 * p ^ 16 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 36, 80, 148, 236, 340, 456, 517, 582, 517, 456, 340, 236, 148, 80, 36, 2] x =
      (fun p q : ℝ => 2 * q ^ 16 + 68 * p ^ 1 * q ^ 15 + 860 * p ^ 2 * q ^ 14 + 6168 * p ^ 3 * q ^ 13 + 29460 * p ^ 4 * q ^ 12 + 101712 * p ^ 5 * q ^ 11 + 266304 * p ^ 6 * q ^ 10 + 544737 * p ^ 7 * q ^ 9 + 886791 * p ^ 8 * q ^ 8 + 1160685 * p ^ 9 * q ^ 7 + 1225107 * p ^ 10 * q ^ 6 + 1038051 * p ^ 11 * q ^ 5 + 696213 * p ^ 12 * q ^ 4 + 359415 * p ^ 13 * q ^ 3 + 135585 * p ^ 14 * q ^ 2 + 33696 * p ^ 15 * q ^ 1 + 4212 * p ^ 16) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 36, 80, 148, 236, 340, 456, 517, 582, 517, 456, 340, 236, 148, 80, 36, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N9_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 36, 80, 148, 236, 340, 456, 517, 582, 517, 456, 340, 236, 148, 80, 36, 2] x := N9_posA2_0 x hu hv

theorem N9_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 9 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 9 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N9_posA2 x h0 h1)
  rw [← N9_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N9_linkB3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 9 3 x * (2 * (1 - x ^ 9) ^ 2) = (1 - x) ^ 4 * peval [18, 72, 180, 333, 522, 738, 918, 999, 918, 738, 522, 333, 180, 72, 18] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 9 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (9 - 3 : ℕ) = 6 by norm_num, show (2 * 9 - 3 : ℕ) = 15 by norm_num, show (9 + 3 : ℕ) = 12 by norm_num, show (9 - 3 + 9 : ℕ) = 15 by norm_num, show (3 + 9 : ℕ) = 12 by norm_num, show (9 + 9 : ℕ) = 18 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N9_posB3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [18, 72, 180, 333, 522, 738, 918, 999, 918, 738, 522, 333, 180, 72, 18] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 18 * q ^ 14 + 324 * p ^ 1 * q ^ 13 + 2754 * p ^ 2 * q ^ 12 + 14661 * p ^ 3 * q ^ 11 + 54675 * p ^ 4 * q ^ 10 + 151389 * p ^ 5 * q ^ 9 + 321813 * p ^ 6 * q ^ 8 + 535329 * p ^ 7 * q ^ 7 + 702999 * p ^ 8 * q ^ 6 + 728271 * p ^ 9 * q ^ 5 + 588303 * p ^ 10 * q ^ 4 + 360855 * p ^ 11 * q ^ 3 + 159651 * p ^ 12 * q ^ 2 + 45927 * p ^ 13 * q ^ 1 + 6561 * p ^ 14 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [18, 72, 180, 333, 522, 738, 918, 999, 918, 738, 522, 333, 180, 72, 18] x =
      (fun p q : ℝ => 18 * q ^ 14 + 324 * p ^ 1 * q ^ 13 + 2754 * p ^ 2 * q ^ 12 + 14661 * p ^ 3 * q ^ 11 + 54675 * p ^ 4 * q ^ 10 + 151389 * p ^ 5 * q ^ 9 + 321813 * p ^ 6 * q ^ 8 + 535329 * p ^ 7 * q ^ 7 + 702999 * p ^ 8 * q ^ 6 + 728271 * p ^ 9 * q ^ 5 + 588303 * p ^ 10 * q ^ 4 + 360855 * p ^ 11 * q ^ 3 + 159651 * p ^ 12 * q ^ 2 + 45927 * p ^ 13 * q ^ 1 + 6561 * p ^ 14) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [18, 72, 180, 333, 522, 738, 918, 999, 918, 738, 522, 333, 180, 72, 18] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N9_posB3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [18, 72, 180, 333, 522, 738, 918, 999, 918, 738, 522, 333, 180, 72, 18] x := N9_posB3_0 x hu hv

theorem N9_B3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 9 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 9 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N9_posB3 x h0 h1)
  rw [← N9_linkB3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N9_linkA3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 9 3 x * (2 * (1 - x) ^ 2 * (1 - x ^ 9) ^ 2) = (1 - x) ^ 4 * peval [-2, 28, 86, 141, 216, 307, 356, 413, 474, 413, 356, 307, 216, 141, 86, 28, -2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 9 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (9 - 3 : ℕ) = 6 by norm_num, show (2 * 9 - 3 : ℕ) = 15 by norm_num, show (9 + 3 : ℕ) = 12 by norm_num, show (9 - 3 + 9 : ℕ) = 15 by norm_num, show (3 + 9 : ℕ) = 12 by norm_num, show (9 + 9 : ℕ) = 18 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N9_posA3_0 (x : ℝ) (hu : ((1 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-2, 28, 86, 141, 216, 307, 356, 413, 474, 413, 356, 307, 216, 141, 86, 28, -2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 28233283612961374 * q ^ 16 + 9128305482979529664 * p ^ 1 * q ^ 15 + 157223543733606586880 * p ^ 2 * q ^ 14 + 1297683004240078295040 * p ^ 3 * q ^ 13 + 6723796729239139450880 * p ^ 4 * q ^ 12 + 24438295484029476012032 * p ^ 5 * q ^ 11 + 66085954231785154936832 * p ^ 6 * q ^ 10 + 137646978865539118530560 * p ^ 7 * q ^ 9 + 225453668057076549550080 * p ^ 8 * q ^ 8 + 293702685774555601960960 * p ^ 9 * q ^ 7 + 305434800032904005025792 * p ^ 10 * q ^ 6 + 252566211792348246966272 * p ^ 11 * q ^ 5 + 163882937843889922375680 * p ^ 12 * q ^ 4 + 81238699792279839703040 * p ^ 13 * q ^ 3 + 29248537707965138862080 * p ^ 14 * q ^ 2 + 6899082283567372304384 * p ^ 15 * q ^ 1 + 811656739243220271104 * p ^ 16 := by
    intro p q hp hq; positivity
  have key : (81091461181640625 : ℝ) * peval [-2, 28, 86, 141, 216, 307, 356, 413, 474, 413, 356, 307, 216, 141, 86, 28, -2] x =
      (fun p q : ℝ => 28233283612961374 * q ^ 16 + 9128305482979529664 * p ^ 1 * q ^ 15 + 157223543733606586880 * p ^ 2 * q ^ 14 + 1297683004240078295040 * p ^ 3 * q ^ 13 + 6723796729239139450880 * p ^ 4 * q ^ 12 + 24438295484029476012032 * p ^ 5 * q ^ 11 + 66085954231785154936832 * p ^ 6 * q ^ 10 + 137646978865539118530560 * p ^ 7 * q ^ 9 + 225453668057076549550080 * p ^ 8 * q ^ 8 + 293702685774555601960960 * p ^ 9 * q ^ 7 + 305434800032904005025792 * p ^ 10 * q ^ 6 + 252566211792348246966272 * p ^ 11 * q ^ 5 + 163882937843889922375680 * p ^ 12 * q ^ 4 + 81238699792279839703040 * p ^ 13 * q ^ 3 + 29248537707965138862080 * p ^ 14 * q ^ 2 + 6899082283567372304384 * p ^ 15 * q ^ 1 + 811656739243220271104 * p ^ 16) (x - ((1 : ℝ) / 16)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (81091461181640625 : ℝ) * peval [-2, 28, 86, 141, 216, 307, 356, 413, 474, 413, 356, 307, 216, 141, 86, 28, -2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N9_posA3 (x : ℝ) (hu : ((1 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-2, 28, 86, 141, 216, 307, 356, 413, 474, 413, 356, 307, 216, 141, 86, 28, -2] x := N9_posA3_0 x hu hv

theorem N9_A3_pos (x : ℝ) (h0 : 0 < x) (hlo : ((1 : ℝ) / 16) < x) (h1 : x < 1) :
    0 < aBlock 9 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 9 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N9_posA3 x hlo h1)
  rw [← N9_linkA3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N9_linkB4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 9 4 x * (2 * (1 - x ^ 9) ^ 2) = (1 - x) ^ 4 * peval [20, 80, 200, 400, 664, 931, 1140, 1230, 1140, 931, 664, 400, 200, 80, 20] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 9 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (9 - 4 : ℕ) = 5 by norm_num, show (2 * 9 - 4 : ℕ) = 14 by norm_num, show (9 + 4 : ℕ) = 13 by norm_num, show (9 - 4 + 9 : ℕ) = 14 by norm_num, show (4 + 9 : ℕ) = 13 by norm_num, show (9 + 9 : ℕ) = 18 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N9_posB4_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [20, 80, 200, 400, 664, 931, 1140, 1230, 1140, 931, 664, 400, 200, 80, 20] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 20 * q ^ 14 + 360 * p ^ 1 * q ^ 13 + 3060 * p ^ 2 * q ^ 12 + 16320 * p ^ 3 * q ^ 11 + 61164 * p ^ 4 * q ^ 10 + 170811 * p ^ 5 * q ^ 9 + 367419 * p ^ 6 * q ^ 8 + 619866 * p ^ 7 * q ^ 7 + 826254 * p ^ 8 * q ^ 6 + 868275 * p ^ 9 * q ^ 5 + 710235 * p ^ 10 * q ^ 4 + 440100 * p ^ 11 * q ^ 3 + 196200 * p ^ 12 * q ^ 2 + 56700 * p ^ 13 * q ^ 1 + 8100 * p ^ 14 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [20, 80, 200, 400, 664, 931, 1140, 1230, 1140, 931, 664, 400, 200, 80, 20] x =
      (fun p q : ℝ => 20 * q ^ 14 + 360 * p ^ 1 * q ^ 13 + 3060 * p ^ 2 * q ^ 12 + 16320 * p ^ 3 * q ^ 11 + 61164 * p ^ 4 * q ^ 10 + 170811 * p ^ 5 * q ^ 9 + 367419 * p ^ 6 * q ^ 8 + 619866 * p ^ 7 * q ^ 7 + 826254 * p ^ 8 * q ^ 6 + 868275 * p ^ 9 * q ^ 5 + 710235 * p ^ 10 * q ^ 4 + 440100 * p ^ 11 * q ^ 3 + 196200 * p ^ 12 * q ^ 2 + 56700 * p ^ 13 * q ^ 1 + 8100 * p ^ 14) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [20, 80, 200, 400, 664, 931, 1140, 1230, 1140, 931, 664, 400, 200, 80, 20] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N9_posB4 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [20, 80, 200, 400, 664, 931, 1140, 1230, 1140, 931, 664, 400, 200, 80, 20] x := N9_posB4_0 x hu hv

theorem N9_B4_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 9 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 9 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N9_posB4 x h0 h1)
  rw [← N9_linkB4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N9_linkA4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 9 4 x * (2 * (1 - x) ^ 2 * (1 - x ^ 9) ^ 2) = (1 - x) ^ 4 * peval [-4, 24, 80, 160, 224, 259, 306, 361, 420, 361, 306, 259, 224, 160, 80, 24, -4] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 9 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (9 - 4 : ℕ) = 5 by norm_num, show (2 * 9 - 4 : ℕ) = 14 by norm_num, show (9 + 4 : ℕ) = 13 by norm_num, show (9 - 4 + 9 : ℕ) = 14 by norm_num, show (4 + 9 : ℕ) = 13 by norm_num, show (9 + 9 : ℕ) = 18 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

/-- The minimal polynomial `f_9` of `ξ*(9)` (constant term first). -/
def f9 : List ℤ := [4, -24, -80, -160, -224, -259, -306, -361, -420, -361, -306, -259, -224, -160, -80, -24, 4]

theorem N9_crit_factor (x : ℝ) : peval [-4, 24, 80, 160, 224, 259, 306, 361, 420, 361, 306, 259, 224, 160, 80, 24, -4] x = -(peval [1] x * peval f9 x) := by
  simp only [f9, peval]; ring

theorem N9_posH (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1] x := by norm_num [peval]

theorem N9_dneg_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [24, 160, 480, 896, 1295, 1836, 2527, 3360, 3249, 3060, 2849, 2688, 2080, 1120, 360, -64] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 24 * q ^ 15 + 520 * p ^ 1 * q ^ 14 + 5240 * p ^ 2 * q ^ 13 + 32616 * p ^ 3 * q ^ 12 + 140487 * p ^ 4 * q ^ 11 + 444729 * p ^ 5 * q ^ 10 + 1072872 * p ^ 6 * q ^ 9 + 2018598 * p ^ 7 * q ^ 8 + 3005643 * p ^ 8 * q ^ 7 + 3568185 * p ^ 9 * q ^ 6 + 3376746 * p ^ 10 * q ^ 5 + 2522700 * p ^ 11 * q ^ 4 + 1454220 * p ^ 12 * q ^ 3 + 618660 * p ^ 13 * q ^ 2 + 176760 * p ^ 14 * q ^ 1 + 25920 * p ^ 15 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [24, 160, 480, 896, 1295, 1836, 2527, 3360, 3249, 3060, 2849, 2688, 2080, 1120, 360, -64] x =
      (fun p q : ℝ => 24 * q ^ 15 + 520 * p ^ 1 * q ^ 14 + 5240 * p ^ 2 * q ^ 13 + 32616 * p ^ 3 * q ^ 12 + 140487 * p ^ 4 * q ^ 11 + 444729 * p ^ 5 * q ^ 10 + 1072872 * p ^ 6 * q ^ 9 + 2018598 * p ^ 7 * q ^ 8 + 3005643 * p ^ 8 * q ^ 7 + 3568185 * p ^ 9 * q ^ 6 + 3376746 * p ^ 10 * q ^ 5 + 2522700 * p ^ 11 * q ^ 4 + 1454220 * p ^ 12 * q ^ 3 + 618660 * p ^ 13 * q ^ 2 + 176760 * p ^ 14 * q ^ 1 + 25920 * p ^ 15) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [24, 160, 480, 896, 1295, 1836, 2527, 3360, 3249, 3060, 2849, 2688, 2080, 1120, 360, -64] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N9_dneg (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [24, 160, 480, 896, 1295, 1836, 2527, 3360, 3249, 3060, 2849, 2688, 2080, 1120, 360, -64] x := N9_dneg_0 x hu hv

theorem N9_dlink (x : ℝ) : dpeval f9 x = -peval [24, 160, 480, 896, 1295, 1836, 2527, 3360, 3249, 3060, 2849, 2688, 2080, 1120, 360, -64] x := by
  simp only [f9, dpeval, peval]; ring

/-- `f_9` changes sign exactly once in `(0,1)`, at a point of `(11287106323660274043/100000000000000000000, 2821776580915068511/25000000000000000000)`. -/
theorem N9_root : ∃ ξs, ((11287106323660274043 : ℝ) / 100000000000000000000) < ξs ∧ ξs < ((2821776580915068511 : ℝ) / 25000000000000000000) ∧ peval f9 ξs = 0 ∧
    ∀ x, 0 < x → x < 1 → (x < ξs → 0 < peval f9 x) ∧ (ξs < x → peval f9 x < 0) := by
  apply sign_change (peval f9) (0 : ℝ) _ _ (by norm_num) (by norm_num) (by norm_num)
    (continuous_peval _).continuousOn
  · apply peval_strictAntiOn
    intro x hu hv
    rw [N9_dlink]
    linarith [N9_dneg x hu hv]
  · exact (fun x hx0 hx => absurd hx (not_le.mpr hx0))
  · norm_num [f9, peval]
  · norm_num [f9, peval]

theorem N9_A4_neg (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : 0 < peval f9 x) :
    aBlock 9 4 x < 0 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 9 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : (1 - x) ^ 4 * peval [-4, 24, 80, 160, 224, 259, 306, 361, 420, 361, 306, 259, 224, 160, 80, 24, -4] x < 0 := by
    rw [N9_crit_factor]
    have := mul_pos (N9_posH x h0 h1) hf
    nlinarith [pow_pos h1x 4]
  rw [← N9_linkA4 x h0 h1] at hn
  have hd := (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))
  by_contra hc
  push Not at hc
  nlinarith [mul_nonneg hc hd.le]

theorem N9_A4_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : peval f9 x < 0) :
    0 < aBlock 9 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 9 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : 0 < (1 - x) ^ 4 * peval [-4, 24, 80, 160, 224, 259, 306, 361, 420, 361, 306, 259, 224, 160, 80, 24, -4] x := by
    rw [N9_crit_factor]
    have := mul_pos (N9_posH x h0 h1) (neg_pos.mpr hf)
    have h4 := pow_pos h1x 4
    nlinarith
  rw [← N9_linkA4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

/-- Note Thm 6(b) for `N = 9`: there is a unique `ξ* ∈ (0,1)` (the unique root of `f_9` in
`(0,1)`, enclosed in `(11287106323660274043/100000000000000000000, 2821776580915068511/25000000000000000000)`) such that, for every `ξ ∈ (0,1)`: all `b_m > 0`;
`a_4 < 0` for `ξ < ξ*` (spectrally unstable, since `b_4 > 0`); and all `a_m > 0` for `ξ > ξ*`. -/
theorem cert_N9 : ∃ ξs : ℝ, ((11287106323660274043 : ℝ) / 100000000000000000000) < ξs ∧ ξs < ((2821776580915068511 : ℝ) / 25000000000000000000) ∧ peval f9 ξs = 0 ∧
    (∀ y, 0 < y → y < 1 → peval f9 y = 0 → y = ξs) ∧
    ∀ x, 0 < x → x < 1 →
      (∀ m, 1 ≤ m → m ≤ 8 → 0 < bBlock 9 m x) ∧
      (x < ξs → aBlock 9 4 x < 0) ∧
      (ξs < x → ∀ m, 1 ≤ m → m ≤ 8 → 0 < aBlock 9 m x) := by
  obtain ⟨ξs, hta, htb, hf0, hs⟩ := N9_root
  refine ⟨ξs, hta, htb, hf0, ?_, ?_⟩
  · intro y hy0 hy1 hy
    rcases lt_trichotomy y ξs with h | h | h
    · linarith [(hs y hy0 hy1).1 h]
    · exact h
    · linarith [(hs y hy0 hy1).2 h]
  intro x h0 h1
  refine ⟨?_, fun hlt => N9_A4_neg x h0 h1 ((hs x h0 h1).1 hlt), ?_⟩
  · intro m hm1 hm2
    interval_cases m
    · exact N9_B1_pos x h0 h1
    · exact N9_B2_pos x h0 h1
    · exact N9_B3_pos x h0 h1
    · exact N9_B4_pos x h0 h1
    · rw [show bBlock 9 5 x = bBlock 9 4 x from bBlock_symm 9 4 (by norm_num) x]
      exact N9_B4_pos x h0 h1
    · rw [show bBlock 9 6 x = bBlock 9 3 x from bBlock_symm 9 3 (by norm_num) x]
      exact N9_B3_pos x h0 h1
    · rw [show bBlock 9 7 x = bBlock 9 2 x from bBlock_symm 9 2 (by norm_num) x]
      exact N9_B2_pos x h0 h1
    · rw [show bBlock 9 8 x = bBlock 9 1 x from bBlock_symm 9 1 (by norm_num) x]
      exact N9_B1_pos x h0 h1
  · intro hgt m hm1 hm2
    have hlo : ((1 : ℝ) / 16) < x := lt_trans (lt_trans (by norm_num) hta) hgt
    have hfn := (hs x h0 h1).2 hgt
    interval_cases m
    · exact N9_A1_pos x h0 h1
    · exact N9_A2_pos x h0 h1
    · exact N9_A3_pos x h0 hlo h1
    · exact N9_A4_pos x h0 h1 hfn
    · rw [show aBlock 9 5 x = aBlock 9 4 x from aBlock_symm 9 4 (by norm_num) x]
      exact N9_A4_pos x h0 h1 hfn
    · rw [show aBlock 9 6 x = aBlock 9 3 x from aBlock_symm 9 3 (by norm_num) x]
      exact N9_A3_pos x h0 hlo h1
    · rw [show aBlock 9 7 x = aBlock 9 2 x from aBlock_symm 9 2 (by norm_num) x]
      exact N9_A2_pos x h0 h1
    · rw [show aBlock 9 8 x = aBlock 9 1 x from aBlock_symm 9 1 (by norm_num) x]
      exact N9_A1_pos x h0 h1

/-! ### N = 10 -/

theorem N10_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 10 1 x * (2 * (1 - x ^ 10) ^ 2) = (1 - x) ^ 4 * peval [9, 26, 50, 80, 115, 154, 196, 240, 285, 240, 196, 154, 115, 80, 50, 26, 9] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 10 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (10 - 1 : ℕ) = 9 by norm_num, show (2 * 10 - 1 : ℕ) = 19 by norm_num, show (10 + 1 : ℕ) = 11 by norm_num, show (10 - 1 + 10 : ℕ) = 19 by norm_num, show (1 + 10 : ℕ) = 11 by norm_num, show (10 + 10 : ℕ) = 20 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N10_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [9, 26, 50, 80, 115, 154, 196, 240, 285, 240, 196, 154, 115, 80, 50, 26, 9] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 9 * q ^ 16 + 170 * p ^ 1 * q ^ 15 + 1520 * p ^ 2 * q ^ 14 + 8550 * p ^ 3 * q ^ 13 + 33915 * p ^ 4 * q ^ 12 + 100776 * p ^ 5 * q ^ 11 + 232560 * p ^ 6 * q ^ 10 + 426360 * p ^ 7 * q ^ 9 + 629850 * p ^ 8 * q ^ 8 + 755730 * p ^ 9 * q ^ 7 + 738216 * p ^ 10 * q ^ 6 + 584640 * p ^ 11 * q ^ 5 + 370440 * p ^ 12 * q ^ 4 + 182700 * p ^ 13 * q ^ 3 + 66600 * p ^ 14 * q ^ 2 + 16200 * p ^ 15 * q ^ 1 + 2025 * p ^ 16 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [9, 26, 50, 80, 115, 154, 196, 240, 285, 240, 196, 154, 115, 80, 50, 26, 9] x =
      (fun p q : ℝ => 9 * q ^ 16 + 170 * p ^ 1 * q ^ 15 + 1520 * p ^ 2 * q ^ 14 + 8550 * p ^ 3 * q ^ 13 + 33915 * p ^ 4 * q ^ 12 + 100776 * p ^ 5 * q ^ 11 + 232560 * p ^ 6 * q ^ 10 + 426360 * p ^ 7 * q ^ 9 + 629850 * p ^ 8 * q ^ 8 + 755730 * p ^ 9 * q ^ 7 + 738216 * p ^ 10 * q ^ 6 + 584640 * p ^ 11 * q ^ 5 + 370440 * p ^ 12 * q ^ 4 + 182700 * p ^ 13 * q ^ 3 + 66600 * p ^ 14 * q ^ 2 + 16200 * p ^ 15 * q ^ 1 + 2025 * p ^ 16) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [9, 26, 50, 80, 115, 154, 196, 240, 285, 240, 196, 154, 115, 80, 50, 26, 9] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N10_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [9, 26, 50, 80, 115, 154, 196, 240, 285, 240, 196, 154, 115, 80, 50, 26, 9] x := N10_posB1_0 x hu hv

theorem N10_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 10 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 10 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N10_posB1 x h0 h1)
  rw [← N10_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N10_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 10 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 10) ^ 2) = (1 - x) ^ 4 * peval [9, 44, 111, 206, 325, 464, 619, 786, 961, 1050, 961, 786, 619, 464, 325, 206, 111, 44, 9] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 10 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (10 - 1 : ℕ) = 9 by norm_num, show (2 * 10 - 1 : ℕ) = 19 by norm_num, show (10 + 1 : ℕ) = 11 by norm_num, show (10 - 1 + 10 : ℕ) = 19 by norm_num, show (1 + 10 : ℕ) = 11 by norm_num, show (10 + 10 : ℕ) = 20 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N10_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [9, 44, 111, 206, 325, 464, 619, 786, 961, 1050, 961, 786, 619, 464, 325, 206, 111, 44, 9] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 9 * q ^ 18 + 206 * p ^ 1 * q ^ 17 + 2236 * p ^ 2 * q ^ 16 + 15310 * p ^ 3 * q ^ 15 + 74195 * p ^ 4 * q ^ 14 + 270636 * p ^ 5 * q ^ 13 + 771324 * p ^ 6 * q ^ 12 + 1759704 * p ^ 7 * q ^ 11 + 3265530 * p ^ 8 * q ^ 10 + 4980570 * p ^ 9 * q ^ 9 + 6280536 * p ^ 10 * q ^ 8 + 6560424 * p ^ 11 * q ^ 7 + 5661864 * p ^ 12 * q ^ 6 + 4003020 * p ^ 13 * q ^ 5 + 2279160 * p ^ 14 * q ^ 4 + 1013400 * p ^ 15 * q ^ 3 + 333225 * p ^ 16 * q ^ 2 + 72900 * p ^ 17 * q ^ 1 + 8100 * p ^ 18 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [9, 44, 111, 206, 325, 464, 619, 786, 961, 1050, 961, 786, 619, 464, 325, 206, 111, 44, 9] x =
      (fun p q : ℝ => 9 * q ^ 18 + 206 * p ^ 1 * q ^ 17 + 2236 * p ^ 2 * q ^ 16 + 15310 * p ^ 3 * q ^ 15 + 74195 * p ^ 4 * q ^ 14 + 270636 * p ^ 5 * q ^ 13 + 771324 * p ^ 6 * q ^ 12 + 1759704 * p ^ 7 * q ^ 11 + 3265530 * p ^ 8 * q ^ 10 + 4980570 * p ^ 9 * q ^ 9 + 6280536 * p ^ 10 * q ^ 8 + 6560424 * p ^ 11 * q ^ 7 + 5661864 * p ^ 12 * q ^ 6 + 4003020 * p ^ 13 * q ^ 5 + 2279160 * p ^ 14 * q ^ 4 + 1013400 * p ^ 15 * q ^ 3 + 333225 * p ^ 16 * q ^ 2 + 72900 * p ^ 17 * q ^ 1 + 8100 * p ^ 18) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [9, 44, 111, 206, 325, 464, 619, 786, 961, 1050, 961, 786, 619, 464, 325, 206, 111, 44, 9] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N10_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [9, 44, 111, 206, 325, 464, 619, 786, 961, 1050, 961, 786, 619, 464, 325, 206, 111, 44, 9] x := N10_posA1_0 x hu hv

theorem N10_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 10 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 10 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N10_posA1 x h0 h1)
  rw [← N10_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N10_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 10 2 x * (2 * (1 - x ^ 10) ^ 2) = (1 - x) ^ 4 * peval [16, 64, 140, 240, 360, 496, 644, 800, 880, 800, 644, 496, 360, 240, 140, 64, 16] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 10 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (10 - 2 : ℕ) = 8 by norm_num, show (2 * 10 - 2 : ℕ) = 18 by norm_num, show (10 + 2 : ℕ) = 12 by norm_num, show (10 - 2 + 10 : ℕ) = 18 by norm_num, show (2 + 10 : ℕ) = 12 by norm_num, show (10 + 10 : ℕ) = 20 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N10_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [16, 64, 140, 240, 360, 496, 644, 800, 880, 800, 644, 496, 360, 240, 140, 64, 16] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 16 * q ^ 16 + 320 * p ^ 1 * q ^ 15 + 3020 * p ^ 2 * q ^ 14 + 17880 * p ^ 3 * q ^ 13 + 74460 * p ^ 4 * q ^ 12 + 231744 * p ^ 5 * q ^ 11 + 558960 * p ^ 6 * q ^ 10 + 1068960 * p ^ 7 * q ^ 9 + 1644160 * p ^ 8 * q ^ 8 + 2049920 * p ^ 9 * q ^ 7 + 2075824 * p ^ 10 * q ^ 6 + 1699040 * p ^ 11 * q ^ 5 + 1108240 * p ^ 12 * q ^ 4 + 560000 * p ^ 13 * q ^ 3 + 208000 * p ^ 14 * q ^ 2 + 51200 * p ^ 15 * q ^ 1 + 6400 * p ^ 16 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [16, 64, 140, 240, 360, 496, 644, 800, 880, 800, 644, 496, 360, 240, 140, 64, 16] x =
      (fun p q : ℝ => 16 * q ^ 16 + 320 * p ^ 1 * q ^ 15 + 3020 * p ^ 2 * q ^ 14 + 17880 * p ^ 3 * q ^ 13 + 74460 * p ^ 4 * q ^ 12 + 231744 * p ^ 5 * q ^ 11 + 558960 * p ^ 6 * q ^ 10 + 1068960 * p ^ 7 * q ^ 9 + 1644160 * p ^ 8 * q ^ 8 + 2049920 * p ^ 9 * q ^ 7 + 2075824 * p ^ 10 * q ^ 6 + 1699040 * p ^ 11 * q ^ 5 + 1108240 * p ^ 12 * q ^ 4 + 560000 * p ^ 13 * q ^ 3 + 208000 * p ^ 14 * q ^ 2 + 51200 * p ^ 15 * q ^ 1 + 6400 * p ^ 16) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [16, 64, 140, 240, 360, 496, 644, 800, 880, 800, 644, 496, 360, 240, 140, 64, 16] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N10_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [16, 64, 140, 240, 360, 496, 644, 800, 880, 800, 644, 496, 360, 240, 140, 64, 16] x := N10_posB2_0 x hu hv

theorem N10_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 10 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 10 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N10_posB2 x h0 h1)
  rw [← N10_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N10_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 10 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 10) ^ 2) = (1 - x) ^ 4 * peval [2, 40, 90, 168, 270, 392, 530, 680, 758, 840, 758, 680, 530, 392, 270, 168, 90, 40, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 10 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (10 - 2 : ℕ) = 8 by norm_num, show (2 * 10 - 2 : ℕ) = 18 by norm_num, show (10 + 2 : ℕ) = 12 by norm_num, show (10 - 2 + 10 : ℕ) = 18 by norm_num, show (2 + 10 : ℕ) = 12 by norm_num, show (10 + 10 : ℕ) = 20 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N10_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 40, 90, 168, 270, 392, 530, 680, 758, 840, 758, 680, 530, 392, 270, 168, 90, 40, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 18 + 76 * p ^ 1 * q ^ 17 + 1076 * p ^ 2 * q ^ 16 + 8680 * p ^ 3 * q ^ 15 + 46910 * p ^ 4 * q ^ 14 + 184548 * p ^ 5 * q ^ 13 + 555084 * p ^ 6 * q ^ 12 + 1317024 * p ^ 7 * q ^ 11 + 2516260 * p ^ 8 * q ^ 10 + 3923320 * p ^ 9 * q ^ 9 + 5031588 * p ^ 10 * q ^ 8 + 5323664 * p ^ 11 * q ^ 7 + 4637344 * p ^ 12 * q ^ 6 + 3298400 * p ^ 13 * q ^ 5 + 1883800 * p ^ 14 * q ^ 4 + 838400 * p ^ 15 * q ^ 3 + 275650 * p ^ 16 * q ^ 2 + 60300 * p ^ 17 * q ^ 1 + 6700 * p ^ 18 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 40, 90, 168, 270, 392, 530, 680, 758, 840, 758, 680, 530, 392, 270, 168, 90, 40, 2] x =
      (fun p q : ℝ => 2 * q ^ 18 + 76 * p ^ 1 * q ^ 17 + 1076 * p ^ 2 * q ^ 16 + 8680 * p ^ 3 * q ^ 15 + 46910 * p ^ 4 * q ^ 14 + 184548 * p ^ 5 * q ^ 13 + 555084 * p ^ 6 * q ^ 12 + 1317024 * p ^ 7 * q ^ 11 + 2516260 * p ^ 8 * q ^ 10 + 3923320 * p ^ 9 * q ^ 9 + 5031588 * p ^ 10 * q ^ 8 + 5323664 * p ^ 11 * q ^ 7 + 4637344 * p ^ 12 * q ^ 6 + 3298400 * p ^ 13 * q ^ 5 + 1883800 * p ^ 14 * q ^ 4 + 838400 * p ^ 15 * q ^ 3 + 275650 * p ^ 16 * q ^ 2 + 60300 * p ^ 17 * q ^ 1 + 6700 * p ^ 18) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 40, 90, 168, 270, 392, 530, 680, 758, 840, 758, 680, 530, 392, 270, 168, 90, 40, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N10_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 40, 90, 168, 270, 392, 530, 680, 758, 840, 758, 680, 530, 392, 270, 168, 90, 40, 2] x := N10_posA2_0 x hu hv

theorem N10_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 10 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 10 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N10_posA2 x h0 h1)
  rw [← N10_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N10_linkB3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 10 3 x * (2 * (1 - x ^ 10) ^ 2) = (1 - x) ^ 4 * peval [21, 84, 210, 390, 615, 876, 1164, 1400, 1505, 1400, 1164, 876, 615, 390, 210, 84, 21] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 10 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (10 - 3 : ℕ) = 7 by norm_num, show (2 * 10 - 3 : ℕ) = 17 by norm_num, show (10 + 3 : ℕ) = 13 by norm_num, show (10 - 3 + 10 : ℕ) = 17 by norm_num, show (3 + 10 : ℕ) = 13 by norm_num, show (10 + 10 : ℕ) = 20 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N10_posB3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [21, 84, 210, 390, 615, 876, 1164, 1400, 1505, 1400, 1164, 876, 615, 390, 210, 84, 21] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 21 * q ^ 16 + 420 * p ^ 1 * q ^ 15 + 3990 * p ^ 2 * q ^ 14 + 23910 * p ^ 3 * q ^ 13 + 101235 * p ^ 4 * q ^ 12 + 321504 * p ^ 5 * q ^ 11 + 793560 * p ^ 6 * q ^ 10 + 1556450 * p ^ 7 * q ^ 9 + 2458820 * p ^ 8 * q ^ 8 + 3150420 * p ^ 9 * q ^ 7 + 3276574 * p ^ 10 * q ^ 6 + 2749390 * p ^ 11 * q ^ 5 + 1833090 * p ^ 12 * q ^ 4 + 943250 * p ^ 13 * q ^ 3 + 355250 * p ^ 14 * q ^ 2 + 88200 * p ^ 15 * q ^ 1 + 11025 * p ^ 16 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [21, 84, 210, 390, 615, 876, 1164, 1400, 1505, 1400, 1164, 876, 615, 390, 210, 84, 21] x =
      (fun p q : ℝ => 21 * q ^ 16 + 420 * p ^ 1 * q ^ 15 + 3990 * p ^ 2 * q ^ 14 + 23910 * p ^ 3 * q ^ 13 + 101235 * p ^ 4 * q ^ 12 + 321504 * p ^ 5 * q ^ 11 + 793560 * p ^ 6 * q ^ 10 + 1556450 * p ^ 7 * q ^ 9 + 2458820 * p ^ 8 * q ^ 8 + 3150420 * p ^ 9 * q ^ 7 + 3276574 * p ^ 10 * q ^ 6 + 2749390 * p ^ 11 * q ^ 5 + 1833090 * p ^ 12 * q ^ 4 + 943250 * p ^ 13 * q ^ 3 + 355250 * p ^ 14 * q ^ 2 + 88200 * p ^ 15 * q ^ 1 + 11025 * p ^ 16) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [21, 84, 210, 390, 615, 876, 1164, 1400, 1505, 1400, 1164, 876, 615, 390, 210, 84, 21] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N10_posB3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [21, 84, 210, 390, 615, 876, 1164, 1400, 1505, 1400, 1164, 876, 615, 390, 210, 84, 21] x := N10_posB3_0 x hu hv

theorem N10_B3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 10 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 10 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N10_posB3 x h0 h1)
  rw [← N10_linkB3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N10_linkA3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 10 3 x * (2 * (1 - x) ^ 2 * (1 - x ^ 10) ^ 2) = (1 - x) ^ 4 * peval [-3, 30, 95, 158, 245, 352, 475, 540, 613, 690, 613, 540, 475, 352, 245, 158, 95, 30, -3] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 10 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (10 - 3 : ℕ) = 7 by norm_num, show (2 * 10 - 3 : ℕ) = 17 by norm_num, show (10 + 3 : ℕ) = 13 by norm_num, show (10 - 3 + 10 : ℕ) = 17 by norm_num, show (3 + 10 : ℕ) = 13 by norm_num, show (10 + 10 : ℕ) = 20 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N10_posA3_0 (x : ℝ) (hu : ((5 : ℝ) / 32) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-3, 30, 95, 158, 245, 352, 475, 540, 613, 690, 613, 540, 475, 352, 245, 158, 95, 30, -3] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 1979372973532663971586203751 * q ^ 18 + 62222474538222592405938801408 * p ^ 1 * q ^ 17 + 820742213924601629665770805248 * p ^ 2 * q ^ 16 + 6396033568752914597294512865280 * p ^ 3 * q ^ 15 + 33873588976740990197197968506880 * p ^ 4 * q ^ 14 + 131357303276466294883104004767744 * p ^ 5 * q ^ 13 + 390150437236793173966683868495872 * p ^ 6 * q ^ 12 + 913559044289036987895837131538432 * p ^ 7 * q ^ 11 + 1718517233453970709588328432271360 * p ^ 8 * q ^ 10 + 2627992697699621003809940842741760 * p ^ 9 * q ^ 9 + 3287653268954474852223905100201984 * p ^ 10 * q ^ 8 + 3369113864448966599892628156710912 * p ^ 11 * q ^ 7 + 2817371232142007989401089759772672 * p ^ 12 * q ^ 6 + 1903335270058767496163709734092800 * p ^ 13 * q ^ 5 + 1019757185774115350781371822899200 * p ^ 14 * q ^ 4 + 419771360668514823014110671667200 * p ^ 15 * q ^ 3 + 125566379849529544957219779379200 * p ^ 16 * q ^ 2 + 24476395714245628122771121766400 * p ^ 17 * q ^ 1 + 2352086074642222522308336025600 * p ^ 18 := by
    intro p q hp hq; positivity
  have key : (19383245667680019896796723 : ℝ) * peval [-3, 30, 95, 158, 245, 352, 475, 540, 613, 690, 613, 540, 475, 352, 245, 158, 95, 30, -3] x =
      (fun p q : ℝ => 1979372973532663971586203751 * q ^ 18 + 62222474538222592405938801408 * p ^ 1 * q ^ 17 + 820742213924601629665770805248 * p ^ 2 * q ^ 16 + 6396033568752914597294512865280 * p ^ 3 * q ^ 15 + 33873588976740990197197968506880 * p ^ 4 * q ^ 14 + 131357303276466294883104004767744 * p ^ 5 * q ^ 13 + 390150437236793173966683868495872 * p ^ 6 * q ^ 12 + 913559044289036987895837131538432 * p ^ 7 * q ^ 11 + 1718517233453970709588328432271360 * p ^ 8 * q ^ 10 + 2627992697699621003809940842741760 * p ^ 9 * q ^ 9 + 3287653268954474852223905100201984 * p ^ 10 * q ^ 8 + 3369113864448966599892628156710912 * p ^ 11 * q ^ 7 + 2817371232142007989401089759772672 * p ^ 12 * q ^ 6 + 1903335270058767496163709734092800 * p ^ 13 * q ^ 5 + 1019757185774115350781371822899200 * p ^ 14 * q ^ 4 + 419771360668514823014110671667200 * p ^ 15 * q ^ 3 + 125566379849529544957219779379200 * p ^ 16 * q ^ 2 + 24476395714245628122771121766400 * p ^ 17 * q ^ 1 + 2352086074642222522308336025600 * p ^ 18) (x - ((5 : ℝ) / 32)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (19383245667680019896796723 : ℝ) * peval [-3, 30, 95, 158, 245, 352, 475, 540, 613, 690, 613, 540, 475, 352, 245, 158, 95, 30, -3] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N10_posA3 (x : ℝ) (hu : ((5 : ℝ) / 32) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-3, 30, 95, 158, 245, 352, 475, 540, 613, 690, 613, 540, 475, 352, 245, 158, 95, 30, -3] x := N10_posA3_0 x hu hv

theorem N10_A3_pos (x : ℝ) (h0 : 0 < x) (hlo : ((5 : ℝ) / 32) < x) (h1 : x < 1) :
    0 < aBlock 10 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 10 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N10_posA3 x hlo h1)
  rw [← N10_linkA3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N10_linkB4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 10 4 x * (2 * (1 - x ^ 10) ^ 2) = (1 - x) ^ 4 * peval [24, 96, 240, 480, 800, 1184, 1556, 1840, 1960, 1840, 1556, 1184, 800, 480, 240, 96, 24] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 10 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (10 - 4 : ℕ) = 6 by norm_num, show (2 * 10 - 4 : ℕ) = 16 by norm_num, show (10 + 4 : ℕ) = 14 by norm_num, show (10 - 4 + 10 : ℕ) = 16 by norm_num, show (4 + 10 : ℕ) = 14 by norm_num, show (10 + 10 : ℕ) = 20 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N10_posB4_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [24, 96, 240, 480, 800, 1184, 1556, 1840, 1960, 1840, 1556, 1184, 800, 480, 240, 96, 24] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 24 * q ^ 16 + 480 * p ^ 1 * q ^ 15 + 4560 * p ^ 2 * q ^ 14 + 27360 * p ^ 3 * q ^ 13 + 116240 * p ^ 4 * q ^ 12 + 371456 * p ^ 5 * q ^ 11 + 925380 * p ^ 6 * q ^ 10 + 1837240 * p ^ 7 * q ^ 9 + 2945020 * p ^ 8 * q ^ 8 + 3834480 * p ^ 9 * q ^ 7 + 4053916 * p ^ 10 * q ^ 6 + 3454840 * p ^ 11 * q ^ 5 + 2335140 * p ^ 12 * q ^ 4 + 1215200 * p ^ 13 * q ^ 3 + 461600 * p ^ 14 * q ^ 2 + 115200 * p ^ 15 * q ^ 1 + 14400 * p ^ 16 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [24, 96, 240, 480, 800, 1184, 1556, 1840, 1960, 1840, 1556, 1184, 800, 480, 240, 96, 24] x =
      (fun p q : ℝ => 24 * q ^ 16 + 480 * p ^ 1 * q ^ 15 + 4560 * p ^ 2 * q ^ 14 + 27360 * p ^ 3 * q ^ 13 + 116240 * p ^ 4 * q ^ 12 + 371456 * p ^ 5 * q ^ 11 + 925380 * p ^ 6 * q ^ 10 + 1837240 * p ^ 7 * q ^ 9 + 2945020 * p ^ 8 * q ^ 8 + 3834480 * p ^ 9 * q ^ 7 + 4053916 * p ^ 10 * q ^ 6 + 3454840 * p ^ 11 * q ^ 5 + 2335140 * p ^ 12 * q ^ 4 + 1215200 * p ^ 13 * q ^ 3 + 461600 * p ^ 14 * q ^ 2 + 115200 * p ^ 15 * q ^ 1 + 14400 * p ^ 16) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [24, 96, 240, 480, 800, 1184, 1556, 1840, 1960, 1840, 1556, 1184, 800, 480, 240, 96, 24] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N10_posB4 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [24, 96, 240, 480, 800, 1184, 1556, 1840, 1960, 1840, 1556, 1184, 800, 480, 240, 96, 24] x := N10_posB4_0 x hu hv

theorem N10_B4_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 10 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 10 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N10_posB4 x h0 h1)
  rw [← N10_linkB4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N10_linkA4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 10 4 x * (2 * (1 - x) ^ 2 * (1 - x ^ 10) ^ 2) = (1 - x) ^ 4 * peval [-6, 24, 86, 176, 250, 344, 394, 456, 526, 600, 526, 456, 394, 344, 250, 176, 86, 24, -6] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 10 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (10 - 4 : ℕ) = 6 by norm_num, show (2 * 10 - 4 : ℕ) = 16 by norm_num, show (10 + 4 : ℕ) = 14 by norm_num, show (10 - 4 + 10 : ℕ) = 16 by norm_num, show (4 + 10 : ℕ) = 14 by norm_num, show (10 + 10 : ℕ) = 20 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N10_posA4_0 (x : ℝ) (hu : ((5 : ℝ) / 32) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-6, 24, 86, 176, 250, 344, 394, 456, 526, 600, 526, 456, 394, 344, 250, 176, 86, 24, -6] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 292583388395686089870850382 * q ^ 18 + 29249875196858034470252838528 * p ^ 1 * q ^ 17 + 517973014284352116128783474688 * p ^ 2 * q ^ 16 + 4662218448339863161145941032960 * p ^ 3 * q ^ 15 + 26929012222194893509469782671360 * p ^ 4 * q ^ 14 + 110628068214984030767157438578688 * p ^ 5 * q ^ 13 + 342018834232739977175607753572352 * p ^ 6 * q ^ 12 + 823542481026919332859239771144192 * p ^ 7 * q ^ 11 + 1578428349211420528956850726502400 * p ^ 8 * q ^ 10 + 2441003290173038941039130061045760 * p ^ 9 * q ^ 9 + 3068955894167525470530609975656448 * p ^ 10 * q ^ 8 + 3144280116186682938438068880801792 * p ^ 11 * q ^ 7 + 2617834566525794333544047425093632 * p ^ 12 * q ^ 6 + 1755476863456152713734410479861760 * p ^ 13 * q ^ 5 + 932015483788481425452347298938880 * p ^ 14 * q ^ 4 + 380122182783681967014110507827200 * p ^ 15 * q ^ 3 + 112839383060391045954771143884800 * p ^ 16 * q ^ 2 + 21899933007482930425637319475200 * p ^ 17 * q ^ 1 + 2104498066785146467328511180800 * p ^ 18 := by
    intro p q hp hq; positivity
  have key : (19383245667680019896796723 : ℝ) * peval [-6, 24, 86, 176, 250, 344, 394, 456, 526, 600, 526, 456, 394, 344, 250, 176, 86, 24, -6] x =
      (fun p q : ℝ => 292583388395686089870850382 * q ^ 18 + 29249875196858034470252838528 * p ^ 1 * q ^ 17 + 517973014284352116128783474688 * p ^ 2 * q ^ 16 + 4662218448339863161145941032960 * p ^ 3 * q ^ 15 + 26929012222194893509469782671360 * p ^ 4 * q ^ 14 + 110628068214984030767157438578688 * p ^ 5 * q ^ 13 + 342018834232739977175607753572352 * p ^ 6 * q ^ 12 + 823542481026919332859239771144192 * p ^ 7 * q ^ 11 + 1578428349211420528956850726502400 * p ^ 8 * q ^ 10 + 2441003290173038941039130061045760 * p ^ 9 * q ^ 9 + 3068955894167525470530609975656448 * p ^ 10 * q ^ 8 + 3144280116186682938438068880801792 * p ^ 11 * q ^ 7 + 2617834566525794333544047425093632 * p ^ 12 * q ^ 6 + 1755476863456152713734410479861760 * p ^ 13 * q ^ 5 + 932015483788481425452347298938880 * p ^ 14 * q ^ 4 + 380122182783681967014110507827200 * p ^ 15 * q ^ 3 + 112839383060391045954771143884800 * p ^ 16 * q ^ 2 + 21899933007482930425637319475200 * p ^ 17 * q ^ 1 + 2104498066785146467328511180800 * p ^ 18) (x - ((5 : ℝ) / 32)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (19383245667680019896796723 : ℝ) * peval [-6, 24, 86, 176, 250, 344, 394, 456, 526, 600, 526, 456, 394, 344, 250, 176, 86, 24, -6] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N10_posA4 (x : ℝ) (hu : ((5 : ℝ) / 32) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-6, 24, 86, 176, 250, 344, 394, 456, 526, 600, 526, 456, 394, 344, 250, 176, 86, 24, -6] x := N10_posA4_0 x hu hv

theorem N10_A4_pos (x : ℝ) (h0 : 0 < x) (hlo : ((5 : ℝ) / 32) < x) (h1 : x < 1) :
    0 < aBlock 10 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 10 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N10_posA4 x hlo h1)
  rw [← N10_linkA4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N10_linkB5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 10 5 x * (2 * (1 - x ^ 10) ^ 2) = (1 - x) ^ 4 * peval [25, 100, 250, 500, 875, 1300, 1700, 2000, 2125, 2000, 1700, 1300, 875, 500, 250, 100, 25] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 10 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (10 - 5 : ℕ) = 5 by norm_num, show (2 * 10 - 5 : ℕ) = 15 by norm_num, show (10 + 5 : ℕ) = 15 by norm_num, show (10 - 5 + 10 : ℕ) = 15 by norm_num, show (5 + 10 : ℕ) = 15 by norm_num, show (10 + 10 : ℕ) = 20 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N10_posB5_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [25, 100, 250, 500, 875, 1300, 1700, 2000, 2125, 2000, 1700, 1300, 875, 500, 250, 100, 25] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 25 * q ^ 16 + 500 * p ^ 1 * q ^ 15 + 4750 * p ^ 2 * q ^ 14 + 28500 * p ^ 3 * q ^ 13 + 121125 * p ^ 4 * q ^ 12 + 387500 * p ^ 5 * q ^ 11 + 967500 * p ^ 6 * q ^ 10 + 1927500 * p ^ 7 * q ^ 9 + 3103750 * p ^ 8 * q ^ 8 + 4062500 * p ^ 9 * q ^ 7 + 4318750 * p ^ 10 * q ^ 6 + 3700000 * p ^ 11 * q ^ 5 + 2512500 * p ^ 12 * q ^ 4 + 1312500 * p ^ 13 * q ^ 3 + 500000 * p ^ 14 * q ^ 2 + 125000 * p ^ 15 * q ^ 1 + 15625 * p ^ 16 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [25, 100, 250, 500, 875, 1300, 1700, 2000, 2125, 2000, 1700, 1300, 875, 500, 250, 100, 25] x =
      (fun p q : ℝ => 25 * q ^ 16 + 500 * p ^ 1 * q ^ 15 + 4750 * p ^ 2 * q ^ 14 + 28500 * p ^ 3 * q ^ 13 + 121125 * p ^ 4 * q ^ 12 + 387500 * p ^ 5 * q ^ 11 + 967500 * p ^ 6 * q ^ 10 + 1927500 * p ^ 7 * q ^ 9 + 3103750 * p ^ 8 * q ^ 8 + 4062500 * p ^ 9 * q ^ 7 + 4318750 * p ^ 10 * q ^ 6 + 3700000 * p ^ 11 * q ^ 5 + 2512500 * p ^ 12 * q ^ 4 + 1312500 * p ^ 13 * q ^ 3 + 500000 * p ^ 14 * q ^ 2 + 125000 * p ^ 15 * q ^ 1 + 15625 * p ^ 16) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [25, 100, 250, 500, 875, 1300, 1700, 2000, 2125, 2000, 1700, 1300, 875, 500, 250, 100, 25] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N10_posB5 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [25, 100, 250, 500, 875, 1300, 1700, 2000, 2125, 2000, 1700, 1300, 875, 500, 250, 100, 25] x := N10_posB5_0 x hu hv

theorem N10_B5_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 10 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 10 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N10_posB5 x h0 h1)
  rw [← N10_linkB5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N10_linkA5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 10 5 x * (2 * (1 - x) ^ 2 * (1 - x ^ 10) ^ 2) = (1 - x) ^ 4 * peval [-7, 22, 83, 172, 285, 318, 367, 428, 497, 570, 497, 428, 367, 318, 285, 172, 83, 22, -7] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 10 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (10 - 5 : ℕ) = 5 by norm_num, show (2 * 10 - 5 : ℕ) = 15 by norm_num, show (10 + 5 : ℕ) = 15 by norm_num, show (10 - 5 + 10 : ℕ) = 15 by norm_num, show (5 + 10 : ℕ) = 15 by norm_num, show (10 + 10 : ℕ) = 20 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

/-- The minimal polynomial `f_10` of `ξ*(10)` (constant term first). -/
def f10 : List ℤ := [7, -29, -54, -118, -167, -158, -187, -158, -167, -118, -54, -29, 7]

theorem N10_crit_factor (x : ℝ) : peval [-7, 22, 83, 172, 285, 318, 367, 428, 497, 570, 497, 428, 367, 318, 285, 172, 83, 22, -7] x = -(peval [1, 1, 0, 0, 0, 1, 1] x * peval f10 x) := by
  simp only [f10, peval]; ring

theorem N10_posH_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1, 1, 0, 0, 0, 1, 1] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 1 * q ^ 6 + 7 * p ^ 1 * q ^ 5 + 20 * p ^ 2 * q ^ 4 + 30 * p ^ 3 * q ^ 3 + 25 * p ^ 4 * q ^ 2 + 12 * p ^ 5 * q ^ 1 + 4 * p ^ 6 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [1, 1, 0, 0, 0, 1, 1] x =
      (fun p q : ℝ => 1 * q ^ 6 + 7 * p ^ 1 * q ^ 5 + 20 * p ^ 2 * q ^ 4 + 30 * p ^ 3 * q ^ 3 + 25 * p ^ 4 * q ^ 2 + 12 * p ^ 5 * q ^ 1 + 4 * p ^ 6) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [1, 1, 0, 0, 0, 1, 1] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N10_posH (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1, 1, 0, 0, 0, 1, 1] x := N10_posH_0 x hu hv

theorem N10_dneg_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [29, 108, 354, 668, 790, 1122, 1106, 1336, 1062, 540, 319, -84] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 29 * q ^ 11 + 427 * p ^ 1 * q ^ 10 + 3029 * p ^ 2 * q ^ 9 + 13499 * p ^ 3 * q ^ 8 + 41408 * p ^ 4 * q ^ 7 + 91170 * p ^ 5 * q ^ 6 + 147054 * p ^ 6 * q ^ 5 + 174960 * p ^ 7 * q ^ 4 + 152445 * p ^ 8 * q ^ 3 + 94125 * p ^ 9 * q ^ 2 + 37650 * p ^ 10 * q ^ 1 + 7350 * p ^ 11 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [29, 108, 354, 668, 790, 1122, 1106, 1336, 1062, 540, 319, -84] x =
      (fun p q : ℝ => 29 * q ^ 11 + 427 * p ^ 1 * q ^ 10 + 3029 * p ^ 2 * q ^ 9 + 13499 * p ^ 3 * q ^ 8 + 41408 * p ^ 4 * q ^ 7 + 91170 * p ^ 5 * q ^ 6 + 147054 * p ^ 6 * q ^ 5 + 174960 * p ^ 7 * q ^ 4 + 152445 * p ^ 8 * q ^ 3 + 94125 * p ^ 9 * q ^ 2 + 37650 * p ^ 10 * q ^ 1 + 7350 * p ^ 11) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [29, 108, 354, 668, 790, 1122, 1106, 1336, 1062, 540, 319, -84] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N10_dneg (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [29, 108, 354, 668, 790, 1122, 1106, 1336, 1062, 540, 319, -84] x := N10_dneg_0 x hu hv

theorem N10_dlink (x : ℝ) : dpeval f10 x = -peval [29, 108, 354, 668, 790, 1122, 1106, 1336, 1062, 540, 319, -84] x := by
  simp only [f10, dpeval, peval]; ring

/-- `f_10` changes sign exactly once in `(0,1)`, at a point of `(129776751192374311/781250000000000000, 16611424152623911809/100000000000000000000)`. -/
theorem N10_root : ∃ ξs, ((129776751192374311 : ℝ) / 781250000000000000) < ξs ∧ ξs < ((16611424152623911809 : ℝ) / 100000000000000000000) ∧ peval f10 ξs = 0 ∧
    ∀ x, 0 < x → x < 1 → (x < ξs → 0 < peval f10 x) ∧ (ξs < x → peval f10 x < 0) := by
  apply sign_change (peval f10) (0 : ℝ) _ _ (by norm_num) (by norm_num) (by norm_num)
    (continuous_peval _).continuousOn
  · apply peval_strictAntiOn
    intro x hu hv
    rw [N10_dlink]
    linarith [N10_dneg x hu hv]
  · exact (fun x hx0 hx => absurd hx (not_le.mpr hx0))
  · norm_num [f10, peval]
  · norm_num [f10, peval]

theorem N10_A5_neg (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : 0 < peval f10 x) :
    aBlock 10 5 x < 0 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 10 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : (1 - x) ^ 4 * peval [-7, 22, 83, 172, 285, 318, 367, 428, 497, 570, 497, 428, 367, 318, 285, 172, 83, 22, -7] x < 0 := by
    rw [N10_crit_factor]
    have := mul_pos (N10_posH x h0 h1) hf
    nlinarith [pow_pos h1x 4]
  rw [← N10_linkA5 x h0 h1] at hn
  have hd := (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))
  by_contra hc
  push Not at hc
  nlinarith [mul_nonneg hc hd.le]

theorem N10_A5_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : peval f10 x < 0) :
    0 < aBlock 10 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 10 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : 0 < (1 - x) ^ 4 * peval [-7, 22, 83, 172, 285, 318, 367, 428, 497, 570, 497, 428, 367, 318, 285, 172, 83, 22, -7] x := by
    rw [N10_crit_factor]
    have := mul_pos (N10_posH x h0 h1) (neg_pos.mpr hf)
    have h4 := pow_pos h1x 4
    nlinarith
  rw [← N10_linkA5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

/-- Note Thm 6(b) for `N = 10`: there is a unique `ξ* ∈ (0,1)` (the unique root of `f_10` in
`(0,1)`, enclosed in `(129776751192374311/781250000000000000, 16611424152623911809/100000000000000000000)`) such that, for every `ξ ∈ (0,1)`: all `b_m > 0`;
`a_5 < 0` for `ξ < ξ*` (spectrally unstable, since `b_5 > 0`); and all `a_m > 0` for `ξ > ξ*`. -/
theorem cert_N10 : ∃ ξs : ℝ, ((129776751192374311 : ℝ) / 781250000000000000) < ξs ∧ ξs < ((16611424152623911809 : ℝ) / 100000000000000000000) ∧ peval f10 ξs = 0 ∧
    (∀ y, 0 < y → y < 1 → peval f10 y = 0 → y = ξs) ∧
    ∀ x, 0 < x → x < 1 →
      (∀ m, 1 ≤ m → m ≤ 9 → 0 < bBlock 10 m x) ∧
      (x < ξs → aBlock 10 5 x < 0) ∧
      (ξs < x → ∀ m, 1 ≤ m → m ≤ 9 → 0 < aBlock 10 m x) := by
  obtain ⟨ξs, hta, htb, hf0, hs⟩ := N10_root
  refine ⟨ξs, hta, htb, hf0, ?_, ?_⟩
  · intro y hy0 hy1 hy
    rcases lt_trichotomy y ξs with h | h | h
    · linarith [(hs y hy0 hy1).1 h]
    · exact h
    · linarith [(hs y hy0 hy1).2 h]
  intro x h0 h1
  refine ⟨?_, fun hlt => N10_A5_neg x h0 h1 ((hs x h0 h1).1 hlt), ?_⟩
  · intro m hm1 hm2
    interval_cases m
    · exact N10_B1_pos x h0 h1
    · exact N10_B2_pos x h0 h1
    · exact N10_B3_pos x h0 h1
    · exact N10_B4_pos x h0 h1
    · exact N10_B5_pos x h0 h1
    · rw [show bBlock 10 6 x = bBlock 10 4 x from bBlock_symm 10 4 (by norm_num) x]
      exact N10_B4_pos x h0 h1
    · rw [show bBlock 10 7 x = bBlock 10 3 x from bBlock_symm 10 3 (by norm_num) x]
      exact N10_B3_pos x h0 h1
    · rw [show bBlock 10 8 x = bBlock 10 2 x from bBlock_symm 10 2 (by norm_num) x]
      exact N10_B2_pos x h0 h1
    · rw [show bBlock 10 9 x = bBlock 10 1 x from bBlock_symm 10 1 (by norm_num) x]
      exact N10_B1_pos x h0 h1
  · intro hgt m hm1 hm2
    have hlo : ((5 : ℝ) / 32) < x := lt_trans (lt_trans (by norm_num) hta) hgt
    have hfn := (hs x h0 h1).2 hgt
    interval_cases m
    · exact N10_A1_pos x h0 h1
    · exact N10_A2_pos x h0 h1
    · exact N10_A3_pos x h0 hlo h1
    · exact N10_A4_pos x h0 hlo h1
    · exact N10_A5_pos x h0 h1 hfn
    · rw [show aBlock 10 6 x = aBlock 10 4 x from aBlock_symm 10 4 (by norm_num) x]
      exact N10_A4_pos x h0 hlo h1
    · rw [show aBlock 10 7 x = aBlock 10 3 x from aBlock_symm 10 3 (by norm_num) x]
      exact N10_A3_pos x h0 hlo h1
    · rw [show aBlock 10 8 x = aBlock 10 2 x from aBlock_symm 10 2 (by norm_num) x]
      exact N10_A2_pos x h0 h1
    · rw [show aBlock 10 9 x = aBlock 10 1 x from aBlock_symm 10 1 (by norm_num) x]
      exact N10_A1_pos x h0 h1

/-! ### N = 11 -/

theorem N11_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 11 1 x * (2 * (1 - x ^ 11) ^ 2) = (1 - x) ^ 4 * peval [10, 29, 56, 90, 130, 175, 224, 276, 330, 385, 330, 276, 224, 175, 130, 90, 56, 29, 10] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 11 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (11 - 1 : ℕ) = 10 by norm_num, show (2 * 11 - 1 : ℕ) = 21 by norm_num, show (11 + 1 : ℕ) = 12 by norm_num, show (11 - 1 + 11 : ℕ) = 21 by norm_num, show (1 + 11 : ℕ) = 12 by norm_num, show (11 + 11 : ℕ) = 22 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N11_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [10, 29, 56, 90, 130, 175, 224, 276, 330, 385, 330, 276, 224, 175, 130, 90, 56, 29, 10] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 10 * q ^ 18 + 209 * p ^ 1 * q ^ 17 + 2079 * p ^ 2 * q ^ 16 + 13090 * p ^ 3 * q ^ 15 + 58520 * p ^ 4 * q ^ 14 + 197505 * p ^ 5 * q ^ 13 + 522291 * p ^ 6 * q ^ 12 + 1108536 * p ^ 7 * q ^ 11 + 1918620 * p ^ 8 * q ^ 10 + 2735810 * p ^ 9 * q ^ 9 + 3233120 * p ^ 10 * q ^ 8 + 3173346 * p ^ 11 * q ^ 7 + 2581656 * p ^ 12 * q ^ 6 + 1727880 * p ^ 13 * q ^ 5 + 936540 * p ^ 14 * q ^ 4 + 399300 * p ^ 15 * q ^ 3 + 127050 * p ^ 16 * q ^ 2 + 27225 * p ^ 17 * q ^ 1 + 3025 * p ^ 18 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [10, 29, 56, 90, 130, 175, 224, 276, 330, 385, 330, 276, 224, 175, 130, 90, 56, 29, 10] x =
      (fun p q : ℝ => 10 * q ^ 18 + 209 * p ^ 1 * q ^ 17 + 2079 * p ^ 2 * q ^ 16 + 13090 * p ^ 3 * q ^ 15 + 58520 * p ^ 4 * q ^ 14 + 197505 * p ^ 5 * q ^ 13 + 522291 * p ^ 6 * q ^ 12 + 1108536 * p ^ 7 * q ^ 11 + 1918620 * p ^ 8 * q ^ 10 + 2735810 * p ^ 9 * q ^ 9 + 3233120 * p ^ 10 * q ^ 8 + 3173346 * p ^ 11 * q ^ 7 + 2581656 * p ^ 12 * q ^ 6 + 1727880 * p ^ 13 * q ^ 5 + 936540 * p ^ 14 * q ^ 4 + 399300 * p ^ 15 * q ^ 3 + 127050 * p ^ 16 * q ^ 2 + 27225 * p ^ 17 * q ^ 1 + 3025 * p ^ 18) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [10, 29, 56, 90, 130, 175, 224, 276, 330, 385, 330, 276, 224, 175, 130, 90, 56, 29, 10] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N11_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [10, 29, 56, 90, 130, 175, 224, 276, 330, 385, 330, 276, 224, 175, 130, 90, 56, 29, 10] x := N11_posB1_0 x hu hv

theorem N11_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 11 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 11 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N11_posB1 x h0 h1)
  rw [← N11_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N11_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 11 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 11) ^ 2) = (1 - x) ^ 4 * peval [10, 49, 124, 231, 366, 525, 704, 899, 1106, 1321, 1430, 1321, 1106, 899, 704, 525, 366, 231, 124, 49, 10] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 11 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (11 - 1 : ℕ) = 10 by norm_num, show (2 * 11 - 1 : ℕ) = 21 by norm_num, show (11 + 1 : ℕ) = 12 by norm_num, show (11 - 1 + 11 : ℕ) = 21 by norm_num, show (1 + 11 : ℕ) = 12 by norm_num, show (11 + 11 : ℕ) = 22 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N11_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [10, 49, 124, 231, 366, 525, 704, 899, 1106, 1321, 1430, 1321, 1106, 899, 704, 525, 366, 231, 124, 49, 10] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 10 * q ^ 20 + 249 * p ^ 1 * q ^ 19 + 2955 * p ^ 2 * q ^ 18 + 22242 * p ^ 3 * q ^ 17 + 119196 * p ^ 4 * q ^ 16 + 483945 * p ^ 5 * q ^ 15 + 1546391 * p ^ 6 * q ^ 14 + 3987720 * p ^ 7 * q ^ 13 + 8441928 * p ^ 8 * q ^ 12 + 14844434 * p ^ 9 * q ^ 11 + 21850840 * p ^ 10 * q ^ 10 + 27049066 * p ^ 11 * q ^ 9 + 28207520 * p ^ 12 * q ^ 8 + 24747888 * p ^ 13 * q ^ 7 + 18174684 * p ^ 14 * q ^ 6 + 11056980 * p ^ 15 * q ^ 5 + 5470410 * p ^ 16 * q ^ 4 + 2132625 * p ^ 17 * q ^ 3 + 620125 * p ^ 18 * q ^ 2 + 121000 * p ^ 19 * q ^ 1 + 12100 * p ^ 20 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [10, 49, 124, 231, 366, 525, 704, 899, 1106, 1321, 1430, 1321, 1106, 899, 704, 525, 366, 231, 124, 49, 10] x =
      (fun p q : ℝ => 10 * q ^ 20 + 249 * p ^ 1 * q ^ 19 + 2955 * p ^ 2 * q ^ 18 + 22242 * p ^ 3 * q ^ 17 + 119196 * p ^ 4 * q ^ 16 + 483945 * p ^ 5 * q ^ 15 + 1546391 * p ^ 6 * q ^ 14 + 3987720 * p ^ 7 * q ^ 13 + 8441928 * p ^ 8 * q ^ 12 + 14844434 * p ^ 9 * q ^ 11 + 21850840 * p ^ 10 * q ^ 10 + 27049066 * p ^ 11 * q ^ 9 + 28207520 * p ^ 12 * q ^ 8 + 24747888 * p ^ 13 * q ^ 7 + 18174684 * p ^ 14 * q ^ 6 + 11056980 * p ^ 15 * q ^ 5 + 5470410 * p ^ 16 * q ^ 4 + 2132625 * p ^ 17 * q ^ 3 + 620125 * p ^ 18 * q ^ 2 + 121000 * p ^ 19 * q ^ 1 + 12100 * p ^ 20) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [10, 49, 124, 231, 366, 525, 704, 899, 1106, 1321, 1430, 1321, 1106, 899, 704, 525, 366, 231, 124, 49, 10] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N11_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [10, 49, 124, 231, 366, 525, 704, 899, 1106, 1321, 1430, 1321, 1106, 899, 704, 525, 366, 231, 124, 49, 10] x := N11_posA1_0 x hu hv

theorem N11_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 11 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 11 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N11_posA1 x h0 h1)
  rw [← N11_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N11_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 11 2 x * (2 * (1 - x ^ 11) ^ 2) = (1 - x) ^ 4 * peval [18, 72, 158, 272, 410, 568, 742, 928, 1122, 1221, 1122, 928, 742, 568, 410, 272, 158, 72, 18] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 11 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (11 - 2 : ℕ) = 9 by norm_num, show (2 * 11 - 2 : ℕ) = 20 by norm_num, show (11 + 2 : ℕ) = 13 by norm_num, show (11 - 2 + 11 : ℕ) = 20 by norm_num, show (2 + 11 : ℕ) = 13 by norm_num, show (11 + 11 : ℕ) = 22 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N11_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [18, 72, 158, 272, 410, 568, 742, 928, 1122, 1221, 1122, 928, 742, 568, 410, 272, 158, 72, 18] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 18 * q ^ 18 + 396 * p ^ 1 * q ^ 17 + 4136 * p ^ 2 * q ^ 16 + 27280 * p ^ 3 * q ^ 15 + 127490 * p ^ 4 * q ^ 14 + 448932 * p ^ 5 * q ^ 13 + 1236444 * p ^ 6 * q ^ 12 + 2728704 * p ^ 7 * q ^ 11 + 4903140 * p ^ 8 * q ^ 10 + 7248021 * p ^ 9 * q ^ 9 + 8867001 * p ^ 10 * q ^ 8 + 8995140 * p ^ 11 * q ^ 7 + 7548948 * p ^ 12 * q ^ 6 + 5198886 * p ^ 13 * q ^ 5 + 2890206 * p ^ 14 * q ^ 4 + 1258884 * p ^ 15 * q ^ 3 + 407286 * p ^ 16 * q ^ 2 + 88209 * p ^ 17 * q ^ 1 + 9801 * p ^ 18 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [18, 72, 158, 272, 410, 568, 742, 928, 1122, 1221, 1122, 928, 742, 568, 410, 272, 158, 72, 18] x =
      (fun p q : ℝ => 18 * q ^ 18 + 396 * p ^ 1 * q ^ 17 + 4136 * p ^ 2 * q ^ 16 + 27280 * p ^ 3 * q ^ 15 + 127490 * p ^ 4 * q ^ 14 + 448932 * p ^ 5 * q ^ 13 + 1236444 * p ^ 6 * q ^ 12 + 2728704 * p ^ 7 * q ^ 11 + 4903140 * p ^ 8 * q ^ 10 + 7248021 * p ^ 9 * q ^ 9 + 8867001 * p ^ 10 * q ^ 8 + 8995140 * p ^ 11 * q ^ 7 + 7548948 * p ^ 12 * q ^ 6 + 5198886 * p ^ 13 * q ^ 5 + 2890206 * p ^ 14 * q ^ 4 + 1258884 * p ^ 15 * q ^ 3 + 407286 * p ^ 16 * q ^ 2 + 88209 * p ^ 17 * q ^ 1 + 9801 * p ^ 18) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [18, 72, 158, 272, 410, 568, 742, 928, 1122, 1221, 1122, 928, 742, 568, 410, 272, 158, 72, 18] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N11_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [18, 72, 158, 272, 410, 568, 742, 928, 1122, 1221, 1122, 928, 742, 568, 410, 272, 158, 72, 18] x := N11_posB2_0 x hu hv

theorem N11_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 11 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 11 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N11_posB2 x h0 h1)
  rw [← N11_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N11_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 11 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 11) ^ 2) = (1 - x) ^ 4 * peval [2, 44, 100, 188, 304, 444, 604, 780, 968, 1065, 1166, 1065, 968, 780, 604, 444, 304, 188, 100, 44, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 11 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (11 - 2 : ℕ) = 9 by norm_num, show (2 * 11 - 2 : ℕ) = 20 by norm_num, show (11 + 2 : ℕ) = 13 by norm_num, show (11 - 2 + 11 : ℕ) = 20 by norm_num, show (2 + 11 : ℕ) = 13 by norm_num, show (11 + 11 : ℕ) = 22 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N11_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 44, 100, 188, 304, 444, 604, 780, 968, 1065, 1166, 1065, 968, 780, 604, 444, 304, 188, 100, 44, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 20 + 84 * p ^ 1 * q ^ 19 + 1316 * p ^ 2 * q ^ 18 + 11792 * p ^ 3 * q ^ 17 + 71126 * p ^ 4 * q ^ 16 + 314028 * p ^ 5 * q ^ 15 + 1066736 * p ^ 6 * q ^ 14 + 2879184 * p ^ 7 * q ^ 13 + 6310128 * p ^ 8 * q ^ 12 + 11397925 * p ^ 9 * q ^ 11 + 17138385 * p ^ 10 * q ^ 10 + 21583980 * p ^ 11 * q ^ 9 + 22828828 * p ^ 12 * q ^ 8 + 20261934 * p ^ 13 * q ^ 7 + 15017310 * p ^ 14 * q ^ 6 + 9198420 * p ^ 15 * q ^ 5 + 4571622 * p ^ 16 * q ^ 4 + 1787049 * p ^ 17 * q ^ 3 + 520421 * p ^ 18 * q ^ 2 + 101640 * p ^ 19 * q ^ 1 + 10164 * p ^ 20 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 44, 100, 188, 304, 444, 604, 780, 968, 1065, 1166, 1065, 968, 780, 604, 444, 304, 188, 100, 44, 2] x =
      (fun p q : ℝ => 2 * q ^ 20 + 84 * p ^ 1 * q ^ 19 + 1316 * p ^ 2 * q ^ 18 + 11792 * p ^ 3 * q ^ 17 + 71126 * p ^ 4 * q ^ 16 + 314028 * p ^ 5 * q ^ 15 + 1066736 * p ^ 6 * q ^ 14 + 2879184 * p ^ 7 * q ^ 13 + 6310128 * p ^ 8 * q ^ 12 + 11397925 * p ^ 9 * q ^ 11 + 17138385 * p ^ 10 * q ^ 10 + 21583980 * p ^ 11 * q ^ 9 + 22828828 * p ^ 12 * q ^ 8 + 20261934 * p ^ 13 * q ^ 7 + 15017310 * p ^ 14 * q ^ 6 + 9198420 * p ^ 15 * q ^ 5 + 4571622 * p ^ 16 * q ^ 4 + 1787049 * p ^ 17 * q ^ 3 + 520421 * p ^ 18 * q ^ 2 + 101640 * p ^ 19 * q ^ 1 + 10164 * p ^ 20) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 44, 100, 188, 304, 444, 604, 780, 968, 1065, 1166, 1065, 968, 780, 604, 444, 304, 188, 100, 44, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N11_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 44, 100, 188, 304, 444, 604, 780, 968, 1065, 1166, 1065, 968, 780, 604, 444, 304, 188, 100, 44, 2] x := N11_posA2_0 x hu hv

theorem N11_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 11 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 11 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N11_posA2 x h0 h1)
  rw [← N11_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N11_linkB3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 11 3 x * (2 * (1 - x ^ 11) ^ 2) = (1 - x) ^ 4 * peval [24, 96, 240, 447, 708, 1014, 1356, 1725, 2024, 2156, 2024, 1725, 1356, 1014, 708, 447, 240, 96, 24] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 11 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (11 - 3 : ℕ) = 8 by norm_num, show (2 * 11 - 3 : ℕ) = 19 by norm_num, show (11 + 3 : ℕ) = 14 by norm_num, show (11 - 3 + 11 : ℕ) = 19 by norm_num, show (3 + 11 : ℕ) = 14 by norm_num, show (11 + 11 : ℕ) = 22 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N11_posB3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [24, 96, 240, 447, 708, 1014, 1356, 1725, 2024, 2156, 2024, 1725, 1356, 1014, 708, 447, 240, 96, 24] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 24 * q ^ 18 + 528 * p ^ 1 * q ^ 17 + 5544 * p ^ 2 * q ^ 16 + 36927 * p ^ 3 * q ^ 15 + 174933 * p ^ 4 * q ^ 14 + 626373 * p ^ 5 * q ^ 13 + 1758735 * p ^ 6 * q ^ 12 + 3965148 * p ^ 7 * q ^ 11 + 7290668 * p ^ 8 * q ^ 10 + 11041492 * p ^ 9 * q ^ 9 + 13848692 * p ^ 10 * q ^ 8 + 14404324 * p ^ 11 * q ^ 7 + 12385076 * p ^ 12 * q ^ 6 + 8724100 * p ^ 13 * q ^ 5 + 4947932 * p ^ 14 * q ^ 4 + 2191552 * p ^ 15 * q ^ 3 + 718256 * p ^ 16 * q ^ 2 + 156816 * p ^ 17 * q ^ 1 + 17424 * p ^ 18 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [24, 96, 240, 447, 708, 1014, 1356, 1725, 2024, 2156, 2024, 1725, 1356, 1014, 708, 447, 240, 96, 24] x =
      (fun p q : ℝ => 24 * q ^ 18 + 528 * p ^ 1 * q ^ 17 + 5544 * p ^ 2 * q ^ 16 + 36927 * p ^ 3 * q ^ 15 + 174933 * p ^ 4 * q ^ 14 + 626373 * p ^ 5 * q ^ 13 + 1758735 * p ^ 6 * q ^ 12 + 3965148 * p ^ 7 * q ^ 11 + 7290668 * p ^ 8 * q ^ 10 + 11041492 * p ^ 9 * q ^ 9 + 13848692 * p ^ 10 * q ^ 8 + 14404324 * p ^ 11 * q ^ 7 + 12385076 * p ^ 12 * q ^ 6 + 8724100 * p ^ 13 * q ^ 5 + 4947932 * p ^ 14 * q ^ 4 + 2191552 * p ^ 15 * q ^ 3 + 718256 * p ^ 16 * q ^ 2 + 156816 * p ^ 17 * q ^ 1 + 17424 * p ^ 18) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [24, 96, 240, 447, 708, 1014, 1356, 1725, 2024, 2156, 2024, 1725, 1356, 1014, 708, 447, 240, 96, 24] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N11_posB3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [24, 96, 240, 447, 708, 1014, 1356, 1725, 2024, 2156, 2024, 1725, 1356, 1014, 708, 447, 240, 96, 24] x := N11_posB3_0 x hu hv

theorem N11_B3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 11 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 11 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N11_posB3 x h0 h1)
  rw [← N11_linkB3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N11_linkA3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 11 3 x * (2 * (1 - x) ^ 2 * (1 - x ^ 11) ^ 2) = (1 - x) ^ 4 * peval [-4, 32, 104, 175, 274, 397, 540, 699, 782, 873, 968, 873, 782, 699, 540, 397, 274, 175, 104, 32, -4] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 11 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (11 - 3 : ℕ) = 8 by norm_num, show (2 * 11 - 3 : ℕ) = 19 by norm_num, show (11 + 3 : ℕ) = 14 by norm_num, show (11 - 3 + 11 : ℕ) = 19 by norm_num, show (3 + 11 : ℕ) = 14 by norm_num, show (11 + 11 : ℕ) = 22 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N11_posA3_0 (x : ℝ) (hu : ((3 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-4, 32, 104, 175, 274, 397, 540, 699, 782, 873, 968, 873, 782, 699, 540, 397, 274, 175, 104, 32, -4] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 8790264015839923102007228 * q ^ 20 + 274192433233135202847526656 * p ^ 1 * q ^ 19 + 3779814361792232413435379712 * p ^ 2 * q ^ 18 + 31583102219292470569863467008 * p ^ 3 * q ^ 17 + 182335595632098140062389436416 * p ^ 4 * q ^ 16 + 780777821604891538667206606848 * p ^ 5 * q ^ 15 + 2590350485545945966580766605312 * p ^ 6 * q ^ 14 + 6852086653586649348445424320512 * p ^ 7 * q ^ 13 + 14736086282723605517169984536576 * p ^ 8 * q ^ 12 + 26107960120682512913414714556416 * p ^ 9 * q ^ 11 + 38431569122686030200072889696256 * p ^ 10 * q ^ 10 + 47220073734731913242884792385536 * p ^ 11 * q ^ 9 + 48475382585219540787790570586112 * p ^ 12 * q ^ 8 + 41464786078640860348137506603008 * p ^ 13 * q ^ 7 + 29344001966625874305737623076864 * p ^ 14 * q ^ 6 + 16963385136332489564994578415616 * p ^ 15 * q ^ 5 + 7845397296945301476284970303488 * p ^ 16 * q ^ 4 + 2806342836728340468057407225856 * p ^ 17 * q ^ 3 + 732894919863872276714598957056 * p ^ 18 * q ^ 2 + 125069420668231461269227438080 * p ^ 19 * q ^ 1 + 10532161740482649370040205312 * p ^ 20 := by
    intro p q hp hq; positivity
  have key : (19004963774880799438801 : ℝ) * peval [-4, 32, 104, 175, 274, 397, 540, 699, 782, 873, 968, 873, 782, 699, 540, 397, 274, 175, 104, 32, -4] x =
      (fun p q : ℝ => 8790264015839923102007228 * q ^ 20 + 274192433233135202847526656 * p ^ 1 * q ^ 19 + 3779814361792232413435379712 * p ^ 2 * q ^ 18 + 31583102219292470569863467008 * p ^ 3 * q ^ 17 + 182335595632098140062389436416 * p ^ 4 * q ^ 16 + 780777821604891538667206606848 * p ^ 5 * q ^ 15 + 2590350485545945966580766605312 * p ^ 6 * q ^ 14 + 6852086653586649348445424320512 * p ^ 7 * q ^ 13 + 14736086282723605517169984536576 * p ^ 8 * q ^ 12 + 26107960120682512913414714556416 * p ^ 9 * q ^ 11 + 38431569122686030200072889696256 * p ^ 10 * q ^ 10 + 47220073734731913242884792385536 * p ^ 11 * q ^ 9 + 48475382585219540787790570586112 * p ^ 12 * q ^ 8 + 41464786078640860348137506603008 * p ^ 13 * q ^ 7 + 29344001966625874305737623076864 * p ^ 14 * q ^ 6 + 16963385136332489564994578415616 * p ^ 15 * q ^ 5 + 7845397296945301476284970303488 * p ^ 16 * q ^ 4 + 2806342836728340468057407225856 * p ^ 17 * q ^ 3 + 732894919863872276714598957056 * p ^ 18 * q ^ 2 + 125069420668231461269227438080 * p ^ 19 * q ^ 1 + 10532161740482649370040205312 * p ^ 20) (x - ((3 : ℝ) / 16)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (19004963774880799438801 : ℝ) * peval [-4, 32, 104, 175, 274, 397, 540, 699, 782, 873, 968, 873, 782, 699, 540, 397, 274, 175, 104, 32, -4] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N11_posA3 (x : ℝ) (hu : ((3 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-4, 32, 104, 175, 274, 397, 540, 699, 782, 873, 968, 873, 782, 699, 540, 397, 274, 175, 104, 32, -4] x := N11_posA3_0 x hu hv

theorem N11_A3_pos (x : ℝ) (h0 : 0 < x) (hlo : ((3 : ℝ) / 16) < x) (h1 : x < 1) :
    0 < aBlock 11 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 11 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N11_posA3 x hlo h1)
  rw [← N11_linkA3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N11_linkB4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 11 4 x * (2 * (1 - x ^ 11) ^ 2) = (1 - x) ^ 4 * peval [28, 112, 280, 560, 936, 1392, 1912, 2403, 2772, 2926, 2772, 2403, 1912, 1392, 936, 560, 280, 112, 28] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 11 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (11 - 4 : ℕ) = 7 by norm_num, show (2 * 11 - 4 : ℕ) = 18 by norm_num, show (11 + 4 : ℕ) = 15 by norm_num, show (11 - 4 + 11 : ℕ) = 18 by norm_num, show (4 + 11 : ℕ) = 15 by norm_num, show (11 + 11 : ℕ) = 22 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N11_posB4_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [28, 112, 280, 560, 936, 1392, 1912, 2403, 2772, 2926, 2772, 2403, 1912, 1392, 936, 560, 280, 112, 28] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 28 * q ^ 18 + 616 * p ^ 1 * q ^ 17 + 6468 * p ^ 2 * q ^ 16 + 43120 * p ^ 3 * q ^ 15 + 204776 * p ^ 4 * q ^ 14 + 736560 * p ^ 5 * q ^ 13 + 2082432 * p ^ 6 * q ^ 12 + 4739251 * p ^ 7 * q ^ 11 + 8817765 * p ^ 8 * q ^ 10 + 13542683 * p ^ 9 * q ^ 9 + 17254237 * p ^ 10 * q ^ 8 + 18246921 * p ^ 11 * q ^ 7 + 15951551 * p ^ 12 * q ^ 6 + 11413325 * p ^ 13 * q ^ 5 + 6563403 * p ^ 14 * q ^ 4 + 2940784 * p ^ 15 * q ^ 3 + 972356 * p ^ 16 * q ^ 2 + 213444 * p ^ 17 * q ^ 1 + 23716 * p ^ 18 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [28, 112, 280, 560, 936, 1392, 1912, 2403, 2772, 2926, 2772, 2403, 1912, 1392, 936, 560, 280, 112, 28] x =
      (fun p q : ℝ => 28 * q ^ 18 + 616 * p ^ 1 * q ^ 17 + 6468 * p ^ 2 * q ^ 16 + 43120 * p ^ 3 * q ^ 15 + 204776 * p ^ 4 * q ^ 14 + 736560 * p ^ 5 * q ^ 13 + 2082432 * p ^ 6 * q ^ 12 + 4739251 * p ^ 7 * q ^ 11 + 8817765 * p ^ 8 * q ^ 10 + 13542683 * p ^ 9 * q ^ 9 + 17254237 * p ^ 10 * q ^ 8 + 18246921 * p ^ 11 * q ^ 7 + 15951551 * p ^ 12 * q ^ 6 + 11413325 * p ^ 13 * q ^ 5 + 6563403 * p ^ 14 * q ^ 4 + 2940784 * p ^ 15 * q ^ 3 + 972356 * p ^ 16 * q ^ 2 + 213444 * p ^ 17 * q ^ 1 + 23716 * p ^ 18) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [28, 112, 280, 560, 936, 1392, 1912, 2403, 2772, 2926, 2772, 2403, 1912, 1392, 936, 560, 280, 112, 28] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N11_posB4 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [28, 112, 280, 560, 936, 1392, 1912, 2403, 2772, 2926, 2772, 2403, 1912, 1392, 936, 560, 280, 112, 28] x := N11_posB4_0 x hu hv

theorem N11_B4_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 11 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 11 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N11_posB4 x h0 h1)
  rw [← N11_linkB4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N11_linkA4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 11 4 x * (2 * (1 - x) ^ 2 * (1 - x ^ 11) ^ 2) = (1 - x) ^ 4 * peval [-8, 24, 92, 192, 276, 384, 512, 579, 658, 745, 836, 745, 658, 579, 512, 384, 276, 192, 92, 24, -8] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 11 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (11 - 4 : ℕ) = 7 by norm_num, show (2 * 11 - 4 : ℕ) = 18 by norm_num, show (11 + 4 : ℕ) = 15 by norm_num, show (11 - 4 + 11 : ℕ) = 18 by norm_num, show (4 + 11 : ℕ) = 15 by norm_num, show (11 + 11 : ℕ) = 22 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N11_posA4_0 (x : ℝ) (hu : ((3 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-8, 24, 92, 192, 276, 384, 512, 579, 658, 745, 836, 745, 658, 579, 512, 384, 276, 192, 92, 24, -8] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 1763039060286064354338808 * q ^ 20 + 123019591452955288502549120 * p ^ 1 * q ^ 19 + 2239370574350346067311940608 * p ^ 2 * q ^ 18 + 21696449004217316090703642624 * p ^ 3 * q ^ 17 + 137474578112821798866713051136 * p ^ 4 * q ^ 16 + 627490135038641744204884082688 * p ^ 5 * q ^ 15 + 2179577855985847209892939563008 * p ^ 6 * q ^ 14 + 5962704737888917319579783397376 * p ^ 7 * q ^ 13 + 13141511332528963618302628724736 * p ^ 8 * q ^ 12 + 23687998350142699835314763988992 * p ^ 9 * q ^ 11 + 35262456720158808159710569562112 * p ^ 10 * q ^ 10 + 43589268250596112897045663055872 * p ^ 11 * q ^ 9 + 44821233400219792214918491537408 * p ^ 12 * q ^ 8 + 38260128199862431981422654783488 * p ^ 13 * q ^ 7 + 26941758546383580885108746354688 * p ^ 14 * q ^ 6 + 15466221193357912801972043382784 * p ^ 15 * q ^ 5 + 7096288534565599866250915217408 * p ^ 16 * q ^ 4 + 2518837877459161778867618709504 * p ^ 17 * q ^ 3 + 653736856157679656134607110144 * p ^ 18 * q ^ 2 + 111172818371761298905979944960 * p ^ 19 * q ^ 1 + 9361921547095688328924626944 * p ^ 20 := by
    intro p q hp hq; positivity
  have key : (19004963774880799438801 : ℝ) * peval [-8, 24, 92, 192, 276, 384, 512, 579, 658, 745, 836, 745, 658, 579, 512, 384, 276, 192, 92, 24, -8] x =
      (fun p q : ℝ => 1763039060286064354338808 * q ^ 20 + 123019591452955288502549120 * p ^ 1 * q ^ 19 + 2239370574350346067311940608 * p ^ 2 * q ^ 18 + 21696449004217316090703642624 * p ^ 3 * q ^ 17 + 137474578112821798866713051136 * p ^ 4 * q ^ 16 + 627490135038641744204884082688 * p ^ 5 * q ^ 15 + 2179577855985847209892939563008 * p ^ 6 * q ^ 14 + 5962704737888917319579783397376 * p ^ 7 * q ^ 13 + 13141511332528963618302628724736 * p ^ 8 * q ^ 12 + 23687998350142699835314763988992 * p ^ 9 * q ^ 11 + 35262456720158808159710569562112 * p ^ 10 * q ^ 10 + 43589268250596112897045663055872 * p ^ 11 * q ^ 9 + 44821233400219792214918491537408 * p ^ 12 * q ^ 8 + 38260128199862431981422654783488 * p ^ 13 * q ^ 7 + 26941758546383580885108746354688 * p ^ 14 * q ^ 6 + 15466221193357912801972043382784 * p ^ 15 * q ^ 5 + 7096288534565599866250915217408 * p ^ 16 * q ^ 4 + 2518837877459161778867618709504 * p ^ 17 * q ^ 3 + 653736856157679656134607110144 * p ^ 18 * q ^ 2 + 111172818371761298905979944960 * p ^ 19 * q ^ 1 + 9361921547095688328924626944 * p ^ 20) (x - ((3 : ℝ) / 16)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (19004963774880799438801 : ℝ) * peval [-8, 24, 92, 192, 276, 384, 512, 579, 658, 745, 836, 745, 658, 579, 512, 384, 276, 192, 92, 24, -8] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N11_posA4 (x : ℝ) (hu : ((3 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-8, 24, 92, 192, 276, 384, 512, 579, 658, 745, 836, 745, 658, 579, 512, 384, 276, 192, 92, 24, -8] x := N11_posA4_0 x hu hv

theorem N11_A4_pos (x : ℝ) (h0 : 0 < x) (hlo : ((3 : ℝ) / 16) < x) (h1 : x < 1) :
    0 < aBlock 11 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 11 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N11_posA4 x hlo h1)
  rw [← N11_linkA4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N11_linkB5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 11 5 x * (2 * (1 - x ^ 11) ^ 2) = (1 - x) ^ 4 * peval [30, 120, 300, 600, 1050, 1625, 2234, 2786, 3190, 3355, 3190, 2786, 2234, 1625, 1050, 600, 300, 120, 30] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 11 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (11 - 5 : ℕ) = 6 by norm_num, show (2 * 11 - 5 : ℕ) = 17 by norm_num, show (11 + 5 : ℕ) = 16 by norm_num, show (11 - 5 + 11 : ℕ) = 17 by norm_num, show (5 + 11 : ℕ) = 16 by norm_num, show (11 + 11 : ℕ) = 22 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N11_posB5_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [30, 120, 300, 600, 1050, 1625, 2234, 2786, 3190, 3355, 3190, 2786, 2234, 1625, 1050, 600, 300, 120, 30] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 30 * q ^ 18 + 660 * p ^ 1 * q ^ 17 + 6930 * p ^ 2 * q ^ 16 + 46200 * p ^ 3 * q ^ 15 + 219450 * p ^ 4 * q ^ 14 + 789965 * p ^ 5 * q ^ 13 + 2237389 * p ^ 6 * q ^ 12 + 5107784 * p ^ 7 * q ^ 11 + 9547780 * p ^ 8 * q ^ 10 + 14754740 * p ^ 9 * q ^ 9 + 18938920 * p ^ 10 * q ^ 8 + 20194174 * p ^ 11 * q ^ 7 + 17803214 * p ^ 12 * q ^ 6 + 12840520 * p ^ 13 * q ^ 5 + 7436660 * p ^ 14 * q ^ 4 + 3351700 * p ^ 15 * q ^ 3 + 1113200 * p ^ 16 * q ^ 2 + 245025 * p ^ 17 * q ^ 1 + 27225 * p ^ 18 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [30, 120, 300, 600, 1050, 1625, 2234, 2786, 3190, 3355, 3190, 2786, 2234, 1625, 1050, 600, 300, 120, 30] x =
      (fun p q : ℝ => 30 * q ^ 18 + 660 * p ^ 1 * q ^ 17 + 6930 * p ^ 2 * q ^ 16 + 46200 * p ^ 3 * q ^ 15 + 219450 * p ^ 4 * q ^ 14 + 789965 * p ^ 5 * q ^ 13 + 2237389 * p ^ 6 * q ^ 12 + 5107784 * p ^ 7 * q ^ 11 + 9547780 * p ^ 8 * q ^ 10 + 14754740 * p ^ 9 * q ^ 9 + 18938920 * p ^ 10 * q ^ 8 + 20194174 * p ^ 11 * q ^ 7 + 17803214 * p ^ 12 * q ^ 6 + 12840520 * p ^ 13 * q ^ 5 + 7436660 * p ^ 14 * q ^ 4 + 3351700 * p ^ 15 * q ^ 3 + 1113200 * p ^ 16 * q ^ 2 + 245025 * p ^ 17 * q ^ 1 + 27225 * p ^ 18) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [30, 120, 300, 600, 1050, 1625, 2234, 2786, 3190, 3355, 3190, 2786, 2234, 1625, 1050, 600, 300, 120, 30] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N11_posB5 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [30, 120, 300, 600, 1050, 1625, 2234, 2786, 3190, 3355, 3190, 2786, 2234, 1625, 1050, 600, 300, 120, 30] x := N11_posB5_0 x hu hv

theorem N11_B5_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 11 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 11 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N11_posB5 x h0 h1)
  rw [← N11_linkB5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N11_linkA5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 11 5 x * (2 * (1 - x) ^ 2 * (1 - x ^ 11) ^ 2) = (1 - x) ^ 4 * peval [-10, 20, 86, 184, 310, 405, 454, 519, 596, 681, 770, 681, 596, 519, 454, 405, 310, 184, 86, 20, -10] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 11 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (11 - 5 : ℕ) = 6 by norm_num, show (2 * 11 - 5 : ℕ) = 17 by norm_num, show (11 + 5 : ℕ) = 16 by norm_num, show (11 - 5 + 11 : ℕ) = 17 by norm_num, show (5 + 11 : ℕ) = 16 by norm_num, show (11 + 11 : ℕ) = 22 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

/-- The minimal polynomial `f_11` of `ξ*(11)` (constant term first). -/
def f11 : List ℤ := [10, -20, -86, -184, -310, -405, -454, -519, -596, -681, -770, -681, -596, -519, -454, -405, -310, -184, -86, -20, 10]

theorem N11_crit_factor (x : ℝ) : peval [-10, 20, 86, 184, 310, 405, 454, 519, 596, 681, 770, 681, 596, 519, 454, 405, 310, 184, 86, 20, -10] x = -(peval [1] x * peval f11 x) := by
  simp only [f11, peval]; ring

theorem N11_posH (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1] x := by norm_num [peval]

theorem N11_dneg_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [20, 172, 552, 1240, 2025, 2724, 3633, 4768, 6129, 7700, 7491, 7152, 6747, 6356, 6075, 4960, 3128, 1548, 380, -200] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 20 * q ^ 19 + 552 * p ^ 1 * q ^ 18 + 7068 * p ^ 2 * q ^ 17 + 56320 * p ^ 3 * q ^ 16 + 314809 * p ^ 4 * q ^ 15 + 1316139 * p ^ 5 * q ^ 14 + 4278890 * p ^ 6 * q ^ 13 + 11094600 * p ^ 7 * q ^ 12 + 23335620 * p ^ 8 * q ^ 11 + 40275796 * p ^ 9 * q ^ 10 + 57471634 * p ^ 10 * q ^ 9 + 68095170 * p ^ 11 * q ^ 8 + 67072720 * p ^ 12 * q ^ 7 + 54783960 * p ^ 13 * q ^ 6 + 36827076 * p ^ 14 * q ^ 5 + 20073900 * p ^ 15 * q ^ 4 + 8642425 * p ^ 16 * q ^ 3 + 2804175 * p ^ 17 * q ^ 2 + 623150 * p ^ 18 * q ^ 1 + 72600 * p ^ 19 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [20, 172, 552, 1240, 2025, 2724, 3633, 4768, 6129, 7700, 7491, 7152, 6747, 6356, 6075, 4960, 3128, 1548, 380, -200] x =
      (fun p q : ℝ => 20 * q ^ 19 + 552 * p ^ 1 * q ^ 18 + 7068 * p ^ 2 * q ^ 17 + 56320 * p ^ 3 * q ^ 16 + 314809 * p ^ 4 * q ^ 15 + 1316139 * p ^ 5 * q ^ 14 + 4278890 * p ^ 6 * q ^ 13 + 11094600 * p ^ 7 * q ^ 12 + 23335620 * p ^ 8 * q ^ 11 + 40275796 * p ^ 9 * q ^ 10 + 57471634 * p ^ 10 * q ^ 9 + 68095170 * p ^ 11 * q ^ 8 + 67072720 * p ^ 12 * q ^ 7 + 54783960 * p ^ 13 * q ^ 6 + 36827076 * p ^ 14 * q ^ 5 + 20073900 * p ^ 15 * q ^ 4 + 8642425 * p ^ 16 * q ^ 3 + 2804175 * p ^ 17 * q ^ 2 + 623150 * p ^ 18 * q ^ 1 + 72600 * p ^ 19) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [20, 172, 552, 1240, 2025, 2724, 3633, 4768, 6129, 7700, 7491, 7152, 6747, 6356, 6075, 4960, 3128, 1548, 380, -200] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N11_dneg (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [20, 172, 552, 1240, 2025, 2724, 3633, 4768, 6129, 7700, 7491, 7152, 6747, 6356, 6075, 4960, 3128, 1548, 380, -200] x := N11_dneg_0 x hu hv

theorem N11_dlink (x : ℝ) : dpeval f11 x = -peval [20, 172, 552, 1240, 2025, 2724, 3633, 4768, 6129, 7700, 7491, 7152, 6747, 6356, 6075, 4960, 3128, 1548, 380, -200] x := by
  simp only [f11, dpeval, peval]; ring

/-- `f_11` changes sign exactly once in `(0,1)`, at a point of `(20460117688266929063/100000000000000000000, 2557514711033366133/12500000000000000000)`. -/
theorem N11_root : ∃ ξs, ((20460117688266929063 : ℝ) / 100000000000000000000) < ξs ∧ ξs < ((2557514711033366133 : ℝ) / 12500000000000000000) ∧ peval f11 ξs = 0 ∧
    ∀ x, 0 < x → x < 1 → (x < ξs → 0 < peval f11 x) ∧ (ξs < x → peval f11 x < 0) := by
  apply sign_change (peval f11) (0 : ℝ) _ _ (by norm_num) (by norm_num) (by norm_num)
    (continuous_peval _).continuousOn
  · apply peval_strictAntiOn
    intro x hu hv
    rw [N11_dlink]
    linarith [N11_dneg x hu hv]
  · exact (fun x hx0 hx => absurd hx (not_le.mpr hx0))
  · norm_num [f11, peval]
  · norm_num [f11, peval]

theorem N11_A5_neg (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : 0 < peval f11 x) :
    aBlock 11 5 x < 0 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 11 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : (1 - x) ^ 4 * peval [-10, 20, 86, 184, 310, 405, 454, 519, 596, 681, 770, 681, 596, 519, 454, 405, 310, 184, 86, 20, -10] x < 0 := by
    rw [N11_crit_factor]
    have := mul_pos (N11_posH x h0 h1) hf
    nlinarith [pow_pos h1x 4]
  rw [← N11_linkA5 x h0 h1] at hn
  have hd := (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))
  by_contra hc
  push Not at hc
  nlinarith [mul_nonneg hc hd.le]

theorem N11_A5_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : peval f11 x < 0) :
    0 < aBlock 11 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 11 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : 0 < (1 - x) ^ 4 * peval [-10, 20, 86, 184, 310, 405, 454, 519, 596, 681, 770, 681, 596, 519, 454, 405, 310, 184, 86, 20, -10] x := by
    rw [N11_crit_factor]
    have := mul_pos (N11_posH x h0 h1) (neg_pos.mpr hf)
    have h4 := pow_pos h1x 4
    nlinarith
  rw [← N11_linkA5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

/-- Note Thm 6(b) for `N = 11`: there is a unique `ξ* ∈ (0,1)` (the unique root of `f_11` in
`(0,1)`, enclosed in `(20460117688266929063/100000000000000000000, 2557514711033366133/12500000000000000000)`) such that, for every `ξ ∈ (0,1)`: all `b_m > 0`;
`a_5 < 0` for `ξ < ξ*` (spectrally unstable, since `b_5 > 0`); and all `a_m > 0` for `ξ > ξ*`. -/
theorem cert_N11 : ∃ ξs : ℝ, ((20460117688266929063 : ℝ) / 100000000000000000000) < ξs ∧ ξs < ((2557514711033366133 : ℝ) / 12500000000000000000) ∧ peval f11 ξs = 0 ∧
    (∀ y, 0 < y → y < 1 → peval f11 y = 0 → y = ξs) ∧
    ∀ x, 0 < x → x < 1 →
      (∀ m, 1 ≤ m → m ≤ 10 → 0 < bBlock 11 m x) ∧
      (x < ξs → aBlock 11 5 x < 0) ∧
      (ξs < x → ∀ m, 1 ≤ m → m ≤ 10 → 0 < aBlock 11 m x) := by
  obtain ⟨ξs, hta, htb, hf0, hs⟩ := N11_root
  refine ⟨ξs, hta, htb, hf0, ?_, ?_⟩
  · intro y hy0 hy1 hy
    rcases lt_trichotomy y ξs with h | h | h
    · linarith [(hs y hy0 hy1).1 h]
    · exact h
    · linarith [(hs y hy0 hy1).2 h]
  intro x h0 h1
  refine ⟨?_, fun hlt => N11_A5_neg x h0 h1 ((hs x h0 h1).1 hlt), ?_⟩
  · intro m hm1 hm2
    interval_cases m
    · exact N11_B1_pos x h0 h1
    · exact N11_B2_pos x h0 h1
    · exact N11_B3_pos x h0 h1
    · exact N11_B4_pos x h0 h1
    · exact N11_B5_pos x h0 h1
    · rw [show bBlock 11 6 x = bBlock 11 5 x from bBlock_symm 11 5 (by norm_num) x]
      exact N11_B5_pos x h0 h1
    · rw [show bBlock 11 7 x = bBlock 11 4 x from bBlock_symm 11 4 (by norm_num) x]
      exact N11_B4_pos x h0 h1
    · rw [show bBlock 11 8 x = bBlock 11 3 x from bBlock_symm 11 3 (by norm_num) x]
      exact N11_B3_pos x h0 h1
    · rw [show bBlock 11 9 x = bBlock 11 2 x from bBlock_symm 11 2 (by norm_num) x]
      exact N11_B2_pos x h0 h1
    · rw [show bBlock 11 10 x = bBlock 11 1 x from bBlock_symm 11 1 (by norm_num) x]
      exact N11_B1_pos x h0 h1
  · intro hgt m hm1 hm2
    have hlo : ((3 : ℝ) / 16) < x := lt_trans (lt_trans (by norm_num) hta) hgt
    have hfn := (hs x h0 h1).2 hgt
    interval_cases m
    · exact N11_A1_pos x h0 h1
    · exact N11_A2_pos x h0 h1
    · exact N11_A3_pos x h0 hlo h1
    · exact N11_A4_pos x h0 hlo h1
    · exact N11_A5_pos x h0 h1 hfn
    · rw [show aBlock 11 6 x = aBlock 11 5 x from aBlock_symm 11 5 (by norm_num) x]
      exact N11_A5_pos x h0 h1 hfn
    · rw [show aBlock 11 7 x = aBlock 11 4 x from aBlock_symm 11 4 (by norm_num) x]
      exact N11_A4_pos x h0 hlo h1
    · rw [show aBlock 11 8 x = aBlock 11 3 x from aBlock_symm 11 3 (by norm_num) x]
      exact N11_A3_pos x h0 hlo h1
    · rw [show aBlock 11 9 x = aBlock 11 2 x from aBlock_symm 11 2 (by norm_num) x]
      exact N11_A2_pos x h0 h1
    · rw [show aBlock 11 10 x = aBlock 11 1 x from aBlock_symm 11 1 (by norm_num) x]
      exact N11_A1_pos x h0 h1

/-! ### N = 12 -/

theorem N12_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 12 1 x * (2 * (1 - x ^ 12) ^ 2) = (1 - x) ^ 4 * peval [11, 32, 62, 100, 145, 196, 252, 312, 375, 440, 506, 440, 375, 312, 252, 196, 145, 100, 62, 32, 11] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 12 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (12 - 1 : ℕ) = 11 by norm_num, show (2 * 12 - 1 : ℕ) = 23 by norm_num, show (12 + 1 : ℕ) = 13 by norm_num, show (12 - 1 + 12 : ℕ) = 23 by norm_num, show (1 + 12 : ℕ) = 13 by norm_num, show (12 + 12 : ℕ) = 24 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N12_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [11, 32, 62, 100, 145, 196, 252, 312, 375, 440, 506, 440, 375, 312, 252, 196, 145, 100, 62, 32, 11] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 11 * q ^ 20 + 252 * p ^ 1 * q ^ 19 + 2760 * p ^ 2 * q ^ 18 + 19228 * p ^ 3 * q ^ 17 + 95634 * p ^ 4 * q ^ 16 + 361284 * p ^ 5 * q ^ 15 + 1076768 * p ^ 6 * q ^ 14 + 2595780 * p ^ 7 * q ^ 13 + 5148297 * p ^ 8 * q ^ 12 + 8498776 * p ^ 9 * q ^ 11 + 11767536 * p ^ 10 * q ^ 10 + 13728660 * p ^ 11 * q ^ 9 + 13519330 * p ^ 12 * q ^ 8 + 11225412 * p ^ 13 * q ^ 7 + 7823376 * p ^ 14 * q ^ 6 + 4533144 * p ^ 15 * q ^ 5 + 2146419 * p ^ 16 * q ^ 4 + 805860 * p ^ 17 * q ^ 3 + 227480 * p ^ 18 * q ^ 2 + 43560 * p ^ 19 * q ^ 1 + 4356 * p ^ 20 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [11, 32, 62, 100, 145, 196, 252, 312, 375, 440, 506, 440, 375, 312, 252, 196, 145, 100, 62, 32, 11] x =
      (fun p q : ℝ => 11 * q ^ 20 + 252 * p ^ 1 * q ^ 19 + 2760 * p ^ 2 * q ^ 18 + 19228 * p ^ 3 * q ^ 17 + 95634 * p ^ 4 * q ^ 16 + 361284 * p ^ 5 * q ^ 15 + 1076768 * p ^ 6 * q ^ 14 + 2595780 * p ^ 7 * q ^ 13 + 5148297 * p ^ 8 * q ^ 12 + 8498776 * p ^ 9 * q ^ 11 + 11767536 * p ^ 10 * q ^ 10 + 13728660 * p ^ 11 * q ^ 9 + 13519330 * p ^ 12 * q ^ 8 + 11225412 * p ^ 13 * q ^ 7 + 7823376 * p ^ 14 * q ^ 6 + 4533144 * p ^ 15 * q ^ 5 + 2146419 * p ^ 16 * q ^ 4 + 805860 * p ^ 17 * q ^ 3 + 227480 * p ^ 18 * q ^ 2 + 43560 * p ^ 19 * q ^ 1 + 4356 * p ^ 20) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [11, 32, 62, 100, 145, 196, 252, 312, 375, 440, 506, 440, 375, 312, 252, 196, 145, 100, 62, 32, 11] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [11, 32, 62, 100, 145, 196, 252, 312, 375, 440, 506, 440, 375, 312, 252, 196, 145, 100, 62, 32, 11] x := N12_posB1_0 x hu hv

theorem N12_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 12 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N12_posB1 x h0 h1)
  rw [← N12_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N12_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 12 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 12) ^ 2) = (1 - x) ^ 4 * peval [11, 54, 137, 256, 407, 586, 789, 1012, 1251, 1502, 1761, 1892, 1761, 1502, 1251, 1012, 789, 586, 407, 256, 137, 54, 11] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 12 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (12 - 1 : ℕ) = 11 by norm_num, show (2 * 12 - 1 : ℕ) = 23 by norm_num, show (12 + 1 : ℕ) = 13 by norm_num, show (12 - 1 + 12 : ℕ) = 23 by norm_num, show (1 + 12 : ℕ) = 13 by norm_num, show (12 + 12 : ℕ) = 24 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N12_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [11, 54, 137, 256, 407, 586, 789, 1012, 1251, 1502, 1761, 1892, 1761, 1502, 1251, 1012, 789, 586, 407, 256, 137, 54, 11] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 11 * q ^ 22 + 296 * p ^ 1 * q ^ 21 + 3812 * p ^ 2 * q ^ 20 + 31276 * p ^ 3 * q ^ 19 + 183586 * p ^ 4 * q ^ 18 + 820732 * p ^ 5 * q ^ 17 + 2904440 * p ^ 6 * q ^ 16 + 8347988 * p ^ 7 * q ^ 15 + 19838489 * p ^ 8 * q ^ 14 + 39475084 * p ^ 9 * q ^ 13 + 66355828 * p ^ 10 * q ^ 12 + 94793908 * p ^ 11 * q ^ 11 + 115504114 * p ^ 12 * q ^ 10 + 120217372 * p ^ 13 * q ^ 9 + 106802344 * p ^ 14 * q ^ 8 + 80728296 * p ^ 15 * q ^ 7 + 51572499 * p ^ 16 * q ^ 6 + 27524112 * p ^ 17 * q ^ 5 + 12036596 * p ^ 18 * q ^ 4 + 4176920 * p ^ 19 * q ^ 3 + 1088516 * p ^ 20 * q ^ 2 + 191664 * p ^ 21 * q ^ 1 + 17424 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [11, 54, 137, 256, 407, 586, 789, 1012, 1251, 1502, 1761, 1892, 1761, 1502, 1251, 1012, 789, 586, 407, 256, 137, 54, 11] x =
      (fun p q : ℝ => 11 * q ^ 22 + 296 * p ^ 1 * q ^ 21 + 3812 * p ^ 2 * q ^ 20 + 31276 * p ^ 3 * q ^ 19 + 183586 * p ^ 4 * q ^ 18 + 820732 * p ^ 5 * q ^ 17 + 2904440 * p ^ 6 * q ^ 16 + 8347988 * p ^ 7 * q ^ 15 + 19838489 * p ^ 8 * q ^ 14 + 39475084 * p ^ 9 * q ^ 13 + 66355828 * p ^ 10 * q ^ 12 + 94793908 * p ^ 11 * q ^ 11 + 115504114 * p ^ 12 * q ^ 10 + 120217372 * p ^ 13 * q ^ 9 + 106802344 * p ^ 14 * q ^ 8 + 80728296 * p ^ 15 * q ^ 7 + 51572499 * p ^ 16 * q ^ 6 + 27524112 * p ^ 17 * q ^ 5 + 12036596 * p ^ 18 * q ^ 4 + 4176920 * p ^ 19 * q ^ 3 + 1088516 * p ^ 20 * q ^ 2 + 191664 * p ^ 21 * q ^ 1 + 17424 * p ^ 22) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [11, 54, 137, 256, 407, 586, 789, 1012, 1251, 1502, 1761, 1892, 1761, 1502, 1251, 1012, 789, 586, 407, 256, 137, 54, 11] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [11, 54, 137, 256, 407, 586, 789, 1012, 1251, 1502, 1761, 1892, 1761, 1502, 1251, 1012, 789, 586, 407, 256, 137, 54, 11] x := N12_posA1_0 x hu hv

theorem N12_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 12 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N12_posA1 x h0 h1)
  rw [← N12_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N12_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 12 2 x * (2 * (1 - x ^ 12) ^ 2) = (1 - x) ^ 4 * peval [20, 80, 176, 304, 460, 640, 840, 1056, 1284, 1520, 1640, 1520, 1284, 1056, 840, 640, 460, 304, 176, 80, 20] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 12 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (12 - 2 : ℕ) = 10 by norm_num, show (2 * 12 - 2 : ℕ) = 22 by norm_num, show (12 + 2 : ℕ) = 14 by norm_num, show (12 - 2 + 12 : ℕ) = 22 by norm_num, show (2 + 12 : ℕ) = 14 by norm_num, show (12 + 12 : ℕ) = 24 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N12_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [20, 80, 176, 304, 460, 640, 840, 1056, 1284, 1520, 1640, 1520, 1284, 1056, 840, 640, 460, 304, 176, 80, 20] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 20 * q ^ 20 + 480 * p ^ 1 * q ^ 19 + 5496 * p ^ 2 * q ^ 18 + 39952 * p ^ 3 * q ^ 17 + 206976 * p ^ 4 * q ^ 16 + 813120 * p ^ 5 * q ^ 15 + 2516360 * p ^ 6 * q ^ 14 + 6290064 * p ^ 7 * q ^ 13 + 12918708 * p ^ 8 * q ^ 12 + 22057024 * p ^ 9 * q ^ 11 + 31550520 * p ^ 10 * q ^ 10 + 37983120 * p ^ 11 * q ^ 9 + 38552944 * p ^ 12 * q ^ 8 + 32951808 * p ^ 13 * q ^ 7 + 23601744 * p ^ 14 * q ^ 6 + 14025120 * p ^ 15 * q ^ 5 + 6791940 * p ^ 16 * q ^ 4 + 2599200 * p ^ 17 * q ^ 3 + 744800 * p ^ 18 * q ^ 2 + 144000 * p ^ 19 * q ^ 1 + 14400 * p ^ 20 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [20, 80, 176, 304, 460, 640, 840, 1056, 1284, 1520, 1640, 1520, 1284, 1056, 840, 640, 460, 304, 176, 80, 20] x =
      (fun p q : ℝ => 20 * q ^ 20 + 480 * p ^ 1 * q ^ 19 + 5496 * p ^ 2 * q ^ 18 + 39952 * p ^ 3 * q ^ 17 + 206976 * p ^ 4 * q ^ 16 + 813120 * p ^ 5 * q ^ 15 + 2516360 * p ^ 6 * q ^ 14 + 6290064 * p ^ 7 * q ^ 13 + 12918708 * p ^ 8 * q ^ 12 + 22057024 * p ^ 9 * q ^ 11 + 31550520 * p ^ 10 * q ^ 10 + 37983120 * p ^ 11 * q ^ 9 + 38552944 * p ^ 12 * q ^ 8 + 32951808 * p ^ 13 * q ^ 7 + 23601744 * p ^ 14 * q ^ 6 + 14025120 * p ^ 15 * q ^ 5 + 6791940 * p ^ 16 * q ^ 4 + 2599200 * p ^ 17 * q ^ 3 + 744800 * p ^ 18 * q ^ 2 + 144000 * p ^ 19 * q ^ 1 + 14400 * p ^ 20) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [20, 80, 176, 304, 460, 640, 840, 1056, 1284, 1520, 1640, 1520, 1284, 1056, 840, 640, 460, 304, 176, 80, 20] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [20, 80, 176, 304, 460, 640, 840, 1056, 1284, 1520, 1640, 1520, 1284, 1056, 840, 640, 460, 304, 176, 80, 20] x := N12_posB2_0 x hu hv

theorem N12_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 12 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N12_posB2 x h0 h1)
  rw [← N12_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N12_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 12 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 12) ^ 2) = (1 - x) ^ 4 * peval [2, 48, 110, 208, 338, 496, 678, 880, 1098, 1328, 1446, 1568, 1446, 1328, 1098, 880, 678, 496, 338, 208, 110, 48, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 12 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (12 - 2 : ℕ) = 10 by norm_num, show (2 * 12 - 2 : ℕ) = 22 by norm_num, show (12 + 2 : ℕ) = 14 by norm_num, show (12 - 2 + 12 : ℕ) = 22 by norm_num, show (2 + 12 : ℕ) = 14 by norm_num, show (12 + 12 : ℕ) = 24 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N12_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 48, 110, 208, 338, 496, 678, 880, 1098, 1328, 1446, 1568, 1446, 1328, 1098, 880, 678, 496, 338, 208, 110, 48, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 22 + 92 * p ^ 1 * q ^ 21 + 1580 * p ^ 2 * q ^ 20 + 15568 * p ^ 3 * q ^ 19 + 103660 * p ^ 4 * q ^ 18 + 507496 * p ^ 5 * q ^ 17 + 1921304 * p ^ 6 * q ^ 16 + 5812400 * p ^ 7 * q ^ 15 + 14370422 * p ^ 8 * q ^ 14 + 29498260 * p ^ 9 * q ^ 13 + 50836204 * p ^ 10 * q ^ 12 + 74117776 * p ^ 11 * q ^ 11 + 91862956 * p ^ 12 * q ^ 10 + 97013608 * p ^ 13 * q ^ 9 + 87280480 * p ^ 14 * q ^ 8 + 66693120 * p ^ 15 * q ^ 7 + 42997362 * p ^ 16 * q ^ 6 + 23116092 * p ^ 17 * q ^ 5 + 10164452 * p ^ 18 * q ^ 4 + 3540800 * p ^ 19 * q ^ 3 + 925112 * p ^ 20 * q ^ 2 + 163152 * p ^ 21 * q ^ 1 + 14832 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 48, 110, 208, 338, 496, 678, 880, 1098, 1328, 1446, 1568, 1446, 1328, 1098, 880, 678, 496, 338, 208, 110, 48, 2] x =
      (fun p q : ℝ => 2 * q ^ 22 + 92 * p ^ 1 * q ^ 21 + 1580 * p ^ 2 * q ^ 20 + 15568 * p ^ 3 * q ^ 19 + 103660 * p ^ 4 * q ^ 18 + 507496 * p ^ 5 * q ^ 17 + 1921304 * p ^ 6 * q ^ 16 + 5812400 * p ^ 7 * q ^ 15 + 14370422 * p ^ 8 * q ^ 14 + 29498260 * p ^ 9 * q ^ 13 + 50836204 * p ^ 10 * q ^ 12 + 74117776 * p ^ 11 * q ^ 11 + 91862956 * p ^ 12 * q ^ 10 + 97013608 * p ^ 13 * q ^ 9 + 87280480 * p ^ 14 * q ^ 8 + 66693120 * p ^ 15 * q ^ 7 + 42997362 * p ^ 16 * q ^ 6 + 23116092 * p ^ 17 * q ^ 5 + 10164452 * p ^ 18 * q ^ 4 + 3540800 * p ^ 19 * q ^ 3 + 925112 * p ^ 20 * q ^ 2 + 163152 * p ^ 21 * q ^ 1 + 14832 * p ^ 22) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 48, 110, 208, 338, 496, 678, 880, 1098, 1328, 1446, 1568, 1446, 1328, 1098, 880, 678, 496, 338, 208, 110, 48, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 48, 110, 208, 338, 496, 678, 880, 1098, 1328, 1446, 1568, 1446, 1328, 1098, 880, 678, 496, 338, 208, 110, 48, 2] x := N12_posA2_0 x hu hv

theorem N12_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 12 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N12_posA2 x h0 h1)
  rw [← N12_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N12_linkB3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 12 3 x * (2 * (1 - x ^ 12) ^ 2) = (1 - x) ^ 4 * peval [27, 108, 270, 504, 801, 1152, 1548, 1980, 2439, 2808, 2970, 2808, 2439, 1980, 1548, 1152, 801, 504, 270, 108, 27] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 12 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (12 - 3 : ℕ) = 9 by norm_num, show (2 * 12 - 3 : ℕ) = 21 by norm_num, show (12 + 3 : ℕ) = 15 by norm_num, show (12 - 3 + 12 : ℕ) = 21 by norm_num, show (3 + 12 : ℕ) = 15 by norm_num, show (12 + 12 : ℕ) = 24 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N12_posB3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [27, 108, 270, 504, 801, 1152, 1548, 1980, 2439, 2808, 2970, 2808, 2439, 1980, 1548, 1152, 801, 504, 270, 108, 27] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 27 * q ^ 20 + 648 * p ^ 1 * q ^ 19 + 7452 * p ^ 2 * q ^ 18 + 54612 * p ^ 3 * q ^ 17 + 286146 * p ^ 4 * q ^ 16 + 1140048 * p ^ 5 * q ^ 15 + 3586212 * p ^ 6 * q ^ 14 + 9129348 * p ^ 7 * q ^ 13 + 19125153 * p ^ 8 * q ^ 12 + 33348996 * p ^ 9 * q ^ 11 + 48766212 * p ^ 10 * q ^ 10 + 60058908 * p ^ 11 * q ^ 9 + 62381826 * p ^ 12 * q ^ 8 + 54553500 * p ^ 13 * q ^ 7 + 39947256 * p ^ 14 * q ^ 6 + 24231960 * p ^ 15 * q ^ 5 + 11951955 * p ^ 16 * q ^ 4 + 4645188 * p ^ 17 * q ^ 3 + 1347192 * p ^ 18 * q ^ 2 + 262440 * p ^ 19 * q ^ 1 + 26244 * p ^ 20 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [27, 108, 270, 504, 801, 1152, 1548, 1980, 2439, 2808, 2970, 2808, 2439, 1980, 1548, 1152, 801, 504, 270, 108, 27] x =
      (fun p q : ℝ => 27 * q ^ 20 + 648 * p ^ 1 * q ^ 19 + 7452 * p ^ 2 * q ^ 18 + 54612 * p ^ 3 * q ^ 17 + 286146 * p ^ 4 * q ^ 16 + 1140048 * p ^ 5 * q ^ 15 + 3586212 * p ^ 6 * q ^ 14 + 9129348 * p ^ 7 * q ^ 13 + 19125153 * p ^ 8 * q ^ 12 + 33348996 * p ^ 9 * q ^ 11 + 48766212 * p ^ 10 * q ^ 10 + 60058908 * p ^ 11 * q ^ 9 + 62381826 * p ^ 12 * q ^ 8 + 54553500 * p ^ 13 * q ^ 7 + 39947256 * p ^ 14 * q ^ 6 + 24231960 * p ^ 15 * q ^ 5 + 11951955 * p ^ 16 * q ^ 4 + 4645188 * p ^ 17 * q ^ 3 + 1347192 * p ^ 18 * q ^ 2 + 262440 * p ^ 19 * q ^ 1 + 26244 * p ^ 20) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [27, 108, 270, 504, 801, 1152, 1548, 1980, 2439, 2808, 2970, 2808, 2439, 1980, 1548, 1152, 801, 504, 270, 108, 27] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_posB3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [27, 108, 270, 504, 801, 1152, 1548, 1980, 2439, 2808, 2970, 2808, 2439, 1980, 1548, 1152, 801, 504, 270, 108, 27] x := N12_posB3_0 x hu hv

theorem N12_B3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 12 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N12_posB3 x h0 h1)
  rw [← N12_linkB3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N12_linkA3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 12 3 x * (2 * (1 - x) ^ 2 * (1 - x ^ 12) ^ 2) = (1 - x) ^ 4 * peval [-5, 34, 113, 192, 303, 442, 605, 788, 987, 1090, 1201, 1316, 1201, 1090, 987, 788, 605, 442, 303, 192, 113, 34, -5] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 12 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (12 - 3 : ℕ) = 9 by norm_num, show (2 * 12 - 3 : ℕ) = 21 by norm_num, show (12 + 3 : ℕ) = 15 by norm_num, show (12 - 3 + 12 : ℕ) = 21 by norm_num, show (3 + 12 : ℕ) = 15 by norm_num, show (12 + 12 : ℕ) = 24 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N12_posA3_0 (x : ℝ) (hu : ((15 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-5, 34, 113, 192, 303, 442, 605, 788, 987, 1090, 1201, 1316, 1201, 1090, 987, 788, 605, 442, 303, 192, 113, 34, -5] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 70875887302052247796727157135026462847835 * q ^ 22 + 2162840274172507152233519178061167476597504 * p ^ 1 * q ^ 21 + 30485595566057662161480907340075686726172672 * p ^ 2 * q ^ 20 + 267157841529678605811577751241446149971771392 * p ^ 3 * q ^ 19 + 1645684981458220968836631200145791003977056256 * p ^ 4 * q ^ 18 + 7619061245323831719606078425526206941194354688 * p ^ 5 * q ^ 17 + 27637519635627513107657851729686790137909870592 * p ^ 6 * q ^ 16 + 80762577111387293528706503037808589246119804928 * p ^ 7 * q ^ 15 + 193830802884026507626193170121121216112967548928 * p ^ 8 * q ^ 14 + 387298246491346378358612932886266550774310895616 * p ^ 9 * q ^ 13 + 650393989959635801409070301130679607747613294592 * p ^ 10 * q ^ 12 + 923594126212595024004205974835762717465633619968 * p ^ 11 * q ^ 11 + 1112730107216481257863563947332019374687489359872 * p ^ 12 * q ^ 10 + 1138062455262644185627610323268054364866121564160 * p ^ 13 * q ^ 9 + 985964770193829499990601413716341091436746244096 * p ^ 14 * q ^ 8 + 719689983322572078435778261087120446229466578944 * p ^ 15 * q ^ 7 + 438498413860512949453375043345403710065499701248 * p ^ 16 * q ^ 6 + 219760567595932365538277973139588396161234370560 * p ^ 17 * q ^ 5 + 88567245929603806486967046105296903200395231232 * p ^ 18 * q ^ 4 + 27702587440985633753889080766762172178141020160 * p ^ 19 * q ^ 3 + 6337497241938045202997380488271492433850662912 * p ^ 20 * q ^ 2 + 947440027441162861297950798481693868719865856 * p ^ 21 * q ^ 1 + 69776941031339957563945743501528662368321536 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (15286700631942576193765185769276826401 : ℝ) * peval [-5, 34, 113, 192, 303, 442, 605, 788, 987, 1090, 1201, 1316, 1201, 1090, 987, 788, 605, 442, 303, 192, 113, 34, -5] x =
      (fun p q : ℝ => 70875887302052247796727157135026462847835 * q ^ 22 + 2162840274172507152233519178061167476597504 * p ^ 1 * q ^ 21 + 30485595566057662161480907340075686726172672 * p ^ 2 * q ^ 20 + 267157841529678605811577751241446149971771392 * p ^ 3 * q ^ 19 + 1645684981458220968836631200145791003977056256 * p ^ 4 * q ^ 18 + 7619061245323831719606078425526206941194354688 * p ^ 5 * q ^ 17 + 27637519635627513107657851729686790137909870592 * p ^ 6 * q ^ 16 + 80762577111387293528706503037808589246119804928 * p ^ 7 * q ^ 15 + 193830802884026507626193170121121216112967548928 * p ^ 8 * q ^ 14 + 387298246491346378358612932886266550774310895616 * p ^ 9 * q ^ 13 + 650393989959635801409070301130679607747613294592 * p ^ 10 * q ^ 12 + 923594126212595024004205974835762717465633619968 * p ^ 11 * q ^ 11 + 1112730107216481257863563947332019374687489359872 * p ^ 12 * q ^ 10 + 1138062455262644185627610323268054364866121564160 * p ^ 13 * q ^ 9 + 985964770193829499990601413716341091436746244096 * p ^ 14 * q ^ 8 + 719689983322572078435778261087120446229466578944 * p ^ 15 * q ^ 7 + 438498413860512949453375043345403710065499701248 * p ^ 16 * q ^ 6 + 219760567595932365538277973139588396161234370560 * p ^ 17 * q ^ 5 + 88567245929603806486967046105296903200395231232 * p ^ 18 * q ^ 4 + 27702587440985633753889080766762172178141020160 * p ^ 19 * q ^ 3 + 6337497241938045202997380488271492433850662912 * p ^ 20 * q ^ 2 + 947440027441162861297950798481693868719865856 * p ^ 21 * q ^ 1 + 69776941031339957563945743501528662368321536 * p ^ 22) (x - ((15 : ℝ) / 64)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (15286700631942576193765185769276826401 : ℝ) * peval [-5, 34, 113, 192, 303, 442, 605, 788, 987, 1090, 1201, 1316, 1201, 1090, 987, 788, 605, 442, 303, 192, 113, 34, -5] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_posA3 (x : ℝ) (hu : ((15 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-5, 34, 113, 192, 303, 442, 605, 788, 987, 1090, 1201, 1316, 1201, 1090, 987, 788, 605, 442, 303, 192, 113, 34, -5] x := N12_posA3_0 x hu hv

theorem N12_A3_pos (x : ℝ) (h0 : 0 < x) (hlo : ((15 : ℝ) / 64) < x) (h1 : x < 1) :
    0 < aBlock 12 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N12_posA3 x hlo h1)
  rw [← N12_linkA3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N12_linkB4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 12 4 x * (2 * (1 - x ^ 12) ^ 2) = (1 - x) ^ 4 * peval [32, 128, 320, 640, 1072, 1600, 2208, 2880, 3504, 3968, 4160, 3968, 3504, 2880, 2208, 1600, 1072, 640, 320, 128, 32] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 12 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (12 - 4 : ℕ) = 8 by norm_num, show (2 * 12 - 4 : ℕ) = 20 by norm_num, show (12 + 4 : ℕ) = 16 by norm_num, show (12 - 4 + 12 : ℕ) = 20 by norm_num, show (4 + 12 : ℕ) = 16 by norm_num, show (12 + 12 : ℕ) = 24 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N12_posB4_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [32, 128, 320, 640, 1072, 1600, 2208, 2880, 3504, 3968, 4160, 3968, 3504, 2880, 2208, 1600, 1072, 640, 320, 128, 32] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 32 * q ^ 20 + 768 * p ^ 1 * q ^ 19 + 8832 * p ^ 2 * q ^ 18 + 64768 * p ^ 3 * q ^ 17 + 339984 * p ^ 4 * q ^ 16 + 1359168 * p ^ 5 * q ^ 15 + 4297952 * p ^ 6 * q ^ 14 + 11020608 * p ^ 7 * q ^ 13 + 23302416 * p ^ 8 * q ^ 12 + 41094400 * p ^ 9 * q ^ 11 + 60888192 * p ^ 10 * q ^ 10 + 76101888 * p ^ 11 * q ^ 9 + 80311936 * p ^ 12 * q ^ 8 + 71397888 * p ^ 13 * q ^ 7 + 53137920 * p ^ 14 * q ^ 6 + 32729088 * p ^ 15 * q ^ 5 + 16363776 * p ^ 16 * q ^ 4 + 6432768 * p ^ 17 * q ^ 3 + 1882112 * p ^ 18 * q ^ 2 + 368640 * p ^ 19 * q ^ 1 + 36864 * p ^ 20 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [32, 128, 320, 640, 1072, 1600, 2208, 2880, 3504, 3968, 4160, 3968, 3504, 2880, 2208, 1600, 1072, 640, 320, 128, 32] x =
      (fun p q : ℝ => 32 * q ^ 20 + 768 * p ^ 1 * q ^ 19 + 8832 * p ^ 2 * q ^ 18 + 64768 * p ^ 3 * q ^ 17 + 339984 * p ^ 4 * q ^ 16 + 1359168 * p ^ 5 * q ^ 15 + 4297952 * p ^ 6 * q ^ 14 + 11020608 * p ^ 7 * q ^ 13 + 23302416 * p ^ 8 * q ^ 12 + 41094400 * p ^ 9 * q ^ 11 + 60888192 * p ^ 10 * q ^ 10 + 76101888 * p ^ 11 * q ^ 9 + 80311936 * p ^ 12 * q ^ 8 + 71397888 * p ^ 13 * q ^ 7 + 53137920 * p ^ 14 * q ^ 6 + 32729088 * p ^ 15 * q ^ 5 + 16363776 * p ^ 16 * q ^ 4 + 6432768 * p ^ 17 * q ^ 3 + 1882112 * p ^ 18 * q ^ 2 + 368640 * p ^ 19 * q ^ 1 + 36864 * p ^ 20) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [32, 128, 320, 640, 1072, 1600, 2208, 2880, 3504, 3968, 4160, 3968, 3504, 2880, 2208, 1600, 1072, 640, 320, 128, 32] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_posB4 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [32, 128, 320, 640, 1072, 1600, 2208, 2880, 3504, 3968, 4160, 3968, 3504, 2880, 2208, 1600, 1072, 640, 320, 128, 32] x := N12_posB4_0 x hu hv

theorem N12_B4_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 12 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N12_posB4 x h0 h1)
  rw [← N12_linkB4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N12_linkA4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 12 4 x * (2 * (1 - x) ^ 2 * (1 - x ^ 12) ^ 2) = (1 - x) ^ 4 * peval [-10, 24, 98, 208, 302, 424, 570, 736, 822, 920, 1026, 1136, 1026, 920, 822, 736, 570, 424, 302, 208, 98, 24, -10] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 12 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (12 - 4 : ℕ) = 8 by norm_num, show (2 * 12 - 4 : ℕ) = 20 by norm_num, show (12 + 4 : ℕ) = 16 by norm_num, show (12 - 4 + 12 : ℕ) = 20 by norm_num, show (4 + 12 : ℕ) = 16 by norm_num, show (12 + 12 : ℕ) = 24 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N12_posA4_0 (x : ℝ) (hu : ((15 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-10, 24, 98, 208, 302, 424, 570, 736, 822, 920, 1026, 1136, 1026, 920, 822, 736, 570, 424, 302, 208, 98, 24, -10] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 27388998654440839177681998018562563080630 * q ^ 22 + 1143615533394856283910622990792825714352384 * p ^ 1 * q ^ 21 + 19095886836638244330607702737719893120335872 * p ^ 2 * q ^ 20 + 186383573535379101007801399993685070428766208 * p ^ 3 * q ^ 19 + 1237305724945086068902664442718781107710459904 * p ^ 4 * q ^ 18 + 6050721554133361827401750112645843825687789568 * p ^ 5 * q ^ 17 + 22874281697341358061603320335974514046346985472 * p ^ 6 * q ^ 16 + 68992188848122428891793909391520089382362021888 * p ^ 7 * q ^ 15 + 169644145795160642446582474597468145043170656256 * p ^ 8 * q ^ 14 + 345224454204937336833602512834313087249903779840 * p ^ 9 * q ^ 13 + 587500102075439353400740941717103122751959334912 * p ^ 10 * q ^ 12 + 841831659473828076581365029678233014879500894208 * p ^ 11 * q ^ 11 + 1019576357409293062948478642790987255776523845632 * p ^ 12 * q ^ 10 + 1044868065394381085659387313199265652587047682048 * p ^ 13 * q ^ 9 + 904492009337595442626188043158381112927641927680 * p ^ 14 * q ^ 8 + 658169101402943304499885753782818195371764744192 * p ^ 15 * q ^ 7 + 399076474110929278887761037824357193038633631744 * p ^ 16 * q ^ 6 + 198821945975094611596414681346447408363722506240 * p ^ 17 * q ^ 5 + 79627937402040677619285598559316190969695567872 * p ^ 18 * q ^ 4 + 24762415511833957396296859715627211619995484160 * p ^ 19 * q ^ 3 + 5639391829371871116323570428907508144665526272 * p ^ 20 * q ^ 2 + 840986091773616472388068686292739501447970816 * p ^ 21 * q ^ 1 + 61936835297481535365749592546300722776375296 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (15286700631942576193765185769276826401 : ℝ) * peval [-10, 24, 98, 208, 302, 424, 570, 736, 822, 920, 1026, 1136, 1026, 920, 822, 736, 570, 424, 302, 208, 98, 24, -10] x =
      (fun p q : ℝ => 27388998654440839177681998018562563080630 * q ^ 22 + 1143615533394856283910622990792825714352384 * p ^ 1 * q ^ 21 + 19095886836638244330607702737719893120335872 * p ^ 2 * q ^ 20 + 186383573535379101007801399993685070428766208 * p ^ 3 * q ^ 19 + 1237305724945086068902664442718781107710459904 * p ^ 4 * q ^ 18 + 6050721554133361827401750112645843825687789568 * p ^ 5 * q ^ 17 + 22874281697341358061603320335974514046346985472 * p ^ 6 * q ^ 16 + 68992188848122428891793909391520089382362021888 * p ^ 7 * q ^ 15 + 169644145795160642446582474597468145043170656256 * p ^ 8 * q ^ 14 + 345224454204937336833602512834313087249903779840 * p ^ 9 * q ^ 13 + 587500102075439353400740941717103122751959334912 * p ^ 10 * q ^ 12 + 841831659473828076581365029678233014879500894208 * p ^ 11 * q ^ 11 + 1019576357409293062948478642790987255776523845632 * p ^ 12 * q ^ 10 + 1044868065394381085659387313199265652587047682048 * p ^ 13 * q ^ 9 + 904492009337595442626188043158381112927641927680 * p ^ 14 * q ^ 8 + 658169101402943304499885753782818195371764744192 * p ^ 15 * q ^ 7 + 399076474110929278887761037824357193038633631744 * p ^ 16 * q ^ 6 + 198821945975094611596414681346447408363722506240 * p ^ 17 * q ^ 5 + 79627937402040677619285598559316190969695567872 * p ^ 18 * q ^ 4 + 24762415511833957396296859715627211619995484160 * p ^ 19 * q ^ 3 + 5639391829371871116323570428907508144665526272 * p ^ 20 * q ^ 2 + 840986091773616472388068686292739501447970816 * p ^ 21 * q ^ 1 + 61936835297481535365749592546300722776375296 * p ^ 22) (x - ((15 : ℝ) / 64)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (15286700631942576193765185769276826401 : ℝ) * peval [-10, 24, 98, 208, 302, 424, 570, 736, 822, 920, 1026, 1136, 1026, 920, 822, 736, 570, 424, 302, 208, 98, 24, -10] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_posA4 (x : ℝ) (hu : ((15 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-10, 24, 98, 208, 302, 424, 570, 736, 822, 920, 1026, 1136, 1026, 920, 822, 736, 570, 424, 302, 208, 98, 24, -10] x := N12_posA4_0 x hu hv

theorem N12_A4_pos (x : ℝ) (h0 : 0 < x) (hlo : ((15 : ℝ) / 64) < x) (h1 : x < 1) :
    0 < aBlock 12 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N12_posA4 x hlo h1)
  rw [← N12_linkA4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N12_linkB5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 12 5 x * (2 * (1 - x ^ 12) ^ 2) = (1 - x) ^ 4 * peval [35, 140, 350, 700, 1225, 1900, 2700, 3516, 4239, 4760, 4970, 4760, 4239, 3516, 2700, 1900, 1225, 700, 350, 140, 35] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 12 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (12 - 5 : ℕ) = 7 by norm_num, show (2 * 12 - 5 : ℕ) = 19 by norm_num, show (12 + 5 : ℕ) = 17 by norm_num, show (12 - 5 + 12 : ℕ) = 19 by norm_num, show (5 + 12 : ℕ) = 17 by norm_num, show (12 + 12 : ℕ) = 24 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N12_posB5_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [35, 140, 350, 700, 1225, 1900, 2700, 3516, 4239, 4760, 4970, 4760, 4239, 3516, 2700, 1900, 1225, 700, 350, 140, 35] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 35 * q ^ 20 + 840 * p ^ 1 * q ^ 19 + 9660 * p ^ 2 * q ^ 18 + 70840 * p ^ 3 * q ^ 17 + 371910 * p ^ 4 * q ^ 16 + 1487580 * p ^ 5 * q ^ 15 + 4709720 * p ^ 6 * q ^ 14 + 12103296 * p ^ 7 * q ^ 13 + 25681917 * p ^ 8 * q ^ 12 + 45518656 * p ^ 9 * q ^ 11 + 67889160 * p ^ 10 * q ^ 10 + 85537200 * p ^ 11 * q ^ 9 + 91102606 * p ^ 12 * q ^ 8 + 81793152 * p ^ 13 * q ^ 7 + 61482036 * p ^ 14 * q ^ 6 + 38225880 * p ^ 15 * q ^ 5 + 19272435 * p ^ 16 * q ^ 4 + 7629300 * p ^ 17 * q ^ 3 + 2244200 * p ^ 18 * q ^ 2 + 441000 * p ^ 19 * q ^ 1 + 44100 * p ^ 20 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [35, 140, 350, 700, 1225, 1900, 2700, 3516, 4239, 4760, 4970, 4760, 4239, 3516, 2700, 1900, 1225, 700, 350, 140, 35] x =
      (fun p q : ℝ => 35 * q ^ 20 + 840 * p ^ 1 * q ^ 19 + 9660 * p ^ 2 * q ^ 18 + 70840 * p ^ 3 * q ^ 17 + 371910 * p ^ 4 * q ^ 16 + 1487580 * p ^ 5 * q ^ 15 + 4709720 * p ^ 6 * q ^ 14 + 12103296 * p ^ 7 * q ^ 13 + 25681917 * p ^ 8 * q ^ 12 + 45518656 * p ^ 9 * q ^ 11 + 67889160 * p ^ 10 * q ^ 10 + 85537200 * p ^ 11 * q ^ 9 + 91102606 * p ^ 12 * q ^ 8 + 81793152 * p ^ 13 * q ^ 7 + 61482036 * p ^ 14 * q ^ 6 + 38225880 * p ^ 15 * q ^ 5 + 19272435 * p ^ 16 * q ^ 4 + 7629300 * p ^ 17 * q ^ 3 + 2244200 * p ^ 18 * q ^ 2 + 441000 * p ^ 19 * q ^ 1 + 44100 * p ^ 20) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [35, 140, 350, 700, 1225, 1900, 2700, 3516, 4239, 4760, 4970, 4760, 4239, 3516, 2700, 1900, 1225, 700, 350, 140, 35] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_posB5 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [35, 140, 350, 700, 1225, 1900, 2700, 3516, 4239, 4760, 4970, 4760, 4239, 3516, 2700, 1900, 1225, 700, 350, 140, 35] x := N12_posB5_0 x hu hv

theorem N12_B5_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 12 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N12_posB5 x h0 h1)
  rw [← N12_linkB5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N12_linkA5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 12 5 x * (2 * (1 - x) ^ 2 * (1 - x ^ 12) ^ 2) = (1 - x) ^ 4 * peval [-13, 18, 89, 196, 335, 442, 573, 640, 723, 818, 921, 1028, 921, 818, 723, 640, 573, 442, 335, 196, 89, 18, -13] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 12 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (12 - 5 : ℕ) = 7 by norm_num, show (2 * 12 - 5 : ℕ) = 19 by norm_num, show (12 + 5 : ℕ) = 17 by norm_num, show (12 - 5 + 12 : ℕ) = 19 by norm_num, show (5 + 12 : ℕ) = 17 by norm_num, show (12 + 12 : ℕ) = 24 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N12_posA5_0 (x : ℝ) (hu : ((15 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-13, 18, 89, 196, 335, 442, 573, 640, 723, 818, 921, 1028, 921, 818, 723, 640, 573, 442, 335, 196, 89, 18, -13] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 453644639665389971287316886481443162707 * q ^ 22 + 507827155487258508119653729233198180799232 * p ^ 1 * q ^ 21 + 11946622918399831022007035303490120707047424 * p ^ 2 * q ^ 20 + 135466144227040809730610116837239749424971776 * p ^ 3 * q ^ 19 + 979626330093277450900000609187692626023809024 * p ^ 4 * q ^ 18 + 5064759454099279360102318211174114054348210176 * p ^ 5 * q ^ 17 + 19908539792803183107927632758904643876514955264 * p ^ 6 * q ^ 16 + 61784581965636985793008772943584681212805906432 * p ^ 7 * q ^ 15 + 155184315697406274275471561904934483415464411136 * p ^ 8 * q ^ 14 + 320821287374997840692528166077767381553829117952 * p ^ 9 * q ^ 13 + 552224518047859219280620531830082161109865857024 * p ^ 10 * q ^ 12 + 797396130298229682748489288602550727159427629056 * p ^ 11 * q ^ 11 + 970103773544944851963714238711867871064771002368 * p ^ 12 * q ^ 10 + 995842694547117665891534524826319540510200430592 * p ^ 13 * q ^ 9 + 861408285346480850481902668072260574129665605632 * p ^ 14 * q ^ 8 + 625079126512246547362211598019748653438787911680 * p ^ 15 * q ^ 7 + 377368188927639378259780973188407239952052518912 * p ^ 16 * q ^ 6 + 186999797395021659954233384946101097949500014592 * p ^ 17 * q ^ 5 + 74464100947613352205870520530671267870326915072 * p ^ 18 * q ^ 4 + 23032349450339637334552157449563468874750361600 * p ^ 19 * q ^ 3 + 5223286017457214624040651916498614975884951552 * p ^ 20 * q ^ 2 + 777113730373088639042139418979366881084833792 * p ^ 21 * q ^ 1 + 57232771857166482046831901973163959021207552 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (15286700631942576193765185769276826401 : ℝ) * peval [-13, 18, 89, 196, 335, 442, 573, 640, 723, 818, 921, 1028, 921, 818, 723, 640, 573, 442, 335, 196, 89, 18, -13] x =
      (fun p q : ℝ => 453644639665389971287316886481443162707 * q ^ 22 + 507827155487258508119653729233198180799232 * p ^ 1 * q ^ 21 + 11946622918399831022007035303490120707047424 * p ^ 2 * q ^ 20 + 135466144227040809730610116837239749424971776 * p ^ 3 * q ^ 19 + 979626330093277450900000609187692626023809024 * p ^ 4 * q ^ 18 + 5064759454099279360102318211174114054348210176 * p ^ 5 * q ^ 17 + 19908539792803183107927632758904643876514955264 * p ^ 6 * q ^ 16 + 61784581965636985793008772943584681212805906432 * p ^ 7 * q ^ 15 + 155184315697406274275471561904934483415464411136 * p ^ 8 * q ^ 14 + 320821287374997840692528166077767381553829117952 * p ^ 9 * q ^ 13 + 552224518047859219280620531830082161109865857024 * p ^ 10 * q ^ 12 + 797396130298229682748489288602550727159427629056 * p ^ 11 * q ^ 11 + 970103773544944851963714238711867871064771002368 * p ^ 12 * q ^ 10 + 995842694547117665891534524826319540510200430592 * p ^ 13 * q ^ 9 + 861408285346480850481902668072260574129665605632 * p ^ 14 * q ^ 8 + 625079126512246547362211598019748653438787911680 * p ^ 15 * q ^ 7 + 377368188927639378259780973188407239952052518912 * p ^ 16 * q ^ 6 + 186999797395021659954233384946101097949500014592 * p ^ 17 * q ^ 5 + 74464100947613352205870520530671267870326915072 * p ^ 18 * q ^ 4 + 23032349450339637334552157449563468874750361600 * p ^ 19 * q ^ 3 + 5223286017457214624040651916498614975884951552 * p ^ 20 * q ^ 2 + 777113730373088639042139418979366881084833792 * p ^ 21 * q ^ 1 + 57232771857166482046831901973163959021207552 * p ^ 22) (x - ((15 : ℝ) / 64)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (15286700631942576193765185769276826401 : ℝ) * peval [-13, 18, 89, 196, 335, 442, 573, 640, 723, 818, 921, 1028, 921, 818, 723, 640, 573, 442, 335, 196, 89, 18, -13] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_posA5 (x : ℝ) (hu : ((15 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-13, 18, 89, 196, 335, 442, 573, 640, 723, 818, 921, 1028, 921, 818, 723, 640, 573, 442, 335, 196, 89, 18, -13] x := N12_posA5_0 x hu hv

theorem N12_A5_pos (x : ℝ) (h0 : 0 < x) (hlo : ((15 : ℝ) / 64) < x) (h1 : x < 1) :
    0 < aBlock 12 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N12_posA5 x hlo h1)
  rw [← N12_linkA5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N12_linkB6 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 12 6 x * (2 * (1 - x ^ 12) ^ 2) = (1 - x) ^ 4 * peval [36, 144, 360, 720, 1260, 2016, 2880, 3744, 4500, 5040, 5256, 5040, 4500, 3744, 2880, 2016, 1260, 720, 360, 144, 36] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 12 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (12 - 6 : ℕ) = 6 by norm_num, show (2 * 12 - 6 : ℕ) = 18 by norm_num, show (12 + 6 : ℕ) = 18 by norm_num, show (12 - 6 + 12 : ℕ) = 18 by norm_num, show (6 + 12 : ℕ) = 18 by norm_num, show (12 + 12 : ℕ) = 24 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N12_posB6_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [36, 144, 360, 720, 1260, 2016, 2880, 3744, 4500, 5040, 5256, 5040, 4500, 3744, 2880, 2016, 1260, 720, 360, 144, 36] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 36 * q ^ 20 + 864 * p ^ 1 * q ^ 19 + 9936 * p ^ 2 * q ^ 18 + 72864 * p ^ 3 * q ^ 17 + 382536 * p ^ 4 * q ^ 16 + 1530144 * p ^ 5 * q ^ 15 + 4845312 * p ^ 6 * q ^ 14 + 12457152 * p ^ 7 * q ^ 13 + 26454924 * p ^ 8 * q ^ 12 + 46952640 * p ^ 9 * q ^ 11 + 70164576 * p ^ 10 * q ^ 10 + 88627392 * p ^ 11 * q ^ 9 + 94676616 * p ^ 12 * q ^ 8 + 85281120 * p ^ 13 * q ^ 7 + 64318320 * p ^ 14 * q ^ 6 + 40116384 * p ^ 15 * q ^ 5 + 20282724 * p ^ 16 * q ^ 4 + 8048160 * p ^ 17 * q ^ 3 + 2371680 * p ^ 18 * q ^ 2 + 466560 * p ^ 19 * q ^ 1 + 46656 * p ^ 20 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [36, 144, 360, 720, 1260, 2016, 2880, 3744, 4500, 5040, 5256, 5040, 4500, 3744, 2880, 2016, 1260, 720, 360, 144, 36] x =
      (fun p q : ℝ => 36 * q ^ 20 + 864 * p ^ 1 * q ^ 19 + 9936 * p ^ 2 * q ^ 18 + 72864 * p ^ 3 * q ^ 17 + 382536 * p ^ 4 * q ^ 16 + 1530144 * p ^ 5 * q ^ 15 + 4845312 * p ^ 6 * q ^ 14 + 12457152 * p ^ 7 * q ^ 13 + 26454924 * p ^ 8 * q ^ 12 + 46952640 * p ^ 9 * q ^ 11 + 70164576 * p ^ 10 * q ^ 10 + 88627392 * p ^ 11 * q ^ 9 + 94676616 * p ^ 12 * q ^ 8 + 85281120 * p ^ 13 * q ^ 7 + 64318320 * p ^ 14 * q ^ 6 + 40116384 * p ^ 15 * q ^ 5 + 20282724 * p ^ 16 * q ^ 4 + 8048160 * p ^ 17 * q ^ 3 + 2371680 * p ^ 18 * q ^ 2 + 466560 * p ^ 19 * q ^ 1 + 46656 * p ^ 20) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [36, 144, 360, 720, 1260, 2016, 2880, 3744, 4500, 5040, 5256, 5040, 4500, 3744, 2880, 2016, 1260, 720, 360, 144, 36] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_posB6 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [36, 144, 360, 720, 1260, 2016, 2880, 3744, 4500, 5040, 5256, 5040, 4500, 3744, 2880, 2016, 1260, 720, 360, 144, 36] x := N12_posB6_0 x hu hv

theorem N12_B6_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 12 6 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N12_posB6 x h0 h1)
  rw [← N12_linkB6 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N12_linkA6 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 12 6 x * (2 * (1 - x) ^ 2 * (1 - x ^ 12) ^ 2) = (1 - x) ^ 4 * peval [-14, 16, 86, 192, 330, 496, 542, 608, 690, 784, 886, 992, 886, 784, 690, 608, 542, 496, 330, 192, 86, 16, -14] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 12 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (12 - 6 : ℕ) = 6 by norm_num, show (2 * 12 - 6 : ℕ) = 18 by norm_num, show (12 + 6 : ℕ) = 18 by norm_num, show (12 - 6 + 12 : ℕ) = 18 by norm_num, show (6 + 12 : ℕ) = 18 by norm_num, show (12 + 12 : ℕ) = 24 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

/-- The minimal polynomial `f_12` of `ξ*(12)` (constant term first). -/
def f12 : List ℤ := [7, -8, -43, -96, -165, -248, -278, -296, -302, -296, -278, -248, -165, -96, -43, -8, 7]

theorem N12_crit_factor (x : ℝ) : peval [-14, 16, 86, 192, 330, 496, 542, 608, 690, 784, 886, 992, 886, 784, 690, 608, 542, 496, 330, 192, 86, 16, -14] x = -(peval [2, 0, 0, 0, 0, 0, 2] x * peval f12 x) := by
  simp only [f12, peval]; ring

theorem N12_posH_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 0, 0, 0, 0, 0, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 6 + 12 * p ^ 1 * q ^ 5 + 30 * p ^ 2 * q ^ 4 + 40 * p ^ 3 * q ^ 3 + 30 * p ^ 4 * q ^ 2 + 12 * p ^ 5 * q ^ 1 + 4 * p ^ 6 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 0, 0, 0, 0, 0, 2] x =
      (fun p q : ℝ => 2 * q ^ 6 + 12 * p ^ 1 * q ^ 5 + 30 * p ^ 2 * q ^ 4 + 40 * p ^ 3 * q ^ 3 + 30 * p ^ 4 * q ^ 2 + 12 * p ^ 5 * q ^ 1 + 4 * p ^ 6) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 0, 0, 0, 0, 0, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_posH (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 0, 0, 0, 0, 0, 2] x := N12_posH_0 x hu hv

theorem N12_dneg_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [8, 86, 288, 660, 1240, 1668, 2072, 2416, 2664, 2780, 2728, 1980, 1248, 602, 120, -112] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 8 * q ^ 15 + 206 * p ^ 1 * q ^ 14 + 2332 * p ^ 2 * q ^ 13 + 15870 * p ^ 3 * q ^ 12 + 73848 * p ^ 4 * q ^ 11 + 251346 * p ^ 5 * q ^ 10 + 650284 * p ^ 6 * q ^ 9 + 1307818 * p ^ 7 * q ^ 8 + 2069504 * p ^ 8 * q ^ 7 + 2588630 * p ^ 9 * q ^ 6 + 2554508 * p ^ 10 * q ^ 5 + 1969858 * p ^ 11 * q ^ 4 + 1161904 * p ^ 12 * q ^ 3 + 501340 * p ^ 13 * q ^ 2 + 143240 * p ^ 14 * q ^ 1 + 20448 * p ^ 15 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [8, 86, 288, 660, 1240, 1668, 2072, 2416, 2664, 2780, 2728, 1980, 1248, 602, 120, -112] x =
      (fun p q : ℝ => 8 * q ^ 15 + 206 * p ^ 1 * q ^ 14 + 2332 * p ^ 2 * q ^ 13 + 15870 * p ^ 3 * q ^ 12 + 73848 * p ^ 4 * q ^ 11 + 251346 * p ^ 5 * q ^ 10 + 650284 * p ^ 6 * q ^ 9 + 1307818 * p ^ 7 * q ^ 8 + 2069504 * p ^ 8 * q ^ 7 + 2588630 * p ^ 9 * q ^ 6 + 2554508 * p ^ 10 * q ^ 5 + 1969858 * p ^ 11 * q ^ 4 + 1161904 * p ^ 12 * q ^ 3 + 501340 * p ^ 13 * q ^ 2 + 143240 * p ^ 14 * q ^ 1 + 20448 * p ^ 15) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [8, 86, 288, 660, 1240, 1668, 2072, 2416, 2664, 2780, 2728, 1980, 1248, 602, 120, -112] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N12_dneg (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [8, 86, 288, 660, 1240, 1668, 2072, 2416, 2664, 2780, 2728, 1980, 1248, 602, 120, -112] x := N12_dneg_0 x hu hv

theorem N12_dlink (x : ℝ) : dpeval f12 x = -peval [8, 86, 288, 660, 1240, 1668, 2072, 2416, 2664, 2780, 2728, 1980, 1248, 602, 120, -112] x := by
  simp only [f12, dpeval, peval]; ring

/-- `f_12` changes sign exactly once in `(0,1)`, at a point of `(24741777680406770917/100000000000000000000, 12370888840203385459/50000000000000000000)`. -/
theorem N12_root : ∃ ξs, ((24741777680406770917 : ℝ) / 100000000000000000000) < ξs ∧ ξs < ((12370888840203385459 : ℝ) / 50000000000000000000) ∧ peval f12 ξs = 0 ∧
    ∀ x, 0 < x → x < 1 → (x < ξs → 0 < peval f12 x) ∧ (ξs < x → peval f12 x < 0) := by
  apply sign_change (peval f12) (0 : ℝ) _ _ (by norm_num) (by norm_num) (by norm_num)
    (continuous_peval _).continuousOn
  · apply peval_strictAntiOn
    intro x hu hv
    rw [N12_dlink]
    linarith [N12_dneg x hu hv]
  · exact (fun x hx0 hx => absurd hx (not_le.mpr hx0))
  · norm_num [f12, peval]
  · norm_num [f12, peval]

theorem N12_A6_neg (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : 0 < peval f12 x) :
    aBlock 12 6 x < 0 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : (1 - x) ^ 4 * peval [-14, 16, 86, 192, 330, 496, 542, 608, 690, 784, 886, 992, 886, 784, 690, 608, 542, 496, 330, 192, 86, 16, -14] x < 0 := by
    rw [N12_crit_factor]
    have := mul_pos (N12_posH x h0 h1) hf
    nlinarith [pow_pos h1x 4]
  rw [← N12_linkA6 x h0 h1] at hn
  have hd := (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))
  by_contra hc
  push Not at hc
  nlinarith [mul_nonneg hc hd.le]

theorem N12_A6_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : peval f12 x < 0) :
    0 < aBlock 12 6 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 12 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : 0 < (1 - x) ^ 4 * peval [-14, 16, 86, 192, 330, 496, 542, 608, 690, 784, 886, 992, 886, 784, 690, 608, 542, 496, 330, 192, 86, 16, -14] x := by
    rw [N12_crit_factor]
    have := mul_pos (N12_posH x h0 h1) (neg_pos.mpr hf)
    have h4 := pow_pos h1x 4
    nlinarith
  rw [← N12_linkA6 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

/-- Note Thm 6(b) for `N = 12`: there is a unique `ξ* ∈ (0,1)` (the unique root of `f_12` in
`(0,1)`, enclosed in `(24741777680406770917/100000000000000000000, 12370888840203385459/50000000000000000000)`) such that, for every `ξ ∈ (0,1)`: all `b_m > 0`;
`a_6 < 0` for `ξ < ξ*` (spectrally unstable, since `b_6 > 0`); and all `a_m > 0` for `ξ > ξ*`. -/
theorem cert_N12 : ∃ ξs : ℝ, ((24741777680406770917 : ℝ) / 100000000000000000000) < ξs ∧ ξs < ((12370888840203385459 : ℝ) / 50000000000000000000) ∧ peval f12 ξs = 0 ∧
    (∀ y, 0 < y → y < 1 → peval f12 y = 0 → y = ξs) ∧
    ∀ x, 0 < x → x < 1 →
      (∀ m, 1 ≤ m → m ≤ 11 → 0 < bBlock 12 m x) ∧
      (x < ξs → aBlock 12 6 x < 0) ∧
      (ξs < x → ∀ m, 1 ≤ m → m ≤ 11 → 0 < aBlock 12 m x) := by
  obtain ⟨ξs, hta, htb, hf0, hs⟩ := N12_root
  refine ⟨ξs, hta, htb, hf0, ?_, ?_⟩
  · intro y hy0 hy1 hy
    rcases lt_trichotomy y ξs with h | h | h
    · linarith [(hs y hy0 hy1).1 h]
    · exact h
    · linarith [(hs y hy0 hy1).2 h]
  intro x h0 h1
  refine ⟨?_, fun hlt => N12_A6_neg x h0 h1 ((hs x h0 h1).1 hlt), ?_⟩
  · intro m hm1 hm2
    interval_cases m
    · exact N12_B1_pos x h0 h1
    · exact N12_B2_pos x h0 h1
    · exact N12_B3_pos x h0 h1
    · exact N12_B4_pos x h0 h1
    · exact N12_B5_pos x h0 h1
    · exact N12_B6_pos x h0 h1
    · rw [show bBlock 12 7 x = bBlock 12 5 x from bBlock_symm 12 5 (by norm_num) x]
      exact N12_B5_pos x h0 h1
    · rw [show bBlock 12 8 x = bBlock 12 4 x from bBlock_symm 12 4 (by norm_num) x]
      exact N12_B4_pos x h0 h1
    · rw [show bBlock 12 9 x = bBlock 12 3 x from bBlock_symm 12 3 (by norm_num) x]
      exact N12_B3_pos x h0 h1
    · rw [show bBlock 12 10 x = bBlock 12 2 x from bBlock_symm 12 2 (by norm_num) x]
      exact N12_B2_pos x h0 h1
    · rw [show bBlock 12 11 x = bBlock 12 1 x from bBlock_symm 12 1 (by norm_num) x]
      exact N12_B1_pos x h0 h1
  · intro hgt m hm1 hm2
    have hlo : ((15 : ℝ) / 64) < x := lt_trans (lt_trans (by norm_num) hta) hgt
    have hfn := (hs x h0 h1).2 hgt
    interval_cases m
    · exact N12_A1_pos x h0 h1
    · exact N12_A2_pos x h0 h1
    · exact N12_A3_pos x h0 hlo h1
    · exact N12_A4_pos x h0 hlo h1
    · exact N12_A5_pos x h0 hlo h1
    · exact N12_A6_pos x h0 h1 hfn
    · rw [show aBlock 12 7 x = aBlock 12 5 x from aBlock_symm 12 5 (by norm_num) x]
      exact N12_A5_pos x h0 hlo h1
    · rw [show aBlock 12 8 x = aBlock 12 4 x from aBlock_symm 12 4 (by norm_num) x]
      exact N12_A4_pos x h0 hlo h1
    · rw [show aBlock 12 9 x = aBlock 12 3 x from aBlock_symm 12 3 (by norm_num) x]
      exact N12_A3_pos x h0 hlo h1
    · rw [show aBlock 12 10 x = aBlock 12 2 x from aBlock_symm 12 2 (by norm_num) x]
      exact N12_A2_pos x h0 h1
    · rw [show aBlock 12 11 x = aBlock 12 1 x from aBlock_symm 12 1 (by norm_num) x]
      exact N12_A1_pos x h0 h1

/-! ### N = 13 -/

theorem N13_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 13 1 x * (2 * (1 - x ^ 13) ^ 2) = (1 - x) ^ 4 * peval [12, 35, 68, 110, 160, 217, 280, 348, 420, 495, 572, 650, 572, 495, 420, 348, 280, 217, 160, 110, 68, 35, 12] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 13 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (13 - 1 : ℕ) = 12 by norm_num, show (2 * 13 - 1 : ℕ) = 25 by norm_num, show (13 + 1 : ℕ) = 14 by norm_num, show (13 - 1 + 13 : ℕ) = 25 by norm_num, show (1 + 13 : ℕ) = 14 by norm_num, show (13 + 13 : ℕ) = 26 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N13_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 35, 68, 110, 160, 217, 280, 348, 420, 495, 572, 650, 572, 495, 420, 348, 280, 217, 160, 110, 68, 35, 12] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 12 * q ^ 22 + 299 * p ^ 1 * q ^ 21 + 3575 * p ^ 2 * q ^ 20 + 27300 * p ^ 3 * q ^ 19 + 149500 * p ^ 4 * q ^ 18 + 624910 * p ^ 5 * q ^ 17 + 2072070 * p ^ 6 * q ^ 16 + 5591300 * p ^ 7 * q ^ 15 + 12498200 * p ^ 8 * q ^ 14 + 23434125 * p ^ 9 * q ^ 13 + 37182145 * p ^ 10 * q ^ 12 + 50220040 * p ^ 11 * q ^ 11 + 57946044 * p ^ 12 * q ^ 10 + 57201430 * p ^ 13 * q ^ 9 + 48278230 * p ^ 14 * q ^ 8 + 34733556 * p ^ 15 * q ^ 7 + 21170292 * p ^ 16 * q ^ 6 + 10813803 * p ^ 17 * q ^ 5 + 4545255 * p ^ 18 * q ^ 4 + 1524380 * p ^ 19 * q ^ 3 + 386672 * p ^ 20 * q ^ 2 + 66924 * p ^ 21 * q ^ 1 + 6084 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [12, 35, 68, 110, 160, 217, 280, 348, 420, 495, 572, 650, 572, 495, 420, 348, 280, 217, 160, 110, 68, 35, 12] x =
      (fun p q : ℝ => 12 * q ^ 22 + 299 * p ^ 1 * q ^ 21 + 3575 * p ^ 2 * q ^ 20 + 27300 * p ^ 3 * q ^ 19 + 149500 * p ^ 4 * q ^ 18 + 624910 * p ^ 5 * q ^ 17 + 2072070 * p ^ 6 * q ^ 16 + 5591300 * p ^ 7 * q ^ 15 + 12498200 * p ^ 8 * q ^ 14 + 23434125 * p ^ 9 * q ^ 13 + 37182145 * p ^ 10 * q ^ 12 + 50220040 * p ^ 11 * q ^ 11 + 57946044 * p ^ 12 * q ^ 10 + 57201430 * p ^ 13 * q ^ 9 + 48278230 * p ^ 14 * q ^ 8 + 34733556 * p ^ 15 * q ^ 7 + 21170292 * p ^ 16 * q ^ 6 + 10813803 * p ^ 17 * q ^ 5 + 4545255 * p ^ 18 * q ^ 4 + 1524380 * p ^ 19 * q ^ 3 + 386672 * p ^ 20 * q ^ 2 + 66924 * p ^ 21 * q ^ 1 + 6084 * p ^ 22) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [12, 35, 68, 110, 160, 217, 280, 348, 420, 495, 572, 650, 572, 495, 420, 348, 280, 217, 160, 110, 68, 35, 12] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N13_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 35, 68, 110, 160, 217, 280, 348, 420, 495, 572, 650, 572, 495, 420, 348, 280, 217, 160, 110, 68, 35, 12] x := N13_posB1_0 x hu hv

theorem N13_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 13 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N13_posB1 x h0 h1)
  rw [← N13_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N13_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 13 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 13) ^ 2) = (1 - x) ^ 4 * peval [12, 59, 150, 281, 448, 647, 874, 1125, 1396, 1683, 1982, 2289, 2444, 2289, 1982, 1683, 1396, 1125, 874, 647, 448, 281, 150, 59, 12] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 13 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (13 - 1 : ℕ) = 12 by norm_num, show (2 * 13 - 1 : ℕ) = 25 by norm_num, show (13 + 1 : ℕ) = 14 by norm_num, show (13 - 1 + 13 : ℕ) = 25 by norm_num, show (1 + 13 : ℕ) = 14 by norm_num, show (13 + 13 : ℕ) = 26 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N13_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 59, 150, 281, 448, 647, 874, 1125, 1396, 1683, 1982, 2289, 2444, 2289, 1982, 1683, 1396, 1125, 874, 647, 448, 281, 150, 59, 12] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 12 * q ^ 24 + 347 * p ^ 1 * q ^ 23 + 4819 * p ^ 2 * q ^ 22 + 42796 * p ^ 3 * q ^ 21 + 273000 * p ^ 4 * q ^ 20 + 1332110 * p ^ 5 * q ^ 19 + 5169710 * p ^ 6 * q ^ 18 + 16379220 * p ^ 7 * q ^ 17 + 43151680 * p ^ 8 * q ^ 16 + 95792125 * p ^ 9 * q ^ 15 + 180911445 * p ^ 10 * q ^ 14 + 292685120 * p ^ 11 * q ^ 13 + 407554784 * p ^ 12 * q ^ 12 + 489865766 * p ^ 13 * q ^ 11 + 508868126 * p ^ 14 * q ^ 10 + 456652196 * p ^ 15 * q ^ 9 + 353217436 * p ^ 16 * q ^ 8 + 234429195 * p ^ 17 * q ^ 7 + 132481635 * p ^ 18 * q ^ 6 + 62960612 * p ^ 19 * q ^ 5 + 24665212 * p ^ 20 * q ^ 4 + 7711132 * p ^ 21 * q ^ 3 + 1820468 * p ^ 22 * q ^ 2 + 292032 * p ^ 23 * q ^ 1 + 24336 * p ^ 24 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [12, 59, 150, 281, 448, 647, 874, 1125, 1396, 1683, 1982, 2289, 2444, 2289, 1982, 1683, 1396, 1125, 874, 647, 448, 281, 150, 59, 12] x =
      (fun p q : ℝ => 12 * q ^ 24 + 347 * p ^ 1 * q ^ 23 + 4819 * p ^ 2 * q ^ 22 + 42796 * p ^ 3 * q ^ 21 + 273000 * p ^ 4 * q ^ 20 + 1332110 * p ^ 5 * q ^ 19 + 5169710 * p ^ 6 * q ^ 18 + 16379220 * p ^ 7 * q ^ 17 + 43151680 * p ^ 8 * q ^ 16 + 95792125 * p ^ 9 * q ^ 15 + 180911445 * p ^ 10 * q ^ 14 + 292685120 * p ^ 11 * q ^ 13 + 407554784 * p ^ 12 * q ^ 12 + 489865766 * p ^ 13 * q ^ 11 + 508868126 * p ^ 14 * q ^ 10 + 456652196 * p ^ 15 * q ^ 9 + 353217436 * p ^ 16 * q ^ 8 + 234429195 * p ^ 17 * q ^ 7 + 132481635 * p ^ 18 * q ^ 6 + 62960612 * p ^ 19 * q ^ 5 + 24665212 * p ^ 20 * q ^ 4 + 7711132 * p ^ 21 * q ^ 3 + 1820468 * p ^ 22 * q ^ 2 + 292032 * p ^ 23 * q ^ 1 + 24336 * p ^ 24) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [12, 59, 150, 281, 448, 647, 874, 1125, 1396, 1683, 1982, 2289, 2444, 2289, 1982, 1683, 1396, 1125, 874, 647, 448, 281, 150, 59, 12] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N13_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 59, 150, 281, 448, 647, 874, 1125, 1396, 1683, 1982, 2289, 2444, 2289, 1982, 1683, 1396, 1125, 874, 647, 448, 281, 150, 59, 12] x := N13_posA1_0 x hu hv

theorem N13_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 13 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N13_posA1 x h0 h1)
  rw [← N13_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N13_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 13 2 x * (2 * (1 - x ^ 13) ^ 2) = (1 - x) ^ 4 * peval [22, 88, 194, 336, 510, 712, 938, 1184, 1446, 1720, 2002, 2145, 2002, 1720, 1446, 1184, 938, 712, 510, 336, 194, 88, 22] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 13 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (13 - 2 : ℕ) = 11 by norm_num, show (2 * 13 - 2 : ℕ) = 24 by norm_num, show (13 + 2 : ℕ) = 15 by norm_num, show (13 - 2 + 13 : ℕ) = 24 by norm_num, show (2 + 13 : ℕ) = 15 by norm_num, show (13 + 13 : ℕ) = 26 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N13_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [22, 88, 194, 336, 510, 712, 938, 1184, 1446, 1720, 2002, 2145, 2002, 1720, 1446, 1184, 938, 712, 510, 336, 194, 88, 22] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 22 * q ^ 22 + 572 * p ^ 1 * q ^ 21 + 7124 * p ^ 2 * q ^ 20 + 56576 * p ^ 3 * q ^ 19 + 321724 * p ^ 4 * q ^ 18 + 1394536 * p ^ 5 * q ^ 17 + 4788784 * p ^ 6 * q ^ 16 + 13366496 * p ^ 7 * q ^ 15 + 30870554 * p ^ 8 * q ^ 14 + 59741396 * p ^ 9 * q ^ 13 + 97735924 * p ^ 10 * q ^ 12 + 135980273 * p ^ 11 * q ^ 11 + 161474599 * p ^ 12 * q ^ 10 + 163898735 * p ^ 13 * q ^ 9 + 142100101 * p ^ 14 * q ^ 8 + 104903370 * p ^ 15 * q ^ 7 + 65518596 * p ^ 16 * q ^ 6 + 34231626 * p ^ 17 * q ^ 5 + 14682382 * p ^ 18 * q ^ 4 + 5010005 * p ^ 19 * q ^ 3 + 1288287 * p ^ 20 * q ^ 2 + 224939 * p ^ 21 * q ^ 1 + 20449 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [22, 88, 194, 336, 510, 712, 938, 1184, 1446, 1720, 2002, 2145, 2002, 1720, 1446, 1184, 938, 712, 510, 336, 194, 88, 22] x =
      (fun p q : ℝ => 22 * q ^ 22 + 572 * p ^ 1 * q ^ 21 + 7124 * p ^ 2 * q ^ 20 + 56576 * p ^ 3 * q ^ 19 + 321724 * p ^ 4 * q ^ 18 + 1394536 * p ^ 5 * q ^ 17 + 4788784 * p ^ 6 * q ^ 16 + 13366496 * p ^ 7 * q ^ 15 + 30870554 * p ^ 8 * q ^ 14 + 59741396 * p ^ 9 * q ^ 13 + 97735924 * p ^ 10 * q ^ 12 + 135980273 * p ^ 11 * q ^ 11 + 161474599 * p ^ 12 * q ^ 10 + 163898735 * p ^ 13 * q ^ 9 + 142100101 * p ^ 14 * q ^ 8 + 104903370 * p ^ 15 * q ^ 7 + 65518596 * p ^ 16 * q ^ 6 + 34231626 * p ^ 17 * q ^ 5 + 14682382 * p ^ 18 * q ^ 4 + 5010005 * p ^ 19 * q ^ 3 + 1288287 * p ^ 20 * q ^ 2 + 224939 * p ^ 21 * q ^ 1 + 20449 * p ^ 22) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [22, 88, 194, 336, 510, 712, 938, 1184, 1446, 1720, 2002, 2145, 2002, 1720, 1446, 1184, 938, 712, 510, 336, 194, 88, 22] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N13_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [22, 88, 194, 336, 510, 712, 938, 1184, 1446, 1720, 2002, 2145, 2002, 1720, 1446, 1184, 938, 712, 510, 336, 194, 88, 22] x := N13_posB2_0 x hu hv

theorem N13_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 13 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N13_posB2 x h0 h1)
  rw [← N13_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N13_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 13 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 13) ^ 2) = (1 - x) ^ 4 * peval [2, 52, 120, 228, 372, 548, 752, 980, 1228, 1492, 1768, 1909, 2054, 1909, 1768, 1492, 1228, 980, 752, 548, 372, 228, 120, 52, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 13 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (13 - 2 : ℕ) = 11 by norm_num, show (2 * 13 - 2 : ℕ) = 24 by norm_num, show (13 + 2 : ℕ) = 15 by norm_num, show (13 - 2 + 13 : ℕ) = 24 by norm_num, show (2 + 13 : ℕ) = 15 by norm_num, show (13 + 13 : ℕ) = 26 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N13_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 52, 120, 228, 372, 548, 752, 980, 1228, 1492, 1768, 1909, 2054, 1909, 1768, 1492, 1228, 980, 752, 548, 372, 228, 120, 52, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 24 + 100 * p ^ 1 * q ^ 23 + 1868 * p ^ 2 * q ^ 22 + 20072 * p ^ 3 * q ^ 21 + 146224 * p ^ 4 * q ^ 20 + 786136 * p ^ 5 * q ^ 19 + 3281824 * p ^ 6 * q ^ 18 + 10998416 * p ^ 7 * q ^ 17 + 30278534 * p ^ 8 * q ^ 16 + 69608396 * p ^ 9 * q ^ 15 + 135230524 * p ^ 10 * q ^ 14 + 223922153 * p ^ 11 * q ^ 13 + 317929339 * p ^ 12 * q ^ 12 + 388551111 * p ^ 13 * q ^ 11 + 409536517 * p ^ 14 * q ^ 10 + 372301930 * p ^ 15 * q ^ 9 + 291342480 * p ^ 16 * q ^ 8 + 195384618 * p ^ 17 * q ^ 7 + 111424742 * p ^ 18 * q ^ 6 + 53358877 * p ^ 19 * q ^ 5 + 21030867 * p ^ 20 * q ^ 4 + 6605027 * p ^ 21 * q ^ 3 + 1564433 * p ^ 22 * q ^ 2 + 251472 * p ^ 23 * q ^ 1 + 20956 * p ^ 24 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 52, 120, 228, 372, 548, 752, 980, 1228, 1492, 1768, 1909, 2054, 1909, 1768, 1492, 1228, 980, 752, 548, 372, 228, 120, 52, 2] x =
      (fun p q : ℝ => 2 * q ^ 24 + 100 * p ^ 1 * q ^ 23 + 1868 * p ^ 2 * q ^ 22 + 20072 * p ^ 3 * q ^ 21 + 146224 * p ^ 4 * q ^ 20 + 786136 * p ^ 5 * q ^ 19 + 3281824 * p ^ 6 * q ^ 18 + 10998416 * p ^ 7 * q ^ 17 + 30278534 * p ^ 8 * q ^ 16 + 69608396 * p ^ 9 * q ^ 15 + 135230524 * p ^ 10 * q ^ 14 + 223922153 * p ^ 11 * q ^ 13 + 317929339 * p ^ 12 * q ^ 12 + 388551111 * p ^ 13 * q ^ 11 + 409536517 * p ^ 14 * q ^ 10 + 372301930 * p ^ 15 * q ^ 9 + 291342480 * p ^ 16 * q ^ 8 + 195384618 * p ^ 17 * q ^ 7 + 111424742 * p ^ 18 * q ^ 6 + 53358877 * p ^ 19 * q ^ 5 + 21030867 * p ^ 20 * q ^ 4 + 6605027 * p ^ 21 * q ^ 3 + 1564433 * p ^ 22 * q ^ 2 + 251472 * p ^ 23 * q ^ 1 + 20956 * p ^ 24) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 52, 120, 228, 372, 548, 752, 980, 1228, 1492, 1768, 1909, 2054, 1909, 1768, 1492, 1228, 980, 752, 548, 372, 228, 120, 52, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N13_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 52, 120, 228, 372, 548, 752, 980, 1228, 1492, 1768, 1909, 2054, 1909, 1768, 1492, 1228, 980, 752, 548, 372, 228, 120, 52, 2] x := N13_posA2_0 x hu hv

theorem N13_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 13 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N13_posA2 x h0 h1)
  rw [← N13_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N13_linkB3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 13 3 x * (2 * (1 - x ^ 13) ^ 2) = (1 - x) ^ 4 * peval [30, 120, 300, 561, 894, 1290, 1740, 2235, 2766, 3324, 3770, 3965, 3770, 3324, 2766, 2235, 1740, 1290, 894, 561, 300, 120, 30] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 13 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (13 - 3 : ℕ) = 10 by norm_num, show (2 * 13 - 3 : ℕ) = 23 by norm_num, show (13 + 3 : ℕ) = 16 by norm_num, show (13 - 3 + 13 : ℕ) = 23 by norm_num, show (3 + 13 : ℕ) = 16 by norm_num, show (13 + 13 : ℕ) = 26 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N13_posB3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [30, 120, 300, 561, 894, 1290, 1740, 2235, 2766, 3324, 3770, 3965, 3770, 3324, 2766, 2235, 1740, 1290, 894, 561, 300, 120, 30] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 30 * q ^ 22 + 780 * p ^ 1 * q ^ 21 + 9750 * p ^ 2 * q ^ 20 + 77961 * p ^ 3 * q ^ 19 + 447603 * p ^ 4 * q ^ 18 + 1963533 * p ^ 5 * q ^ 17 + 6837831 * p ^ 6 * q ^ 16 + 19388655 * p ^ 7 * q ^ 15 + 45555939 * p ^ 8 * q ^ 14 + 89799567 * p ^ 9 * q ^ 13 + 149790797 * p ^ 10 * q ^ 12 + 212660474 * p ^ 11 * q ^ 11 + 257844990 * p ^ 12 * q ^ 10 + 267326904 * p ^ 13 * q ^ 9 + 236766972 * p ^ 14 * q ^ 8 + 178507602 * p ^ 15 * q ^ 7 + 113771814 * p ^ 16 * q ^ 6 + 60576360 * p ^ 17 * q ^ 5 + 26424840 * p ^ 18 * q ^ 4 + 9147125 * p ^ 19 * q ^ 3 + 2378675 * p ^ 20 * q ^ 2 + 418275 * p ^ 21 * q ^ 1 + 38025 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [30, 120, 300, 561, 894, 1290, 1740, 2235, 2766, 3324, 3770, 3965, 3770, 3324, 2766, 2235, 1740, 1290, 894, 561, 300, 120, 30] x =
      (fun p q : ℝ => 30 * q ^ 22 + 780 * p ^ 1 * q ^ 21 + 9750 * p ^ 2 * q ^ 20 + 77961 * p ^ 3 * q ^ 19 + 447603 * p ^ 4 * q ^ 18 + 1963533 * p ^ 5 * q ^ 17 + 6837831 * p ^ 6 * q ^ 16 + 19388655 * p ^ 7 * q ^ 15 + 45555939 * p ^ 8 * q ^ 14 + 89799567 * p ^ 9 * q ^ 13 + 149790797 * p ^ 10 * q ^ 12 + 212660474 * p ^ 11 * q ^ 11 + 257844990 * p ^ 12 * q ^ 10 + 267326904 * p ^ 13 * q ^ 9 + 236766972 * p ^ 14 * q ^ 8 + 178507602 * p ^ 15 * q ^ 7 + 113771814 * p ^ 16 * q ^ 6 + 60576360 * p ^ 17 * q ^ 5 + 26424840 * p ^ 18 * q ^ 4 + 9147125 * p ^ 19 * q ^ 3 + 2378675 * p ^ 20 * q ^ 2 + 418275 * p ^ 21 * q ^ 1 + 38025 * p ^ 22) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [30, 120, 300, 561, 894, 1290, 1740, 2235, 2766, 3324, 3770, 3965, 3770, 3324, 2766, 2235, 1740, 1290, 894, 561, 300, 120, 30] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N13_posB3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [30, 120, 300, 561, 894, 1290, 1740, 2235, 2766, 3324, 3770, 3965, 3770, 3324, 2766, 2235, 1740, 1290, 894, 561, 300, 120, 30] x := N13_posB3_0 x hu hv

theorem N13_B3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 13 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N13_posB3 x h0 h1)
  rw [← N13_linkB3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N13_linkA3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 13 3 x * (2 * (1 - x) ^ 2 * (1 - x ^ 13) ^ 2) = (1 - x) ^ 4 * peval [-6, 36, 122, 209, 332, 487, 670, 877, 1104, 1347, 1472, 1605, 1742, 1605, 1472, 1347, 1104, 877, 670, 487, 332, 209, 122, 36, -6] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 13 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (13 - 3 : ℕ) = 10 by norm_num, show (2 * 13 - 3 : ℕ) = 23 by norm_num, show (13 + 3 : ℕ) = 16 by norm_num, show (13 - 3 + 13 : ℕ) = 23 by norm_num, show (3 + 13 : ℕ) = 16 by norm_num, show (13 + 13 : ℕ) = 26 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N13_posA3_0 (x : ℝ) (hu : ((17 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-6, 36, 122, 209, 332, 487, 670, 877, 1104, 1347, 1472, 1605, 1742, 1605, 1472, 1347, 1104, 877, 670, 487, 332, 209, 122, 36, -6] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 417941727343338095952394679429950683433932794 * q ^ 24 + 13155353159325192540557918936984360175774799104 * p ^ 1 * q ^ 23 + 194589377641702221798597183012406830022850240512 * p ^ 2 * q ^ 22 + 1811449261770880820859203352582916864205593509888 * p ^ 3 * q ^ 21 + 11967237188918141566695064020778724942662965657600 * p ^ 4 * q ^ 20 + 59909298032498763243775650247877102473620957954048 * p ^ 5 * q ^ 19 + 236754413502338509444393160676711020906854750355456 * p ^ 6 * q ^ 18 + 759232189467163217522460938249583769593779218022400 * p ^ 7 * q ^ 17 + 2014439371508493476826085268996431369440833061257216 * p ^ 8 * q ^ 16 + 4484479805786523319390545156047805939264693920071680 * p ^ 9 * q ^ 15 + 8461309458075397442696088380555655328333140076789760 * p ^ 10 * q ^ 14 + 13627882364836184944328781249017954488157783967399936 * p ^ 11 * q ^ 13 + 18824414503817040110719643892983651901851520647299072 * p ^ 12 * q ^ 12 + 22357142407524349543117449747533613835938395659960320 * p ^ 13 * q ^ 11 + 22841034416507147108603121340735546880854681999900672 * p ^ 14 * q ^ 10 + 20039826443845822780371229898234537623783094272131072 * p ^ 15 * q ^ 9 + 15037783056984390115550348182931038820725907009830912 * p ^ 16 * q ^ 8 + 9584327443236202411162954610305408475046001476894720 * p ^ 17 * q ^ 7 + 5133060070314982690192805661773271167080955193065472 * p ^ 18 * q ^ 6 + 2273610901702301121167851893976136379184669802889216 * p ^ 19 * q ^ 5 + 813354555947886645813817813783852527597197253935104 * p ^ 20 * q ^ 4 + 226567529861114475355727466910812579405523750748160 * p ^ 21 * q ^ 3 + 46265681204470534714097018801273763891137365409792 * p ^ 22 * q ^ 2 + 6181816745709385428735769560544239268271159771136 * p ^ 23 * q ^ 1 + 407033201363580933579309929912377894207154552832 * p ^ 24 := by
    intro p q hp hq; positivity
  have key : (13500460747057082764996435506735298654081 : ℝ) * peval [-6, 36, 122, 209, 332, 487, 670, 877, 1104, 1347, 1472, 1605, 1742, 1605, 1472, 1347, 1104, 877, 670, 487, 332, 209, 122, 36, -6] x =
      (fun p q : ℝ => 417941727343338095952394679429950683433932794 * q ^ 24 + 13155353159325192540557918936984360175774799104 * p ^ 1 * q ^ 23 + 194589377641702221798597183012406830022850240512 * p ^ 2 * q ^ 22 + 1811449261770880820859203352582916864205593509888 * p ^ 3 * q ^ 21 + 11967237188918141566695064020778724942662965657600 * p ^ 4 * q ^ 20 + 59909298032498763243775650247877102473620957954048 * p ^ 5 * q ^ 19 + 236754413502338509444393160676711020906854750355456 * p ^ 6 * q ^ 18 + 759232189467163217522460938249583769593779218022400 * p ^ 7 * q ^ 17 + 2014439371508493476826085268996431369440833061257216 * p ^ 8 * q ^ 16 + 4484479805786523319390545156047805939264693920071680 * p ^ 9 * q ^ 15 + 8461309458075397442696088380555655328333140076789760 * p ^ 10 * q ^ 14 + 13627882364836184944328781249017954488157783967399936 * p ^ 11 * q ^ 13 + 18824414503817040110719643892983651901851520647299072 * p ^ 12 * q ^ 12 + 22357142407524349543117449747533613835938395659960320 * p ^ 13 * q ^ 11 + 22841034416507147108603121340735546880854681999900672 * p ^ 14 * q ^ 10 + 20039826443845822780371229898234537623783094272131072 * p ^ 15 * q ^ 9 + 15037783056984390115550348182931038820725907009830912 * p ^ 16 * q ^ 8 + 9584327443236202411162954610305408475046001476894720 * p ^ 17 * q ^ 7 + 5133060070314982690192805661773271167080955193065472 * p ^ 18 * q ^ 6 + 2273610901702301121167851893976136379184669802889216 * p ^ 19 * q ^ 5 + 813354555947886645813817813783852527597197253935104 * p ^ 20 * q ^ 4 + 226567529861114475355727466910812579405523750748160 * p ^ 21 * q ^ 3 + 46265681204470534714097018801273763891137365409792 * p ^ 22 * q ^ 2 + 6181816745709385428735769560544239268271159771136 * p ^ 23 * q ^ 1 + 407033201363580933579309929912377894207154552832 * p ^ 24) (x - ((17 : ℝ) / 64)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (13500460747057082764996435506735298654081 : ℝ) * peval [-6, 36, 122, 209, 332, 487, 670, 877, 1104, 1347, 1472, 1605, 1742, 1605, 1472, 1347, 1104, 877, 670, 487, 332, 209, 122, 36, -6] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N13_posA3 (x : ℝ) (hu : ((17 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-6, 36, 122, 209, 332, 487, 670, 877, 1104, 1347, 1472, 1605, 1742, 1605, 1472, 1347, 1104, 877, 670, 487, 332, 209, 122, 36, -6] x := N13_posA3_0 x hu hv

theorem N13_A3_pos (x : ℝ) (h0 : 0 < x) (hlo : ((17 : ℝ) / 64) < x) (h1 : x < 1) :
    0 < aBlock 13 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N13_posA3 x hlo h1)
  rw [← N13_linkA3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N13_linkB4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 13 4 x * (2 * (1 - x ^ 13) ^ 2) = (1 - x) ^ 4 * peval [36, 144, 360, 720, 1208, 1808, 2504, 3280, 4120, 4891, 5460, 5694, 5460, 4891, 4120, 3280, 2504, 1808, 1208, 720, 360, 144, 36] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 13 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (13 - 4 : ℕ) = 9 by norm_num, show (2 * 13 - 4 : ℕ) = 22 by norm_num, show (13 + 4 : ℕ) = 17 by norm_num, show (13 - 4 + 13 : ℕ) = 22 by norm_num, show (4 + 13 : ℕ) = 17 by norm_num, show (13 + 13 : ℕ) = 26 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N13_posB4_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [36, 144, 360, 720, 1208, 1808, 2504, 3280, 4120, 4891, 5460, 5694, 5460, 4891, 4120, 3280, 2504, 1808, 1208, 720, 360, 144, 36] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 36 * q ^ 22 + 936 * p ^ 1 * q ^ 21 + 11700 * p ^ 2 * q ^ 20 + 93600 * p ^ 3 * q ^ 19 + 538148 * p ^ 4 * q ^ 18 + 2366936 * p ^ 5 * q ^ 17 + 8276268 * p ^ 6 * q ^ 16 + 23600720 * p ^ 7 * q ^ 15 + 55861520 * p ^ 8 * q ^ 14 + 111114315 * p ^ 9 * q ^ 13 + 187340595 * p ^ 10 * q ^ 12 + 269257560 * p ^ 11 * q ^ 11 + 330969600 * p ^ 12 * q ^ 10 + 348277566 * p ^ 13 * q ^ 9 + 313331070 * p ^ 14 * q ^ 8 + 240032052 * p ^ 15 * q ^ 7 + 155397528 * p ^ 16 * q ^ 6 + 83963763 * p ^ 17 * q ^ 5 + 37110879 * p ^ 18 * q ^ 4 + 12989340 * p ^ 19 * q ^ 3 + 3407040 * p ^ 20 * q ^ 2 + 602316 * p ^ 21 * q ^ 1 + 54756 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [36, 144, 360, 720, 1208, 1808, 2504, 3280, 4120, 4891, 5460, 5694, 5460, 4891, 4120, 3280, 2504, 1808, 1208, 720, 360, 144, 36] x =
      (fun p q : ℝ => 36 * q ^ 22 + 936 * p ^ 1 * q ^ 21 + 11700 * p ^ 2 * q ^ 20 + 93600 * p ^ 3 * q ^ 19 + 538148 * p ^ 4 * q ^ 18 + 2366936 * p ^ 5 * q ^ 17 + 8276268 * p ^ 6 * q ^ 16 + 23600720 * p ^ 7 * q ^ 15 + 55861520 * p ^ 8 * q ^ 14 + 111114315 * p ^ 9 * q ^ 13 + 187340595 * p ^ 10 * q ^ 12 + 269257560 * p ^ 11 * q ^ 11 + 330969600 * p ^ 12 * q ^ 10 + 348277566 * p ^ 13 * q ^ 9 + 313331070 * p ^ 14 * q ^ 8 + 240032052 * p ^ 15 * q ^ 7 + 155397528 * p ^ 16 * q ^ 6 + 83963763 * p ^ 17 * q ^ 5 + 37110879 * p ^ 18 * q ^ 4 + 12989340 * p ^ 19 * q ^ 3 + 3407040 * p ^ 20 * q ^ 2 + 602316 * p ^ 21 * q ^ 1 + 54756 * p ^ 22) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [36, 144, 360, 720, 1208, 1808, 2504, 3280, 4120, 4891, 5460, 5694, 5460, 4891, 4120, 3280, 2504, 1808, 1208, 720, 360, 144, 36] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N13_posB4 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [36, 144, 360, 720, 1208, 1808, 2504, 3280, 4120, 4891, 5460, 5694, 5460, 4891, 4120, 3280, 2504, 1808, 1208, 720, 360, 144, 36] x := N13_posB4_0 x hu hv

theorem N13_B4_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 13 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N13_posB4 x h0 h1)
  rw [← N13_linkB4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N13_linkA4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 13 4 x * (2 * (1 - x) ^ 2 * (1 - x ^ 13) ^ 2) = (1 - x) ^ 4 * peval [-12, 24, 104, 224, 328, 464, 628, 816, 1024, 1131, 1250, 1377, 1508, 1377, 1250, 1131, 1024, 816, 628, 464, 328, 224, 104, 24, -12] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 13 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (13 - 4 : ℕ) = 9 by norm_num, show (2 * 13 - 4 : ℕ) = 22 by norm_num, show (13 + 4 : ℕ) = 17 by norm_num, show (13 - 4 + 13 : ℕ) = 22 by norm_num, show (4 + 13 : ℕ) = 17 by norm_num, show (13 + 13 : ℕ) = 26 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N13_posA4_0 (x : ℝ) (hu : ((17 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-12, 24, 104, 224, 328, 464, 628, 816, 1024, 1131, 1250, 1377, 1508, 1377, 1250, 1131, 1024, 816, 628, 464, 328, 224, 104, 24, -12] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 189334480217393087211752291118838758942319604 * q ^ 24 + 7343337863292954691290321392756503449008758272 * p ^ 1 * q ^ 23 + 123782975918028097390699044205197909101393870848 * p ^ 2 * q ^ 22 + 1261022666431270800526056106502394406640101097472 * p ^ 3 * q ^ 21 + 8899291928954210499784433612124918224217496879104 * p ^ 4 * q ^ 20 + 46844221964339708929778686964012791457568459325440 * p ^ 5 * q ^ 19 + 192506076017913085292845207434328470448435944751104 * p ^ 6 * q ^ 18 + 636698017510994067317128378313944690980895177310208 * p ^ 7 * q ^ 17 + 1731208822445231177118599660299017335778745513934848 * p ^ 8 * q ^ 16 + 3929149584456332187114015961096083434806905925533696 * p ^ 9 * q ^ 15 + 7525541025031527297174213132855239841907008327385088 * p ^ 10 * q ^ 14 + 12258266543338548900572288215328310577172170030448640 * p ^ 11 * q ^ 13 + 17069052973281708988175084031218734782039920742498304 * p ^ 12 * q ^ 12 + 20376980676436042463972766999882092456907590414106624 * p ^ 13 * q ^ 11 + 20872084638313457191235683990119999426399340758827008 * p ^ 14 * q ^ 10 + 18318920135024794744089265452165632451035293059907584 * p ^ 15 * q ^ 9 + 13725265475540194702464808236868186214441209657556992 * p ^ 16 * q ^ 8 + 8721012244732171401864503853904837286332053288124416 * p ^ 17 * q ^ 7 + 4651270383231307894463278070680635048148046830370816 * p ^ 18 * q ^ 6 + 2050369850022238400820243416988279736959499147149312 * p ^ 19 * q ^ 5 + 729927296382035908943040288819439339478070239166464 * p ^ 20 * q ^ 4 + 202437202802666313953132815056564627443032378048512 * p ^ 21 * q ^ 3 + 41200028575995904339313578804055938866141118595072 * p ^ 22 * q ^ 2 + 5494948218408342603320684053817101571796586463232 * p ^ 23 * q ^ 1 + 361807290100960829848275493255447017073026269184 * p ^ 24 := by
    intro p q hp hq; positivity
  have key : (13500460747057082764996435506735298654081 : ℝ) * peval [-12, 24, 104, 224, 328, 464, 628, 816, 1024, 1131, 1250, 1377, 1508, 1377, 1250, 1131, 1024, 816, 628, 464, 328, 224, 104, 24, -12] x =
      (fun p q : ℝ => 189334480217393087211752291118838758942319604 * q ^ 24 + 7343337863292954691290321392756503449008758272 * p ^ 1 * q ^ 23 + 123782975918028097390699044205197909101393870848 * p ^ 2 * q ^ 22 + 1261022666431270800526056106502394406640101097472 * p ^ 3 * q ^ 21 + 8899291928954210499784433612124918224217496879104 * p ^ 4 * q ^ 20 + 46844221964339708929778686964012791457568459325440 * p ^ 5 * q ^ 19 + 192506076017913085292845207434328470448435944751104 * p ^ 6 * q ^ 18 + 636698017510994067317128378313944690980895177310208 * p ^ 7 * q ^ 17 + 1731208822445231177118599660299017335778745513934848 * p ^ 8 * q ^ 16 + 3929149584456332187114015961096083434806905925533696 * p ^ 9 * q ^ 15 + 7525541025031527297174213132855239841907008327385088 * p ^ 10 * q ^ 14 + 12258266543338548900572288215328310577172170030448640 * p ^ 11 * q ^ 13 + 17069052973281708988175084031218734782039920742498304 * p ^ 12 * q ^ 12 + 20376980676436042463972766999882092456907590414106624 * p ^ 13 * q ^ 11 + 20872084638313457191235683990119999426399340758827008 * p ^ 14 * q ^ 10 + 18318920135024794744089265452165632451035293059907584 * p ^ 15 * q ^ 9 + 13725265475540194702464808236868186214441209657556992 * p ^ 16 * q ^ 8 + 8721012244732171401864503853904837286332053288124416 * p ^ 17 * q ^ 7 + 4651270383231307894463278070680635048148046830370816 * p ^ 18 * q ^ 6 + 2050369850022238400820243416988279736959499147149312 * p ^ 19 * q ^ 5 + 729927296382035908943040288819439339478070239166464 * p ^ 20 * q ^ 4 + 202437202802666313953132815056564627443032378048512 * p ^ 21 * q ^ 3 + 41200028575995904339313578804055938866141118595072 * p ^ 22 * q ^ 2 + 5494948218408342603320684053817101571796586463232 * p ^ 23 * q ^ 1 + 361807290100960829848275493255447017073026269184 * p ^ 24) (x - ((17 : ℝ) / 64)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (13500460747057082764996435506735298654081 : ℝ) * peval [-12, 24, 104, 224, 328, 464, 628, 816, 1024, 1131, 1250, 1377, 1508, 1377, 1250, 1131, 1024, 816, 628, 464, 328, 224, 104, 24, -12] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N13_posA4 (x : ℝ) (hu : ((17 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-12, 24, 104, 224, 328, 464, 628, 816, 1024, 1131, 1250, 1377, 1508, 1377, 1250, 1131, 1024, 816, 628, 464, 328, 224, 104, 24, -12] x := N13_posA4_0 x hu hv

theorem N13_A4_pos (x : ℝ) (h0 : 0 < x) (hlo : ((17 : ℝ) / 64) < x) (h1 : x < 1) :
    0 < aBlock 13 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N13_posA4 x hlo h1)
  rw [← N13_linkA4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N13_linkB5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 13 5 x * (2 * (1 - x ^ 13) ^ 2) = (1 - x) ^ 4 * peval [40, 160, 400, 800, 1400, 2175, 3100, 4150, 5196, 6109, 6760, 7020, 6760, 6109, 5196, 4150, 3100, 2175, 1400, 800, 400, 160, 40] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 13 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (13 - 5 : ℕ) = 8 by norm_num, show (2 * 13 - 5 : ℕ) = 21 by norm_num, show (13 + 5 : ℕ) = 18 by norm_num, show (13 - 5 + 13 : ℕ) = 21 by norm_num, show (5 + 13 : ℕ) = 18 by norm_num, show (13 + 13 : ℕ) = 26 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N13_posB5_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [40, 160, 400, 800, 1400, 2175, 3100, 4150, 5196, 6109, 6760, 7020, 6760, 6109, 5196, 4150, 3100, 2175, 1400, 800, 400, 160, 40] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 40 * q ^ 22 + 1040 * p ^ 1 * q ^ 21 + 13000 * p ^ 2 * q ^ 20 + 104000 * p ^ 3 * q ^ 19 + 598000 * p ^ 4 * q ^ 18 + 2631135 * p ^ 5 * q ^ 17 + 9207835 * p ^ 6 * q ^ 16 + 26298350 * p ^ 7 * q ^ 15 + 62404446 * p ^ 8 * q ^ 14 + 124591103 * p ^ 9 * q ^ 13 + 211130803 * p ^ 10 * q ^ 12 + 305434376 * p ^ 11 * q ^ 11 + 378431560 * p ^ 12 * q ^ 10 + 401906336 * p ^ 13 * q ^ 9 + 365275248 * p ^ 14 * q ^ 8 + 282830288 * p ^ 15 * q ^ 7 + 185065816 * p ^ 16 * q ^ 6 + 101004540 * p ^ 17 * q ^ 5 + 45045260 * p ^ 18 * q ^ 4 + 15886000 * p ^ 19 * q ^ 3 + 4191200 * p ^ 20 * q ^ 2 + 743600 * p ^ 21 * q ^ 1 + 67600 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [40, 160, 400, 800, 1400, 2175, 3100, 4150, 5196, 6109, 6760, 7020, 6760, 6109, 5196, 4150, 3100, 2175, 1400, 800, 400, 160, 40] x =
      (fun p q : ℝ => 40 * q ^ 22 + 1040 * p ^ 1 * q ^ 21 + 13000 * p ^ 2 * q ^ 20 + 104000 * p ^ 3 * q ^ 19 + 598000 * p ^ 4 * q ^ 18 + 2631135 * p ^ 5 * q ^ 17 + 9207835 * p ^ 6 * q ^ 16 + 26298350 * p ^ 7 * q ^ 15 + 62404446 * p ^ 8 * q ^ 14 + 124591103 * p ^ 9 * q ^ 13 + 211130803 * p ^ 10 * q ^ 12 + 305434376 * p ^ 11 * q ^ 11 + 378431560 * p ^ 12 * q ^ 10 + 401906336 * p ^ 13 * q ^ 9 + 365275248 * p ^ 14 * q ^ 8 + 282830288 * p ^ 15 * q ^ 7 + 185065816 * p ^ 16 * q ^ 6 + 101004540 * p ^ 17 * q ^ 5 + 45045260 * p ^ 18 * q ^ 4 + 15886000 * p ^ 19 * q ^ 3 + 4191200 * p ^ 20 * q ^ 2 + 743600 * p ^ 21 * q ^ 1 + 67600 * p ^ 22) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [40, 160, 400, 800, 1400, 2175, 3100, 4150, 5196, 6109, 6760, 7020, 6760, 6109, 5196, 4150, 3100, 2175, 1400, 800, 400, 160, 40] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N13_posB5 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [40, 160, 400, 800, 1400, 2175, 3100, 4150, 5196, 6109, 6760, 7020, 6760, 6109, 5196, 4150, 3100, 2175, 1400, 800, 400, 160, 40] x := N13_posB5_0 x hu hv

theorem N13_B5_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 13 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N13_posB5 x h0 h1)
  rw [← N13_linkB5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N13_linkA5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 13 5 x * (2 * (1 - x) ^ 2 * (1 - x ^ 13) ^ 2) = (1 - x) ^ 4 * peval [-16, 16, 92, 208, 360, 479, 626, 797, 884, 987, 1102, 1225, 1352, 1225, 1102, 987, 884, 797, 626, 479, 360, 208, 92, 16, -16] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 13 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (13 - 5 : ℕ) = 8 by norm_num, show (2 * 13 - 5 : ℕ) = 21 by norm_num, show (13 + 5 : ℕ) = 18 by norm_num, show (13 - 5 + 13 : ℕ) = 21 by norm_num, show (5 + 13 : ℕ) = 18 by norm_num, show (13 + 13 : ℕ) = 26 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N13_posA5_0 (x : ℝ) (hu : ((17 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-16, 16, 92, 208, 360, 479, 626, 797, 884, 987, 1102, 1225, 1352, 1225, 1102, 987, 884, 797, 626, 479, 360, 208, 92, 16, -16] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 31006890920922213026088250843083669448305648 * q ^ 24 + 3294502779388019845693320187988240643135759360 * p ^ 1 * q ^ 23 + 74223778814007813393883449133462438380764708864 * p ^ 2 * q ^ 22 + 874684963048205686035907625558403916541491937280 * p ^ 3 * q ^ 21 + 6745926163205825394081864607731197519544465752064 * p ^ 4 * q ^ 20 + 37707886944356887715723847632538305803418282754048 * p ^ 5 * q ^ 19 + 161819514234085442779303835603060747534255670689792 * p ^ 6 * q ^ 18 + 552871308919339306346250161471477396878657132167168 * p ^ 7 * q ^ 17 + 1541162318908712909217366574005728691478231466901504 * p ^ 8 * q ^ 16 + 3565642866842441294876113108352943761471618087387136 * p ^ 9 * q ^ 15 + 6930437145511981417113178729706989510155499575705600 * p ^ 10 * q ^ 14 + 11413254589656745443019393910925064892003742614814720 * p ^ 11 * q ^ 13 + 16015742029564371248207099074713087248784645144182784 * p ^ 12 * q ^ 12 + 19213251401259047318491538372778205146977603904602112 * p ^ 13 * q ^ 11 + 19726660670859379665839983353255949775030903547887616 * p ^ 14 * q ^ 10 + 17315937435854352166765474432437111149064171055742976 * p ^ 15 * q ^ 9 + 12950651626282883920250023254243322254046645750595584 * p ^ 16 * q ^ 8 + 8201269447097344984169252083027384707966271204884480 * p ^ 17 * q ^ 7 + 4354374628659328688513188035223709162545832664236032 * p ^ 18 * q ^ 6 + 1909555501880399760933156406828814545265888304562176 * p ^ 19 * q ^ 5 + 676195717131193747431502332607873880976492628803584 * p ^ 20 * q ^ 4 + 186633287564958333542137251979714015243432749957120 * p ^ 21 * q ^ 3 + 37843252352239690756356929814735421138618655178752 * p ^ 22 * q ^ 2 + 5037035866874314053043960382665676440813537591296 * p ^ 23 * q ^ 1 + 331656682592547427360919202150826432316940746752 * p ^ 24 := by
    intro p q hp hq; positivity
  have key : (13500460747057082764996435506735298654081 : ℝ) * peval [-16, 16, 92, 208, 360, 479, 626, 797, 884, 987, 1102, 1225, 1352, 1225, 1102, 987, 884, 797, 626, 479, 360, 208, 92, 16, -16] x =
      (fun p q : ℝ => 31006890920922213026088250843083669448305648 * q ^ 24 + 3294502779388019845693320187988240643135759360 * p ^ 1 * q ^ 23 + 74223778814007813393883449133462438380764708864 * p ^ 2 * q ^ 22 + 874684963048205686035907625558403916541491937280 * p ^ 3 * q ^ 21 + 6745926163205825394081864607731197519544465752064 * p ^ 4 * q ^ 20 + 37707886944356887715723847632538305803418282754048 * p ^ 5 * q ^ 19 + 161819514234085442779303835603060747534255670689792 * p ^ 6 * q ^ 18 + 552871308919339306346250161471477396878657132167168 * p ^ 7 * q ^ 17 + 1541162318908712909217366574005728691478231466901504 * p ^ 8 * q ^ 16 + 3565642866842441294876113108352943761471618087387136 * p ^ 9 * q ^ 15 + 6930437145511981417113178729706989510155499575705600 * p ^ 10 * q ^ 14 + 11413254589656745443019393910925064892003742614814720 * p ^ 11 * q ^ 13 + 16015742029564371248207099074713087248784645144182784 * p ^ 12 * q ^ 12 + 19213251401259047318491538372778205146977603904602112 * p ^ 13 * q ^ 11 + 19726660670859379665839983353255949775030903547887616 * p ^ 14 * q ^ 10 + 17315937435854352166765474432437111149064171055742976 * p ^ 15 * q ^ 9 + 12950651626282883920250023254243322254046645750595584 * p ^ 16 * q ^ 8 + 8201269447097344984169252083027384707966271204884480 * p ^ 17 * q ^ 7 + 4354374628659328688513188035223709162545832664236032 * p ^ 18 * q ^ 6 + 1909555501880399760933156406828814545265888304562176 * p ^ 19 * q ^ 5 + 676195717131193747431502332607873880976492628803584 * p ^ 20 * q ^ 4 + 186633287564958333542137251979714015243432749957120 * p ^ 21 * q ^ 3 + 37843252352239690756356929814735421138618655178752 * p ^ 22 * q ^ 2 + 5037035866874314053043960382665676440813537591296 * p ^ 23 * q ^ 1 + 331656682592547427360919202150826432316940746752 * p ^ 24) (x - ((17 : ℝ) / 64)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (13500460747057082764996435506735298654081 : ℝ) * peval [-16, 16, 92, 208, 360, 479, 626, 797, 884, 987, 1102, 1225, 1352, 1225, 1102, 987, 884, 797, 626, 479, 360, 208, 92, 16, -16] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N13_posA5 (x : ℝ) (hu : ((17 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-16, 16, 92, 208, 360, 479, 626, 797, 884, 987, 1102, 1225, 1352, 1225, 1102, 987, 884, 797, 626, 479, 360, 208, 92, 16, -16] x := N13_posA5_0 x hu hv

theorem N13_A5_pos (x : ℝ) (h0 : 0 < x) (hlo : ((17 : ℝ) / 64) < x) (h1 : x < 1) :
    0 < aBlock 13 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N13_posA5 x hlo h1)
  rw [← N13_linkA5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N13_linkB6 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 13 6 x * (2 * (1 - x ^ 13) ^ 2) = (1 - x) ^ 4 * peval [42, 168, 420, 840, 1470, 2352, 3450, 4637, 5786, 6770, 7462, 7735, 7462, 6770, 5786, 4637, 3450, 2352, 1470, 840, 420, 168, 42] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 13 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (13 - 6 : ℕ) = 7 by norm_num, show (2 * 13 - 6 : ℕ) = 20 by norm_num, show (13 + 6 : ℕ) = 19 by norm_num, show (13 - 6 + 13 : ℕ) = 20 by norm_num, show (6 + 13 : ℕ) = 19 by norm_num, show (13 + 13 : ℕ) = 26 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N13_posB6_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [42, 168, 420, 840, 1470, 2352, 3450, 4637, 5786, 6770, 7462, 7735, 7462, 6770, 5786, 4637, 3450, 2352, 1470, 840, 420, 168, 42] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 42 * q ^ 22 + 1092 * p ^ 1 * q ^ 21 + 13650 * p ^ 2 * q ^ 20 + 109200 * p ^ 3 * q ^ 19 + 627900 * p ^ 4 * q ^ 18 + 2762760 * p ^ 5 * q ^ 17 + 9669582 * p ^ 6 * q ^ 16 + 27625949 * p ^ 7 * q ^ 15 + 65599001 * p ^ 8 * q ^ 14 + 131126619 * p ^ 9 * q ^ 13 + 222626781 * p ^ 10 * q ^ 12 + 322936692 * p ^ 11 * q ^ 11 + 401541972 * p ^ 12 * q ^ 10 + 428309882 * p ^ 13 * q ^ 9 + 391215734 * p ^ 14 * q ^ 8 + 304539690 * p ^ 15 * q ^ 7 + 200348148 * p ^ 16 * q ^ 6 + 109905432 * p ^ 17 * q ^ 5 + 49238826 * p ^ 18 * q ^ 4 + 17431505 * p ^ 19 * q ^ 3 + 4612517 * p ^ 20 * q ^ 2 + 819819 * p ^ 21 * q ^ 1 + 74529 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [42, 168, 420, 840, 1470, 2352, 3450, 4637, 5786, 6770, 7462, 7735, 7462, 6770, 5786, 4637, 3450, 2352, 1470, 840, 420, 168, 42] x =
      (fun p q : ℝ => 42 * q ^ 22 + 1092 * p ^ 1 * q ^ 21 + 13650 * p ^ 2 * q ^ 20 + 109200 * p ^ 3 * q ^ 19 + 627900 * p ^ 4 * q ^ 18 + 2762760 * p ^ 5 * q ^ 17 + 9669582 * p ^ 6 * q ^ 16 + 27625949 * p ^ 7 * q ^ 15 + 65599001 * p ^ 8 * q ^ 14 + 131126619 * p ^ 9 * q ^ 13 + 222626781 * p ^ 10 * q ^ 12 + 322936692 * p ^ 11 * q ^ 11 + 401541972 * p ^ 12 * q ^ 10 + 428309882 * p ^ 13 * q ^ 9 + 391215734 * p ^ 14 * q ^ 8 + 304539690 * p ^ 15 * q ^ 7 + 200348148 * p ^ 16 * q ^ 6 + 109905432 * p ^ 17 * q ^ 5 + 49238826 * p ^ 18 * q ^ 4 + 17431505 * p ^ 19 * q ^ 3 + 4612517 * p ^ 20 * q ^ 2 + 819819 * p ^ 21 * q ^ 1 + 74529 * p ^ 22) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [42, 168, 420, 840, 1470, 2352, 3450, 4637, 5786, 6770, 7462, 7735, 7462, 6770, 5786, 4637, 3450, 2352, 1470, 840, 420, 168, 42] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N13_posB6 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [42, 168, 420, 840, 1470, 2352, 3450, 4637, 5786, 6770, 7462, 7735, 7462, 6770, 5786, 4637, 3450, 2352, 1470, 840, 420, 168, 42] x := N13_posB6_0 x hu hv

theorem N13_B6_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 13 6 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N13_posB6 x h0 h1)
  rw [← N13_linkB6 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N13_linkA6 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 13 6 x * (2 * (1 - x) ^ 2 * (1 - x ^ 13) ^ 2) = (1 - x) ^ 4 * peval [-18, 12, 86, 200, 350, 532, 664, 729, 814, 915, 1028, 1149, 1274, 1149, 1028, 915, 814, 729, 664, 532, 350, 200, 86, 12, -18] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 13 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (13 - 6 : ℕ) = 7 by norm_num, show (2 * 13 - 6 : ℕ) = 20 by norm_num, show (13 + 6 : ℕ) = 19 by norm_num, show (13 - 6 + 13 : ℕ) = 20 by norm_num, show (6 + 13 : ℕ) = 19 by norm_num, show (13 + 13 : ℕ) = 26 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

/-- The minimal polynomial `f_13` of `ξ*(13)` (constant term first). -/
def f13 : List ℤ := [18, -12, -86, -200, -350, -532, -664, -729, -814, -915, -1028, -1149, -1274, -1149, -1028, -915, -814, -729, -664, -532, -350, -200, -86, -12, 18]

theorem N13_crit_factor (x : ℝ) : peval [-18, 12, 86, 200, 350, 532, 664, 729, 814, 915, 1028, 1149, 1274, 1149, 1028, 915, 814, 729, 664, 532, 350, 200, 86, 12, -18] x = -(peval [1] x * peval f13 x) := by
  simp only [f13, peval]; ring

theorem N13_posH (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1] x := by norm_num [peval]

theorem N13_dneg_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 172, 600, 1400, 2660, 3984, 5103, 6512, 8235, 10280, 12639, 15288, 14937, 14392, 13725, 13024, 12393, 11952, 10108, 7000, 4200, 1892, 276, -432] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 12 * q ^ 23 + 448 * p ^ 1 * q ^ 22 + 7420 * p ^ 2 * q ^ 21 + 74984 * p ^ 3 * q ^ 20 + 527800 * p ^ 4 * q ^ 19 + 2780492 * p ^ 5 * q ^ 18 + 11459487 * p ^ 6 * q ^ 17 + 38048075 * p ^ 7 * q ^ 16 + 103848875 * p ^ 8 * q ^ 15 + 236345525 * p ^ 9 * q ^ 14 + 453026158 * p ^ 10 * q ^ 13 + 736446308 * p ^ 11 * q ^ 12 + 1019933850 * p ^ 12 * q ^ 11 + 1206510942 * p ^ 13 * q ^ 10 + 1220019450 * p ^ 14 * q ^ 9 + 1053525382 * p ^ 15 * q ^ 8 + 774467512 * p ^ 16 * q ^ 7 + 481741260 * p ^ 17 * q ^ 6 + 250948607 * p ^ 18 * q ^ 5 + 107635255 * p ^ 19 * q ^ 4 + 36972299 * p ^ 20 * q ^ 3 + 9694685 * p ^ 21 * q ^ 2 + 1762670 * p ^ 22 * q ^ 1 + 170352 * p ^ 23 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [12, 172, 600, 1400, 2660, 3984, 5103, 6512, 8235, 10280, 12639, 15288, 14937, 14392, 13725, 13024, 12393, 11952, 10108, 7000, 4200, 1892, 276, -432] x =
      (fun p q : ℝ => 12 * q ^ 23 + 448 * p ^ 1 * q ^ 22 + 7420 * p ^ 2 * q ^ 21 + 74984 * p ^ 3 * q ^ 20 + 527800 * p ^ 4 * q ^ 19 + 2780492 * p ^ 5 * q ^ 18 + 11459487 * p ^ 6 * q ^ 17 + 38048075 * p ^ 7 * q ^ 16 + 103848875 * p ^ 8 * q ^ 15 + 236345525 * p ^ 9 * q ^ 14 + 453026158 * p ^ 10 * q ^ 13 + 736446308 * p ^ 11 * q ^ 12 + 1019933850 * p ^ 12 * q ^ 11 + 1206510942 * p ^ 13 * q ^ 10 + 1220019450 * p ^ 14 * q ^ 9 + 1053525382 * p ^ 15 * q ^ 8 + 774467512 * p ^ 16 * q ^ 7 + 481741260 * p ^ 17 * q ^ 6 + 250948607 * p ^ 18 * q ^ 5 + 107635255 * p ^ 19 * q ^ 4 + 36972299 * p ^ 20 * q ^ 3 + 9694685 * p ^ 21 * q ^ 2 + 1762670 * p ^ 22 * q ^ 1 + 170352 * p ^ 23) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [12, 172, 600, 1400, 2660, 3984, 5103, 6512, 8235, 10280, 12639, 15288, 14937, 14392, 13725, 13024, 12393, 11952, 10108, 7000, 4200, 1892, 276, -432] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N13_dneg (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [12, 172, 600, 1400, 2660, 3984, 5103, 6512, 8235, 10280, 12639, 15288, 14937, 14392, 13725, 13024, 12393, 11952, 10108, 7000, 4200, 1892, 276, -432] x := N13_dneg_0 x hu hv

theorem N13_dlink (x : ℝ) : dpeval f13 x = -peval [12, 172, 600, 1400, 2660, 3984, 5103, 6512, 8235, 10280, 12639, 15288, 14937, 14392, 13725, 13024, 12393, 11952, 10108, 7000, 4200, 1892, 276, -432] x := by
  simp only [f13, dpeval, peval]; ring

/-- `f_13` changes sign exactly once in `(0,1)`, at a point of `(13993069432195334351/50000000000000000000, 27986138864390668703/100000000000000000000)`. -/
theorem N13_root : ∃ ξs, ((13993069432195334351 : ℝ) / 50000000000000000000) < ξs ∧ ξs < ((27986138864390668703 : ℝ) / 100000000000000000000) ∧ peval f13 ξs = 0 ∧
    ∀ x, 0 < x → x < 1 → (x < ξs → 0 < peval f13 x) ∧ (ξs < x → peval f13 x < 0) := by
  apply sign_change (peval f13) (0 : ℝ) _ _ (by norm_num) (by norm_num) (by norm_num)
    (continuous_peval _).continuousOn
  · apply peval_strictAntiOn
    intro x hu hv
    rw [N13_dlink]
    linarith [N13_dneg x hu hv]
  · exact (fun x hx0 hx => absurd hx (not_le.mpr hx0))
  · norm_num [f13, peval]
  · norm_num [f13, peval]

theorem N13_A6_neg (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : 0 < peval f13 x) :
    aBlock 13 6 x < 0 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : (1 - x) ^ 4 * peval [-18, 12, 86, 200, 350, 532, 664, 729, 814, 915, 1028, 1149, 1274, 1149, 1028, 915, 814, 729, 664, 532, 350, 200, 86, 12, -18] x < 0 := by
    rw [N13_crit_factor]
    have := mul_pos (N13_posH x h0 h1) hf
    nlinarith [pow_pos h1x 4]
  rw [← N13_linkA6 x h0 h1] at hn
  have hd := (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))
  by_contra hc
  push Not at hc
  nlinarith [mul_nonneg hc hd.le]

theorem N13_A6_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : peval f13 x < 0) :
    0 < aBlock 13 6 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 13 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : 0 < (1 - x) ^ 4 * peval [-18, 12, 86, 200, 350, 532, 664, 729, 814, 915, 1028, 1149, 1274, 1149, 1028, 915, 814, 729, 664, 532, 350, 200, 86, 12, -18] x := by
    rw [N13_crit_factor]
    have := mul_pos (N13_posH x h0 h1) (neg_pos.mpr hf)
    have h4 := pow_pos h1x 4
    nlinarith
  rw [← N13_linkA6 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

/-- Note Thm 6(b) for `N = 13`: there is a unique `ξ* ∈ (0,1)` (the unique root of `f_13` in
`(0,1)`, enclosed in `(13993069432195334351/50000000000000000000, 27986138864390668703/100000000000000000000)`) such that, for every `ξ ∈ (0,1)`: all `b_m > 0`;
`a_6 < 0` for `ξ < ξ*` (spectrally unstable, since `b_6 > 0`); and all `a_m > 0` for `ξ > ξ*`. -/
theorem cert_N13 : ∃ ξs : ℝ, ((13993069432195334351 : ℝ) / 50000000000000000000) < ξs ∧ ξs < ((27986138864390668703 : ℝ) / 100000000000000000000) ∧ peval f13 ξs = 0 ∧
    (∀ y, 0 < y → y < 1 → peval f13 y = 0 → y = ξs) ∧
    ∀ x, 0 < x → x < 1 →
      (∀ m, 1 ≤ m → m ≤ 12 → 0 < bBlock 13 m x) ∧
      (x < ξs → aBlock 13 6 x < 0) ∧
      (ξs < x → ∀ m, 1 ≤ m → m ≤ 12 → 0 < aBlock 13 m x) := by
  obtain ⟨ξs, hta, htb, hf0, hs⟩ := N13_root
  refine ⟨ξs, hta, htb, hf0, ?_, ?_⟩
  · intro y hy0 hy1 hy
    rcases lt_trichotomy y ξs with h | h | h
    · linarith [(hs y hy0 hy1).1 h]
    · exact h
    · linarith [(hs y hy0 hy1).2 h]
  intro x h0 h1
  refine ⟨?_, fun hlt => N13_A6_neg x h0 h1 ((hs x h0 h1).1 hlt), ?_⟩
  · intro m hm1 hm2
    interval_cases m
    · exact N13_B1_pos x h0 h1
    · exact N13_B2_pos x h0 h1
    · exact N13_B3_pos x h0 h1
    · exact N13_B4_pos x h0 h1
    · exact N13_B5_pos x h0 h1
    · exact N13_B6_pos x h0 h1
    · rw [show bBlock 13 7 x = bBlock 13 6 x from bBlock_symm 13 6 (by norm_num) x]
      exact N13_B6_pos x h0 h1
    · rw [show bBlock 13 8 x = bBlock 13 5 x from bBlock_symm 13 5 (by norm_num) x]
      exact N13_B5_pos x h0 h1
    · rw [show bBlock 13 9 x = bBlock 13 4 x from bBlock_symm 13 4 (by norm_num) x]
      exact N13_B4_pos x h0 h1
    · rw [show bBlock 13 10 x = bBlock 13 3 x from bBlock_symm 13 3 (by norm_num) x]
      exact N13_B3_pos x h0 h1
    · rw [show bBlock 13 11 x = bBlock 13 2 x from bBlock_symm 13 2 (by norm_num) x]
      exact N13_B2_pos x h0 h1
    · rw [show bBlock 13 12 x = bBlock 13 1 x from bBlock_symm 13 1 (by norm_num) x]
      exact N13_B1_pos x h0 h1
  · intro hgt m hm1 hm2
    have hlo : ((17 : ℝ) / 64) < x := lt_trans (lt_trans (by norm_num) hta) hgt
    have hfn := (hs x h0 h1).2 hgt
    interval_cases m
    · exact N13_A1_pos x h0 h1
    · exact N13_A2_pos x h0 h1
    · exact N13_A3_pos x h0 hlo h1
    · exact N13_A4_pos x h0 hlo h1
    · exact N13_A5_pos x h0 hlo h1
    · exact N13_A6_pos x h0 h1 hfn
    · rw [show aBlock 13 7 x = aBlock 13 6 x from aBlock_symm 13 6 (by norm_num) x]
      exact N13_A6_pos x h0 h1 hfn
    · rw [show aBlock 13 8 x = aBlock 13 5 x from aBlock_symm 13 5 (by norm_num) x]
      exact N13_A5_pos x h0 hlo h1
    · rw [show aBlock 13 9 x = aBlock 13 4 x from aBlock_symm 13 4 (by norm_num) x]
      exact N13_A4_pos x h0 hlo h1
    · rw [show aBlock 13 10 x = aBlock 13 3 x from aBlock_symm 13 3 (by norm_num) x]
      exact N13_A3_pos x h0 hlo h1
    · rw [show aBlock 13 11 x = aBlock 13 2 x from aBlock_symm 13 2 (by norm_num) x]
      exact N13_A2_pos x h0 h1
    · rw [show aBlock 13 12 x = aBlock 13 1 x from aBlock_symm 13 1 (by norm_num) x]
      exact N13_A1_pos x h0 h1

/-! ### N = 14 -/

theorem N14_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 14 1 x * (2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [13, 38, 74, 120, 175, 238, 308, 384, 465, 550, 638, 728, 819, 728, 638, 550, 465, 384, 308, 238, 175, 120, 74, 38, 13] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (14 - 1 : ℕ) = 13 by norm_num, show (2 * 14 - 1 : ℕ) = 27 by norm_num, show (14 + 1 : ℕ) = 15 by norm_num, show (14 - 1 + 14 : ℕ) = 27 by norm_num, show (1 + 14 : ℕ) = 15 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [13, 38, 74, 120, 175, 238, 308, 384, 465, 550, 638, 728, 819, 728, 638, 550, 465, 384, 308, 238, 175, 120, 74, 38, 13] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 13 * q ^ 24 + 350 * p ^ 1 * q ^ 23 + 4536 * p ^ 2 * q ^ 22 + 37674 * p ^ 3 * q ^ 21 + 225225 * p ^ 4 * q ^ 20 + 1031940 * p ^ 5 * q ^ 19 + 3767400 * p ^ 6 * q ^ 18 + 11248380 * p ^ 7 * q ^ 17 + 27972945 * p ^ 8 * q ^ 16 + 58708650 * p ^ 9 * q ^ 15 + 104984880 * p ^ 10 * q ^ 14 + 161056350 * p ^ 11 * q ^ 13 + 212952285 * p ^ 12 * q ^ 12 + 243373858 * p ^ 13 * q ^ 11 + 240697236 * p ^ 14 * q ^ 10 + 205917712 * p ^ 15 * q ^ 9 + 152056905 * p ^ 16 * q ^ 8 + 96504408 * p ^ 17 * q ^ 7 + 52260208 * p ^ 18 * q ^ 6 + 23865842 * p ^ 19 * q ^ 5 + 9018009 * p ^ 20 * q ^ 4 + 2732730 * p ^ 21 * q ^ 3 + 629356 * p ^ 22 * q ^ 2 + 99372 * p ^ 23 * q ^ 1 + 8281 * p ^ 24 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [13, 38, 74, 120, 175, 238, 308, 384, 465, 550, 638, 728, 819, 728, 638, 550, 465, 384, 308, 238, 175, 120, 74, 38, 13] x =
      (fun p q : ℝ => 13 * q ^ 24 + 350 * p ^ 1 * q ^ 23 + 4536 * p ^ 2 * q ^ 22 + 37674 * p ^ 3 * q ^ 21 + 225225 * p ^ 4 * q ^ 20 + 1031940 * p ^ 5 * q ^ 19 + 3767400 * p ^ 6 * q ^ 18 + 11248380 * p ^ 7 * q ^ 17 + 27972945 * p ^ 8 * q ^ 16 + 58708650 * p ^ 9 * q ^ 15 + 104984880 * p ^ 10 * q ^ 14 + 161056350 * p ^ 11 * q ^ 13 + 212952285 * p ^ 12 * q ^ 12 + 243373858 * p ^ 13 * q ^ 11 + 240697236 * p ^ 14 * q ^ 10 + 205917712 * p ^ 15 * q ^ 9 + 152056905 * p ^ 16 * q ^ 8 + 96504408 * p ^ 17 * q ^ 7 + 52260208 * p ^ 18 * q ^ 6 + 23865842 * p ^ 19 * q ^ 5 + 9018009 * p ^ 20 * q ^ 4 + 2732730 * p ^ 21 * q ^ 3 + 629356 * p ^ 22 * q ^ 2 + 99372 * p ^ 23 * q ^ 1 + 8281 * p ^ 24) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [13, 38, 74, 120, 175, 238, 308, 384, 465, 550, 638, 728, 819, 728, 638, 550, 465, 384, 308, 238, 175, 120, 74, 38, 13] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [13, 38, 74, 120, 175, 238, 308, 384, 465, 550, 638, 728, 819, 728, 638, 550, 465, 384, 308, 238, 175, 120, 74, 38, 13] x := N14_posB1_0 x hu hv

theorem N14_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 14 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posB1 x h0 h1)
  rw [← N14_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N14_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 14 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [13, 64, 163, 306, 489, 708, 959, 1238, 1541, 1864, 2203, 2554, 2913, 3094, 2913, 2554, 2203, 1864, 1541, 1238, 959, 708, 489, 306, 163, 64, 13] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (14 - 1 : ℕ) = 13 by norm_num, show (2 * 14 - 1 : ℕ) = 27 by norm_num, show (14 + 1 : ℕ) = 15 by norm_num, show (14 - 1 + 14 : ℕ) = 27 by norm_num, show (1 + 14 : ℕ) = 15 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [13, 64, 163, 306, 489, 708, 959, 1238, 1541, 1864, 2203, 2554, 2913, 3094, 2913, 2554, 2203, 1864, 1541, 1238, 959, 708, 489, 306, 163, 64, 13] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 13 * q ^ 26 + 402 * p ^ 1 * q ^ 25 + 5988 * p ^ 2 * q ^ 24 + 57218 * p ^ 3 * q ^ 23 + 394065 * p ^ 4 * q ^ 22 + 2083536 * p ^ 5 * q ^ 21 + 8796060 * p ^ 6 * q ^ 20 + 30445740 * p ^ 7 * q ^ 19 + 88036065 * p ^ 8 * q ^ 18 + 215593950 * p ^ 9 * q ^ 17 + 451711260 * p ^ 10 * q ^ 16 + 815830470 * p ^ 11 * q ^ 15 + 1277117205 * p ^ 12 * q ^ 14 + 1739408398 * p ^ 13 * q ^ 13 + 2066001808 * p ^ 14 * q ^ 12 + 2142202088 * p ^ 15 * q ^ 11 + 1938516697 * p ^ 16 * q ^ 10 + 1528402876 * p ^ 17 * q ^ 9 + 1046505460 * p ^ 18 * q ^ 8 + 618924306 * p ^ 19 * q ^ 7 + 313522209 * p ^ 20 * q ^ 6 + 134268134 * p ^ 21 * q ^ 5 + 47632312 * p ^ 22 * q ^ 4 + 13547716 * p ^ 23 * q ^ 3 + 2923193 * p ^ 24 * q ^ 2 + 430612 * p ^ 25 * q ^ 1 + 33124 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [13, 64, 163, 306, 489, 708, 959, 1238, 1541, 1864, 2203, 2554, 2913, 3094, 2913, 2554, 2203, 1864, 1541, 1238, 959, 708, 489, 306, 163, 64, 13] x =
      (fun p q : ℝ => 13 * q ^ 26 + 402 * p ^ 1 * q ^ 25 + 5988 * p ^ 2 * q ^ 24 + 57218 * p ^ 3 * q ^ 23 + 394065 * p ^ 4 * q ^ 22 + 2083536 * p ^ 5 * q ^ 21 + 8796060 * p ^ 6 * q ^ 20 + 30445740 * p ^ 7 * q ^ 19 + 88036065 * p ^ 8 * q ^ 18 + 215593950 * p ^ 9 * q ^ 17 + 451711260 * p ^ 10 * q ^ 16 + 815830470 * p ^ 11 * q ^ 15 + 1277117205 * p ^ 12 * q ^ 14 + 1739408398 * p ^ 13 * q ^ 13 + 2066001808 * p ^ 14 * q ^ 12 + 2142202088 * p ^ 15 * q ^ 11 + 1938516697 * p ^ 16 * q ^ 10 + 1528402876 * p ^ 17 * q ^ 9 + 1046505460 * p ^ 18 * q ^ 8 + 618924306 * p ^ 19 * q ^ 7 + 313522209 * p ^ 20 * q ^ 6 + 134268134 * p ^ 21 * q ^ 5 + 47632312 * p ^ 22 * q ^ 4 + 13547716 * p ^ 23 * q ^ 3 + 2923193 * p ^ 24 * q ^ 2 + 430612 * p ^ 25 * q ^ 1 + 33124 * p ^ 26) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [13, 64, 163, 306, 489, 708, 959, 1238, 1541, 1864, 2203, 2554, 2913, 3094, 2913, 2554, 2203, 1864, 1541, 1238, 959, 708, 489, 306, 163, 64, 13] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [13, 64, 163, 306, 489, 708, 959, 1238, 1541, 1864, 2203, 2554, 2913, 3094, 2913, 2554, 2203, 1864, 1541, 1238, 959, 708, 489, 306, 163, 64, 13] x := N14_posA1_0 x hu hv

theorem N14_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 14 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posA1 x h0 h1)
  rw [← N14_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N14_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 14 2 x * (2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [24, 96, 212, 368, 560, 784, 1036, 1312, 1608, 1920, 2244, 2576, 2744, 2576, 2244, 1920, 1608, 1312, 1036, 784, 560, 368, 212, 96, 24] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (14 - 2 : ℕ) = 12 by norm_num, show (2 * 14 - 2 : ℕ) = 26 by norm_num, show (14 + 2 : ℕ) = 16 by norm_num, show (14 - 2 + 14 : ℕ) = 26 by norm_num, show (2 + 14 : ℕ) = 16 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [24, 96, 212, 368, 560, 784, 1036, 1312, 1608, 1920, 2244, 2576, 2744, 2576, 2244, 1920, 1608, 1312, 1036, 784, 560, 368, 212, 96, 24] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 24 * q ^ 24 + 672 * p ^ 1 * q ^ 23 + 9044 * p ^ 2 * q ^ 22 + 77896 * p ^ 3 * q ^ 21 + 482300 * p ^ 4 * q ^ 20 + 2285920 * p ^ 5 * q ^ 19 + 8623160 * p ^ 6 * q ^ 18 + 26575120 * p ^ 7 * q ^ 17 + 68148080 * p ^ 8 * q ^ 16 + 147347200 * p ^ 9 * q ^ 15 + 271210940 * p ^ 10 * q ^ 14 + 427892920 * p ^ 11 * q ^ 13 + 581393372 * p ^ 12 * q ^ 12 + 682276672 * p ^ 13 * q ^ 11 + 692362984 * p ^ 14 * q ^ 10 + 607305776 * p ^ 15 * q ^ 9 + 459431896 * p ^ 16 * q ^ 8 + 298437216 * p ^ 17 * q ^ 7 + 165213972 * p ^ 18 * q ^ 6 + 77008008 * p ^ 19 * q ^ 5 + 29638532 * p ^ 20 * q ^ 4 + 9124192 * p ^ 21 * q ^ 3 + 2127776 * p ^ 22 * q ^ 2 + 338688 * p ^ 23 * q ^ 1 + 28224 * p ^ 24 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [24, 96, 212, 368, 560, 784, 1036, 1312, 1608, 1920, 2244, 2576, 2744, 2576, 2244, 1920, 1608, 1312, 1036, 784, 560, 368, 212, 96, 24] x =
      (fun p q : ℝ => 24 * q ^ 24 + 672 * p ^ 1 * q ^ 23 + 9044 * p ^ 2 * q ^ 22 + 77896 * p ^ 3 * q ^ 21 + 482300 * p ^ 4 * q ^ 20 + 2285920 * p ^ 5 * q ^ 19 + 8623160 * p ^ 6 * q ^ 18 + 26575120 * p ^ 7 * q ^ 17 + 68148080 * p ^ 8 * q ^ 16 + 147347200 * p ^ 9 * q ^ 15 + 271210940 * p ^ 10 * q ^ 14 + 427892920 * p ^ 11 * q ^ 13 + 581393372 * p ^ 12 * q ^ 12 + 682276672 * p ^ 13 * q ^ 11 + 692362984 * p ^ 14 * q ^ 10 + 607305776 * p ^ 15 * q ^ 9 + 459431896 * p ^ 16 * q ^ 8 + 298437216 * p ^ 17 * q ^ 7 + 165213972 * p ^ 18 * q ^ 6 + 77008008 * p ^ 19 * q ^ 5 + 29638532 * p ^ 20 * q ^ 4 + 9124192 * p ^ 21 * q ^ 3 + 2127776 * p ^ 22 * q ^ 2 + 338688 * p ^ 23 * q ^ 1 + 28224 * p ^ 24) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [24, 96, 212, 368, 560, 784, 1036, 1312, 1608, 1920, 2244, 2576, 2744, 2576, 2244, 1920, 1608, 1312, 1036, 784, 560, 368, 212, 96, 24] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [24, 96, 212, 368, 560, 784, 1036, 1312, 1608, 1920, 2244, 2576, 2744, 2576, 2244, 1920, 1608, 1312, 1036, 784, 560, 368, 212, 96, 24] x := N14_posB2_0 x hu hv

theorem N14_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 14 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posB2 x h0 h1)
  rw [← N14_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N14_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 14 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [2, 56, 130, 248, 406, 600, 826, 1080, 1358, 1656, 1970, 2296, 2462, 2632, 2462, 2296, 1970, 1656, 1358, 1080, 826, 600, 406, 248, 130, 56, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (14 - 2 : ℕ) = 12 by norm_num, show (2 * 14 - 2 : ℕ) = 26 by norm_num, show (14 + 2 : ℕ) = 16 by norm_num, show (14 - 2 + 14 : ℕ) = 26 by norm_num, show (2 + 14 : ℕ) = 16 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 56, 130, 248, 406, 600, 826, 1080, 1358, 1656, 1970, 2296, 2462, 2632, 2462, 2296, 1970, 1656, 1358, 1080, 826, 600, 406, 248, 130, 56, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 26 + 108 * p ^ 1 * q ^ 25 + 2180 * p ^ 2 * q ^ 24 + 25368 * p ^ 3 * q ^ 23 + 200690 * p ^ 4 * q ^ 22 + 1175356 * p ^ 5 * q ^ 21 + 5363540 * p ^ 6 * q ^ 20 + 19723600 * p ^ 7 * q ^ 19 + 59832890 * p ^ 8 * q ^ 18 + 152280700 * p ^ 9 * q ^ 17 + 329228900 * p ^ 10 * q ^ 16 + 610235080 * p ^ 11 * q ^ 15 + 976279682 * p ^ 12 * q ^ 14 + 1354583692 * p ^ 13 * q ^ 13 + 1635102400 * p ^ 14 * q ^ 12 + 1719863248 * p ^ 15 * q ^ 11 + 1576617082 * p ^ 16 * q ^ 10 + 1257919740 * p ^ 17 * q ^ 9 + 870794848 * p ^ 18 * q ^ 8 + 520202760 * p ^ 19 * q ^ 7 + 265896554 * p ^ 20 * q ^ 6 + 114761724 * p ^ 21 * q ^ 5 + 40974584 * p ^ 22 * q ^ 4 + 11712960 * p ^ 23 * q ^ 3 + 2536730 * p ^ 24 * q ^ 2 + 374556 * p ^ 25 * q ^ 1 + 28812 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 56, 130, 248, 406, 600, 826, 1080, 1358, 1656, 1970, 2296, 2462, 2632, 2462, 2296, 1970, 1656, 1358, 1080, 826, 600, 406, 248, 130, 56, 2] x =
      (fun p q : ℝ => 2 * q ^ 26 + 108 * p ^ 1 * q ^ 25 + 2180 * p ^ 2 * q ^ 24 + 25368 * p ^ 3 * q ^ 23 + 200690 * p ^ 4 * q ^ 22 + 1175356 * p ^ 5 * q ^ 21 + 5363540 * p ^ 6 * q ^ 20 + 19723600 * p ^ 7 * q ^ 19 + 59832890 * p ^ 8 * q ^ 18 + 152280700 * p ^ 9 * q ^ 17 + 329228900 * p ^ 10 * q ^ 16 + 610235080 * p ^ 11 * q ^ 15 + 976279682 * p ^ 12 * q ^ 14 + 1354583692 * p ^ 13 * q ^ 13 + 1635102400 * p ^ 14 * q ^ 12 + 1719863248 * p ^ 15 * q ^ 11 + 1576617082 * p ^ 16 * q ^ 10 + 1257919740 * p ^ 17 * q ^ 9 + 870794848 * p ^ 18 * q ^ 8 + 520202760 * p ^ 19 * q ^ 7 + 265896554 * p ^ 20 * q ^ 6 + 114761724 * p ^ 21 * q ^ 5 + 40974584 * p ^ 22 * q ^ 4 + 11712960 * p ^ 23 * q ^ 3 + 2536730 * p ^ 24 * q ^ 2 + 374556 * p ^ 25 * q ^ 1 + 28812 * p ^ 26) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 56, 130, 248, 406, 600, 826, 1080, 1358, 1656, 1970, 2296, 2462, 2632, 2462, 2296, 1970, 1656, 1358, 1080, 826, 600, 406, 248, 130, 56, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 56, 130, 248, 406, 600, 826, 1080, 1358, 1656, 1970, 2296, 2462, 2632, 2462, 2296, 1970, 1656, 1358, 1080, 826, 600, 406, 248, 130, 56, 2] x := N14_posA2_0 x hu hv

theorem N14_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 14 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posA2 x h0 h1)
  rw [← N14_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N14_linkB3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 14 3 x * (2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [33, 132, 330, 618, 987, 1428, 1932, 2490, 3093, 3732, 4398, 4928, 5159, 4928, 4398, 3732, 3093, 2490, 1932, 1428, 987, 618, 330, 132, 33] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (14 - 3 : ℕ) = 11 by norm_num, show (2 * 14 - 3 : ℕ) = 25 by norm_num, show (14 + 3 : ℕ) = 17 by norm_num, show (14 - 3 + 14 : ℕ) = 25 by norm_num, show (3 + 14 : ℕ) = 17 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posB3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [33, 132, 330, 618, 987, 1428, 1932, 2490, 3093, 3732, 4398, 4928, 5159, 4928, 4398, 3732, 3093, 2490, 1932, 1428, 987, 618, 330, 132, 33] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 33 * q ^ 24 + 924 * p ^ 1 * q ^ 23 + 12474 * p ^ 2 * q ^ 22 + 108066 * p ^ 3 * q ^ 21 + 674625 * p ^ 4 * q ^ 20 + 3230640 * p ^ 5 * q ^ 19 + 12335820 * p ^ 6 * q ^ 18 + 38542020 * p ^ 7 * q ^ 17 + 100336005 * p ^ 8 * q ^ 16 + 220489500 * p ^ 9 * q ^ 15 + 412873230 * p ^ 10 * q ^ 14 + 663221636 * p ^ 11 * q ^ 13 + 918110347 * p ^ 12 * q ^ 12 + 1098282416 * p ^ 13 * q ^ 11 + 1136532606 * p ^ 14 * q ^ 10 + 1016816724 * p ^ 15 * q ^ 9 + 784582029 * p ^ 16 * q ^ 8 + 519649746 * p ^ 17 * q ^ 7 + 293104350 * p ^ 18 * q ^ 6 + 139023192 * p ^ 19 * q ^ 5 + 54351143 * p ^ 20 * q ^ 4 + 16956940 * p ^ 21 * q ^ 3 + 3996146 * p ^ 22 * q ^ 2 + 640332 * p ^ 23 * q ^ 1 + 53361 * p ^ 24 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [33, 132, 330, 618, 987, 1428, 1932, 2490, 3093, 3732, 4398, 4928, 5159, 4928, 4398, 3732, 3093, 2490, 1932, 1428, 987, 618, 330, 132, 33] x =
      (fun p q : ℝ => 33 * q ^ 24 + 924 * p ^ 1 * q ^ 23 + 12474 * p ^ 2 * q ^ 22 + 108066 * p ^ 3 * q ^ 21 + 674625 * p ^ 4 * q ^ 20 + 3230640 * p ^ 5 * q ^ 19 + 12335820 * p ^ 6 * q ^ 18 + 38542020 * p ^ 7 * q ^ 17 + 100336005 * p ^ 8 * q ^ 16 + 220489500 * p ^ 9 * q ^ 15 + 412873230 * p ^ 10 * q ^ 14 + 663221636 * p ^ 11 * q ^ 13 + 918110347 * p ^ 12 * q ^ 12 + 1098282416 * p ^ 13 * q ^ 11 + 1136532606 * p ^ 14 * q ^ 10 + 1016816724 * p ^ 15 * q ^ 9 + 784582029 * p ^ 16 * q ^ 8 + 519649746 * p ^ 17 * q ^ 7 + 293104350 * p ^ 18 * q ^ 6 + 139023192 * p ^ 19 * q ^ 5 + 54351143 * p ^ 20 * q ^ 4 + 16956940 * p ^ 21 * q ^ 3 + 3996146 * p ^ 22 * q ^ 2 + 640332 * p ^ 23 * q ^ 1 + 53361 * p ^ 24) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [33, 132, 330, 618, 987, 1428, 1932, 2490, 3093, 3732, 4398, 4928, 5159, 4928, 4398, 3732, 3093, 2490, 1932, 1428, 987, 618, 330, 132, 33] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posB3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [33, 132, 330, 618, 987, 1428, 1932, 2490, 3093, 3732, 4398, 4928, 5159, 4928, 4398, 3732, 3093, 2490, 1932, 1428, 987, 618, 330, 132, 33] x := N14_posB3_0 x hu hv

theorem N14_B3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 14 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posB3 x h0 h1)
  rw [← N14_linkB3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N14_linkA3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 14 3 x * (2 * (1 - x) ^ 2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [-7, 38, 131, 226, 361, 532, 735, 966, 1221, 1496, 1787, 1936, 2093, 2254, 2093, 1936, 1787, 1496, 1221, 966, 735, 532, 361, 226, 131, 38, -7] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (14 - 3 : ℕ) = 11 by norm_num, show (2 * 14 - 3 : ℕ) = 25 by norm_num, show (14 + 3 : ℕ) = 17 by norm_num, show (14 - 3 + 14 : ℕ) = 25 by norm_num, show (3 + 14 : ℕ) = 17 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posA3_0 (x : ℝ) (hu : ((5 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-7, 38, 131, 226, 361, 532, 735, 966, 1221, 1496, 1787, 1936, 2093, 2254, 2093, 1936, 1787, 1496, 1221, 966, 735, 532, 361, 226, 131, 38, -7] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 623407089545721843881328086182353 * q ^ 26 + 20110612599147351031285578983821056 * p ^ 1 * q ^ 25 + 309190192967220147567677962630821376 * p ^ 2 * q ^ 24 + 3024071593530067335765388694888325120 * p ^ 3 * q ^ 23 + 21176018708849542321989766570355261440 * p ^ 4 * q ^ 22 + 113225678012237074393848334577037213696 * p ^ 5 * q ^ 21 + 481247917931881134580305958370415738880 * p ^ 6 * q ^ 20 + 1670839563909685238631546775447409786880 * p ^ 7 * q ^ 19 + 4830937847956399542947783448967564492800 * p ^ 8 * q ^ 18 + 11797330407753287761407371496393592012800 * p ^ 9 * q ^ 17 + 24587596643165623542928530822158289469440 * p ^ 10 * q ^ 16 + 44070805419155013697348027072050433622016 * p ^ 11 * q ^ 15 + 68305181550126150254473944580676147740672 * p ^ 12 * q ^ 14 + 91871132569857987351151106371064394940416 * p ^ 13 * q ^ 13 + 107436317255063829987451253103898392526848 * p ^ 14 * q ^ 12 + 109266683932839049705566113109143401267200 * p ^ 15 * q ^ 11 + 96511075510217525474373704953398722297856 * p ^ 16 * q ^ 10 + 73790531402563775701391144042610026348544 * p ^ 17 * q ^ 9 + 48572388809576114994221588930958513405952 * p ^ 18 * q ^ 8 + 27301525503124140354688088719619518365696 * p ^ 19 * q ^ 7 + 12949733541515866633164254427302863044608 * p ^ 20 * q ^ 6 + 5096318812259079709732052994714117341184 * p ^ 21 * q ^ 5 + 1623512802748844516645709948977332355072 * p ^ 22 * q ^ 4 + 403258319492965810951881272396240388096 * p ^ 23 * q ^ 3 + 73455258942342339305263074591603949568 * p ^ 24 * q ^ 2 + 8749998832894560747171841769610412032 * p ^ 25 * q ^ 1 + 512820444418728834999082301515628544 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (1191817653772720942460132761 : ℝ) * peval [-7, 38, 131, 226, 361, 532, 735, 966, 1221, 1496, 1787, 1936, 2093, 2254, 2093, 1936, 1787, 1496, 1221, 966, 735, 532, 361, 226, 131, 38, -7] x =
      (fun p q : ℝ => 623407089545721843881328086182353 * q ^ 26 + 20110612599147351031285578983821056 * p ^ 1 * q ^ 25 + 309190192967220147567677962630821376 * p ^ 2 * q ^ 24 + 3024071593530067335765388694888325120 * p ^ 3 * q ^ 23 + 21176018708849542321989766570355261440 * p ^ 4 * q ^ 22 + 113225678012237074393848334577037213696 * p ^ 5 * q ^ 21 + 481247917931881134580305958370415738880 * p ^ 6 * q ^ 20 + 1670839563909685238631546775447409786880 * p ^ 7 * q ^ 19 + 4830937847956399542947783448967564492800 * p ^ 8 * q ^ 18 + 11797330407753287761407371496393592012800 * p ^ 9 * q ^ 17 + 24587596643165623542928530822158289469440 * p ^ 10 * q ^ 16 + 44070805419155013697348027072050433622016 * p ^ 11 * q ^ 15 + 68305181550126150254473944580676147740672 * p ^ 12 * q ^ 14 + 91871132569857987351151106371064394940416 * p ^ 13 * q ^ 13 + 107436317255063829987451253103898392526848 * p ^ 14 * q ^ 12 + 109266683932839049705566113109143401267200 * p ^ 15 * q ^ 11 + 96511075510217525474373704953398722297856 * p ^ 16 * q ^ 10 + 73790531402563775701391144042610026348544 * p ^ 17 * q ^ 9 + 48572388809576114994221588930958513405952 * p ^ 18 * q ^ 8 + 27301525503124140354688088719619518365696 * p ^ 19 * q ^ 7 + 12949733541515866633164254427302863044608 * p ^ 20 * q ^ 6 + 5096318812259079709732052994714117341184 * p ^ 21 * q ^ 5 + 1623512802748844516645709948977332355072 * p ^ 22 * q ^ 4 + 403258319492965810951881272396240388096 * p ^ 23 * q ^ 3 + 73455258942342339305263074591603949568 * p ^ 24 * q ^ 2 + 8749998832894560747171841769610412032 * p ^ 25 * q ^ 1 + 512820444418728834999082301515628544 * p ^ 26) (x - ((5 : ℝ) / 16)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1191817653772720942460132761 : ℝ) * peval [-7, 38, 131, 226, 361, 532, 735, 966, 1221, 1496, 1787, 1936, 2093, 2254, 2093, 1936, 1787, 1496, 1221, 966, 735, 532, 361, 226, 131, 38, -7] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posA3 (x : ℝ) (hu : ((5 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-7, 38, 131, 226, 361, 532, 735, 966, 1221, 1496, 1787, 1936, 2093, 2254, 2093, 1936, 1787, 1496, 1221, 966, 735, 532, 361, 226, 131, 38, -7] x := N14_posA3_0 x hu hv

theorem N14_A3_pos (x : ℝ) (h0 : 0 < x) (hlo : ((5 : ℝ) / 16) < x) (h1 : x < 1) :
    0 < aBlock 14 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posA3 x hlo h1)
  rw [← N14_linkA3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N14_linkB4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 14 4 x * (2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [40, 160, 400, 800, 1344, 2016, 2800, 3680, 4640, 5664, 6596, 7280, 7560, 7280, 6596, 5664, 4640, 3680, 2800, 2016, 1344, 800, 400, 160, 40] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (14 - 4 : ℕ) = 10 by norm_num, show (2 * 14 - 4 : ℕ) = 24 by norm_num, show (14 + 4 : ℕ) = 18 by norm_num, show (14 - 4 + 14 : ℕ) = 24 by norm_num, show (4 + 14 : ℕ) = 18 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posB4_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [40, 160, 400, 800, 1344, 2016, 2800, 3680, 4640, 5664, 6596, 7280, 7560, 7280, 6596, 5664, 4640, 3680, 2800, 2016, 1344, 800, 400, 160, 40] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 40 * q ^ 24 + 1120 * p ^ 1 * q ^ 23 + 15120 * p ^ 2 * q ^ 22 + 131040 * p ^ 3 * q ^ 21 + 818944 * p ^ 4 * q ^ 20 + 3929856 * p ^ 5 * q ^ 19 + 15054144 * p ^ 6 * q ^ 18 + 47248256 * p ^ 7 * q ^ 17 + 123729144 * p ^ 8 * q ^ 16 + 273895776 * p ^ 9 * q ^ 15 + 517386884 * p ^ 10 * q ^ 14 + 839582856 * p ^ 11 * q ^ 13 + 1175662404 * p ^ 12 * q ^ 12 + 1424351936 * p ^ 13 * q ^ 11 + 1494405576 * p ^ 14 * q ^ 10 + 1356707184 * p ^ 15 * q ^ 9 + 1062866896 * p ^ 16 * q ^ 8 + 714841344 * p ^ 17 * q ^ 7 + 409280116 * p ^ 18 * q ^ 6 + 196872200 * p ^ 19 * q ^ 5 + 77942340 * p ^ 20 * q ^ 4 + 24578400 * p ^ 21 * q ^ 3 + 5840800 * p ^ 22 * q ^ 2 + 940800 * p ^ 23 * q ^ 1 + 78400 * p ^ 24 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [40, 160, 400, 800, 1344, 2016, 2800, 3680, 4640, 5664, 6596, 7280, 7560, 7280, 6596, 5664, 4640, 3680, 2800, 2016, 1344, 800, 400, 160, 40] x =
      (fun p q : ℝ => 40 * q ^ 24 + 1120 * p ^ 1 * q ^ 23 + 15120 * p ^ 2 * q ^ 22 + 131040 * p ^ 3 * q ^ 21 + 818944 * p ^ 4 * q ^ 20 + 3929856 * p ^ 5 * q ^ 19 + 15054144 * p ^ 6 * q ^ 18 + 47248256 * p ^ 7 * q ^ 17 + 123729144 * p ^ 8 * q ^ 16 + 273895776 * p ^ 9 * q ^ 15 + 517386884 * p ^ 10 * q ^ 14 + 839582856 * p ^ 11 * q ^ 13 + 1175662404 * p ^ 12 * q ^ 12 + 1424351936 * p ^ 13 * q ^ 11 + 1494405576 * p ^ 14 * q ^ 10 + 1356707184 * p ^ 15 * q ^ 9 + 1062866896 * p ^ 16 * q ^ 8 + 714841344 * p ^ 17 * q ^ 7 + 409280116 * p ^ 18 * q ^ 6 + 196872200 * p ^ 19 * q ^ 5 + 77942340 * p ^ 20 * q ^ 4 + 24578400 * p ^ 21 * q ^ 3 + 5840800 * p ^ 22 * q ^ 2 + 940800 * p ^ 23 * q ^ 1 + 78400 * p ^ 24) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [40, 160, 400, 800, 1344, 2016, 2800, 3680, 4640, 5664, 6596, 7280, 7560, 7280, 6596, 5664, 4640, 3680, 2800, 2016, 1344, 800, 400, 160, 40] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posB4 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [40, 160, 400, 800, 1344, 2016, 2800, 3680, 4640, 5664, 6596, 7280, 7560, 7280, 6596, 5664, 4640, 3680, 2800, 2016, 1344, 800, 400, 160, 40] x := N14_posB4_0 x hu hv

theorem N14_B4_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 14 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posB4 x h0 h1)
  rw [← N14_linkB4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N14_linkA4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 14 4 x * (2 * (1 - x) ^ 2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [-14, 24, 110, 240, 354, 504, 686, 896, 1130, 1384, 1514, 1656, 1806, 1960, 1806, 1656, 1514, 1384, 1130, 896, 686, 504, 354, 240, 110, 24, -14] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (14 - 4 : ℕ) = 10 by norm_num, show (2 * 14 - 4 : ℕ) = 24 by norm_num, show (14 + 4 : ℕ) = 18 by norm_num, show (14 - 4 + 14 : ℕ) = 24 by norm_num, show (4 + 14 : ℕ) = 18 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posA4_0 (x : ℝ) (hu : ((5 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-14, 24, 110, 240, 354, 504, 686, 896, 1130, 1384, 1514, 1656, 1806, 1960, 1806, 1656, 1514, 1384, 1130, 896, 686, 504, 354, 240, 110, 24, -14] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 355075458171949566908290788066146 * q ^ 26 + 12757936274799742496316948038132160 * p ^ 1 * q ^ 25 + 212223683417136200727647567531868160 * p ^ 2 * q ^ 24 + 2204251199172908594211699941856706560 * p ^ 3 * q ^ 23 + 16181480606330614695451829385285730304 * p ^ 4 * q ^ 22 + 89857546086159026941668632165112348672 * p ^ 5 * q ^ 21 + 393845231208402916841720314512397041664 * p ^ 6 * q ^ 20 + 1402194038621006227013044196991760859136 * p ^ 7 * q ^ 19 + 4138550894502158554883751235495934296064 * p ^ 8 * q ^ 18 + 10277750707113523014395230380488713568256 * p ^ 9 * q ^ 17 + 21713116784646352640038685396231671250944 * p ^ 10 * q ^ 16 + 39339091999590498584403108556758671425536 * p ^ 11 * q ^ 15 + 61476745863251520185975640707799474765824 * p ^ 12 * q ^ 14 + 83186209849893865050639860551884478611456 * p ^ 13 * q ^ 13 + 97671500988262154154209641221974850011136 * p ^ 14 * q ^ 12 + 99556276488555334450182656308182811410432 * p ^ 15 * q ^ 11 + 87989490830056996210914452690356588773376 * p ^ 16 * q ^ 10 + 67224738161918686195620769626641896833024 * p ^ 17 * q ^ 9 + 44166935275432709523457070591520786087936 * p ^ 18 * q ^ 8 + 24756722542122384506629934724092610478080 * p ^ 19 * q ^ 7 + 11703424589472264264398985807098003062784 * p ^ 20 * q ^ 6 + 4589269018251611566395029089406681088000 * p ^ 21 * q ^ 5 + 1456858994407181122020178616499956613120 * p ^ 22 * q ^ 4 + 360756464627016295653326533067160944640 * p ^ 23 * q ^ 3 + 65563536337891766185047629036120965120 * p ^ 24 * q ^ 2 + 7800386556456391363757843438024785920 * p ^ 25 * q ^ 1 + 457165512466308651355771043986800640 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (1191817653772720942460132761 : ℝ) * peval [-14, 24, 110, 240, 354, 504, 686, 896, 1130, 1384, 1514, 1656, 1806, 1960, 1806, 1656, 1514, 1384, 1130, 896, 686, 504, 354, 240, 110, 24, -14] x =
      (fun p q : ℝ => 355075458171949566908290788066146 * q ^ 26 + 12757936274799742496316948038132160 * p ^ 1 * q ^ 25 + 212223683417136200727647567531868160 * p ^ 2 * q ^ 24 + 2204251199172908594211699941856706560 * p ^ 3 * q ^ 23 + 16181480606330614695451829385285730304 * p ^ 4 * q ^ 22 + 89857546086159026941668632165112348672 * p ^ 5 * q ^ 21 + 393845231208402916841720314512397041664 * p ^ 6 * q ^ 20 + 1402194038621006227013044196991760859136 * p ^ 7 * q ^ 19 + 4138550894502158554883751235495934296064 * p ^ 8 * q ^ 18 + 10277750707113523014395230380488713568256 * p ^ 9 * q ^ 17 + 21713116784646352640038685396231671250944 * p ^ 10 * q ^ 16 + 39339091999590498584403108556758671425536 * p ^ 11 * q ^ 15 + 61476745863251520185975640707799474765824 * p ^ 12 * q ^ 14 + 83186209849893865050639860551884478611456 * p ^ 13 * q ^ 13 + 97671500988262154154209641221974850011136 * p ^ 14 * q ^ 12 + 99556276488555334450182656308182811410432 * p ^ 15 * q ^ 11 + 87989490830056996210914452690356588773376 * p ^ 16 * q ^ 10 + 67224738161918686195620769626641896833024 * p ^ 17 * q ^ 9 + 44166935275432709523457070591520786087936 * p ^ 18 * q ^ 8 + 24756722542122384506629934724092610478080 * p ^ 19 * q ^ 7 + 11703424589472264264398985807098003062784 * p ^ 20 * q ^ 6 + 4589269018251611566395029089406681088000 * p ^ 21 * q ^ 5 + 1456858994407181122020178616499956613120 * p ^ 22 * q ^ 4 + 360756464627016295653326533067160944640 * p ^ 23 * q ^ 3 + 65563536337891766185047629036120965120 * p ^ 24 * q ^ 2 + 7800386556456391363757843438024785920 * p ^ 25 * q ^ 1 + 457165512466308651355771043986800640 * p ^ 26) (x - ((5 : ℝ) / 16)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1191817653772720942460132761 : ℝ) * peval [-14, 24, 110, 240, 354, 504, 686, 896, 1130, 1384, 1514, 1656, 1806, 1960, 1806, 1656, 1514, 1384, 1130, 896, 686, 504, 354, 240, 110, 24, -14] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posA4 (x : ℝ) (hu : ((5 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-14, 24, 110, 240, 354, 504, 686, 896, 1130, 1384, 1514, 1656, 1806, 1960, 1806, 1656, 1514, 1384, 1130, 896, 686, 504, 354, 240, 110, 24, -14] x := N14_posA4_0 x hu hv

theorem N14_A4_pos (x : ℝ) (h0 : 0 < x) (hlo : ((5 : ℝ) / 16) < x) (h1 : x < 1) :
    0 < aBlock 14 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posA4 x hlo h1)
  rw [← N14_linkA4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N14_linkB5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 14 5 x * (2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [45, 180, 450, 900, 1575, 2450, 3500, 4700, 6025, 7324, 8446, 9240, 9555, 9240, 8446, 7324, 6025, 4700, 3500, 2450, 1575, 900, 450, 180, 45] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (14 - 5 : ℕ) = 9 by norm_num, show (2 * 14 - 5 : ℕ) = 23 by norm_num, show (14 + 5 : ℕ) = 19 by norm_num, show (14 - 5 + 14 : ℕ) = 23 by norm_num, show (5 + 14 : ℕ) = 19 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posB5_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [45, 180, 450, 900, 1575, 2450, 3500, 4700, 6025, 7324, 8446, 9240, 9555, 9240, 8446, 7324, 6025, 4700, 3500, 2450, 1575, 900, 450, 180, 45] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 45 * q ^ 24 + 1260 * p ^ 1 * q ^ 23 + 17010 * p ^ 2 * q ^ 22 + 147420 * p ^ 3 * q ^ 21 + 921375 * p ^ 4 * q ^ 20 + 4422530 * p ^ 5 * q ^ 19 + 16951690 * p ^ 6 * q ^ 18 + 53264090 * p ^ 7 * q ^ 17 + 139740755 * p ^ 8 * q ^ 16 + 310190524 * p ^ 9 * q ^ 15 + 588182126 * p ^ 10 * q ^ 14 + 959250264 * p ^ 11 * q ^ 13 + 1351695891 * p ^ 12 * q ^ 12 + 1650086844 * p ^ 13 * q ^ 11 + 1746578874 * p ^ 14 * q ^ 10 + 1601398176 * p ^ 15 * q ^ 9 + 1268012109 * p ^ 16 * q ^ 8 + 862279236 * p ^ 17 * q ^ 7 + 499118004 * p ^ 18 * q ^ 6 + 242572050 * p ^ 19 * q ^ 5 + 96925185 * p ^ 20 * q ^ 4 + 30803850 * p ^ 21 * q ^ 3 + 7364700 * p ^ 22 * q ^ 2 + 1190700 * p ^ 23 * q ^ 1 + 99225 * p ^ 24 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [45, 180, 450, 900, 1575, 2450, 3500, 4700, 6025, 7324, 8446, 9240, 9555, 9240, 8446, 7324, 6025, 4700, 3500, 2450, 1575, 900, 450, 180, 45] x =
      (fun p q : ℝ => 45 * q ^ 24 + 1260 * p ^ 1 * q ^ 23 + 17010 * p ^ 2 * q ^ 22 + 147420 * p ^ 3 * q ^ 21 + 921375 * p ^ 4 * q ^ 20 + 4422530 * p ^ 5 * q ^ 19 + 16951690 * p ^ 6 * q ^ 18 + 53264090 * p ^ 7 * q ^ 17 + 139740755 * p ^ 8 * q ^ 16 + 310190524 * p ^ 9 * q ^ 15 + 588182126 * p ^ 10 * q ^ 14 + 959250264 * p ^ 11 * q ^ 13 + 1351695891 * p ^ 12 * q ^ 12 + 1650086844 * p ^ 13 * q ^ 11 + 1746578874 * p ^ 14 * q ^ 10 + 1601398176 * p ^ 15 * q ^ 9 + 1268012109 * p ^ 16 * q ^ 8 + 862279236 * p ^ 17 * q ^ 7 + 499118004 * p ^ 18 * q ^ 6 + 242572050 * p ^ 19 * q ^ 5 + 96925185 * p ^ 20 * q ^ 4 + 30803850 * p ^ 21 * q ^ 3 + 7364700 * p ^ 22 * q ^ 2 + 1190700 * p ^ 23 * q ^ 1 + 99225 * p ^ 24) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [45, 180, 450, 900, 1575, 2450, 3500, 4700, 6025, 7324, 8446, 9240, 9555, 9240, 8446, 7324, 6025, 4700, 3500, 2450, 1575, 900, 450, 180, 45] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posB5 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [45, 180, 450, 900, 1575, 2450, 3500, 4700, 6025, 7324, 8446, 9240, 9555, 9240, 8446, 7324, 6025, 4700, 3500, 2450, 1575, 900, 450, 180, 45] x := N14_posB5_0 x hu hv

theorem N14_B5_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 14 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posB5 x h0 h1)
  rw [← N14_linkB5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N14_linkA5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 14 5 x * (2 * (1 - x) ^ 2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [-19, 14, 95, 220, 385, 516, 679, 870, 1085, 1194, 1319, 1456, 1601, 1750, 1601, 1456, 1319, 1194, 1085, 870, 679, 516, 385, 220, 95, 14, -19] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (14 - 5 : ℕ) = 9 by norm_num, show (2 * 14 - 5 : ℕ) = 23 by norm_num, show (14 + 5 : ℕ) = 19 by norm_num, show (14 - 5 + 14 : ℕ) = 23 by norm_num, show (5 + 14 : ℕ) = 19 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posA5_0 (x : ℝ) (hu : ((5 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-19, 14, 95, 220, 385, 516, 679, 870, 1085, 1194, 1319, 1456, 1601, 1750, 1601, 1456, 1319, 1194, 1085, 870, 679, 516, 385, 220, 95, 14, -19] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 154382654228080309776940141983141 * q ^ 26 + 7239869379129550550973498782123520 * p ^ 1 * q ^ 25 + 139314858773497327366846844016401920 * p ^ 2 * q ^ 24 + 1587947389061467972219639989371289600 * p ^ 3 * q ^ 23 + 12437391151094482020956463837086679040 * p ^ 4 * q ^ 22 + 72444822170771404830768591100094447616 * p ^ 5 * q ^ 21 + 329344727663619506500898141686191882240 * p ^ 6 * q ^ 20 + 1206640710295657986449917532751255306240 * p ^ 7 * q ^ 19 + 3643451121505329835667124475990177218560 * p ^ 8 * q ^ 18 + 9214408632777135076364133158030036434944 * p ^ 9 * q ^ 17 + 19750441105623805963685628668667647492096 * p ^ 10 * q ^ 16 + 36190794534118513942655936322616462147584 * p ^ 11 * q ^ 15 + 57044979662322557109991616315199582109696 * p ^ 12 * q ^ 14 + 77667285992337637336522336623694524186624 * p ^ 13 * q ^ 13 + 91556526628121108661335073256635308703744 * p ^ 14 * q ^ 12 + 93513204537929922364534025503328031473664 * p ^ 15 * q ^ 11 + 82671740567338964566493343182473106817024 * p ^ 16 * q ^ 10 + 63083049091391183182152640009110558867456 * p ^ 17 * q ^ 9 + 41341006466893918701421142305056902610944 * p ^ 18 * q ^ 8 + 23091057656844375722280246533464963153920 * p ^ 19 * q ^ 7 + 10870257034139886280549142797195111563264 * p ^ 20 * q ^ 6 + 4243437419311638400122165439496863088640 * p ^ 21 * q ^ 5 + 1341186218981117637398513559653559828480 * p ^ 22 * q ^ 4 + 330841904733515768954211988913447239680 * p ^ 23 * q ^ 3 + 59954776246959895494634980964239933440 * p ^ 24 * q ^ 2 + 7122092073286270375604987486892195840 * p ^ 25 * q ^ 1 + 417411989643151377324834431466209280 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (1191817653772720942460132761 : ℝ) * peval [-19, 14, 95, 220, 385, 516, 679, 870, 1085, 1194, 1319, 1456, 1601, 1750, 1601, 1456, 1319, 1194, 1085, 870, 679, 516, 385, 220, 95, 14, -19] x =
      (fun p q : ℝ => 154382654228080309776940141983141 * q ^ 26 + 7239869379129550550973498782123520 * p ^ 1 * q ^ 25 + 139314858773497327366846844016401920 * p ^ 2 * q ^ 24 + 1587947389061467972219639989371289600 * p ^ 3 * q ^ 23 + 12437391151094482020956463837086679040 * p ^ 4 * q ^ 22 + 72444822170771404830768591100094447616 * p ^ 5 * q ^ 21 + 329344727663619506500898141686191882240 * p ^ 6 * q ^ 20 + 1206640710295657986449917532751255306240 * p ^ 7 * q ^ 19 + 3643451121505329835667124475990177218560 * p ^ 8 * q ^ 18 + 9214408632777135076364133158030036434944 * p ^ 9 * q ^ 17 + 19750441105623805963685628668667647492096 * p ^ 10 * q ^ 16 + 36190794534118513942655936322616462147584 * p ^ 11 * q ^ 15 + 57044979662322557109991616315199582109696 * p ^ 12 * q ^ 14 + 77667285992337637336522336623694524186624 * p ^ 13 * q ^ 13 + 91556526628121108661335073256635308703744 * p ^ 14 * q ^ 12 + 93513204537929922364534025503328031473664 * p ^ 15 * q ^ 11 + 82671740567338964566493343182473106817024 * p ^ 16 * q ^ 10 + 63083049091391183182152640009110558867456 * p ^ 17 * q ^ 9 + 41341006466893918701421142305056902610944 * p ^ 18 * q ^ 8 + 23091057656844375722280246533464963153920 * p ^ 19 * q ^ 7 + 10870257034139886280549142797195111563264 * p ^ 20 * q ^ 6 + 4243437419311638400122165439496863088640 * p ^ 21 * q ^ 5 + 1341186218981117637398513559653559828480 * p ^ 22 * q ^ 4 + 330841904733515768954211988913447239680 * p ^ 23 * q ^ 3 + 59954776246959895494634980964239933440 * p ^ 24 * q ^ 2 + 7122092073286270375604987486892195840 * p ^ 25 * q ^ 1 + 417411989643151377324834431466209280 * p ^ 26) (x - ((5 : ℝ) / 16)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1191817653772720942460132761 : ℝ) * peval [-19, 14, 95, 220, 385, 516, 679, 870, 1085, 1194, 1319, 1456, 1601, 1750, 1601, 1456, 1319, 1194, 1085, 870, 679, 516, 385, 220, 95, 14, -19] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posA5 (x : ℝ) (hu : ((5 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-19, 14, 95, 220, 385, 516, 679, 870, 1085, 1194, 1319, 1456, 1601, 1750, 1601, 1456, 1319, 1194, 1085, 870, 679, 516, 385, 220, 95, 14, -19] x := N14_posA5_0 x hu hv

theorem N14_A5_pos (x : ℝ) (h0 : 0 < x) (hlo : ((5 : ℝ) / 16) < x) (h1 : x < 1) :
    0 < aBlock 14 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posA5 x hlo h1)
  rw [← N14_linkA5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N14_linkB6 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 14 6 x * (2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [48, 192, 480, 960, 1680, 2688, 3948, 5424, 6968, 8432, 9668, 10528, 10864, 10528, 9668, 8432, 6968, 5424, 3948, 2688, 1680, 960, 480, 192, 48] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (14 - 6 : ℕ) = 8 by norm_num, show (2 * 14 - 6 : ℕ) = 22 by norm_num, show (14 + 6 : ℕ) = 20 by norm_num, show (14 - 6 + 14 : ℕ) = 22 by norm_num, show (6 + 14 : ℕ) = 20 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posB6_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [48, 192, 480, 960, 1680, 2688, 3948, 5424, 6968, 8432, 9668, 10528, 10864, 10528, 9668, 8432, 6968, 5424, 3948, 2688, 1680, 960, 480, 192, 48] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 48 * q ^ 24 + 1344 * p ^ 1 * q ^ 23 + 18144 * p ^ 2 * q ^ 22 + 157248 * p ^ 3 * q ^ 21 + 982800 * p ^ 4 * q ^ 20 + 4717440 * p ^ 5 * q ^ 19 + 18083436 * p ^ 6 * q ^ 18 + 56832072 * p ^ 7 * q ^ 17 + 149169524 * p ^ 8 * q ^ 16 + 331399600 * p ^ 9 * q ^ 15 + 629273540 * p ^ 10 * q ^ 14 + 1028420904 * p ^ 11 * q ^ 13 + 1453434108 * p ^ 12 * q ^ 12 + 1781161536 * p ^ 13 * q ^ 11 + 1894395296 * p ^ 14 * q ^ 10 + 1746762304 * p ^ 15 * q ^ 9 + 1391844272 * p ^ 16 * q ^ 8 + 952800576 * p ^ 17 * q ^ 7 + 555194640 * p ^ 18 * q ^ 6 + 271528992 * p ^ 19 * q ^ 5 + 109108496 * p ^ 20 * q ^ 4 + 34840960 * p ^ 21 * q ^ 3 + 8360576 * p ^ 22 * q ^ 2 + 1354752 * p ^ 23 * q ^ 1 + 112896 * p ^ 24 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [48, 192, 480, 960, 1680, 2688, 3948, 5424, 6968, 8432, 9668, 10528, 10864, 10528, 9668, 8432, 6968, 5424, 3948, 2688, 1680, 960, 480, 192, 48] x =
      (fun p q : ℝ => 48 * q ^ 24 + 1344 * p ^ 1 * q ^ 23 + 18144 * p ^ 2 * q ^ 22 + 157248 * p ^ 3 * q ^ 21 + 982800 * p ^ 4 * q ^ 20 + 4717440 * p ^ 5 * q ^ 19 + 18083436 * p ^ 6 * q ^ 18 + 56832072 * p ^ 7 * q ^ 17 + 149169524 * p ^ 8 * q ^ 16 + 331399600 * p ^ 9 * q ^ 15 + 629273540 * p ^ 10 * q ^ 14 + 1028420904 * p ^ 11 * q ^ 13 + 1453434108 * p ^ 12 * q ^ 12 + 1781161536 * p ^ 13 * q ^ 11 + 1894395296 * p ^ 14 * q ^ 10 + 1746762304 * p ^ 15 * q ^ 9 + 1391844272 * p ^ 16 * q ^ 8 + 952800576 * p ^ 17 * q ^ 7 + 555194640 * p ^ 18 * q ^ 6 + 271528992 * p ^ 19 * q ^ 5 + 109108496 * p ^ 20 * q ^ 4 + 34840960 * p ^ 21 * q ^ 3 + 8360576 * p ^ 22 * q ^ 2 + 1354752 * p ^ 23 * q ^ 1 + 112896 * p ^ 24) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [48, 192, 480, 960, 1680, 2688, 3948, 5424, 6968, 8432, 9668, 10528, 10864, 10528, 9668, 8432, 6968, 5424, 3948, 2688, 1680, 960, 480, 192, 48] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posB6 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [48, 192, 480, 960, 1680, 2688, 3948, 5424, 6968, 8432, 9668, 10528, 10864, 10528, 9668, 8432, 6968, 5424, 3948, 2688, 1680, 960, 480, 192, 48] x := N14_posB6_0 x hu hv

theorem N14_B6_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 14 6 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posB6 x h0 h1)
  rw [← N14_linkB6 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N14_linkA6 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 14 6 x * (2 * (1 - x) ^ 2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [-22, 8, 86, 208, 370, 568, 714, 888, 974, 1080, 1202, 1336, 1478, 1624, 1478, 1336, 1202, 1080, 974, 888, 714, 568, 370, 208, 86, 8, -22] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (14 - 6 : ℕ) = 8 by norm_num, show (2 * 14 - 6 : ℕ) = 22 by norm_num, show (14 + 6 : ℕ) = 20 by norm_num, show (14 - 6 + 14 : ℕ) = 22 by norm_num, show (6 + 14 : ℕ) = 20 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posA6_0 (x : ℝ) (hu : ((5 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-22, 8, 86, 208, 370, 568, 714, 888, 974, 1080, 1202, 1336, 1478, 1624, 1478, 1336, 1202, 1080, 974, 888, 714, 568, 370, 208, 86, 8, -22] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2814512106697145785559959199758 * q ^ 26 + 348503174517992588449221195475776 * p ^ 1 * q ^ 25 + 8559530792076646304741612276453376 * p ^ 2 * q ^ 24 + 109562877064442505347992013944061952 * p ^ 3 * q ^ 23 + 918969866509672196900785282396323840 * p ^ 4 * q ^ 22 + 5601497019342551206812957517824393216 * p ^ 5 * q ^ 21 + 26302619501986372789763837776097181696 * p ^ 6 * q ^ 20 + 98721423376903484579328150849938522112 * p ^ 7 * q ^ 19 + 303670440637769657749161591173759893504 * p ^ 8 * q ^ 18 + 779193037551665013701761390698980966400 * p ^ 9 * q ^ 17 + 1689192524549467221089695471978456023040 * p ^ 10 * q ^ 16 + 3122628695345398923302479905515869569024 * p ^ 11 * q ^ 15 + 4954798313704104336028442830961689755648 * p ^ 12 * q ^ 14 + 6778262980430574561671338067706538622976 * p ^ 13 * q ^ 13 + 8015213416746600282593904376396977799168 * p ^ 14 * q ^ 12 + 8199555193464736323073915628005003100160 * p ^ 15 * q ^ 11 + 7250680215367466353820385327597028900864 * p ^ 16 * q ^ 10 + 5527392217380239030957786656199329972224 * p ^ 17 * q ^ 9 + 3615209292892578143143311458750097784832 * p ^ 18 * q ^ 8 + 2013679420747576167812263545286170771456 * p ^ 19 * q ^ 7 + 944793675290214551741742165771037966336 * p ^ 20 * q ^ 6 + 367492098154194527257815248030121590784 * p ^ 21 * q ^ 5 + 115738579710939069216893311115604262912 * p ^ 22 * q ^ 4 + 28460975631324749693283007688003813376 * p ^ 23 * q ^ 3 + 5145526731183082076233600277369847808 * p ^ 24 * q ^ 2 + 610465034853108889337570356019331072 * p ^ 25 * q ^ 1 + 35778170540841546627842951268532224 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (108347059433883722041830251 : ℝ) * peval [-22, 8, 86, 208, 370, 568, 714, 888, 974, 1080, 1202, 1336, 1478, 1624, 1478, 1336, 1202, 1080, 974, 888, 714, 568, 370, 208, 86, 8, -22] x =
      (fun p q : ℝ => 2814512106697145785559959199758 * q ^ 26 + 348503174517992588449221195475776 * p ^ 1 * q ^ 25 + 8559530792076646304741612276453376 * p ^ 2 * q ^ 24 + 109562877064442505347992013944061952 * p ^ 3 * q ^ 23 + 918969866509672196900785282396323840 * p ^ 4 * q ^ 22 + 5601497019342551206812957517824393216 * p ^ 5 * q ^ 21 + 26302619501986372789763837776097181696 * p ^ 6 * q ^ 20 + 98721423376903484579328150849938522112 * p ^ 7 * q ^ 19 + 303670440637769657749161591173759893504 * p ^ 8 * q ^ 18 + 779193037551665013701761390698980966400 * p ^ 9 * q ^ 17 + 1689192524549467221089695471978456023040 * p ^ 10 * q ^ 16 + 3122628695345398923302479905515869569024 * p ^ 11 * q ^ 15 + 4954798313704104336028442830961689755648 * p ^ 12 * q ^ 14 + 6778262980430574561671338067706538622976 * p ^ 13 * q ^ 13 + 8015213416746600282593904376396977799168 * p ^ 14 * q ^ 12 + 8199555193464736323073915628005003100160 * p ^ 15 * q ^ 11 + 7250680215367466353820385327597028900864 * p ^ 16 * q ^ 10 + 5527392217380239030957786656199329972224 * p ^ 17 * q ^ 9 + 3615209292892578143143311458750097784832 * p ^ 18 * q ^ 8 + 2013679420747576167812263545286170771456 * p ^ 19 * q ^ 7 + 944793675290214551741742165771037966336 * p ^ 20 * q ^ 6 + 367492098154194527257815248030121590784 * p ^ 21 * q ^ 5 + 115738579710939069216893311115604262912 * p ^ 22 * q ^ 4 + 28460975631324749693283007688003813376 * p ^ 23 * q ^ 3 + 5145526731183082076233600277369847808 * p ^ 24 * q ^ 2 + 610465034853108889337570356019331072 * p ^ 25 * q ^ 1 + 35778170540841546627842951268532224 * p ^ 26) (x - ((5 : ℝ) / 16)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (108347059433883722041830251 : ℝ) * peval [-22, 8, 86, 208, 370, 568, 714, 888, 974, 1080, 1202, 1336, 1478, 1624, 1478, 1336, 1202, 1080, 974, 888, 714, 568, 370, 208, 86, 8, -22] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posA6 (x : ℝ) (hu : ((5 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-22, 8, 86, 208, 370, 568, 714, 888, 974, 1080, 1202, 1336, 1478, 1624, 1478, 1336, 1202, 1080, 974, 888, 714, 568, 370, 208, 86, 8, -22] x := N14_posA6_0 x hu hv

theorem N14_A6_pos (x : ℝ) (h0 : 0 < x) (hlo : ((5 : ℝ) / 16) < x) (h1 : x < 1) :
    0 < aBlock 14 6 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posA6 x hlo h1)
  rw [← N14_linkA6 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N14_linkB7 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 14 7 x * (2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [49, 196, 490, 980, 1715, 2744, 4116, 5684, 7301, 8820, 10094, 10976, 11319, 10976, 10094, 8820, 7301, 5684, 4116, 2744, 1715, 980, 490, 196, 49] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (14 - 7 : ℕ) = 7 by norm_num, show (2 * 14 - 7 : ℕ) = 21 by norm_num, show (14 + 7 : ℕ) = 21 by norm_num, show (14 - 7 + 14 : ℕ) = 21 by norm_num, show (7 + 14 : ℕ) = 21 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N14_posB7_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [49, 196, 490, 980, 1715, 2744, 4116, 5684, 7301, 8820, 10094, 10976, 11319, 10976, 10094, 8820, 7301, 5684, 4116, 2744, 1715, 980, 490, 196, 49] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 49 * q ^ 24 + 1372 * p ^ 1 * q ^ 23 + 18522 * p ^ 2 * q ^ 22 + 160524 * p ^ 3 * q ^ 21 + 1003275 * p ^ 4 * q ^ 20 + 4815720 * p ^ 5 * q ^ 19 + 18460260 * p ^ 6 * q ^ 18 + 58017764 * p ^ 7 * q ^ 17 + 152293029 * p ^ 8 * q ^ 16 + 338396940 * p ^ 9 * q ^ 15 + 642771710 * p ^ 10 * q ^ 14 + 1051061760 * p ^ 11 * q ^ 13 + 1486677591 * p ^ 12 * q ^ 12 + 1824030096 * p ^ 13 * q ^ 11 + 1942922814 * p ^ 14 * q ^ 10 + 1794785916 * p ^ 15 * q ^ 9 + 1433082469 * p ^ 16 * q ^ 8 + 983209500 * p ^ 17 * q ^ 7 + 574194348 * p ^ 18 * q ^ 6 + 281416408 * p ^ 19 * q ^ 5 + 113295987 * p ^ 20 * q ^ 4 + 36235892 * p ^ 21 * q ^ 3 + 8706026 * p ^ 22 * q ^ 2 + 1411788 * p ^ 23 * q ^ 1 + 117649 * p ^ 24 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [49, 196, 490, 980, 1715, 2744, 4116, 5684, 7301, 8820, 10094, 10976, 11319, 10976, 10094, 8820, 7301, 5684, 4116, 2744, 1715, 980, 490, 196, 49] x =
      (fun p q : ℝ => 49 * q ^ 24 + 1372 * p ^ 1 * q ^ 23 + 18522 * p ^ 2 * q ^ 22 + 160524 * p ^ 3 * q ^ 21 + 1003275 * p ^ 4 * q ^ 20 + 4815720 * p ^ 5 * q ^ 19 + 18460260 * p ^ 6 * q ^ 18 + 58017764 * p ^ 7 * q ^ 17 + 152293029 * p ^ 8 * q ^ 16 + 338396940 * p ^ 9 * q ^ 15 + 642771710 * p ^ 10 * q ^ 14 + 1051061760 * p ^ 11 * q ^ 13 + 1486677591 * p ^ 12 * q ^ 12 + 1824030096 * p ^ 13 * q ^ 11 + 1942922814 * p ^ 14 * q ^ 10 + 1794785916 * p ^ 15 * q ^ 9 + 1433082469 * p ^ 16 * q ^ 8 + 983209500 * p ^ 17 * q ^ 7 + 574194348 * p ^ 18 * q ^ 6 + 281416408 * p ^ 19 * q ^ 5 + 113295987 * p ^ 20 * q ^ 4 + 36235892 * p ^ 21 * q ^ 3 + 8706026 * p ^ 22 * q ^ 2 + 1411788 * p ^ 23 * q ^ 1 + 117649 * p ^ 24) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [49, 196, 490, 980, 1715, 2744, 4116, 5684, 7301, 8820, 10094, 10976, 11319, 10976, 10094, 8820, 7301, 5684, 4116, 2744, 1715, 980, 490, 196, 49] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posB7 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [49, 196, 490, 980, 1715, 2744, 4116, 5684, 7301, 8820, 10094, 10976, 11319, 10976, 10094, 8820, 7301, 5684, 4116, 2744, 1715, 980, 490, 196, 49] x := N14_posB7_0 x hu hv

theorem N14_B7_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 14 7 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N14_posB7 x h0 h1)
  rw [← N14_linkB7 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N14_linkA7 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 14 7 x * (2 * (1 - x) ^ 2 * (1 - x ^ 14) ^ 2) = (1 - x) ^ 4 * peval [-23, 6, 83, 204, 365, 562, 791, 852, 937, 1042, 1163, 1296, 1437, 1582, 1437, 1296, 1163, 1042, 937, 852, 791, 562, 365, 204, 83, 6, -23] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 14 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (14 - 7 : ℕ) = 7 by norm_num, show (2 * 14 - 7 : ℕ) = 21 by norm_num, show (14 + 7 : ℕ) = 21 by norm_num, show (14 - 7 + 14 : ℕ) = 21 by norm_num, show (7 + 14 : ℕ) = 21 by norm_num, show (14 + 14 : ℕ) = 28 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

/-- The minimal polynomial `f_14` of `ξ*(14)` (constant term first). -/
def f14 : List ℤ := [23, -29, -54, -150, -215, -347, -444, -431, -500, -459, -500, -431, -444, -347, -215, -150, -54, -29, 23]

theorem N14_crit_factor (x : ℝ) : peval [-23, 6, 83, 204, 365, 562, 791, 852, 937, 1042, 1163, 1296, 1437, 1582, 1437, 1296, 1163, 1042, 937, 852, 791, 562, 365, 204, 83, 6, -23] x = -(peval [1, 1, 0, 0, 0, 0, 0, 1, 1] x * peval f14 x) := by
  simp only [f14, peval]; ring

theorem N14_posH_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1, 1, 0, 0, 0, 0, 0, 1, 1] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 1 * q ^ 8 + 9 * p ^ 1 * q ^ 7 + 35 * p ^ 2 * q ^ 6 + 77 * p ^ 3 * q ^ 5 + 105 * p ^ 4 * q ^ 4 + 91 * p ^ 5 * q ^ 3 + 49 * p ^ 6 * q ^ 2 + 16 * p ^ 7 * q ^ 1 + 4 * p ^ 8 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [1, 1, 0, 0, 0, 0, 0, 1, 1] x =
      (fun p q : ℝ => 1 * q ^ 8 + 9 * p ^ 1 * q ^ 7 + 35 * p ^ 2 * q ^ 6 + 77 * p ^ 3 * q ^ 5 + 105 * p ^ 4 * q ^ 4 + 91 * p ^ 5 * q ^ 3 + 49 * p ^ 6 * q ^ 2 + 16 * p ^ 7 * q ^ 1 + 4 * p ^ 8) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [1, 1, 0, 0, 0, 0, 0, 1, 1] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_posH (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1, 1, 0, 0, 0, 0, 0, 1, 1] x := N14_posH_0 x hu hv

theorem N14_dneg_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [29, 108, 450, 860, 1735, 2664, 3017, 4000, 4131, 5000, 4741, 5328, 4511, 3010, 2250, 864, 493, -414] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 29 * q ^ 17 + 601 * p ^ 1 * q ^ 16 + 6122 * p ^ 2 * q ^ 15 + 40290 * p ^ 3 * q ^ 14 + 190525 * p ^ 4 * q ^ 13 + 684241 * p ^ 5 * q ^ 12 + 1928253 * p ^ 6 * q ^ 11 + 4350287 * p ^ 7 * q ^ 10 + 7951151 * p ^ 8 * q ^ 9 + 11844889 * p ^ 9 * q ^ 8 + 14402997 * p ^ 10 * q ^ 7 + 14256767 * p ^ 11 * q ^ 6 + 11405499 * p ^ 12 * q ^ 5 + 7277431 * p ^ 13 * q ^ 4 + 3616984 * p ^ 14 * q ^ 3 + 1338484 * p ^ 15 * q ^ 2 + 334621 * p ^ 16 * q ^ 1 + 42777 * p ^ 17 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [29, 108, 450, 860, 1735, 2664, 3017, 4000, 4131, 5000, 4741, 5328, 4511, 3010, 2250, 864, 493, -414] x =
      (fun p q : ℝ => 29 * q ^ 17 + 601 * p ^ 1 * q ^ 16 + 6122 * p ^ 2 * q ^ 15 + 40290 * p ^ 3 * q ^ 14 + 190525 * p ^ 4 * q ^ 13 + 684241 * p ^ 5 * q ^ 12 + 1928253 * p ^ 6 * q ^ 11 + 4350287 * p ^ 7 * q ^ 10 + 7951151 * p ^ 8 * q ^ 9 + 11844889 * p ^ 9 * q ^ 8 + 14402997 * p ^ 10 * q ^ 7 + 14256767 * p ^ 11 * q ^ 6 + 11405499 * p ^ 12 * q ^ 5 + 7277431 * p ^ 13 * q ^ 4 + 3616984 * p ^ 14 * q ^ 3 + 1338484 * p ^ 15 * q ^ 2 + 334621 * p ^ 16 * q ^ 1 + 42777 * p ^ 17) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [29, 108, 450, 860, 1735, 2664, 3017, 4000, 4131, 5000, 4741, 5328, 4511, 3010, 2250, 864, 493, -414] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N14_dneg (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [29, 108, 450, 860, 1735, 2664, 3017, 4000, 4131, 5000, 4741, 5328, 4511, 3010, 2250, 864, 493, -414] x := N14_dneg_0 x hu hv

theorem N14_dlink (x : ℝ) : dpeval f14 x = -peval [29, 108, 450, 860, 1735, 2664, 3017, 4000, 4131, 5000, 4741, 5328, 4511, 3010, 2250, 864, 493, -414] x := by
  simp only [f14, dpeval, peval]; ring

/-- `f_14` changes sign exactly once in `(0,1)`, at a point of `(31493896311550662339/100000000000000000000, 1574694815577533117/5000000000000000000)`. -/
theorem N14_root : ∃ ξs, ((31493896311550662339 : ℝ) / 100000000000000000000) < ξs ∧ ξs < ((1574694815577533117 : ℝ) / 5000000000000000000) ∧ peval f14 ξs = 0 ∧
    ∀ x, 0 < x → x < 1 → (x < ξs → 0 < peval f14 x) ∧ (ξs < x → peval f14 x < 0) := by
  apply sign_change (peval f14) (0 : ℝ) _ _ (by norm_num) (by norm_num) (by norm_num)
    (continuous_peval _).continuousOn
  · apply peval_strictAntiOn
    intro x hu hv
    rw [N14_dlink]
    linarith [N14_dneg x hu hv]
  · exact (fun x hx0 hx => absurd hx (not_le.mpr hx0))
  · norm_num [f14, peval]
  · norm_num [f14, peval]

theorem N14_A7_neg (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : 0 < peval f14 x) :
    aBlock 14 7 x < 0 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : (1 - x) ^ 4 * peval [-23, 6, 83, 204, 365, 562, 791, 852, 937, 1042, 1163, 1296, 1437, 1582, 1437, 1296, 1163, 1042, 937, 852, 791, 562, 365, 204, 83, 6, -23] x < 0 := by
    rw [N14_crit_factor]
    have := mul_pos (N14_posH x h0 h1) hf
    nlinarith [pow_pos h1x 4]
  rw [← N14_linkA7 x h0 h1] at hn
  have hd := (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))
  by_contra hc
  push Not at hc
  nlinarith [mul_nonneg hc hd.le]

theorem N14_A7_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : peval f14 x < 0) :
    0 < aBlock 14 7 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 14 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : 0 < (1 - x) ^ 4 * peval [-23, 6, 83, 204, 365, 562, 791, 852, 937, 1042, 1163, 1296, 1437, 1582, 1437, 1296, 1163, 1042, 937, 852, 791, 562, 365, 204, 83, 6, -23] x := by
    rw [N14_crit_factor]
    have := mul_pos (N14_posH x h0 h1) (neg_pos.mpr hf)
    have h4 := pow_pos h1x 4
    nlinarith
  rw [← N14_linkA7 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

/-- Note Thm 6(b) for `N = 14`: there is a unique `ξ* ∈ (0,1)` (the unique root of `f_14` in
`(0,1)`, enclosed in `(31493896311550662339/100000000000000000000, 1574694815577533117/5000000000000000000)`) such that, for every `ξ ∈ (0,1)`: all `b_m > 0`;
`a_7 < 0` for `ξ < ξ*` (spectrally unstable, since `b_7 > 0`); and all `a_m > 0` for `ξ > ξ*`. -/
theorem cert_N14 : ∃ ξs : ℝ, ((31493896311550662339 : ℝ) / 100000000000000000000) < ξs ∧ ξs < ((1574694815577533117 : ℝ) / 5000000000000000000) ∧ peval f14 ξs = 0 ∧
    (∀ y, 0 < y → y < 1 → peval f14 y = 0 → y = ξs) ∧
    ∀ x, 0 < x → x < 1 →
      (∀ m, 1 ≤ m → m ≤ 13 → 0 < bBlock 14 m x) ∧
      (x < ξs → aBlock 14 7 x < 0) ∧
      (ξs < x → ∀ m, 1 ≤ m → m ≤ 13 → 0 < aBlock 14 m x) := by
  obtain ⟨ξs, hta, htb, hf0, hs⟩ := N14_root
  refine ⟨ξs, hta, htb, hf0, ?_, ?_⟩
  · intro y hy0 hy1 hy
    rcases lt_trichotomy y ξs with h | h | h
    · linarith [(hs y hy0 hy1).1 h]
    · exact h
    · linarith [(hs y hy0 hy1).2 h]
  intro x h0 h1
  refine ⟨?_, fun hlt => N14_A7_neg x h0 h1 ((hs x h0 h1).1 hlt), ?_⟩
  · intro m hm1 hm2
    interval_cases m
    · exact N14_B1_pos x h0 h1
    · exact N14_B2_pos x h0 h1
    · exact N14_B3_pos x h0 h1
    · exact N14_B4_pos x h0 h1
    · exact N14_B5_pos x h0 h1
    · exact N14_B6_pos x h0 h1
    · exact N14_B7_pos x h0 h1
    · rw [show bBlock 14 8 x = bBlock 14 6 x from bBlock_symm 14 6 (by norm_num) x]
      exact N14_B6_pos x h0 h1
    · rw [show bBlock 14 9 x = bBlock 14 5 x from bBlock_symm 14 5 (by norm_num) x]
      exact N14_B5_pos x h0 h1
    · rw [show bBlock 14 10 x = bBlock 14 4 x from bBlock_symm 14 4 (by norm_num) x]
      exact N14_B4_pos x h0 h1
    · rw [show bBlock 14 11 x = bBlock 14 3 x from bBlock_symm 14 3 (by norm_num) x]
      exact N14_B3_pos x h0 h1
    · rw [show bBlock 14 12 x = bBlock 14 2 x from bBlock_symm 14 2 (by norm_num) x]
      exact N14_B2_pos x h0 h1
    · rw [show bBlock 14 13 x = bBlock 14 1 x from bBlock_symm 14 1 (by norm_num) x]
      exact N14_B1_pos x h0 h1
  · intro hgt m hm1 hm2
    have hlo : ((5 : ℝ) / 16) < x := lt_trans (lt_trans (by norm_num) hta) hgt
    have hfn := (hs x h0 h1).2 hgt
    interval_cases m
    · exact N14_A1_pos x h0 h1
    · exact N14_A2_pos x h0 h1
    · exact N14_A3_pos x h0 hlo h1
    · exact N14_A4_pos x h0 hlo h1
    · exact N14_A5_pos x h0 hlo h1
    · exact N14_A6_pos x h0 hlo h1
    · exact N14_A7_pos x h0 h1 hfn
    · rw [show aBlock 14 8 x = aBlock 14 6 x from aBlock_symm 14 6 (by norm_num) x]
      exact N14_A6_pos x h0 hlo h1
    · rw [show aBlock 14 9 x = aBlock 14 5 x from aBlock_symm 14 5 (by norm_num) x]
      exact N14_A5_pos x h0 hlo h1
    · rw [show aBlock 14 10 x = aBlock 14 4 x from aBlock_symm 14 4 (by norm_num) x]
      exact N14_A4_pos x h0 hlo h1
    · rw [show aBlock 14 11 x = aBlock 14 3 x from aBlock_symm 14 3 (by norm_num) x]
      exact N14_A3_pos x h0 hlo h1
    · rw [show aBlock 14 12 x = aBlock 14 2 x from aBlock_symm 14 2 (by norm_num) x]
      exact N14_A2_pos x h0 h1
    · rw [show aBlock 14 13 x = aBlock 14 1 x from aBlock_symm 14 1 (by norm_num) x]
      exact N14_A1_pos x h0 h1

/-! ### N = 15 -/

theorem N15_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 15 1 x * (2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [14, 41, 80, 130, 190, 259, 336, 420, 510, 605, 704, 806, 910, 1015, 910, 806, 704, 605, 510, 420, 336, 259, 190, 130, 80, 41, 14] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (15 - 1 : ℕ) = 14 by norm_num, show (2 * 15 - 1 : ℕ) = 29 by norm_num, show (15 + 1 : ℕ) = 16 by norm_num, show (15 - 1 + 15 : ℕ) = 29 by norm_num, show (1 + 15 : ℕ) = 16 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [14, 41, 80, 130, 190, 259, 336, 420, 510, 605, 704, 806, 910, 1015, 910, 806, 704, 605, 510, 420, 336, 259, 190, 130, 80, 41, 14] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 14 * q ^ 26 + 405 * p ^ 1 * q ^ 25 + 5655 * p ^ 2 * q ^ 24 + 50750 * p ^ 3 * q ^ 23 + 328860 * p ^ 4 * q ^ 22 + 1638819 * p ^ 5 * q ^ 21 + 6531525 * p ^ 6 * q ^ 20 + 21375900 * p ^ 7 * q ^ 19 + 58529250 * p ^ 8 * q ^ 18 + 135917925 * p ^ 9 * q ^ 17 + 270405135 * p ^ 10 * q ^ 16 + 464332050 * p ^ 11 * q ^ 15 + 691945800 * p ^ 12 * q ^ 14 + 898198875 * p ^ 13 * q ^ 13 + 1017958515 * p ^ 14 * q ^ 12 + 1008260942 * p ^ 15 * q ^ 11 + 872516970 * p ^ 16 * q ^ 10 + 658602945 * p ^ 17 * q ^ 9 + 432256825 * p ^ 18 * q ^ 8 + 245405160 * p ^ 19 * q ^ 7 + 119555436 * p ^ 20 * q ^ 6 + 49364315 * p ^ 21 * q ^ 5 + 16941015 * p ^ 22 * q ^ 4 + 4681950 * p ^ 23 * q ^ 3 + 987350 * p ^ 24 * q ^ 2 + 143325 * p ^ 25 * q ^ 1 + 11025 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [14, 41, 80, 130, 190, 259, 336, 420, 510, 605, 704, 806, 910, 1015, 910, 806, 704, 605, 510, 420, 336, 259, 190, 130, 80, 41, 14] x =
      (fun p q : ℝ => 14 * q ^ 26 + 405 * p ^ 1 * q ^ 25 + 5655 * p ^ 2 * q ^ 24 + 50750 * p ^ 3 * q ^ 23 + 328860 * p ^ 4 * q ^ 22 + 1638819 * p ^ 5 * q ^ 21 + 6531525 * p ^ 6 * q ^ 20 + 21375900 * p ^ 7 * q ^ 19 + 58529250 * p ^ 8 * q ^ 18 + 135917925 * p ^ 9 * q ^ 17 + 270405135 * p ^ 10 * q ^ 16 + 464332050 * p ^ 11 * q ^ 15 + 691945800 * p ^ 12 * q ^ 14 + 898198875 * p ^ 13 * q ^ 13 + 1017958515 * p ^ 14 * q ^ 12 + 1008260942 * p ^ 15 * q ^ 11 + 872516970 * p ^ 16 * q ^ 10 + 658602945 * p ^ 17 * q ^ 9 + 432256825 * p ^ 18 * q ^ 8 + 245405160 * p ^ 19 * q ^ 7 + 119555436 * p ^ 20 * q ^ 6 + 49364315 * p ^ 21 * q ^ 5 + 16941015 * p ^ 22 * q ^ 4 + 4681950 * p ^ 23 * q ^ 3 + 987350 * p ^ 24 * q ^ 2 + 143325 * p ^ 25 * q ^ 1 + 11025 * p ^ 26) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [14, 41, 80, 130, 190, 259, 336, 420, 510, 605, 704, 806, 910, 1015, 910, 806, 704, 605, 510, 420, 336, 259, 190, 130, 80, 41, 14] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [14, 41, 80, 130, 190, 259, 336, 420, 510, 605, 704, 806, 910, 1015, 910, 806, 704, 605, 510, 420, 336, 259, 190, 130, 80, 41, 14] x := N15_posB1_0 x hu hv

theorem N15_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 15 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posB1 x h0 h1)
  rw [← N15_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N15_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 15 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [14, 69, 176, 331, 530, 769, 1044, 1351, 1686, 2045, 2424, 2819, 3226, 3641, 3850, 3641, 3226, 2819, 2424, 2045, 1686, 1351, 1044, 769, 530, 331, 176, 69, 14] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (15 - 1 : ℕ) = 14 by norm_num, show (2 * 15 - 1 : ℕ) = 29 by norm_num, show (15 + 1 : ℕ) = 16 by norm_num, show (15 - 1 + 15 : ℕ) = 29 by norm_num, show (1 + 15 : ℕ) = 16 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [14, 69, 176, 331, 530, 769, 1044, 1351, 1686, 2045, 2424, 2819, 3226, 3641, 3850, 3641, 3226, 2819, 2424, 2045, 1686, 1351, 1044, 769, 530, 331, 176, 69, 14] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 14 * q ^ 28 + 461 * p ^ 1 * q ^ 27 + 7331 * p ^ 2 * q ^ 26 + 74990 * p ^ 3 * q ^ 25 + 554480 * p ^ 4 * q ^ 24 + 3157259 * p ^ 5 * q ^ 23 + 14402241 * p ^ 6 * q ^ 22 + 54057276 * p ^ 7 * q ^ 21 + 170158950 * p ^ 8 * q ^ 20 + 455538525 * p ^ 9 * q ^ 19 + 1048193835 * p ^ 10 * q ^ 18 + 2089624290 * p ^ 11 * q ^ 17 + 3630894540 * p ^ 12 * q ^ 16 + 5523310275 * p ^ 13 * q ^ 15 + 7378537215 * p ^ 14 * q ^ 14 + 8672890502 * p ^ 15 * q ^ 13 + 8977394798 * p ^ 16 * q ^ 12 + 8181714593 * p ^ 17 * q ^ 11 + 6556736485 * p ^ 18 * q ^ 10 + 4608844240 * p ^ 19 * q ^ 9 + 2830203376 * p ^ 20 * q ^ 8 + 1509206699 * p ^ 21 * q ^ 7 + 692620019 * p ^ 22 * q ^ 6 + 269903270 * p ^ 23 * q ^ 5 + 87479210 * p ^ 24 * q ^ 4 + 22820525 * p ^ 25 * q ^ 3 + 4533725 * p ^ 26 * q ^ 2 + 617400 * p ^ 27 * q ^ 1 + 44100 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [14, 69, 176, 331, 530, 769, 1044, 1351, 1686, 2045, 2424, 2819, 3226, 3641, 3850, 3641, 3226, 2819, 2424, 2045, 1686, 1351, 1044, 769, 530, 331, 176, 69, 14] x =
      (fun p q : ℝ => 14 * q ^ 28 + 461 * p ^ 1 * q ^ 27 + 7331 * p ^ 2 * q ^ 26 + 74990 * p ^ 3 * q ^ 25 + 554480 * p ^ 4 * q ^ 24 + 3157259 * p ^ 5 * q ^ 23 + 14402241 * p ^ 6 * q ^ 22 + 54057276 * p ^ 7 * q ^ 21 + 170158950 * p ^ 8 * q ^ 20 + 455538525 * p ^ 9 * q ^ 19 + 1048193835 * p ^ 10 * q ^ 18 + 2089624290 * p ^ 11 * q ^ 17 + 3630894540 * p ^ 12 * q ^ 16 + 5523310275 * p ^ 13 * q ^ 15 + 7378537215 * p ^ 14 * q ^ 14 + 8672890502 * p ^ 15 * q ^ 13 + 8977394798 * p ^ 16 * q ^ 12 + 8181714593 * p ^ 17 * q ^ 11 + 6556736485 * p ^ 18 * q ^ 10 + 4608844240 * p ^ 19 * q ^ 9 + 2830203376 * p ^ 20 * q ^ 8 + 1509206699 * p ^ 21 * q ^ 7 + 692620019 * p ^ 22 * q ^ 6 + 269903270 * p ^ 23 * q ^ 5 + 87479210 * p ^ 24 * q ^ 4 + 22820525 * p ^ 25 * q ^ 3 + 4533725 * p ^ 26 * q ^ 2 + 617400 * p ^ 27 * q ^ 1 + 44100 * p ^ 28) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [14, 69, 176, 331, 530, 769, 1044, 1351, 1686, 2045, 2424, 2819, 3226, 3641, 3850, 3641, 3226, 2819, 2424, 2045, 1686, 1351, 1044, 769, 530, 331, 176, 69, 14] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [14, 69, 176, 331, 530, 769, 1044, 1351, 1686, 2045, 2424, 2819, 3226, 3641, 3850, 3641, 3226, 2819, 2424, 2045, 1686, 1351, 1044, 769, 530, 331, 176, 69, 14] x := N15_posA1_0 x hu hv

theorem N15_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 15 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posA1 x h0 h1)
  rw [← N15_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N15_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 15 2 x * (2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [26, 104, 230, 400, 610, 856, 1134, 1440, 1770, 2120, 2486, 2864, 3250, 3445, 3250, 2864, 2486, 2120, 1770, 1440, 1134, 856, 610, 400, 230, 104, 26] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (15 - 2 : ℕ) = 13 by norm_num, show (2 * 15 - 2 : ℕ) = 28 by norm_num, show (15 + 2 : ℕ) = 17 by norm_num, show (15 - 2 + 15 : ℕ) = 28 by norm_num, show (2 + 15 : ℕ) = 17 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [26, 104, 230, 400, 610, 856, 1134, 1440, 1770, 2120, 2486, 2864, 3250, 3445, 3250, 2864, 2486, 2120, 1770, 1440, 1134, 856, 610, 400, 230, 104, 26] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 26 * q ^ 26 + 780 * p ^ 1 * q ^ 25 + 11280 * p ^ 2 * q ^ 24 + 104720 * p ^ 3 * q ^ 23 + 701190 * p ^ 4 * q ^ 22 + 3606876 * p ^ 5 * q ^ 21 + 14823900 * p ^ 6 * q ^ 20 + 49982400 * p ^ 7 * q ^ 19 + 140873850 * p ^ 8 * q ^ 18 + 336464700 * p ^ 9 * q ^ 17 + 687927240 * p ^ 10 * q ^ 16 + 1213102800 * p ^ 11 * q ^ 15 + 1855130550 * p ^ 12 * q ^ 14 + 2469530505 * p ^ 13 * q ^ 13 + 2868333585 * p ^ 14 * q ^ 12 + 2909764598 * p ^ 15 * q ^ 11 + 2577364920 * p ^ 16 * q ^ 10 + 1990068795 * p ^ 17 * q ^ 9 + 1335143095 * p ^ 18 * q ^ 8 + 774199140 * p ^ 19 * q ^ 7 + 384824154 * p ^ 20 * q ^ 6 + 161891015 * p ^ 21 * q ^ 5 + 56502615 * p ^ 22 * q ^ 4 + 15843750 * p ^ 23 * q ^ 3 + 3380000 * p ^ 24 * q ^ 2 + 494325 * p ^ 25 * q ^ 1 + 38025 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [26, 104, 230, 400, 610, 856, 1134, 1440, 1770, 2120, 2486, 2864, 3250, 3445, 3250, 2864, 2486, 2120, 1770, 1440, 1134, 856, 610, 400, 230, 104, 26] x =
      (fun p q : ℝ => 26 * q ^ 26 + 780 * p ^ 1 * q ^ 25 + 11280 * p ^ 2 * q ^ 24 + 104720 * p ^ 3 * q ^ 23 + 701190 * p ^ 4 * q ^ 22 + 3606876 * p ^ 5 * q ^ 21 + 14823900 * p ^ 6 * q ^ 20 + 49982400 * p ^ 7 * q ^ 19 + 140873850 * p ^ 8 * q ^ 18 + 336464700 * p ^ 9 * q ^ 17 + 687927240 * p ^ 10 * q ^ 16 + 1213102800 * p ^ 11 * q ^ 15 + 1855130550 * p ^ 12 * q ^ 14 + 2469530505 * p ^ 13 * q ^ 13 + 2868333585 * p ^ 14 * q ^ 12 + 2909764598 * p ^ 15 * q ^ 11 + 2577364920 * p ^ 16 * q ^ 10 + 1990068795 * p ^ 17 * q ^ 9 + 1335143095 * p ^ 18 * q ^ 8 + 774199140 * p ^ 19 * q ^ 7 + 384824154 * p ^ 20 * q ^ 6 + 161891015 * p ^ 21 * q ^ 5 + 56502615 * p ^ 22 * q ^ 4 + 15843750 * p ^ 23 * q ^ 3 + 3380000 * p ^ 24 * q ^ 2 + 494325 * p ^ 25 * q ^ 1 + 38025 * p ^ 26) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [26, 104, 230, 400, 610, 856, 1134, 1440, 1770, 2120, 2486, 2864, 3250, 3445, 3250, 2864, 2486, 2120, 1770, 1440, 1134, 856, 610, 400, 230, 104, 26] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [26, 104, 230, 400, 610, 856, 1134, 1440, 1770, 2120, 2486, 2864, 3250, 3445, 3250, 2864, 2486, 2120, 1770, 1440, 1134, 856, 610, 400, 230, 104, 26] x := N15_posB2_0 x hu hv

theorem N15_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 15 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posB2 x h0 h1)
  rw [← N15_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N15_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 15 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [2, 60, 140, 268, 440, 652, 900, 1180, 1488, 1820, 2172, 2540, 2920, 3113, 3310, 3113, 2920, 2540, 2172, 1820, 1488, 1180, 900, 652, 440, 268, 140, 60, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (15 - 2 : ℕ) = 13 by norm_num, show (2 * 15 - 2 : ℕ) = 28 by norm_num, show (15 + 2 : ℕ) = 17 by norm_num, show (15 - 2 + 15 : ℕ) = 28 by norm_num, show (2 + 15 : ℕ) = 17 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 60, 140, 268, 440, 652, 900, 1180, 1488, 1820, 2172, 2540, 2920, 3113, 3310, 3113, 2920, 2540, 2172, 1820, 1488, 1180, 900, 652, 440, 268, 140, 60, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 28 + 116 * p ^ 1 * q ^ 27 + 2516 * p ^ 2 * q ^ 26 + 31520 * p ^ 3 * q ^ 25 + 269090 * p ^ 4 * q ^ 24 + 1705172 * p ^ 5 * q ^ 23 + 8444016 * p ^ 6 * q ^ 22 + 33804576 * p ^ 7 * q ^ 21 + 112033350 * p ^ 8 * q ^ 20 + 312713700 * p ^ 9 * q ^ 19 + 744635580 * p ^ 10 * q ^ 18 + 1527339840 * p ^ 11 * q ^ 17 + 2718241890 * p ^ 12 * q ^ 16 + 4220405505 * p ^ 13 * q ^ 15 + 5738768085 * p ^ 14 * q ^ 14 + 6851573726 * p ^ 15 * q ^ 13 + 7192099268 * p ^ 16 * q ^ 12 + 6638949083 * p ^ 17 * q ^ 11 + 5383807195 * p ^ 18 * q ^ 10 + 3826648540 * p ^ 19 * q ^ 9 + 2374535878 * p ^ 20 * q ^ 8 + 1278602039 * p ^ 21 * q ^ 7 + 592020299 * p ^ 22 * q ^ 6 + 232514750 * p ^ 23 * q ^ 5 + 75861500 * p ^ 24 * q ^ 4 + 19895525 * p ^ 25 * q ^ 3 + 3968525 * p ^ 26 * q ^ 2 + 541800 * p ^ 27 * q ^ 1 + 38700 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 60, 140, 268, 440, 652, 900, 1180, 1488, 1820, 2172, 2540, 2920, 3113, 3310, 3113, 2920, 2540, 2172, 1820, 1488, 1180, 900, 652, 440, 268, 140, 60, 2] x =
      (fun p q : ℝ => 2 * q ^ 28 + 116 * p ^ 1 * q ^ 27 + 2516 * p ^ 2 * q ^ 26 + 31520 * p ^ 3 * q ^ 25 + 269090 * p ^ 4 * q ^ 24 + 1705172 * p ^ 5 * q ^ 23 + 8444016 * p ^ 6 * q ^ 22 + 33804576 * p ^ 7 * q ^ 21 + 112033350 * p ^ 8 * q ^ 20 + 312713700 * p ^ 9 * q ^ 19 + 744635580 * p ^ 10 * q ^ 18 + 1527339840 * p ^ 11 * q ^ 17 + 2718241890 * p ^ 12 * q ^ 16 + 4220405505 * p ^ 13 * q ^ 15 + 5738768085 * p ^ 14 * q ^ 14 + 6851573726 * p ^ 15 * q ^ 13 + 7192099268 * p ^ 16 * q ^ 12 + 6638949083 * p ^ 17 * q ^ 11 + 5383807195 * p ^ 18 * q ^ 10 + 3826648540 * p ^ 19 * q ^ 9 + 2374535878 * p ^ 20 * q ^ 8 + 1278602039 * p ^ 21 * q ^ 7 + 592020299 * p ^ 22 * q ^ 6 + 232514750 * p ^ 23 * q ^ 5 + 75861500 * p ^ 24 * q ^ 4 + 19895525 * p ^ 25 * q ^ 3 + 3968525 * p ^ 26 * q ^ 2 + 541800 * p ^ 27 * q ^ 1 + 38700 * p ^ 28) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 60, 140, 268, 440, 652, 900, 1180, 1488, 1820, 2172, 2540, 2920, 3113, 3310, 3113, 2920, 2540, 2172, 1820, 1488, 1180, 900, 652, 440, 268, 140, 60, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 60, 140, 268, 440, 652, 900, 1180, 1488, 1820, 2172, 2540, 2920, 3113, 3310, 3113, 2920, 2540, 2172, 1820, 1488, 1180, 900, 652, 440, 268, 140, 60, 2] x := N15_posA2_0 x hu hv

theorem N15_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 15 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posA2 x h0 h1)
  rw [← N15_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N15_linkB3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 15 3 x * (2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [36, 144, 360, 675, 1080, 1566, 2124, 2745, 3420, 4140, 4896, 5679, 6300, 6570, 6300, 5679, 4896, 4140, 3420, 2745, 2124, 1566, 1080, 675, 360, 144, 36] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (15 - 3 : ℕ) = 12 by norm_num, show (2 * 15 - 3 : ℕ) = 27 by norm_num, show (15 + 3 : ℕ) = 18 by norm_num, show (15 - 3 + 15 : ℕ) = 27 by norm_num, show (3 + 15 : ℕ) = 18 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posB3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [36, 144, 360, 675, 1080, 1566, 2124, 2745, 3420, 4140, 4896, 5679, 6300, 6570, 6300, 5679, 4896, 4140, 3420, 2745, 2124, 1566, 1080, 675, 360, 144, 36] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 36 * q ^ 26 + 1080 * p ^ 1 * q ^ 25 + 15660 * p ^ 2 * q ^ 24 + 146115 * p ^ 3 * q ^ 23 + 985365 * p ^ 4 * q ^ 22 + 5114421 * p ^ 5 * q ^ 21 + 21244275 * p ^ 6 * q ^ 20 + 72499050 * p ^ 7 * q ^ 19 + 207072450 * p ^ 8 * q ^ 18 + 501736950 * p ^ 9 * q ^ 17 + 1041659190 * p ^ 10 * q ^ 16 + 1866679425 * p ^ 11 * q ^ 15 + 2902848795 * p ^ 12 * q ^ 14 + 3931718535 * p ^ 13 * q ^ 13 + 4648483485 * p ^ 14 * q ^ 12 + 4801810518 * p ^ 15 * q ^ 11 + 4332047670 * p ^ 16 * q ^ 10 + 3407228550 * p ^ 17 * q ^ 9 + 2328312690 * p ^ 18 * q ^ 8 + 1374663015 * p ^ 19 * q ^ 7 + 695234529 * p ^ 20 * q ^ 6 + 297250965 * p ^ 21 * q ^ 5 + 105268815 * p ^ 22 * q ^ 4 + 29889000 * p ^ 23 * q ^ 3 + 6439500 * p ^ 24 * q ^ 2 + 947700 * p ^ 25 * q ^ 1 + 72900 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [36, 144, 360, 675, 1080, 1566, 2124, 2745, 3420, 4140, 4896, 5679, 6300, 6570, 6300, 5679, 4896, 4140, 3420, 2745, 2124, 1566, 1080, 675, 360, 144, 36] x =
      (fun p q : ℝ => 36 * q ^ 26 + 1080 * p ^ 1 * q ^ 25 + 15660 * p ^ 2 * q ^ 24 + 146115 * p ^ 3 * q ^ 23 + 985365 * p ^ 4 * q ^ 22 + 5114421 * p ^ 5 * q ^ 21 + 21244275 * p ^ 6 * q ^ 20 + 72499050 * p ^ 7 * q ^ 19 + 207072450 * p ^ 8 * q ^ 18 + 501736950 * p ^ 9 * q ^ 17 + 1041659190 * p ^ 10 * q ^ 16 + 1866679425 * p ^ 11 * q ^ 15 + 2902848795 * p ^ 12 * q ^ 14 + 3931718535 * p ^ 13 * q ^ 13 + 4648483485 * p ^ 14 * q ^ 12 + 4801810518 * p ^ 15 * q ^ 11 + 4332047670 * p ^ 16 * q ^ 10 + 3407228550 * p ^ 17 * q ^ 9 + 2328312690 * p ^ 18 * q ^ 8 + 1374663015 * p ^ 19 * q ^ 7 + 695234529 * p ^ 20 * q ^ 6 + 297250965 * p ^ 21 * q ^ 5 + 105268815 * p ^ 22 * q ^ 4 + 29889000 * p ^ 23 * q ^ 3 + 6439500 * p ^ 24 * q ^ 2 + 947700 * p ^ 25 * q ^ 1 + 72900 * p ^ 26) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [36, 144, 360, 675, 1080, 1566, 2124, 2745, 3420, 4140, 4896, 5679, 6300, 6570, 6300, 5679, 4896, 4140, 3420, 2745, 2124, 1566, 1080, 675, 360, 144, 36] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posB3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [36, 144, 360, 675, 1080, 1566, 2124, 2745, 3420, 4140, 4896, 5679, 6300, 6570, 6300, 5679, 4896, 4140, 3420, 2745, 2124, 1566, 1080, 675, 360, 144, 36] x := N15_posB3_0 x hu hv

theorem N15_B3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 15 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posB3 x h0 h1)
  rw [← N15_linkB3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N15_linkA3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 15 3 x * (2 * (1 - x) ^ 2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [-8, 40, 140, 243, 390, 577, 800, 1055, 1338, 1645, 1972, 2315, 2490, 2673, 2860, 2673, 2490, 2315, 1972, 1645, 1338, 1055, 800, 577, 390, 243, 140, 40, -8] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (15 - 3 : ℕ) = 12 by norm_num, show (2 * 15 - 3 : ℕ) = 27 by norm_num, show (15 + 3 : ℕ) = 18 by norm_num, show (15 - 3 + 15 : ℕ) = 27 by norm_num, show (3 + 15 : ℕ) = 18 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posA3_0 (x : ℝ) (hu : ((21 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-8, 40, 140, 243, 390, 577, 800, 1055, 1338, 1645, 1972, 2315, 2490, 2673, 2860, 2673, 2490, 2315, 1972, 1645, 1338, 1055, 800, 577, 390, 243, 140, 40, -8] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 13928180626317474863283629755696497911249096745134968 * q ^ 28 + 474007057870846457503024140783085915751962172806146560 * p ^ 1 * q ^ 27 + 7727718487926799403100034133309907150720002726805192704 * p ^ 2 * q ^ 26 + 80515420926473057156078695458526940864377479566295040000 * p ^ 3 * q ^ 25 + 603198569637673088258271882986728267234897867852218368000 * p ^ 4 * q ^ 24 + 3464989898591880210904015441700225669442066335406227456000 * p ^ 5 * q ^ 23 + 15888290504767630401652312357054614590861041529699224780800 * p ^ 6 * q ^ 22 + 59764952215911840085825590754521940584561636037096164556800 * p ^ 7 * q ^ 21 + 188053221890165096953814890564878145759704610634225449697280 * p ^ 8 * q ^ 20 + 502141158728947739290755418481360301954823040804381628825600 * p ^ 9 * q ^ 19 + 1150191513215927462092308229048970483264056902673113292996608 * p ^ 10 * q ^ 18 + 2278502874262442257385519008182215012032994574458130750504960 * p ^ 11 * q ^ 17 + 3927371533486257904828926378002879000641480914185758697324544 * p ^ 12 * q ^ 16 + 5916015448691103256729132229849364687629841683245072121856000 * p ^ 13 * q ^ 15 + 7810830837724082430362832838609009657909589914709096870707200 * p ^ 14 * q ^ 14 + 9052940008390508283851435053851675316603156219032695209984000 * p ^ 15 * q ^ 13 + 9213524658077960483947112682593945853154282986180229411635200 * p ^ 16 * q ^ 12 + 8225139617158323939523160650137774429461853244149257889382400 * p ^ 17 * q ^ 11 + 6424571044091627454060082400318024236802422235668725094154240 * p ^ 18 * q ^ 10 + 4372165269029722514256229406980509169977428630822376492236800 * p ^ 19 * q ^ 9 + 2576178178152955016811868036176625417068974097948826572685312 * p ^ 20 * q ^ 8 + 1302611596903250658286660117312991420791232196271463559331840 * p ^ 21 * q ^ 7 + 558222723995450111045332970084916406485064370036699505360896 * p ^ 22 * q ^ 6 + 199229724544469193485303100592289991101572345159605026816000 * p ^ 23 * q ^ 5 + 57745019922766115847433556285430579086609819122229523251200 * p ^ 24 * q ^ 4 + 13086082309977767216196840510202385917555926017024983040000 * p ^ 25 * q ^ 3 + 2179875318333282126626273599282207485277343271492780032000 * p ^ 26 * q ^ 2 + 237920774544372347860214261287549145335828570536148992000 * p ^ 27 * q ^ 1 + 12795739135159521229456901447397601093691620600263475200 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (5459046029871041534743397034507955551606171601 : ℝ) * peval [-8, 40, 140, 243, 390, 577, 800, 1055, 1338, 1645, 1972, 2315, 2490, 2673, 2860, 2673, 2490, 2315, 1972, 1645, 1338, 1055, 800, 577, 390, 243, 140, 40, -8] x =
      (fun p q : ℝ => 13928180626317474863283629755696497911249096745134968 * q ^ 28 + 474007057870846457503024140783085915751962172806146560 * p ^ 1 * q ^ 27 + 7727718487926799403100034133309907150720002726805192704 * p ^ 2 * q ^ 26 + 80515420926473057156078695458526940864377479566295040000 * p ^ 3 * q ^ 25 + 603198569637673088258271882986728267234897867852218368000 * p ^ 4 * q ^ 24 + 3464989898591880210904015441700225669442066335406227456000 * p ^ 5 * q ^ 23 + 15888290504767630401652312357054614590861041529699224780800 * p ^ 6 * q ^ 22 + 59764952215911840085825590754521940584561636037096164556800 * p ^ 7 * q ^ 21 + 188053221890165096953814890564878145759704610634225449697280 * p ^ 8 * q ^ 20 + 502141158728947739290755418481360301954823040804381628825600 * p ^ 9 * q ^ 19 + 1150191513215927462092308229048970483264056902673113292996608 * p ^ 10 * q ^ 18 + 2278502874262442257385519008182215012032994574458130750504960 * p ^ 11 * q ^ 17 + 3927371533486257904828926378002879000641480914185758697324544 * p ^ 12 * q ^ 16 + 5916015448691103256729132229849364687629841683245072121856000 * p ^ 13 * q ^ 15 + 7810830837724082430362832838609009657909589914709096870707200 * p ^ 14 * q ^ 14 + 9052940008390508283851435053851675316603156219032695209984000 * p ^ 15 * q ^ 13 + 9213524658077960483947112682593945853154282986180229411635200 * p ^ 16 * q ^ 12 + 8225139617158323939523160650137774429461853244149257889382400 * p ^ 17 * q ^ 11 + 6424571044091627454060082400318024236802422235668725094154240 * p ^ 18 * q ^ 10 + 4372165269029722514256229406980509169977428630822376492236800 * p ^ 19 * q ^ 9 + 2576178178152955016811868036176625417068974097948826572685312 * p ^ 20 * q ^ 8 + 1302611596903250658286660117312991420791232196271463559331840 * p ^ 21 * q ^ 7 + 558222723995450111045332970084916406485064370036699505360896 * p ^ 22 * q ^ 6 + 199229724544469193485303100592289991101572345159605026816000 * p ^ 23 * q ^ 5 + 57745019922766115847433556285430579086609819122229523251200 * p ^ 24 * q ^ 4 + 13086082309977767216196840510202385917555926017024983040000 * p ^ 25 * q ^ 3 + 2179875318333282126626273599282207485277343271492780032000 * p ^ 26 * q ^ 2 + 237920774544372347860214261287549145335828570536148992000 * p ^ 27 * q ^ 1 + 12795739135159521229456901447397601093691620600263475200 * p ^ 28) (x - ((21 : ℝ) / 64)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (5459046029871041534743397034507955551606171601 : ℝ) * peval [-8, 40, 140, 243, 390, 577, 800, 1055, 1338, 1645, 1972, 2315, 2490, 2673, 2860, 2673, 2490, 2315, 1972, 1645, 1338, 1055, 800, 577, 390, 243, 140, 40, -8] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posA3 (x : ℝ) (hu : ((21 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-8, 40, 140, 243, 390, 577, 800, 1055, 1338, 1645, 1972, 2315, 2490, 2673, 2860, 2673, 2490, 2315, 1972, 1645, 1338, 1055, 800, 577, 390, 243, 140, 40, -8] x := N15_posA3_0 x hu hv

theorem N15_A3_pos (x : ℝ) (h0 : 0 < x) (hlo : ((21 : ℝ) / 64) < x) (h1 : x < 1) :
    0 < aBlock 15 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posA3 x hlo h1)
  rw [← N15_linkA3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N15_linkB4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 15 4 x * (2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [44, 176, 440, 880, 1480, 2224, 3096, 4080, 5160, 6320, 7544, 8651, 9460, 9790, 9460, 8651, 7544, 6320, 5160, 4080, 3096, 2224, 1480, 880, 440, 176, 44] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (15 - 4 : ℕ) = 11 by norm_num, show (2 * 15 - 4 : ℕ) = 26 by norm_num, show (15 + 4 : ℕ) = 19 by norm_num, show (15 - 4 + 15 : ℕ) = 26 by norm_num, show (4 + 15 : ℕ) = 19 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posB4_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [44, 176, 440, 880, 1480, 2224, 3096, 4080, 5160, 6320, 7544, 8651, 9460, 9790, 9460, 8651, 7544, 6320, 5160, 4080, 3096, 2224, 1480, 880, 440, 176, 44] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 44 * q ^ 26 + 1320 * p ^ 1 * q ^ 25 + 19140 * p ^ 2 * q ^ 24 + 178640 * p ^ 3 * q ^ 23 + 1205760 * p ^ 4 * q ^ 22 + 6268704 * p ^ 5 * q ^ 21 + 26106600 * p ^ 6 * q ^ 20 + 89419200 * p ^ 7 * q ^ 19 + 256631700 * p ^ 8 * q ^ 18 + 625567800 * p ^ 9 * q ^ 17 + 1308166860 * p ^ 10 * q ^ 16 + 2364133035 * p ^ 11 * q ^ 15 + 3711962265 * p ^ 12 * q ^ 14 + 5081932185 * p ^ 13 * q ^ 13 + 6079733715 * p ^ 14 * q ^ 12 + 6360962102 * p ^ 15 * q ^ 11 + 5817222510 * p ^ 16 * q ^ 10 + 4640958630 * p ^ 17 * q ^ 9 + 3218090590 * p ^ 18 * q ^ 8 + 1928054535 * p ^ 19 * q ^ 7 + 989119461 * p ^ 20 * q ^ 6 + 428601965 * p ^ 21 * q ^ 5 + 153623415 * p ^ 22 * q ^ 4 + 44068200 * p ^ 23 * q ^ 3 + 9571100 * p ^ 24 * q ^ 2 + 1415700 * p ^ 25 * q ^ 1 + 108900 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [44, 176, 440, 880, 1480, 2224, 3096, 4080, 5160, 6320, 7544, 8651, 9460, 9790, 9460, 8651, 7544, 6320, 5160, 4080, 3096, 2224, 1480, 880, 440, 176, 44] x =
      (fun p q : ℝ => 44 * q ^ 26 + 1320 * p ^ 1 * q ^ 25 + 19140 * p ^ 2 * q ^ 24 + 178640 * p ^ 3 * q ^ 23 + 1205760 * p ^ 4 * q ^ 22 + 6268704 * p ^ 5 * q ^ 21 + 26106600 * p ^ 6 * q ^ 20 + 89419200 * p ^ 7 * q ^ 19 + 256631700 * p ^ 8 * q ^ 18 + 625567800 * p ^ 9 * q ^ 17 + 1308166860 * p ^ 10 * q ^ 16 + 2364133035 * p ^ 11 * q ^ 15 + 3711962265 * p ^ 12 * q ^ 14 + 5081932185 * p ^ 13 * q ^ 13 + 6079733715 * p ^ 14 * q ^ 12 + 6360962102 * p ^ 15 * q ^ 11 + 5817222510 * p ^ 16 * q ^ 10 + 4640958630 * p ^ 17 * q ^ 9 + 3218090590 * p ^ 18 * q ^ 8 + 1928054535 * p ^ 19 * q ^ 7 + 989119461 * p ^ 20 * q ^ 6 + 428601965 * p ^ 21 * q ^ 5 + 153623415 * p ^ 22 * q ^ 4 + 44068200 * p ^ 23 * q ^ 3 + 9571100 * p ^ 24 * q ^ 2 + 1415700 * p ^ 25 * q ^ 1 + 108900 * p ^ 26) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [44, 176, 440, 880, 1480, 2224, 3096, 4080, 5160, 6320, 7544, 8651, 9460, 9790, 9460, 8651, 7544, 6320, 5160, 4080, 3096, 2224, 1480, 880, 440, 176, 44] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posB4 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [44, 176, 440, 880, 1480, 2224, 3096, 4080, 5160, 6320, 7544, 8651, 9460, 9790, 9460, 8651, 7544, 6320, 5160, 4080, 3096, 2224, 1480, 880, 440, 176, 44] x := N15_posB4_0 x hu hv

theorem N15_B4_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 15 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posB4 x h0 h1)
  rw [← N15_linkB4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N15_linkA4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 15 4 x * (2 * (1 - x) ^ 2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [-16, 24, 116, 256, 380, 544, 744, 976, 1236, 1520, 1824, 1979, 2146, 2321, 2500, 2321, 2146, 1979, 1824, 1520, 1236, 976, 744, 544, 380, 256, 116, 24, -16] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (15 - 4 : ℕ) = 11 by norm_num, show (2 * 15 - 4 : ℕ) = 26 by norm_num, show (15 + 4 : ℕ) = 19 by norm_num, show (15 - 4 + 15 : ℕ) = 26 by norm_num, show (4 + 15 : ℕ) = 19 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posA4_0 (x : ℝ) (hu : ((21 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-16, 24, 116, 256, 380, 544, 744, 976, 1236, 1520, 1824, 1979, 2146, 2321, 2500, 2321, 2146, 1979, 1824, 1520, 1236, 976, 744, 544, 380, 256, 116, 24, -16] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 8038338592185349774187846186097317455889041768571120 * q ^ 28 + 300675414101385632321723545766476038306858675432988160 * p ^ 1 * q ^ 27 + 5265168826203597820926885156226066860359320695789830144 * p ^ 2 * q ^ 26 + 58014958537035597683769638986437345322589897362610585600 * p ^ 3 * q ^ 25 + 454575301393469744969476529113096322390173719423635947520 * p ^ 4 * q ^ 24 + 2708592533959442478492617517081560168383179606150506086400 * p ^ 5 * q ^ 23 + 12801006212211830581091275032252022842293224072760178769920 * p ^ 6 * q ^ 22 + 49378029221634105057120559437391739117571810402911649792000 * p ^ 7 * q ^ 21 + 158667256723269121282963938409532304382734711233176613683200 * p ^ 8 * q ^ 20 + 431168392966644223612734260274565571459462100074039869440000 * p ^ 9 * q ^ 19 + 1002127153851587839269826619499167563755480659072942447001600 * p ^ 10 * q ^ 18 + 2009187295041553053493558051723497120757086307532267700879360 * p ^ 11 * q ^ 17 + 3497115310815129619262451933716073135353944300938753636040704 * p ^ 12 * q ^ 16 + 5308833542522491962109529379960371061009228000101535213158400 * p ^ 13 * q ^ 15 + 7050838056135477084727588646389856869157501559364812621742080 * p ^ 14 * q ^ 14 + 8207170555069842117502265225200748071671424280005576962867200 * p ^ 15 * q ^ 13 + 8376143814403673175206710218330758717092410624825436468674560 * p ^ 16 * q ^ 12 + 7488507152380421568129737497394730483799298269995836899328000 * p ^ 17 * q ^ 11 + 5850799113472760817306251830905389645713288426735096614092800 * p ^ 18 * q ^ 10 + 3978718759469036901007932082014809707336364805964397805568000 * p ^ 19 * q ^ 9 + 2340633663399514616366976960467799752557598594344162133278720 * p ^ 20 * q ^ 8 + 1180882895771988387669618676563665003525734113838392556912640 * p ^ 21 * q ^ 7 + 504725842402620832892755001535019788983279073673516119228416 * p ^ 22 * q ^ 6 + 179633599532971743005467163633964777389090624626295098572800 * p ^ 23 * q ^ 5 + 51924970691534282672399567376224092880099608602243718184960 * p ^ 24 * q ^ 4 + 11739694518916383264390859767012477403013725892076634112000 * p ^ 25 * q ^ 3 + 1952215710247069114882194454113038545115443511942551961600 * p ^ 26 * q ^ 2 + 212876482487069995453875917994122919511004510479712256000 * p ^ 27 * q ^ 1 + 11448819226195361100040385505566274662776713168656793600 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (5459046029871041534743397034507955551606171601 : ℝ) * peval [-16, 24, 116, 256, 380, 544, 744, 976, 1236, 1520, 1824, 1979, 2146, 2321, 2500, 2321, 2146, 1979, 1824, 1520, 1236, 976, 744, 544, 380, 256, 116, 24, -16] x =
      (fun p q : ℝ => 8038338592185349774187846186097317455889041768571120 * q ^ 28 + 300675414101385632321723545766476038306858675432988160 * p ^ 1 * q ^ 27 + 5265168826203597820926885156226066860359320695789830144 * p ^ 2 * q ^ 26 + 58014958537035597683769638986437345322589897362610585600 * p ^ 3 * q ^ 25 + 454575301393469744969476529113096322390173719423635947520 * p ^ 4 * q ^ 24 + 2708592533959442478492617517081560168383179606150506086400 * p ^ 5 * q ^ 23 + 12801006212211830581091275032252022842293224072760178769920 * p ^ 6 * q ^ 22 + 49378029221634105057120559437391739117571810402911649792000 * p ^ 7 * q ^ 21 + 158667256723269121282963938409532304382734711233176613683200 * p ^ 8 * q ^ 20 + 431168392966644223612734260274565571459462100074039869440000 * p ^ 9 * q ^ 19 + 1002127153851587839269826619499167563755480659072942447001600 * p ^ 10 * q ^ 18 + 2009187295041553053493558051723497120757086307532267700879360 * p ^ 11 * q ^ 17 + 3497115310815129619262451933716073135353944300938753636040704 * p ^ 12 * q ^ 16 + 5308833542522491962109529379960371061009228000101535213158400 * p ^ 13 * q ^ 15 + 7050838056135477084727588646389856869157501559364812621742080 * p ^ 14 * q ^ 14 + 8207170555069842117502265225200748071671424280005576962867200 * p ^ 15 * q ^ 13 + 8376143814403673175206710218330758717092410624825436468674560 * p ^ 16 * q ^ 12 + 7488507152380421568129737497394730483799298269995836899328000 * p ^ 17 * q ^ 11 + 5850799113472760817306251830905389645713288426735096614092800 * p ^ 18 * q ^ 10 + 3978718759469036901007932082014809707336364805964397805568000 * p ^ 19 * q ^ 9 + 2340633663399514616366976960467799752557598594344162133278720 * p ^ 20 * q ^ 8 + 1180882895771988387669618676563665003525734113838392556912640 * p ^ 21 * q ^ 7 + 504725842402620832892755001535019788983279073673516119228416 * p ^ 22 * q ^ 6 + 179633599532971743005467163633964777389090624626295098572800 * p ^ 23 * q ^ 5 + 51924970691534282672399567376224092880099608602243718184960 * p ^ 24 * q ^ 4 + 11739694518916383264390859767012477403013725892076634112000 * p ^ 25 * q ^ 3 + 1952215710247069114882194454113038545115443511942551961600 * p ^ 26 * q ^ 2 + 212876482487069995453875917994122919511004510479712256000 * p ^ 27 * q ^ 1 + 11448819226195361100040385505566274662776713168656793600 * p ^ 28) (x - ((21 : ℝ) / 64)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (5459046029871041534743397034507955551606171601 : ℝ) * peval [-16, 24, 116, 256, 380, 544, 744, 976, 1236, 1520, 1824, 1979, 2146, 2321, 2500, 2321, 2146, 1979, 1824, 1520, 1236, 976, 744, 544, 380, 256, 116, 24, -16] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posA4 (x : ℝ) (hu : ((21 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-16, 24, 116, 256, 380, 544, 744, 976, 1236, 1520, 1824, 1979, 2146, 2321, 2500, 2321, 2146, 1979, 1824, 1520, 1236, 976, 744, 544, 380, 256, 116, 24, -16] x := N15_posA4_0 x hu hv

theorem N15_A4_pos (x : ℝ) (h0 : 0 < x) (hlo : ((21 : ℝ) / 64) < x) (h1 : x < 1) :
    0 < aBlock 15 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posA4 x hlo h1)
  rw [← N15_linkA4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N15_linkB5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 15 5 x * (2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [50, 200, 500, 1000, 1750, 2725, 3900, 5250, 6750, 8375, 9950, 11300, 12250, 12625, 12250, 11300, 9950, 8375, 6750, 5250, 3900, 2725, 1750, 1000, 500, 200, 50] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (15 - 5 : ℕ) = 10 by norm_num, show (2 * 15 - 5 : ℕ) = 25 by norm_num, show (15 + 5 : ℕ) = 20 by norm_num, show (15 - 5 + 15 : ℕ) = 25 by norm_num, show (5 + 15 : ℕ) = 20 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posB5_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [50, 200, 500, 1000, 1750, 2725, 3900, 5250, 6750, 8375, 9950, 11300, 12250, 12625, 12250, 11300, 9950, 8375, 6750, 5250, 3900, 2725, 1750, 1000, 500, 200, 50] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 50 * q ^ 26 + 1500 * p ^ 1 * q ^ 25 + 21750 * p ^ 2 * q ^ 24 + 203000 * p ^ 3 * q ^ 23 + 1370250 * p ^ 4 * q ^ 22 + 7125225 * p ^ 5 * q ^ 21 + 29686875 * p ^ 6 * q ^ 20 + 101767500 * p ^ 7 * q ^ 19 + 292473750 * p ^ 8 * q ^ 18 + 714408750 * p ^ 9 * q ^ 17 + 1498265850 * p ^ 10 * q ^ 16 + 2718079500 * p ^ 11 * q ^ 15 + 4288580250 * p ^ 12 * q ^ 14 + 5906703375 * p ^ 13 * q ^ 13 + 7117183875 * p ^ 14 * q ^ 12 + 7508393750 * p ^ 15 * q ^ 11 + 6931020000 * p ^ 16 * q ^ 10 + 5586378750 * p ^ 17 * q ^ 9 + 3915902500 * p ^ 18 * q ^ 8 + 2372343750 * p ^ 19 * q ^ 7 + 1230431250 * p ^ 20 * q ^ 6 + 538690625 * p ^ 21 * q ^ 5 + 194878125 * p ^ 22 * q ^ 4 + 56343750 * p ^ 23 * q ^ 3 + 12312500 * p ^ 24 * q ^ 2 + 1828125 * p ^ 25 * q ^ 1 + 140625 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [50, 200, 500, 1000, 1750, 2725, 3900, 5250, 6750, 8375, 9950, 11300, 12250, 12625, 12250, 11300, 9950, 8375, 6750, 5250, 3900, 2725, 1750, 1000, 500, 200, 50] x =
      (fun p q : ℝ => 50 * q ^ 26 + 1500 * p ^ 1 * q ^ 25 + 21750 * p ^ 2 * q ^ 24 + 203000 * p ^ 3 * q ^ 23 + 1370250 * p ^ 4 * q ^ 22 + 7125225 * p ^ 5 * q ^ 21 + 29686875 * p ^ 6 * q ^ 20 + 101767500 * p ^ 7 * q ^ 19 + 292473750 * p ^ 8 * q ^ 18 + 714408750 * p ^ 9 * q ^ 17 + 1498265850 * p ^ 10 * q ^ 16 + 2718079500 * p ^ 11 * q ^ 15 + 4288580250 * p ^ 12 * q ^ 14 + 5906703375 * p ^ 13 * q ^ 13 + 7117183875 * p ^ 14 * q ^ 12 + 7508393750 * p ^ 15 * q ^ 11 + 6931020000 * p ^ 16 * q ^ 10 + 5586378750 * p ^ 17 * q ^ 9 + 3915902500 * p ^ 18 * q ^ 8 + 2372343750 * p ^ 19 * q ^ 7 + 1230431250 * p ^ 20 * q ^ 6 + 538690625 * p ^ 21 * q ^ 5 + 194878125 * p ^ 22 * q ^ 4 + 56343750 * p ^ 23 * q ^ 3 + 12312500 * p ^ 24 * q ^ 2 + 1828125 * p ^ 25 * q ^ 1 + 140625 * p ^ 26) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [50, 200, 500, 1000, 1750, 2725, 3900, 5250, 6750, 8375, 9950, 11300, 12250, 12625, 12250, 11300, 9950, 8375, 6750, 5250, 3900, 2725, 1750, 1000, 500, 200, 50] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posB5 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [50, 200, 500, 1000, 1750, 2725, 3900, 5250, 6750, 8375, 9950, 11300, 12250, 12625, 12250, 11300, 9950, 8375, 6750, 5250, 3900, 2725, 1750, 1000, 500, 200, 50] x := N15_posB5_0 x hu hv

theorem N15_B5_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 15 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posB5 x h0 h1)
  rw [← N15_linkB5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N15_linkA5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 15 5 x * (2 * (1 - x) ^ 2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [-22, 12, 98, 232, 410, 553, 732, 943, 1182, 1445, 1578, 1727, 1888, 2057, 2230, 2057, 1888, 1727, 1578, 1445, 1182, 943, 732, 553, 410, 232, 98, 12, -22] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (15 - 5 : ℕ) = 10 by norm_num, show (2 * 15 - 5 : ℕ) = 25 by norm_num, show (15 + 5 : ℕ) = 20 by norm_num, show (15 - 5 + 15 : ℕ) = 25 by norm_num, show (5 + 15 : ℕ) = 20 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posA5_0 (x : ℝ) (hu : ((21 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-22, 12, 98, 232, 410, 553, 732, 943, 1182, 1445, 1578, 1727, 1888, 2057, 2230, 2057, 1888, 1727, 1578, 1445, 1182, 943, 732, 553, 410, 232, 98, 12, -22] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 3404268495802221971025222109217653233562517380466954 * q ^ 28 + 163931932384612691454773085366813801106935298626786560 * p ^ 1 * q ^ 27 + 3319827423973317988598359656698566310311909323487944704 * p ^ 2 * q ^ 26 + 40247494858184182467968194790964432416106039374564556800 * p ^ 3 * q ^ 25 + 337515335631043991367446825641823026303059583367181762560 * p ^ 4 * q ^ 24 + 2115872434045192322761644943189416673277655164886450176000 * p ^ 5 * q ^ 23 + 10401143905763557590185751854333448970170926507846585548800 * p ^ 6 * q ^ 22 + 41394136840786775840135625107788168353639422681452589875200 * p ^ 7 * q ^ 21 + 136406947489045008170066498042144005447144544889963722833920 * p ^ 8 * q ^ 20 + 378357190254340791269712902419135999887446766410575695052800 * p ^ 9 * q ^ 19 + 894212434714384200212113809394020753925774881605969409736704 * p ^ 10 * q ^ 18 + 1817323010479502258374808350673007854736174465621890653224960 * p ^ 11 * q ^ 17 + 3197729273521713454525806305681756866076761699986001618796544 * p ^ 12 * q ^ 16 + 4895760265565096543092231404770224642882943798733058552627200 * p ^ 13 * q ^ 15 + 6543789041022101982262184626652171861914599734666980900208640 * p ^ 14 * q ^ 14 + 7650951883862329133998289660085361120968905559454389370880000 * p ^ 15 * q ^ 13 + 7829637804776033245091295281722467119961264333061653266432000 * p ^ 16 * q ^ 12 + 7007812861967083130092544520260176572769383163060643535257600 * p ^ 17 * q ^ 11 + 5473693482960155916128385111337359156936925408690028947701760 * p ^ 18 * q ^ 10 + 3716688272257649880890420461673624715174319043580701075046400 * p ^ 19 * q ^ 9 + 2181000186549074378376666806478868899782091253369496431230976 * p ^ 20 * q ^ 8 + 1096746051456122456515502072166809605380960350208685600931840 * p ^ 21 * q ^ 7 + 467003686441146199903467553872832012870498089461524144848896 * p ^ 22 * q ^ 6 + 165554791791318000452859753381989304354063368703585969766400 * p ^ 23 * q ^ 5 + 47674976115530341366159284679547609429756745429404855828480 * p ^ 24 * q ^ 4 + 10743682125209457619981972539341361466405582632736784384000 * p ^ 25 * q ^ 3 + 1782269032031914503987726310713233504025787413908632371200 * p ^ 26 * q ^ 2 + 194093263444093231149122160524053250142386465437384704000 * p ^ 27 * q ^ 1 + 10438629294472241002977998549192779839590532594951782400 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (5459046029871041534743397034507955551606171601 : ℝ) * peval [-22, 12, 98, 232, 410, 553, 732, 943, 1182, 1445, 1578, 1727, 1888, 2057, 2230, 2057, 1888, 1727, 1578, 1445, 1182, 943, 732, 553, 410, 232, 98, 12, -22] x =
      (fun p q : ℝ => 3404268495802221971025222109217653233562517380466954 * q ^ 28 + 163931932384612691454773085366813801106935298626786560 * p ^ 1 * q ^ 27 + 3319827423973317988598359656698566310311909323487944704 * p ^ 2 * q ^ 26 + 40247494858184182467968194790964432416106039374564556800 * p ^ 3 * q ^ 25 + 337515335631043991367446825641823026303059583367181762560 * p ^ 4 * q ^ 24 + 2115872434045192322761644943189416673277655164886450176000 * p ^ 5 * q ^ 23 + 10401143905763557590185751854333448970170926507846585548800 * p ^ 6 * q ^ 22 + 41394136840786775840135625107788168353639422681452589875200 * p ^ 7 * q ^ 21 + 136406947489045008170066498042144005447144544889963722833920 * p ^ 8 * q ^ 20 + 378357190254340791269712902419135999887446766410575695052800 * p ^ 9 * q ^ 19 + 894212434714384200212113809394020753925774881605969409736704 * p ^ 10 * q ^ 18 + 1817323010479502258374808350673007854736174465621890653224960 * p ^ 11 * q ^ 17 + 3197729273521713454525806305681756866076761699986001618796544 * p ^ 12 * q ^ 16 + 4895760265565096543092231404770224642882943798733058552627200 * p ^ 13 * q ^ 15 + 6543789041022101982262184626652171861914599734666980900208640 * p ^ 14 * q ^ 14 + 7650951883862329133998289660085361120968905559454389370880000 * p ^ 15 * q ^ 13 + 7829637804776033245091295281722467119961264333061653266432000 * p ^ 16 * q ^ 12 + 7007812861967083130092544520260176572769383163060643535257600 * p ^ 17 * q ^ 11 + 5473693482960155916128385111337359156936925408690028947701760 * p ^ 18 * q ^ 10 + 3716688272257649880890420461673624715174319043580701075046400 * p ^ 19 * q ^ 9 + 2181000186549074378376666806478868899782091253369496431230976 * p ^ 20 * q ^ 8 + 1096746051456122456515502072166809605380960350208685600931840 * p ^ 21 * q ^ 7 + 467003686441146199903467553872832012870498089461524144848896 * p ^ 22 * q ^ 6 + 165554791791318000452859753381989304354063368703585969766400 * p ^ 23 * q ^ 5 + 47674976115530341366159284679547609429756745429404855828480 * p ^ 24 * q ^ 4 + 10743682125209457619981972539341361466405582632736784384000 * p ^ 25 * q ^ 3 + 1782269032031914503987726310713233504025787413908632371200 * p ^ 26 * q ^ 2 + 194093263444093231149122160524053250142386465437384704000 * p ^ 27 * q ^ 1 + 10438629294472241002977998549192779839590532594951782400 * p ^ 28) (x - ((21 : ℝ) / 64)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (5459046029871041534743397034507955551606171601 : ℝ) * peval [-22, 12, 98, 232, 410, 553, 732, 943, 1182, 1445, 1578, 1727, 1888, 2057, 2230, 2057, 1888, 1727, 1578, 1445, 1182, 943, 732, 553, 410, 232, 98, 12, -22] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posA5 (x : ℝ) (hu : ((21 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-22, 12, 98, 232, 410, 553, 732, 943, 1182, 1445, 1578, 1727, 1888, 2057, 2230, 2057, 1888, 1727, 1578, 1445, 1182, 943, 732, 553, 410, 232, 98, 12, -22] x := N15_posA5_0 x hu hv

theorem N15_A5_pos (x : ℝ) (h0 : 0 < x) (hlo : ((21 : ℝ) / 64) < x) (h1 : x < 1) :
    0 < aBlock 15 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posA5 x hlo h1)
  rw [← N15_linkA5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N15_linkB6 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 15 6 x * (2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [54, 216, 540, 1080, 1890, 3024, 4446, 6120, 8010, 9945, 11754, 13266, 14310, 14715, 14310, 13266, 11754, 9945, 8010, 6120, 4446, 3024, 1890, 1080, 540, 216, 54] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (15 - 6 : ℕ) = 9 by norm_num, show (2 * 15 - 6 : ℕ) = 24 by norm_num, show (15 + 6 : ℕ) = 21 by norm_num, show (15 - 6 + 15 : ℕ) = 24 by norm_num, show (6 + 15 : ℕ) = 21 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posB6_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [54, 216, 540, 1080, 1890, 3024, 4446, 6120, 8010, 9945, 11754, 13266, 14310, 14715, 14310, 13266, 11754, 9945, 8010, 6120, 4446, 3024, 1890, 1080, 540, 216, 54] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 54 * q ^ 26 + 1620 * p ^ 1 * q ^ 25 + 23490 * p ^ 2 * q ^ 24 + 219240 * p ^ 3 * q ^ 23 + 1479870 * p ^ 4 * q ^ 22 + 7695324 * p ^ 5 * q ^ 21 + 32063760 * p ^ 6 * q ^ 20 + 109931040 * p ^ 7 * q ^ 19 + 316033110 * p ^ 8 * q ^ 18 + 772403805 * p ^ 9 * q ^ 17 + 1621471635 * p ^ 10 * q ^ 16 + 2946020490 * p ^ 11 * q ^ 15 + 4658340960 * p ^ 12 * q ^ 14 + 6435074565 * p ^ 13 * q ^ 13 + 7783884945 * p ^ 14 * q ^ 12 + 8251345422 * p ^ 15 * q ^ 11 + 7660618740 * p ^ 16 * q ^ 10 + 6214943700 * p ^ 17 * q ^ 9 + 4387735170 * p ^ 18 * q ^ 8 + 2678071410 * p ^ 19 * q ^ 7 + 1399328136 * p ^ 20 * q ^ 6 + 616941765 * p ^ 21 * q ^ 5 + 224593965 * p ^ 22 * q ^ 4 + 65281950 * p ^ 23 * q ^ 3 + 14324850 * p ^ 24 * q ^ 2 + 2132325 * p ^ 25 * q ^ 1 + 164025 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [54, 216, 540, 1080, 1890, 3024, 4446, 6120, 8010, 9945, 11754, 13266, 14310, 14715, 14310, 13266, 11754, 9945, 8010, 6120, 4446, 3024, 1890, 1080, 540, 216, 54] x =
      (fun p q : ℝ => 54 * q ^ 26 + 1620 * p ^ 1 * q ^ 25 + 23490 * p ^ 2 * q ^ 24 + 219240 * p ^ 3 * q ^ 23 + 1479870 * p ^ 4 * q ^ 22 + 7695324 * p ^ 5 * q ^ 21 + 32063760 * p ^ 6 * q ^ 20 + 109931040 * p ^ 7 * q ^ 19 + 316033110 * p ^ 8 * q ^ 18 + 772403805 * p ^ 9 * q ^ 17 + 1621471635 * p ^ 10 * q ^ 16 + 2946020490 * p ^ 11 * q ^ 15 + 4658340960 * p ^ 12 * q ^ 14 + 6435074565 * p ^ 13 * q ^ 13 + 7783884945 * p ^ 14 * q ^ 12 + 8251345422 * p ^ 15 * q ^ 11 + 7660618740 * p ^ 16 * q ^ 10 + 6214943700 * p ^ 17 * q ^ 9 + 4387735170 * p ^ 18 * q ^ 8 + 2678071410 * p ^ 19 * q ^ 7 + 1399328136 * p ^ 20 * q ^ 6 + 616941765 * p ^ 21 * q ^ 5 + 224593965 * p ^ 22 * q ^ 4 + 65281950 * p ^ 23 * q ^ 3 + 14324850 * p ^ 24 * q ^ 2 + 2132325 * p ^ 25 * q ^ 1 + 164025 * p ^ 26) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [54, 216, 540, 1080, 1890, 3024, 4446, 6120, 8010, 9945, 11754, 13266, 14310, 14715, 14310, 13266, 11754, 9945, 8010, 6120, 4446, 3024, 1890, 1080, 540, 216, 54] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posB6 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [54, 216, 540, 1080, 1890, 3024, 4446, 6120, 8010, 9945, 11754, 13266, 14310, 14715, 14310, 13266, 11754, 9945, 8010, 6120, 4446, 3024, 1890, 1080, 540, 216, 54] x := N15_posB6_0 x hu hv

theorem N15_B6_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 15 6 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posB6 x h0 h1)
  rw [← N15_linkB6 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N15_linkA6 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 15 6 x * (2 * (1 - x) ^ 2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [-26, 4, 86, 216, 390, 604, 764, 956, 1176, 1285, 1414, 1559, 1716, 1881, 2050, 1881, 1716, 1559, 1414, 1285, 1176, 956, 764, 604, 390, 216, 86, 4, -26] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (15 - 6 : ℕ) = 9 by norm_num, show (2 * 15 - 6 : ℕ) = 24 by norm_num, show (15 + 6 : ℕ) = 21 by norm_num, show (15 - 6 + 15 : ℕ) = 24 by norm_num, show (6 + 15 : ℕ) = 21 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posA6_0 (x : ℝ) (hu : ((21 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-26, 4, 86, 216, 390, 604, 764, 956, 1176, 1285, 1414, 1559, 1716, 1881, 2050, 1881, 1716, 1559, 1414, 1285, 1176, 956, 764, 604, 390, 216, 86, 4, -26] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 229179658221124786001506259795856870860434484569030 * q ^ 28 + 69902106886291174001957742498538760625399142631509760 * p ^ 1 * q ^ 27 + 1977779272860173102432192649407004672928364813752131584 * p ^ 2 * q ^ 26 + 27958005396591579366027458196904897971735676472315084800 * p ^ 3 * q ^ 25 + 256413604429268081811268921885211516839716367460634460160 * p ^ 4 * q ^ 24 + 1705104048656473989108089158990931793060026025388487475200 * p ^ 5 * q ^ 23 + 8740443188872131823836380833673811033740383867208529346560 * p ^ 6 * q ^ 22 + 35889590027938062085576812010047217629624077773892380262400 * p ^ 7 * q ^ 21 + 121155798259429588733760023614684505122120018900894655447040 * p ^ 8 * q ^ 20 + 342508463028221488617551521651742352070136881884413558784000 * p ^ 9 * q ^ 19 + 821862973613311156660684384137194064893641823808279564779520 * p ^ 10 * q ^ 18 + 1690665478376643815545292063439062188979423462832031829852160 * p ^ 11 * q ^ 17 + 3003611022247116335803595772836776092338334390255733255962624 * p ^ 12 * q ^ 16 + 4633050699874866476740426091238937374167888456178248712192000 * p ^ 13 * q ^ 15 + 6227337116382366068121960630708756628853163464019799336550400 * p ^ 14 * q ^ 14 + 7309377652026535295349720481227187510988570146143173450137600 * p ^ 15 * q ^ 13 + 7497775041780267220308583379517958653766261131947327073812480 * p ^ 16 * q ^ 12 + 6717291560654702041845005063776257760391945652619151881011200 * p ^ 17 * q ^ 11 + 5245288783815350481202793872904079379857899880271525833605120 * p ^ 18 * q ^ 10 + 3556682604340018530823178087293046686723670650396465954816000 * p ^ 19 * q ^ 9 + 2082303784674524094914374337716767917680191580721346814935040 * p ^ 20 * q ^ 8 + 1043958672847630990417156843243346879692910066974272929136640 * p ^ 21 * q ^ 7 + 442977538244338226249524241955243455442652402287249264214016 * p ^ 22 * q ^ 6 + 156460607892736436002953246971684273590794992386151599308800 * p ^ 23 * q ^ 5 + 44896256924321881530587729242216963644013877151041490780160 * p ^ 24 * q ^ 4 + 10086235029209179723445380258856482008332776095258509312000 * p ^ 25 * q ^ 3 + 1669351259816813881445505270102445221409716073137805721600 * p ^ 26 * q ^ 2 + 181571117415442054945952988877340137229974435409166336000 * p ^ 27 * q ^ 1 + 9765169339990160938269740578277116624133078879148441600 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (5459046029871041534743397034507955551606171601 : ℝ) * peval [-26, 4, 86, 216, 390, 604, 764, 956, 1176, 1285, 1414, 1559, 1716, 1881, 2050, 1881, 1716, 1559, 1414, 1285, 1176, 956, 764, 604, 390, 216, 86, 4, -26] x =
      (fun p q : ℝ => 229179658221124786001506259795856870860434484569030 * q ^ 28 + 69902106886291174001957742498538760625399142631509760 * p ^ 1 * q ^ 27 + 1977779272860173102432192649407004672928364813752131584 * p ^ 2 * q ^ 26 + 27958005396591579366027458196904897971735676472315084800 * p ^ 3 * q ^ 25 + 256413604429268081811268921885211516839716367460634460160 * p ^ 4 * q ^ 24 + 1705104048656473989108089158990931793060026025388487475200 * p ^ 5 * q ^ 23 + 8740443188872131823836380833673811033740383867208529346560 * p ^ 6 * q ^ 22 + 35889590027938062085576812010047217629624077773892380262400 * p ^ 7 * q ^ 21 + 121155798259429588733760023614684505122120018900894655447040 * p ^ 8 * q ^ 20 + 342508463028221488617551521651742352070136881884413558784000 * p ^ 9 * q ^ 19 + 821862973613311156660684384137194064893641823808279564779520 * p ^ 10 * q ^ 18 + 1690665478376643815545292063439062188979423462832031829852160 * p ^ 11 * q ^ 17 + 3003611022247116335803595772836776092338334390255733255962624 * p ^ 12 * q ^ 16 + 4633050699874866476740426091238937374167888456178248712192000 * p ^ 13 * q ^ 15 + 6227337116382366068121960630708756628853163464019799336550400 * p ^ 14 * q ^ 14 + 7309377652026535295349720481227187510988570146143173450137600 * p ^ 15 * q ^ 13 + 7497775041780267220308583379517958653766261131947327073812480 * p ^ 16 * q ^ 12 + 6717291560654702041845005063776257760391945652619151881011200 * p ^ 17 * q ^ 11 + 5245288783815350481202793872904079379857899880271525833605120 * p ^ 18 * q ^ 10 + 3556682604340018530823178087293046686723670650396465954816000 * p ^ 19 * q ^ 9 + 2082303784674524094914374337716767917680191580721346814935040 * p ^ 20 * q ^ 8 + 1043958672847630990417156843243346879692910066974272929136640 * p ^ 21 * q ^ 7 + 442977538244338226249524241955243455442652402287249264214016 * p ^ 22 * q ^ 6 + 156460607892736436002953246971684273590794992386151599308800 * p ^ 23 * q ^ 5 + 44896256924321881530587729242216963644013877151041490780160 * p ^ 24 * q ^ 4 + 10086235029209179723445380258856482008332776095258509312000 * p ^ 25 * q ^ 3 + 1669351259816813881445505270102445221409716073137805721600 * p ^ 26 * q ^ 2 + 181571117415442054945952988877340137229974435409166336000 * p ^ 27 * q ^ 1 + 9765169339990160938269740578277116624133078879148441600 * p ^ 28) (x - ((21 : ℝ) / 64)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (5459046029871041534743397034507955551606171601 : ℝ) * peval [-26, 4, 86, 216, 390, 604, 764, 956, 1176, 1285, 1414, 1559, 1716, 1881, 2050, 1881, 1716, 1559, 1414, 1285, 1176, 956, 764, 604, 390, 216, 86, 4, -26] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posA6 (x : ℝ) (hu : ((21 : ℝ) / 64) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-26, 4, 86, 216, 390, 604, 764, 956, 1176, 1285, 1414, 1559, 1716, 1881, 2050, 1881, 1716, 1559, 1414, 1285, 1176, 956, 764, 604, 390, 216, 86, 4, -26] x := N15_posA6_0 x hu hv

theorem N15_A6_pos (x : ℝ) (h0 : 0 < x) (hlo : ((21 : ℝ) / 64) < x) (h1 : x < 1) :
    0 < aBlock 15 6 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posA6 x hlo h1)
  rw [← N15_linkA6 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N15_linkB7 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 15 7 x * (2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [56, 224, 560, 1120, 1960, 3136, 4704, 6615, 8700, 10790, 12716, 14309, 15400, 15820, 15400, 14309, 12716, 10790, 8700, 6615, 4704, 3136, 1960, 1120, 560, 224, 56] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (15 - 7 : ℕ) = 8 by norm_num, show (2 * 15 - 7 : ℕ) = 23 by norm_num, show (15 + 7 : ℕ) = 22 by norm_num, show (15 - 7 + 15 : ℕ) = 23 by norm_num, show (7 + 15 : ℕ) = 22 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N15_posB7_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [56, 224, 560, 1120, 1960, 3136, 4704, 6615, 8700, 10790, 12716, 14309, 15400, 15820, 15400, 14309, 12716, 10790, 8700, 6615, 4704, 3136, 1960, 1120, 560, 224, 56] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 56 * q ^ 26 + 1680 * p ^ 1 * q ^ 25 + 24360 * p ^ 2 * q ^ 24 + 227360 * p ^ 3 * q ^ 23 + 1534680 * p ^ 4 * q ^ 22 + 7980336 * p ^ 5 * q ^ 21 + 33251400 * p ^ 6 * q ^ 20 + 114004695 * p ^ 7 * q ^ 19 + 327761265 * p ^ 8 * q ^ 18 + 801171195 * p ^ 9 * q ^ 17 + 1682307165 * p ^ 10 * q ^ 16 + 3058014225 * p ^ 11 * q ^ 15 + 4839209655 * p ^ 12 * q ^ 14 + 6692792085 * p ^ 13 * q ^ 13 + 8108974755 * p ^ 14 * q ^ 12 + 8614633208 * p ^ 15 * q ^ 11 + 8019497520 * p ^ 16 * q ^ 10 + 6526769760 * p ^ 17 * q ^ 9 + 4624208680 * p ^ 18 * q ^ 8 + 2832970140 * p ^ 19 * q ^ 7 + 1485811404 * p ^ 20 * q ^ 6 + 657396740 * p ^ 21 * q ^ 5 + 240083340 * p ^ 22 * q ^ 4 + 69972000 * p ^ 23 * q ^ 3 + 15386000 * p ^ 24 * q ^ 2 + 2293200 * p ^ 25 * q ^ 1 + 176400 * p ^ 26 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [56, 224, 560, 1120, 1960, 3136, 4704, 6615, 8700, 10790, 12716, 14309, 15400, 15820, 15400, 14309, 12716, 10790, 8700, 6615, 4704, 3136, 1960, 1120, 560, 224, 56] x =
      (fun p q : ℝ => 56 * q ^ 26 + 1680 * p ^ 1 * q ^ 25 + 24360 * p ^ 2 * q ^ 24 + 227360 * p ^ 3 * q ^ 23 + 1534680 * p ^ 4 * q ^ 22 + 7980336 * p ^ 5 * q ^ 21 + 33251400 * p ^ 6 * q ^ 20 + 114004695 * p ^ 7 * q ^ 19 + 327761265 * p ^ 8 * q ^ 18 + 801171195 * p ^ 9 * q ^ 17 + 1682307165 * p ^ 10 * q ^ 16 + 3058014225 * p ^ 11 * q ^ 15 + 4839209655 * p ^ 12 * q ^ 14 + 6692792085 * p ^ 13 * q ^ 13 + 8108974755 * p ^ 14 * q ^ 12 + 8614633208 * p ^ 15 * q ^ 11 + 8019497520 * p ^ 16 * q ^ 10 + 6526769760 * p ^ 17 * q ^ 9 + 4624208680 * p ^ 18 * q ^ 8 + 2832970140 * p ^ 19 * q ^ 7 + 1485811404 * p ^ 20 * q ^ 6 + 657396740 * p ^ 21 * q ^ 5 + 240083340 * p ^ 22 * q ^ 4 + 69972000 * p ^ 23 * q ^ 3 + 15386000 * p ^ 24 * q ^ 2 + 2293200 * p ^ 25 * q ^ 1 + 176400 * p ^ 26) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [56, 224, 560, 1120, 1960, 3136, 4704, 6615, 8700, 10790, 12716, 14309, 15400, 15820, 15400, 14309, 12716, 10790, 8700, 6615, 4704, 3136, 1960, 1120, 560, 224, 56] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_posB7 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [56, 224, 560, 1120, 1960, 3136, 4704, 6615, 8700, 10790, 12716, 14309, 15400, 15820, 15400, 14309, 12716, 10790, 8700, 6615, 4704, 3136, 1960, 1120, 560, 224, 56] x := N15_posB7_0 x hu hv

theorem N15_B7_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 15 7 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N15_posB7 x h0 h1)
  rw [← N15_linkB7 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N15_linkA7 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 15 7 x * (2 * (1 - x) ^ 2 * (1 - x ^ 15) ^ 2) = (1 - x) ^ 4 * peval [-28, 0, 80, 208, 380, 592, 840, 1015, 1098, 1205, 1332, 1475, 1630, 1793, 1960, 1793, 1630, 1475, 1332, 1205, 1098, 1015, 840, 592, 380, 208, 80, 0, -28] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 15 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (15 - 7 : ℕ) = 8 by norm_num, show (2 * 15 - 7 : ℕ) = 23 by norm_num, show (15 + 7 : ℕ) = 22 by norm_num, show (15 - 7 + 15 : ℕ) = 23 by norm_num, show (7 + 15 : ℕ) = 22 by norm_num, show (15 + 15 : ℕ) = 30 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

/-- The minimal polynomial `f_15` of `ξ*(15)` (constant term first). -/
def f15 : List ℤ := [28, 0, -80, -208, -380, -592, -840, -1015, -1098, -1205, -1332, -1475, -1630, -1793, -1960, -1793, -1630, -1475, -1332, -1205, -1098, -1015, -840, -592, -380, -208, -80, 0, 28]

theorem N15_crit_factor (x : ℝ) : peval [-28, 0, 80, 208, 380, 592, 840, 1015, 1098, 1205, 1332, 1475, 1630, 1793, 1960, 1793, 1630, 1475, 1332, 1205, 1098, 1015, 840, 592, 380, 208, 80, 0, -28] x = -(peval [1] x * peval f15 x) := by
  simp only [f15, peval]; ring

theorem N15_posH (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [1] x := by norm_num [peval]

theorem N15_dneg_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [0, 160, 624, 1520, 2960, 5040, 7105, 8784, 10845, 13320, 16225, 19560, 23309, 27440, 26895, 26080, 25075, 23976, 22895, 21960, 21315, 18480, 13616, 9120, 5200, 2080, 0, -784] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 160 * p ^ 1 * q ^ 26 + 4784 * p ^ 2 * q ^ 25 + 69120 * p ^ 3 * q ^ 24 + 642640 * p ^ 4 * q ^ 23 + 4319840 * p ^ 5 * q ^ 22 + 22361745 * p ^ 6 * q ^ 21 + 92705829 * p ^ 7 * q ^ 20 + 316015455 * p ^ 8 * q ^ 19 + 902313345 * p ^ 9 * q ^ 18 + 2187083025 * p ^ 10 * q ^ 17 + 4544304435 * p ^ 11 * q ^ 16 + 8151349185 * p ^ 12 * q ^ 15 + 12685541295 * p ^ 13 * q ^ 14 + 17184186270 * p ^ 14 * q ^ 13 + 20299884360 * p ^ 15 * q ^ 12 + 20924740000 * p ^ 16 * q ^ 11 + 18809221712 * p ^ 17 * q ^ 10 + 14717622660 * p ^ 18 * q ^ 9 + 9992143700 * p ^ 19 * q ^ 8 + 5856405940 * p ^ 20 * q ^ 7 + 2940394380 * p ^ 21 * q ^ 6 + 1249959256 * p ^ 22 * q ^ 5 + 441862400 * p ^ 23 * q ^ 4 + 126243600 * p ^ 24 * q ^ 3 + 27773200 * p ^ 25 * q ^ 2 + 4272800 * p ^ 26 * q ^ 1 + 352800 * p ^ 27 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [0, 160, 624, 1520, 2960, 5040, 7105, 8784, 10845, 13320, 16225, 19560, 23309, 27440, 26895, 26080, 25075, 23976, 22895, 21960, 21315, 18480, 13616, 9120, 5200, 2080, 0, -784] x =
      (fun p q : ℝ => 160 * p ^ 1 * q ^ 26 + 4784 * p ^ 2 * q ^ 25 + 69120 * p ^ 3 * q ^ 24 + 642640 * p ^ 4 * q ^ 23 + 4319840 * p ^ 5 * q ^ 22 + 22361745 * p ^ 6 * q ^ 21 + 92705829 * p ^ 7 * q ^ 20 + 316015455 * p ^ 8 * q ^ 19 + 902313345 * p ^ 9 * q ^ 18 + 2187083025 * p ^ 10 * q ^ 17 + 4544304435 * p ^ 11 * q ^ 16 + 8151349185 * p ^ 12 * q ^ 15 + 12685541295 * p ^ 13 * q ^ 14 + 17184186270 * p ^ 14 * q ^ 13 + 20299884360 * p ^ 15 * q ^ 12 + 20924740000 * p ^ 16 * q ^ 11 + 18809221712 * p ^ 17 * q ^ 10 + 14717622660 * p ^ 18 * q ^ 9 + 9992143700 * p ^ 19 * q ^ 8 + 5856405940 * p ^ 20 * q ^ 7 + 2940394380 * p ^ 21 * q ^ 6 + 1249959256 * p ^ 22 * q ^ 5 + 441862400 * p ^ 23 * q ^ 4 + 126243600 * p ^ 24 * q ^ 3 + 27773200 * p ^ 25 * q ^ 2 + 4272800 * p ^ 26 * q ^ 1 + 352800 * p ^ 27) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [0, 160, 624, 1520, 2960, 5040, 7105, 8784, 10845, 13320, 16225, 19560, 23309, 27440, 26895, 26080, 25075, 23976, 22895, 21960, 21315, 18480, 13616, 9120, 5200, 2080, 0, -784] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N15_dneg (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [0, 160, 624, 1520, 2960, 5040, 7105, 8784, 10845, 13320, 16225, 19560, 23309, 27440, 26895, 26080, 25075, 23976, 22895, 21960, 21315, 18480, 13616, 9120, 5200, 2080, 0, -784] x := N15_dneg_0 x hu hv

theorem N15_dlink (x : ℝ) : dpeval f15 x = -peval [0, 160, 624, 1520, 2960, 5040, 7105, 8784, 10845, 13320, 16225, 19560, 23309, 27440, 26895, 26080, 25075, 23976, 22895, 21960, 21315, 18480, 13616, 9120, 5200, 2080, 0, -784] x := by
  simp only [f15, dpeval, peval]; ring

/-- `f_15` changes sign exactly once in `(0,1)`, at a point of `(1712292917560851283/5000000000000000000, 34245858351217025661/100000000000000000000)`. -/
theorem N15_root : ∃ ξs, ((1712292917560851283 : ℝ) / 5000000000000000000) < ξs ∧ ξs < ((34245858351217025661 : ℝ) / 100000000000000000000) ∧ peval f15 ξs = 0 ∧
    ∀ x, 0 < x → x < 1 → (x < ξs → 0 < peval f15 x) ∧ (ξs < x → peval f15 x < 0) := by
  apply sign_change (peval f15) (0 : ℝ) _ _ (by norm_num) (by norm_num) (by norm_num)
    (continuous_peval _).continuousOn
  · apply peval_strictAntiOn
    intro x hu hv
    rw [N15_dlink]
    linarith [N15_dneg x hu hv]
  · exact (fun x hx0 hx => absurd hx (not_le.mpr hx0))
  · norm_num [f15, peval]
  · norm_num [f15, peval]

theorem N15_A7_neg (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : 0 < peval f15 x) :
    aBlock 15 7 x < 0 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : (1 - x) ^ 4 * peval [-28, 0, 80, 208, 380, 592, 840, 1015, 1098, 1205, 1332, 1475, 1630, 1793, 1960, 1793, 1630, 1475, 1332, 1205, 1098, 1015, 840, 592, 380, 208, 80, 0, -28] x < 0 := by
    rw [N15_crit_factor]
    have := mul_pos (N15_posH x h0 h1) hf
    nlinarith [pow_pos h1x 4]
  rw [← N15_linkA7 x h0 h1] at hn
  have hd := (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))
  by_contra hc
  push Not at hc
  nlinarith [mul_nonneg hc hd.le]

theorem N15_A7_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : peval f15 x < 0) :
    0 < aBlock 15 7 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 15 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : 0 < (1 - x) ^ 4 * peval [-28, 0, 80, 208, 380, 592, 840, 1015, 1098, 1205, 1332, 1475, 1630, 1793, 1960, 1793, 1630, 1475, 1332, 1205, 1098, 1015, 840, 592, 380, 208, 80, 0, -28] x := by
    rw [N15_crit_factor]
    have := mul_pos (N15_posH x h0 h1) (neg_pos.mpr hf)
    have h4 := pow_pos h1x 4
    nlinarith
  rw [← N15_linkA7 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

/-- Note Thm 6(b) for `N = 15`: there is a unique `ξ* ∈ (0,1)` (the unique root of `f_15` in
`(0,1)`, enclosed in `(1712292917560851283/5000000000000000000, 34245858351217025661/100000000000000000000)`) such that, for every `ξ ∈ (0,1)`: all `b_m > 0`;
`a_7 < 0` for `ξ < ξ*` (spectrally unstable, since `b_7 > 0`); and all `a_m > 0` for `ξ > ξ*`. -/
theorem cert_N15 : ∃ ξs : ℝ, ((1712292917560851283 : ℝ) / 5000000000000000000) < ξs ∧ ξs < ((34245858351217025661 : ℝ) / 100000000000000000000) ∧ peval f15 ξs = 0 ∧
    (∀ y, 0 < y → y < 1 → peval f15 y = 0 → y = ξs) ∧
    ∀ x, 0 < x → x < 1 →
      (∀ m, 1 ≤ m → m ≤ 14 → 0 < bBlock 15 m x) ∧
      (x < ξs → aBlock 15 7 x < 0) ∧
      (ξs < x → ∀ m, 1 ≤ m → m ≤ 14 → 0 < aBlock 15 m x) := by
  obtain ⟨ξs, hta, htb, hf0, hs⟩ := N15_root
  refine ⟨ξs, hta, htb, hf0, ?_, ?_⟩
  · intro y hy0 hy1 hy
    rcases lt_trichotomy y ξs with h | h | h
    · linarith [(hs y hy0 hy1).1 h]
    · exact h
    · linarith [(hs y hy0 hy1).2 h]
  intro x h0 h1
  refine ⟨?_, fun hlt => N15_A7_neg x h0 h1 ((hs x h0 h1).1 hlt), ?_⟩
  · intro m hm1 hm2
    interval_cases m
    · exact N15_B1_pos x h0 h1
    · exact N15_B2_pos x h0 h1
    · exact N15_B3_pos x h0 h1
    · exact N15_B4_pos x h0 h1
    · exact N15_B5_pos x h0 h1
    · exact N15_B6_pos x h0 h1
    · exact N15_B7_pos x h0 h1
    · rw [show bBlock 15 8 x = bBlock 15 7 x from bBlock_symm 15 7 (by norm_num) x]
      exact N15_B7_pos x h0 h1
    · rw [show bBlock 15 9 x = bBlock 15 6 x from bBlock_symm 15 6 (by norm_num) x]
      exact N15_B6_pos x h0 h1
    · rw [show bBlock 15 10 x = bBlock 15 5 x from bBlock_symm 15 5 (by norm_num) x]
      exact N15_B5_pos x h0 h1
    · rw [show bBlock 15 11 x = bBlock 15 4 x from bBlock_symm 15 4 (by norm_num) x]
      exact N15_B4_pos x h0 h1
    · rw [show bBlock 15 12 x = bBlock 15 3 x from bBlock_symm 15 3 (by norm_num) x]
      exact N15_B3_pos x h0 h1
    · rw [show bBlock 15 13 x = bBlock 15 2 x from bBlock_symm 15 2 (by norm_num) x]
      exact N15_B2_pos x h0 h1
    · rw [show bBlock 15 14 x = bBlock 15 1 x from bBlock_symm 15 1 (by norm_num) x]
      exact N15_B1_pos x h0 h1
  · intro hgt m hm1 hm2
    have hlo : ((21 : ℝ) / 64) < x := lt_trans (lt_trans (by norm_num) hta) hgt
    have hfn := (hs x h0 h1).2 hgt
    interval_cases m
    · exact N15_A1_pos x h0 h1
    · exact N15_A2_pos x h0 h1
    · exact N15_A3_pos x h0 hlo h1
    · exact N15_A4_pos x h0 hlo h1
    · exact N15_A5_pos x h0 hlo h1
    · exact N15_A6_pos x h0 hlo h1
    · exact N15_A7_pos x h0 h1 hfn
    · rw [show aBlock 15 8 x = aBlock 15 7 x from aBlock_symm 15 7 (by norm_num) x]
      exact N15_A7_pos x h0 h1 hfn
    · rw [show aBlock 15 9 x = aBlock 15 6 x from aBlock_symm 15 6 (by norm_num) x]
      exact N15_A6_pos x h0 hlo h1
    · rw [show aBlock 15 10 x = aBlock 15 5 x from aBlock_symm 15 5 (by norm_num) x]
      exact N15_A5_pos x h0 hlo h1
    · rw [show aBlock 15 11 x = aBlock 15 4 x from aBlock_symm 15 4 (by norm_num) x]
      exact N15_A4_pos x h0 hlo h1
    · rw [show aBlock 15 12 x = aBlock 15 3 x from aBlock_symm 15 3 (by norm_num) x]
      exact N15_A3_pos x h0 hlo h1
    · rw [show aBlock 15 13 x = aBlock 15 2 x from aBlock_symm 15 2 (by norm_num) x]
      exact N15_A2_pos x h0 h1
    · rw [show aBlock 15 14 x = aBlock 15 1 x from aBlock_symm 15 1 (by norm_num) x]
      exact N15_A1_pos x h0 h1

/-! ### N = 16 -/

theorem N16_linkB1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 16 1 x * (2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [15, 44, 86, 140, 205, 280, 364, 456, 555, 660, 770, 884, 1001, 1120, 1240, 1120, 1001, 884, 770, 660, 555, 456, 364, 280, 205, 140, 86, 44, 15] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (16 - 1 : ℕ) = 15 by norm_num, show (2 * 16 - 1 : ℕ) = 31 by norm_num, show (16 + 1 : ℕ) = 17 by norm_num, show (16 - 1 + 16 : ℕ) = 31 by norm_num, show (1 + 16 : ℕ) = 17 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posB1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [15, 44, 86, 140, 205, 280, 364, 456, 555, 660, 770, 884, 1001, 1120, 1240, 1120, 1001, 884, 770, 660, 555, 456, 364, 280, 205, 140, 86, 44, 15] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 15 * q ^ 28 + 464 * p ^ 1 * q ^ 27 + 6944 * p ^ 2 * q ^ 26 + 66960 * p ^ 3 * q ^ 25 + 467480 * p ^ 4 * q ^ 24 + 2517200 * p ^ 5 * q ^ 23 + 10874304 * p ^ 6 * q ^ 22 + 38707344 * p ^ 7 * q ^ 21 + 115701300 * p ^ 8 * q ^ 20 + 294512400 * p ^ 9 * q ^ 19 + 645122400 * p ^ 10 * q ^ 18 + 1225732560 * p ^ 11 * q ^ 17 + 2032135560 * p ^ 12 * q ^ 16 + 2952675600 * p ^ 13 * q ^ 15 + 3771484800 * p ^ 14 * q ^ 14 + 4242920160 * p ^ 15 * q ^ 13 + 4207559132 * p ^ 16 * q ^ 12 + 3677172512 * p ^ 17 * q ^ 11 + 2828504640 * p ^ 18 * q ^ 10 + 1910228320 * p ^ 19 * q ^ 9 + 1128247120 * p ^ 20 * q ^ 8 + 579417696 * p ^ 21 * q ^ 7 + 256520576 * p ^ 22 * q ^ 6 + 96649280 * p ^ 23 * q ^ 5 + 30379440 * p ^ 24 * q ^ 4 + 7716800 * p ^ 25 * q ^ 3 + 1500800 * p ^ 26 * q ^ 2 + 201600 * p ^ 27 * q ^ 1 + 14400 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [15, 44, 86, 140, 205, 280, 364, 456, 555, 660, 770, 884, 1001, 1120, 1240, 1120, 1001, 884, 770, 660, 555, 456, 364, 280, 205, 140, 86, 44, 15] x =
      (fun p q : ℝ => 15 * q ^ 28 + 464 * p ^ 1 * q ^ 27 + 6944 * p ^ 2 * q ^ 26 + 66960 * p ^ 3 * q ^ 25 + 467480 * p ^ 4 * q ^ 24 + 2517200 * p ^ 5 * q ^ 23 + 10874304 * p ^ 6 * q ^ 22 + 38707344 * p ^ 7 * q ^ 21 + 115701300 * p ^ 8 * q ^ 20 + 294512400 * p ^ 9 * q ^ 19 + 645122400 * p ^ 10 * q ^ 18 + 1225732560 * p ^ 11 * q ^ 17 + 2032135560 * p ^ 12 * q ^ 16 + 2952675600 * p ^ 13 * q ^ 15 + 3771484800 * p ^ 14 * q ^ 14 + 4242920160 * p ^ 15 * q ^ 13 + 4207559132 * p ^ 16 * q ^ 12 + 3677172512 * p ^ 17 * q ^ 11 + 2828504640 * p ^ 18 * q ^ 10 + 1910228320 * p ^ 19 * q ^ 9 + 1128247120 * p ^ 20 * q ^ 8 + 579417696 * p ^ 21 * q ^ 7 + 256520576 * p ^ 22 * q ^ 6 + 96649280 * p ^ 23 * q ^ 5 + 30379440 * p ^ 24 * q ^ 4 + 7716800 * p ^ 25 * q ^ 3 + 1500800 * p ^ 26 * q ^ 2 + 201600 * p ^ 27 * q ^ 1 + 14400 * p ^ 28) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [15, 44, 86, 140, 205, 280, 364, 456, 555, 660, 770, 884, 1001, 1120, 1240, 1120, 1001, 884, 770, 660, 555, 456, 364, 280, 205, 140, 86, 44, 15] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posB1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [15, 44, 86, 140, 205, 280, 364, 456, 555, 660, 770, 884, 1001, 1120, 1240, 1120, 1001, 884, 770, 660, 555, 456, 364, 280, 205, 140, 86, 44, 15] x := N16_posB1_0 x hu hv

theorem N16_B1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 16 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posB1 x h0 h1)
  rw [← N16_linkB1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N16_linkA1 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 16 1 x * (2 * (1 - x) ^ 2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [15, 74, 189, 356, 571, 830, 1129, 1464, 1831, 2226, 2645, 3084, 3539, 4006, 4481, 4720, 4481, 4006, 3539, 3084, 2645, 2226, 1831, 1464, 1129, 830, 571, 356, 189, 74, 15] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (16 - 1 : ℕ) = 15 by norm_num, show (2 * 16 - 1 : ℕ) = 31 by norm_num, show (16 + 1 : ℕ) = 17 by norm_num, show (16 - 1 + 16 : ℕ) = 31 by norm_num, show (1 + 16 : ℕ) = 17 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posA1_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [15, 74, 189, 356, 571, 830, 1129, 1464, 1831, 2226, 2645, 3084, 3539, 4006, 4481, 4720, 4481, 4006, 3539, 3084, 2645, 2226, 1831, 1464, 1129, 830, 571, 356, 189, 74, 15] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 15 * q ^ 30 + 524 * p ^ 1 * q ^ 29 + 8860 * p ^ 2 * q ^ 28 + 96592 * p ^ 3 * q ^ 27 + 763096 * p ^ 4 * q ^ 26 + 4654960 * p ^ 5 * q ^ 25 + 22813024 * p ^ 6 * q ^ 24 + 92273360 * p ^ 7 * q ^ 23 + 314027892 * p ^ 8 * q ^ 22 + 912146976 * p ^ 9 * q ^ 21 + 2285977200 * p ^ 10 * q ^ 20 + 4984271760 * p ^ 11 * q ^ 19 + 9515555400 * p ^ 12 * q ^ 18 + 15984148080 * p ^ 13 * q ^ 17 + 23710729440 * p ^ 14 * q ^ 16 + 31139561760 * p ^ 15 * q ^ 15 + 36265178972 * p ^ 16 * q ^ 14 + 37479089680 * p ^ 17 * q ^ 13 + 34367431216 * p ^ 18 * q ^ 12 + 27932936928 * p ^ 19 * q ^ 11 + 20083178960 * p ^ 20 * q ^ 10 + 12733319456 * p ^ 21 * q ^ 9 + 7087179840 * p ^ 22 * q ^ 8 + 3440402368 * p ^ 23 * q ^ 7 + 1443058864 * p ^ 24 * q ^ 6 + 515831680 * p ^ 25 * q ^ 5 + 153885760 * p ^ 26 * q ^ 4 + 37072000 * p ^ 27 * q ^ 3 + 6824000 * p ^ 28 * q ^ 2 + 864000 * p ^ 29 * q ^ 1 + 57600 * p ^ 30 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [15, 74, 189, 356, 571, 830, 1129, 1464, 1831, 2226, 2645, 3084, 3539, 4006, 4481, 4720, 4481, 4006, 3539, 3084, 2645, 2226, 1831, 1464, 1129, 830, 571, 356, 189, 74, 15] x =
      (fun p q : ℝ => 15 * q ^ 30 + 524 * p ^ 1 * q ^ 29 + 8860 * p ^ 2 * q ^ 28 + 96592 * p ^ 3 * q ^ 27 + 763096 * p ^ 4 * q ^ 26 + 4654960 * p ^ 5 * q ^ 25 + 22813024 * p ^ 6 * q ^ 24 + 92273360 * p ^ 7 * q ^ 23 + 314027892 * p ^ 8 * q ^ 22 + 912146976 * p ^ 9 * q ^ 21 + 2285977200 * p ^ 10 * q ^ 20 + 4984271760 * p ^ 11 * q ^ 19 + 9515555400 * p ^ 12 * q ^ 18 + 15984148080 * p ^ 13 * q ^ 17 + 23710729440 * p ^ 14 * q ^ 16 + 31139561760 * p ^ 15 * q ^ 15 + 36265178972 * p ^ 16 * q ^ 14 + 37479089680 * p ^ 17 * q ^ 13 + 34367431216 * p ^ 18 * q ^ 12 + 27932936928 * p ^ 19 * q ^ 11 + 20083178960 * p ^ 20 * q ^ 10 + 12733319456 * p ^ 21 * q ^ 9 + 7087179840 * p ^ 22 * q ^ 8 + 3440402368 * p ^ 23 * q ^ 7 + 1443058864 * p ^ 24 * q ^ 6 + 515831680 * p ^ 25 * q ^ 5 + 153885760 * p ^ 26 * q ^ 4 + 37072000 * p ^ 27 * q ^ 3 + 6824000 * p ^ 28 * q ^ 2 + 864000 * p ^ 29 * q ^ 1 + 57600 * p ^ 30) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [15, 74, 189, 356, 571, 830, 1129, 1464, 1831, 2226, 2645, 3084, 3539, 4006, 4481, 4720, 4481, 4006, 3539, 3084, 2645, 2226, 1831, 1464, 1129, 830, 571, 356, 189, 74, 15] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posA1 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [15, 74, 189, 356, 571, 830, 1129, 1464, 1831, 2226, 2645, 3084, 3539, 4006, 4481, 4720, 4481, 4006, 3539, 3084, 2645, 2226, 1831, 1464, 1129, 830, 571, 356, 189, 74, 15] x := N16_posA1_0 x hu hv

theorem N16_A1_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 16 1 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posA1 x h0 h1)
  rw [← N16_linkA1 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N16_linkB2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 16 2 x * (2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [28, 112, 248, 432, 660, 928, 1232, 1568, 1932, 2320, 2728, 3152, 3588, 4032, 4256, 4032, 3588, 3152, 2728, 2320, 1932, 1568, 1232, 928, 660, 432, 248, 112, 28] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (16 - 2 : ℕ) = 14 by norm_num, show (2 * 16 - 2 : ℕ) = 30 by norm_num, show (16 + 2 : ℕ) = 18 by norm_num, show (16 - 2 + 16 : ℕ) = 30 by norm_num, show (2 + 16 : ℕ) = 18 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posB2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [28, 112, 248, 432, 660, 928, 1232, 1568, 1932, 2320, 2728, 3152, 3588, 4032, 4256, 4032, 3588, 3152, 2728, 2320, 1932, 1568, 1232, 928, 660, 432, 248, 112, 28] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 28 * q ^ 28 + 896 * p ^ 1 * q ^ 27 + 13856 * p ^ 2 * q ^ 26 + 137920 * p ^ 3 * q ^ 25 + 992960 * p ^ 4 * q ^ 24 + 5508608 * p ^ 5 * q ^ 23 + 24496416 * p ^ 6 * q ^ 22 + 89683776 * p ^ 7 * q ^ 21 + 275511600 * p ^ 8 * q ^ 20 + 720220800 * p ^ 9 * q ^ 19 + 1619049120 * p ^ 10 * q ^ 18 + 3154856640 * p ^ 11 * q ^ 17 + 5360759040 * p ^ 12 * q ^ 16 + 7978387200 * p ^ 13 * q ^ 15 + 10432413376 * p ^ 14 * q ^ 14 + 12007916928 * p ^ 15 * q ^ 13 + 12176691504 * p ^ 16 * q ^ 12 + 10876300032 * p ^ 17 * q ^ 11 + 8546040256 * p ^ 18 * q ^ 10 + 5892478592 * p ^ 19 * q ^ 9 + 3551067520 * p ^ 20 * q ^ 8 + 1859393536 * p ^ 21 * q ^ 7 + 838527872 * p ^ 22 * q ^ 6 + 321414912 * p ^ 23 * q ^ 5 + 102613056 * p ^ 24 * q ^ 4 + 26417664 * p ^ 25 * q ^ 3 + 5193216 * p ^ 26 * q ^ 2 + 702464 * p ^ 27 * q ^ 1 + 50176 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [28, 112, 248, 432, 660, 928, 1232, 1568, 1932, 2320, 2728, 3152, 3588, 4032, 4256, 4032, 3588, 3152, 2728, 2320, 1932, 1568, 1232, 928, 660, 432, 248, 112, 28] x =
      (fun p q : ℝ => 28 * q ^ 28 + 896 * p ^ 1 * q ^ 27 + 13856 * p ^ 2 * q ^ 26 + 137920 * p ^ 3 * q ^ 25 + 992960 * p ^ 4 * q ^ 24 + 5508608 * p ^ 5 * q ^ 23 + 24496416 * p ^ 6 * q ^ 22 + 89683776 * p ^ 7 * q ^ 21 + 275511600 * p ^ 8 * q ^ 20 + 720220800 * p ^ 9 * q ^ 19 + 1619049120 * p ^ 10 * q ^ 18 + 3154856640 * p ^ 11 * q ^ 17 + 5360759040 * p ^ 12 * q ^ 16 + 7978387200 * p ^ 13 * q ^ 15 + 10432413376 * p ^ 14 * q ^ 14 + 12007916928 * p ^ 15 * q ^ 13 + 12176691504 * p ^ 16 * q ^ 12 + 10876300032 * p ^ 17 * q ^ 11 + 8546040256 * p ^ 18 * q ^ 10 + 5892478592 * p ^ 19 * q ^ 9 + 3551067520 * p ^ 20 * q ^ 8 + 1859393536 * p ^ 21 * q ^ 7 + 838527872 * p ^ 22 * q ^ 6 + 321414912 * p ^ 23 * q ^ 5 + 102613056 * p ^ 24 * q ^ 4 + 26417664 * p ^ 25 * q ^ 3 + 5193216 * p ^ 26 * q ^ 2 + 702464 * p ^ 27 * q ^ 1 + 50176 * p ^ 28) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [28, 112, 248, 432, 660, 928, 1232, 1568, 1932, 2320, 2728, 3152, 3588, 4032, 4256, 4032, 3588, 3152, 2728, 2320, 1932, 1568, 1232, 928, 660, 432, 248, 112, 28] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posB2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [28, 112, 248, 432, 660, 928, 1232, 1568, 1932, 2320, 2728, 3152, 3588, 4032, 4256, 4032, 3588, 3152, 2728, 2320, 1932, 1568, 1232, 928, 660, 432, 248, 112, 28] x := N16_posB2_0 x hu hv

theorem N16_B2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 16 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posB2 x h0 h1)
  rw [← N16_linkB2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N16_linkA2 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 16 2 x * (2 * (1 - x) ^ 2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [2, 64, 150, 288, 474, 704, 974, 1280, 1618, 1984, 2374, 2784, 3210, 3648, 3870, 4096, 3870, 3648, 3210, 2784, 2374, 1984, 1618, 1280, 974, 704, 474, 288, 150, 64, 2] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (16 - 2 : ℕ) = 14 by norm_num, show (2 * 16 - 2 : ℕ) = 30 by norm_num, show (16 + 2 : ℕ) = 18 by norm_num, show (16 - 2 + 16 : ℕ) = 30 by norm_num, show (2 + 16 : ℕ) = 18 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posA2_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 64, 150, 288, 474, 704, 974, 1280, 1618, 1984, 2374, 2784, 3210, 3648, 3870, 4096, 3870, 3648, 3210, 2784, 2374, 1984, 1618, 1280, 974, 704, 474, 288, 150, 64, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 30 + 124 * p ^ 1 * q ^ 29 + 2876 * p ^ 2 * q ^ 28 + 38592 * p ^ 3 * q ^ 27 + 353616 * p ^ 4 * q ^ 26 + 2410592 * p ^ 5 * q ^ 25 + 12874144 * p ^ 6 * q ^ 24 + 55737536 * p ^ 7 * q ^ 23 + 200362392 * p ^ 8 * q ^ 22 + 608586576 * p ^ 9 * q ^ 21 + 1582585680 * p ^ 10 * q ^ 20 + 3558759360 * p ^ 11 * q ^ 19 + 6973565040 * p ^ 12 * q ^ 18 + 11978146080 * p ^ 13 * q ^ 17 + 18114332416 * p ^ 14 * q ^ 16 + 24195767808 * p ^ 15 * q ^ 15 + 28606221256 * p ^ 16 * q ^ 14 + 29969427312 * p ^ 17 * q ^ 13 + 27827647472 * p ^ 18 * q ^ 12 + 22883502720 * p ^ 19 * q ^ 11 + 16635480160 * p ^ 20 * q ^ 10 + 10658885952 * p ^ 21 * q ^ 9 + 5992285312 * p ^ 22 * q ^ 8 + 2936494080 * p ^ 23 * q ^ 7 + 1242485920 * p ^ 24 * q ^ 6 + 447615168 * p ^ 25 * q ^ 5 + 134433600 * p ^ 26 * q ^ 4 + 32564224 * p ^ 27 * q ^ 3 + 6019456 * p ^ 28 * q ^ 2 + 764160 * p ^ 29 * q ^ 1 + 50944 * p ^ 30 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 64, 150, 288, 474, 704, 974, 1280, 1618, 1984, 2374, 2784, 3210, 3648, 3870, 4096, 3870, 3648, 3210, 2784, 2374, 1984, 1618, 1280, 974, 704, 474, 288, 150, 64, 2] x =
      (fun p q : ℝ => 2 * q ^ 30 + 124 * p ^ 1 * q ^ 29 + 2876 * p ^ 2 * q ^ 28 + 38592 * p ^ 3 * q ^ 27 + 353616 * p ^ 4 * q ^ 26 + 2410592 * p ^ 5 * q ^ 25 + 12874144 * p ^ 6 * q ^ 24 + 55737536 * p ^ 7 * q ^ 23 + 200362392 * p ^ 8 * q ^ 22 + 608586576 * p ^ 9 * q ^ 21 + 1582585680 * p ^ 10 * q ^ 20 + 3558759360 * p ^ 11 * q ^ 19 + 6973565040 * p ^ 12 * q ^ 18 + 11978146080 * p ^ 13 * q ^ 17 + 18114332416 * p ^ 14 * q ^ 16 + 24195767808 * p ^ 15 * q ^ 15 + 28606221256 * p ^ 16 * q ^ 14 + 29969427312 * p ^ 17 * q ^ 13 + 27827647472 * p ^ 18 * q ^ 12 + 22883502720 * p ^ 19 * q ^ 11 + 16635480160 * p ^ 20 * q ^ 10 + 10658885952 * p ^ 21 * q ^ 9 + 5992285312 * p ^ 22 * q ^ 8 + 2936494080 * p ^ 23 * q ^ 7 + 1242485920 * p ^ 24 * q ^ 6 + 447615168 * p ^ 25 * q ^ 5 + 134433600 * p ^ 26 * q ^ 4 + 32564224 * p ^ 27 * q ^ 3 + 6019456 * p ^ 28 * q ^ 2 + 764160 * p ^ 29 * q ^ 1 + 50944 * p ^ 30) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 64, 150, 288, 474, 704, 974, 1280, 1618, 1984, 2374, 2784, 3210, 3648, 3870, 4096, 3870, 3648, 3210, 2784, 2374, 1984, 1618, 1280, 974, 704, 474, 288, 150, 64, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posA2 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 64, 150, 288, 474, 704, 974, 1280, 1618, 1984, 2374, 2784, 3210, 3648, 3870, 4096, 3870, 3648, 3210, 2784, 2374, 1984, 1618, 1280, 974, 704, 474, 288, 150, 64, 2] x := N16_posA2_0 x hu hv

theorem N16_A2_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < aBlock 16 2 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posA2 x h0 h1)
  rw [← N16_linkA2 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N16_linkB3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 16 3 x * (2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [39, 156, 390, 732, 1173, 1704, 2316, 3000, 3747, 4548, 5394, 6276, 7185, 7904, 8216, 7904, 7185, 6276, 5394, 4548, 3747, 3000, 2316, 1704, 1173, 732, 390, 156, 39] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (16 - 3 : ℕ) = 13 by norm_num, show (2 * 16 - 3 : ℕ) = 29 by norm_num, show (16 + 3 : ℕ) = 19 by norm_num, show (16 - 3 + 16 : ℕ) = 29 by norm_num, show (3 + 16 : ℕ) = 19 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posB3_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [39, 156, 390, 732, 1173, 1704, 2316, 3000, 3747, 4548, 5394, 6276, 7185, 7904, 8216, 7904, 7185, 6276, 5394, 4548, 3747, 3000, 2316, 1704, 1173, 732, 390, 156, 39] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 39 * q ^ 28 + 1248 * p ^ 1 * q ^ 27 + 19344 * p ^ 2 * q ^ 26 + 193392 * p ^ 3 * q ^ 25 + 1401048 * p ^ 4 * q ^ 24 + 7834176 * p ^ 5 * q ^ 23 + 35166096 * p ^ 6 * q ^ 22 + 130128336 * p ^ 7 * q ^ 21 + 404513460 * p ^ 8 * q ^ 20 + 1071102240 * p ^ 9 * q ^ 19 + 2441059920 * p ^ 10 * q ^ 18 + 4825931760 * p ^ 11 * q ^ 17 + 8325200520 * p ^ 12 * q ^ 16 + 12586129712 * p ^ 13 * q ^ 15 + 16725314528 * p ^ 14 * q ^ 14 + 19572145632 * p ^ 15 * q ^ 13 + 20184466172 * p ^ 16 * q ^ 12 + 18339566336 * p ^ 17 * q ^ 11 + 14660801376 * p ^ 18 * q ^ 10 + 10284701856 * p ^ 19 * q ^ 9 + 6305222352 * p ^ 20 * q ^ 8 + 3357443232 * p ^ 21 * q ^ 7 + 1538738240 * p ^ 22 * q ^ 6 + 598784576 * p ^ 23 * q ^ 5 + 193787568 * p ^ 24 * q ^ 4 + 50478272 * p ^ 25 * q ^ 3 + 10015616 * p ^ 26 * q ^ 2 + 1362816 * p ^ 27 * q ^ 1 + 97344 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [39, 156, 390, 732, 1173, 1704, 2316, 3000, 3747, 4548, 5394, 6276, 7185, 7904, 8216, 7904, 7185, 6276, 5394, 4548, 3747, 3000, 2316, 1704, 1173, 732, 390, 156, 39] x =
      (fun p q : ℝ => 39 * q ^ 28 + 1248 * p ^ 1 * q ^ 27 + 19344 * p ^ 2 * q ^ 26 + 193392 * p ^ 3 * q ^ 25 + 1401048 * p ^ 4 * q ^ 24 + 7834176 * p ^ 5 * q ^ 23 + 35166096 * p ^ 6 * q ^ 22 + 130128336 * p ^ 7 * q ^ 21 + 404513460 * p ^ 8 * q ^ 20 + 1071102240 * p ^ 9 * q ^ 19 + 2441059920 * p ^ 10 * q ^ 18 + 4825931760 * p ^ 11 * q ^ 17 + 8325200520 * p ^ 12 * q ^ 16 + 12586129712 * p ^ 13 * q ^ 15 + 16725314528 * p ^ 14 * q ^ 14 + 19572145632 * p ^ 15 * q ^ 13 + 20184466172 * p ^ 16 * q ^ 12 + 18339566336 * p ^ 17 * q ^ 11 + 14660801376 * p ^ 18 * q ^ 10 + 10284701856 * p ^ 19 * q ^ 9 + 6305222352 * p ^ 20 * q ^ 8 + 3357443232 * p ^ 21 * q ^ 7 + 1538738240 * p ^ 22 * q ^ 6 + 598784576 * p ^ 23 * q ^ 5 + 193787568 * p ^ 24 * q ^ 4 + 50478272 * p ^ 25 * q ^ 3 + 10015616 * p ^ 26 * q ^ 2 + 1362816 * p ^ 27 * q ^ 1 + 97344 * p ^ 28) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [39, 156, 390, 732, 1173, 1704, 2316, 3000, 3747, 4548, 5394, 6276, 7185, 7904, 8216, 7904, 7185, 6276, 5394, 4548, 3747, 3000, 2316, 1704, 1173, 732, 390, 156, 39] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posB3 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [39, 156, 390, 732, 1173, 1704, 2316, 3000, 3747, 4548, 5394, 6276, 7185, 7904, 8216, 7904, 7185, 6276, 5394, 4548, 3747, 3000, 2316, 1704, 1173, 732, 390, 156, 39] x := N16_posB3_0 x hu hv

theorem N16_B3_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 16 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posB3 x h0 h1)
  rw [← N16_linkB3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N16_linkA3 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 16 3 x * (2 * (1 - x) ^ 2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [-9, 42, 149, 260, 419, 622, 865, 1144, 1455, 1794, 2157, 2540, 2939, 3142, 3353, 3568, 3353, 3142, 2939, 2540, 2157, 1794, 1455, 1144, 865, 622, 419, 260, 149, 42, -9] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (16 - 3 : ℕ) = 13 by norm_num, show (2 * 16 - 3 : ℕ) = 29 by norm_num, show (16 + 3 : ℕ) = 19 by norm_num, show (16 - 3 + 16 : ℕ) = 29 by norm_num, show (3 + 16 : ℕ) = 19 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posA3_0 (x : ℝ) (hu : ((47 : ℝ) / 128) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-9, 42, 149, 260, 419, 622, 865, 1144, 1455, 1794, 2157, 2540, 2939, 3142, 3353, 3568, 3353, 3142, 2939, 2540, 2157, 1794, 1455, 1144, 865, 622, 419, 260, 149, 42, -9] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 30258534353326863991696891348541148990167321794658991445380348125 * q ^ 30 + 1071249991393303239590857292290868977973208105267274661464865932800 * p ^ 1 * q ^ 29 + 18277485370420656033487952433573772602895307596224233411546879950848 * p ^ 2 * q ^ 28 + 200350335487956277845928498678126768798601699891059141928886783180800 * p ^ 3 * q ^ 27 + 1586699168514353714492583945858979976649067495351264866647539149963264 * p ^ 4 * q ^ 26 + 9678441167455249946707957962618836215533966531622534137521864756428800 * p ^ 5 * q ^ 25 + 47328389863366587576676957490911317949826206259961342347587050050420736 * p ^ 6 * q ^ 24 + 190667083646839769971098262080970349665294166466728363977635933034905600 * p ^ 7 * q ^ 23 + 645279287588897200878981767111030955432323836979633347086876217152897024 * p ^ 8 * q ^ 22 + 1861350434262087141170872452767657668367145017631601958687017902538752000 * p ^ 9 * q ^ 21 + 4626830975539693237837417479455981007419721441973004987028312550344228864 * p ^ 10 * q ^ 20 + 9994559500025538353311177545650104090778783254220593887262374312384921600 * p ^ 11 * q ^ 19 + 18882495752606730891091037845530007667029697507921452652979984211195723776 * p ^ 12 * q ^ 18 + 31352525122103107548874684016552938406488703599742347493718851737092096000 * p ^ 13 * q ^ 17 + 45911896453555162484914291346947259768659301694346944662781283088707092480 * p ^ 14 * q ^ 16 + 59433204469182432584224506638600407432683688377331502080545774681194496000 * p ^ 15 * q ^ 15 + 68096235587172355605356720335678545948765549060264703886138325755410513920 * p ^ 16 * q ^ 14 + 69069267724122055470048374236029333117925624096718213675279905341623500800 * p ^ 17 * q ^ 13 + 61960983279676031741602493122282423855911394750532982136781115309730824192 * p ^ 18 * q ^ 12 + 49058572134726542215257290646048857271491956208850969964076028634595328000 * p ^ 19 * q ^ 11 + 34165223141744606576651696314467144773430638499189306483699001836435406848 * p ^ 20 * q ^ 10 + 20823049668174988576141878630222786363742196277325604696568225540027187200 * p ^ 21 * q ^ 9 + 11029451654803811029819114845937826006537290000120024085614470536764588032 * p ^ 22 * q ^ 8 + 5028708349341342711090219508507373841263657282687731656408892231699660800 * p ^ 23 * q ^ 7 + 1947950673447229514802909651609720485717041431525923219463034248515026944 * p ^ 24 * q ^ 6 + 629604315171859343858913963516904576196637980976331622439868428189696000 * p ^ 25 * q ^ 5 + 165479945752022702588770466131960524394559207659933500649395060442398720 * p ^ 26 * q ^ 4 + 34032785172739653757639294270873530031205117859735230258013464113971200 * p ^ 27 * q ^ 3 + 5146141355932339461545525084932242507859945883586178793797627441840128 * p ^ 28 * q ^ 2 + 509695036630243571557501704297139318918402434380001362672231645184000 * p ^ 29 * q ^ 1 + 24853700833779496060708654533346222027259242514529590256017390698496 * p ^ 30 := by
    intro p q hp hq; positivity
  have key : (599003433304810403471059943169868346577158542512617035467 : ℝ) * peval [-9, 42, 149, 260, 419, 622, 865, 1144, 1455, 1794, 2157, 2540, 2939, 3142, 3353, 3568, 3353, 3142, 2939, 2540, 2157, 1794, 1455, 1144, 865, 622, 419, 260, 149, 42, -9] x =
      (fun p q : ℝ => 30258534353326863991696891348541148990167321794658991445380348125 * q ^ 30 + 1071249991393303239590857292290868977973208105267274661464865932800 * p ^ 1 * q ^ 29 + 18277485370420656033487952433573772602895307596224233411546879950848 * p ^ 2 * q ^ 28 + 200350335487956277845928498678126768798601699891059141928886783180800 * p ^ 3 * q ^ 27 + 1586699168514353714492583945858979976649067495351264866647539149963264 * p ^ 4 * q ^ 26 + 9678441167455249946707957962618836215533966531622534137521864756428800 * p ^ 5 * q ^ 25 + 47328389863366587576676957490911317949826206259961342347587050050420736 * p ^ 6 * q ^ 24 + 190667083646839769971098262080970349665294166466728363977635933034905600 * p ^ 7 * q ^ 23 + 645279287588897200878981767111030955432323836979633347086876217152897024 * p ^ 8 * q ^ 22 + 1861350434262087141170872452767657668367145017631601958687017902538752000 * p ^ 9 * q ^ 21 + 4626830975539693237837417479455981007419721441973004987028312550344228864 * p ^ 10 * q ^ 20 + 9994559500025538353311177545650104090778783254220593887262374312384921600 * p ^ 11 * q ^ 19 + 18882495752606730891091037845530007667029697507921452652979984211195723776 * p ^ 12 * q ^ 18 + 31352525122103107548874684016552938406488703599742347493718851737092096000 * p ^ 13 * q ^ 17 + 45911896453555162484914291346947259768659301694346944662781283088707092480 * p ^ 14 * q ^ 16 + 59433204469182432584224506638600407432683688377331502080545774681194496000 * p ^ 15 * q ^ 15 + 68096235587172355605356720335678545948765549060264703886138325755410513920 * p ^ 16 * q ^ 14 + 69069267724122055470048374236029333117925624096718213675279905341623500800 * p ^ 17 * q ^ 13 + 61960983279676031741602493122282423855911394750532982136781115309730824192 * p ^ 18 * q ^ 12 + 49058572134726542215257290646048857271491956208850969964076028634595328000 * p ^ 19 * q ^ 11 + 34165223141744606576651696314467144773430638499189306483699001836435406848 * p ^ 20 * q ^ 10 + 20823049668174988576141878630222786363742196277325604696568225540027187200 * p ^ 21 * q ^ 9 + 11029451654803811029819114845937826006537290000120024085614470536764588032 * p ^ 22 * q ^ 8 + 5028708349341342711090219508507373841263657282687731656408892231699660800 * p ^ 23 * q ^ 7 + 1947950673447229514802909651609720485717041431525923219463034248515026944 * p ^ 24 * q ^ 6 + 629604315171859343858913963516904576196637980976331622439868428189696000 * p ^ 25 * q ^ 5 + 165479945752022702588770466131960524394559207659933500649395060442398720 * p ^ 26 * q ^ 4 + 34032785172739653757639294270873530031205117859735230258013464113971200 * p ^ 27 * q ^ 3 + 5146141355932339461545525084932242507859945883586178793797627441840128 * p ^ 28 * q ^ 2 + 509695036630243571557501704297139318918402434380001362672231645184000 * p ^ 29 * q ^ 1 + 24853700833779496060708654533346222027259242514529590256017390698496 * p ^ 30) (x - ((47 : ℝ) / 128)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (599003433304810403471059943169868346577158542512617035467 : ℝ) * peval [-9, 42, 149, 260, 419, 622, 865, 1144, 1455, 1794, 2157, 2540, 2939, 3142, 3353, 3568, 3353, 3142, 2939, 2540, 2157, 1794, 1455, 1144, 865, 622, 419, 260, 149, 42, -9] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posA3 (x : ℝ) (hu : ((47 : ℝ) / 128) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-9, 42, 149, 260, 419, 622, 865, 1144, 1455, 1794, 2157, 2540, 2939, 3142, 3353, 3568, 3353, 3142, 2939, 2540, 2157, 1794, 1455, 1144, 865, 622, 419, 260, 149, 42, -9] x := N16_posA3_0 x hu hv

theorem N16_A3_pos (x : ℝ) (h0 : 0 < x) (hlo : ((47 : ℝ) / 128) < x) (h1 : x < 1) :
    0 < aBlock 16 3 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posA3 x hlo h1)
  rw [← N16_linkA3 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N16_linkB4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 16 4 x * (2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [48, 192, 480, 960, 1616, 2432, 3392, 4480, 5680, 6976, 8352, 9792, 11088, 12032, 12416, 12032, 11088, 9792, 8352, 6976, 5680, 4480, 3392, 2432, 1616, 960, 480, 192, 48] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (16 - 4 : ℕ) = 12 by norm_num, show (2 * 16 - 4 : ℕ) = 28 by norm_num, show (16 + 4 : ℕ) = 20 by norm_num, show (16 - 4 + 16 : ℕ) = 28 by norm_num, show (4 + 16 : ℕ) = 20 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posB4_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [48, 192, 480, 960, 1616, 2432, 3392, 4480, 5680, 6976, 8352, 9792, 11088, 12032, 12416, 12032, 11088, 9792, 8352, 6976, 5680, 4480, 3392, 2432, 1616, 960, 480, 192, 48] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 48 * q ^ 28 + 1536 * p ^ 1 * q ^ 27 + 23808 * p ^ 2 * q ^ 26 + 238080 * p ^ 3 * q ^ 25 + 1726016 * p ^ 4 * q ^ 24 + 9664256 * p ^ 5 * q ^ 23 + 43473024 * p ^ 6 * q ^ 22 + 161351424 * p ^ 7 * q ^ 21 + 503568000 * p ^ 8 * q ^ 20 + 1340052480 * p ^ 9 * q ^ 19 + 3072476160 * p ^ 10 * q ^ 18 + 6117396480 * p ^ 11 * q ^ 17 + 10639137408 * p ^ 12 * q ^ 16 + 16231887360 * p ^ 13 * q ^ 15 + 21788993280 * p ^ 14 * q ^ 14 + 25780124160 * p ^ 15 * q ^ 13 + 26903936576 * p ^ 16 * q ^ 12 + 24755422208 * p ^ 17 * q ^ 11 + 20054054400 * p ^ 18 * q ^ 10 + 14262983680 * p ^ 19 * q ^ 9 + 8867634688 * p ^ 20 * q ^ 8 + 4788393984 * p ^ 21 * q ^ 7 + 2224553984 * p ^ 22 * q ^ 6 + 876769280 * p ^ 23 * q ^ 5 + 287032320 * p ^ 24 * q ^ 4 + 75505664 * p ^ 25 * q ^ 3 + 15097856 * p ^ 26 * q ^ 2 + 2064384 * p ^ 27 * q ^ 1 + 147456 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [48, 192, 480, 960, 1616, 2432, 3392, 4480, 5680, 6976, 8352, 9792, 11088, 12032, 12416, 12032, 11088, 9792, 8352, 6976, 5680, 4480, 3392, 2432, 1616, 960, 480, 192, 48] x =
      (fun p q : ℝ => 48 * q ^ 28 + 1536 * p ^ 1 * q ^ 27 + 23808 * p ^ 2 * q ^ 26 + 238080 * p ^ 3 * q ^ 25 + 1726016 * p ^ 4 * q ^ 24 + 9664256 * p ^ 5 * q ^ 23 + 43473024 * p ^ 6 * q ^ 22 + 161351424 * p ^ 7 * q ^ 21 + 503568000 * p ^ 8 * q ^ 20 + 1340052480 * p ^ 9 * q ^ 19 + 3072476160 * p ^ 10 * q ^ 18 + 6117396480 * p ^ 11 * q ^ 17 + 10639137408 * p ^ 12 * q ^ 16 + 16231887360 * p ^ 13 * q ^ 15 + 21788993280 * p ^ 14 * q ^ 14 + 25780124160 * p ^ 15 * q ^ 13 + 26903936576 * p ^ 16 * q ^ 12 + 24755422208 * p ^ 17 * q ^ 11 + 20054054400 * p ^ 18 * q ^ 10 + 14262983680 * p ^ 19 * q ^ 9 + 8867634688 * p ^ 20 * q ^ 8 + 4788393984 * p ^ 21 * q ^ 7 + 2224553984 * p ^ 22 * q ^ 6 + 876769280 * p ^ 23 * q ^ 5 + 287032320 * p ^ 24 * q ^ 4 + 75505664 * p ^ 25 * q ^ 3 + 15097856 * p ^ 26 * q ^ 2 + 2064384 * p ^ 27 * q ^ 1 + 147456 * p ^ 28) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [48, 192, 480, 960, 1616, 2432, 3392, 4480, 5680, 6976, 8352, 9792, 11088, 12032, 12416, 12032, 11088, 9792, 8352, 6976, 5680, 4480, 3392, 2432, 1616, 960, 480, 192, 48] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posB4 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [48, 192, 480, 960, 1616, 2432, 3392, 4480, 5680, 6976, 8352, 9792, 11088, 12032, 12416, 12032, 11088, 9792, 8352, 6976, 5680, 4480, 3392, 2432, 1616, 960, 480, 192, 48] x := N16_posB4_0 x hu hv

theorem N16_B4_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 16 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posB4 x h0 h1)
  rw [← N16_linkB4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N16_linkA4 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 16 4 x * (2 * (1 - x) ^ 2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [-18, 24, 122, 272, 406, 584, 802, 1056, 1342, 1656, 1994, 2352, 2534, 2728, 2930, 3136, 2930, 2728, 2534, 2352, 1994, 1656, 1342, 1056, 802, 584, 406, 272, 122, 24, -18] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (16 - 4 : ℕ) = 12 by norm_num, show (2 * 16 - 4 : ℕ) = 28 by norm_num, show (16 + 4 : ℕ) = 20 by norm_num, show (16 - 4 + 16 : ℕ) = 28 by norm_num, show (4 + 16 : ℕ) = 20 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posA4_0 (x : ℝ) (hu : ((47 : ℝ) / 128) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-18, 24, 122, 272, 406, 584, 802, 1056, 1342, 1656, 1994, 2352, 2534, 2728, 2930, 3136, 2930, 2728, 2534, 2352, 1994, 1656, 1342, 1056, 802, 584, 406, 272, 122, 24, -18] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 19591558032085489297748451414214794519935486350967134566694198202 * q ^ 30 + 735731434919545127319767729046267600427414585430251815342791308800 * p ^ 1 * q ^ 29 + 13167671939221670962265141905356981685889540860404572822228739883008 * p ^ 2 * q ^ 28 + 150147045508870878320593911593486365928753445129184022000762145996800 * p ^ 3 * q ^ 27 + 1228984240287157297273605358793097193853407292508298689824947955564544 * p ^ 4 * q ^ 26 + 7708128335722203114651626366449400878932925390355070065889101689651200 * p ^ 5 * q ^ 25 + 38595828069869118825881121058122856763224906689553867554391083517476864 * p ^ 6 * q ^ 24 + 158657870837014432387122857985313457753011157698664733801119093882880000 * p ^ 7 * q ^ 23 + 546294699119816593221684840083704634196070822786865466838745837941555200 * p ^ 8 * q ^ 22 + 1599216928282432673157129822199143956334402197570037463706446526939136000 * p ^ 9 * q ^ 21 + 4025418410432956747264697261320777773148120809376718796349641989661130752 * p ^ 10 * q ^ 20 + 8788204157496276413884138905954384614457383411216163487546626911974195200 * p ^ 11 * q ^ 19 + 16751627267865181946486044363300917489144226599584329195435152367285174272 * p ^ 12 * q ^ 18 + 28019355186425989371106464684280041461075062169641338792016634251208294400 * p ^ 13 * q ^ 17 + 41275126112707031357664256576304021396867883827931195930022106328069046272 * p ^ 14 * q ^ 16 + 53680094618105932116937644721641549482765625367182980539947492725017804800 * p ^ 15 * q ^ 15 + 61719065895453000480700722005884057084143480116757096192968523761420599296 * p ^ 16 * q ^ 14 + 62752147083999749686537627159585609624823326276024240926016670533681152000 * p ^ 17 * q ^ 13 + 56375168472918013599297984962495891780813578799671660469736562722582036480 * p ^ 18 * q ^ 12 + 44661398932805864438656307889884199330822870139553773853883587566239744000 * p ^ 19 * q ^ 11 + 31096967018942817103653549638030256054566323184999776467922142097209032704 * p ^ 20 * q ^ 10 + 18937107975791458922698068212904189311107927656655478468351648759965286400 * p ^ 21 * q ^ 9 + 10016917922151519886850604746133610990237966630884810723764196905827958784 * p ^ 22 * q ^ 8 + 4559160764359401680674429003096952456221558976578127871027030550747545600 * p ^ 23 * q ^ 7 + 1762624408321126422437593830238530415289644561859822057707278791457374208 * p ^ 24 * q ^ 6 + 568564910296277442963693779805778736691954810071527307471254715577139200 * p ^ 25 * q ^ 5 + 149155536578763900677842472758727560643648593160941381782114654952095744 * p ^ 26 * q ^ 4 + 30626426778676588256122284366625764532360999457908568362514752785612800 * p ^ 27 * q ^ 3 + 4625657665088302058112050857707391460825380292835330491718920759672832 * p ^ 28 * q ^ 2 + 457861643074625581229620175046582778011446254612543596976750460928000 * p ^ 29 * q ^ 1 + 22326205833734123579958621868938131651605760224916411585913927237632 * p ^ 30 := by
    intro p q hp hq; positivity
  have key : (599003433304810403471059943169868346577158542512617035467 : ℝ) * peval [-18, 24, 122, 272, 406, 584, 802, 1056, 1342, 1656, 1994, 2352, 2534, 2728, 2930, 3136, 2930, 2728, 2534, 2352, 1994, 1656, 1342, 1056, 802, 584, 406, 272, 122, 24, -18] x =
      (fun p q : ℝ => 19591558032085489297748451414214794519935486350967134566694198202 * q ^ 30 + 735731434919545127319767729046267600427414585430251815342791308800 * p ^ 1 * q ^ 29 + 13167671939221670962265141905356981685889540860404572822228739883008 * p ^ 2 * q ^ 28 + 150147045508870878320593911593486365928753445129184022000762145996800 * p ^ 3 * q ^ 27 + 1228984240287157297273605358793097193853407292508298689824947955564544 * p ^ 4 * q ^ 26 + 7708128335722203114651626366449400878932925390355070065889101689651200 * p ^ 5 * q ^ 25 + 38595828069869118825881121058122856763224906689553867554391083517476864 * p ^ 6 * q ^ 24 + 158657870837014432387122857985313457753011157698664733801119093882880000 * p ^ 7 * q ^ 23 + 546294699119816593221684840083704634196070822786865466838745837941555200 * p ^ 8 * q ^ 22 + 1599216928282432673157129822199143956334402197570037463706446526939136000 * p ^ 9 * q ^ 21 + 4025418410432956747264697261320777773148120809376718796349641989661130752 * p ^ 10 * q ^ 20 + 8788204157496276413884138905954384614457383411216163487546626911974195200 * p ^ 11 * q ^ 19 + 16751627267865181946486044363300917489144226599584329195435152367285174272 * p ^ 12 * q ^ 18 + 28019355186425989371106464684280041461075062169641338792016634251208294400 * p ^ 13 * q ^ 17 + 41275126112707031357664256576304021396867883827931195930022106328069046272 * p ^ 14 * q ^ 16 + 53680094618105932116937644721641549482765625367182980539947492725017804800 * p ^ 15 * q ^ 15 + 61719065895453000480700722005884057084143480116757096192968523761420599296 * p ^ 16 * q ^ 14 + 62752147083999749686537627159585609624823326276024240926016670533681152000 * p ^ 17 * q ^ 13 + 56375168472918013599297984962495891780813578799671660469736562722582036480 * p ^ 18 * q ^ 12 + 44661398932805864438656307889884199330822870139553773853883587566239744000 * p ^ 19 * q ^ 11 + 31096967018942817103653549638030256054566323184999776467922142097209032704 * p ^ 20 * q ^ 10 + 18937107975791458922698068212904189311107927656655478468351648759965286400 * p ^ 21 * q ^ 9 + 10016917922151519886850604746133610990237966630884810723764196905827958784 * p ^ 22 * q ^ 8 + 4559160764359401680674429003096952456221558976578127871027030550747545600 * p ^ 23 * q ^ 7 + 1762624408321126422437593830238530415289644561859822057707278791457374208 * p ^ 24 * q ^ 6 + 568564910296277442963693779805778736691954810071527307471254715577139200 * p ^ 25 * q ^ 5 + 149155536578763900677842472758727560643648593160941381782114654952095744 * p ^ 26 * q ^ 4 + 30626426778676588256122284366625764532360999457908568362514752785612800 * p ^ 27 * q ^ 3 + 4625657665088302058112050857707391460825380292835330491718920759672832 * p ^ 28 * q ^ 2 + 457861643074625581229620175046582778011446254612543596976750460928000 * p ^ 29 * q ^ 1 + 22326205833734123579958621868938131651605760224916411585913927237632 * p ^ 30) (x - ((47 : ℝ) / 128)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (599003433304810403471059943169868346577158542512617035467 : ℝ) * peval [-18, 24, 122, 272, 406, 584, 802, 1056, 1342, 1656, 1994, 2352, 2534, 2728, 2930, 3136, 2930, 2728, 2534, 2352, 1994, 1656, 1342, 1056, 802, 584, 406, 272, 122, 24, -18] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posA4 (x : ℝ) (hu : ((47 : ℝ) / 128) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-18, 24, 122, 272, 406, 584, 802, 1056, 1342, 1656, 1994, 2352, 2534, 2728, 2930, 3136, 2930, 2728, 2534, 2352, 1994, 1656, 1342, 1056, 802, 584, 406, 272, 122, 24, -18] x := N16_posA4_0 x hu hv

theorem N16_A4_pos (x : ℝ) (h0 : 0 < x) (hlo : ((47 : ℝ) / 128) < x) (h1 : x < 1) :
    0 < aBlock 16 4 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posA4 x hlo h1)
  rw [← N16_linkA4 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N16_linkB5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 16 5 x * (2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [55, 220, 550, 1100, 1925, 3000, 4300, 5800, 7475, 9300, 11250, 13124, 14721, 15840, 16280, 15840, 14721, 13124, 11250, 9300, 7475, 5800, 4300, 3000, 1925, 1100, 550, 220, 55] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (16 - 5 : ℕ) = 11 by norm_num, show (2 * 16 - 5 : ℕ) = 27 by norm_num, show (16 + 5 : ℕ) = 21 by norm_num, show (16 - 5 + 16 : ℕ) = 27 by norm_num, show (5 + 16 : ℕ) = 21 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posB5_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [55, 220, 550, 1100, 1925, 3000, 4300, 5800, 7475, 9300, 11250, 13124, 14721, 15840, 16280, 15840, 14721, 13124, 11250, 9300, 7475, 5800, 4300, 3000, 1925, 1100, 550, 220, 55] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 55 * q ^ 28 + 1760 * p ^ 1 * q ^ 27 + 27280 * p ^ 2 * q ^ 26 + 272800 * p ^ 3 * q ^ 25 + 1977800 * p ^ 4 * q ^ 24 + 11075600 * p ^ 5 * q ^ 23 + 49838400 * p ^ 6 * q ^ 22 + 185094000 * p ^ 7 * q ^ 21 + 578272500 * p ^ 8 * q ^ 20 + 1541280000 * p ^ 9 * q ^ 19 + 3541714800 * p ^ 10 * q ^ 18 + 7072665424 * p ^ 11 * q ^ 17 + 12347560104 * p ^ 12 * q ^ 16 + 18927905040 * p ^ 13 * q ^ 15 + 25553777920 * p ^ 14 * q ^ 14 + 30438793440 * p ^ 15 * q ^ 13 + 32012808828 * p ^ 16 * q ^ 12 + 29714496768 * p ^ 17 * q ^ 11 + 24303876960 * p ^ 18 * q ^ 10 + 17465294880 * p ^ 19 * q ^ 9 + 10976897360 * p ^ 20 * q ^ 8 + 5992992544 * p ^ 21 * q ^ 7 + 2814409664 * p ^ 22 * q ^ 6 + 1120595520 * p ^ 23 * q ^ 5 + 370230960 * p ^ 24 * q ^ 4 + 98155200 * p ^ 25 * q ^ 3 + 19747200 * p ^ 26 * q ^ 2 + 2710400 * p ^ 27 * q ^ 1 + 193600 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [55, 220, 550, 1100, 1925, 3000, 4300, 5800, 7475, 9300, 11250, 13124, 14721, 15840, 16280, 15840, 14721, 13124, 11250, 9300, 7475, 5800, 4300, 3000, 1925, 1100, 550, 220, 55] x =
      (fun p q : ℝ => 55 * q ^ 28 + 1760 * p ^ 1 * q ^ 27 + 27280 * p ^ 2 * q ^ 26 + 272800 * p ^ 3 * q ^ 25 + 1977800 * p ^ 4 * q ^ 24 + 11075600 * p ^ 5 * q ^ 23 + 49838400 * p ^ 6 * q ^ 22 + 185094000 * p ^ 7 * q ^ 21 + 578272500 * p ^ 8 * q ^ 20 + 1541280000 * p ^ 9 * q ^ 19 + 3541714800 * p ^ 10 * q ^ 18 + 7072665424 * p ^ 11 * q ^ 17 + 12347560104 * p ^ 12 * q ^ 16 + 18927905040 * p ^ 13 * q ^ 15 + 25553777920 * p ^ 14 * q ^ 14 + 30438793440 * p ^ 15 * q ^ 13 + 32012808828 * p ^ 16 * q ^ 12 + 29714496768 * p ^ 17 * q ^ 11 + 24303876960 * p ^ 18 * q ^ 10 + 17465294880 * p ^ 19 * q ^ 9 + 10976897360 * p ^ 20 * q ^ 8 + 5992992544 * p ^ 21 * q ^ 7 + 2814409664 * p ^ 22 * q ^ 6 + 1120595520 * p ^ 23 * q ^ 5 + 370230960 * p ^ 24 * q ^ 4 + 98155200 * p ^ 25 * q ^ 3 + 19747200 * p ^ 26 * q ^ 2 + 2710400 * p ^ 27 * q ^ 1 + 193600 * p ^ 28) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [55, 220, 550, 1100, 1925, 3000, 4300, 5800, 7475, 9300, 11250, 13124, 14721, 15840, 16280, 15840, 14721, 13124, 11250, 9300, 7475, 5800, 4300, 3000, 1925, 1100, 550, 220, 55] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posB5 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [55, 220, 550, 1100, 1925, 3000, 4300, 5800, 7475, 9300, 11250, 13124, 14721, 15840, 16280, 15840, 14721, 13124, 11250, 9300, 7475, 5800, 4300, 3000, 1925, 1100, 550, 220, 55] x := N16_posB5_0 x hu hv

theorem N16_B5_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 16 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posB5 x h0 h1)
  rw [← N16_linkB5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N16_linkA5 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 16 5 x * (2 * (1 - x) ^ 2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [-25, 10, 101, 244, 435, 590, 785, 1016, 1279, 1570, 1885, 2044, 2219, 2406, 2601, 2800, 2601, 2406, 2219, 2044, 1885, 1570, 1279, 1016, 785, 590, 435, 244, 101, 10, -25] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (16 - 5 : ℕ) = 11 by norm_num, show (2 * 16 - 5 : ℕ) = 27 by norm_num, show (16 + 5 : ℕ) = 21 by norm_num, show (16 - 5 + 16 : ℕ) = 27 by norm_num, show (5 + 16 : ℕ) = 21 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posA5_0 (x : ℝ) (hu : ((47 : ℝ) / 128) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-25, 10, 101, 244, 435, 590, 785, 1016, 1279, 1570, 1885, 2044, 2219, 2406, 2601, 2800, 2601, 2406, 2219, 2044, 1885, 1570, 1279, 1016, 785, 590, 435, 244, 101, 10, -25] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 32590106432366307653812699477249084578263009745389828892134794375 * q ^ 30 + 1383209987055804446909068767821928810075613493439998056973669798400 * p ^ 1 * q ^ 29 + 26966069309113839878856795254594317071479880345544481907057346412544 * p ^ 2 * q ^ 28 + 327590927602793872993365677421240966350724352060216043711438730035200 * p ^ 3 * q ^ 27 + 2815587720710904984567944043136919790781636266188218987838945693270016 * p ^ 4 * q ^ 26 + 18357011287092384917029059630045867401193346164906452670364017229824000 * p ^ 5 * q ^ 25 + 94848939494591127948846941389414972065087054217504070704333017348833280 * p ^ 6 * q ^ 24 + 400103610552342206312347645251479518721609826588108386316375696605184000 * p ^ 7 * q ^ 23 + 1407528799719568001472148776909893095060366929572160188285457119167119360 * p ^ 8 * q ^ 22 + 4194895145113606853929552944729026394702713597026530092968435708329984000 * p ^ 9 * q ^ 21 + 10718444491923856403243351484285551399647759116200901034204925064105689088 * p ^ 10 * q ^ 20 + 23694138149249848707658598070370050591402505197811864169969496780020121600 * p ^ 11 * q ^ 19 + 45632248902730168035325907782290947808583853274757148061613234613466955776 * p ^ 12 * q ^ 18 + 76967668369252635129494212906635810894542819185410276801549513367171891200 * p ^ 13 * q ^ 17 + 114134994452998436582317723659009729341097176198219073202929328974364409856 * p ^ 14 * q ^ 16 + 149188915242235231866054916434547373457088987241798853964697072195076096000 * p ^ 15 * q ^ 15 + 172149533333913654228160809236949207335329739385341013838010892022626385920 * p ^ 16 * q ^ 14 + 175428853742431558905786424603184224584544267537166787415565716936458240000 * p ^ 17 * q ^ 13 + 157768339329834660231374078500804316717774136160853791587012014723668377600 * p ^ 18 * q ^ 12 + 124982820651831854824862364062415646316131468779804512863949758018879488000 * p ^ 19 * q ^ 11 + 86937249391757831214432940185367693791324493566685822965582355488370065408 * p ^ 20 * q ^ 10 + 52846901503595898730423304600359071378499788414808600618405070382091468800 * p ^ 21 * q ^ 9 + 27885527264216505472140805623052441514294574265208752495272173017736675328 * p ^ 22 * q ^ 8 + 12655267339862644509075636282033816157346458926435632734531278410140876800 * p ^ 23 * q ^ 7 + 4877324439147995124419208959875064595784345939942064526743456883696205824 * p ^ 24 * q ^ 6 + 1568302040068231952999853789179102766891416099680631146732713154707456000 * p ^ 25 * q ^ 5 + 410215690356483219089868107183682364675135592462395895131103861586001920 * p ^ 26 * q ^ 4 + 84021518008241998173315229213909157670633117693493768681400534302720000 * p ^ 27 * q ^ 3 + 12667234367949907922360007580617145570852300881968655025872114076876800 * p ^ 28 * q ^ 2 + 1252640344260768099590470290221783071918107677713562670974128619520000 * p ^ 29 * q ^ 1 + 61081129167763168284792456056528850744959155332318484527500366970880 * p ^ 30 := by
    intro p q hp hq; positivity
  have key : (1797010299914431210413179829509605039731475627537851106401 : ℝ) * peval [-25, 10, 101, 244, 435, 590, 785, 1016, 1279, 1570, 1885, 2044, 2219, 2406, 2601, 2800, 2601, 2406, 2219, 2044, 1885, 1570, 1279, 1016, 785, 590, 435, 244, 101, 10, -25] x =
      (fun p q : ℝ => 32590106432366307653812699477249084578263009745389828892134794375 * q ^ 30 + 1383209987055804446909068767821928810075613493439998056973669798400 * p ^ 1 * q ^ 29 + 26966069309113839878856795254594317071479880345544481907057346412544 * p ^ 2 * q ^ 28 + 327590927602793872993365677421240966350724352060216043711438730035200 * p ^ 3 * q ^ 27 + 2815587720710904984567944043136919790781636266188218987838945693270016 * p ^ 4 * q ^ 26 + 18357011287092384917029059630045867401193346164906452670364017229824000 * p ^ 5 * q ^ 25 + 94848939494591127948846941389414972065087054217504070704333017348833280 * p ^ 6 * q ^ 24 + 400103610552342206312347645251479518721609826588108386316375696605184000 * p ^ 7 * q ^ 23 + 1407528799719568001472148776909893095060366929572160188285457119167119360 * p ^ 8 * q ^ 22 + 4194895145113606853929552944729026394702713597026530092968435708329984000 * p ^ 9 * q ^ 21 + 10718444491923856403243351484285551399647759116200901034204925064105689088 * p ^ 10 * q ^ 20 + 23694138149249848707658598070370050591402505197811864169969496780020121600 * p ^ 11 * q ^ 19 + 45632248902730168035325907782290947808583853274757148061613234613466955776 * p ^ 12 * q ^ 18 + 76967668369252635129494212906635810894542819185410276801549513367171891200 * p ^ 13 * q ^ 17 + 114134994452998436582317723659009729341097176198219073202929328974364409856 * p ^ 14 * q ^ 16 + 149188915242235231866054916434547373457088987241798853964697072195076096000 * p ^ 15 * q ^ 15 + 172149533333913654228160809236949207335329739385341013838010892022626385920 * p ^ 16 * q ^ 14 + 175428853742431558905786424603184224584544267537166787415565716936458240000 * p ^ 17 * q ^ 13 + 157768339329834660231374078500804316717774136160853791587012014723668377600 * p ^ 18 * q ^ 12 + 124982820651831854824862364062415646316131468779804512863949758018879488000 * p ^ 19 * q ^ 11 + 86937249391757831214432940185367693791324493566685822965582355488370065408 * p ^ 20 * q ^ 10 + 52846901503595898730423304600359071378499788414808600618405070382091468800 * p ^ 21 * q ^ 9 + 27885527264216505472140805623052441514294574265208752495272173017736675328 * p ^ 22 * q ^ 8 + 12655267339862644509075636282033816157346458926435632734531278410140876800 * p ^ 23 * q ^ 7 + 4877324439147995124419208959875064595784345939942064526743456883696205824 * p ^ 24 * q ^ 6 + 1568302040068231952999853789179102766891416099680631146732713154707456000 * p ^ 25 * q ^ 5 + 410215690356483219089868107183682364675135592462395895131103861586001920 * p ^ 26 * q ^ 4 + 84021518008241998173315229213909157670633117693493768681400534302720000 * p ^ 27 * q ^ 3 + 12667234367949907922360007580617145570852300881968655025872114076876800 * p ^ 28 * q ^ 2 + 1252640344260768099590470290221783071918107677713562670974128619520000 * p ^ 29 * q ^ 1 + 61081129167763168284792456056528850744959155332318484527500366970880 * p ^ 30) (x - ((47 : ℝ) / 128)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1797010299914431210413179829509605039731475627537851106401 : ℝ) * peval [-25, 10, 101, 244, 435, 590, 785, 1016, 1279, 1570, 1885, 2044, 2219, 2406, 2601, 2800, 2601, 2406, 2219, 2044, 1885, 1570, 1279, 1016, 785, 590, 435, 244, 101, 10, -25] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posA5 (x : ℝ) (hu : ((47 : ℝ) / 128) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-25, 10, 101, 244, 435, 590, 785, 1016, 1279, 1570, 1885, 2044, 2219, 2406, 2601, 2800, 2601, 2406, 2219, 2044, 1885, 1570, 1279, 1016, 785, 590, 435, 244, 101, 10, -25] x := N16_posA5_0 x hu hv

theorem N16_A5_pos (x : ℝ) (h0 : 0 < x) (hlo : ((47 : ℝ) / 128) < x) (h1 : x < 1) :
    0 < aBlock 16 5 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posA5 x hlo h1)
  rw [← N16_linkA5 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N16_linkB6 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 16 6 x * (2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [60, 240, 600, 1200, 2100, 3360, 4944, 6816, 8940, 11280, 13640, 15824, 17636, 18880, 19360, 18880, 17636, 15824, 13640, 11280, 8940, 6816, 4944, 3360, 2100, 1200, 600, 240, 60] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (16 - 6 : ℕ) = 10 by norm_num, show (2 * 16 - 6 : ℕ) = 26 by norm_num, show (16 + 6 : ℕ) = 22 by norm_num, show (16 - 6 + 16 : ℕ) = 26 by norm_num, show (6 + 16 : ℕ) = 22 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posB6_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [60, 240, 600, 1200, 2100, 3360, 4944, 6816, 8940, 11280, 13640, 15824, 17636, 18880, 19360, 18880, 17636, 15824, 13640, 11280, 8940, 6816, 4944, 3360, 2100, 1200, 600, 240, 60] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 60 * q ^ 28 + 1920 * p ^ 1 * q ^ 27 + 29760 * p ^ 2 * q ^ 26 + 297600 * p ^ 3 * q ^ 25 + 2157600 * p ^ 4 * q ^ 24 + 12082560 * p ^ 5 * q ^ 23 + 54371424 * p ^ 6 * q ^ 22 + 201948864 * p ^ 7 * q ^ 21 + 631066800 * p ^ 8 * q ^ 20 + 1682678400 * p ^ 9 * q ^ 19 + 3869299040 * p ^ 10 * q ^ 18 + 7735150400 * p ^ 11 * q ^ 17 + 13525431360 * p ^ 12 * q ^ 16 + 20779020800 * p ^ 13 * q ^ 15 + 28134987200 * p ^ 14 * q ^ 14 + 33639192960 * p ^ 15 * q ^ 13 + 35542959152 * p ^ 16 * q ^ 12 + 33174371072 * p ^ 17 * q ^ 11 + 27307880640 * p ^ 18 * q ^ 10 + 19764590720 * p ^ 19 * q ^ 9 + 12517681280 * p ^ 20 * q ^ 8 + 6888600576 * p ^ 21 * q ^ 7 + 3260483456 * p ^ 22 * q ^ 6 + 1307866880 * p ^ 23 * q ^ 5 + 434994240 * p ^ 24 * q ^ 4 + 115980800 * p ^ 25 * q ^ 3 + 23436800 * p ^ 26 * q ^ 2 + 3225600 * p ^ 27 * q ^ 1 + 230400 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [60, 240, 600, 1200, 2100, 3360, 4944, 6816, 8940, 11280, 13640, 15824, 17636, 18880, 19360, 18880, 17636, 15824, 13640, 11280, 8940, 6816, 4944, 3360, 2100, 1200, 600, 240, 60] x =
      (fun p q : ℝ => 60 * q ^ 28 + 1920 * p ^ 1 * q ^ 27 + 29760 * p ^ 2 * q ^ 26 + 297600 * p ^ 3 * q ^ 25 + 2157600 * p ^ 4 * q ^ 24 + 12082560 * p ^ 5 * q ^ 23 + 54371424 * p ^ 6 * q ^ 22 + 201948864 * p ^ 7 * q ^ 21 + 631066800 * p ^ 8 * q ^ 20 + 1682678400 * p ^ 9 * q ^ 19 + 3869299040 * p ^ 10 * q ^ 18 + 7735150400 * p ^ 11 * q ^ 17 + 13525431360 * p ^ 12 * q ^ 16 + 20779020800 * p ^ 13 * q ^ 15 + 28134987200 * p ^ 14 * q ^ 14 + 33639192960 * p ^ 15 * q ^ 13 + 35542959152 * p ^ 16 * q ^ 12 + 33174371072 * p ^ 17 * q ^ 11 + 27307880640 * p ^ 18 * q ^ 10 + 19764590720 * p ^ 19 * q ^ 9 + 12517681280 * p ^ 20 * q ^ 8 + 6888600576 * p ^ 21 * q ^ 7 + 3260483456 * p ^ 22 * q ^ 6 + 1307866880 * p ^ 23 * q ^ 5 + 434994240 * p ^ 24 * q ^ 4 + 115980800 * p ^ 25 * q ^ 3 + 23436800 * p ^ 26 * q ^ 2 + 3225600 * p ^ 27 * q ^ 1 + 230400 * p ^ 28) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [60, 240, 600, 1200, 2100, 3360, 4944, 6816, 8940, 11280, 13640, 15824, 17636, 18880, 19360, 18880, 17636, 15824, 13640, 11280, 8940, 6816, 4944, 3360, 2100, 1200, 600, 240, 60] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posB6 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [60, 240, 600, 1200, 2100, 3360, 4944, 6816, 8940, 11280, 13640, 15824, 17636, 18880, 19360, 18880, 17636, 15824, 13640, 11280, 8940, 6816, 4944, 3360, 2100, 1200, 600, 240, 60] x := N16_posB6_0 x hu hv

theorem N16_B6_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 16 6 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posB6 x h0 h1)
  rw [← N16_linkB6 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N16_linkA6 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 16 6 x * (2 * (1 - x) ^ 2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [-30, 0, 86, 224, 410, 640, 814, 1024, 1266, 1536, 1670, 1824, 1994, 2176, 2366, 2560, 2366, 2176, 1994, 1824, 1670, 1536, 1266, 1024, 814, 640, 410, 224, 86, 0, -30] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (16 - 6 : ℕ) = 10 by norm_num, show (2 * 16 - 6 : ℕ) = 26 by norm_num, show (16 + 6 : ℕ) = 22 by norm_num, show (16 - 6 + 16 : ℕ) = 26 by norm_num, show (6 + 16 : ℕ) = 22 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posA6_0 (x : ℝ) (hu : ((47 : ℝ) / 128) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-30, 0, 86, 224, 410, 640, 814, 1024, 1266, 1536, 1670, 1824, 1994, 2176, 2366, 2560, 2366, 2176, 1994, 1824, 1670, 1536, 1266, 1024, 814, 640, 410, 224, 86, 0, -30] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 490904797029579978512147356365522781290084128559636240281795078 * q ^ 30 + 28633263796909548690183871813472731253967057260581687033880102400 * p ^ 1 * q ^ 29 + 654193099988220707147549955484109599082831403741012643772187869184 * p ^ 2 * q ^ 28 + 8752135977569476613321476822618962674229182246882883001971743129600 * p ^ 3 * q ^ 27 + 80293433541610536874925265380914300042499435221712636120180056915968 * p ^ 4 * q ^ 26 + 548783114353123750614498063841255630213158502893991958880574242816000 * p ^ 5 * q ^ 25 + 2938619161843024090715324918474020141441752733186898469281945113067520 * p ^ 6 * q ^ 24 + 12747155574543552897524335135622209422027371068883084637425281007616000 * p ^ 7 * q ^ 23 + 45856443144759702888928538938630199664370165226266593451906175745392640 * p ^ 8 * q ^ 22 + 139168793643074675082300045568830841968946153963250327252763420917760000 * p ^ 9 * q ^ 21 + 360912975006446280642783455006275933241039517477404828644499614146232320 * p ^ 10 * q ^ 20 + 807613022759769031683085203739135582767026862931964481857959308243763200 * p ^ 11 * q ^ 19 + 1570923772328951049559541982997545249366034956162668281773247158395338752 * p ^ 12 * q ^ 18 + 2671007554024262514501451604068204857934301750956149862677320713214361600 * p ^ 13 * q ^ 17 + 3985927816059353150976356863618934135626536595482897642826701708101419008 * p ^ 14 * q ^ 16 + 5235101242466759358392491821315787871385647513483036701332610882732032000 * p ^ 15 * q ^ 15 + 6061291791157119221641792893079025292959729510737668096361611332828528640 * p ^ 16 * q ^ 14 + 6189767847368568104602873471297206580933658004976729096714148703633408000 * p ^ 17 * q ^ 13 + 5571822448947863598135748722637337391563004524206732324589738859083857920 * p ^ 18 * q ^ 12 + 4413362327591268778258323131184043292779046709159404882959089952358400000 * p ^ 19 * q ^ 11 + 3066596729092859782509257623009739127173379471702971963594467457787494400 * p ^ 20 * q ^ 10 + 1860588827240694092227653378870315915338667353983115184305897684651212800 * p ^ 21 * q ^ 9 + 979279368006022209584663374750668678898880342885793197663008011445075968 * p ^ 22 * q ^ 8 + 443089565616187829567641190837292095291007696962198513442722410817126400 * p ^ 23 * q ^ 7 + 170207953265049607139743519666102792494979283260619448111067025482186752 * p ^ 24 * q ^ 6 + 54549608205176547555063550069268235431809232408961907492483353280512000 * p ^ 25 * q ^ 5 + 14224310016128336709129829576679103413681224020855216314671974674595840 * p ^ 26 * q ^ 4 + 2905824558929115622459146730991062764300959809092117612707195125760000 * p ^ 27 * q ^ 3 + 437246877680044698580802560148774939448729921286501660809803228774400 * p ^ 28 * q ^ 2 + 43194494629681658606567941042130450755796816472881471412900986880000 * p ^ 29 * q ^ 1 + 2106245833371143733958360553673408646377901908010982225086219550720 * p ^ 30 := by
    intro p q hp hq; positivity
  have key : (66555937033867822607895549241096482953017615834735226163 : ℝ) * peval [-30, 0, 86, 224, 410, 640, 814, 1024, 1266, 1536, 1670, 1824, 1994, 2176, 2366, 2560, 2366, 2176, 1994, 1824, 1670, 1536, 1266, 1024, 814, 640, 410, 224, 86, 0, -30] x =
      (fun p q : ℝ => 490904797029579978512147356365522781290084128559636240281795078 * q ^ 30 + 28633263796909548690183871813472731253967057260581687033880102400 * p ^ 1 * q ^ 29 + 654193099988220707147549955484109599082831403741012643772187869184 * p ^ 2 * q ^ 28 + 8752135977569476613321476822618962674229182246882883001971743129600 * p ^ 3 * q ^ 27 + 80293433541610536874925265380914300042499435221712636120180056915968 * p ^ 4 * q ^ 26 + 548783114353123750614498063841255630213158502893991958880574242816000 * p ^ 5 * q ^ 25 + 2938619161843024090715324918474020141441752733186898469281945113067520 * p ^ 6 * q ^ 24 + 12747155574543552897524335135622209422027371068883084637425281007616000 * p ^ 7 * q ^ 23 + 45856443144759702888928538938630199664370165226266593451906175745392640 * p ^ 8 * q ^ 22 + 139168793643074675082300045568830841968946153963250327252763420917760000 * p ^ 9 * q ^ 21 + 360912975006446280642783455006275933241039517477404828644499614146232320 * p ^ 10 * q ^ 20 + 807613022759769031683085203739135582767026862931964481857959308243763200 * p ^ 11 * q ^ 19 + 1570923772328951049559541982997545249366034956162668281773247158395338752 * p ^ 12 * q ^ 18 + 2671007554024262514501451604068204857934301750956149862677320713214361600 * p ^ 13 * q ^ 17 + 3985927816059353150976356863618934135626536595482897642826701708101419008 * p ^ 14 * q ^ 16 + 5235101242466759358392491821315787871385647513483036701332610882732032000 * p ^ 15 * q ^ 15 + 6061291791157119221641792893079025292959729510737668096361611332828528640 * p ^ 16 * q ^ 14 + 6189767847368568104602873471297206580933658004976729096714148703633408000 * p ^ 17 * q ^ 13 + 5571822448947863598135748722637337391563004524206732324589738859083857920 * p ^ 18 * q ^ 12 + 4413362327591268778258323131184043292779046709159404882959089952358400000 * p ^ 19 * q ^ 11 + 3066596729092859782509257623009739127173379471702971963594467457787494400 * p ^ 20 * q ^ 10 + 1860588827240694092227653378870315915338667353983115184305897684651212800 * p ^ 21 * q ^ 9 + 979279368006022209584663374750668678898880342885793197663008011445075968 * p ^ 22 * q ^ 8 + 443089565616187829567641190837292095291007696962198513442722410817126400 * p ^ 23 * q ^ 7 + 170207953265049607139743519666102792494979283260619448111067025482186752 * p ^ 24 * q ^ 6 + 54549608205176547555063550069268235431809232408961907492483353280512000 * p ^ 25 * q ^ 5 + 14224310016128336709129829576679103413681224020855216314671974674595840 * p ^ 26 * q ^ 4 + 2905824558929115622459146730991062764300959809092117612707195125760000 * p ^ 27 * q ^ 3 + 437246877680044698580802560148774939448729921286501660809803228774400 * p ^ 28 * q ^ 2 + 43194494629681658606567941042130450755796816472881471412900986880000 * p ^ 29 * q ^ 1 + 2106245833371143733958360553673408646377901908010982225086219550720 * p ^ 30) (x - ((47 : ℝ) / 128)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (66555937033867822607895549241096482953017615834735226163 : ℝ) * peval [-30, 0, 86, 224, 410, 640, 814, 1024, 1266, 1536, 1670, 1824, 1994, 2176, 2366, 2560, 2366, 2176, 1994, 1824, 1670, 1536, 1266, 1024, 814, 640, 410, 224, 86, 0, -30] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posA6 (x : ℝ) (hu : ((47 : ℝ) / 128) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-30, 0, 86, 224, 410, 640, 814, 1024, 1266, 1536, 1670, 1824, 1994, 2176, 2366, 2560, 2366, 2176, 1994, 1824, 1670, 1536, 1266, 1024, 814, 640, 410, 224, 86, 0, -30] x := N16_posA6_0 x hu hv

theorem N16_A6_pos (x : ℝ) (h0 : 0 < x) (hlo : ((47 : ℝ) / 128) < x) (h1 : x < 1) :
    0 < aBlock 16 6 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posA6 x hlo h1)
  rw [← N16_linkA6 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N16_linkB7 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 16 7 x * (2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [63, 252, 630, 1260, 2205, 3528, 5292, 7448, 9947, 12596, 15202, 17572, 19513, 20832, 21336, 20832, 19513, 17572, 15202, 12596, 9947, 7448, 5292, 3528, 2205, 1260, 630, 252, 63] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (16 - 7 : ℕ) = 9 by norm_num, show (2 * 16 - 7 : ℕ) = 25 by norm_num, show (16 + 7 : ℕ) = 23 by norm_num, show (16 - 7 + 16 : ℕ) = 25 by norm_num, show (7 + 16 : ℕ) = 23 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posB7_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [63, 252, 630, 1260, 2205, 3528, 5292, 7448, 9947, 12596, 15202, 17572, 19513, 20832, 21336, 20832, 19513, 17572, 15202, 12596, 9947, 7448, 5292, 3528, 2205, 1260, 630, 252, 63] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 63 * q ^ 28 + 2016 * p ^ 1 * q ^ 27 + 31248 * p ^ 2 * q ^ 26 + 312480 * p ^ 3 * q ^ 25 + 2265480 * p ^ 4 * q ^ 24 + 12686688 * p ^ 5 * q ^ 23 + 57090096 * p ^ 6 * q ^ 22 + 212048816 * p ^ 7 * q ^ 21 + 662650100 * p ^ 8 * q ^ 20 + 1767040656 * p ^ 9 * q ^ 19 + 4064010208 * p ^ 10 * q ^ 18 + 8127089008 * p ^ 11 * q ^ 17 + 14218743336 * p ^ 12 * q ^ 16 + 21863426480 * p ^ 13 * q ^ 15 + 29641758944 * p ^ 14 * q ^ 14 + 35504858592 * p ^ 15 * q ^ 13 + 37603949148 * p ^ 16 * q ^ 12 + 35203831200 * p ^ 17 * q ^ 11 + 29083584960 * p ^ 18 * q ^ 10 + 21137573856 * p ^ 19 * q ^ 9 + 13448518608 * p ^ 20 * q ^ 8 + 7436300256 * p ^ 21 * q ^ 7 + 3536515584 * p ^ 22 * q ^ 6 + 1425001536 * p ^ 23 * q ^ 5 + 475877808 * p ^ 24 * q ^ 4 + 127318464 * p ^ 25 * q ^ 3 + 25796736 * p ^ 26 * q ^ 2 + 3556224 * p ^ 27 * q ^ 1 + 254016 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [63, 252, 630, 1260, 2205, 3528, 5292, 7448, 9947, 12596, 15202, 17572, 19513, 20832, 21336, 20832, 19513, 17572, 15202, 12596, 9947, 7448, 5292, 3528, 2205, 1260, 630, 252, 63] x =
      (fun p q : ℝ => 63 * q ^ 28 + 2016 * p ^ 1 * q ^ 27 + 31248 * p ^ 2 * q ^ 26 + 312480 * p ^ 3 * q ^ 25 + 2265480 * p ^ 4 * q ^ 24 + 12686688 * p ^ 5 * q ^ 23 + 57090096 * p ^ 6 * q ^ 22 + 212048816 * p ^ 7 * q ^ 21 + 662650100 * p ^ 8 * q ^ 20 + 1767040656 * p ^ 9 * q ^ 19 + 4064010208 * p ^ 10 * q ^ 18 + 8127089008 * p ^ 11 * q ^ 17 + 14218743336 * p ^ 12 * q ^ 16 + 21863426480 * p ^ 13 * q ^ 15 + 29641758944 * p ^ 14 * q ^ 14 + 35504858592 * p ^ 15 * q ^ 13 + 37603949148 * p ^ 16 * q ^ 12 + 35203831200 * p ^ 17 * q ^ 11 + 29083584960 * p ^ 18 * q ^ 10 + 21137573856 * p ^ 19 * q ^ 9 + 13448518608 * p ^ 20 * q ^ 8 + 7436300256 * p ^ 21 * q ^ 7 + 3536515584 * p ^ 22 * q ^ 6 + 1425001536 * p ^ 23 * q ^ 5 + 475877808 * p ^ 24 * q ^ 4 + 127318464 * p ^ 25 * q ^ 3 + 25796736 * p ^ 26 * q ^ 2 + 3556224 * p ^ 27 * q ^ 1 + 254016 * p ^ 28) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [63, 252, 630, 1260, 2205, 3528, 5292, 7448, 9947, 12596, 15202, 17572, 19513, 20832, 21336, 20832, 19513, 17572, 15202, 12596, 9947, 7448, 5292, 3528, 2205, 1260, 630, 252, 63] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posB7 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [63, 252, 630, 1260, 2205, 3528, 5292, 7448, 9947, 12596, 15202, 17572, 19513, 20832, 21336, 20832, 19513, 17572, 15202, 12596, 9947, 7448, 5292, 3528, 2205, 1260, 630, 252, 63] x := N16_posB7_0 x hu hv

theorem N16_B7_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 16 7 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posB7 x h0 h1)
  rw [← N16_linkB7 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N16_linkA7 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 16 7 x * (2 * (1 - x) ^ 2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [-33, -6, 77, 212, 395, 622, 889, 1080, 1303, 1410, 1541, 1692, 1859, 2038, 2225, 2416, 2225, 2038, 1859, 1692, 1541, 1410, 1303, 1080, 889, 622, 395, 212, 77, -6, -33] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (16 - 7 : ℕ) = 9 by norm_num, show (2 * 16 - 7 : ℕ) = 25 by norm_num, show (16 + 7 : ℕ) = 23 by norm_num, show (16 - 7 + 16 : ℕ) = 25 by norm_num, show (7 + 16 : ℕ) = 23 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posA7_0 (x : ℝ) (hu : ((47 : ℝ) / 128) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-33, -6, 77, 212, 395, 622, 889, 1080, 1303, 1410, 1541, 1692, 1859, 2038, 2225, 2416, 2225, 2038, 1859, 1692, 1541, 1410, 1303, 1080, 889, 622, 395, 212, 77, -6, -33] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 479071337305915349826841979738464411347008383321982281777223125 * q ^ 30 + 133113124089747095953687476180120289598594450261782877707553932800 * p ^ 1 * q ^ 29 + 3983837114900636689357415777440710490818251052551130260039178190848 * p ^ 2 * q ^ 28 + 60050979317599267962481724511748286113408333278504457337948864512000 * p ^ 3 * q ^ 27 + 589634004573501639176014693964511515784039771446618700775592778792960 * p ^ 4 * q ^ 26 + 4211473761355947616332945682037851735219310346525943150706457732710400 * p ^ 5 * q ^ 25 + 23260575554168776438233812473219730958703452380641597955776328906047488 * p ^ 6 * q ^ 24 + 103242341214824151910164242824141126108516748635996723177582822529433600 * p ^ 7 * q ^ 23 + 378019031557883389340960131737015373976433669506712436471685759467782144 * p ^ 8 * q ^ 22 + 1163340170359504031204048821524649350641613650953474314895499099897856000 * p ^ 9 * q ^ 21 + 3050845501437372743227866738930291799010237092602131610891326011622817792 * p ^ 10 * q ^ 20 + 6888801219548718046387749440352485177077598191436726056652188117735833600 * p ^ 11 * q ^ 19 + 13497890248597791476881340024238278599240118575913049536889803871507972096 * p ^ 12 * q ^ 18 + 23084767448272222914771946721193371666182017756108657426322468428683673600 * p ^ 13 * q ^ 17 + 34607758182121980980512121954958494898182560679007306760860087955193069568 * p ^ 14 * q ^ 16 + 45611892912563223825367531962912501933262023857149231882440824114380800000 * p ^ 15 * q ^ 15 + 52940595441791283415575159603216940012287952837707168802718669221134336000 * p ^ 16 * q ^ 14 + 54146172091400514943917596810248119034721877443336107472152763393743257600 * p ^ 17 * q ^ 13 + 48774504706552661330136366128698547910609373169207335697189906666250829824 * p ^ 18 * q ^ 12 + 38630693863313106908643305767101495700837789376994238125246783500058624000 * p ^ 19 * q ^ 11 + 26821751301236855718444937181236669121197593848142398615640850183470710784 * p ^ 20 * q ^ 10 + 16251318626462377942400957023230525576462111751374115391067088646622412800 * p ^ 21 * q ^ 9 + 8537676235748029665471637666960773379229881217628521797259238616458067968 * p ^ 22 * q ^ 8 + 3854447330213746094063692553830568926355078244404910090384655677967564800 * p ^ 23 * q ^ 7 + 1477050966647729856584124296816005949121650236297117214717672392872689664 * p ^ 24 * q ^ 6 + 472208618120282281241011573528866128575578717318411019966742162256691200 * p ^ 25 * q ^ 5 + 122846610858813507381162872807400019567044654771063677169463702673424384 * p ^ 26 * q ^ 4 + 25046027708242250060307338456628280622106225968706092840276355501260800 * p ^ 27 * q ^ 3 + 3763245544858930338617885936044141225661975431694843595408100873469952 * p ^ 28 * q ^ 2 + 371472653815262264016484292962321876499852621666780654150948487168000 * p ^ 29 * q ^ 1 + 18113714166991836112041900761591314358849956408894447135741488136192 * p ^ 30 := by
    intro p q hp hq; positivity
  have key : (599003433304810403471059943169868346577158542512617035467 : ℝ) * peval [-33, -6, 77, 212, 395, 622, 889, 1080, 1303, 1410, 1541, 1692, 1859, 2038, 2225, 2416, 2225, 2038, 1859, 1692, 1541, 1410, 1303, 1080, 889, 622, 395, 212, 77, -6, -33] x =
      (fun p q : ℝ => 479071337305915349826841979738464411347008383321982281777223125 * q ^ 30 + 133113124089747095953687476180120289598594450261782877707553932800 * p ^ 1 * q ^ 29 + 3983837114900636689357415777440710490818251052551130260039178190848 * p ^ 2 * q ^ 28 + 60050979317599267962481724511748286113408333278504457337948864512000 * p ^ 3 * q ^ 27 + 589634004573501639176014693964511515784039771446618700775592778792960 * p ^ 4 * q ^ 26 + 4211473761355947616332945682037851735219310346525943150706457732710400 * p ^ 5 * q ^ 25 + 23260575554168776438233812473219730958703452380641597955776328906047488 * p ^ 6 * q ^ 24 + 103242341214824151910164242824141126108516748635996723177582822529433600 * p ^ 7 * q ^ 23 + 378019031557883389340960131737015373976433669506712436471685759467782144 * p ^ 8 * q ^ 22 + 1163340170359504031204048821524649350641613650953474314895499099897856000 * p ^ 9 * q ^ 21 + 3050845501437372743227866738930291799010237092602131610891326011622817792 * p ^ 10 * q ^ 20 + 6888801219548718046387749440352485177077598191436726056652188117735833600 * p ^ 11 * q ^ 19 + 13497890248597791476881340024238278599240118575913049536889803871507972096 * p ^ 12 * q ^ 18 + 23084767448272222914771946721193371666182017756108657426322468428683673600 * p ^ 13 * q ^ 17 + 34607758182121980980512121954958494898182560679007306760860087955193069568 * p ^ 14 * q ^ 16 + 45611892912563223825367531962912501933262023857149231882440824114380800000 * p ^ 15 * q ^ 15 + 52940595441791283415575159603216940012287952837707168802718669221134336000 * p ^ 16 * q ^ 14 + 54146172091400514943917596810248119034721877443336107472152763393743257600 * p ^ 17 * q ^ 13 + 48774504706552661330136366128698547910609373169207335697189906666250829824 * p ^ 18 * q ^ 12 + 38630693863313106908643305767101495700837789376994238125246783500058624000 * p ^ 19 * q ^ 11 + 26821751301236855718444937181236669121197593848142398615640850183470710784 * p ^ 20 * q ^ 10 + 16251318626462377942400957023230525576462111751374115391067088646622412800 * p ^ 21 * q ^ 9 + 8537676235748029665471637666960773379229881217628521797259238616458067968 * p ^ 22 * q ^ 8 + 3854447330213746094063692553830568926355078244404910090384655677967564800 * p ^ 23 * q ^ 7 + 1477050966647729856584124296816005949121650236297117214717672392872689664 * p ^ 24 * q ^ 6 + 472208618120282281241011573528866128575578717318411019966742162256691200 * p ^ 25 * q ^ 5 + 122846610858813507381162872807400019567044654771063677169463702673424384 * p ^ 26 * q ^ 4 + 25046027708242250060307338456628280622106225968706092840276355501260800 * p ^ 27 * q ^ 3 + 3763245544858930338617885936044141225661975431694843595408100873469952 * p ^ 28 * q ^ 2 + 371472653815262264016484292962321876499852621666780654150948487168000 * p ^ 29 * q ^ 1 + 18113714166991836112041900761591314358849956408894447135741488136192 * p ^ 30) (x - ((47 : ℝ) / 128)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (599003433304810403471059943169868346577158542512617035467 : ℝ) * peval [-33, -6, 77, 212, 395, 622, 889, 1080, 1303, 1410, 1541, 1692, 1859, 2038, 2225, 2416, 2225, 2038, 1859, 1692, 1541, 1410, 1303, 1080, 889, 622, 395, 212, 77, -6, -33] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posA7 (x : ℝ) (hu : ((47 : ℝ) / 128) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-33, -6, 77, 212, 395, 622, 889, 1080, 1303, 1410, 1541, 1692, 1859, 2038, 2225, 2416, 2225, 2038, 1859, 1692, 1541, 1410, 1303, 1080, 889, 622, 395, 212, 77, -6, -33] x := N16_posA7_0 x hu hv

theorem N16_A7_pos (x : ℝ) (h0 : 0 < x) (hlo : ((47 : ℝ) / 128) < x) (h1 : x < 1) :
    0 < aBlock 16 7 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posA7 x hlo h1)
  rw [← N16_linkA7 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

theorem N16_linkB8 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    bBlock 16 8 x * (2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [64, 256, 640, 1280, 2240, 3584, 5376, 7680, 10304, 13056, 15744, 18176, 20160, 21504, 22016, 21504, 20160, 18176, 15744, 13056, 10304, 7680, 5376, 3584, 2240, 1280, 640, 256, 64] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold bBlock W
  simp only [show (16 - 8 : ℕ) = 8 by norm_num, show (2 * 16 - 8 : ℕ) = 24 by norm_num, show (16 + 8 : ℕ) = 24 by norm_num, show (16 - 8 + 16 : ℕ) = 24 by norm_num, show (8 + 16 : ℕ) = 24 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

theorem N16_posB8_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [64, 256, 640, 1280, 2240, 3584, 5376, 7680, 10304, 13056, 15744, 18176, 20160, 21504, 22016, 21504, 20160, 18176, 15744, 13056, 10304, 7680, 5376, 3584, 2240, 1280, 640, 256, 64] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 64 * q ^ 28 + 2048 * p ^ 1 * q ^ 27 + 31744 * p ^ 2 * q ^ 26 + 317440 * p ^ 3 * q ^ 25 + 2301440 * p ^ 4 * q ^ 24 + 12888064 * p ^ 5 * q ^ 23 + 57996288 * p ^ 6 * q ^ 22 + 215414784 * p ^ 7 * q ^ 21 + 673170944 * p ^ 8 * q ^ 20 + 1795117056 * p ^ 9 * q ^ 19 + 4128712704 * p ^ 10 * q ^ 18 + 8257048576 * p ^ 11 * q ^ 17 + 14448021504 * p ^ 12 * q ^ 16 + 22221029376 * p ^ 13 * q ^ 15 + 30137421824 * p ^ 14 * q ^ 14 + 36117651456 * p ^ 15 * q ^ 13 + 38280864768 * p ^ 16 * q ^ 12 + 35871539200 * p ^ 17 * q ^ 11 + 29669842944 * p ^ 18 * q ^ 10 + 21593112576 * p ^ 19 * q ^ 9 + 13759176704 * p ^ 20 * q ^ 8 + 7620231168 * p ^ 21 * q ^ 7 + 3629776896 * p ^ 22 * q ^ 6 + 1464795136 * p ^ 23 * q ^ 5 + 489832448 * p ^ 24 * q ^ 4 + 131203072 * p ^ 25 * q ^ 3 + 26607616 * p ^ 26 * q ^ 2 + 3670016 * p ^ 27 * q ^ 1 + 262144 * p ^ 28 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [64, 256, 640, 1280, 2240, 3584, 5376, 7680, 10304, 13056, 15744, 18176, 20160, 21504, 22016, 21504, 20160, 18176, 15744, 13056, 10304, 7680, 5376, 3584, 2240, 1280, 640, 256, 64] x =
      (fun p q : ℝ => 64 * q ^ 28 + 2048 * p ^ 1 * q ^ 27 + 31744 * p ^ 2 * q ^ 26 + 317440 * p ^ 3 * q ^ 25 + 2301440 * p ^ 4 * q ^ 24 + 12888064 * p ^ 5 * q ^ 23 + 57996288 * p ^ 6 * q ^ 22 + 215414784 * p ^ 7 * q ^ 21 + 673170944 * p ^ 8 * q ^ 20 + 1795117056 * p ^ 9 * q ^ 19 + 4128712704 * p ^ 10 * q ^ 18 + 8257048576 * p ^ 11 * q ^ 17 + 14448021504 * p ^ 12 * q ^ 16 + 22221029376 * p ^ 13 * q ^ 15 + 30137421824 * p ^ 14 * q ^ 14 + 36117651456 * p ^ 15 * q ^ 13 + 38280864768 * p ^ 16 * q ^ 12 + 35871539200 * p ^ 17 * q ^ 11 + 29669842944 * p ^ 18 * q ^ 10 + 21593112576 * p ^ 19 * q ^ 9 + 13759176704 * p ^ 20 * q ^ 8 + 7620231168 * p ^ 21 * q ^ 7 + 3629776896 * p ^ 22 * q ^ 6 + 1464795136 * p ^ 23 * q ^ 5 + 489832448 * p ^ 24 * q ^ 4 + 131203072 * p ^ 25 * q ^ 3 + 26607616 * p ^ 26 * q ^ 2 + 3670016 * p ^ 27 * q ^ 1 + 262144 * p ^ 28) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [64, 256, 640, 1280, 2240, 3584, 5376, 7680, 10304, 13056, 15744, 18176, 20160, 21504, 22016, 21504, 20160, 18176, 15744, 13056, 10304, 7680, 5376, 3584, 2240, 1280, 640, 256, 64] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posB8 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [64, 256, 640, 1280, 2240, 3584, 5376, 7680, 10304, 13056, 15744, 18176, 20160, 21504, 22016, 21504, 20160, 18176, 15744, 13056, 10304, 7680, 5376, 3584, 2240, 1280, 640, 256, 64] x := N16_posB8_0 x hu hv

theorem N16_B8_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    0 < bBlock 16 8 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := mul_pos (pow_pos h1x 4) (N16_posB8 x h0 h1)
  rw [← N16_linkB8 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos two_pos (pow_pos hX 2))).mp hn

theorem N16_linkA8 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 16 8 x * (2 * (1 - x) ^ 2 * (1 - x ^ 16) ^ 2) = (1 - x) ^ 4 * peval [-34, -8, 74, 208, 390, 616, 882, 1184, 1262, 1368, 1498, 1648, 1814, 1992, 2178, 2368, 2178, 1992, 1814, 1648, 1498, 1368, 1262, 1184, 882, 616, 390, 208, 74, -8, -34] x := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ 16 ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold aBlock cR
  simp only [show (16 - 8 : ℕ) = 8 by norm_num, show (2 * 16 - 8 : ℕ) = 24 by norm_num, show (16 + 8 : ℕ) = 24 by norm_num, show (16 - 8 + 16 : ℕ) = 24 by norm_num, show (8 + 16 : ℕ) = 24 by norm_num, show (16 + 16 : ℕ) = 32 by norm_num]
  push_cast
  field_simp
  simp only [peval]
  ring

/-- The minimal polynomial `f_16` of `ξ*(16)` (constant term first). -/
def f16 : List ℤ := [17, 4, -37, -104, -195, -308, -441, -592, -648, -688, -712, -720, -712, -688, -648, -592, -441, -308, -195, -104, -37, 4, 17]

theorem N16_crit_factor (x : ℝ) : peval [-34, -8, 74, 208, 390, 616, 882, 1184, 1262, 1368, 1498, 1648, 1814, 1992, 2178, 2368, 2178, 1992, 1814, 1648, 1498, 1368, 1262, 1184, 882, 616, 390, 208, 74, -8, -34] x = -(peval [2, 0, 0, 0, 0, 0, 0, 0, 2] x * peval f16 x) := by
  simp only [f16, peval]; ring

theorem N16_posH_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 0, 0, 0, 0, 0, 0, 0, 2] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 2 * q ^ 8 + 16 * p ^ 1 * q ^ 7 + 56 * p ^ 2 * q ^ 6 + 112 * p ^ 3 * q ^ 5 + 140 * p ^ 4 * q ^ 4 + 112 * p ^ 5 * q ^ 3 + 56 * p ^ 6 * q ^ 2 + 16 * p ^ 7 * q ^ 1 + 4 * p ^ 8 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [2, 0, 0, 0, 0, 0, 0, 0, 2] x =
      (fun p q : ℝ => 2 * q ^ 8 + 16 * p ^ 1 * q ^ 7 + 56 * p ^ 2 * q ^ 6 + 112 * p ^ 3 * q ^ 5 + 140 * p ^ 4 * q ^ 4 + 112 * p ^ 5 * q ^ 3 + 56 * p ^ 6 * q ^ 2 + 16 * p ^ 7 * q ^ 1 + 4 * p ^ 8) (x - (0 : ℝ)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [2, 0, 0, 0, 0, 0, 0, 0, 2] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_posH (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < (1 : ℝ)) :
    0 < peval [2, 0, 0, 0, 0, 0, 0, 0, 2] x := N16_posH_0 x hu hv

theorem N16_dneg_0 (x : ℝ) (hu : ((1 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-4, 74, 312, 780, 1540, 2646, 4144, 5184, 6192, 7120, 7920, 8544, 8944, 9072, 8880, 7056, 5236, 3510, 1976, 740, -84, -374] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 39855276687118872574451018 * q ^ 21 + 3083250416495700768654264288 * p ^ 1 * q ^ 20 + 61827406619493274068779796480 * p ^ 2 * q ^ 19 + 662498280927754997073965219840 * p ^ 3 * q ^ 18 + 4668622899575814471222924410880 * p ^ 4 * q ^ 17 + 23853974556710831884699258847232 * p ^ 5 * q ^ 16 + 93421862606894698957166862139392 * p ^ 6 * q ^ 15 + 290444244599568395897051692400640 * p ^ 7 * q ^ 14 + 733653264517543216889518903787520 * p ^ 8 * q ^ 13 + 1529434308486461512688904636989440 * p ^ 9 * q ^ 12 + 2658686256113383661353899537727488 * p ^ 10 * q ^ 11 + 3878114025712421547490137077710848 * p ^ 11 * q ^ 10 + 4760565918054095210522901049507840 * p ^ 12 * q ^ 9 + 4917477114199022970068539036139520 * p ^ 13 * q ^ 8 + 4261196346982110896715511093002240 * p ^ 14 * q ^ 7 + 3077608942760213555584840973156352 * p ^ 15 * q ^ 6 + 1832513655508657747709567195676672 * p ^ 16 * q ^ 5 + 883815242661683010209262002503680 * p ^ 17 * q ^ 4 + 335263081655243434890828279971840 * p ^ 18 * q ^ 3 + 94900289983486113533098910023680 * p ^ 19 * q ^ 2 + 18014116595418577078996149731328 * p ^ 20 * q ^ 1 + 1729402234881676244034076540928 * p ^ 21 := by
    intro p q hp hq; positivity
  have key : (4987885095119476318359375 : ℝ) * peval [-4, 74, 312, 780, 1540, 2646, 4144, 5184, 6192, 7120, 7920, 8544, 8944, 9072, 8880, 7056, 5236, 3510, 1976, 740, -84, -374] x =
      (fun p q : ℝ => 39855276687118872574451018 * q ^ 21 + 3083250416495700768654264288 * p ^ 1 * q ^ 20 + 61827406619493274068779796480 * p ^ 2 * q ^ 19 + 662498280927754997073965219840 * p ^ 3 * q ^ 18 + 4668622899575814471222924410880 * p ^ 4 * q ^ 17 + 23853974556710831884699258847232 * p ^ 5 * q ^ 16 + 93421862606894698957166862139392 * p ^ 6 * q ^ 15 + 290444244599568395897051692400640 * p ^ 7 * q ^ 14 + 733653264517543216889518903787520 * p ^ 8 * q ^ 13 + 1529434308486461512688904636989440 * p ^ 9 * q ^ 12 + 2658686256113383661353899537727488 * p ^ 10 * q ^ 11 + 3878114025712421547490137077710848 * p ^ 11 * q ^ 10 + 4760565918054095210522901049507840 * p ^ 12 * q ^ 9 + 4917477114199022970068539036139520 * p ^ 13 * q ^ 8 + 4261196346982110896715511093002240 * p ^ 14 * q ^ 7 + 3077608942760213555584840973156352 * p ^ 15 * q ^ 6 + 1832513655508657747709567195676672 * p ^ 16 * q ^ 5 + 883815242661683010209262002503680 * p ^ 17 * q ^ 4 + 335263081655243434890828279971840 * p ^ 18 * q ^ 3 + 94900289983486113533098910023680 * p ^ 19 * q ^ 2 + 18014116595418577078996149731328 * p ^ 20 * q ^ 1 + 1729402234881676244034076540928 * p ^ 21) (x - ((1 : ℝ) / 16)) ((1 : ℝ) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (4987885095119476318359375 : ℝ) * peval [-4, 74, 312, 780, 1540, 2646, 4144, 5184, 6192, 7120, 7920, 8544, 8944, 9072, 8880, 7056, 5236, 3510, 1976, 740, -84, -374] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_dneg (x : ℝ) (hu : ((1 : ℝ) / 16) < x) (hv : x < (1 : ℝ)) :
    0 < peval [-4, 74, 312, 780, 1540, 2646, 4144, 5184, 6192, 7120, 7920, 8544, 8944, 9072, 8880, 7056, 5236, 3510, 1976, 740, -84, -374] x := N16_dneg_0 x hu hv

theorem N16_fleft_0 (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < ((1 : ℝ) / 16)) :
    0 < peval [17, 4, -37, -104, -195, -308, -441, -592, -648, -688, -712, -720, -712, -688, -648, -592, -441, -308, -195, -104, -37, 4, 17] x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < 5261245166962866168321277952 * q ^ 22 + 115824764925638391970249310208 * p ^ 1 * q ^ 21 + 1216927699614658405213556178944 * p ^ 2 * q ^ 20 + 8117663057014092194643900891136 * p ^ 3 * q ^ 19 + 38580263190386885847239060619264 * p ^ 4 * q ^ 18 + 138960344294213548849039824388096 * p ^ 5 * q ^ 17 + 393907238307152449805502805704704 * p ^ 6 * q ^ 16 + 900747549925128789162765588627456 * p ^ 7 * q ^ 15 + 1689557097909154449154102463037440 * p ^ 8 * q ^ 14 + 2629104043056082551598458540654592 * p ^ 9 * q ^ 13 + 3418856581888031124602506782441472 * p ^ 10 * q ^ 12 + 3730603493889019168823801884442624 * p ^ 11 * q ^ 11 + 3420420999663427661451048343568384 * p ^ 12 * q ^ 10 + 2631504422754699263721799600308224 * p ^ 13 * q ^ 9 + 1691861434031026383781488819699712 * p ^ 14 * q ^ 8 + 902375914746042474560800486850560 * p ^ 15 * q ^ 7 + 394790522223366833224069060493312 * p ^ 16 * q ^ 6 + 139330696787536880456996417961984 * p ^ 17 * q ^ 5 + 38698849600387708753329728061440 * p ^ 18 * q ^ 4 + 8145796289722544048581246812160 * p ^ 19 * q ^ 3 + 1221598819692150411816586599168 * p ^ 20 * q ^ 2 + 116310020534013130520892814912 * p ^ 21 * q ^ 1 + 5285007511696636893105380177 * p ^ 22 := by
    intro p q hp hq; positivity
  have key : (1 : ℝ) * peval [17, 4, -37, -104, -195, -308, -441, -592, -648, -688, -712, -720, -712, -688, -648, -592, -441, -308, -195, -104, -37, 4, 17] x =
      (fun p q : ℝ => 5261245166962866168321277952 * q ^ 22 + 115824764925638391970249310208 * p ^ 1 * q ^ 21 + 1216927699614658405213556178944 * p ^ 2 * q ^ 20 + 8117663057014092194643900891136 * p ^ 3 * q ^ 19 + 38580263190386885847239060619264 * p ^ 4 * q ^ 18 + 138960344294213548849039824388096 * p ^ 5 * q ^ 17 + 393907238307152449805502805704704 * p ^ 6 * q ^ 16 + 900747549925128789162765588627456 * p ^ 7 * q ^ 15 + 1689557097909154449154102463037440 * p ^ 8 * q ^ 14 + 2629104043056082551598458540654592 * p ^ 9 * q ^ 13 + 3418856581888031124602506782441472 * p ^ 10 * q ^ 12 + 3730603493889019168823801884442624 * p ^ 11 * q ^ 11 + 3420420999663427661451048343568384 * p ^ 12 * q ^ 10 + 2631504422754699263721799600308224 * p ^ 13 * q ^ 9 + 1691861434031026383781488819699712 * p ^ 14 * q ^ 8 + 902375914746042474560800486850560 * p ^ 15 * q ^ 7 + 394790522223366833224069060493312 * p ^ 16 * q ^ 6 + 139330696787536880456996417961984 * p ^ 17 * q ^ 5 + 38698849600387708753329728061440 * p ^ 18 * q ^ 4 + 8145796289722544048581246812160 * p ^ 19 * q ^ 3 + 1221598819692150411816586599168 * p ^ 20 * q ^ 2 + 116310020534013130520892814912 * p ^ 21 * q ^ 1 + 5285007511696636893105380177 * p ^ 22) (x - (0 : ℝ)) (((1 : ℝ) / 16) - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < (1 : ℝ) * peval [17, 4, -37, -104, -195, -308, -441, -592, -648, -688, -712, -720, -712, -688, -648, -592, -441, -308, -195, -104, -37, 4, 17] x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)

theorem N16_fleft (x : ℝ) (hu : (0 : ℝ) < x) (hv : x < ((1 : ℝ) / 16)) :
    0 < peval [17, 4, -37, -104, -195, -308, -441, -592, -648, -688, -712, -720, -712, -688, -648, -592, -441, -308, -195, -104, -37, 4, 17] x := N16_fleft_0 x hu hv

theorem N16_dlink (x : ℝ) : dpeval f16 x = -peval [-4, 74, 312, 780, 1540, 2646, 4144, 5184, 6192, 7120, 7920, 8544, 8944, 9072, 8880, 7056, 5236, 3510, 1976, 740, -84, -374] x := by
  simp only [f16, dpeval, peval]; ring

/-- `f_16` changes sign exactly once in `(0,1)`, at a point of `(37168181390604545663/100000000000000000000, 290376417114098013/781250000000000000)`. -/
theorem N16_root : ∃ ξs, ((37168181390604545663 : ℝ) / 100000000000000000000) < ξs ∧ ξs < ((290376417114098013 : ℝ) / 781250000000000000) ∧ peval f16 ξs = 0 ∧
    ∀ x, 0 < x → x < 1 → (x < ξs → 0 < peval f16 x) ∧ (ξs < x → peval f16 x < 0) := by
  apply sign_change (peval f16) ((1 : ℝ) / 16) _ _ (by norm_num) (by norm_num) (by norm_num)
    (continuous_peval _).continuousOn
  · apply peval_strictAntiOn
    intro x hu hv
    rw [N16_dlink]
    linarith [N16_dneg x hu hv]
  · exact (fun x hx0 hx => by
      rcases hx.lt_or_eq with h | h
      · exact N16_fleft x hx0 h
      · rw [h]; norm_num [f16, peval])
  · norm_num [f16, peval]
  · norm_num [f16, peval]

theorem N16_A8_neg (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : 0 < peval f16 x) :
    aBlock 16 8 x < 0 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : (1 - x) ^ 4 * peval [-34, -8, 74, 208, 390, 616, 882, 1184, 1262, 1368, 1498, 1648, 1814, 1992, 2178, 2368, 2178, 1992, 1814, 1648, 1498, 1368, 1262, 1184, 882, 616, 390, 208, 74, -8, -34] x < 0 := by
    rw [N16_crit_factor]
    have := mul_pos (N16_posH x h0 h1) hf
    nlinarith [pow_pos h1x 4]
  rw [← N16_linkA8 x h0 h1] at hn
  have hd := (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))
  by_contra hc
  push Not at hc
  nlinarith [mul_nonneg hc hd.le]

theorem N16_A8_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : peval f16 x < 0) :
    0 < aBlock 16 8 x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ 16 := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : 0 < (1 - x) ^ 4 * peval [-34, -8, 74, 208, 390, 616, 882, 1184, 1262, 1368, 1498, 1648, 1814, 1992, 2178, 2368, 2178, 1992, 1814, 1648, 1498, 1368, 1262, 1184, 882, 616, 390, 208, 74, -8, -34] x := by
    rw [N16_crit_factor]
    have := mul_pos (N16_posH x h0 h1) (neg_pos.mpr hf)
    have h4 := pow_pos h1x 4
    nlinarith
  rw [← N16_linkA8 x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right (mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))).mp hn

/-- Note Thm 6(b) for `N = 16`: there is a unique `ξ* ∈ (0,1)` (the unique root of `f_16` in
`(0,1)`, enclosed in `(37168181390604545663/100000000000000000000, 290376417114098013/781250000000000000)`) such that, for every `ξ ∈ (0,1)`: all `b_m > 0`;
`a_8 < 0` for `ξ < ξ*` (spectrally unstable, since `b_8 > 0`); and all `a_m > 0` for `ξ > ξ*`. -/
theorem cert_N16 : ∃ ξs : ℝ, ((37168181390604545663 : ℝ) / 100000000000000000000) < ξs ∧ ξs < ((290376417114098013 : ℝ) / 781250000000000000) ∧ peval f16 ξs = 0 ∧
    (∀ y, 0 < y → y < 1 → peval f16 y = 0 → y = ξs) ∧
    ∀ x, 0 < x → x < 1 →
      (∀ m, 1 ≤ m → m ≤ 15 → 0 < bBlock 16 m x) ∧
      (x < ξs → aBlock 16 8 x < 0) ∧
      (ξs < x → ∀ m, 1 ≤ m → m ≤ 15 → 0 < aBlock 16 m x) := by
  obtain ⟨ξs, hta, htb, hf0, hs⟩ := N16_root
  refine ⟨ξs, hta, htb, hf0, ?_, ?_⟩
  · intro y hy0 hy1 hy
    rcases lt_trichotomy y ξs with h | h | h
    · linarith [(hs y hy0 hy1).1 h]
    · exact h
    · linarith [(hs y hy0 hy1).2 h]
  intro x h0 h1
  refine ⟨?_, fun hlt => N16_A8_neg x h0 h1 ((hs x h0 h1).1 hlt), ?_⟩
  · intro m hm1 hm2
    interval_cases m
    · exact N16_B1_pos x h0 h1
    · exact N16_B2_pos x h0 h1
    · exact N16_B3_pos x h0 h1
    · exact N16_B4_pos x h0 h1
    · exact N16_B5_pos x h0 h1
    · exact N16_B6_pos x h0 h1
    · exact N16_B7_pos x h0 h1
    · exact N16_B8_pos x h0 h1
    · rw [show bBlock 16 9 x = bBlock 16 7 x from bBlock_symm 16 7 (by norm_num) x]
      exact N16_B7_pos x h0 h1
    · rw [show bBlock 16 10 x = bBlock 16 6 x from bBlock_symm 16 6 (by norm_num) x]
      exact N16_B6_pos x h0 h1
    · rw [show bBlock 16 11 x = bBlock 16 5 x from bBlock_symm 16 5 (by norm_num) x]
      exact N16_B5_pos x h0 h1
    · rw [show bBlock 16 12 x = bBlock 16 4 x from bBlock_symm 16 4 (by norm_num) x]
      exact N16_B4_pos x h0 h1
    · rw [show bBlock 16 13 x = bBlock 16 3 x from bBlock_symm 16 3 (by norm_num) x]
      exact N16_B3_pos x h0 h1
    · rw [show bBlock 16 14 x = bBlock 16 2 x from bBlock_symm 16 2 (by norm_num) x]
      exact N16_B2_pos x h0 h1
    · rw [show bBlock 16 15 x = bBlock 16 1 x from bBlock_symm 16 1 (by norm_num) x]
      exact N16_B1_pos x h0 h1
  · intro hgt m hm1 hm2
    have hlo : ((47 : ℝ) / 128) < x := lt_trans (lt_trans (by norm_num) hta) hgt
    have hfn := (hs x h0 h1).2 hgt
    interval_cases m
    · exact N16_A1_pos x h0 h1
    · exact N16_A2_pos x h0 h1
    · exact N16_A3_pos x h0 hlo h1
    · exact N16_A4_pos x h0 hlo h1
    · exact N16_A5_pos x h0 hlo h1
    · exact N16_A6_pos x h0 hlo h1
    · exact N16_A7_pos x h0 hlo h1
    · exact N16_A8_pos x h0 h1 hfn
    · rw [show aBlock 16 9 x = aBlock 16 7 x from aBlock_symm 16 7 (by norm_num) x]
      exact N16_A7_pos x h0 hlo h1
    · rw [show aBlock 16 10 x = aBlock 16 6 x from aBlock_symm 16 6 (by norm_num) x]
      exact N16_A6_pos x h0 hlo h1
    · rw [show aBlock 16 11 x = aBlock 16 5 x from aBlock_symm 16 5 (by norm_num) x]
      exact N16_A5_pos x h0 hlo h1
    · rw [show aBlock 16 12 x = aBlock 16 4 x from aBlock_symm 16 4 (by norm_num) x]
      exact N16_A4_pos x h0 hlo h1
    · rw [show aBlock 16 13 x = aBlock 16 3 x from aBlock_symm 16 3 (by norm_num) x]
      exact N16_A3_pos x h0 hlo h1
    · rw [show aBlock 16 14 x = aBlock 16 2 x from aBlock_symm 16 2 (by norm_num) x]
      exact N16_A2_pos x h0 h1
    · rw [show aBlock 16 15 x = aBlock 16 1 x from aBlock_symm 16 1 (by norm_num) x]
      exact N16_A1_pos x h0 h1

/-! ### Factorisations, palindromy, Krein window -/

/-- Factorisation of the critical numerator for `N = 8` (note §6, "Structure of the numerator"):
`A_4 = 2(1−ξ)²(1−ξ^8)² a_4 = −2·(ξ−1)⁴ (ξ^4+1) f_8(ξ)`. -/
theorem A_crit_factor_N8 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 8 4 x * (2 * (1 - x) ^ 2 * (1 - x ^ 8) ^ 2)
      = -(2 * (x ^ 4 + 1)) * (x - 1) ^ 4 * peval f8 x := by
  rw [N8_linkA4 x h0 h1, N8_crit_factor]
  simp only [peval]; ring

/-- Factorisation of the critical numerator for `N = 9`:
`A_4 = −1·(ξ−1)⁴ f_9(ξ)`. -/
theorem A_crit_factor_N9 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 9 4 x * (2 * (1 - x) ^ 2 * (1 - x ^ 9) ^ 2) = -1 * (x - 1) ^ 4 * peval f9 x := by
  rw [N9_linkA4 x h0 h1, N9_crit_factor]
  simp only [peval]; ring

/-- Factorisation of the critical numerator for `N = 10` (note §6, "Structure of the numerator"):
`A_5 = 2(1−ξ)²(1−ξ^10)² a_5 = −1·(ξ−1)⁴ (ξ^5+1)(ξ+1) f_10(ξ)`. -/
theorem A_crit_factor_N10 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 10 5 x * (2 * (1 - x) ^ 2 * (1 - x ^ 10) ^ 2)
      = -(1 * (x ^ 5 + 1) * (x + 1)) * (x - 1) ^ 4 * peval f10 x := by
  rw [N10_linkA5 x h0 h1, N10_crit_factor]
  simp only [peval]; ring

/-- Factorisation of the critical numerator for `N = 11`:
`A_5 = −1·(ξ−1)⁴ f_11(ξ)`. -/
theorem A_crit_factor_N11 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 11 5 x * (2 * (1 - x) ^ 2 * (1 - x ^ 11) ^ 2) = -1 * (x - 1) ^ 4 * peval f11 x := by
  rw [N11_linkA5 x h0 h1, N11_crit_factor]
  simp only [peval]; ring

/-- Factorisation of the critical numerator for `N = 12` (note §6, "Structure of the numerator"):
`A_6 = 2(1−ξ)²(1−ξ^12)² a_6 = −2·(ξ−1)⁴ (ξ^6+1) f_12(ξ)`. -/
theorem A_crit_factor_N12 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 12 6 x * (2 * (1 - x) ^ 2 * (1 - x ^ 12) ^ 2)
      = -(2 * (x ^ 6 + 1)) * (x - 1) ^ 4 * peval f12 x := by
  rw [N12_linkA6 x h0 h1, N12_crit_factor]
  simp only [peval]; ring

/-- Factorisation of the critical numerator for `N = 13`:
`A_6 = −1·(ξ−1)⁴ f_13(ξ)`. -/
theorem A_crit_factor_N13 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 13 6 x * (2 * (1 - x) ^ 2 * (1 - x ^ 13) ^ 2) = -1 * (x - 1) ^ 4 * peval f13 x := by
  rw [N13_linkA6 x h0 h1, N13_crit_factor]
  simp only [peval]; ring

/-- Factorisation of the critical numerator for `N = 14` (note §6, "Structure of the numerator"):
`A_7 = 2(1−ξ)²(1−ξ^14)² a_7 = −1·(ξ−1)⁴ (ξ^7+1)(ξ+1) f_14(ξ)`. -/
theorem A_crit_factor_N14 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 14 7 x * (2 * (1 - x) ^ 2 * (1 - x ^ 14) ^ 2)
      = -(1 * (x ^ 7 + 1) * (x + 1)) * (x - 1) ^ 4 * peval f14 x := by
  rw [N14_linkA7 x h0 h1, N14_crit_factor]
  simp only [peval]; ring

/-- Factorisation of the critical numerator for `N = 15`:
`A_7 = −1·(ξ−1)⁴ f_15(ξ)`. -/
theorem A_crit_factor_N15 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 15 7 x * (2 * (1 - x) ^ 2 * (1 - x ^ 15) ^ 2) = -1 * (x - 1) ^ 4 * peval f15 x := by
  rw [N15_linkA7 x h0 h1, N15_crit_factor]
  simp only [peval]; ring

/-- Factorisation of the critical numerator for `N = 16` (note §6, "Structure of the numerator"):
`A_8 = 2(1−ξ)²(1−ξ^16)² a_8 = −2·(ξ−1)⁴ (ξ^8+1) f_16(ξ)`. -/
theorem A_crit_factor_N16 (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    aBlock 16 8 x * (2 * (1 - x) ^ 2 * (1 - x ^ 16) ^ 2)
      = -(2 * (x ^ 8 + 1)) * (x - 1) ^ 4 * peval f16 x := by
  rw [N16_linkA8 x h0 h1, N16_crit_factor]
  simp only [peval]; ring

/-- `f_8` is palindromic, of degree 10, with leading coefficient 1 (note Table 1). -/
theorem f8_data : f8.reverse = f8 ∧ f8.length = 11 ∧ f8.getLast? = some 1 := by
  refine ⟨?_, ?_, ?_⟩ <;> rfl

/-- `f_9` is palindromic, of degree 16, with leading coefficient 4 (note Table 1). -/
theorem f9_data : f9.reverse = f9 ∧ f9.length = 17 ∧ f9.getLast? = some 4 := by
  refine ⟨?_, ?_, ?_⟩ <;> rfl

/-- `f_10` is palindromic, of degree 12, with leading coefficient 7 (note Table 1). -/
theorem f10_data : f10.reverse = f10 ∧ f10.length = 13 ∧ f10.getLast? = some 7 := by
  refine ⟨?_, ?_, ?_⟩ <;> rfl

/-- `f_11` is palindromic, of degree 20, with leading coefficient 10 (note Table 1). -/
theorem f11_data : f11.reverse = f11 ∧ f11.length = 21 ∧ f11.getLast? = some 10 := by
  refine ⟨?_, ?_, ?_⟩ <;> rfl

/-- `f_12` is palindromic, of degree 16, with leading coefficient 7 (note Table 1). -/
theorem f12_data : f12.reverse = f12 ∧ f12.length = 17 ∧ f12.getLast? = some 7 := by
  refine ⟨?_, ?_, ?_⟩ <;> rfl

/-- `f_13` is palindromic, of degree 24, with leading coefficient 18 (note Table 1). -/
theorem f13_data : f13.reverse = f13 ∧ f13.length = 25 ∧ f13.getLast? = some 18 := by
  refine ⟨?_, ?_, ?_⟩ <;> rfl

/-- `f_14` is palindromic, of degree 18, with leading coefficient 23 (note Table 1). -/
theorem f14_data : f14.reverse = f14 ∧ f14.length = 19 ∧ f14.getLast? = some 23 := by
  refine ⟨?_, ?_, ?_⟩ <;> rfl

/-- `f_15` is palindromic, of degree 28, with leading coefficient 28 (note Table 1). -/
theorem f15_data : f15.reverse = f15 ∧ f15.length = 29 ∧ f15.getLast? = some 28 := by
  refine ⟨?_, ?_, ?_⟩ <;> rfl

/-- `f_16` is palindromic, of degree 22, with leading coefficient 17 (note Table 1). -/
theorem f16_data : f16.reverse = f16 ∧ f16.length = 23 ∧ f16.getLast? = some 17 := by
  refine ⟨?_, ?_, ?_⟩ <;> rfl

theorem krein_D_N9 : aBlock 9 4 ((1128710633 : ℝ) / 10000000000) * bBlock 9 4 ((1128710633 : ℝ) / 10000000000) - gBlock 9 4 ((1128710633 : ℝ) / 10000000000) ^ 2 < 0 := by
  norm_num [aBlock, bBlock, gBlock, W, cR]

/-- Odd-`N` Krein window for `N = 9` (note Prop 9, qualitative part): just above the spectral
threshold `ξ*` there is a point where the ring is spectrally stable (all `a_m, b_m > 0`) but the
critical block is indefinite (`a b − g² < 0`); hence `ξ_EM > ξ*`. -/
theorem krein_N9 : ∃ ξs : ℝ, peval f9 ξs = 0 ∧ ((11287106323660274043 : ℝ) / 100000000000000000000) < ξs ∧ ξs < ((2821776580915068511 : ℝ) / 25000000000000000000) ∧
    (∀ y, 0 < y → y < 1 → peval f9 y = 0 → y = ξs) ∧
    ∃ x : ℝ, ξs < x ∧ x < 1 ∧
      (∀ m, 1 ≤ m → m ≤ 8 → 0 < aBlock 9 m x ∧ 0 < bBlock 9 m x) ∧
      aBlock 9 4 x * bBlock 9 4 x - gBlock 9 4 x ^ 2 < 0 := by
  obtain ⟨ξs, hta, htb, hf0, huniq, hs⟩ := cert_N9
  have hw : ξs < ((1128710633 : ℝ) / 10000000000) := lt_trans htb (by norm_num)
  refine ⟨ξs, hf0, hta, htb, huniq, ((1128710633 : ℝ) / 10000000000), hw, by norm_num, fun m hm1 hm2 => ?_, krein_D_N9⟩
  have h := hs ((1128710633 : ℝ) / 10000000000) (by norm_num) (by norm_num)
  exact ⟨h.2.2 hw m hm1 hm2, h.1 m hm1 hm2⟩

theorem krein_D_N11 : aBlock 11 5 ((2046011769 : ℝ) / 10000000000) * bBlock 11 5 ((2046011769 : ℝ) / 10000000000) - gBlock 11 5 ((2046011769 : ℝ) / 10000000000) ^ 2 < 0 := by
  norm_num [aBlock, bBlock, gBlock, W, cR]

/-- Odd-`N` Krein window for `N = 11` (note Prop 9, qualitative part): just above the spectral
threshold `ξ*` there is a point where the ring is spectrally stable (all `a_m, b_m > 0`) but the
critical block is indefinite (`a b − g² < 0`); hence `ξ_EM > ξ*`. -/
theorem krein_N11 : ∃ ξs : ℝ, peval f11 ξs = 0 ∧ ((20460117688266929063 : ℝ) / 100000000000000000000) < ξs ∧ ξs < ((2557514711033366133 : ℝ) / 12500000000000000000) ∧
    (∀ y, 0 < y → y < 1 → peval f11 y = 0 → y = ξs) ∧
    ∃ x : ℝ, ξs < x ∧ x < 1 ∧
      (∀ m, 1 ≤ m → m ≤ 10 → 0 < aBlock 11 m x ∧ 0 < bBlock 11 m x) ∧
      aBlock 11 5 x * bBlock 11 5 x - gBlock 11 5 x ^ 2 < 0 := by
  obtain ⟨ξs, hta, htb, hf0, huniq, hs⟩ := cert_N11
  have hw : ξs < ((2046011769 : ℝ) / 10000000000) := lt_trans htb (by norm_num)
  refine ⟨ξs, hf0, hta, htb, huniq, ((2046011769 : ℝ) / 10000000000), hw, by norm_num, fun m hm1 hm2 => ?_, krein_D_N11⟩
  have h := hs ((2046011769 : ℝ) / 10000000000) (by norm_num) (by norm_num)
  exact ⟨h.2.2 hw m hm1 hm2, h.1 m hm1 hm2⟩

theorem krein_D_N13 : aBlock 13 6 ((2798613887 : ℝ) / 10000000000) * bBlock 13 6 ((2798613887 : ℝ) / 10000000000) - gBlock 13 6 ((2798613887 : ℝ) / 10000000000) ^ 2 < 0 := by
  norm_num [aBlock, bBlock, gBlock, W, cR]

/-- Odd-`N` Krein window for `N = 13` (note Prop 9, qualitative part): just above the spectral
threshold `ξ*` there is a point where the ring is spectrally stable (all `a_m, b_m > 0`) but the
critical block is indefinite (`a b − g² < 0`); hence `ξ_EM > ξ*`. -/
theorem krein_N13 : ∃ ξs : ℝ, peval f13 ξs = 0 ∧ ((13993069432195334351 : ℝ) / 50000000000000000000) < ξs ∧ ξs < ((27986138864390668703 : ℝ) / 100000000000000000000) ∧
    (∀ y, 0 < y → y < 1 → peval f13 y = 0 → y = ξs) ∧
    ∃ x : ℝ, ξs < x ∧ x < 1 ∧
      (∀ m, 1 ≤ m → m ≤ 12 → 0 < aBlock 13 m x ∧ 0 < bBlock 13 m x) ∧
      aBlock 13 6 x * bBlock 13 6 x - gBlock 13 6 x ^ 2 < 0 := by
  obtain ⟨ξs, hta, htb, hf0, huniq, hs⟩ := cert_N13
  have hw : ξs < ((2798613887 : ℝ) / 10000000000) := lt_trans htb (by norm_num)
  refine ⟨ξs, hf0, hta, htb, huniq, ((2798613887 : ℝ) / 10000000000), hw, by norm_num, fun m hm1 hm2 => ?_, krein_D_N13⟩
  have h := hs ((2798613887 : ℝ) / 10000000000) (by norm_num) (by norm_num)
  exact ⟨h.2.2 hw m hm1 hm2, h.1 m hm1 hm2⟩

theorem krein_D_N15 : aBlock 15 7 ((856146459 : ℝ) / 2500000000) * bBlock 15 7 ((856146459 : ℝ) / 2500000000) - gBlock 15 7 ((856146459 : ℝ) / 2500000000) ^ 2 < 0 := by
  norm_num [aBlock, bBlock, gBlock, W, cR]

/-- Odd-`N` Krein window for `N = 15` (note Prop 9, qualitative part): just above the spectral
threshold `ξ*` there is a point where the ring is spectrally stable (all `a_m, b_m > 0`) but the
critical block is indefinite (`a b − g² < 0`); hence `ξ_EM > ξ*`. -/
theorem krein_N15 : ∃ ξs : ℝ, peval f15 ξs = 0 ∧ ((1712292917560851283 : ℝ) / 5000000000000000000) < ξs ∧ ξs < ((34245858351217025661 : ℝ) / 100000000000000000000) ∧
    (∀ y, 0 < y → y < 1 → peval f15 y = 0 → y = ξs) ∧
    ∃ x : ℝ, ξs < x ∧ x < 1 ∧
      (∀ m, 1 ≤ m → m ≤ 14 → 0 < aBlock 15 m x ∧ 0 < bBlock 15 m x) ∧
      aBlock 15 7 x * bBlock 15 7 x - gBlock 15 7 x ^ 2 < 0 := by
  obtain ⟨ξs, hta, htb, hf0, huniq, hs⟩ := cert_N15
  have hw : ξs < ((856146459 : ℝ) / 2500000000) := lt_trans htb (by norm_num)
  refine ⟨ξs, hf0, hta, htb, huniq, ((856146459 : ℝ) / 2500000000), hw, by norm_num, fun m hm1 hm2 => ?_, krein_D_N15⟩
  have h := hs ((856146459 : ℝ) / 2500000000) (by norm_num) (by norm_num)
  exact ⟨h.2.2 hw m hm1 hm2, h.1 m hm1 hm2⟩


/-! ### `ξ*(8)`: an algebraic unit (note Cor 7) -/

open Polynomial in
/-- The polynomial of note Cor 7, as an element of `ℤ[X]`. -/
def P8 : ℤ[X] :=
  X ^ 10 - 12 * X ^ 9 - 37 * X ^ 8 - 72 * X ^ 7 - 84 * X ^ 6 - 88 * X ^ 5 - 84 * X ^ 4
    - 72 * X ^ 3 - 37 * X ^ 2 - 12 * X + 1

theorem P8_monic : P8.Monic := by unfold P8; monicity!

theorem P8_natDegree : P8.natDegree = 10 := by unfold P8; compute_degree!

theorem aeval_P8 (ξ : ℝ) : Polynomial.aeval ξ P8 = peval f8 ξ := by
  simp [P8, f8, peval, map_ofNat] <;> ring

theorem peval_f8_inv (ξ : ℝ) (h0 : ξ ≠ 0) : peval f8 ξ⁻¹ = peval f8 ξ / ξ ^ 10 := by
  simp only [f8, peval]; field_simp; push_cast; ring

/-- Note Cor 7: any real root of `P_8` (in particular `ξ*(8)`) is a nonzero algebraic integer whose
inverse is also an algebraic integer, i.e. an algebraic unit. -/
theorem f8_root_unit (ξ : ℝ) (h : peval f8 ξ = 0) :
    ξ ≠ 0 ∧ IsIntegral ℤ ξ ∧ IsIntegral ℤ ξ⁻¹ ∧ ξ * ξ⁻¹ = 1 := by
  have h0 : ξ ≠ 0 := by
    intro h0; rw [h0] at h; norm_num [f8, peval] at h
  refine ⟨h0, ⟨P8, P8_monic, ?_⟩, ⟨P8, P8_monic, ?_⟩, mul_inv_cancel₀ h0⟩
  · rw [← Polynomial.aeval_def, aeval_P8, h]
  · rw [← Polynomial.aeval_def, aeval_P8, peval_f8_inv ξ h0, h, zero_div]

/-- Note Cor 7: `C = ξ + 1/ξ` satisfies `C⁵ − 12C⁴ − 42C³ − 24C² + 32C + 32 = 0`. -/
theorem f8_C_poly (ξ : ℝ) (h : peval f8 ξ = 0) :
    (ξ + ξ⁻¹) ^ 5 - 12 * (ξ + ξ⁻¹) ^ 4 - 42 * (ξ + ξ⁻¹) ^ 3 - 24 * (ξ + ξ⁻¹) ^ 2
      + 32 * (ξ + ξ⁻¹) + 32 = 0 := by
  have h0 : ξ ≠ 0 := (f8_root_unit ξ h).1
  have key : (ξ + ξ⁻¹) ^ 5 - 12 * (ξ + ξ⁻¹) ^ 4 - 42 * (ξ + ξ⁻¹) ^ 3 - 24 * (ξ + ξ⁻¹) ^ 2
      + 32 * (ξ + ξ⁻¹) + 32 = peval f8 ξ / ξ ^ 5 := by
    simp only [f8, peval]; field_simp; push_cast; ring
  rw [key, h, zero_div]

/-- Note Cor 7 / Thm 6 for `N = 8`, packaged: the spectral threshold `ξ*(8)` is an algebraic unit,
root of the monic palindromic degree-10 `P_8`, and `ξ* + 1/ξ*` is a root of the quintic. -/
theorem xi8_unit : ∃ ξs : ℝ, ((3367811150204063113 : ℝ) / 50000000000000000000) < ξs ∧
    ξs < ((6735622300408126227 : ℝ) / 100000000000000000000) ∧
    Polynomial.aeval ξs P8 = 0 ∧ IsIntegral ℤ ξs ∧ IsIntegral ℤ ξs⁻¹ ∧
    (ξs + ξs⁻¹) ^ 5 - 12 * (ξs + ξs⁻¹) ^ 4 - 42 * (ξs + ξs⁻¹) ^ 3 - 24 * (ξs + ξs⁻¹) ^ 2
      + 32 * (ξs + ξs⁻¹) + 32 = 0 ∧
    ∀ x, 0 < x → x < 1 →
      (x < ξs → aBlock 8 4 x < 0) ∧ (ξs < x → ∀ m, 1 ≤ m → m ≤ 7 → 0 < aBlock 8 m x) := by
  obtain ⟨ξs, hta, htb, hf0, _, hs⟩ := cert_N8
  obtain ⟨_, hi, hi', _⟩ := f8_root_unit ξs hf0
  exact ⟨ξs, hta, htb, by rw [aeval_P8, hf0], hi, hi', f8_C_poly ξs hf0,
    fun x h0 h1 => ⟨(hs x h0 h1).2.1, (hs x h0 h1).2.2⟩⟩

/-! ### Relation to the hyperbolic radius -/

/-- With `ξ = tanh²(r/2)` (`r > 0`): `c = (1+ξ)/(1−ξ) = cosh r` and `ξ + 1/ξ = 2 + 4/sinh² r`
(note §3.1 and §5, "The natural invariant"). -/
theorem tanh_relations (r : ℝ) (hr : 0 < r) :
    cR (Real.tanh (r / 2) ^ 2) = Real.cosh r ∧
      Real.tanh (r / 2) ^ 2 + (Real.tanh (r / 2) ^ 2)⁻¹ = 2 + 4 / Real.sinh r ^ 2 := by
  have hs : 0 < Real.sinh (r / 2) := Real.sinh_pos_iff.mpr (by linarith)
  have hc : 0 < Real.cosh (r / 2) := Real.cosh_pos _
  have hid : Real.cosh (r / 2) ^ 2 - Real.sinh (r / 2) ^ 2 = 1 := Real.cosh_sq_sub_sinh_sq _
  have e1 : Real.cosh r = Real.cosh (r / 2) ^ 2 + Real.sinh (r / 2) ^ 2 := by
    rw [← Real.cosh_two_mul]; ring_nf
  have e2 : Real.sinh r = 2 * Real.sinh (r / 2) * Real.cosh (r / 2) := by
    rw [← Real.sinh_two_mul]; ring_nf
  rw [Real.tanh_eq_sinh_div_cosh, e1, e2]
  have hsc : Real.cosh (r / 2) ^ 2 - Real.sinh (r / 2) ^ 2 ≠ 0 := by rw [hid]; norm_num
  constructor
  · unfold cR
    rw [div_pow]
    field_simp
    nlinarith [hid]
  · rw [div_pow]
    field_simp
    nlinarith [hid]

end Havelock.Tanh

#print axioms Havelock.Tanh.W_hasSum
#print axioms Havelock.Tanh.aSum_eq_aBlock
#print axioms Havelock.Tanh.aBlock_recip
#print axioms Havelock.Tanh.aBlock_even
#print axioms Havelock.Tanh.numerA_recip
#print axioms Havelock.Tanh.charpoly_block
#print axioms Havelock.Tanh.tanh_relations
#print axioms Havelock.Tanh.stable_N2
#print axioms Havelock.Tanh.stable_N3
#print axioms Havelock.Tanh.stable_N4
#print axioms Havelock.Tanh.stable_N5
#print axioms Havelock.Tanh.stable_N6
#print axioms Havelock.Tanh.stable_N7
#print axioms Havelock.Tanh.cert_N8
#print axioms Havelock.Tanh.cert_N9
#print axioms Havelock.Tanh.cert_N10
#print axioms Havelock.Tanh.cert_N11
#print axioms Havelock.Tanh.cert_N12
#print axioms Havelock.Tanh.cert_N13
#print axioms Havelock.Tanh.cert_N14
#print axioms Havelock.Tanh.cert_N15
#print axioms Havelock.Tanh.cert_N16
#print axioms Havelock.Tanh.A_crit_factor_N8
#print axioms Havelock.Tanh.A_crit_factor_N9
#print axioms Havelock.Tanh.A_crit_factor_N10
#print axioms Havelock.Tanh.A_crit_factor_N11
#print axioms Havelock.Tanh.A_crit_factor_N12
#print axioms Havelock.Tanh.A_crit_factor_N13
#print axioms Havelock.Tanh.A_crit_factor_N14
#print axioms Havelock.Tanh.A_crit_factor_N15
#print axioms Havelock.Tanh.A_crit_factor_N16
#print axioms Havelock.Tanh.krein_N9
#print axioms Havelock.Tanh.krein_N11
#print axioms Havelock.Tanh.krein_N13
#print axioms Havelock.Tanh.krein_N15
#print axioms Havelock.Tanh.f8_root_unit
#print axioms Havelock.Tanh.xi8_unit
