import Havelock.Tanh.General
import Mathlib.FieldTheory.Finite.Extension
import Mathlib.RingTheory.Polynomial.GaussLemma

/-! # Irreducibility of integer polynomials via a reduction mod `p` (note §"Irreducibility")

A kernel-checkable certificate format for "`f ∈ ℤ[X]` is irreducible over `ℤ` and over `ℚ`".

* `toPoly l` is the polynomial with integer coefficient list `l` (constant term first), the same
  convention as `peval` (`aeval_toPoly`).
* Lifting (`irreducible_of_irreducible_map_zmod`): if `f ∈ ℤ[X]` is primitive, `p ∤ lead f` and
  `f mod p` is irreducible in `(ZMod p)[X]`, then `f` is irreducible in `ℤ[X]`, hence (Gauss) in `ℚ[X]`.
* Finite fields: an irreducible `q ∈ 𝔽_p[X]` of degree `d` divides `X^{p^d} − X`. Hence `h̄` of degree
  `n` is irreducible as soon as `gcd(h̄, X^{p^d} − X) = 1` for `1 ≤ d ≤ n/2`.
* Certificate (`Cert`, `checkAll`): `X^p ≡ m₁ (mod h̄)`, a table `m_j ≡ X^{jp}` (`j < n`) verified by
  explicit quotients, the Frobenius chain `r_d ≡ r_{d−1}^p = Σ_j r_{d−1,j} m_j (mod h̄)` (`r_0 = X`),
  and Bezout identities `u_d h̄ + v_d (r_d − X) = 1`; all identities are in `ℤ[X]` modulo `p`
  coefficientwise and are checked by a computable `Bool` function (`decide +kernel`).
-/

namespace Havelock.Tanh.Irred

open Polynomial

/-! ### Integer coefficient lists as polynomials -/

/-- The polynomial `Σ l[i] X^i` (constant term first) over a commutative ring `R`. -/
noncomputable def toPoly {R : Type*} [CommRing R] : List ℤ → R[X]
  | [] => 0
  | c :: l => C (c : R) + X * toPoly l

/-- Coefficientwise sum. -/
def padd : List ℤ → List ℤ → List ℤ
  | [], l => l
  | a :: l, [] => a :: l
  | a :: l, b :: m => (a + b) :: padd l m

/-- Scalar multiple. -/
def pscale (c : ℤ) : List ℤ → List ℤ
  | [] => []
  | a :: l => (c * a) :: pscale c l

/-- Schoolbook product. -/
def pmul : List ℤ → List ℤ → List ℤ
  | [], _ => []
  | a :: l, m => padd (pscale a m) (0 :: pmul l m)

/-- Difference. -/
def psub (a b : List ℤ) : List ℤ := padd a (pscale (-1) b)

/-- All entries divisible by `p`. -/
def allDvd (p : ℤ) : List ℤ → Bool
  | [] => true
  | a :: l => decide (a % p = 0) && allDvd p l

/-- The coefficient list of `X^n`. -/
def xpow (n : ℕ) : List ℤ := List.replicate n 0 ++ [1]

section Lists
variable {R : Type*} [CommRing R]

@[simp] theorem toPoly_nil : (toPoly [] : R[X]) = 0 := rfl

@[simp] theorem toPoly_cons (c : ℤ) (l : List ℤ) :
    (toPoly (c :: l) : R[X]) = C (c : R) + X * toPoly l := rfl

theorem toPoly_padd : ∀ a b : List ℤ, (toPoly (padd a b) : R[X]) = toPoly a + toPoly b
  | [], b => by simp [padd]
  | a :: l, [] => by simp [padd]
  | a :: l, b :: m => by
    simp only [padd, toPoly_cons, toPoly_padd l m]; push_cast; rw [C_add]; ring

theorem toPoly_pscale (c : ℤ) : ∀ a : List ℤ, (toPoly (pscale c a) : R[X]) = C (c : R) * toPoly a
  | [] => by simp [pscale]
  | a :: l => by
    simp only [pscale, toPoly_cons, toPoly_pscale c l]; push_cast; rw [C_mul]; ring

theorem toPoly_pmul : ∀ a b : List ℤ, (toPoly (pmul a b) : R[X]) = toPoly a * toPoly b
  | [], b => by simp [pmul]
  | a :: l, b => by
    simp only [pmul, toPoly_padd, toPoly_pscale, toPoly_cons, toPoly_pmul l b]; push_cast
    simp only [C_0]; ring

theorem toPoly_psub (a b : List ℤ) : (toPoly (psub a b) : R[X]) = toPoly a - toPoly b := by
  simp only [psub, toPoly_padd, toPoly_pscale]; push_cast; simp only [C_neg, C_1]; ring

theorem toPoly_xpow (n : ℕ) : (toPoly (xpow n) : R[X]) = X ^ n := by
  unfold xpow
  induction n with
  | zero => simp
  | succ n ih => rw [List.replicate_succ, List.cons_append, toPoly_cons, ih]; simp [pow_succ]; ring

theorem map_toPoly {S : Type*} [CommRing S] (φ : R →+* S) :
    ∀ l : List ℤ, (toPoly l : R[X]).map φ = toPoly l
  | [] => by simp
  | c :: l => by simp [Polynomial.map_add, Polynomial.map_mul, map_toPoly φ l]

theorem coeff_toPoly : ∀ (l : List ℤ) (i : ℕ), (toPoly l : R[X]).coeff i = ((l.getD i 0 : ℤ) : R)
  | [], i => by simp
  | c :: l, 0 => by simp
  | c :: l, i + 1 => by
    rw [toPoly_cons, coeff_add, coeff_C, coeff_X_mul, coeff_toPoly l i]; simp

theorem natDegree_toPoly_le : ∀ l : List ℤ, (toPoly l : R[X]).natDegree ≤ l.length - 1
  | [] => by simp
  | c :: l => by
    rw [natDegree_le_iff_coeff_eq_zero]
    intro i hi
    rw [coeff_toPoly, List.getD_eq_getElem?_getD, List.getElem?_eq_none (by simp at hi ⊢; omega)]
    simp

theorem aeval_toPoly_real (ξ : ℝ) : ∀ l : List ℤ, aeval ξ (toPoly l : ℤ[X]) = peval l ξ
  | [] => by simp [peval]
  | c :: l => by simp [peval, aeval_toPoly_real ξ l]

end Lists

theorem toPoly_eq_zero_of_allDvd (p : ℕ) : ∀ l : List ℤ, allDvd (p : ℤ) l = true →
    (toPoly l : (ZMod p)[X]) = 0
  | [] => fun _ => rfl
  | c :: l => by
    intro h
    simp only [allDvd, Bool.and_eq_true, decide_eq_true_eq] at h
    rw [toPoly_cons, toPoly_eq_zero_of_allDvd p l h.2,
      (ZMod.intCast_zmod_eq_zero_iff_dvd c p).mpr (Int.dvd_of_emod_eq_zero h.1)]
    simp

theorem eq_of_allDvd_psub (p : ℕ) {a b : List ℤ} (h : allDvd (p : ℤ) (psub a b) = true) :
    (toPoly a : (ZMod p)[X]) = toPoly b := by
  have := toPoly_eq_zero_of_allDvd p _ h
  rw [toPoly_psub, sub_eq_zero] at this
  exact this

/-! ### Lifting irreducibility from `ZMod p` to `ℤ` and `ℚ` -/

/-- If `f ∈ ℤ[X]` is primitive, its leading coefficient is nonzero mod `p`, and `f mod p` is
irreducible, then `f` is irreducible in `ℤ[X]`. -/
theorem irreducible_of_irreducible_map_zmod (p : ℕ) [Fact p.Prime] (f : ℤ[X])
    (hprim : f.IsPrimitive) (hlead : ((f.leadingCoeff : ℤ) : ZMod p) ≠ 0)
    (hirr : Irreducible (f.map (Int.castRingHom (ZMod p)))) : Irreducible f := by
  set φ := Int.castRingHom (ZMod p)
  refine ⟨fun hu => hirr.not_isUnit (hu.map (mapRingHom φ)), fun a b hab => ?_⟩
  have hf0 : f ≠ 0 := by rintro rfl; simp at hlead
  have ha0 : a ≠ 0 := by rintro rfl; simp at hab; exact hf0 hab
  have hb0 : b ≠ 0 := by rintro rfl; simp at hab; exact hf0 hab
  have hdeg : (f.map φ).natDegree = f.natDegree := natDegree_map_of_leadingCoeff_ne_zero φ hlead
  have hmul : f.map φ = a.map φ * b.map φ := by rw [hab, Polynomial.map_mul]
  have hfm0 : f.map φ ≠ 0 := hirr.ne_zero
  have ham0 : a.map φ ≠ 0 := by rintro h; rw [h, zero_mul] at hmul; exact hfm0 hmul
  have hbm0 : b.map φ ≠ 0 := by rintro h; rw [h, mul_zero] at hmul; exact hfm0 hmul
  have e1 : (f.map φ).natDegree = (a.map φ).natDegree + (b.map φ).natDegree := by
    rw [hmul, natDegree_mul ham0 hbm0]
  have e2 : f.natDegree = a.natDegree + b.natDegree := by rw [hab, natDegree_mul ha0 hb0]
  have la : (a.map φ).natDegree ≤ a.natDegree := natDegree_map_le
  have lb : (b.map φ).natDegree ≤ b.natDegree := natDegree_map_le
  have key : ∀ c : ℤ[X], c ∣ f → c.natDegree = 0 → IsUnit c := by
    intro c hc hc0
    rw [natDegree_eq_zero] at hc0
    obtain ⟨r, rfl⟩ := hc0
    exact isUnit_C.mpr (hprim r hc)
  rcases hirr.isUnit_or_isUnit hmul with hu | hu
  · left
    have : (a.map φ).natDegree = 0 := natDegree_eq_zero_of_isUnit hu
    exact key a ⟨b, hab⟩ (by omega)
  · right
    have : (b.map φ).natDegree = 0 := natDegree_eq_zero_of_isUnit hu
    exact key b ⟨a, by rw [hab, mul_comm]⟩ (by omega)

/-- A polynomial with two coprime coefficients is primitive. -/
theorem isPrimitive_of_coprime_coeff (f : ℤ[X]) (i j : ℕ) (h : Int.gcd (f.coeff i) (f.coeff j) = 1) :
    f.IsPrimitive := by
  intro r hr
  rw [C_dvd_iff_dvd_coeff] at hr
  exact (Int.isCoprime_iff_gcd_eq_one.mpr h).isUnit_of_dvd' (hr i) (hr j)

/-! ### Finite fields -/

/-- An irreducible polynomial of degree `d` over `𝔽_p` divides `X^{p^d} − X`. -/
theorem dvd_X_pow_sub_X (p : ℕ) [Fact p.Prime] {q : (ZMod p)[X]} (hq : Irreducible q) :
    q ∣ X ^ p ^ q.natDegree - X := by
  have := (hq.natDegree_dvd_iff_dvd_X_pow_card_pow_sub_X (n := q.natDegree)).mp dvd_rfl
  rwa [Nat.card_zmod] at this

/-- Frobenius on a `ZMod p`-algebra: `(aeval x P)^p = aeval (x^p) P`. -/
theorem aeval_pow_card (p : ℕ) [Fact p.Prime] {A : Type*} [CommRing A] [Algebra (ZMod p) A]
    (x : A) (P : (ZMod p)[X]) : (aeval x P) ^ p = aeval (x ^ p) P := by
  rw [← expand_aeval, ZMod.expand_card, map_pow]

/-- A polynomial of degree `n` over `𝔽_p` with no irreducible factor of degree `d ≤ D` dividing
`X^{p^d} − X`, where `n ≤ 2D + 1`, is irreducible. -/
theorem irreducible_of_no_small_factor (p : ℕ) [Fact p.Prime] (f : (ZMod p)[X]) (D : ℕ)
    (hpos : 0 < f.natDegree) (hD : f.natDegree / 2 ≤ D)
    (hno : ∀ d, 1 ≤ d → d ≤ D → ∀ q : (ZMod p)[X], q ∣ f → q ∣ X ^ p ^ d - X → IsUnit q) :
    Irreducible f := by
  have hf0 : f ≠ 0 := by rintro rfl; simp at hpos
  have hnu : ¬IsUnit f := fun hu => by
    rw [natDegree_eq_zero_of_isUnit hu] at hpos; exact lt_irrefl 0 hpos
  refine ⟨hnu, fun a b hab => ?_⟩
  by_contra hcon
  push Not at hcon
  obtain ⟨hua, hub⟩ := hcon
  have ha0 : a ≠ 0 := by rintro rfl; simp at hab; exact hf0 hab
  have hb0 : b ≠ 0 := by rintro rfl; simp at hab; exact hf0 hab
  have e : f.natDegree = a.natDegree + b.natDegree := by rw [hab, natDegree_mul ha0 hb0]
  -- one of the two factors has degree ≤ n/2
  have small : ∀ c : (ZMod p)[X], c ∣ f → c ≠ 0 → ¬IsUnit c → c.natDegree ≤ D → False := by
    intro c hc hc0 hcu hcd
    obtain ⟨q, hqi, hqc⟩ := WfDvdMonoid.exists_irreducible_factor hcu hc0
    have hqd : q.natDegree ≤ c.natDegree := natDegree_le_of_dvd hqc hc0
    exact hqi.not_isUnit (hno q.natDegree hqi.natDegree_pos (by omega) q (hqc.trans hc)
      (dvd_X_pow_sub_X p hqi))
  rcases le_total a.natDegree b.natDegree with h | h
  · exact small a ⟨b, hab⟩ ha0 hua (by omega)
  · exact small b ⟨a, by rw [hab, mul_comm]⟩ hb0 hub (by omega)

/-! ### The certificate checker -/

/-- `Σ_j r_j · ms_j` (requires `r.length ≤ ms.length`). -/
def lincomb : List (List ℤ) → List ℤ → List ℤ
  | _, [] => []
  | [], _ :: _ => []
  | m :: ms, c :: r => padd (pscale c m) (lincomb ms r)

section Quot
variable (p : ℕ) [Fact p.Prime] (f : (ZMod p)[X])

/-- Value of a coefficient list in `AdjoinRoot f`. -/
noncomputable abbrev φ (l : List ℤ) : AdjoinRoot f := AdjoinRoot.mk f (toPoly l)

/-- `TableOK y z ms` : `ms_j ≡ z * y^j (mod f)` for all `j`. -/
def TableOK (y : AdjoinRoot f) : AdjoinRoot f → List (List ℤ) → Prop
  | _, [] => True
  | z, m :: ms => φ p f m = z ∧ TableOK y (z * y) ms

theorem φ_cons (c : ℤ) (l : List ℤ) :
    φ p f (c :: l) = (c : AdjoinRoot f) + AdjoinRoot.root f * φ p f l := by
  simp [φ, map_add, map_mul, AdjoinRoot.mk_X]

theorem aeval_toPoly_cons (y : AdjoinRoot f) (c : ℤ) (l : List ℤ) :
    aeval y (toPoly (c :: l) : (ZMod p)[X]) =
      (c : AdjoinRoot f) + y * aeval y (toPoly l : (ZMod p)[X]) := by
  simp

theorem lincomb_spec (y : AdjoinRoot f) :
    ∀ (r : List ℤ) (ms : List (List ℤ)) (z : AdjoinRoot f), TableOK p f y z ms →
      r.length ≤ ms.length → φ p f (lincomb ms r) = z * aeval y (toPoly r : (ZMod p)[X])
  | [], ms, z, _, _ => by cases ms <;> simp [lincomb, φ]
  | c :: r, [], z, _, hl => by simp at hl
  | c :: r, m :: ms, z, ht, hl => by
    obtain ⟨hm, ht'⟩ := ht
    have ih := lincomb_spec y r ms (z * y) ht' (by simpa using hl)
    simp only [φ, lincomb, toPoly_padd, toPoly_pscale, map_add, map_mul] at ih hm ⊢
    rw [ih, hm, aeval_toPoly_cons p f]
    simp only [map_intCast]
    ring

theorem mul_of_check {a b c q h : List ℤ} (hf : f = toPoly h)
    (hc : allDvd (p : ℤ) (psub (pmul a b) (padd c (pmul q h))) = true) :
    φ p f a * φ p f b = φ p f c := by
  have e := eq_of_allDvd_psub p hc
  rw [toPoly_pmul, toPoly_padd, toPoly_pmul, ← hf] at e
  show AdjoinRoot.mk f (toPoly a) * AdjoinRoot.mk f (toPoly b) = AdjoinRoot.mk f (toPoly c)
  rw [← map_mul, e, map_add, map_mul, AdjoinRoot.mk_self, mul_zero, add_zero]

/-- Table check: `cur * m1 ≡ nxt (mod h)` with explicit quotients. -/
def checkTable (p : ℤ) (h m1 : List ℤ) : List ℤ → List (List ℤ × List ℤ) → Bool
  | _, [] => true
  | cur, (nxt, q) :: rest =>
    allDvd p (psub (pmul cur m1) (padd nxt (pmul q h))) && checkTable p h m1 nxt rest

theorem checkTable_spec {h m1 : List ℤ} (hf : f = toPoly h) (y : AdjoinRoot f)
    (hy : φ p f m1 = y) :
    ∀ (tab : List (List ℤ × List ℤ)) (cur : List ℤ) (z : AdjoinRoot f),
      checkTable p h m1 cur tab = true → φ p f cur = z → TableOK p f y (z * y) (tab.map Prod.fst)
  | [], _, _, _, _ => trivial
  | (nxt, q) :: rest, cur, z, hc, hz => by
    simp only [checkTable, Bool.and_eq_true] at hc
    have e := mul_of_check p f hf hc.1
    rw [hz, hy] at e
    exact ⟨e.symm, checkTable_spec hf y hy rest nxt (z * y) hc.2 e.symm⟩

/-- Frobenius chain + Bezout check. -/
def checkChain (p : ℤ) (h : List ℤ) (ms : List (List ℤ)) :
    List ℤ → List (List ℤ × List ℤ × List ℤ) → Bool
  | _, [] => true
  | cur, (nxt, u, v) :: rest =>
    decide (cur.length ≤ ms.length) && allDvd p (psub (lincomb ms cur) nxt) &&
    allDvd p (psub (padd (pmul u h) (pmul v (psub nxt [0, 1]))) [1]) &&
    checkChain p h ms nxt rest

theorem checkChain_spec {h : List ℤ} (hf : f = toPoly h) (ms : List (List ℤ))
    (hms : TableOK p f (AdjoinRoot.root f ^ p) 1 ms) :
    ∀ (ch : List (List ℤ × List ℤ × List ℤ)) (cur : List ℤ) (k : ℕ),
      checkChain p h ms cur ch = true → φ p f cur = AdjoinRoot.root f ^ p ^ k →
      ∀ i, i < ch.length → ∀ q : (ZMod p)[X], q ∣ f → q ∣ X ^ p ^ (k + i + 1) - X → IsUnit q
  | [], _, _, _, _, i, hi => by simp at hi
  | (nxt, u, v) :: rest, cur, k, hc, hk, i, hi => by
    simp only [checkChain, Bool.and_eq_true, decide_eq_true_eq] at hc
    obtain ⟨⟨⟨hlen, hfrob⟩, hbez⟩, hrest⟩ := hc
    -- the Frobenius step
    have hn : φ p f nxt = AdjoinRoot.root f ^ p ^ (k + 1) := by
      have e1 := eq_of_allDvd_psub p hfrob
      have e2 := lincomb_spec p f (AdjoinRoot.root f ^ p) cur ms 1 hms hlen
      simp only [φ] at e2 ⊢
      rw [← e1, e2, one_mul, ← aeval_pow_card, AdjoinRoot.aeval_eq]
      change φ p f cur ^ p = _
      rw [hk, ← pow_mul, pow_succ]
    cases i with
    | zero =>
      intro q hqf hqX
      have hb := eq_of_allDvd_psub p hbez
      rw [toPoly_padd, toPoly_pmul, toPoly_pmul, toPoly_psub, ← hf] at hb
      have hd : f ∣ X ^ p ^ (k + 1) - toPoly nxt := by
        rw [← AdjoinRoot.mk_eq_mk, map_pow, AdjoinRoot.mk_X]; exact hn.symm
      have h1 : q ∣ toPoly nxt - X := by
        have := dvd_sub hqX (hqf.trans hd)
        rwa [show (X ^ p ^ (k + 0 + 1) - X : (ZMod p)[X]) - (X ^ p ^ (k + 1) - toPoly nxt)
          = toPoly nxt - X by ring] at this
      have h2 : q ∣ (toPoly [1] : (ZMod p)[X]) := by
        rw [← hb]
        exact dvd_add (dvd_mul_of_dvd_right hqf _)
          (dvd_mul_of_dvd_right (by simpa using h1) _)
      exact isUnit_of_dvd_one (by simpa using h2)
    | succ i =>
      intro q hqf hqX
      exact checkChain_spec hf ms hms rest nxt (k + 1) hrest hn i (by simpa using hi) q hqf
        (by rwa [show k + 1 + i + 1 = k + (i + 1) + 1 by ring])

end Quot

/-- The certificate data for one polynomial and one prime. -/
structure Cert where
  /-- quotient in `X^p = m₁ + q₀ h (mod p)` -/
  q0 : List ℤ
  /-- `m₁ ≡ X^p (mod h)` -/
  m1 : List ℤ
  /-- `(m_j, q_j)`, `j = 2, …, n−1`, with `m_{j−1} m₁ = m_j + q_j h (mod p)` -/
  tab : List (List ℤ × List ℤ)
  /-- `(r_d, u_d, v_d)`, `d = 1, …, D`, with `r_d ≡ r_{d−1}^p (mod h)` and
  `u_d h + v_d (r_d − X) = 1 (mod p)` -/
  chain : List (List ℤ × List ℤ × List ℤ)
  /-- two coefficient indices with coprime coefficients (primitivity) -/
  i : ℕ
  /-- two coefficient indices with coprime coefficients (primitivity) -/
  j : ℕ

/-- The full certificate check. -/
def checkAll (p : ℕ) (h : List ℤ) (c : Cert) : Bool :=
  allDvd p (psub (xpow p) (padd c.m1 (pmul c.q0 h))) &&
  checkTable p h c.m1 c.m1 c.tab &&
  checkChain p h ([1] :: c.m1 :: c.tab.map Prod.fst) [0, 1] c.chain &&
  decide (2 ≤ h.length) &&
  decide ((h.length - 1) / 2 ≤ c.chain.length) &&
  decide (h.getLastD 0 % (p : ℤ) ≠ 0) &&
  decide (Int.gcd (h.getD c.i 0) (h.getD c.j 0) = 1)

theorem natDegree_toPoly_eq {R : Type*} [CommRing R] (h : List ℤ)
    (hlead : ((h.getLastD 0 : ℤ) : R) ≠ 0) :
    (toPoly h : R[X]).natDegree = h.length - 1 := by
  apply le_antisymm (natDegree_toPoly_le h)
  apply le_natDegree_of_ne_zero
  rw [coeff_toPoly]
  convert hlead using 3
  rw [List.getLastD_eq_getLast?, List.getLast?_eq_getElem?, List.getD_eq_getElem?_getD]

theorem lead_toPoly (h : List ℤ) (h0 : h.getLastD 0 ≠ 0) :
    (toPoly h : ℤ[X]).leadingCoeff = h.getLastD 0 := by
  rw [leadingCoeff, natDegree_toPoly_eq h (by simpa using h0), coeff_toPoly]
  simp [List.getLastD_eq_getLast?, List.getLast?_eq_getElem?, List.getD_eq_getElem?_getD]

/-- Soundness of the certificate: `h mod p` is irreducible over `𝔽_p`. -/
theorem irreducible_zmod_of_checkAll (p : ℕ) [Fact p.Prime] (h : List ℤ) (c : Cert)
    (hc : checkAll p h c = true) : Irreducible (toPoly h : (ZMod p)[X]) := by
  simp only [checkAll, Bool.and_eq_true, decide_eq_true_eq] at hc
  obtain ⟨⟨⟨⟨⟨⟨hX, htab⟩, hch⟩, hl⟩, hD⟩, hlead⟩, -⟩ := hc
  set f : (ZMod p)[X] := toPoly h with hf
  have hlead' : ((h.getLastD 0 : ℤ) : ZMod p) ≠ 0 := by
    rw [ne_eq, ZMod.intCast_zmod_eq_zero_iff_dvd]
    exact fun hd => hlead (Int.emod_eq_zero_of_dvd hd)
  have hdeg : f.natDegree = h.length - 1 := natDegree_toPoly_eq (R := ZMod p) h hlead'
  set y := AdjoinRoot.root f ^ p
  have hm1 : φ p f c.m1 = y := by
    have e := eq_of_allDvd_psub p hX
    rw [toPoly_xpow, toPoly_padd, toPoly_pmul, ← hf] at e
    simp only [φ, y, ← AdjoinRoot.mk_X, ← map_pow, e, map_add, map_mul, AdjoinRoot.mk_self,
      mul_zero, add_zero]
  have hT : TableOK p f y 1 ([1] :: c.m1 :: c.tab.map Prod.fst) := by
    refine ⟨by simp [φ], by rw [hm1, one_mul], ?_⟩
    rw [one_mul]
    exact checkTable_spec p f hf y hm1 c.tab c.m1 y htab hm1
  refine irreducible_of_no_small_factor p f c.chain.length (by omega) (by omega) ?_
  intro d hd1 hd2 q hqf hqX
  have := checkChain_spec p f hf _ hT c.chain [0, 1] 0 hch (by simp [φ, AdjoinRoot.mk_X])
    (d - 1) (by omega) q hqf
  rw [show 0 + (d - 1) + 1 = d by omega] at this
  exact this hqX

/-- Main criterion: a certificate proves `toPoly h` irreducible over `ℤ` and over `ℚ`. -/
theorem irreducible_of_checkAll (p : ℕ) (hp : p.Prime) (h : List ℤ) (c : Cert)
    (hc : checkAll p h c = true) :
    Irreducible (toPoly h : ℤ[X]) ∧ Irreducible (toPoly h : ℚ[X]) := by
  have := Fact.mk hp
  have hirr := irreducible_zmod_of_checkAll p h c hc
  simp only [checkAll, Bool.and_eq_true, decide_eq_true_eq] at hc
  obtain ⟨⟨⟨⟨⟨⟨-, -⟩, -⟩, hl⟩, -⟩, hlead⟩, hgcd⟩ := hc
  have h0 : h.getLastD 0 ≠ 0 := fun e => hlead (by rw [e]; rfl)
  have hprim : (toPoly h : ℤ[X]).IsPrimitive :=
    isPrimitive_of_coprime_coeff _ c.i c.j (by rw [coeff_toPoly, coeff_toPoly]; simpa using hgcd)
  have hZ : Irreducible (toPoly h : ℤ[X]) := by
    refine irreducible_of_irreducible_map_zmod p _ hprim ?_ (by rwa [map_toPoly])
    rw [lead_toPoly h h0, ne_eq, ZMod.intCast_zmod_eq_zero_iff_dvd]
    exact fun hd => hlead (Int.emod_eq_zero_of_dvd hd)
  refine ⟨hZ, ?_⟩
  have := (IsPrimitive.Int.irreducible_iff_irreducible_map_cast hprim).mp hZ
  rwa [map_toPoly] at this

end Havelock.Tanh.Irred
