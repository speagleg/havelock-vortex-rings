/-
Havelock — Conjecture 13: the abstract Lyapunov argument of Section "Nonlinear stability for a
degenerate minimum" (paper `conj13/paper/main.tex`, `sec:nonlinear`), in a metric-space setting.

This module formalises the *topological* half of the section:

* `lyapunov_on_section` (paper Lemma B4): a function `V` which is conserved along trajectories in a
  cross-section `S`, is constant on a compact set `O` (the relative-equilibrium orbit `K·x₀`) and is
  strictly larger than that constant on a punctured neighbourhood of `O` in `S` (paper Lemma B3,
  phrased intrinsically) makes `O` Lyapunov stable for the trajectories in `S`.
* `stable_modulo_orbit` (paper Thm `thm:abstract`, "Theorem 17"): a continuous map `σ` into a
  topological group `G` acting continuously on the phase space, with `σ = 1` on `O` and
  `σ x • x ∈ S` (paper Lemma A), together with `G`-equivariance of the trajectories, upgrades this to
  Lyapunov stability of `O` in the whole phase space, with no restriction on the perturbation.
* `stable_modulo_orbit_global`: the same statement for globally defined trajectories.

Abstraction (see the report for the gap to the vortex ring): the phase space is any metric space
`P`; the Hamiltonian flow is replaced by an arbitrary family of trajectories `c x : ℝ → P` defined on
preconnected time sets `I x ∋ 0` (e.g. the maximal existence interval); "symplectic + `H`-invariant"
is replaced by the consequence actually used in the paper, namely that `G` maps trajectories to
trajectories. `infDist x O` is the paper's `dist(x, K·x₀)`.
-/
import Mathlib.Topology.MetricSpace.Thickening
import Mathlib.Topology.Algebra.MulAction
import Mathlib.Topology.Algebra.Group.Defs
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Topology.Order.Compact

open Metric Set Filter Topology

namespace Havelock.Conj13.Lyapunov

variable {P : Type*} [MetricSpace P]

/-- **Paper Lemma B4 (Lyapunov stability on the cross-section `S`).**

`O` compact and nonempty (the orbit `K·x₀`), `V` continuous on `S`, equal to `V₀` on `O`, and
`V > V₀` at every point of `S` at positive distance `≤ r` from `O` (the strict minimum modulo `K`
of paper Lemma B3); the `r`-neighbourhood of `O` in `S` is compact (local compactness of `S`).
Trajectories `c x` starting in `S` are continuous on a preconnected time set `I x ∋ 0`, stay in `S`,
and conserve `V`. Then for every `ε > 0` there is `δ > 0` such that trajectories starting in `S`
within `δ` of `O` stay within `ε` of `O` for all `t ∈ I x`. -/
theorem lyapunov_on_section (S O : Set P) (hO : IsCompact O) (hOne : O.Nonempty)
    (V : P → ℝ) (V₀ : ℝ) (hVc : ContinuousOn V S) (hVO : ∀ x ∈ O, V x = V₀)
    (r : ℝ) (hr : 0 < r) (hK : IsCompact {x | x ∈ S ∧ infDist x O ≤ r})
    (hmin : ∀ x ∈ S, infDist x O ≤ r → 0 < infDist x O → V₀ < V x)
    (I : P → Set ℝ) (c : P → ℝ → P)
    (hI0 : ∀ x ∈ S, (0 : ℝ) ∈ I x) (hIc : ∀ x ∈ S, IsPreconnected (I x))
    (hc0 : ∀ x ∈ S, c x 0 = x) (hcc : ∀ x ∈ S, ContinuousOn (c x) (I x))
    (hcS : ∀ x ∈ S, ∀ t ∈ I x, c x t ∈ S) (hcV : ∀ x ∈ S, ∀ t ∈ I x, V (c x t) = V x) :
    ∀ ε > 0, ∃ δ > 0, ∀ x ∈ S, infDist x O < δ → ∀ t ∈ I x, infDist (c x t) O < ε := by
  intro ε hε
  set Kr := {x | x ∈ S ∧ infDist x O ≤ r} with hKr
  have hKrc : IsClosed Kr := hK.isClosed
  set e := min ε r / 2 with he
  have he0 : 0 < e := by positivity
  have heε : e < ε := by
    have := min_le_left ε r; linarith
  have her : e < r := by
    have := min_le_right ε r; linarith
  -- the sphere `Σ = {x ∈ S | dist(x, O) = e}` and the positive minimum `m` of `V - V₀` on it
  set Sig := Kr ∩ {x | infDist x O = e} with hSig
  have hSigc : IsCompact Sig :=
    hK.inter_right (isClosed_eq (continuous_infDist_pt O) continuous_const)
  have hm : ∃ m > 0, ∀ x ∈ Sig, V₀ + m ≤ V x := by
    rcases Sig.eq_empty_or_nonempty with h | h
    · exact ⟨1, one_pos, by simp [h]⟩
    · obtain ⟨y, hy, hymin⟩ := hSigc.exists_isMinOn h (hVc.mono fun x hx => hx.1.1)
      have hy' : V₀ < V y := hmin y hy.1.1 hy.1.2 (by rw [hy.2]; exact he0)
      refine ⟨V y - V₀, by linarith, fun x hx => ?_⟩
      have := hymin hx
      simp only [mem_ofPred_eq] at this
      linarith
  obtain ⟨m, hm0, hmSig⟩ := hm
  -- the set `C` of points of `Kr` where `V ≥ V₀ + m` is compact and disjoint from `O`
  set C := Kr ∩ V ⁻¹' Ici (V₀ + m) with hC
  have hCc : IsCompact C :=
    hK.of_isClosed_subset ((hVc.mono fun x hx => hx.1).preimage_isClosed_of_isClosed hKrc
      isClosed_Ici) inter_subset_left
  have hδ0 : ∃ δ₀ > 0, ∀ x ∈ C, δ₀ ≤ infDist x O := by
    rcases C.eq_empty_or_nonempty with h | h
    · exact ⟨1, one_pos, by simp [h]⟩
    · obtain ⟨y, hy, hymin⟩ :=
        hCc.exists_isMinOn h (continuous_infDist_pt O).continuousOn
      refine ⟨infDist y O, ?_, fun x hx => hymin hx⟩
      rcases (infDist_nonneg (x := y) (s := O)).lt_or_eq with h0 | h0
      · exact h0
      · exfalso
        have hyO : y ∈ O := (hO.isClosed.mem_iff_infDist_zero hOne).2 h0.symm
        have h1 : V₀ + m ≤ V y := hy.2
        rw [hVO y hyO] at h1
        linarith
  obtain ⟨δ₀, hδ₀, hδ₀C⟩ := hδ0
  refine ⟨min δ₀ e, lt_min hδ₀ he0, fun x hxS hx t ht => ?_⟩
  have hxe : infDist x O < e := lt_of_lt_of_le hx (min_le_right _ _)
  have hVx : V x < V₀ + m := by
    by_contra hcon
    have hxC : x ∈ C := ⟨⟨hxS, by linarith⟩, le_of_not_gt hcon⟩
    have := hδ₀C x hxC
    have := min_le_left δ₀ e
    linarith
  by_contra hcon
  replace hcon := le_of_not_gt hcon
  -- intermediate value theorem on the preconnected time set
  have hcont : ContinuousOn (fun s => infDist (c x s) O) (I x) :=
    (continuous_infDist_pt O).comp_continuousOn (hcc x hxS)
  have hmem : e ∈ Icc (infDist (c x 0) O) (infDist (c x t) O) := by
    rw [hc0 x hxS]; exact ⟨hxe.le, by linarith⟩
  obtain ⟨s, hs, hse⟩ := (hIc x hxS).intermediate_value (hI0 x hxS) ht hcont hmem
  have hsS : c x s ∈ S := hcS x hxS s hs
  have hsSig : c x s ∈ Sig := ⟨⟨hsS, by simp only at hse; rw [hse]; exact her.le⟩, hse⟩
  have := hmSig _ hsSig
  rw [hcV x hxS s hs] at this
  linarith

variable {G : Type*} [Group G] [TopologicalSpace G] [IsTopologicalGroup G] [MulAction G P]
  [ContinuousSMul G P]

/-- **Paper Theorem `thm:abstract` ("Theorem 17"), abstract form: Lyapunov stability modulo the
orbit in the full phase space, with no restriction on the momentum of the perturbation.**

Hypotheses of `lyapunov_on_section` on the cross-section `S` (paper Lemmas B3, B4), plus:
* (paper Lemma A) an open `U ⊇ O` and a map `σ : P → G`, continuous on `U`, with `σ = 1` on `O`
  and `σ x • x ∈ S` for `x ∈ U`;
* (`G` commutes with the flow) for every `g` and `x`, `I (g • x) = I x` and
  `c (g • x) t = g • c x t`;
* every trajectory is continuous on its preconnected time set `I x ∋ 0`, starts at `x`, and
  trajectories starting in `S` stay in `S` and conserve `V`.

Conclusion: for every `ε > 0` there is `δ > 0` such that every trajectory with
`dist(x(0), O) < δ` satisfies `dist(x(t), O) < ε` for all `t ∈ I x`. -/
theorem stable_modulo_orbit (S O : Set P) (hO : IsCompact O) (hOne : O.Nonempty)
    (V : P → ℝ) (V₀ : ℝ) (hVc : ContinuousOn V S) (hVO : ∀ x ∈ O, V x = V₀)
    (r : ℝ) (hr : 0 < r) (hK : IsCompact {x | x ∈ S ∧ infDist x O ≤ r})
    (hmin : ∀ x ∈ S, infDist x O ≤ r → 0 < infDist x O → V₀ < V x)
    (U : Set P) (hU : IsOpen U) (hOU : O ⊆ U) (σ : P → G) (hσc : ContinuousOn σ U)
    (hσO : ∀ x ∈ O, σ x = 1) (hσS : ∀ x ∈ U, σ x • x ∈ S)
    (I : P → Set ℝ) (c : P → ℝ → P)
    (hI0 : ∀ x, (0 : ℝ) ∈ I x) (hIc : ∀ x, IsPreconnected (I x))
    (hc0 : ∀ x, c x 0 = x) (hcc : ∀ x, ContinuousOn (c x) (I x))
    (hcS : ∀ x ∈ S, ∀ t ∈ I x, c x t ∈ S) (hcV : ∀ x ∈ S, ∀ t ∈ I x, V (c x t) = V x)
    (hIG : ∀ (g : G) x, I (g • x) = I x) (hcG : ∀ (g : G) x, ∀ t ∈ I x, c (g • x) t = g • c x t) :
    ∀ ε > 0, ∃ δ > 0, ∀ x, infDist x O < δ → ∀ t ∈ I x, infDist (c x t) O < ε := by
  intro ε hε
  set Kr := {x | x ∈ S ∧ infDist x O ≤ r} with hKr
  set ε₁ := min (ε / 2) r with hε₁
  have hε₁0 : 0 < ε₁ := lt_min (by linarith) hr
  -- Lemma B4 on the section with `ε₁`
  obtain ⟨δ₁, hδ₁, hB4⟩ := lyapunov_on_section S O hO hOne V V₀ hVc hVO r hr hK hmin I c
    (fun x _ => hI0 x) (fun x _ => hIc x) (fun x _ => hc0 x) (fun x _ => hcc x) hcS hcV
    ε₁ hε₁0
  -- `ω(g) = sup_{y ∈ Kr} dist(g⁻¹ • y, y) < ε/2` for `g` near `1` (tube lemma)
  have hn1 : IsOpen {p : G × P | dist (p.1⁻¹ • p.2) p.2 < ε / 2} :=
    isOpen_lt ((continuous_fst.inv.smul continuous_snd).dist continuous_snd) continuous_const
  obtain ⟨u₁, v₁, hu₁, -, h1u₁, hKv₁, huv₁⟩ := generalized_tube_lemma (isCompact_singleton (x := (1 : G))) hK hn1
    (by
      rintro ⟨g, y⟩ ⟨hg, -⟩
      rw [mem_singleton_iff] at hg
      subst hg
      show dist ((1 : G)⁻¹ • y) y < ε / 2
      rw [inv_one, one_smul, dist_self]; positivity)
  -- `dist(g • x, x) < δ₁/2` for `g` near `1` and `x` near `O` (tube lemma)
  have hn2 : IsOpen {p : G × P | dist (p.1 • p.2) p.2 < δ₁ / 2} :=
    isOpen_lt ((continuous_fst.smul continuous_snd).dist continuous_snd) continuous_const
  obtain ⟨u₂, v₂, hu₂, hv₂, h1u₂, hOv₂, huv₂⟩ := generalized_tube_lemma (isCompact_singleton (x := (1 : G))) hO hn2
    (by
      rintro ⟨g, y⟩ ⟨hg, -⟩
      rw [mem_singleton_iff] at hg
      subst hg
      show dist ((1 : G) • y) y < δ₁ / 2
      rw [one_smul, dist_self]; positivity)
  have h1u₁' : (1 : G) ∈ u₁ := h1u₁ rfl
  have h1u₂' : (1 : G) ∈ u₂ := h1u₂ rfl
  -- an open neighbourhood `W` of `O` on which `σ ∈ u₁ ∩ u₂`
  set W := v₂ ∩ (U ∩ σ ⁻¹' (u₁ ∩ u₂)) with hW
  have hWo : IsOpen W := hv₂.inter (hσc.isOpen_inter_preimage hU (hu₁.inter hu₂))
  have hOW : O ⊆ W := fun y hy => ⟨hOv₂ hy, hOU hy, by
    simp only [mem_preimage, hσO y hy]; exact ⟨h1u₁', h1u₂'⟩⟩
  obtain ⟨δ₂, hδ₂, hδ₂W⟩ := hO.exists_cthickening_subset_open hWo hOW
  refine ⟨min δ₂ (δ₁ / 2), lt_min hδ₂ (by linarith), fun x hx t ht => ?_⟩
  obtain ⟨z, hzO, hxz⟩ := (infDist_lt_iff hOne).1 hx
  have hxW : x ∈ W := hδ₂W (mem_cthickening_of_dist_le x z δ₂ O hzO
    (hxz.le.trans (min_le_left _ _)))
  obtain ⟨hxv₂, hxU, hσu⟩ := hxW
  set g := σ x with hg
  set y₀ := g • x with hy₀
  have hy₀S : y₀ ∈ S := hσS x hxU
  have hdy₀ : dist y₀ x < δ₁ / 2 := huv₂ (mk_mem_prod hσu.2 hxv₂)
  have hy₀O : infDist y₀ O < δ₁ := by
    have h1 := infDist_le_infDist_add_dist (s := O) (x := y₀) (y := x)
    have h2 : infDist x O < δ₁ / 2 := lt_of_lt_of_le hx (min_le_right _ _)
    linarith
  -- the trajectory of `y₀ = σ(x)•x` stays near `O`, and `x(t) = σ(x)⁻¹ • y(t)`
  have hty : t ∈ I y₀ := by rw [hy₀, hIG]; exact ht
  have hyt := hB4 y₀ hy₀S hy₀O t hty
  have hytK : c y₀ t ∈ Kr := ⟨hcS y₀ hy₀S t hty, hyt.le.trans (min_le_right _ _)⟩
  have hxt : c x t = g⁻¹ • c y₀ t := by
    rw [hy₀, hcG g x t ht, inv_smul_smul]
  have hω : dist (g⁻¹ • c y₀ t) (c y₀ t) < ε / 2 := huv₁ (mk_mem_prod hσu.1 (hKv₁ hytK))
  have h1 := infDist_le_infDist_add_dist (s := O) (x := c x t) (y := c y₀ t)
  rw [hxt] at h1 ⊢
  have : ε₁ ≤ ε / 2 := min_le_left _ _
  linarith

/-- **Theorem 17 with global existence.** If every trajectory is defined for all `t ∈ ℝ`
(`I x = univ`), the conclusion of `stable_modulo_orbit` holds for all times — the exact form of
paper Thm `thm:abstract` / Thm `thm:main` (minus the instantiation). -/
theorem stable_modulo_orbit_global (S O : Set P) (hO : IsCompact O) (hOne : O.Nonempty)
    (V : P → ℝ) (V₀ : ℝ) (hVc : ContinuousOn V S) (hVO : ∀ x ∈ O, V x = V₀)
    (r : ℝ) (hr : 0 < r) (hK : IsCompact {x | x ∈ S ∧ infDist x O ≤ r})
    (hmin : ∀ x ∈ S, infDist x O ≤ r → 0 < infDist x O → V₀ < V x)
    (U : Set P) (hU : IsOpen U) (hOU : O ⊆ U) (σ : P → G) (hσc : ContinuousOn σ U)
    (hσO : ∀ x ∈ O, σ x = 1) (hσS : ∀ x ∈ U, σ x • x ∈ S)
    (c : P → ℝ → P) (hc0 : ∀ x, c x 0 = x) (hcc : ∀ x, Continuous (c x))
    (hcS : ∀ x ∈ S, ∀ t, c x t ∈ S) (hcV : ∀ x ∈ S, ∀ t, V (c x t) = V x)
    (hcG : ∀ (g : G) x t, c (g • x) t = g • c x t) :
    ∀ ε > 0, ∃ δ > 0, ∀ x, infDist x O < δ → ∀ t : ℝ, infDist (c x t) O < ε := by
  intro ε hε
  obtain ⟨δ, hδ, h⟩ := stable_modulo_orbit S O hO hOne V V₀ hVc hVO r hr hK hmin U hU hOU σ
    hσc hσO hσS (fun _ => univ) c (fun _ => mem_univ _) (fun _ => isPreconnected_univ) hc0
    (fun x => (hcc x).continuousOn) (fun x hx t _ => hcS x hx t) (fun x hx t _ => hcV x hx t)
    (fun _ _ => rfl) (fun g x t _ => hcG g x t) ε hε
  exact ⟨δ, hδ, fun x hx t => h x hx t (mem_univ t)⟩

end Havelock.Conj13.Lyapunov
