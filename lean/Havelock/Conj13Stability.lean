/-
Havelock — Conjecture 13 paper, §"Nonlinear stability for a degenerate minimum": the analytic and
metric-topological core, in an abstract setting.

* `Conj13.LyapunovAbstract`: `lyapunov_on_section` (Lemma B4 on a metric space: a conserved function with a
  strict minimum modulo a compact orbit), `stable_modulo_orbit`, `stable_modulo_orbit_global` (Theorem 17 for
  a topological group acting continuously, given a slice map σ as in Lemma A).
* `Conj13.LyapunovAnalytic`: `lemma_B1` (quartic lower bound in chart coordinates), `lemma_B2` (absorption,
  AM–GM), `lemma_B3` (penalty function with a strict minimum), with `quadratic_bound`, `taylor_lower`,
  `quartic_form_coercive`.

Not formalised: the instantiation to the vortex ring — existence and global existence of the vortex flow,
conservation of H and J, Lemma A (the slice), the manifold chart on S/K and the slaving map (hypotheses
here), and the verification of (A1)–(A4) for the ring.
-/
import Havelock.Conj13.LyapunovAbstract
import Havelock.Conj13.LyapunovAnalytic
