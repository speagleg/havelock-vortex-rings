# integrable/ — is the large-spacing vortex row an integrable lattice?

Follow-up to ../quantum/ section 11, which found finite-time Lyapunov exponents of the Toda-limit row decaying like
those of the integrable Toda lattice. **Verdict so far: not established.**
- The row's motion is regular in every test at small and moderate amplitude.
- But a control with no special structure behaves the same way.
- No extra integral of low polynomial or rational degree exists.
- Exact integrability is proved only for a different member of the family, where the chain is the Volterra lattice.

| Script | What it does |
|---|---|
| `find_integral.py` | sympy: degree-1 integrals sum g(Y_j, Y_{j+1}) m e^{-D_j}. Only H itself (m constant) |
| `search_deg2.py`, `search_deg2b.py` | Linear-algebra search (random points, complex-step brackets) for a local degree-2 integral on an open 6-site chain, polynomial x s-power ansatz, function-level dependencies removed. None found |
| `reduced_check.py` | The closed reduced dynamics in bond variables (checked to 3e-17) |
| `search_reduced.py` | The same search in reduced variables with a rational ansatz (denominators (1+Y^2)^a). None found |
| `many_ics.py` | Lyapunov exponents for 20 random initial conditions (N = 3, 4) -> `many_ics.out` |
| `poincare.py` | N = 3 Poincare sections: Toda-limit row vs full tanh row -> `poincare.png` |
| `phi_family.py`, `phi_family2.py` | Lyapunov tests of H = sum phi(Y_j) e^{-D_j} for several phi -> `phi_family.out`, `phi_family2.out` |
| `poincare_phi.py`, `poincare_strong.py` | Poincare sections across the phi family; large-amplitude attempt -> `poincare_phi.png`, `poincare_strong.png` |

**1. Reduced dynamics** (DERIVED, checked). With c_j = s_j s_{j+1} e^{-(u_{j+1}-u_j)}:
    dc_j/dt = 4 c_j (f(Y_{j+1}) c_{j+1} - f(Y_j) c_{j-1}),   dY_j/dt = 2 (c_{j-1} - c_j),   f(Y) = Y/(1+Y^2).
- The Y equation is exactly Toda's (Flaschka variables with b = -Y, a^2 = c); the c equation is deformed by f.
- H has degree 1 in the c's and the flow raises degree by exactly 1, so every integral is homogeneous in c.
- The only degree-0 integral is the impulse sum Y_j (exact). A periodic ring can mix degrees through
  prod c_j = e^{-L}/prod(1 + Y_j^2).

**2. Gauge family and an exact Volterra member** (DERIVED).
- *Gauge family.* The canonical shift u_j -> u_j + g(Y_j) maps the row to H = sum phi(Y_j) e^{-D_j} with
  phi = 2/(1+Y^2).
- *Volterra member.* For phi = e^{cY}, the bond variables a_j = e^{c Y_j - D_j} satisfy {a_j, a_{j+1}} = -c a_j a_{j+1}.
  This is the Volterra (Kac-van Moerbeke) Poisson algebra, so H = sum a_j is the integrable Volterra lattice.
- *Other phi.* The bracket picks up a factor (ln phi)'(Y_{j+1}), constant only for exponentials, so the vortex row is
  not Volterra.

**3. Integral searches** (NUMERICAL, negative).
- *Degree 1:* exactly ruled out (sympy). The adjacent-bond condition forces (1+Y^2)M'' - 2Y M' = 0, so M is constant.
- *Degree 2:* no local integral, for polynomial x s-power ansatze up to degree 4, or rational ansatze with
  denominators (1+Y^2)^{0..2}, on open and periodic chains.
- *Not ruled out:* higher degree, non-local, non-rational, or ring-specific integrals.

**4. Regularity is real but not specific** (NUMERICAL).
- *The row is regular.* Twenty random initial conditions (N = 3, 4) all give finite-time exponents at the level of
  the periodic Toda control (`many_ics.out`), and N = 3 Poincare sections show clean invariant curves (`poincare.png`).
  The full tanh row is mixed, with chaotic trajectories.
- *So are other decaying phi.* In the phi family:
  - phi = e^Y (Volterra, exactly integrable), 2/(1+Y^2), 1/(Y^2 + 0.6Y + 2), 1/Y^2, 1/(Y^2-1) and the control
    1/(1+Y^2)^2 all look regular;
  - 1 + Y^2, cosh Y and a random phi are chaotic (`phi_family*.out`).
  - The control has no structure that would make it integrable, yet its Poincare section also shows curves
    (`poincare_phi.png`).
  - So regularity here is generic for decaying phi in this energy range and is not evidence of exact integrability.
- *Large amplitude doesn't discriminate.* There the chains come apart rather than returning to the section
  (`poincare_strong.out`).

**Status of the earlier claims.**
- ../quantum/ section 11 said the Toda-limit row "looks integrable". That is downgraded to "regular at small amplitude,
  like other decaying members of its gauge family; exact integrability unproven".
- The -ln cosh(d/2) chain's resemblance to the SU(1,1) Ishimori-Haldane lattice Landau-Lifshitz chain is a structural
  observation. Its regular-looking Lyapunov data need the same control comparison before carrying weight.

**What would settle it.** A Lax pair or an extra integral (ring-specific, or of degree >= 3). Alternatively, a
singularity-analysis (Painleve) test of the reduced equations. Or long, high-precision runs looking for thin chaotic
layers near separatrices, which near-integrable systems have and integrable ones do not.
