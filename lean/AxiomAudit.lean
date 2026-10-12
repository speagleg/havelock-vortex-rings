/-
Axiom audit for CI (`.github/workflows/lean.yml`): `lake env lean AxiomAudit.lean` after `lake build`.
Every main theorem of the three areas must depend on a subset of `propext`, `Classical.choice`,
`Quot.sound` only; `.github/scripts/check_axioms.sh` parses this output and enforces that.
This file is not part of the `Havelock` library.
-/
import Havelock

-- Sinh model (SinhArithmetic.lean)
#print axioms Havelock.Sinh.xi_add_inv
#print axioms Havelock.Sinh.c_int_iff
#print axioms Havelock.Sinh.c_values
#print axioms Havelock.Sinh.isIntegral_of_c_int
#print axioms Havelock.Sinh.xi8_eq
#print axioms Havelock.Sinh.xi8_root
#print axioms Havelock.Sinh.xi10_eq
#print axioms Havelock.Sinh.not_isIntegral_one_seventh

-- Conjecture 13 (Conj13Algebra.lean)
#print axioms Havelock.Conj13.bn_spec
#print axioms Havelock.Conj13.bn_unique
#print axioms Havelock.Conj13.bn_even_paper
#print axioms Havelock.Conj13.bn_odd_paper
#print axioms Havelock.Conj13.bn_sign
#print axioms Havelock.Conj13.phi_crit_and_C1
#print axioms Havelock.Conj13.C1_eq_iff
#print axioms Havelock.Conj13.am_nonneg_and_kernel
#print axioms Havelock.Conj13.am_next
#print axioms Havelock.Conj13.bracket_value
#print axioms Havelock.Conj13.c4even_reduction
#print axioms Havelock.Conj13.c4odd_reduction
#print axioms Havelock.Conj13.c4_at_bn
#print axioms Havelock.Conj13.c4even_pos
#print axioms Havelock.Conj13.c4odd_pos
#print axioms Havelock.Conj13.table_values
#print axioms Havelock.Conj13.table_n7
#print axioms Havelock.Conj13.table_n8
#print axioms Havelock.Conj13.fejer
#print axioms Havelock.Conj13.sum_cos_sq_mul_cos

-- Tanh model (TanhAlgebra.lean and Havelock/Tanh/*.lean)
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
#print axioms Havelock.Tanh.W_self
#print axioms Havelock.Tanh.Vpp_eq
#print axioms Havelock.Tanh.aBlock_symm
#print axioms Havelock.Tanh.bBlock_symm
#print axioms Havelock.Tanh.gBlock_symm
#print axioms Havelock.Tanh.gBlock_half
#print axioms Havelock.Tanh.f8_C_poly
#print axioms Havelock.Tanh.P8_monic
#print axioms Havelock.Tanh.P8_natDegree
#print axioms Havelock.Tanh.aeval_P8
#print axioms Havelock.Tanh.f8_data
#print axioms Havelock.Tanh.f9_data
#print axioms Havelock.Tanh.f10_data
#print axioms Havelock.Tanh.f11_data
#print axioms Havelock.Tanh.f12_data
#print axioms Havelock.Tanh.f13_data
#print axioms Havelock.Tanh.f14_data
#print axioms Havelock.Tanh.f15_data
#print axioms Havelock.Tanh.f16_data

-- Tanh model, large N (TanhLargeN.lean; note Thms 11–18, Prop 19, Cor 20)
#print axioms Havelock.Tanh.LargeN.cR_xiOf
#print axioms Havelock.Tanh.LargeN.aBlock_xiOf
#print axioms Havelock.Tanh.LargeN.bBlock_xiOf
#print axioms Havelock.Tanh.LargeN.gBlock_xiOf
#print axioms Havelock.Tanh.LargeN.beta_pos
#print axioms Havelock.Tanh.LargeN.bBlock_pos
#print axioms Havelock.Tanh.LargeN.K_pos
#print axioms Havelock.Tanh.LargeN.beta_mul_G
#print axioms Havelock.Tanh.LargeN.eight_mul_Phi
#print axioms Havelock.Tanh.LargeN.Phi_strictAntiOn
#print axioms Havelock.Tanh.LargeN.tendsto_Phi_zero
#print axioms Havelock.Tanh.LargeN.tendsto_Phi_atTop
#print axioms Havelock.Tanh.LargeN.existsUnique_Phi_zero
#print axioms Havelock.Tanh.LargeN.tauInf_bracket
#print axioms Havelock.Tanh.LargeN.psi_strictAntiOn
#print axioms Havelock.Tanh.LargeN.aBlock_half_le
#print axioms Havelock.Tanh.LargeN.aBlock_half_lt
#print axioms Havelock.Tanh.LargeN.eight_mul_aBlock_half
#print axioms Havelock.Tanh.LargeN.evenRHS_strictAntiOn
#print axioms Havelock.Tanh.LargeN.tendsto_evenRHS_zero
#print axioms Havelock.Tanh.LargeN.tendsto_evenRHS_atTop
#print axioms Havelock.Tanh.LargeN.aBlock_half_pos_of_le_six
#print axioms Havelock.Tanh.LargeN.aBlock_half_threshold
#print axioms Havelock.Tanh.LargeN.specStable_iff
#print axioms Havelock.Tanh.LargeN.det_xiOf
#print axioms Havelock.Tanh.LargeN.emDefinite_iff_even
#print axioms Havelock.Tanh.LargeN.cR_xiStar_rate
#print axioms Havelock.Tanh.LargeN.tendsto_cR_xiStar
#print axioms Havelock.Tanh.LargeN.cR_xiEM_rate
#print axioms Havelock.Tanh.LargeN.tendsto_cR_xiEM
#print axioms Havelock.Tanh.LargeN.krein_window_odd
#print axioms Havelock.Tanh.LargeN.xiEM_sub_xiStar_le
#print axioms Havelock.Tanh.LargeN.aBlock_half_strictMonoOn
#print axioms Havelock.Tanh.LargeN.xiEM_eq_xiStar_of_even
#print axioms Havelock.Tanh.LargeN.specStable_all_of_even_le_six
#print axioms Havelock.Tanh.LargeN.spacing_rate
#print axioms Havelock.Tanh.LargeN.tendsto_spacing
#print axioms Havelock.Tanh.LargeN.hasSum_rowP
#print axioms Havelock.Tanh.LargeN.hasSum_inv_sq_add_sq'
#print axioms Havelock.Tanh.LargeN.hasSum_note_sum4
#print axioms Havelock.Tanh.LargeN.row_bloch_symbol
#print axioms Havelock.Tanh.LargeN.rowM_real_parts
#print axioms Havelock.Tanh.LargeN.row_stable_iff
#print axioms Havelock.Tanh.LargeN.row_definite_iff
#print axioms Havelock.Tanh.LargeN.rowMyy_half_crit
#print axioms Havelock.Tanh.LargeN.pair_hessian_row
#print axioms Havelock.Tanh.LargeN.row_force_balance_y
#print axioms Havelock.Tanh.LargeN.row_bloch_from_hessian
#print axioms Havelock.Tanh.LargeN.flowBlock_hasEigenvalue_iff
#print axioms Havelock.Tanh.LargeN.flowBlock_re_eigenvalues_iff
#print axioms Havelock.Tanh.LargeN.hessBlock_posDef_iff
#print axioms Havelock.Tanh.LargeN.specStable_iff_spectral
#print axioms Havelock.Tanh.LargeN.emDefinite_iff_posDef
#print axioms Havelock.Tanh.LargeN.row_spectral_iff
#print axioms Havelock.Tanh.LargeN.row_posDef_iff
#print axioms Havelock.Tanh.LargeN.aCrit_strictAntiOn
#print axioms Havelock.Tanh.LargeN.aBlock_half_threshold_of_ge_17
#print axioms Havelock.Tanh.LargeN.xiStar_unique_of_ge_17
#print axioms Havelock.Tanh.LargeN.unique_threshold_of_ge_eight
#print axioms Havelock.Tanh.LargeN.stable_all_radii_of_le_seven

-- Tanh model: irreducibility of f_N (TanhIrred.lean)
#print axioms Havelock.Tanh.Irred.irreducible_of_irreducible_map_zmod
#print axioms Havelock.Tanh.Irred.dvd_X_pow_sub_X
#print axioms Havelock.Tanh.Irred.irreducible_of_checkAll
#print axioms Havelock.Tanh.Irred.f8_irreducible
#print axioms Havelock.Tanh.Irred.f9_irreducible
#print axioms Havelock.Tanh.Irred.f10_irreducible
#print axioms Havelock.Tanh.Irred.f11_irreducible
#print axioms Havelock.Tanh.Irred.f12_irreducible
#print axioms Havelock.Tanh.Irred.f13_irreducible
#print axioms Havelock.Tanh.Irred.f14_irreducible
#print axioms Havelock.Tanh.Irred.f15_irreducible
#print axioms Havelock.Tanh.Irred.f16_irreducible

-- Tanh model: derivation of the blocks (TanhBlocks.lean; note Lemma 1, Thm 2)
#print axioms Havelock.Tanh.fejerC
#print axioms Havelock.Tanh.ring_critical
#print axioms Havelock.Tanh.ring_critical_iff
#print axioms Havelock.Tanh.grad_formula
#print axioms Havelock.Tanh.hessQ_eq
#print axioms Havelock.Tanh.hessSS_eq
#print axioms Havelock.Tanh.hessTT_eq
#print axioms Havelock.Tanh.hessST_eq
#print axioms Havelock.Tanh.hessTS_eq
#print axioms Havelock.Tanh.blocks_eigen
#print axioms Havelock.Tanh.blocks_eigen_zero
#print axioms Havelock.Tanh.blocks_hermitian

-- Sinh model: mode decomposition and threshold (SinhModes.lean)
#print axioms Havelock.Sinh.sinh_ring_critical_iff
#print axioms Havelock.Sinh.sinhD_omega
#print axioms Havelock.Sinh.sinh_blocks_eigen
#print axioms Havelock.Sinh.sinh_blocks_hermitian
#print axioms Havelock.Sinh.threshold
#print axioms Havelock.Sinh.threshold_strict
#print axioms Havelock.Sinh.threshold_cosh
#print axioms Havelock.Sinh.stable_small
#print axioms Havelock.Sinh.threshold_xiStar
#print axioms Havelock.Sinh.xiStar_add_inv

-- Conjecture 13: nonlinear-stability core, abstract (Conj13Stability.lean)
#print axioms Havelock.Conj13.Lyapunov.lyapunov_on_section
#print axioms Havelock.Conj13.Lyapunov.stable_modulo_orbit
#print axioms Havelock.Conj13.Lyapunov.stable_modulo_orbit_global
#print axioms Havelock.Conj13.Lyapunov.lemma_B1
#print axioms Havelock.Conj13.Lyapunov.lemma_B2
#print axioms Havelock.Conj13.Lyapunov.lemma_B3
#print axioms Havelock.Conj13.Lyapunov.quartic_form_coercive

-- Tanh model: the full linearised flow (LinearFlow.lean; note Prop 3)
#print axioms Havelock.Tanh.LargeN.linFlow_mode
#print axioms Havelock.Tanh.LargeN.linFlow_hasEigenvalue_iff
#print axioms Havelock.Tanh.LargeN.linFlow_spectral_iff
#print axioms Havelock.Tanh.LargeN.hessMat_parseval
#print axioms Havelock.Tanh.LargeN.hessMat_posDef_compl_iff
#print axioms Havelock.Tanh.LargeN.hessQ_posDef_compl_iff

-- Conjecture 13: series inputs derived from the Hamiltonian (Conj13Series.lean)
#print axioms Havelock.Conj13.Config.Ham_stag_sub
#print axioms Havelock.Conj13.EvenQuartic.Ham_curve_jet
#print axioms Havelock.Conj13.EvenQuartic.c4even_limit
#print axioms Havelock.Conj13.EvenQuartic.c4even_limit_at_bn
#print axioms Havelock.Conj13.OddSums.sum_inv_one_sub_cos
#print axioms Havelock.Conj13.Odd.Ham_oddLine_jet
#print axioms Havelock.Conj13.Odd.J1_oddAnsatz_jet
#print axioms Havelock.Conj13.Odd.epsOdd_hasDerivAt
#print axioms Havelock.Conj13.Odd.oddNorm_pow4
#print axioms Havelock.Conj13.Taylor.even_iteratedDeriv
#print axioms Havelock.Conj13.Taylor.odd_iteratedDeriv
#print axioms Havelock.Conj13.Taylor.even_iteratedDeriv_at_bn
#print axioms Havelock.Conj13.Taylor.odd_iteratedDeriv_at_bn
#print axioms Havelock.Conj13.Odd4.fejer1
#print axioms Havelock.Conj13.Odd4.fejer2
#print axioms Havelock.Conj13.Odd4.jet_pair
#print axioms Havelock.Conj13.Odd4.jet_logPw
#print axioms Havelock.Conj13.Odd4.Ham_oddAnsatz_jet
#print axioms Havelock.Conj13.Odd4.Ham_oddAnsatz_jet4
