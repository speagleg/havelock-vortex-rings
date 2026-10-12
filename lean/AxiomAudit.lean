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
#print axioms Havelock.Conj13.P1_n7
#print axioms Havelock.Conj13.P1_n8
#print axioms Havelock.Conj13.P1_even_rows
#print axioms Havelock.Conj13.P2_complete_square
#print axioms Havelock.Conj13.P2_n7
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
