"""
    SquarePXPDynamics

Tools for square-lattice PXP dynamics with dense local gates and a minimal
ITensors-backed PEPS/iPEPS prototype stack: finite product states, periodic
Gamma-lambda iPEPS states, link weights, QR-reduced square-star updates,
Trotter evolution, simple/local diagnostics, and ScarFinder-lite orchestration.
"""
module SquarePXPDynamics

import EDKit

include("Internals.jl")
include("Lattice.jl")          # SquareGeometry + SquareUnitCells
include("PXPModel.jl")         # SpinOps + SquarePXP + StarModels
include("SquareIPEPS.jl")
include("Observables.jl")      # Observables + FiniteIPEPSObservables
include("PEPSKitMeasurements.jl")
include("CTMTrust.jl")
include("StarSimpleUpdate.jl")
include("IPEPSEvolution.jl")
include("FinitePXPEEDBenchmark.jl")
include("PXPValidation.jl")
include("ScarFinderSupport.jl") # CandidateSnapshots + IPEPSCompression
include("ScarFinder.jl")
include("ScarFinderAudit.jl")

using .SpinOps:
    pauli_x,
    pauli_y,
    pauli_z,
    identity2,
    projector_up,
    projector_down,
    kron_all,
    embed_one_site
using .SquareGeometry:
    SquareCoord,
    SquareUnitCell,
    OneSiteSquareUC,
    FiveSiteSquareUC,
    square_neighbor,
    square_star_sites,
    square_star_color,
    disjoint_square_stars,
    unit_cell_representatives,
    wrap_square_coord
using .SquarePXP:
    SQUARE_STAR_SITES,
    square_pxp_star_hamiltonian,
    square_star_blockade_projector,
    square_pxp_gate,
    projected_square_pxp_gate,
    square_star_basis_allowed
using .StarModels:
    AbstractStarModel,
    PXPStarModel,
    AbstractModelProtocol,
    StaticModel,
    model_at,
    star_site_order,
    star_hamiltonian,
    star_gate,
    star_gate_itensor
using .SquareUnitCells:
    PeriodicSquareUnitCell,
    wrap,
    neighbor,
    update_centers,
    assert_five_color_compatible,
    stars_are_disjoint_mod_unitcell,
    BondKey,
    bondkey
using .SquareIPEPS: SquareIPEPSState, product_square_ipeps, checkerboard_square_ipeps
using .SquareIPEPS: physical_index, link_index
using .SquareIPEPS: unitcell_reps, physical_dim, simple_weight_dim, copy_state
using .SquareIPEPS: link_weight, set_link_weight!, link_weight_tensor
using .SquareIPEPS: state_version, log_norm
using .SquareIPEPS: absorb_link_weight, deabsorb_link_weight
using .SquareIPEPS: weight_entropy, bond_entropy, all_bond_entropies
using .SquareIPEPS: normalize_link_weights!
using .SquareIPEPS: square_pxp_gate_itensor, projected_square_pxp_gate_itensor
using .Observables: local_density_simple, density_simple, sublattice_densities
using .Observables: sublattice_imbalance_simple, checkerboard_structure_factor_simple
using .Observables: nearest_neighbor_density_simple, blockade_violation_simple
using .Observables: star_expectation_simple, pxp_energy_density_simple
using .Observables: mean_bond_entropy, max_bond_entropy
using .Observables: SimpleObservableSummary, measure_simple
using .Observables: local_x_simple, local_y_simple, local_z_simple
using .Observables: nearest_neighbor_zz_simple
using .FiniteIPEPSObservables:
    dense_state_finite,
    exact_one_site_expectation_finite,
    exact_nearest_neighbor_expectation_finite,
    exact_star_expectation_finite,
    exact_density_finite,
    exact_all_down_return_probability_finite,
    exact_blockade_violation_finite,
    exact_pxp_energy_density_finite
using .PEPSKitMeasurements: PEPSKitCTMRGParams, default_ctmrg_params
using .PEPSKitMeasurements: PEPSKitMeasurementContext, CTMRGDiagnostics
using .PEPSKitMeasurements: CTMObservableSummary, CTMValidationPoint
using .PEPSKitMeasurements: configure_ctm_threading!, configure_ctm_threading_from_env!
using .PEPSKitMeasurements: to_pepskit_infinitepeps
using .PEPSKitMeasurements: pepskit_ctmrg_context, local_density_ctm
using .PEPSKitMeasurements: nearest_neighbor_density_ctm, blockade_violation_ctm
using .PEPSKitMeasurements: star_expectation_ctm, pxp_energy_density_ctm, measure_ctm
using .PEPSKitMeasurements: correlator_ctm, correlation_length_ctm
using .PEPSKitMeasurements: ctm_diagnostics, validate_ctm_sweep, write_ctm_validation_csv
using .PEPSKitMeasurements: ctm_chi_sensitivity_sweep
using .PEPSKitMeasurements: assert_fresh_pepskit_context
using .CTMTrust: CTMTrustPolicy, CTMTrustAssessment, assess_ctm_trust, write_ctm_trust_csv
using .CTMTrust: tight_ctm_trust_policy, calibrated_ctm_trust_policy
using .StarSimpleUpdate: StarUpdateInfo, project_star!
using .IPEPSEvolution: TrotterParams, EvolutionLog, trotter_sequence, evolve!, reverse_evolve!
using .FinitePXPEEDBenchmark:
    PXPSquareSpaceGroupBasis,
    PXPEEDBenchmarkConfig,
    PXPEEDSample,
    PXPEEDBenchmarkResult,
    pxp_ed_space_group_basis,
    pxp_ed_constrained_count,
    pxp_ed_group_order,
    pxp_ed_boundary_condition,
    pxp_ed_symmetry_sector,
    pxp_ed_observable_scope,
    pxp_ed_reference_label,
    pxp_ed_site_density_operator,
    pxp_ed_region_density_operator,
    pxp_ed_initial_state,
    pxp_ed_hamiltonian_operator,
    sparse_pxp_ed_hamiltonian,
    run_pxp_ed_benchmark,
    write_pxp_ed_benchmark_json
using .PXPValidation:
    TrustedCTMMeasurement,
    measure_ctm_trusted,
    PXPValidationConfig,
    PXPValidationMetadata,
    PXPIPEPSSample,
    PXPEDComparisonSample,
    PXPValidationReport,
    validate_pxp_ed_ipeps,
    write_pxp_validation_json,
    PXPConvergenceConfig,
    PXPConvergenceReport,
    validate_pxp_convergence,
    write_pxp_convergence_json,
    PXPReversibilityReport,
    validate_pxp_reversibility,
    PXPLargerDBenchmarkConfig,
    PXPLargerDBenchmarkSummary,
    PXPLargerDBenchmarkRun,
    PXPLargerDBenchmarkReport,
    run_pxp_larger_d_benchmark,
    write_pxp_larger_d_benchmark_json,
    write_pxp_larger_d_benchmark_csv
using .CandidateSnapshots:
    SQUARE_IPEPS_SNAPSHOT_FORMAT_VERSION,
    write_square_ipeps_snapshot,
    load_square_ipeps_snapshot
using .IPEPSCompression:
    IPEPSCompressionInfo,
    compress_to_target_maxdim!
using .ScarFinder:
    ScarFinderParams,
    ScarFinderCandidateScore,
    ScarFinderIteration,
    ScarFinderResult,
    ScarFinderCompressionInfo,
    MeasurementBackend,
    SimpleBackend,
    TrustedCTMBackend,
    measure_scarfinder,
    CandidateStore,
    NoCandidateStore,
    JSONCandidateStore,
    JLD2CandidateStore,
    ScarFinderObjective,
    RevivalObjective,
    TargetEnergyObjective,
    LowVarianceObjective,
    CompositeObjective,
    rank_scarfinder_candidates,
    write_scarfinder_log,
    scarfinder!
using .ScarFinderAudit:
    ScarFinderAuditConfig,
    ScarFinderAuditRow,
    ScarFinderAuditStability,
    ScarFinderAuditReport,
    run_scarfinder_audit,
    write_scarfinder_audit_json,
    write_scarfinder_audit_csv

export pauli_x, pauli_y, pauli_z, identity2, projector_up, projector_down
export kron_all, embed_one_site
export SquareCoord, SquareUnitCell, OneSiteSquareUC, FiveSiteSquareUC
export square_neighbor, square_star_sites, square_star_color, disjoint_square_stars
export unit_cell_representatives, wrap_square_coord
export SQUARE_STAR_SITES, square_pxp_star_hamiltonian, square_star_blockade_projector
export square_pxp_gate, projected_square_pxp_gate, square_star_basis_allowed
export AbstractStarModel, PXPStarModel
export AbstractModelProtocol, StaticModel, model_at
export star_site_order
export star_hamiltonian, star_gate, star_gate_itensor
export physical_index, link_index
export PeriodicSquareUnitCell
export wrap, neighbor, update_centers, assert_five_color_compatible
export stars_are_disjoint_mod_unitcell
export BondKey, bondkey
export SquareIPEPSState, product_square_ipeps, checkerboard_square_ipeps
export unitcell_reps, physical_dim, simple_weight_dim, copy_state
export link_weight, set_link_weight!, link_weight_tensor
export state_version, log_norm
export absorb_link_weight, deabsorb_link_weight
export weight_entropy, bond_entropy, all_bond_entropies, normalize_link_weights!
export square_pxp_gate_itensor, projected_square_pxp_gate_itensor
export local_density_simple, density_simple, sublattice_densities
export sublattice_imbalance_simple, checkerboard_structure_factor_simple
export nearest_neighbor_density_simple, blockade_violation_simple
export star_expectation_simple, pxp_energy_density_simple
export mean_bond_entropy, max_bond_entropy
export SimpleObservableSummary, measure_simple
export local_x_simple, local_y_simple, local_z_simple
export nearest_neighbor_zz_simple
export dense_state_finite
export exact_one_site_expectation_finite, exact_nearest_neighbor_expectation_finite
export exact_star_expectation_finite, exact_density_finite
export exact_all_down_return_probability_finite
export exact_blockade_violation_finite, exact_pxp_energy_density_finite
export PEPSKitCTMRGParams, default_ctmrg_params
export PEPSKitMeasurementContext, CTMRGDiagnostics, CTMObservableSummary
export CTMValidationPoint
export configure_ctm_threading!, configure_ctm_threading_from_env!
export to_pepskit_infinitepeps, pepskit_ctmrg_context
export local_density_ctm, nearest_neighbor_density_ctm
export blockade_violation_ctm, star_expectation_ctm, pxp_energy_density_ctm, measure_ctm
export correlator_ctm, correlation_length_ctm
export ctm_diagnostics, validate_ctm_sweep, write_ctm_validation_csv
export ctm_chi_sensitivity_sweep
export assert_fresh_pepskit_context
export CTMTrustPolicy, CTMTrustAssessment, assess_ctm_trust, write_ctm_trust_csv
export tight_ctm_trust_policy, calibrated_ctm_trust_policy
export StarUpdateInfo, project_star!
export TrotterParams, EvolutionLog, trotter_sequence, evolve!, reverse_evolve!
export PXPSquareSpaceGroupBasis
export PXPEEDBenchmarkConfig, PXPEEDSample, PXPEEDBenchmarkResult
export pxp_ed_space_group_basis, pxp_ed_constrained_count, pxp_ed_group_order
export pxp_ed_boundary_condition, pxp_ed_symmetry_sector, pxp_ed_observable_scope
export pxp_ed_reference_label, pxp_ed_site_density_operator, pxp_ed_region_density_operator
export pxp_ed_initial_state, pxp_ed_hamiltonian_operator, sparse_pxp_ed_hamiltonian
export run_pxp_ed_benchmark, write_pxp_ed_benchmark_json
export TrustedCTMMeasurement, measure_ctm_trusted
export PXPValidationConfig, PXPValidationMetadata, PXPIPEPSSample
export PXPEDComparisonSample, PXPValidationReport, validate_pxp_ed_ipeps
export write_pxp_validation_json
export PXPConvergenceConfig, PXPConvergenceReport, validate_pxp_convergence
export write_pxp_convergence_json
export PXPReversibilityReport, validate_pxp_reversibility
export PXPLargerDBenchmarkConfig, PXPLargerDBenchmarkSummary
export PXPLargerDBenchmarkRun, PXPLargerDBenchmarkReport
export run_pxp_larger_d_benchmark
export write_pxp_larger_d_benchmark_json, write_pxp_larger_d_benchmark_csv
export ScarFinderParams, ScarFinderCandidateScore, ScarFinderIteration, ScarFinderResult
export ScarFinderCompressionInfo
export MeasurementBackend, SimpleBackend, TrustedCTMBackend, measure_scarfinder
export CandidateStore, NoCandidateStore, JSONCandidateStore, JLD2CandidateStore
export SQUARE_IPEPS_SNAPSHOT_FORMAT_VERSION
export write_square_ipeps_snapshot, load_square_ipeps_snapshot
export IPEPSCompressionInfo, compress_to_target_maxdim!
export ScarFinderObjective, RevivalObjective, TargetEnergyObjective
export LowVarianceObjective, CompositeObjective
export rank_scarfinder_candidates, write_scarfinder_log, scarfinder!
export ScarFinderAuditConfig, ScarFinderAuditRow, ScarFinderAuditStability
export ScarFinderAuditReport
export run_scarfinder_audit, write_scarfinder_audit_json, write_scarfinder_audit_csv

end
