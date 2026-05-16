"""
    SquarePXPDynamics

Tools for square-lattice PXP dynamics with dense local gates and a minimal
ITensors-backed PEPS/iPEPS prototype stack: finite product states, periodic
Gamma-lambda iPEPS states, link weights, QR-reduced square-star updates,
Trotter evolution, simple/local diagnostics, and ScarFinder-lite orchestration.
"""
module SquarePXPDynamics

import EDKit

include("SpinOps.jl")
include("SquareGeometry.jl")
include("SquarePXP.jl")
include("SquarePEPS.jl")
include("SquareUnitCells.jl")
include("SquareIPEPS.jl")
include("GaugeDiagnostics.jl")
include("StarModels.jl")
include("Observables.jl")
include("PEPSKitMeasurements.jl")
include("CTMTrust.jl")
include("StarSimpleUpdate.jl")
include("IPEPSEvolution.jl")
include("Benchmarks.jl")
include("FiniteTFIMReference.jl")
include("FiniteMPSTFIMReference.jl")
include("FinitePXPEEDBenchmark.jl")
include("ScarFinder.jl")

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
    TFIMStarModel,
    AbstractModelProtocol,
    StaticModel,
    model_at,
    star_site_order,
    tfim_pauli_convention,
    star_hamiltonian,
    star_gate,
    star_gate_itensor,
    tfim_product_basis_energy
using .SquarePEPS:
    SquarePEPSState, product_square_peps, site_tensor, physical_index, link_index
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
using .SquareIPEPS: unitcell_reps, physical_dim, simple_weight_dim, copy_state
using .SquareIPEPS: link_weight, set_link_weight!, link_weight_tensor
using .SquareIPEPS: state_version, log_norm
using .SquareIPEPS: absorb_link_weight, deabsorb_link_weight
using .SquareIPEPS: weight_entropy, bond_entropy, all_bond_entropies
using .SquareIPEPS: square_pxp_gate_itensor, projected_square_pxp_gate_itensor
using .GaugeDiagnostics: SimpleGaugeDiagnostic
using .GaugeDiagnostics: gauge_diagnostic_simple, gauge_deviation_simple
using .GaugeDiagnostics: all_gauge_deviations_simple
using .Observables: local_density_simple, density_simple, sublattice_densities
using .Observables: nearest_neighbor_density_simple, blockade_violation_simple
using .Observables: star_expectation_simple, pxp_energy_density_simple
using .Observables: mean_bond_entropy, max_bond_entropy
using .Observables: SimpleObservableSummary, measure_simple
using .Observables: local_x_simple, local_y_simple, local_z_simple
using .Observables: nearest_neighbor_zz_simple
using .Observables: tfim_energy_density_star_simple, tfim_energy_density_decomposed_simple
using .Observables: TFIMObservableSummary, measure_tfim_simple
using .PEPSKitMeasurements: PEPSKitCTMRGParams, PEPSKitMeasurementContext, CTMRGDiagnostics
using .PEPSKitMeasurements: CTMObservableSummary, CTMValidationPoint
using .PEPSKitMeasurements: to_pepskit_infinitepeps
using .PEPSKitMeasurements: pepskit_ctmrg_context, local_density_ctm
using .PEPSKitMeasurements: nearest_neighbor_density_ctm, blockade_violation_ctm
using .PEPSKitMeasurements: star_expectation_ctm, pxp_energy_density_ctm, measure_ctm
using .PEPSKitMeasurements: ctm_diagnostics, validate_ctm_sweep, write_ctm_validation_csv
using .CTMTrust: CTMTrustPolicy, CTMTrustAssessment, assess_ctm_trust, write_ctm_trust_csv
using .StarSimpleUpdate: StarUpdateInfo, project_star!
using .IPEPSEvolution: TrotterParams, EvolutionLog, trotter_sequence, evolve!
using .Benchmarks:
    BenchmarkSpec,
    BenchmarkMetadata,
    EvolutionDiagnostics,
    BenchmarkSample,
    BenchmarkResult,
    run_benchmark,
    write_benchmark_json,
    write_benchmark_csv
using .FiniteTFIMReference:
    FiniteTFIMReferenceSample,
    finite_tfim_hamiltonian,
    finite_tfim_product_state,
    measure_finite_tfim,
    run_finite_tfim_reference
using .FiniteMPSTFIMReference:
    FiniteMPSTFIMMetadata,
    FiniteMPSTFIMSample,
    FiniteMPSTFIMResult,
    finite_mps_site_index,
    finite_mps_square_lattice_bonds,
    run_finite_mps_tfim_reference
using .FinitePXPEEDBenchmark:
    PXPSquareSpaceGroupBasis,
    PXPEEDBenchmarkConfig,
    PXPEEDSample,
    PXPEEDBenchmarkResult,
    pxp_ed_space_group_basis,
    pxp_ed_constrained_count,
    pxp_ed_group_order,
    pxp_ed_initial_state,
    pxp_ed_hamiltonian_operator,
    sparse_pxp_ed_hamiltonian,
    run_pxp_ed_benchmark,
    write_pxp_ed_benchmark_json
using .ScarFinder:
    ScarFinderParams,
    ScarFinderCandidateScore,
    ScarFinderIteration,
    ScarFinderResult,
    rank_scarfinder_candidates,
    write_scarfinder_log,
    scarfinder!

export pauli_x, pauli_y, pauli_z, identity2, projector_up, projector_down
export kron_all, embed_one_site
export SquareCoord, SquareUnitCell, OneSiteSquareUC, FiveSiteSquareUC
export square_neighbor, square_star_sites, square_star_color, disjoint_square_stars
export unit_cell_representatives, wrap_square_coord
export SQUARE_STAR_SITES, square_pxp_star_hamiltonian, square_star_blockade_projector
export square_pxp_gate, projected_square_pxp_gate, square_star_basis_allowed
export AbstractStarModel, PXPStarModel, TFIMStarModel
export AbstractModelProtocol, StaticModel, model_at
export star_site_order, tfim_pauli_convention
export star_hamiltonian, star_gate, star_gate_itensor, tfim_product_basis_energy
export SquarePEPSState, product_square_peps, site_tensor, physical_index, link_index
export PeriodicSquareUnitCell
export wrap, neighbor, update_centers, assert_five_color_compatible
export stars_are_disjoint_mod_unitcell
export BondKey, bondkey
export SquareIPEPSState, product_square_ipeps, checkerboard_square_ipeps
export unitcell_reps, physical_dim, simple_weight_dim, copy_state
export link_weight, set_link_weight!, link_weight_tensor
export state_version, log_norm
export absorb_link_weight, deabsorb_link_weight
export weight_entropy, bond_entropy, all_bond_entropies
export square_pxp_gate_itensor, projected_square_pxp_gate_itensor
export SimpleGaugeDiagnostic
export gauge_diagnostic_simple, gauge_deviation_simple, all_gauge_deviations_simple
export local_density_simple, density_simple, sublattice_densities
export nearest_neighbor_density_simple, blockade_violation_simple
export star_expectation_simple, pxp_energy_density_simple
export mean_bond_entropy, max_bond_entropy
export SimpleObservableSummary, measure_simple
export local_x_simple, local_y_simple, local_z_simple
export nearest_neighbor_zz_simple
export tfim_energy_density_star_simple, tfim_energy_density_decomposed_simple
export TFIMObservableSummary, measure_tfim_simple
export PEPSKitCTMRGParams, PEPSKitMeasurementContext, CTMRGDiagnostics, CTMObservableSummary
export CTMValidationPoint
export to_pepskit_infinitepeps, pepskit_ctmrg_context
export local_density_ctm, nearest_neighbor_density_ctm
export blockade_violation_ctm, star_expectation_ctm, pxp_energy_density_ctm, measure_ctm
export ctm_diagnostics, validate_ctm_sweep, write_ctm_validation_csv
export CTMTrustPolicy, CTMTrustAssessment, assess_ctm_trust, write_ctm_trust_csv
export StarUpdateInfo, project_star!
export TrotterParams, EvolutionLog, trotter_sequence, evolve!
export BenchmarkSpec, BenchmarkMetadata, EvolutionDiagnostics, BenchmarkSample, BenchmarkResult
export run_benchmark, write_benchmark_json, write_benchmark_csv
export FiniteTFIMReferenceSample
export finite_tfim_hamiltonian, finite_tfim_product_state
export measure_finite_tfim, run_finite_tfim_reference
export FiniteMPSTFIMMetadata, FiniteMPSTFIMSample, FiniteMPSTFIMResult
export finite_mps_site_index, finite_mps_square_lattice_bonds
export run_finite_mps_tfim_reference
export PXPSquareSpaceGroupBasis
export PXPEEDBenchmarkConfig, PXPEEDSample, PXPEEDBenchmarkResult
export pxp_ed_space_group_basis, pxp_ed_constrained_count, pxp_ed_group_order
export pxp_ed_initial_state, pxp_ed_hamiltonian_operator, sparse_pxp_ed_hamiltonian
export run_pxp_ed_benchmark, write_pxp_ed_benchmark_json
export ScarFinderParams, ScarFinderCandidateScore, ScarFinderIteration, ScarFinderResult
export rank_scarfinder_candidates, write_scarfinder_log, scarfinder!

end
