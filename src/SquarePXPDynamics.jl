"""
    SquarePXPDynamics

Tools for square-lattice PXP dynamics with dense local gates and a minimal
ITensors-backed finite square PEPS state container.
"""
module SquarePXPDynamics

include("SpinOps.jl")
include("SquareGeometry.jl")
include("SquarePXP.jl")
include("SquarePEPS.jl")
include("SquareUnitCells.jl")
include("SquareIPEPS.jl")
include("Observables.jl")
include("StarSimpleUpdate.jl")
include("IPEPSEvolution.jl")
include("ScarFinder.jl")

using .SpinOps: pauli_x, pauli_y, pauli_z, identity2, projector_up, projector_down,
                kron_all, embed_one_site
using .SquareGeometry: SquareCoord, SquareUnitCell, OneSiteSquareUC, FiveSiteSquareUC,
                       square_neighbor, square_star_sites, square_star_color,
                       disjoint_square_stars, unit_cell_representatives, wrap_square_coord
using .SquarePXP: SQUARE_STAR_SITES, square_pxp_star_hamiltonian,
                  square_star_blockade_projector, square_pxp_gate,
                  projected_square_pxp_gate, square_star_basis_allowed
using .SquarePEPS: SquarePEPSState, product_square_peps, site_tensor, physical_index, link_index
using .SquareUnitCells: PeriodicSquareUnitCell, wrap, neighbor, update_centers,
                        assert_five_color_compatible, stars_are_disjoint_mod_unitcell,
                        BondKey, bondkey
using .SquareIPEPS: SquareIPEPSState, product_square_ipeps, checkerboard_square_ipeps
using .SquareIPEPS: link_weight, set_link_weight!, link_weight_tensor
using .SquareIPEPS: absorb_link_weight, deabsorb_link_weight
using .SquareIPEPS: weight_entropy, bond_entropy, all_bond_entropies
using .SquareIPEPS: square_pxp_gate_itensor, projected_square_pxp_gate_itensor
using .Observables: local_density_simple, density_simple, sublattice_densities
using .Observables: nearest_neighbor_density_simple, blockade_violation_simple
using .Observables: star_expectation_simple, pxp_energy_density_simple
using .Observables: mean_bond_entropy, max_bond_entropy
using .Observables: SimpleObservableSummary, measure_simple
using .StarSimpleUpdate: StarUpdateInfo, project_star!
using .IPEPSEvolution: TrotterParams, EvolutionLog, trotter_sequence, evolve!
using .ScarFinder: ScarFinderParams, ScarFinderIteration, ScarFinderResult, scarfinder!

export pauli_x, pauli_y, pauli_z, identity2, projector_up, projector_down
export kron_all, embed_one_site
export SquareCoord, SquareUnitCell, OneSiteSquareUC, FiveSiteSquareUC
export square_neighbor, square_star_sites, square_star_color, disjoint_square_stars
export unit_cell_representatives, wrap_square_coord
export SQUARE_STAR_SITES, square_pxp_star_hamiltonian, square_star_blockade_projector
export square_pxp_gate, projected_square_pxp_gate, square_star_basis_allowed
export SquarePEPSState, product_square_peps, site_tensor, physical_index, link_index
export PeriodicSquareUnitCell
export wrap, neighbor, update_centers, assert_five_color_compatible
export stars_are_disjoint_mod_unitcell
export BondKey, bondkey
export SquareIPEPSState, product_square_ipeps, checkerboard_square_ipeps
export link_weight, set_link_weight!, link_weight_tensor
export absorb_link_weight, deabsorb_link_weight
export weight_entropy, bond_entropy, all_bond_entropies
export square_pxp_gate_itensor, projected_square_pxp_gate_itensor
export local_density_simple, density_simple, sublattice_densities
export nearest_neighbor_density_simple, blockade_violation_simple
export star_expectation_simple, pxp_energy_density_simple
export mean_bond_entropy, max_bond_entropy
export SimpleObservableSummary, measure_simple
export StarUpdateInfo, project_star!
export TrotterParams, EvolutionLog, trotter_sequence, evolve!
export ScarFinderParams, ScarFinderIteration, ScarFinderResult, scarfinder!

end
