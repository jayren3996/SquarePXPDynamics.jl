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

using .SpinOps: pauli_x, pauli_y, pauli_z, identity2, projector_up, projector_down,
                kron_all, embed_one_site
using .SquareGeometry: SquareCoord, SquareUnitCell, OneSiteSquareUC, FiveSiteSquareUC,
                       square_neighbor, square_star_sites, square_star_color,
                       disjoint_square_stars, unit_cell_representatives, wrap_square_coord
using .SquarePXP: SQUARE_STAR_SITES, square_pxp_star_hamiltonian,
                  square_star_blockade_projector, square_pxp_gate,
                  projected_square_pxp_gate, square_star_basis_allowed
using .SquarePEPS: SquarePEPSState, product_square_peps, site_tensor, physical_index, link_index

export pauli_x, pauli_y, pauli_z, identity2, projector_up, projector_down
export kron_all, embed_one_site
export SquareCoord, SquareUnitCell, OneSiteSquareUC, FiveSiteSquareUC
export square_neighbor, square_star_sites, square_star_color, disjoint_square_stars
export unit_cell_representatives, wrap_square_coord
export SQUARE_STAR_SITES, square_pxp_star_hamiltonian, square_star_blockade_projector
export square_pxp_gate, projected_square_pxp_gate, square_star_basis_allowed
export SquarePEPSState, product_square_peps, site_tensor, physical_index, link_index

end
