module SquarePXP

using LinearAlgebra
using ..SpinOps

export SQUARE_STAR_SITES
export square_pxp_star_hamiltonian, square_star_blockade_projector
export square_pxp_gate, projected_square_pxp_gate
export square_star_basis_allowed

"""
    SQUARE_STAR_SITES

Number of sites in the dense square-star convention: center, right, up, left, down.
"""
const SQUARE_STAR_SITES = 5

"""
    square_pxp_star_hamiltonian()

Return the dense 5-site square-star PXP Hamiltonian term
`X_center * P_down_right * P_down_up * P_down_left * P_down_down`.
"""
function square_pxp_star_hamiltonian()
    return kron_all([
        pauli_x(),
        projector_down(),
        projector_down(),
        projector_down(),
        projector_down(),
    ])
end

function _basis_bits(index::Integer, nsites::Integer)
    1 <= index <= 2^nsites || throw(ArgumentError("index out of range"))
    value = index - 1
    return ntuple(site -> (value >> (nsites - site)) & 1, nsites)
end

"""
    square_star_basis_allowed(bits)

Return whether a 5-site square-star computational basis tuple satisfies the local
blockade constraint. The basis order is `(center, right, up, left, down)`, with
`0` representing `|up>` and `1` representing `|down>`.
"""
function square_star_basis_allowed(bits)
    length(bits) == SQUARE_STAR_SITES ||
        throw(ArgumentError("square star basis must have 5 sites"))
    all(bit -> bit === 0 || bit === 1, bits) ||
        throw(ArgumentError("square star basis bits must be integers 0 or 1"))
    center_is_up = bits[1] == 0
    for site = 2:SQUARE_STAR_SITES
        center_is_up && bits[site] == 0 && return false
    end
    return true
end

"""
    square_star_blockade_projector()

Return the dense diagonal projector onto locally blockade-allowed square-star
basis states.
"""
function square_star_blockade_projector()
    dim = 2^SQUARE_STAR_SITES
    projector = zeros(ComplexF64, dim, dim)
    for i = 1:dim
        projector[i, i] =
            square_star_basis_allowed(_basis_bits(i, SQUARE_STAR_SITES)) ? 1 : 0
    end
    return projector
end

function _validate_finite_step(step::Real)
    isfinite(step) || throw(ArgumentError("step must be finite"))
    return step
end

"""
    square_pxp_gate(step; evolution = :real)

Return the dense square-star PXP evolution gate. Use `evolution = :real` for
`exp(-im * step * H)` and `evolution = :imaginary` for `exp(-step * H)`.
"""
function square_pxp_gate(step::Real; evolution::Symbol = :real)
    _validate_finite_step(step)
    H = square_pxp_star_hamiltonian()
    if evolution === :real
        return exp(-im * step * H)
    elseif evolution === :imaginary
        return exp(-step * H)
    else
        throw(ArgumentError("evolution must be :real or :imaginary"))
    end
end

"""
    projected_square_pxp_gate(step; evolution = :real)

Return `square_star_blockade_projector() * square_pxp_gate(step; evolution)`,
keeping local blockade enforcement explicit.
"""
function projected_square_pxp_gate(step::Real; evolution::Symbol = :real)
    _validate_finite_step(step)
    return square_star_blockade_projector() * square_pxp_gate(step; evolution)
end

end
