module SpinOps

export pauli_x, pauli_y, pauli_z, identity2, projector_up, projector_down
export kron_all, embed_one_site

"""
    identity2()

Return the 2x2 identity matrix as `ComplexF64`.
"""
function identity2()
    return ComplexF64[1 0; 0 1]
end

"""
    pauli_x()

Return the spin-1/2 Pauli X matrix as `ComplexF64`.
"""
function pauli_x()
    return ComplexF64[0 1; 1 0]
end

"""
    pauli_y()

Return the spin-1/2 Pauli Y matrix as `ComplexF64`.
"""
function pauli_y()
    return ComplexF64[0 -im; im 0]
end

"""
    pauli_z()

Return the spin-1/2 Pauli Z matrix as `ComplexF64`.
"""
function pauli_z()
    return ComplexF64[1 0; 0 -1]
end

"""
    projector_up()

Return the single-site projector onto `|up> = |0>`.
"""
function projector_up()
    return ComplexF64[1 0; 0 0]
end

"""
    projector_down()

Return the single-site projector onto `|down> = |1>`.
"""
function projector_down()
    return ComplexF64[0 0; 0 1]
end

"""
    kron_all(ops)

Return the left-to-right Kronecker product of a nonempty vector of matrices.
All inputs are converted to `ComplexF64` matrices before multiplication.
"""
function kron_all(ops::AbstractVector{<:AbstractMatrix})
    isempty(ops) && throw(ArgumentError("ops must be nonempty"))
    out = Matrix{ComplexF64}(ops[1])
    for op in ops[2:end]
        out = kron(out, Matrix{ComplexF64}(op))
    end
    return out
end

"""
    embed_one_site(op, site, nsites)

Embed a 2x2 single-site operator into an `nsites` spin-1/2 Hilbert space.
Site numbering is one-based and follows the left-to-right tensor-product order.
"""
function embed_one_site(op::AbstractMatrix, site::Integer, nsites::Integer)
    size(op) == (2, 2) || throw(ArgumentError("op must be 2x2"))
    1 <= site <= nsites || throw(ArgumentError("site must be in 1:nsites"))
    return kron_all([i == site ? op : identity2() for i = 1:nsites])
end

end


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
    projected_square_pxp_gate(step; evolution = :real, projection = :left)

Return a dense square-star PXP evolution gate with local blockade projection.
The default `projection = :left` preserves the historical `P * U` behavior and
assumes the input vector is already in the constrained sector. Use
`projection = :sandwich` to return `P * U * P` for explicit constrained-sector
action on raw local vectors. For the current square-star PXP Hamiltonian `H`,
`H` is supported on states with all four neighbors down (so `U = I` outside
that support) and maps constrained basis states to constrained basis states
inside it; hence `[P, U] = 0`, and `:left` and `:sandwich` produce equivalent
operators. The keyword makes the chosen convention explicit at call sites.
"""
function projected_square_pxp_gate(
    step::Real;
    evolution::Symbol = :real,
    projection::Symbol = :left,
)
    _validate_finite_step(step)
    P = square_star_blockade_projector()
    U = square_pxp_gate(step; evolution)
    projection === :left && return P * U
    projection === :sandwich && return P * U * P
    throw(ArgumentError("projection must be :left or :sandwich"))
end

end


module StarModels

using ITensors
import ..SquarePXP:
    SQUARE_STAR_SITES,
    square_pxp_star_hamiltonian,
    square_pxp_gate,
    projected_square_pxp_gate

export AbstractStarModel, PXPStarModel
export AbstractModelProtocol, StaticModel, model_at
export star_site_order
export star_hamiltonian, star_gate, star_gate_itensor

"""
    AbstractStarModel

Abstract supertype for dense five-site square-star models.
"""
abstract type AbstractStarModel end

"""
    PXPStarModel(projected)

Five-site square-star PXP model wrapper. When `projected` is `true`, dense
gates include the local blockade projector.
"""
struct PXPStarModel <: AbstractStarModel
    projected::Bool
end

"""
    AbstractModelProtocol

Abstract supertype for protocols that select a star model during evolution.
"""
abstract type AbstractModelProtocol end

"""
    StaticModel(model)

Model protocol that returns the same star model for every time and step.
"""
struct StaticModel{M<:AbstractStarModel} <: AbstractModelProtocol
    model::M
end

"""
    model_at(protocol, time, step)

Return the star model selected by `protocol` at `time` and integer `step`.
"""
model_at(protocol::StaticModel, time, step) = protocol.model

"""
    star_site_order()

Return the dense square-star site order `(center, right, up, left, down)`.
"""
star_site_order() = (:center, :right, :up, :left, :down)

"""
    star_hamiltonian(model)

Return the dense 32x32 Hamiltonian for a five-site square-star model in
`star_site_order()`.
"""
star_hamiltonian(model::PXPStarModel) = square_pxp_star_hamiltonian()

"""
    star_gate(model, step; evolution = :real)

Return the dense 32x32 square-star evolution gate for `model`. Use
`evolution = :real` for `exp(-im * step * H)` and `evolution = :imaginary` for
`exp(-step * H)`.
"""
function star_gate(model::PXPStarModel, step::Real; evolution::Symbol = :real)
    return model.projected ? projected_square_pxp_gate(step; evolution) :
           square_pxp_gate(step; evolution)
end

function _square_star_indices(site_indices)
    sites = collect(site_indices)
    length(sites) == SQUARE_STAR_SITES || throw(
        ArgumentError(
            "square-star gate requires 5 physical indices in (center, right, up, left, down) order",
        ),
    )
    all(i -> dim(i) == 2, sites) ||
        throw(ArgumentError("square-star physical indices must all have dimension 2"))
    return sites
end

function _square_star_dense_index(values)
    idx = 1
    for (site, value) in enumerate(values)
        1 <= value <= 2 || throw(ArgumentError("square-star basis values must be 1 or 2"))
        idx += (value - 1) * 2^(SQUARE_STAR_SITES - site)
    end
    return idx
end

function _square_star_gate_itensor(dense_gate::AbstractMatrix, site_indices)
    size(dense_gate) == (2^SQUARE_STAR_SITES, 2^SQUARE_STAR_SITES) ||
        throw(ArgumentError("dense square-star gate must be 32x32"))
    sites = _square_star_indices(site_indices)
    out = prime.(sites)
    data = zeros(ComplexF64, ntuple(Returns(2), 2 * SQUARE_STAR_SITES))

    for out_values in Iterators.product((1:2 for _ = 1:SQUARE_STAR_SITES)...)
        out_idx = _square_star_dense_index(out_values)
        for in_values in Iterators.product((1:2 for _ = 1:SQUARE_STAR_SITES)...)
            in_idx = _square_star_dense_index(in_values)
            data[out_values..., in_values...] = dense_gate[out_idx, in_idx]
        end
    end

    return ITensor(data, out..., sites...)
end

"""
    star_gate_itensor(model, site_indices, step; evolution = :real)

Wrap `star_gate(model, step; evolution)` as an ITensor with primed output
indices followed by unprimed input indices. `site_indices` must contain five
dimension-2 indices in dense square-star order `(center, right, up, left, down)`.
"""
function star_gate_itensor(
    model::AbstractStarModel,
    site_indices,
    step::Real;
    evolution::Symbol = :real,
)
    return _square_star_gate_itensor(star_gate(model, step; evolution), site_indices)
end

function star_gate_itensor(
    model::AbstractStarModel,
    step::Real,
    site_indices::Union{Tuple,AbstractVector};
    evolution::Symbol = :real,
)
    return star_gate_itensor(model, site_indices, step; evolution)
end

end
