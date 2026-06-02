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
