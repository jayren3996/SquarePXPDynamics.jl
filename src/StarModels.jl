module StarModels

using ITensors
using LinearAlgebra
import ..SpinOps: pauli_x, pauli_z, identity2, kron_all, embed_one_site
import ..SquarePXP:
    SQUARE_STAR_SITES,
    square_pxp_star_hamiltonian,
    square_pxp_gate,
    projected_square_pxp_gate

export AbstractStarModel, PXPStarModel, TFIMStarModel
export AbstractModelProtocol, StaticModel, model_at
export star_site_order, tfim_pauli_convention
export star_hamiltonian, star_gate, star_gate_itensor, tfim_product_basis_energy

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
    TFIMStarModel(J, h)

Five-site square-star transverse-field Ising model with finite coupling `J`
and finite transverse field `h`. Values are promoted to a common real type.
The Hamiltonian convention is
`-h * X_center - (J / 2) * Z_center * (Z_right + Z_up + Z_left + Z_down)`.
"""
struct TFIMStarModel{T<:Real} <: AbstractStarModel
    J::T
    h::T

    function TFIMStarModel{T}(J::T, h::T) where {T<:Real}
        isfinite(J) || throw(ArgumentError("J must be finite"))
        isfinite(h) || throw(ArgumentError("h must be finite"))
        return new{T}(J, h)
    end
end

TFIMStarModel(J::Real, h::Real) =
    TFIMStarModel{promote_type(typeof(J), typeof(h))}(promote(J, h)...)

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
    tfim_pauli_convention()

Return the TFIM Pauli convention: `Z |up> = +|up>` and transverse `X` field.
"""
tfim_pauli_convention() = (:Z_up_is_plus_one, :X_field)

function _validate_finite_step(step::Real)
    isfinite(step) || throw(ArgumentError("step must be finite"))
    return step
end

function _validate_evolution(evolution::Symbol)
    evolution === :real || evolution === :imaginary ||
        throw(ArgumentError("evolution must be :real or :imaginary"))
    return evolution
end

"""
    star_hamiltonian(model)

Return the dense 32x32 Hamiltonian for a five-site square-star model in
`star_site_order()`.
"""
star_hamiltonian(model::PXPStarModel) = square_pxp_star_hamiltonian()

function star_hamiltonian(model::TFIMStarModel)
    nsites = SQUARE_STAR_SITES
    z_center = embed_one_site(pauli_z(), 1, nsites)
    H = -model.h * embed_one_site(pauli_x(), 1, nsites)
    for site = 2:nsites
        H .-= (model.J / 2) .* (z_center * embed_one_site(pauli_z(), site, nsites))
    end
    return H
end

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

function star_gate(model::TFIMStarModel, step::Real; evolution::Symbol = :real)
    _validate_finite_step(step)
    _validate_evolution(evolution)
    H = star_hamiltonian(model)
    return evolution === :real ? exp(-im * step * H) : exp(-step * H)
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

function _tfim_z_value(state)
    state === :up && return 1
    state === :z_up && return 1
    state isa Integer && !(state isa Bool) && state == 1 && return 1
    state === :down && return -1
    state === :z_down && return -1
    state isa Integer && !(state isa Bool) && state == 2 && return -1
    state isa Integer && !(state isa Bool) && state == -1 && return -1
    throw(ArgumentError("TFIM product basis states must be :up, :down, :z_up, :z_down, 1, 2, +1, or -1"))
end

"""
    tfim_product_basis_energy(model, states)

Return the diagonal TFIM bond energy
`-(J / 2) * Z_center * (Z_right + Z_up + Z_left + Z_down)` for five product
basis labels in `star_site_order()`. Accepted labels include `:up`, `:down`,
`:z_up`, `:z_down`, `1`, `2`, `+1`, and `-1`.
"""
function tfim_product_basis_energy(model::TFIMStarModel, states)
    length(states) == SQUARE_STAR_SITES ||
        throw(ArgumentError("TFIM square-star product basis must have 5 states"))
    z_center = _tfim_z_value(states[1])
    z_neighbors = sum(_tfim_z_value(states[site]) for site = 2:SQUARE_STAR_SITES)
    return -(model.J / 2) * z_center * z_neighbors
end

end
