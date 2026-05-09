module Evolution

using LinearAlgebra: norm
using ITensors: dim
using ..Geometry: Coord
using ..States: TriangularIPEPS, unit_cell_representatives
using ..Gates: dense_gate, projected_gate
using ..Models: pxp_star_hamiltonian
using ..SpinOps: pauli_x, pauli_z, projector_down, projector_up
using ..Observables: local_expectation, mean_blockade_violation, tensor_norm
using ..Schedules: first_order_colors, second_order_colors, schedule_layers
using ..SimpleUpdate: SimpleUpdateDiagnostics, apply_star_gate_simple_update!

export ProjectedPXPStepDiagnostics
export evolve_step!, color_canonical_center
export projected_pxp_step!, imaginary_projected_pxp_step!, run_projected_pxp!

struct ProjectedPXPStepDiagnostics
    layer_diagnostics::Vector{SimpleUpdateDiagnostics}
    discarded_weights::Vector{Float64}
    max_bond_dim::Int
    mean_bond_dim::Float64
    lambda_summaries::Dict{Tuple{Coord,Int},NamedTuple{(:min, :max, :norm),Tuple{Float64,Float64,Float64}}}
    blockade_violation::Float64
    tensor_norms::Vector{Float64}
    local_z::Dict{Coord,Float64}
    local_x::Dict{Coord,Float64}
    local_projector_up::Dict{Coord,Float64}
end

"""
    color_canonical_center(color) -> Coord

Return a canonical center coordinate whose `star_color` equals `color`.
Used by `evolve_step!` to pick a representative center per color.
"""
function color_canonical_center(color::Integer)
    1 <= color <= 7 || throw(ArgumentError("color must be in 1:7"))
    # star_color(c) = (q + 3r) mod 7 + 1  →  pick (color-1, 0).
    return Coord(color - 1, 0)
end

"""
    evolve_step!(state, gate; order=:second, update=:simple) -> state

Apply one Trotter "step" to `state` by sweeping the star gate over all 7
colors of the triangular partition. `order` selects the schedule
(`:first` or `:second`); `update` selects the bond-update backend
(`:simple` is the only currently supported value).

The same `gate` is applied at each scheduled color; for translationally
invariant states the schedule reduces to repeated applications at
color-canonical centers.
"""
function evolve_step!(state::TriangularIPEPS, gate::AbstractMatrix;
                      order::Symbol = :second,
                      update::Symbol = :simple,
                      maxdim::Union{Nothing,Integer} = nothing,
                      cutoff::Real = 1e-12,
                      chi::Union{Nothing,Integer} = nothing,
                      maxiter::Integer = 4,
                      tol::Real = 1e-10,
                      regularization::Real = 1e-10)
    schedule = if order === :first
        first_order_colors()
    elseif order === :second
        second_order_colors()
    else
        throw(ArgumentError("order must be :first or :second"))
    end

    config = _update_config(update, maxdim, cutoff, chi, maxiter, tol, regularization)
    layer_gate = order === :second ? _half_step_gate(gate) : Matrix{ComplexF64}(gate)

    for color in schedule
        center = color_canonical_center(color)
        _apply_update!(state, layer_gate, center, update, config; cutoff = cutoff, maxdim = maxdim)
    end
    return state
end

function evolve_step!(state::TriangularIPEPS,
                      H::AbstractMatrix,
                      dt::Real;
                      order::Symbol = :second,
                      update::Symbol = :simple,
                      evolution::Symbol = :real,
                      projected::Bool = false,
                      projector::Union{Nothing,AbstractMatrix} = nothing,
                      maxdim::Union{Nothing,Integer} = nothing,
                      cutoff::Real = 1e-12,
                      chi::Union{Nothing,Integer} = nothing,
                      maxiter::Integer = 4,
                      tol::Real = 1e-10,
                      regularization::Real = 1e-10)
    config = _update_config(update, maxdim, cutoff, chi, maxiter, tol, regularization)

    for layer in schedule_layers(order)
        step = dt * layer.scale
        gate = _cached_gate(H, step, evolution, projected, projector)
        _apply_update!(state, gate, color_canonical_center(layer.color), update, config;
                       cutoff = cutoff, maxdim = maxdim)
    end
    return state
end

function _half_step_gate(gate::AbstractMatrix)
    size(gate) == (128, 128) || throw(ArgumentError("gate must be 128x128"))
    return sqrt(Matrix{ComplexF64}(gate))
end

"""
    projected_pxp_step!(state, dt; order=:second, maxdim, cutoff=1e-12, evolution=:real)

Apply one scheduled projected-PXP sweep. For `order=:second`, `dt` is the
full step parameter and each of the 14 symmetric layers uses `dt/2`; over a
complete translationally invariant product-gate sweep this gives the same
per-site total angle as the 7-layer first-order schedule.
"""
function projected_pxp_step!(state::TriangularIPEPS,
                             dt::Real;
                             order::Symbol = :second,
                             maxdim::Union{Nothing,Integer} = nothing,
                             cutoff::Real = 1e-12,
                             evolution::Symbol = :real,
                             update::Symbol = :simple,
                             chi::Union{Nothing,Integer} = nothing,
                             maxiter::Integer = 4,
                             tol::Real = 1e-10,
                             regularization::Real = 1e-10)
    _validate_projected_pxp_options(order, maxdim, cutoff, evolution)
    config = _update_config(update, maxdim, cutoff, chi, maxiter, tol, regularization)
    H = pxp_star_hamiltonian(projector_down(), pauli_x())
    layer_diags = SimpleUpdateDiagnostics[]
    for layer in schedule_layers(order)
        gate = _cached_gate(H, dt * layer.scale, evolution, true, nothing)
        push!(layer_diags, _apply_update!(
            state, gate, color_canonical_center(layer.color), update, config;
            cutoff = cutoff, maxdim = maxdim,
        ))
    end
    return _step_diagnostics(state, layer_diags)
end

function imaginary_projected_pxp_step!(state::TriangularIPEPS,
                                       dτ::Real;
                                       order::Symbol = :second,
                                       maxdim::Union{Nothing,Integer} = nothing,
                                       cutoff::Real = 1e-12,
                                       update::Symbol = :simple,
                                       chi::Union{Nothing,Integer} = nothing,
                                       maxiter::Integer = 4,
                                       tol::Real = 1e-10,
                                       regularization::Real = 1e-10)
    return projected_pxp_step!(
        state, dτ; order = order, maxdim = maxdim, cutoff = cutoff,
        evolution = :imaginary, update = update, chi = chi, maxiter = maxiter,
        tol = tol, regularization = regularization,
    )
end

function run_projected_pxp!(state::TriangularIPEPS,
                            dt::Real,
                            nsteps::Integer;
                            order::Symbol = :second,
                            maxdim::Union{Nothing,Integer} = nothing,
                            cutoff::Real = 1e-12,
                            evolution::Symbol = :real,
                            update::Symbol = :simple,
                            chi::Union{Nothing,Integer} = nothing,
                            maxiter::Integer = 4,
                            tol::Real = 1e-10,
                            regularization::Real = 1e-10)
    nsteps >= 0 || throw(ArgumentError("nsteps must be nonnegative"))
    history = ProjectedPXPStepDiagnostics[]
    for _ in 1:nsteps
        push!(history, projected_pxp_step!(
            state, dt; order = order, maxdim = maxdim, cutoff = cutoff,
            evolution = evolution, update = update, chi = chi, maxiter = maxiter,
            tol = tol, regularization = regularization,
        ))
    end
    return history
end

const _GATE_CACHE = Dict{Tuple{UInt,Float64,Symbol,Bool,UInt},Matrix{ComplexF64}}()

function _cached_gate(H::AbstractMatrix,
                      step::Real,
                      evolution::Symbol,
                      projected::Bool,
                      projector::Union{Nothing,AbstractMatrix})
    key = (objectid(H), Float64(step), evolution, projected,
           projector === nothing ? UInt(0) : objectid(projector))
    return get!(_GATE_CACHE, key) do
        if projected
            projector === nothing ?
                projected_gate(H, step; evolution = evolution) :
                projected_gate(H, step; evolution = evolution, projector = projector)
        else
            dense_gate(H, step; evolution = evolution)
        end
    end
end

function _update_config(update::Symbol,
                        maxdim::Union{Nothing,Integer},
                        cutoff::Real,
                        chi::Union{Nothing,Integer},
                        maxiter::Integer,
                        tol::Real,
                        regularization::Real)
    update === :simple || throw(ArgumentError("update must be :simple"))
    return nothing
end

function _apply_update!(state::TriangularIPEPS,
                        gate::AbstractMatrix,
                        center::Coord,
                        update::Symbol,
                        config;
                        cutoff::Real,
                        maxdim::Union{Nothing,Integer})
    update === :simple || throw(ArgumentError("update must be :simple"))
    return apply_star_gate_simple_update!(
        state, gate, center; maxdim = maxdim, cutoff = cutoff,
    )
end

function _validate_projected_pxp_options(order::Symbol,
                                         maxdim::Union{Nothing,Integer},
                                         cutoff::Real,
                                         evolution::Symbol)
    order in (:first, :second) || throw(ArgumentError("order must be :first or :second"))
    evolution in (:real, :imaginary) || throw(ArgumentError("evolution must be :real or :imaginary"))
    maxdim !== nothing || throw(ArgumentError("maxdim is required"))
    maxdim >= 1 || throw(ArgumentError("maxdim must be >= 1"))
    cutoff >= 0 || throw(ArgumentError("cutoff must be nonnegative"))
    return nothing
end

function _step_diagnostics(state::TriangularIPEPS,
                           layer_diags::Vector{SimpleUpdateDiagnostics})
    reps = collect(unit_cell_representatives(state.unitcell))
    bond_dims = [dim(idx) for idx in values(state.bond_inds)]
    lambda_summaries = Dict{Tuple{Coord,Int},NamedTuple{(:min, :max, :norm),Tuple{Float64,Float64,Float64}}}()
    for (bond, λ) in state.lambdas
        lambda_summaries[bond] = (min = minimum(λ), max = maximum(λ), norm = norm(λ))
    end
    local_z = Dict(c => real(local_expectation(state, c, pauli_z())) for c in reps)
    local_x = Dict(c => real(local_expectation(state, c, pauli_x())) for c in reps)
    local_pup = Dict(c => real(local_expectation(state, c, projector_up())) for c in reps)
    return ProjectedPXPStepDiagnostics(
        layer_diags,
        [d.discarded_weight for d in layer_diags],
        maximum(bond_dims),
        sum(bond_dims) / length(bond_dims),
        lambda_summaries,
        mean_blockade_violation(state, reps),
        [tensor_norm(state, c) for c in reps],
        local_z,
        local_x,
        local_pup,
    )
end

end
