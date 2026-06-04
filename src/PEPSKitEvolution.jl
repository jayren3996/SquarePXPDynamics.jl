module PEPSKitEvolution

# PEPSKit-native Trotter evolution for the square PXP star dynamics. Preserves
# the legacy five-colour scheduling and palindromic second-order sequence; the
# per-star kernel is the PEPSKit-native simple update `project_star_pepskit!`.
# This deliberately does NOT use PEPSKit's generic `time_evolve`, which cannot
# represent the five-site cross gate.

using ..SquareGeometry: SquareCoord
using ..SquareUnitCells:
    PeriodicSquareUnitCell, update_centers, assert_five_color_compatible
using ..PEPSKitBackend: PXPIPEPSState, state_version, log_norm
using ..PEPSKitStarUpdate: StarUpdateInfoPK, project_star_pepskit!

export TrotterParamsPK, trotter_sequence_pk, evolve_pepskit!, EvolutionLogPK

"""
    TrotterParamsPK

Parameters for one PEPSKit-native Trotter evolution. `order` is 1 or 2,
`schedule` is `:five_color` or `:serial`, `evolution` is `:real` or
`:imaginary`.
"""
struct TrotterParamsPK
    dt::Float64
    order::Int
    schedule::Symbol
    maxdim::Int
    cutoff::Float64
    rel_floor::Float64
    evolution::Symbol
    projected::Bool

    function TrotterParamsPK(;
        dt::Real,
        order::Integer = 2,
        schedule::Symbol = :five_color,
        maxdim::Integer,
        cutoff::Real = 1.0e-12,
        rel_floor::Real = 1.0e-3,
        evolution::Symbol = :real,
        projected::Bool = true,
    )
        isfinite(dt) && dt > 0 || throw(ArgumentError("dt must be finite and positive"))
        order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
        schedule in (:five_color, :serial) ||
            throw(ArgumentError("schedule must be :five_color or :serial"))
        maxdim >= 1 || throw(ArgumentError("maxdim must be at least 1"))
        evolution in (:real, :imaginary) ||
            throw(ArgumentError("evolution must be :real or :imaginary"))
        return new(
            Float64(dt), Int(order), schedule, Int(maxdim),
            Float64(cutoff), Float64(rel_floor), evolution, projected,
        )
    end
end

"""
    EvolutionLogPK

Diagnostics accumulated over a PEPSKit-native evolution run.
"""
struct EvolutionLogPK
    nsteps::Int
    max_truncerr::Float64
    max_bond_dim::Int
    log_norm_before::Float64
    log_norm_after::Float64
end

"""
    trotter_sequence_pk(params)::Vector{Tuple{Int,Float64}}

Return the five-colour Trotter layer sequence as `(color, layer_dt)` pairs.
First order sweeps colours `1:5` with full `dt`; second order is the palindromic
`1,2,3,4,5,4,3,2,1` sweep with colour 5 receiving the full `dt` and the others
`dt/2`. Identical schedule to the legacy `trotter_sequence`.
"""
function trotter_sequence_pk(params::TrotterParamsPK)::Vector{Tuple{Int, Float64}}
    if params.order == 1
        return [(color, params.dt) for color in 1:5]
    else
        h = params.dt / 2
        return [
            (1, h), (2, h), (3, h), (4, h), (5, params.dt),
            (4, h), (3, h), (2, h), (1, h),
        ]
    end
end

# A Trotter step is a list of (center, layer_dt) substeps. For :five_color the
# colour sequence maps onto disjoint star layers; for :serial each representative
# is its own layer in representative order, with the order-2 palindrome reversing
# the half sweep (matching the legacy IPEPSEvolution serial schedule).
function _step_substeps(
    cell::PeriodicSquareUnitCell, params::TrotterParamsPK,
)::Vector{Tuple{SquareCoord, Float64}}
    subs = Tuple{SquareCoord, Float64}[]
    if params.schedule === :five_color
        for (color, layer_dt) in trotter_sequence_pk(params)
            for center in update_centers(cell, color)
                push!(subs, (center, layer_dt))
            end
        end
    else # :serial
        if params.order == 1
            for center in cell.reps
                push!(subs, (center, params.dt))
            end
        else
            h = params.dt / 2
            for center in cell.reps
                push!(subs, (center, h))
            end
            for center in reverse(cell.reps)
                push!(subs, (center, h))
            end
        end
    end
    return subs
end

function _nsteps(total_time::Real, dt::Float64)
    total = Float64(total_time)
    isfinite(total) && total >= 0 || throw(ArgumentError("total_time must be finite and >= 0"))
    total == 0 && return 0
    n = round(Int, total / dt)
    isapprox(total / dt, n; rtol = 1e-10, atol = 1e-12) ||
        throw(ArgumentError("total_time must be an integer multiple of dt"))
    return n
end

"""
    evolve_pepskit!(state, total_time; dt, order = 2, schedule = :five_color,
                    maxdim, cutoff = 1e-12, rel_floor = 1e-3,
                    evolution = :real, projected = true)

Evolve `state` for `total_time` with the PEPSKit-native five-site star simple
update, preserving the legacy five-colour palindromic schedule. Returns an
`EvolutionLogPK`. The `:five_color` schedule requires a five-colour-compatible
unit cell; `:serial` applies every representative once per step.
"""
function evolve_pepskit!(
    state::PXPIPEPSState,
    total_time::Real;
    dt::Real,
    order::Integer = 2,
    schedule::Symbol = :five_color,
    maxdim::Integer,
    cutoff::Real = 1.0e-12,
    rel_floor::Real = 1.0e-3,
    evolution::Symbol = :real,
    projected::Bool = true,
    reverse::Bool = false,
)::EvolutionLogPK
    params = TrotterParamsPK(;
        dt, order, schedule, maxdim, cutoff, rel_floor, evolution, projected,
    )
    schedule === :five_color && assert_five_color_compatible(state.unitcell)
    nsteps = _nsteps(total_time, params.dt)
    substeps = _step_substeps(state.unitcell, params)
    if reverse
        # The exact inverse of one forward step: apply the substeps in inverse
        # order with negated layer times (U(dt)^{-1} = U(-dt)). Requires real time.
        evolution === :real || throw(
            ArgumentError("reverse evolution is defined only for real-time evolution"),
        )
        substeps = [(center, -layer_dt) for (center, layer_dt) in Base.reverse(substeps)]
    end

    log_before = log_norm(state)
    max_truncerr = 0.0
    for _ in 1:nsteps
        for (center, layer_dt) in substeps
            info = project_star_pepskit!(
                state, center;
                dt = layer_dt, maxdim = params.maxdim, cutoff = params.cutoff,
                rel_floor = params.rel_floor, evolution = params.evolution,
                projected = params.projected,
            )
            max_truncerr = max(max_truncerr, info.max_truncerr)
        end
    end

    max_bond_dim = maximum(_weight_dim(w) for w in state.weights.data)
    return EvolutionLogPK(
        nsteps, max_truncerr, max_bond_dim, log_before, log_norm(state),
    )
end

_weight_dim(w) = length(w.data)

end
