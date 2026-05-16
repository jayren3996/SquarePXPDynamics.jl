module IPEPSEvolution

using ..SquareUnitCells:
    assert_five_color_compatible, update_centers, stars_are_disjoint_mod_unitcell
using ..SquareIPEPS: SquareIPEPSState, all_bond_entropies, log_norm
using ..StarSimpleUpdate: StarUpdateInfo, project_star!
using ..StarModels: AbstractModelProtocol, PXPStarModel, StaticModel, model_at

export TrotterParams, EvolutionLog, trotter_sequence, evolve!

const _TROTTER_SPLIT_DIRECTIONS = (:right, :up, :left, :down)

"""
    TrotterParams(
        dt,
        order,
        evolution,
        maxdim,
        cutoff,
        split_order = (:right, :up, :left, :down);
        schedule = :five_color,
    )

Parameters for deterministic iPEPS Trotter evolution. `dt` is the positive
full time-step size, `order` must be `1` or `2`, `evolution` must be `:real` or
`:imaginary`, `schedule` must be `:five_color` or `:serial`, and
`maxdim`/`cutoff`/`split_order` are forwarded to [`project_star!`](@ref).
The default `:five_color` schedule applies disjoint star layers. The `:serial`
schedule applies one unit-cell center per layer in representative order, with
second order using the corresponding reversed half sweep.

The legacy six-argument PXP constructor
`TrotterParams(dt, order, evolution, projected, maxdim, cutoff)` returns a
compatibility wrapper accepted by `evolve!`, `trotter_sequence`, and existing
PXP orchestration APIs. New model-aware code should keep `TrotterParams`
model-agnostic and pass an explicit model protocol to `evolve!`.
"""
struct TrotterParams
    dt::Float64
    order::Int
    evolution::Symbol
    maxdim::Int
    cutoff::Float64
    split_order::NTuple{4,Symbol}
    schedule::Symbol

    function TrotterParams(
        dt::Real,
        order::Integer,
        evolution::Symbol,
        maxdim::Integer,
        cutoff::Real,
        split_order = _TROTTER_SPLIT_DIRECTIONS;
        schedule::Symbol = :five_color,
    )
        step = Float64(dt)
        isfinite(step) && step > 0 || throw(ArgumentError("dt must be finite and positive"))
        ord = Int(order)
        ord in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
        evolution in (:real, :imaginary) ||
            throw(ArgumentError("evolution must be :real or :imaginary"))
        dim = Int(maxdim)
        dim >= 1 || throw(ArgumentError("maxdim must be at least 1"))
        trunc_cutoff = Float64(cutoff)
        isfinite(trunc_cutoff) && trunc_cutoff >= 0 ||
            throw(ArgumentError("cutoff must be finite and nonnegative"))
        order_values = _validate_trotter_split_order(split_order)
        schedule in (:five_color, :serial) ||
            throw(ArgumentError("schedule must be :five_color or :serial"))
        return new(step, ord, evolution, dim, trunc_cutoff, order_values, schedule)
    end
end

struct LegacyPXPParams
    trotter::TrotterParams
    protocol::StaticModel{PXPStarModel}
end

"""
    legacy_trotter_params(params)

Return the model-agnostic [`TrotterParams`](@ref) from either a current
`TrotterParams` value or the legacy PXP compatibility wrapper returned by the
old six-argument constructor.
"""
legacy_trotter_params(params::TrotterParams) = params
legacy_trotter_params(params::LegacyPXPParams) = params.trotter
legacy_trotter_params(params) = throw(
    ArgumentError("trotter must be a TrotterParams or legacy PXP TrotterParams"),
)
legacy_trotter_protocol(params::TrotterParams) = nothing
legacy_trotter_protocol(params::LegacyPXPParams) = params.protocol
legacy_trotter_protocol(params) = legacy_trotter_params(params)

function Base.getproperty(params::LegacyPXPParams, name::Symbol)
    if name === :projected
        return getfield(getfield(params, :protocol), :model).projected
    elseif name in (:dt, :order, :evolution, :maxdim, :cutoff, :split_order, :schedule)
        return getproperty(getfield(params, :trotter), name)
    else
        return getfield(params, name)
    end
end

function Base.propertynames(params::LegacyPXPParams; private::Bool = false)
    public = (:dt, :order, :evolution, :maxdim, :cutoff, :split_order, :schedule, :projected)
    return private ? (public..., :trotter, :protocol) : public
end

function TrotterParams(
    dt::Real,
    order::Integer,
    evolution::Symbol,
    projected::Bool,
    maxdim::Integer,
    cutoff::Real,
)
    return LegacyPXPParams(
        TrotterParams(dt, order, evolution, maxdim, cutoff),
        StaticModel(PXPStarModel(projected)),
    )
end

function _validate_trotter_split_order(split_order)
    order_values = try
        collect(split_order)
    catch
        throw(
            ArgumentError(
                "split_order must be a permutation of (:right, :up, :left, :down)",
            ),
        )
    end
    all(dir -> dir isa Symbol, order_values) || throw(
        ArgumentError("split_order must be a permutation of (:right, :up, :left, :down)"),
    )
    order = Tuple(order_values)
    length(order) == length(_TROTTER_SPLIT_DIRECTIONS) ||
        throw(ArgumentError("split_order must contain four directions"))
    all(dir -> dir in _TROTTER_SPLIT_DIRECTIONS, order) || throw(
        ArgumentError("split_order must be a permutation of (:right, :up, :left, :down)"),
    )
    Set(order) == Set(_TROTTER_SPLIT_DIRECTIONS) || throw(
        ArgumentError("split_order must be a permutation of (:right, :up, :left, :down)"),
    )
    return order::NTuple{4,Symbol}
end

"""
    EvolutionLog

Diagnostics returned by [`evolve!`](@ref). `layer_infos` stores every
[`StarUpdateInfo`](@ref) grouped by applied Trotter layer, `total_time` is the
requested applied time, `model_metadata` records the protocol/model identity
used for the call, and the entropy fields summarize the final simple-update
link weights. `log_norm_before`, `log_norm_after`, and `log_norm_delta` record
the normalization ledger around the evolution call.
"""
struct EvolutionLog
    total_time::Float64
    params::TrotterParams
    model_metadata::NamedTuple
    nsteps::Int
    layer_infos::Vector{Vector{StarUpdateInfo}}
    max_truncerr::Float64
    max_bond_entropy::Float64
    mean_bond_entropy::Float64
    log_norm_before::Float64
    log_norm_after::Float64
    log_norm_delta::Float64
end

"""
    trotter_sequence(params::TrotterParams)::Vector{Tuple{Int,Float64}}

Return the deterministic five-color Trotter layer sequence for `params`. First
order sweeps colors `1:5` with full `dt`; second order performs a symmetric
`1,2,3,4,5,4,3,2,1` sweep where color 5 receives the full `dt` and all other
layers receive `dt / 2`.
"""
function trotter_sequence(params::TrotterParams)::Vector{Tuple{Int,Float64}}
    if params.order == 1
        return [(color, params.dt) for color = 1:5]
    elseif params.order == 2
        half_step = params.dt / 2
        return [
            (1, half_step),
            (2, half_step),
            (3, half_step),
            (4, half_step),
            (5, params.dt),
            (4, half_step),
            (3, half_step),
            (2, half_step),
            (1, half_step),
        ]
    else
        throw(ArgumentError("order must be 1 or 2"))
    end
end

trotter_sequence(params::LegacyPXPParams)::Vector{Tuple{Int,Float64}} =
    trotter_sequence(params.trotter)

function _nsteps_for_total_time(total_time::Real, params::TrotterParams)
    total = Float64(total_time)
    isfinite(total) || throw(ArgumentError("total_time must be finite"))
    total == 0 && return 0, total
    total > 0 || throw(ArgumentError("total_time must be nonnegative"))

    nsteps_float = total / params.dt
    nsteps = round(Int, nsteps_float)
    isapprox(nsteps_float, nsteps; rtol = 1e-10, atol = 1e-12) ||
        throw(ArgumentError("total_time must be an integer multiple of params.dt"))
    return nsteps, total
end

function _assert_finite_diagnostics(max_truncerr, max_bond_entropy, mean_bond_entropy)
    all(isfinite, (max_truncerr, max_bond_entropy, mean_bond_entropy)) ||
        throw(ArgumentError("evolution diagnostics must be finite"))
    return nothing
end

function _assert_finite_log_norms(log_norm_before, log_norm_after, log_norm_delta)
    all(isfinite, (log_norm_before, log_norm_after, log_norm_delta)) ||
        throw(ArgumentError("log-norm diagnostics must be finite"))
    return nothing
end

function _assert_finite_star_update_info(info::StarUpdateInfo)
    isfinite(info.max_truncerr) && info.max_truncerr >= 0 ||
        throw(ArgumentError("star update max_truncerr must be finite and nonnegative"))
    all(isfinite, values(info.truncerrs)) ||
        throw(ArgumentError("star update truncation errors must be finite"))
    all(>=(0), values(info.truncerrs)) ||
        throw(ArgumentError("star update truncation errors must be nonnegative"))
    all(>(0), values(info.keptdims)) ||
        throw(ArgumentError("star update kept dimensions must be positive"))
    all(isfinite, values(info.min_lambda)) ||
        throw(ArgumentError("star update minimum link weights must be finite"))
    all(>=(0), values(info.min_lambda)) ||
        throw(ArgumentError("star update minimum link weights must be nonnegative"))
    all(isfinite, values(info.touched_min_lambda)) ||
        throw(ArgumentError("star update touched minimum link weights must be finite"))
    all(>=(0), values(info.touched_min_lambda)) ||
        throw(ArgumentError("star update touched minimum link weights must be nonnegative"))
    all(isfinite, values(info.norm_factors)) ||
        throw(ArgumentError("star update norm factors must be finite"))
    all(>(0), values(info.norm_factors)) ||
        throw(ArgumentError("star update norm factors must be positive"))
    return nothing
end

function _assert_finite_star_update_infos(layer_infos)
    for info in Iterators.flatten(layer_infos)
        _assert_finite_star_update_info(info)
    end
    return nothing
end

function _model_metadata(protocol::AbstractModelProtocol)
    model = model_at(protocol, 0.0, 1)
    return (
        protocol_type = string(nameof(typeof(protocol))),
        model_type = string(nameof(typeof(model))),
        pxp_projected = model isa PXPStarModel ? model.projected : nothing,
    )
end

function _apply_star_layer!(
    psi::SquareIPEPSState,
    center,
    layer_dt::Float64,
    params::TrotterParams,
    model,
)
    return project_star!(
        psi,
        center,
        layer_dt;
        model = model,
        evolution = params.evolution,
        maxdim = params.maxdim,
        cutoff = params.cutoff,
        split_order = params.split_order,
    )
end

function _finish_evolution_log(
    psi::SquareIPEPSState,
    total_time::Float64,
    params::TrotterParams,
    model_metadata::NamedTuple,
    nsteps::Int,
    layer_infos,
    log_norm_before::Float64,
)
    max_truncerr =
        isempty(layer_infos) ? 0.0 :
        maximum(info.max_truncerr for info in Iterators.flatten(layer_infos))
    _assert_finite_star_update_infos(layer_infos)
    entropy_values = collect(values(all_bond_entropies(psi)))
    max_bond_entropy = maximum(entropy_values)
    mean_bond_entropy = sum(entropy_values) / length(entropy_values)
    _assert_finite_diagnostics(max_truncerr, max_bond_entropy, mean_bond_entropy)
    log_norm_after = log_norm(psi)
    log_norm_delta = log_norm_after - log_norm_before
    _assert_finite_log_norms(log_norm_before, log_norm_after, log_norm_delta)

    return EvolutionLog(
        total_time,
        params,
        model_metadata,
        nsteps,
        layer_infos,
        max_truncerr,
        max_bond_entropy,
        mean_bond_entropy,
        log_norm_before,
        log_norm_after,
        log_norm_delta,
    )
end

function _evolve_five_color_step!(
    psi::SquareIPEPSState,
    params::TrotterParams,
    model,
    layer_infos,
)
    for (color, layer_dt) in trotter_sequence(params)
        centers = update_centers(psi.unitcell, color)
        stars_are_disjoint_mod_unitcell(psi.unitcell, centers) ||
            throw(ArgumentError("color layer $color has overlapping wrapped stars"))
        infos = StarUpdateInfo[]
        for center in centers
            push!(infos, _apply_star_layer!(psi, center, layer_dt, params, model))
        end
        push!(layer_infos, infos)
    end
    return nothing
end

function _serial_trotter_sequence(centers, params::TrotterParams)
    if params.order == 1
        return [(center, params.dt) for center in centers]
    elseif params.order == 2
        half_step = params.dt / 2
        forward = [
            (center, idx == length(centers) ? params.dt : half_step) for
            (idx, center) in enumerate(centers)
        ]
        reverse_half = [(center, half_step) for center in reverse(centers[1:(end - 1)])]
        return [forward; reverse_half]
    else
        throw(ArgumentError("order must be 1 or 2"))
    end
end

function _evolve_serial_step!(
    psi::SquareIPEPSState,
    params::TrotterParams,
    model,
    layer_infos,
)
    for (center, layer_dt) in _serial_trotter_sequence(psi.unitcell.reps, params)
        info = _apply_star_layer!(psi, center, layer_dt, params, model)
        push!(layer_infos, [info])
    end
    return nothing
end

function _evolve_with_params!(
    psi::SquareIPEPSState,
    total_time::Real,
    params::TrotterParams,
    protocol::AbstractModelProtocol,
)::EvolutionLog
    nsteps, total = _nsteps_for_total_time(total_time, params)
    log_norm_before = log_norm(psi)
    metadata = _model_metadata(protocol)
    layer_infos = Vector{Vector{StarUpdateInfo}}()
    nsteps == 0 &&
        return _finish_evolution_log(
            psi,
            total,
            params,
            metadata,
            nsteps,
            layer_infos,
            log_norm_before,
        )

    params.schedule === :five_color && assert_five_color_compatible(psi.unitcell)
    for step_index = 1:nsteps
        time_before_step = (step_index - 1) * params.dt
        model = model_at(protocol, time_before_step, step_index)
        if params.schedule === :five_color
            _evolve_five_color_step!(psi, params, model, layer_infos)
        elseif params.schedule === :serial
            _evolve_serial_step!(psi, params, model, layer_infos)
        else
            throw(ArgumentError("schedule must be :five_color or :serial"))
        end
    end

    return _finish_evolution_log(
        psi,
        total,
        params,
        metadata,
        nsteps,
        layer_infos,
        log_norm_before,
    )
end

"""
    evolve!(
        psi::SquareIPEPSState,
        total_time::Real;
        params::TrotterParams,
        protocol = nothing,
    )::EvolutionLog

    evolve!(
        psi::SquareIPEPSState,
        total_time::Real;
        dt::Real,
        order::Integer = 2,
        evolution::Symbol = :real,
        projected::Bool = true,
        maxdim::Integer = psi.maxdim,
        cutoff::Real = 1e-12,
        schedule::Symbol = :five_color,
        protocol = nothing,
    )::EvolutionLog

Apply deterministic iPEPS Trotter evolution by orchestrating square-star
updates and calling [`project_star!`](@ref) for each center. `schedule` may be
`:five_color` for disjoint color layers or `:serial` for one center per layer.
`total_time` must be zero or a positive integer multiple of `dt`. Zero time
returns diagnostics without mutating `psi`. The convenience keyword form
constructs [`TrotterParams`](@ref) and then uses the same evolution path.
"""
function evolve!(
    psi::SquareIPEPSState,
    total_time::Real;
    params::Union{TrotterParams,LegacyPXPParams,Nothing} = nothing,
    protocol::Union{AbstractModelProtocol,Nothing} = nothing,
    dt = nothing,
    order::Integer = 2,
    evolution::Symbol = :real,
    projected::Bool = true,
    maxdim::Integer = psi.maxdim,
    cutoff::Real = 1e-12,
    schedule::Symbol = :five_color,
)::EvolutionLog
    actual_params, actual_protocol = if params === nothing
        dt === nothing && throw(UndefKeywordError(:dt))
        (
            TrotterParams(dt, order, evolution, maxdim, cutoff; schedule = schedule),
            protocol === nothing ? StaticModel(PXPStarModel(projected)) : protocol,
        )
    elseif params isa LegacyPXPParams
        protocol === nothing || throw(
            ArgumentError("cannot pass an explicit protocol with legacy PXP TrotterParams"),
        )
        (params.trotter, params.protocol)
    else
        (
            params,
            protocol === nothing ? StaticModel(PXPStarModel(projected)) : protocol,
        )
    end
    return _evolve_with_params!(psi, total_time, actual_params, actual_protocol)
end

end
