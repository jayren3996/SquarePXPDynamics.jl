module ScarFinder

using ..SquareIPEPS: SquareIPEPSState
using ..IPEPSEvolution: TrotterParams, EvolutionLog, evolve!
using ..Observables: SimpleObservableSummary, measure_simple

export ScarFinderParams, ScarFinderIteration, ScarFinderResult, scarfinder!

"""
    ScarFinderParams(
        projection_time,
        trotter,
        iterations,
        max_truncerr,
        max_blockade_violation,
        max_bond_entropy,
        stop_on_reject,
    )

Parameters for the S6-lite ScarFinder orchestration loop. Each iteration
evolves the supplied iPEPS state with [`evolve!`](@ref), records
[`measure_simple`](@ref) diagnostics, and accepts or rejects the iteration
using simple-update/local observables only. This does not run CTMRG, perform
energy targeting, rank candidates, or apply any imaginary-time correction.
"""
struct ScarFinderParams
    projection_time::Float64
    trotter::TrotterParams
    iterations::Int
    max_truncerr::Float64
    max_blockade_violation::Float64
    max_bond_entropy::Float64
    stop_on_reject::Bool

    function ScarFinderParams(
        projection_time::Real,
        trotter::TrotterParams,
        iterations::Integer,
        max_truncerr::Real,
        max_blockade_violation::Real,
        max_bond_entropy::Real,
        stop_on_reject::Bool,
    )
        time = Float64(projection_time)
        time >= 0 || throw(ArgumentError("projection_time must be nonnegative"))
        niterations = Int(iterations)
        niterations >= 0 || throw(ArgumentError("iterations must be nonnegative"))
        truncerr_limit = _nonnegative_float(max_truncerr, "max_truncerr")
        blockade_limit = _nonnegative_float(max_blockade_violation, "max_blockade_violation")
        entropy_limit = _nonnegative_float(max_bond_entropy, "max_bond_entropy")

        return new(
            time,
            trotter,
            niterations,
            truncerr_limit,
            blockade_limit,
            entropy_limit,
            stop_on_reject,
        )
    end
end

"""
    ScarFinderIteration

Diagnostics for one S6-lite ScarFinder iteration. `evolution` is the
[`EvolutionLog`](@ref) returned by [`evolve!`](@ref), and `observables` is the
[`SimpleObservableSummary`](@ref) returned by [`measure_simple`](@ref).
Rejected iterations carry a short deterministic `reject_reason`.
"""
struct ScarFinderIteration
    iteration::Int
    accepted::Bool
    reject_reason::Union{Nothing,String}
    evolution::EvolutionLog
    observables::SimpleObservableSummary
end

"""
    ScarFinderResult

Result of [`scarfinder!`](@ref). `state` is the same mutably evolved
[`SquareIPEPSState`](@ref) passed to `scarfinder!`; iteration records contain
the local/simple diagnostics used for acceptance. These diagnostics are not
CTMRG-quality environment measurements.
"""
struct ScarFinderResult
    state::SquareIPEPSState
    params::ScarFinderParams
    iterations::Vector{ScarFinderIteration}
    accepted_iterations::Int
    rejected_iterations::Int
end

function _nonnegative_float(value::Real, name::String)
    converted = Float64(value)
    converted >= 0 || throw(ArgumentError("$name must be nonnegative"))
    return converted
end

function _finite_summary(obs::SimpleObservableSummary)
    return all(isfinite, (
        obs.density,
        obs.density_even,
        obs.density_odd,
        obs.blockade_violation,
        obs.pxp_energy_density,
        obs.mean_bond_entropy,
        obs.max_bond_entropy,
    ))
end

function _evaluate_scarfinder_iteration(
    log::EvolutionLog,
    obs::SimpleObservableSummary,
    params::ScarFinderParams,
)
    all(isfinite, (log.max_truncerr, log.max_bond_entropy, log.mean_bond_entropy)) ||
        return false, "non-finite evolution diagnostic"
    _finite_summary(obs) || return false, "non-finite simple observable diagnostic"

    if log.max_truncerr > params.max_truncerr
        return false, "max_truncerr exceeds threshold"
    elseif obs.blockade_violation > params.max_blockade_violation
        return false, "blockade_violation exceeds threshold"
    elseif obs.max_bond_entropy > params.max_bond_entropy
        return false, "max_bond_entropy exceeds threshold"
    else
        return true, nothing
    end
end

_count_accepted(iterations) = count(iteration -> iteration.accepted, iterations)
_count_rejected(iterations) = count(iteration -> !iteration.accepted, iterations)

"""
    scarfinder!(psi::SquareIPEPSState, params::ScarFinderParams)::ScarFinderResult

Run the S6-lite ScarFinder projection scaffold in place on `psi`.
For each requested iteration, this calls [`evolve!`](@ref) for
`params.projection_time`, calls [`measure_simple`](@ref), records diagnostics,
and accepts or rejects the iteration by truncation error, simple blockade
violation, simple bond entropy, and finite-diagnostic checks. `psi` is mutated
like [`evolve!`](@ref). No CTMRG, energy correction, candidate ranking, or
direct star-projection logic is performed here.
"""
function scarfinder!(
    psi::SquareIPEPSState,
    params::ScarFinderParams,
)::ScarFinderResult
    iterations = ScarFinderIteration[]

    for n in 1:params.iterations
        log = evolve!(psi, params.projection_time; params = params.trotter)
        obs = measure_simple(psi)
        accepted, reason = _evaluate_scarfinder_iteration(log, obs, params)
        push!(iterations, ScarFinderIteration(n, accepted, reason, log, obs))

        if !accepted && params.stop_on_reject
            break
        end
    end

    return ScarFinderResult(
        psi,
        params,
        iterations,
        _count_accepted(iterations),
        _count_rejected(iterations),
    )
end

"""
    scarfinder!(
        psi::SquareIPEPSState;
        projection_time,
        trotter,
        iterations,
        max_truncerr = Inf,
        max_blockade_violation = Inf,
        max_bond_entropy = Inf,
        stop_on_reject = false,
    )::ScarFinderResult

Convenience keyword constructor for [`ScarFinderParams`](@ref), then delegates
to [`scarfinder!(psi, params)`](@ref). The loop records only simple/local
diagnostics from [`measure_simple`](@ref).
"""
function scarfinder!(
    psi::SquareIPEPSState;
    projection_time::Real,
    trotter::TrotterParams,
    iterations::Integer,
    max_truncerr::Real = Inf,
    max_blockade_violation::Real = Inf,
    max_bond_entropy::Real = Inf,
    stop_on_reject::Bool = false,
)::ScarFinderResult
    params = ScarFinderParams(
        projection_time,
        trotter,
        iterations,
        max_truncerr,
        max_blockade_violation,
        max_bond_entropy,
        stop_on_reject,
    )
    return scarfinder!(psi, params)
end

end
