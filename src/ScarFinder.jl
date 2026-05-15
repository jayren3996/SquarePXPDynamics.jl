module ScarFinder

using ..SquareIPEPS: SquareIPEPSState
using ..IPEPSEvolution: TrotterParams, EvolutionLog, evolve!
using ..Observables: SimpleObservableSummary, measure_simple
using ..PEPSKitMeasurements: CTMObservableSummary, CTMRGDiagnostics

export ScarFinderParams, ScarFinderCandidateScore, ScarFinderIteration, ScarFinderResult
export rank_scarfinder_candidates, write_scarfinder_log, scarfinder!

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
energy targeting, or apply any imaginary-time correction. Candidate scores are
computed from recorded diagnostics by the orchestration layer.
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
        isfinite(time) || throw(ArgumentError("projection_time must be finite"))
        time >= 0 || throw(ArgumentError("projection_time must be nonnegative"))
        niterations = Int(iterations)
        niterations >= 0 || throw(ArgumentError("iterations must be nonnegative"))
        truncerr_limit = _nonnegative_float(max_truncerr, "max_truncerr")
        blockade_limit =
            _nonnegative_float(max_blockade_violation, "max_blockade_violation")
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
    ScarFinderCandidateScore

One ranked ScarFinder candidate diagnostic record. `diagnostics` is `:simple`
for the default local/simple measurements or `:ctm` for records supplied by an
optional CTM measurement callback. The scalar `score` is an operational sorting
key built from existing diagnostics only; it is not an energy correction or a
physics claim.
"""
struct ScarFinderCandidateScore
    iteration::Int
    diagnostics::Symbol
    accepted::Bool
    reject_reason::Union{Nothing,String}
    density::Float64
    density_even::Float64
    density_odd::Float64
    blockade_violation::Float64
    pxp_energy_density::Float64
    mean_bond_entropy::Union{Nothing,Float64}
    max_bond_entropy::Union{Nothing,Float64}
    max_truncerr::Float64
    score::Float64
    ctm_chi::Union{Nothing,Int}
    ctm_tol::Union{Nothing,Float64}
    ctm_maxiter::Union{Nothing,Int}
    ctm_iterations::Union{Nothing,Int}
    ctm_residual::Union{Nothing,Float64}
    ctm_converged::Union{Nothing,Bool}
    ctm_accepted::Union{Nothing,Bool}
end

"""
    ScarFinderIteration

Diagnostics for one S6-lite ScarFinder iteration. `evolution` is the
[`EvolutionLog`](@ref) returned by [`evolve!`](@ref), and `observables` is the
[`SimpleObservableSummary`](@ref) returned by [`measure_simple`](@ref). The
`simple_score` field is always present. The `ctm_score` field is populated only
when a caller supplies a CTM measurement callback and the callback is scheduled
for that iteration. Rejected iterations carry a short deterministic
`reject_reason`.
"""
struct ScarFinderIteration
    iteration::Int
    accepted::Bool
    reject_reason::Union{Nothing,String}
    evolution::EvolutionLog
    observables::SimpleObservableSummary
    simple_score::ScarFinderCandidateScore
    ctm_score::Union{Nothing,ScarFinderCandidateScore}
end

ScarFinderIteration(
    iteration::Int,
    accepted::Bool,
    reject_reason::Union{Nothing,String},
    evolution::EvolutionLog,
    observables::SimpleObservableSummary,
) = ScarFinderIteration(
    iteration,
    accepted,
    reject_reason,
    evolution,
    observables,
    _candidate_score(iteration, :simple, accepted, reject_reason, evolution, observables),
    nothing,
)

"""
    ScarFinderResult

Result of [`scarfinder!`](@ref). `state` is the same mutably evolved
[`SquareIPEPSState`](@ref) passed to `scarfinder!`; iteration records contain
the local/simple diagnostics used for acceptance plus optional CTM callback
diagnostics. Simple diagnostics are not CTMRG-quality environment
measurements.
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
    return all(
        isfinite,
        (
            obs.density,
            obs.density_even,
            obs.density_odd,
            obs.blockade_violation,
            obs.pxp_energy_density,
            obs.mean_bond_entropy,
            obs.max_bond_entropy,
        ),
    )
end

function _finite_summary(obs::CTMObservableSummary)
    return all(
        isfinite,
        (
            obs.density,
            obs.density_even,
            obs.density_odd,
            obs.blockade_violation,
            obs.pxp_energy_density,
        ),
    )
end

function _score_value(
    obs,
    log::EvolutionLog,
    mean_bond_entropy::Union{Nothing,Float64},
    max_bond_entropy::Union{Nothing,Float64},
)
    entropy_penalty = max_bond_entropy === nothing ? 0.0 : max_bond_entropy
    return abs(obs.pxp_energy_density) +
           obs.blockade_violation +
           entropy_penalty +
           log.max_truncerr
end

function _candidate_score(
    iteration::Int,
    diagnostics::Symbol,
    accepted::Bool,
    reject_reason::Union{Nothing,String},
    log::EvolutionLog,
    obs::SimpleObservableSummary,
)
    return ScarFinderCandidateScore(
        iteration,
        diagnostics,
        accepted,
        reject_reason,
        obs.density,
        obs.density_even,
        obs.density_odd,
        obs.blockade_violation,
        obs.pxp_energy_density,
        obs.mean_bond_entropy,
        obs.max_bond_entropy,
        log.max_truncerr,
        _score_value(obs, log, obs.mean_bond_entropy, obs.max_bond_entropy),
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )
end

function _ctm_score_fields(::Nothing)
    return (nothing, nothing, nothing, nothing, nothing, nothing, nothing)
end

function _ctm_score_fields(diagnostics::CTMRGDiagnostics)
    return (
        diagnostics.chi,
        diagnostics.tol,
        diagnostics.maxiter,
        diagnostics.iterations,
        diagnostics.residual,
        diagnostics.converged,
        diagnostics.accepted,
    )
end

function _candidate_score(
    iteration::Int,
    diagnostics::Symbol,
    accepted::Bool,
    reject_reason::Union{Nothing,String},
    log::EvolutionLog,
    obs::CTMObservableSummary,
)
    ctm_fields = _ctm_score_fields(obs.diagnostics)
    return ScarFinderCandidateScore(
        iteration,
        diagnostics,
        accepted,
        reject_reason,
        obs.density,
        obs.density_even,
        obs.density_odd,
        obs.blockade_violation,
        obs.pxp_energy_density,
        nothing,
        nothing,
        log.max_truncerr,
        _score_value(obs, log, nothing, nothing),
        ctm_fields...,
    )
end

function _candidate_score(
    iteration::Int,
    diagnostics::Symbol,
    accepted::Bool,
    reject_reason::Union{Nothing,String},
    log::EvolutionLog,
    obs,
)
    throw(
        ArgumentError(
            "CTM callback must return SimpleObservableSummary or CTMObservableSummary",
        ),
    )
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

function _nonnegative_int(value::Integer, name::String)
    converted = Int(value)
    converted >= 0 || throw(ArgumentError("$name must be nonnegative"))
    return converted
end

function _should_measure_ctm(iteration::Int, total::Int, ctm_every::Int, ctm_at_end::Bool)
    return (ctm_every > 0 && iteration % ctm_every == 0) ||
           (ctm_at_end && iteration == total)
end

function _iteration_score(iteration::ScarFinderIteration, diagnostics::Symbol)
    if diagnostics === :simple
        return iteration.simple_score
    elseif diagnostics === :ctm
        return iteration.ctm_score
    else
        throw(ArgumentError("diagnostics must be :simple or :ctm"))
    end
end

"""
    rank_scarfinder_candidates(
        result;
        diagnostics = :simple,
        accepted_only = false,
        require_ctm_accepted = false,
        by = :score,
    )

Return candidate scores sorted by a field of [`ScarFinderCandidateScore`](@ref).
Simple ranking is available for every iteration. CTM ranking includes only
iterations where an optional CTM callback supplied diagnostics. Set
`require_ctm_accepted = true` with `diagnostics = :ctm` to exclude CTM records
whose diagnostics are missing or flagged as unaccepted.
"""
function rank_scarfinder_candidates(
    result::ScarFinderResult;
    diagnostics::Symbol = :simple,
    accepted_only::Bool = false,
    require_ctm_accepted::Bool = false,
    by::Symbol = :score,
    rev::Bool = false,
)
    require_ctm_accepted && diagnostics !== :ctm &&
        throw(ArgumentError("require_ctm_accepted is only valid with diagnostics = :ctm"))
    scores = ScarFinderCandidateScore[]
    for iteration in result.iterations
        accepted_only && !iteration.accepted && continue
        score = _iteration_score(iteration, diagnostics)
        score === nothing && continue
        require_ctm_accepted && score.ctm_accepted !== true && continue
        push!(scores, score)
    end
    if !hasfield(ScarFinderCandidateScore, by)
        throw(ArgumentError("unknown ScarFinderCandidateScore field: $by"))
    end
    return sort(scores; by = score -> getfield(score, by), rev)
end

"""
    scarfinder!(psi::SquareIPEPSState, params::ScarFinderParams)::ScarFinderResult

Run the S6-lite ScarFinder projection scaffold in place on `psi`.
For each requested iteration, this calls [`evolve!`](@ref) for
`params.projection_time`, calls [`measure_simple`](@ref), records diagnostics,
and accepts or rejects the iteration by truncation error, simple blockade
violation, simple bond entropy, and finite-diagnostic checks. `psi` is mutated
like [`evolve!`](@ref). CTM diagnostics are recorded only when a caller
supplies `ctm_callback` and schedules it with `ctm_every` or `ctm_at_end`. No
energy correction, new CTMRG algorithm, or direct star-projection logic is
performed here.
"""
function scarfinder!(
    psi::SquareIPEPSState,
    params::ScarFinderParams;
    ctm_callback = nothing,
    ctm_every::Integer = 0,
    ctm_at_end::Bool = false,
    log_path::Union{Nothing,AbstractString} = nothing,
    log_format::Symbol = :csv,
)::ScarFinderResult
    ctm_period = _nonnegative_int(ctm_every, "ctm_every")
    iterations = ScarFinderIteration[]

    for n = 1:params.iterations
        log = evolve!(psi, params.projection_time; params = params.trotter)
        obs = measure_simple(psi)
        accepted, reason = _evaluate_scarfinder_iteration(log, obs, params)
        simple_score = _candidate_score(n, :simple, accepted, reason, log, obs)
        ctm_score = nothing
        if ctm_callback !== nothing &&
           _should_measure_ctm(n, params.iterations, ctm_period, ctm_at_end)
            ctm_obs = ctm_callback(psi, n, simple_score)
            _finite_summary(ctm_obs) || throw(ArgumentError("non-finite CTM diagnostic"))
            ctm_score = _candidate_score(n, :ctm, accepted, reason, log, ctm_obs)
        end
        push!(
            iterations,
            ScarFinderIteration(n, accepted, reason, log, obs, simple_score, ctm_score),
        )

        if !accepted && params.stop_on_reject
            break
        end
    end

    result = ScarFinderResult(
        psi,
        params,
        iterations,
        _count_accepted(iterations),
        _count_rejected(iterations),
    )
    log_path === nothing || write_scarfinder_log(result, log_path; format = log_format)
    return result
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
    ctm_callback = nothing,
    ctm_every::Integer = 0,
    ctm_at_end::Bool = false,
    log_path::Union{Nothing,AbstractString} = nothing,
    log_format::Symbol = :csv,
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
    return scarfinder!(
        psi,
        params;
        ctm_callback,
        ctm_every,
        ctm_at_end,
        log_path,
        log_format,
    )
end

function _csv_value(value::Nothing)
    return ""
end

function _csv_value(value::Bool)
    return string(value)
end

function _csv_value(value::Real)
    return string(value)
end

function _csv_value(value::Symbol)
    return String(value)
end

function _csv_value(value::AbstractString)
    escaped = replace(value, "\"" => "\"\"")
    return any(ch -> ch in escaped, (',', '"', '\n', '\r')) ? "\"$escaped\"" : escaped
end

function _score_rows(result::ScarFinderResult)
    rows = ScarFinderCandidateScore[]
    for iteration in result.iterations
        push!(rows, iteration.simple_score)
        iteration.ctm_score === nothing || push!(rows, iteration.ctm_score)
    end
    return rows
end

function _write_csv_log(io, result::ScarFinderResult)
    header = (
        "iteration",
        "accepted",
        "reject_reason",
        "diagnostics",
        "density",
        "density_even",
        "density_odd",
        "blockade_violation",
        "pxp_energy_density",
        "mean_bond_entropy",
        "max_bond_entropy",
        "max_truncerr",
        "score",
        "ctm_chi",
        "ctm_tol",
        "ctm_maxiter",
        "ctm_iterations",
        "ctm_residual",
        "ctm_converged",
        "ctm_accepted",
    )
    println(io, join(header, ","))
    for score in _score_rows(result)
        row = (
            score.iteration,
            score.accepted,
            score.reject_reason,
            score.diagnostics,
            score.density,
            score.density_even,
            score.density_odd,
            score.blockade_violation,
            score.pxp_energy_density,
            score.mean_bond_entropy,
            score.max_bond_entropy,
            score.max_truncerr,
            score.score,
            score.ctm_chi,
            score.ctm_tol,
            score.ctm_maxiter,
            score.ctm_iterations,
            score.ctm_residual,
            score.ctm_converged,
            score.ctm_accepted,
        )
        println(io, join(_csv_value.(row), ","))
    end
end

function _json_escape(value::AbstractString)
    escaped = replace(value, "\\" => "\\\\")
    escaped = replace(escaped, "\"" => "\\\"")
    escaped = replace(escaped, "\n" => "\\n")
    escaped = replace(escaped, "\r" => "\\r")
    return escaped
end

_json_value(value::Nothing) = "null"
_json_value(value::Bool) = value ? "true" : "false"
_json_value(value::Real) = isfinite(value) ? string(value) : "\"$(value)\""
_json_value(value::Symbol) = _json_value(String(value))
_json_value(value::AbstractString) = "\"$(_json_escape(value))\""

function _score_json(score::ScarFinderCandidateScore)
    fields = (
        :iteration,
        :accepted,
        :reject_reason,
        :diagnostics,
        :density,
        :density_even,
        :density_odd,
        :blockade_violation,
        :pxp_energy_density,
        :mean_bond_entropy,
        :max_bond_entropy,
        :max_truncerr,
        :score,
        :ctm_chi,
        :ctm_tol,
        :ctm_maxiter,
        :ctm_iterations,
        :ctm_residual,
        :ctm_converged,
        :ctm_accepted,
    )
    pairs = ["\"$(field)\":$(_json_value(getfield(score, field)))" for field in fields]
    return "{" * join(pairs, ",") * "}"
end

function _write_json_log(io, result::ScarFinderResult)
    rows = _score_json.(_score_rows(result))
    println(io, "{\"iterations\":[" * join(rows, ",") * "]}")
end

"""
    write_scarfinder_log(result, path; format = :csv)

Write recorded ScarFinder candidate diagnostics to CSV or JSON. This function
serializes diagnostics already stored in `result`; it does not run CTMRG or
refresh any environment.
"""
function write_scarfinder_log(
    result::ScarFinderResult,
    path::AbstractString;
    format::Symbol = :csv,
)
    open(path, "w") do io
        if format === :csv
            _write_csv_log(io, result)
        elseif format === :json
            _write_json_log(io, result)
        else
            throw(ArgumentError("log format must be :csv or :json"))
        end
    end
    return path
end

end
