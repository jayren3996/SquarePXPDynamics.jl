module ScarFinder

using ..SquareIPEPS: SquareIPEPSState, copy_state, state_version, log_norm
using ..IPEPSEvolution:
    TrotterParams, EvolutionLog, legacy_trotter_params, legacy_trotter_protocol, evolve!
using ..StarModels: AbstractModelProtocol
using ..Observables: SimpleObservableSummary, measure_simple
using ..PXPValidation: TrustedCTMMeasurement, measure_ctm_trusted
using ..CTMTrust: CTMTrustPolicy
using ..PEPSKitMeasurements: CTMObservableSummary, CTMRGDiagnostics
using ..PEPSKitMeasurements: PEPSKitCTMRGParams, measure_ctm

export ScarFinderParams, ScarFinderCandidateScore, ScarFinderIteration, ScarFinderResult
export MeasurementBackend, SimpleBackend, TrustedCTMBackend, measure_scarfinder
export ScarFinderObjective, RevivalObjective, TargetEnergyObjective
export LowVarianceObjective, CompositeObjective
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
        target_energy = nothing,
        correction_time = 0,
        correction_attempts = 0,
    )

Parameters for the S6-lite ScarFinder orchestration loop. Each iteration
evolves the supplied iPEPS state with [`evolve!`](@ref), records
[`measure_simple`](@ref) diagnostics, and accepts or rejects the iteration
using simple-update/local observables only. If `target_energy` is supplied,
short imaginary-time correction attempts may be applied using simple/local PXP
energy as the diagnostic objective. This does not run CTMRG or make
physics-quality energy claims. Candidate scores are computed from recorded
diagnostics by the orchestration layer.
"""
struct ScarFinderParams
    projection_time::Float64
    trotter::TrotterParams
    protocol::Union{Nothing,AbstractModelProtocol}
    iterations::Int
    max_truncerr::Float64
    max_blockade_violation::Float64
    max_bond_entropy::Float64
    stop_on_reject::Bool
    target_energy::Union{Nothing,Float64}
    correction_time::Float64
    correction_attempts::Int

    function ScarFinderParams(
        projection_time::Real,
        trotter,
        iterations::Integer,
        max_truncerr::Real,
        max_blockade_violation::Real,
        max_bond_entropy::Real,
        stop_on_reject::Bool,
        ;
        target_energy::Union{Nothing,Real} = nothing,
        correction_time::Real = 0,
        correction_attempts::Integer = 0,
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
        energy_target = _optional_finite_float(target_energy, "target_energy")
        correction_step = _finite_nonnegative_float(correction_time, "correction_time")
        correction_count = _nonnegative_int(correction_attempts, "correction_attempts")

        return new(
            time,
            legacy_trotter_params(trotter),
            legacy_trotter_protocol(trotter),
            niterations,
            truncerr_limit,
            blockade_limit,
            entropy_limit,
            stop_on_reject,
            energy_target,
            correction_step,
            correction_count,
        )
    end
end

"""
    MeasurementBackend

Abstract ScarFinder measurement backend interface.
"""
abstract type MeasurementBackend end

"""Development backend using local/simple-update observables only."""
struct SimpleBackend <: MeasurementBackend end

"""Trusted finite-chi CTM backend for production ScarFinder diagnostics."""
struct TrustedCTMBackend <: MeasurementBackend
    params::Tuple
    policy::CTMTrustPolicy
    measure::Any

    function TrustedCTMBackend(
        params,
        policy::CTMTrustPolicy = CTMTrustPolicy();
        measure = measure_ctm,
    )
        collected = Tuple(params)
        isempty(collected) &&
            throw(ArgumentError("TrustedCTMBackend requires at least one CTMRG parameter point"))
        all(p -> p isa PEPSKitCTMRGParams, collected) ||
            throw(ArgumentError("TrustedCTMBackend params must be PEPSKitCTMRGParams"))
        return new(collected, policy, measure)
    end
end

"""Measure `psi` with a ScarFinder measurement backend."""
measure_scarfinder(psi::SquareIPEPSState, ::SimpleBackend) = measure_simple(psi)

function measure_scarfinder(psi::SquareIPEPSState, backend::TrustedCTMBackend)
    return measure_ctm_trusted(
        psi;
        params = backend.params,
        policy = backend.policy,
        measure = backend.measure,
    )
end

"""
    ScarFinderObjective

Abstract type for ScarFinder candidate ranking objectives.
"""
abstract type ScarFinderObjective end

"""
    RevivalObjective(observable = :sublattice_imbalance, weight = 1.0)

Reward revival strength in ScarFinder ranking. Supported observables are
`:sublattice_imbalance` and `:density`.
"""
struct RevivalObjective <: ScarFinderObjective
    observable::Symbol
    weight::Float64

    function RevivalObjective(observable::Symbol = :sublattice_imbalance, weight::Real = 1.0)
        observable in (:sublattice_imbalance, :density) ||
            throw(ArgumentError("revival observable must be :sublattice_imbalance or :density"))
        w = Float64(weight)
        isfinite(w) && w >= 0 ||
            throw(ArgumentError("revival weight must be finite and nonnegative"))
        return new(observable, w)
    end
end

"""
    TargetEnergyObjective(target, weight = 1.0)

Penalize distance from a target simple/local PXP energy density in ScarFinder
ranking.
"""
struct TargetEnergyObjective <: ScarFinderObjective
    target::Float64
    weight::Float64

    function TargetEnergyObjective(target::Real, weight::Real = 1.0)
        t = Float64(target)
        isfinite(t) || throw(ArgumentError("target energy must be finite"))
        w = Float64(weight)
        isfinite(w) && w >= 0 ||
            throw(ArgumentError("target energy weight must be finite and nonnegative"))
        return new(t, w)
    end
end

"""
    LowVarianceObjective(weight = 1.0)

Placeholder objective component for future energy-variance proxy penalties.
"""
struct LowVarianceObjective <: ScarFinderObjective
    weight::Float64

    function LowVarianceObjective(weight::Real = 1.0)
        w = Float64(weight)
        isfinite(w) && w >= 0 ||
            throw(ArgumentError("low variance weight must be finite and nonnegative"))
        return new(w)
    end
end

"""
    CompositeObjective(; revival = RevivalObjective(), target_energy = nothing,
                         low_variance = nothing, blockade_weight = 100.0,
                         truncation_weight = 10.0, finite_chi_weight = 10.0,
                         entropy_weight = 1.0)

Weighted ScarFinder ranking objective combining revival rewards with blockade,
truncation, finite-`chi`, entropy, and optional target-energy penalties.
"""
struct CompositeObjective <: ScarFinderObjective
    revival::Union{Nothing,RevivalObjective}
    target_energy::Union{Nothing,TargetEnergyObjective}
    low_variance::Union{Nothing,LowVarianceObjective}
    blockade_weight::Float64
    truncation_weight::Float64
    finite_chi_weight::Float64
    entropy_weight::Float64

    function CompositeObjective(;
        revival::Union{Nothing,RevivalObjective} = RevivalObjective(),
        target_energy::Union{Nothing,TargetEnergyObjective} = nothing,
        low_variance::Union{Nothing,LowVarianceObjective} = nothing,
        blockade_weight::Real = 100.0,
        truncation_weight::Real = 10.0,
        finite_chi_weight::Real = 10.0,
        entropy_weight::Real = 1.0,
    )
        weights = Float64.((blockade_weight, truncation_weight, finite_chi_weight, entropy_weight))
        all(w -> isfinite(w) && w >= 0, weights) ||
            throw(ArgumentError("objective weights must be finite and nonnegative"))
        return new(revival, target_energy, low_variance, weights...)
    end
end

"""
    ScarFinderCandidateScore

One ranked ScarFinder candidate diagnostic record. `diagnostics` is `:simple`
for the default local/simple measurements or `:ctm` for records supplied by an
optional CTM measurement callback. The scalar `score` is an operational sorting
key built from existing diagnostics only; it is not an energy correction or a
physics claim. `objective_parameters` records the deterministic objective
settings needed to audit that scalar score. The `log_norm_*` fields mirror the
[`EvolutionLog`](@ref) normalization ledger for the candidate-producing
evolution call.
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
    objective_name::String
    objective_parameters::String
    revival_strength::Union{Nothing,Float64}
    finite_chi_drift::Union{Nothing,Float64}
    energy_variance_proxy::Union{Nothing,Float64}
    log_norm_before::Float64
    log_norm_after::Float64
    log_norm_delta::Float64
    ctm_chi::Union{Nothing,Int}
    ctm_tol::Union{Nothing,Float64}
    ctm_maxiter::Union{Nothing,Int}
    ctm_iterations::Union{Nothing,Int}
    ctm_residual::Union{Nothing,Float64}
    ctm_converged::Union{Nothing,Bool}
    ctm_accepted::Union{Nothing,Bool}
    correction_accepted::Union{Nothing,Bool}
    correction_energy_before::Union{Nothing,Float64}
    correction_energy_after::Union{Nothing,Float64}
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
    correction_accepted::Union{Nothing,Bool}
    correction_energy_before::Union{Nothing,Float64}
    correction_energy_after::Union{Nothing,Float64}
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
    nothing,
    nothing,
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

function _finite_nonnegative_float(value::Real, name::String)
    converted = Float64(value)
    isfinite(converted) || throw(ArgumentError("$name must be finite"))
    converted >= 0 || throw(ArgumentError("$name must be nonnegative"))
    return converted
end

_optional_finite_float(::Nothing, name::String) = nothing

function _optional_finite_float(value::Real, name::String)
    converted = Float64(value)
    isfinite(converted) || throw(ArgumentError("$name must be finite"))
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

_scarfinder_observable(obs) = obs
_scarfinder_observable(obs::TrustedCTMMeasurement) = obs.measurement
_scarfinder_trusted_ctm(obs) = nothing
_scarfinder_trusted_ctm(obs::TrustedCTMMeasurement) = obs

function _imbalance(obs)
    return obs.density_even - obs.density_odd
end

function _revival_strength(obs, objective::RevivalObjective)
    objective.observable === :sublattice_imbalance && return abs(_imbalance(obs))
    objective.observable === :density && return abs(obs.density)
    throw(ArgumentError("unsupported revival observable"))
end

_finite_chi_drift(::Nothing) = nothing

function _finite_chi_drift(trusted::TrustedCTMMeasurement)
    vals = (
        trusted.trust.finite_chi_density_delta,
        trusted.trust.finite_chi_blockade_delta,
        trusted.trust.finite_chi_energy_delta,
    )
    present = Float64[]
    for v in vals
        v === nothing || push!(present, abs(v))
    end
    return isempty(present) ? nothing : maximum(present)
end

function _composite_objective(objective::CompositeObjective)
    return objective
end

function _composite_objective(objective::RevivalObjective)
    return CompositeObjective(; revival = objective)
end

function _composite_objective(objective::TargetEnergyObjective)
    return CompositeObjective(; revival = nothing, target_energy = objective)
end

function _composite_objective(objective::LowVarianceObjective)
    return CompositeObjective(; revival = nothing, low_variance = objective)
end

function _objective_parameter_value(value::Nothing)
    return "nothing"
end

function _objective_parameter_value(value::Symbol)
    return String(value)
end

function _objective_parameter_value(value::Real)
    return string(Float64(value))
end

function _objective_parameters(objective::CompositeObjective)
    revival_observable =
        objective.revival === nothing ? nothing : objective.revival.observable
    revival_weight = objective.revival === nothing ? nothing : objective.revival.weight
    target_energy =
        objective.target_energy === nothing ? nothing : objective.target_energy.target
    target_energy_weight =
        objective.target_energy === nothing ? nothing : objective.target_energy.weight
    low_variance_weight =
        objective.low_variance === nothing ? nothing : objective.low_variance.weight
    fields = (
        :revival_observable => revival_observable,
        :revival_weight => revival_weight,
        :target_energy => target_energy,
        :target_energy_weight => target_energy_weight,
        :low_variance_weight => low_variance_weight,
        :blockade_weight => objective.blockade_weight,
        :truncation_weight => objective.truncation_weight,
        :finite_chi_weight => objective.finite_chi_weight,
        :entropy_weight => objective.entropy_weight,
    )
    return join(
        ("$(key)=$(_objective_parameter_value(value))" for (key, value) in fields),
        ";",
    )
end

function _score_value(
    obs,
    log::EvolutionLog,
    mean_bond_entropy::Union{Nothing,Float64},
    max_bond_entropy::Union{Nothing,Float64},
    objective::CompositeObjective;
    trusted_ctm = nothing,
)
    revival = objective.revival === nothing ? nothing : _revival_strength(obs, objective.revival)
    finite_chi = _finite_chi_drift(trusted_ctm)
    entropy_penalty = max_bond_entropy === nothing ? 0.0 : max_bond_entropy
    score = obs.blockade_violation * objective.blockade_weight +
            log.max_truncerr * objective.truncation_weight +
            entropy_penalty * objective.entropy_weight
    finite_chi === nothing || (score += finite_chi * objective.finite_chi_weight)
    revival === nothing || (score -= objective.revival.weight * revival)
    objective.target_energy === nothing ||
        (score += objective.target_energy.weight * abs(obs.pxp_energy_density - objective.target_energy.target))
    return score, revival, finite_chi
end

function _candidate_score(
    iteration::Int,
    diagnostics::Symbol,
    accepted::Bool,
    reject_reason::Union{Nothing,String},
    log::EvolutionLog,
    obs::SimpleObservableSummary,
    correction_accepted::Union{Nothing,Bool} = nothing,
    correction_energy_before::Union{Nothing,Float64} = nothing,
    correction_energy_after::Union{Nothing,Float64} = nothing,
    objective::ScarFinderObjective = CompositeObjective(),
    trusted_ctm = nothing,
)
    composite = _composite_objective(objective)
    score, revival, finite_chi = _score_value(
        obs,
        log,
        obs.mean_bond_entropy,
        obs.max_bond_entropy,
        composite;
        trusted_ctm,
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
        score,
        String(nameof(typeof(objective))),
        _objective_parameters(composite),
        revival,
        finite_chi,
        nothing,
        log.log_norm_before,
        log.log_norm_after,
        log.log_norm_delta,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        correction_accepted,
        correction_energy_before,
        correction_energy_after,
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
    correction_accepted::Union{Nothing,Bool} = nothing,
    correction_energy_before::Union{Nothing,Float64} = nothing,
    correction_energy_after::Union{Nothing,Float64} = nothing,
    objective::ScarFinderObjective = CompositeObjective(),
    trusted_ctm = nothing,
)
    ctm_fields = _ctm_score_fields(obs.diagnostics)
    composite = _composite_objective(objective)
    score, revival, finite_chi =
        _score_value(obs, log, nothing, nothing, composite; trusted_ctm)
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
        score,
        String(nameof(typeof(objective))),
        _objective_parameters(composite),
        revival,
        finite_chi,
        nothing,
        log.log_norm_before,
        log.log_norm_after,
        log.log_norm_delta,
        ctm_fields...,
        correction_accepted,
        correction_energy_before,
        correction_energy_after,
    )
end

function _candidate_score(
    iteration::Int,
    diagnostics::Symbol,
    accepted::Bool,
    reject_reason::Union{Nothing,String},
    log::EvolutionLog,
    obs,
    correction_accepted::Union{Nothing,Bool} = nothing,
    correction_energy_before::Union{Nothing,Float64} = nothing,
    correction_energy_after::Union{Nothing,Float64} = nothing,
    objective::ScarFinderObjective = CompositeObjective(),
    trusted_ctm = nothing,
)
    throw(
        ArgumentError(
            "CTM callback must return SimpleObservableSummary or CTMObservableSummary",
        ),
    )
end

function _replace_state!(target::SquareIPEPSState, source::SquareIPEPSState)
    empty!(target.tensors)
    for (c, tensor) in source.tensors
        target.tensors[c] = copy(tensor)
    end
    empty!(target.physical_indices)
    for (c, index) in source.physical_indices
        target.physical_indices[c] = index
    end
    empty!(target.link_indices)
    for (key, index) in source.link_indices
        target.link_indices[key] = index
    end
    empty!(target.link_weights)
    for (key, values) in source.link_weights
        target.link_weights[key] = copy(values)
    end
    target.mutation_version[] = state_version(source) + 1
    target.log_norm_value[] = log_norm(source)
    return target
end

function _correction_params(params::ScarFinderParams)
    return TrotterParams(
        params.trotter.dt,
        params.trotter.order,
        :imaginary,
        params.trotter.maxdim,
        params.trotter.cutoff,
        params.trotter.split_order;
        schedule = params.trotter.schedule,
    )
end

function _maybe_apply_energy_correction!(
    psi::SquareIPEPSState,
    obs::SimpleObservableSummary,
    params::ScarFinderParams,
)
    target = params.target_energy
    target === nothing && return (nothing, nothing, nothing, obs)
    before_energy = obs.pxp_energy_density
    best_objective = abs(before_energy - target)
    best_state = nothing
    best_obs = obs
    accepted = false

    if params.correction_time > 0 && params.correction_attempts > 0
        correction_params = _correction_params(params)
        trial_base = copy_state(psi)
        for _ = 1:params.correction_attempts
            trial = copy_state(trial_base)
            evolve!(
                trial,
                params.correction_time;
                params = correction_params,
                protocol = params.protocol,
            )
            trial_obs = measure_simple(trial)
            trial_objective = abs(trial_obs.pxp_energy_density - target)
            if trial_objective < best_objective
                best_objective = trial_objective
                best_state = trial
                best_obs = trial_obs
                trial_base = trial
                accepted = true
            end
        end
    end

    if accepted
        _replace_state!(psi, best_state)
    end
    return (accepted, before_energy, best_obs.pxp_energy_density, best_obs)
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
whose diagnostics are missing or flagged as unaccepted. Nullable sort fields,
such as CTM residuals, are sorted with missing values last.
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
    present = filter(score -> getfield(score, by) !== nothing, scores)
    missing = filter(score -> getfield(score, by) === nothing, scores)
    return vcat(sort(present; by = score -> getfield(score, by), rev), missing)
end

"""
    scarfinder!(psi::SquareIPEPSState, params::ScarFinderParams)::ScarFinderResult

Run the S6-lite ScarFinder projection scaffold in place on `psi`.
For each requested iteration, this calls [`evolve!`](@ref) for
`params.projection_time`, calls [`measure_simple`](@ref), records diagnostics,
and accepts or rejects the iteration by truncation error, simple blockade
violation, simple bond entropy, and finite-diagnostic checks. `psi` is mutated
like [`evolve!`](@ref). If `params.target_energy` is supplied, the loop may
apply guarded simple/local imaginary-time correction attempts and records their
diagnostic outcome. CTM diagnostics are recorded only when a caller supplies
`ctm_callback` and schedules it with `ctm_every` or `ctm_at_end`. No new CTMRG
algorithm or direct star-projection logic is performed here.
"""
function scarfinder!(
    psi::SquareIPEPSState,
    params::ScarFinderParams;
    objective::ScarFinderObjective = CompositeObjective(),
    ctm_callback = nothing,
    ctm_every::Integer = 0,
    ctm_at_end::Bool = false,
    log_path::Union{Nothing,AbstractString} = nothing,
    log_format::Symbol = :csv,
)::ScarFinderResult
    ctm_period = _nonnegative_int(ctm_every, "ctm_every")
    iterations = ScarFinderIteration[]

    for n = 1:params.iterations
        log = evolve!(
            psi,
            params.projection_time;
            params = params.trotter,
            protocol = params.protocol,
        )
        obs_before_correction = measure_simple(psi)
        correction_accepted, correction_energy_before, correction_energy_after, obs =
            _maybe_apply_energy_correction!(psi, obs_before_correction, params)
        accepted, reason = _evaluate_scarfinder_iteration(log, obs, params)
        simple_score = _candidate_score(
            n,
            :simple,
            accepted,
            reason,
            log,
            obs,
            correction_accepted,
            correction_energy_before,
            correction_energy_after,
            objective,
        )
        ctm_score = nothing
        if ctm_callback !== nothing &&
           _should_measure_ctm(n, params.iterations, ctm_period, ctm_at_end)
            raw_ctm = ctm_callback(psi, n, simple_score)
            ctm_obs = _scarfinder_observable(raw_ctm)
            _finite_summary(ctm_obs) || throw(ArgumentError("non-finite CTM diagnostic"))
            ctm_score = _candidate_score(
                n,
                :ctm,
                accepted,
                reason,
                log,
                ctm_obs,
                correction_accepted,
                correction_energy_before,
                correction_energy_after,
                objective,
                _scarfinder_trusted_ctm(raw_ctm),
            )
        end
        push!(
            iterations,
            ScarFinderIteration(
                n,
                accepted,
                reason,
                log,
                obs,
                simple_score,
                ctm_score,
                correction_accepted,
                correction_energy_before,
                correction_energy_after,
            ),
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
        target_energy = nothing,
        correction_time = 0,
        correction_attempts = 0,
        objective = CompositeObjective(),
    )::ScarFinderResult

Convenience keyword constructor for [`ScarFinderParams`](@ref), then delegates
to [`scarfinder!(psi, params)`](@ref). The loop records simple/local
diagnostics from [`measure_simple`](@ref), with optional guarded energy
correction fields when `target_energy` is supplied.
"""
function scarfinder!(
    psi::SquareIPEPSState;
    projection_time::Real,
    trotter,
    iterations::Integer,
    max_truncerr::Real = Inf,
    max_blockade_violation::Real = Inf,
    max_bond_entropy::Real = Inf,
    stop_on_reject::Bool = false,
    target_energy::Union{Nothing,Real} = nothing,
    correction_time::Real = 0,
    correction_attempts::Integer = 0,
    ctm_callback = nothing,
    ctm_every::Integer = 0,
    ctm_at_end::Bool = false,
    log_path::Union{Nothing,AbstractString} = nothing,
    log_format::Symbol = :csv,
    objective::ScarFinderObjective = CompositeObjective(),
)::ScarFinderResult
    params = ScarFinderParams(
        projection_time,
        trotter,
        iterations,
        max_truncerr,
        max_blockade_violation,
        max_bond_entropy,
        stop_on_reject;
        target_energy,
        correction_time,
        correction_attempts,
    )
    return scarfinder!(
        psi,
        params;
        objective,
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
        "objective_name",
        "objective_parameters",
        "revival_strength",
        "finite_chi_drift",
        "energy_variance_proxy",
        "log_norm_before",
        "log_norm_after",
        "log_norm_delta",
        "ctm_chi",
        "ctm_tol",
        "ctm_maxiter",
        "ctm_iterations",
        "ctm_residual",
        "ctm_converged",
        "ctm_accepted",
        "correction_accepted",
        "correction_energy_before",
        "correction_energy_after",
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
            score.objective_name,
            score.objective_parameters,
            score.revival_strength,
            score.finite_chi_drift,
            score.energy_variance_proxy,
            score.log_norm_before,
            score.log_norm_after,
            score.log_norm_delta,
            score.ctm_chi,
            score.ctm_tol,
            score.ctm_maxiter,
            score.ctm_iterations,
            score.ctm_residual,
            score.ctm_converged,
            score.ctm_accepted,
            score.correction_accepted,
            score.correction_energy_before,
            score.correction_energy_after,
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
        :objective_name,
        :objective_parameters,
        :revival_strength,
        :finite_chi_drift,
        :energy_variance_proxy,
        :log_norm_before,
        :log_norm_after,
        :log_norm_delta,
        :ctm_chi,
        :ctm_tol,
        :ctm_maxiter,
        :ctm_iterations,
        :ctm_residual,
        :ctm_converged,
        :ctm_accepted,
        :correction_accepted,
        :correction_energy_before,
        :correction_energy_after,
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
