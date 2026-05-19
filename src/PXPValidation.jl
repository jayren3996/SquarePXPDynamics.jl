module PXPValidation

using JSON3

using ..SquareIPEPS: SquareIPEPSState, copy_state, log_norm, product_square_ipeps
using ..SquareUnitCells: PeriodicSquareUnitCell
using ..Observables: SimpleObservableSummary, measure_simple
using ..FiniteIPEPSObservables: exact_density_finite, exact_all_down_return_probability_finite
using ..PEPSKitMeasurements:
    CTMObservableSummary,
    CTMValidationPoint,
    PEPSKitCTMRGParams,
    default_ctmrg_params,
    measure_ctm,
    validate_ctm_sweep
using ..CTMTrust: CTMTrustAssessment, CTMTrustPolicy, assess_ctm_trust
using ..IPEPSEvolution: EvolutionLog, TrotterParams, evolve!, reverse_evolve!
using ..FinitePXPEEDBenchmark:
    PXPEEDBenchmarkConfig,
    PXPEEDBenchmarkResult,
    PXPEEDSample,
    run_pxp_ed_benchmark

export TrustedCTMMeasurement, measure_ctm_trusted
export PXPValidationConfig, PXPValidationMetadata, PXPIPEPSSample
export PXPEDComparisonSample, PXPValidationReport, validate_pxp_ed_ipeps
export write_pxp_validation_json
export PXPConvergenceConfig, PXPConvergenceReport, validate_pxp_convergence
export write_pxp_convergence_json
export PXPReversibilityReport, validate_pxp_reversibility
export PXPAuditConfig, PXPAuditSummary, PXPAuditRun, PXPAuditReport
export run_pxp_audit_campaign, write_pxp_audit_json, write_pxp_audit_csv
export PXPLargerDBenchmarkConfig, PXPLargerDBenchmarkSummary
export PXPLargerDBenchmarkRun, PXPLargerDBenchmarkReport
export run_pxp_larger_d_benchmark
export write_pxp_larger_d_benchmark_json, write_pxp_larger_d_benchmark_csv

"""
    TrustedCTMMeasurement(measurement, points, trust, policy = CTMTrustPolicy())

Finite-`chi` CTMRG measurement bundle for one iPEPS state. `points` stores the
full validation sweep, `measurement` is the last sweep point's CTM observable
summary, `trust` is the finite-`chi` assessment returned by
[`assess_ctm_trust`](@ref), and `policy` records the trust thresholds used to
produce that assessment.
"""
struct TrustedCTMMeasurement
    measurement::CTMObservableSummary
    points::Vector{CTMValidationPoint}
    trust::CTMTrustAssessment
    policy::CTMTrustPolicy

    function TrustedCTMMeasurement(
        measurement::CTMObservableSummary,
        points::Vector{CTMValidationPoint},
        trust::CTMTrustAssessment,
        policy::CTMTrustPolicy = CTMTrustPolicy(),
    )
        isempty(points) &&
            throw(ArgumentError("trusted CTM measurement requires at least one sweep point"))
        points[end].measurement == measurement ||
            throw(ArgumentError("measurement must match the final CTM validation point"))
        return new(measurement, points, trust, policy)
    end
end

"""
    measure_ctm_trusted(psi; params, policy = CTMTrustPolicy(),
                        reference = measure_simple(psi), measure = measure_ctm)

Run a CTMRG validation sweep for `psi`, assess finite-`chi` trust, and return a
[`TrustedCTMMeasurement`](@ref). The `measure` keyword exists so tests and
benchmark scripts can supply deterministic synthetic CTM summaries without
running PEPSKit CTMRG.
"""
function measure_ctm_trusted(
    psi::SquareIPEPSState;
    params,
    policy::CTMTrustPolicy = CTMTrustPolicy(),
    reference = measure_simple(psi),
    measure = measure_ctm,
)::TrustedCTMMeasurement
    points = validate_ctm_sweep(psi; params, reference, measure)
    assessment = assess_ctm_trust(points; policy)
    return TrustedCTMMeasurement(points[end].measurement, points, assessment, policy)
end

function _finite_nonnegative(value::Real, label::String)
    converted = Float64(value)
    isfinite(converted) && converted >= 0 ||
        throw(ArgumentError("$label must be finite and nonnegative"))
    return converted
end

function _finite_positive(value::Real, label::String)
    converted = Float64(value)
    isfinite(converted) && converted > 0 ||
        throw(ArgumentError("$label must be finite and positive"))
    return converted
end

function _nonnegative_int(value::Integer, label::String)
    converted = Int(value)
    converted >= 0 || throw(ArgumentError("$label must be nonnegative"))
    return converted
end

function _positive_int(value::Integer, label::String)
    converted = Int(value)
    converted >= 1 || throw(ArgumentError("$label must be at least 1"))
    return converted
end

include("PXPValidationTypes.jl")


function _validation_ed_config(config::PXPValidationConfig)
    return PXPEEDBenchmarkConfig(
        config.n;
        total_time = config.total_time,
        dt = config.dt,
        measure_every = config.measure_every,
        initial_state = config.initial_state,
        point_group = config.point_group,
        use_sparse = config.use_sparse,
        tol = config.ed_tol,
        m_init = config.ed_m_init,
        m_max = config.ed_m_max,
        extend_step = config.ed_extend_step,
    )
end

function _copy_config(
    base::PXPValidationConfig;
    dt = base.dt,
    maxdim = base.maxdim,
    cutoff = base.cutoff,
)
    return PXPValidationConfig(
        base.n;
        total_time = base.total_time,
        dt,
        measure_every = max(1, round(Int, base.measure_every * base.dt / dt)),
        initial_state = base.initial_state,
        point_group = base.point_group,
        use_sparse = base.use_sparse,
        ed_tol = base.ed_tol,
        ed_m_init = base.ed_m_init,
        ed_m_max = base.ed_m_max,
        ed_extend_step = base.ed_extend_step,
        order = base.order,
        maxdim,
        cutoff,
        schedule = base.schedule,
        exact_finite_observables = base.exact_finite_observables,
        exact_finite_max_sites = base.exact_finite_max_sites,
    )
end

function _validation_trotter(config::PXPValidationConfig)
    return TrotterParams(
        config.dt,
        config.order,
        :real,
        config.maxdim,
        config.cutoff;
        schedule = config.schedule,
    )
end

"""
    validate_pxp_reversibility(psi, total_time; params, protocol = nothing)

Measure simple PXP observables before evolution, after forward real-time
evolution, and after applying [`reverse_evolve!`](@ref) for the same duration.
The supplied state is copied before mutation, so the returned
[`PXPReversibilityReport`](@ref) is a validation artifact and does not change
the caller's `psi`.
"""
function validate_pxp_reversibility(
    psi::SquareIPEPSState,
    total_time::Real;
    params::TrotterParams,
    protocol = nothing,
)::PXPReversibilityReport
    work = copy_state(psi)
    before = measure_simple(work)
    forward_log = evolve!(work, total_time; params, protocol)
    after_forward = measure_simple(work)
    reverse_log = reverse_evolve!(work, total_time; params, protocol)
    after_reverse = measure_simple(work)

    return PXPReversibilityReport(
        before,
        after_forward,
        after_reverse,
        forward_log,
        reverse_log,
        abs(after_reverse.density - before.density),
        abs(after_reverse.blockade_violation - before.blockade_violation),
        abs(after_reverse.pxp_energy_density - before.pxp_energy_density),
    )
end

function _validation_initial_state(config::PXPValidationConfig)
    cell = PeriodicSquareUnitCell(config.n, config.n)
    state = config.initial_state === :all_down ? :down : config.initial_state
    return product_square_ipeps(cell; state, maxdim = config.maxdim)
end

function _git_commit()
    package_root = abspath(joinpath(@__DIR__, ".."))
    git_dir = joinpath(package_root, ".git")
    for git_command in (
        `git -C $package_root rev-parse HEAD`,
        `git --git-dir $git_dir --work-tree $package_root rev-parse HEAD`,
    )
        try
            command = pipeline(git_command; stderr = devnull)
            commit = chomp(read(command, String))
            isempty(commit) || return commit
        catch
        end
    end
    return nothing
end

function _validation_metadata()
    project = Base.active_project()
    return PXPValidationMetadata(_git_commit(), string(VERSION), project)
end

function _maybe_trusted_ctm(
    psi::SquareIPEPSState,
    ctm_params,
    trust_policy::CTMTrustPolicy,
    ctm_measure,
)
    ctm_params === nothing && return nothing
    return measure_ctm_trusted(
        psi;
        params = ctm_params,
        policy = trust_policy,
        measure = ctm_measure,
    )
end

function _comparison(ed::PXPEEDSample, sample::PXPIPEPSSample)
    ctm = sample.ctm
    ctm_density = ctm === nothing ? nothing : ctm.measurement.density
    exact_density = sample.exact_finite_density
    return PXPEDComparisonSample(
        ed.step,
        ed.time,
        ed.return_probability,
        ed.excitation_density,
        sample.simple.density,
        ctm_density,
        exact_density,
        sample.simple.density - ed.excitation_density,
        ctm_density === nothing ? nothing : ctm_density - ed.excitation_density,
        exact_density === nothing ? nothing : exact_density - ed.excitation_density,
        sample.simple.blockade_violation,
        ctm === nothing ? nothing : ctm.measurement.blockade_violation,
        ctm === nothing ? nothing : ctm.trust.trusted,
        ctm === nothing ? nothing : ctm.trust.reason,
    )
end

function _validate_pxp_ipeps_against_ed(
    config::PXPValidationConfig,
    ed_result::PXPEEDBenchmarkResult;
    ctm_params = nothing,
    trust_policy::CTMTrustPolicy = CTMTrustPolicy(),
    ctm_measure = measure_ctm,
)::PXPValidationReport
    psi = _validation_initial_state(config)
    trotter = _validation_trotter(config)
    samples = PXPIPEPSSample[]
    last_time = 0.0

    for ed_sample in ed_result.samples
        interval = ed_sample.time - last_time
        evolution = if iszero(interval)
            nothing
        else
            evolve!(psi, interval; params = trotter)
        end
        last_time = ed_sample.time

        simple = measure_simple(psi)
        exact_finite_density = config.exact_finite_observables ?
            exact_density_finite(psi; max_sites = config.exact_finite_max_sites) : nothing
        ctm = _maybe_trusted_ctm(psi, ctm_params, trust_policy, ctm_measure)
        push!(
            samples,
            PXPIPEPSSample(
                ed_sample.step,
                ed_sample.time,
                simple,
                evolution,
                ctm,
                log_norm(psi),
                exact_finite_density,
            ),
        )
    end

    comparisons = [
        _comparison(ed_sample, ipeps_sample) for
        (ed_sample, ipeps_sample) in zip(ed_result.samples, samples)
    ]
    return PXPValidationReport(config, ed_result, samples, comparisons, _validation_metadata())
end

"""
    validate_pxp_ed_ipeps(config; ctm_params = nothing,
                          trust_policy = CTMTrustPolicy(),
                          ctm_measure = measure_ctm)

Run a finite periodic PXP ED trajectory and a matched all-down iPEPS trajectory
on an `n x n` unit cell. iPEPS samples are measured at the same times as ED.
When `ctm_params` is supplied, every iPEPS sample also receives a trusted CTM
measurement bundle from [`measure_ctm_trusted`](@ref).

# Example

```julia
using SquarePXPDynamics
config = PXPValidationConfig(3; total_time = 0.02, dt = 0.01, measure_every = 1)
report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
# report.comparisons[i] holds per-time density / energy / return-probability
# residuals between ED and iPEPS.
```
"""
function validate_pxp_ed_ipeps(
    config::PXPValidationConfig;
    ctm_params = nothing,
    trust_policy::CTMTrustPolicy = CTMTrustPolicy(),
    ctm_measure = measure_ctm,
)::PXPValidationReport
    ed_result = run_pxp_ed_benchmark(_validation_ed_config(config))
    return _validate_pxp_ipeps_against_ed(
        config,
        ed_result,
        ctm_params = ctm_params,
        trust_policy = trust_policy,
        ctm_measure = ctm_measure,
    )
end

"""
    validate_pxp_convergence(config; trust_policy = CTMTrustPolicy(),
                             ctm_measure = measure_ctm)

Run a `dt x D x cutoff` validation grid from [`PXPConvergenceConfig`](@ref)
and aggregate simple and CTM density error budgets. `config.chi_values`
provides the finite-`chi` trust sweep evaluated within each grid run. When it
is empty, CTM sweeps are skipped and CTM summary fields are `nothing`.
"""
function validate_pxp_convergence(
    config::PXPConvergenceConfig;
    trust_policy::CTMTrustPolicy = CTMTrustPolicy(),
    ctm_measure = measure_ctm,
)::PXPConvergenceReport
    ctm_params = isempty(config.chi_values) ? nothing :
        Tuple(default_ctmrg_params(; chi = chi) for chi in config.chi_values)
    runs = PXPValidationReport[]

    for dt in config.dt_values, D in config.D_values, cutoff in config.cutoff_values
        run_config = _copy_config(config.base; dt, maxdim = D, cutoff)
        push!(
            runs,
            validate_pxp_ed_ipeps(
                run_config;
                ctm_params,
                trust_policy,
                ctm_measure,
            ),
        )
    end

    simple_errors = [abs(c.density_error_simple) for r in runs for c in r.comparisons]
    exact_errors = [
        abs(c.density_error_exact_finite) for r in runs for c in r.comparisons if
        c.density_error_exact_finite !== nothing
    ]
    ctm_errors = Float64[]
    trust_flags = Bool[]
    for run in runs
        for (ed_sample, ipeps_sample) in zip(run.ed_result.samples, run.ipeps_samples)
            ipeps_sample.ctm === nothing && continue
            for point in ipeps_sample.ctm.points
                push!(ctm_errors, abs(point.measurement.density - ed_sample.excitation_density))
            end
            push!(trust_flags, ipeps_sample.ctm.trust.trusted)
        end
    end

    return PXPConvergenceReport(
        config,
        runs,
        maximum(simple_errors),
        isempty(ctm_errors) ? nothing : maximum(ctm_errors),
        _finite_max_or_nothing(exact_errors),
        isempty(trust_flags) ? nothing : all(==(true), trust_flags),
    )
end


include("PXPValidationCampaigns.jl")

include("PXPValidationSerialization.jl")

end
