module PXPValidation

using JSON3

using ..SquareIPEPS: SquareIPEPSState, log_norm, product_square_ipeps
using ..SquareUnitCells: PeriodicSquareUnitCell
using ..Observables: SimpleObservableSummary, measure_simple
using ..PEPSKitMeasurements:
    CTMObservableSummary,
    CTMValidationPoint,
    measure_ctm,
    validate_ctm_sweep
using ..CTMTrust: CTMTrustAssessment, CTMTrustPolicy, assess_ctm_trust
using ..IPEPSEvolution: EvolutionLog, TrotterParams, evolve!
using ..FinitePXPEEDBenchmark:
    PXPEEDBenchmarkConfig,
    PXPEEDBenchmarkResult,
    PXPEEDSample,
    run_pxp_ed_benchmark

export TrustedCTMMeasurement, measure_ctm_trusted
export PXPValidationConfig, PXPValidationMetadata, PXPIPEPSSample
export PXPEDComparisonSample, PXPValidationReport, validate_pxp_ed_ipeps
export write_pxp_validation_json

"""
    TrustedCTMMeasurement(measurement, points, trust)

Finite-`chi` CTMRG measurement bundle for one iPEPS state. `points` stores the
full validation sweep, `measurement` is the last sweep point's CTM observable
summary, and `trust` is the finite-`chi` assessment returned by
[`assess_ctm_trust`](@ref).
"""
struct TrustedCTMMeasurement
    measurement::CTMObservableSummary
    points::Vector{CTMValidationPoint}
    trust::CTMTrustAssessment

    function TrustedCTMMeasurement(
        measurement::CTMObservableSummary,
        points::Vector{CTMValidationPoint},
        trust::CTMTrustAssessment,
    )
        isempty(points) &&
            throw(ArgumentError("trusted CTM measurement requires at least one sweep point"))
        points[end].measurement == measurement ||
            throw(ArgumentError("measurement must match the final CTM validation point"))
        return new(measurement, points, trust)
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
    return TrustedCTMMeasurement(points[end].measurement, points, assessment)
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

"""
    PXPValidationConfig(n; kwargs...)

Controls for a short-time finite-ED versus iPEPS PXP validation run. The ED
side uses [`PXPEEDBenchmarkConfig`](@ref); the iPEPS side uses a periodic
`n x n` unit cell, all-down product initialization, real-time PXP evolution,
and a serial star schedule by default so `3 x 3` smoke validation is supported.
"""
struct PXPValidationConfig
    n::Int
    total_time::Float64
    dt::Float64
    measure_every::Int
    initial_state::Symbol
    point_group::Bool
    use_sparse::Bool
    ed_tol::Float64
    ed_m_init::Int
    ed_m_max::Int
    ed_extend_step::Int
    order::Int
    maxdim::Int
    cutoff::Float64
    schedule::Symbol

    function PXPValidationConfig(
        n::Integer;
        total_time::Real = 0.02,
        dt::Real = 0.01,
        measure_every::Integer = 1,
        initial_state::Symbol = :down,
        point_group::Bool = true,
        use_sparse::Bool = true,
        ed_tol::Real = 1e-10,
        ed_m_init::Integer = 30,
        ed_m_max::Integer = 60,
        ed_extend_step::Integer = 10,
        order::Integer = 1,
        maxdim::Integer = 1,
        cutoff::Real = 1e-12,
        schedule::Symbol = :serial,
    )
        n_int = _positive_int(n, "n")
        total = _finite_nonnegative(total_time, "total_time")
        step = _finite_positive(dt, "dt")
        cadence = _positive_int(measure_every, "measure_every")
        initial_state in (:down, :all_down) ||
            throw(ArgumentError("initial_state must be :down or :all_down"))
        tol = _finite_positive(ed_tol, "ed_tol")
        m_init = _positive_int(ed_m_init, "ed_m_init")
        m_max = _positive_int(ed_m_max, "ed_m_max")
        m_max >= m_init || throw(ArgumentError("ed_m_max must be at least ed_m_init"))
        extend = _positive_int(ed_extend_step, "ed_extend_step")
        ord = _positive_int(order, "order")
        ord in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
        dim = _positive_int(maxdim, "maxdim")
        trunc_cutoff = _finite_nonnegative(cutoff, "cutoff")
        schedule in (:serial, :five_color) ||
            throw(ArgumentError("schedule must be :serial or :five_color"))
        schedule === :five_color && n_int % 5 != 0 &&
            throw(ArgumentError("five_color validation requires n divisible by 5"))

        PXPEEDBenchmarkConfig(
            n_int;
            total_time = total,
            dt = step,
            measure_every = cadence,
            initial_state,
            point_group,
            use_sparse,
            tol,
            m_init,
            m_max,
            extend_step = extend,
        )
        TrotterParams(step, ord, :real, dim, trunc_cutoff; schedule)

        return new(
            n_int,
            total,
            step,
            cadence,
            initial_state,
            point_group,
            use_sparse,
            tol,
            m_init,
            m_max,
            extend,
            ord,
            dim,
            trunc_cutoff,
            schedule,
        )
    end
end

"""
    PXPValidationMetadata

Reproducibility metadata attached to a PXP validation report. `git_commit` is
`nothing` when the repository commit cannot be read from the current process.
"""
struct PXPValidationMetadata
    git_commit::Union{Nothing,String}
    julia_version::String
    project_path::Union{Nothing,String}
end

"""
    PXPIPEPSSample

One iPEPS validation sample at an ED measurement time. The first sample has
`evolution === nothing`; later samples store the incremental [`EvolutionLog`](@ref)
for the time interval since the previous sample. `ctm` is present only when a
finite-`chi` CTM sweep was requested.
"""
struct PXPIPEPSSample
    step::Int
    time::Float64
    simple::SimpleObservableSummary
    evolution::Union{Nothing,EvolutionLog}
    ctm::Union{Nothing,TrustedCTMMeasurement}
    log_norm::Float64
end

"""
    PXPEDComparisonSample

Observable comparison at one matched ED/iPEPS sample time. Density errors are
reported as `iPEPS - ED`; CTM fields are `nothing` when no CTM sweep was run.
"""
struct PXPEDComparisonSample
    step::Int
    time::Float64
    ed_return_probability::Float64
    ed_excitation_density::Float64
    ipeps_simple_density::Float64
    ipeps_ctm_density::Union{Nothing,Float64}
    density_error_simple::Float64
    density_error_ctm::Union{Nothing,Float64}
    simple_blockade_violation::Float64
    ctm_blockade_violation::Union{Nothing,Float64}
    ctm_trusted::Union{Nothing,Bool}
    ctm_reason::Union{Nothing,Symbol}
end

"""
    PXPValidationReport

Complete short-time finite PXP validation artifact: input config, finite ED
result, matched iPEPS trajectory samples, observable comparisons, and
reproducibility metadata.
"""
struct PXPValidationReport
    config::PXPValidationConfig
    ed_result::PXPEEDBenchmarkResult
    ipeps_samples::Vector{PXPIPEPSSample}
    comparisons::Vector{PXPEDComparisonSample}
    metadata::PXPValidationMetadata
end

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

function _validation_initial_state(config::PXPValidationConfig)
    cell = PeriodicSquareUnitCell(config.n, config.n)
    state = config.initial_state === :all_down ? :down : config.initial_state
    return product_square_ipeps(cell; state, maxdim = config.maxdim)
end

function _git_commit()
    try
        package_root = abspath(joinpath(@__DIR__, ".."))
        command = pipeline(`git -C $package_root rev-parse HEAD`; stderr = devnull)
        commit = chomp(read(command, String))
        isempty(commit) && return nothing
        return commit
    catch
        return nothing
    end
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

function _ctm_density(ctm::Nothing)
    return nothing
end

function _ctm_density(ctm::TrustedCTMMeasurement)
    return ctm.measurement.density
end

function _ctm_blockade(ctm::Nothing)
    return nothing
end

function _ctm_blockade(ctm::TrustedCTMMeasurement)
    return ctm.measurement.blockade_violation
end

function _ctm_trusted(ctm::Nothing)
    return nothing
end

function _ctm_trusted(ctm::TrustedCTMMeasurement)
    return ctm.trust.trusted
end

function _ctm_reason(ctm::Nothing)
    return nothing
end

function _ctm_reason(ctm::TrustedCTMMeasurement)
    return ctm.trust.reason
end

function _comparison(ed::PXPEEDSample, sample::PXPIPEPSSample)
    ctm_density = _ctm_density(sample.ctm)
    return PXPEDComparisonSample(
        ed.step,
        ed.time,
        ed.return_probability,
        ed.excitation_density,
        sample.simple.density,
        ctm_density,
        sample.simple.density - ed.excitation_density,
        ctm_density === nothing ? nothing : ctm_density - ed.excitation_density,
        sample.simple.blockade_violation,
        _ctm_blockade(sample.ctm),
        _ctm_trusted(sample.ctm),
        _ctm_reason(sample.ctm),
    )
end

"""
    validate_pxp_ed_ipeps(config; ctm_params = nothing,
                          trust_policy = CTMTrustPolicy(),
                          ctm_measure = measure_ctm)

Run a finite periodic PXP ED trajectory and a matched all-down iPEPS trajectory
on an `n x n` unit cell. iPEPS samples are measured at the same times as ED.
When `ctm_params` is supplied, every iPEPS sample also receives a trusted CTM
measurement bundle from [`measure_ctm_trusted`](@ref).
"""
function validate_pxp_ed_ipeps(
    config::PXPValidationConfig;
    ctm_params = nothing,
    trust_policy::CTMTrustPolicy = CTMTrustPolicy(),
    ctm_measure = measure_ctm,
)::PXPValidationReport
    ed_result = run_pxp_ed_benchmark(_validation_ed_config(config))
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
            ),
        )
    end

    comparisons = [
        _comparison(ed_sample, ipeps_sample) for
        (ed_sample, ipeps_sample) in zip(ed_result.samples, samples)
    ]
    return PXPValidationReport(
        config,
        ed_result,
        samples,
        comparisons,
        _validation_metadata(),
    )
end

function _json_value(value::Nothing)
    return nothing
end

function _json_value(value::Symbol)
    return String(value)
end

function _config_data(config::PXPValidationConfig)
    return (;
        n = config.n,
        total_time = config.total_time,
        dt = config.dt,
        measure_every = config.measure_every,
        initial_state = String(config.initial_state),
        point_group = config.point_group,
        use_sparse = config.use_sparse,
        ed_tol = config.ed_tol,
        ed_m_init = config.ed_m_init,
        ed_m_max = config.ed_m_max,
        ed_extend_step = config.ed_extend_step,
        order = config.order,
        maxdim = config.maxdim,
        cutoff = config.cutoff,
        schedule = String(config.schedule),
    )
end

function _metadata_data(metadata::PXPValidationMetadata)
    return (;
        git_commit = metadata.git_commit,
        julia_version = metadata.julia_version,
        project_path = metadata.project_path,
    )
end

function _simple_data(summary::SimpleObservableSummary)
    return (;
        density = summary.density,
        density_even = summary.density_even,
        density_odd = summary.density_odd,
        blockade_violation = summary.blockade_violation,
        pxp_energy_density = summary.pxp_energy_density,
        mean_bond_entropy = summary.mean_bond_entropy,
        max_bond_entropy = summary.max_bond_entropy,
    )
end

function _ctm_diagnostics_data(diagnostics::Nothing)
    return nothing
end

function _ctm_diagnostics_data(diagnostics)
    return (;
        chi = diagnostics.chi,
        tol = diagnostics.tol,
        maxiter = diagnostics.maxiter,
        iterations = diagnostics.iterations,
        residual = diagnostics.residual,
        converged = diagnostics.converged,
        accepted = diagnostics.accepted,
    )
end

function _ctm_summary_data(summary::CTMObservableSummary)
    return (;
        density = summary.density,
        density_even = summary.density_even,
        density_odd = summary.density_odd,
        blockade_violation = summary.blockade_violation,
        pxp_energy_density = summary.pxp_energy_density,
        diagnostics = _ctm_diagnostics_data(summary.diagnostics),
    )
end

function _ctm_point_data(point::CTMValidationPoint)
    return (;
        chi = point.params.chi,
        tol = point.params.tol,
        maxiter = point.params.maxiter,
        verbosity = point.params.verbosity,
        measurement = _ctm_summary_data(point.measurement),
        delta_density = point.delta_density,
        delta_density_even = point.delta_density_even,
        delta_density_odd = point.delta_density_odd,
        delta_blockade_violation = point.delta_blockade_violation,
        delta_pxp_energy_density = point.delta_pxp_energy_density,
    )
end

function _trust_data(trust::CTMTrustAssessment)
    return (;
        trusted = trust.trusted,
        reason = String(trust.reason),
        message = trust.message,
        compared_points = trust.compared_points,
        finite_chi_density_delta = trust.finite_chi_density_delta,
        finite_chi_blockade_delta = trust.finite_chi_blockade_delta,
        finite_chi_energy_delta = trust.finite_chi_energy_delta,
        observed_max_residual = trust.observed_max_residual,
    )
end

function _trusted_ctm_data(ctm::Nothing)
    return nothing
end

function _trusted_ctm_data(ctm::TrustedCTMMeasurement)
    return (;
        measurement = _ctm_summary_data(ctm.measurement),
        points = _ctm_point_data.(ctm.points),
        trust = _trust_data(ctm.trust),
    )
end

function _evolution_data(evolution::Nothing)
    return nothing
end

function _evolution_data(evolution::EvolutionLog)
    return (;
        total_time = evolution.total_time,
        nsteps = evolution.nsteps,
        dt = evolution.params.dt,
        order = evolution.params.order,
        evolution = String(evolution.params.evolution),
        maxdim = evolution.params.maxdim,
        cutoff = evolution.params.cutoff,
        schedule = String(evolution.params.schedule),
        max_truncerr = evolution.max_truncerr,
        max_bond_entropy = evolution.max_bond_entropy,
        mean_bond_entropy = evolution.mean_bond_entropy,
        log_norm_before = evolution.log_norm_before,
        log_norm_after = evolution.log_norm_after,
        log_norm_delta = evolution.log_norm_delta,
        model_metadata = evolution.model_metadata,
    )
end

function _ed_sample_data(sample::PXPEEDSample)
    return (;
        step = sample.step,
        time = sample.time,
        norm = sample.norm,
        return_probability = sample.return_probability,
        excitation_density = sample.excitation_density,
    )
end

function _ed_diagnostics_data(diagnostics)
    return (;
        basis_builds = diagnostics.basis_builds,
        basis_extensions = diagnostics.basis_extensions,
        restarts = diagnostics.restarts,
        matvecs = diagnostics.matvecs,
        total_times_served = diagnostics.total_times_served,
        max_dim_used = diagnostics.max_dim_used,
        accepted_intervals = diagnostics.accepted_intervals,
    )
end

function _ed_result_data(result::PXPEEDBenchmarkResult)
    return (;
        lattice_size = result.lattice_size,
        basis_dimension = result.basis_dimension,
        constrained_dimension = result.constrained_dimension,
        group_order = result.group_order,
        point_group = result.point_group,
        hamiltonian_nnz = result.hamiltonian_nnz,
        samples = _ed_sample_data.(result.samples),
        diagnostics = _ed_diagnostics_data(result.diagnostics),
    )
end

function _ipeps_sample_data(sample::PXPIPEPSSample)
    return (;
        step = sample.step,
        time = sample.time,
        simple = _simple_data(sample.simple),
        evolution = _evolution_data(sample.evolution),
        ctm = _trusted_ctm_data(sample.ctm),
        log_norm = sample.log_norm,
    )
end

function _comparison_data(sample::PXPEDComparisonSample)
    return (;
        step = sample.step,
        time = sample.time,
        ed_return_probability = sample.ed_return_probability,
        ed_excitation_density = sample.ed_excitation_density,
        ipeps_simple_density = sample.ipeps_simple_density,
        ipeps_ctm_density = sample.ipeps_ctm_density,
        density_error_simple = sample.density_error_simple,
        density_error_ctm = sample.density_error_ctm,
        simple_blockade_violation = sample.simple_blockade_violation,
        ctm_blockade_violation = sample.ctm_blockade_violation,
        ctm_trusted = sample.ctm_trusted,
        ctm_reason = _json_value(sample.ctm_reason),
    )
end

function _report_data(report::PXPValidationReport)
    return (;
        config = _config_data(report.config),
        metadata = _metadata_data(report.metadata),
        ed_result = _ed_result_data(report.ed_result),
        ipeps_samples = _ipeps_sample_data.(report.ipeps_samples),
        comparisons = _comparison_data.(report.comparisons),
    )
end

"""
    write_pxp_validation_json(report, path)

Write a [`PXPValidationReport`](@ref) to `path` as a JSON artifact containing
configuration, metadata, ED samples, iPEPS samples, CTM trust data when present,
and matched observable comparisons.
"""
function write_pxp_validation_json(report::PXPValidationReport, path::AbstractString)
    open(path, "w") do io
        JSON3.write(io, _report_data(report))
        write(io, '\n')
    end
    return path
end

end
