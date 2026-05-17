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
    exact_finite_observables::Bool
    exact_finite_max_sites::Int

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
        exact_finite_observables::Bool = false,
        exact_finite_max_sites::Integer = 12,
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
        exact_limit = _positive_int(exact_finite_max_sites, "exact_finite_max_sites")
        exact_finite_observables && n_int^2 > exact_limit &&
            throw(ArgumentError("exact finite observables require n^2 <= exact_finite_max_sites"))

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
            exact_finite_observables,
            exact_limit,
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
    exact_finite_density::Union{Nothing,Float64}
end

function PXPIPEPSSample(step, time, simple, evolution, ctm, log_norm)
    return PXPIPEPSSample(step, time, simple, evolution, ctm, log_norm, nothing)
end

"""
    PXPEDComparisonSample

Observable comparison at one matched ED/iPEPS sample time. Density errors are
reported as `iPEPS - ED`; CTM fields are `nothing` when no CTM sweep was run.
For D>1 no-CTM runs, `density_error_simple` is a simple/local-environment
diagnostic offset, not an exact finite-PEPS validation error.
"""
struct PXPEDComparisonSample
    step::Int
    time::Float64
    ed_return_probability::Float64
    ed_excitation_density::Float64
    ipeps_simple_density::Float64
    ipeps_ctm_density::Union{Nothing,Float64}
    ipeps_exact_finite_density::Union{Nothing,Float64}
    density_error_simple::Float64
    density_error_ctm::Union{Nothing,Float64}
    density_error_exact_finite::Union{Nothing,Float64}
    simple_blockade_violation::Float64
    ctm_blockade_violation::Union{Nothing,Float64}
    ctm_trusted::Union{Nothing,Bool}
    ctm_reason::Union{Nothing,Symbol}
end

function PXPEDComparisonSample(
    step,
    time,
    ed_return_probability,
    ed_excitation_density,
    ipeps_simple_density,
    ipeps_ctm_density,
    density_error_simple,
    density_error_ctm,
    simple_blockade_violation,
    ctm_blockade_violation,
    ctm_trusted,
    ctm_reason,
)
    return PXPEDComparisonSample(
        step,
        time,
        ed_return_probability,
        ed_excitation_density,
        ipeps_simple_density,
        ipeps_ctm_density,
        nothing,
        density_error_simple,
        density_error_ctm,
        nothing,
        simple_blockade_violation,
        ctm_blockade_violation,
        ctm_trusted,
        ctm_reason,
    )
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

"""
    PXPConvergenceConfig(base; dt_values, D_values, chi_values = Int[],
                         cutoff_values)

Parameter grid for repeated PXP ED-versus-iPEPS validation runs. `base`
provides all controls not swept explicitly. `dt_values`, `D_values`, and
`cutoff_values` define the run grid and must be nonempty and positive.
`chi_values` are positive finite-`chi` CTM trust-sweep values for each run;
leave them empty to skip CTM comparisons.
"""
struct PXPConvergenceConfig
    base::PXPValidationConfig
    dt_values::Vector{Float64}
    D_values::Vector{Int}
    chi_values::Vector{Int}
    cutoff_values::Vector{Float64}

    function PXPConvergenceConfig(
        base::PXPValidationConfig;
        dt_values,
        D_values,
        chi_values = Int[],
        cutoff_values,
    )
        dt_grid = [_finite_positive(dt, "dt_values entry") for dt in dt_values]
        D_grid = [_positive_int(D, "D_values entry") for D in D_values]
        chi_grid = [_positive_int(chi, "chi_values entry") for chi in chi_values]
        cutoff_grid = [_finite_positive(cutoff, "cutoff_values entry") for cutoff in cutoff_values]

        isempty(dt_grid) && throw(ArgumentError("dt_values must be nonempty"))
        isempty(D_grid) && throw(ArgumentError("D_values must be nonempty"))
        isempty(cutoff_grid) && throw(ArgumentError("cutoff_values must be nonempty"))

        return new(base, dt_grid, D_grid, chi_grid, cutoff_grid)
    end
end

"""
    PXPConvergenceReport

Machine-readable convergence/error-budget report containing every validation
run in a `dt x D x cutoff` sweep, the maximum absolute simple-density error,
and optional CTM density-error and trust summaries when finite-`chi` CTM trust
sweeps were requested. CTM density error is aggregated across all validation
points in each trust sweep, not just the final `chi`. Simple-density errors are
cheap local-environment audit signals; they are not exact finite contractions
for D>1 loopy periodic PEPS.
"""
struct PXPConvergenceReport
    config::PXPConvergenceConfig
    runs::Vector{PXPValidationReport}
    max_abs_density_error_simple::Float64
    max_abs_density_error_ctm::Union{Nothing,Float64}
    max_abs_density_error_exact_finite::Union{Nothing,Float64}
    all_ctm_trusted::Union{Nothing,Bool}
end

function PXPConvergenceReport(
    config,
    runs,
    max_abs_density_error_simple,
    max_abs_density_error_ctm,
    all_ctm_trusted,
)
    return PXPConvergenceReport(
        config,
        runs,
        max_abs_density_error_simple,
        max_abs_density_error_ctm,
        nothing,
        all_ctm_trusted,
    )
end

"""
    PXPReversibilityReport

Simple-observable reversibility diagnostics for a forward real-time PXP
evolution followed by [`reverse_evolve!`](@ref). `before`, `after_forward`, and
`after_reverse` store the simple measurements around the protocol, the log
fields store both evolution calls, and the drift fields are absolute differences
between `before` and `after_reverse`.
"""
struct PXPReversibilityReport
    before::SimpleObservableSummary
    after_forward::SimpleObservableSummary
    after_reverse::SimpleObservableSummary
    forward_log::EvolutionLog
    reverse_log::EvolutionLog
    density_drift::Float64
    blockade_drift::Float64
    energy_drift::Float64
end

"""
    PXPAuditConfig(; kwargs...)

Small PXP audit-campaign grid for ED/iPEPS, optional trusted CTM attachment,
and reversibility diagnostics. Defaults are deliberately short-time and
serial-schedule oriented: `n_values = [3]`, `total_time = 0.02`,
`dt_values = [0.02, 0.01]`, `D_values = [1, 2]`,
`cutoff_values = [1e-12]`, and no CTM sweep unless `chi_values` is supplied.
"""
struct PXPAuditConfig
    n_values::Vector{Int}
    total_time::Float64
    dt_values::Vector{Float64}
    D_values::Vector{Int}
    cutoff_values::Vector{Float64}
    chi_values::Vector{Int}
    measure_every::Int
    order::Int
    schedule::Symbol
    ctm_tol::Float64
    ctm_maxiter::Int
    ctm_verbosity::Int
    ctm_seed::Union{Nothing,Int}
    exact_finite_observables::Bool
    exact_finite_max_sites::Int

    function PXPAuditConfig(;
        n_values = [3],
        total_time::Real = 0.02,
        dt_values = [0.02, 0.01],
        D_values = [1, 2],
        cutoff_values = [1e-12],
        chi_values = Int[],
        measure_every::Integer = 1,
        order::Integer = 1,
        schedule::Symbol = :serial,
        ctm_tol::Real = 1e-8,
        ctm_maxiter::Integer = 100,
        ctm_verbosity::Integer = 0,
        ctm_seed = 0,
        exact_finite_observables::Bool = false,
        exact_finite_max_sites::Integer = 12,
    )
        n_grid = [_positive_int(n, "n_values entry") for n in n_values]
        dt_grid = [_finite_positive(dt, "dt_values entry") for dt in dt_values]
        D_grid = [_positive_int(D, "D_values entry") for D in D_values]
        cutoff_grid = [_finite_positive(cutoff, "cutoff_values entry") for cutoff in cutoff_values]
        chi_grid = [_positive_int(chi, "chi_values entry") for chi in chi_values]
        isempty(n_grid) && throw(ArgumentError("n_values must be nonempty"))
        isempty(dt_grid) && throw(ArgumentError("dt_values must be nonempty"))
        isempty(D_grid) && throw(ArgumentError("D_values must be nonempty"))
        isempty(cutoff_grid) && throw(ArgumentError("cutoff_values must be nonempty"))
        cadence = _positive_int(measure_every, "measure_every")
        ord = _positive_int(order, "order")
        ord in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
        schedule in (:serial, :five_color) ||
            throw(ArgumentError("schedule must be :serial or :five_color"))
        tol = _finite_positive(ctm_tol, "ctm_tol")
        maxiter = _positive_int(ctm_maxiter, "ctm_maxiter")
        verbosity = _nonnegative_int(ctm_verbosity, "ctm_verbosity")
        seed = ctm_seed === nothing ? nothing : _nonnegative_int(ctm_seed, "ctm_seed")
        total = _finite_nonnegative(total_time, "total_time")
        exact_limit = _positive_int(exact_finite_max_sites, "exact_finite_max_sites")

        for n in n_grid, dt in dt_grid, D in D_grid, cutoff in cutoff_grid
            PXPValidationConfig(
                n;
                total_time = total,
                dt,
                measure_every = cadence,
                order = ord,
                maxdim = D,
                cutoff,
                schedule,
                exact_finite_observables,
                exact_finite_max_sites = exact_limit,
            )
        end

        return new(
            n_grid,
            total,
            dt_grid,
            D_grid,
            cutoff_grid,
            chi_grid,
            cadence,
            ord,
            schedule,
            tol,
            maxiter,
            verbosity,
            seed,
            exact_finite_observables,
            exact_limit,
        )
    end
end

"""
    PXPAuditSummary

Flat per-run audit row intended for JSON and CSV bottleneck triage. CTM fields
are `nothing` or `:not_run` when the audit configuration does not request a
finite-`chi` CTM sweep. The simple fields summarize no-CTM local diagnostics,
not exact finite-PEPS observables for D>1.
"""
struct PXPAuditSummary
    n::Int
    total_time::Float64
    dt::Float64
    D::Int
    cutoff::Float64
    schedule::Symbol
    chi_values::Vector{Int}
    max_abs_density_error_simple::Float64
    max_abs_density_error_ctm::Union{Nothing,Float64}
    max_abs_density_error_exact_finite::Union{Nothing,Float64}
    max_blockade_violation_simple::Float64
    max_blockade_violation_ctm::Union{Nothing,Float64}
    pxp_energy_drift_simple::Float64
    pxp_energy_drift_ctm::Union{Nothing,Float64}
    ctm_trust_status::Symbol
    ctm_trust_reason::Symbol
    finite_chi_density_delta::Union{Nothing,Float64}
    finite_chi_blockade_delta::Union{Nothing,Float64}
    finite_chi_energy_delta::Union{Nothing,Float64}
    finite_chi_max_residual::Union{Nothing,Float64}
    max_truncerr::Float64
    log_norm_initial::Float64
    log_norm_final::Float64
    log_norm_delta::Float64
    log_norm_delta_abs::Float64
    reversibility_density_drift::Float64
    reversibility_blockade_drift::Float64
    reversibility_energy_drift::Float64
end

function PXPAuditSummary(
    n,
    total_time,
    dt,
    D,
    cutoff,
    schedule,
    chi_values,
    max_abs_density_error_simple,
    max_abs_density_error_ctm,
    max_blockade_violation_simple,
    max_blockade_violation_ctm,
    pxp_energy_drift_simple,
    pxp_energy_drift_ctm,
    ctm_trust_status,
    ctm_trust_reason,
    finite_chi_density_delta,
    finite_chi_blockade_delta,
    finite_chi_energy_delta,
    finite_chi_max_residual,
    max_truncerr,
    log_norm_initial,
    log_norm_final,
    log_norm_delta,
    log_norm_delta_abs,
    reversibility_density_drift,
    reversibility_blockade_drift,
    reversibility_energy_drift,
)
    return PXPAuditSummary(
        n,
        total_time,
        dt,
        D,
        cutoff,
        schedule,
        chi_values,
        max_abs_density_error_simple,
        max_abs_density_error_ctm,
        nothing,
        max_blockade_violation_simple,
        max_blockade_violation_ctm,
        pxp_energy_drift_simple,
        pxp_energy_drift_ctm,
        ctm_trust_status,
        ctm_trust_reason,
        finite_chi_density_delta,
        finite_chi_blockade_delta,
        finite_chi_energy_delta,
        finite_chi_max_residual,
        max_truncerr,
        log_norm_initial,
        log_norm_final,
        log_norm_delta,
        log_norm_delta_abs,
        reversibility_density_drift,
        reversibility_blockade_drift,
        reversibility_energy_drift,
    )
end

"""
    PXPAuditRun

One audit-grid point containing the full ED/iPEPS validation report, the
matched reversibility report, and the flattened [`PXPAuditSummary`](@ref).
"""
struct PXPAuditRun
    validation::PXPValidationReport
    reversibility::PXPReversibilityReport
    summary::PXPAuditSummary
end

"""
    PXPAuditReport

Machine-readable PXP audit-campaign artifact containing the campaign
configuration, per-run validation and reversibility reports, and reproducibility
metadata. Use [`write_pxp_audit_json`](@ref) for nested JSON and
[`write_pxp_audit_csv`](@ref) for the flat summary table.
"""
struct PXPAuditReport
    config::PXPAuditConfig
    runs::Vector{PXPAuditRun}
    metadata::PXPValidationMetadata
end

"""
    PXPLargerDBenchmarkConfig(; kwargs...)

Controls the M3 larger-D PXP dynamics benchmark. `ed_mode = :symmetric_pbc`
uses the current finite PBC ED path, and `observable_mode = :auto` selects
`:exact_finite` for tiny cells with exact finite observables enabled and
`:symmetric_pbc_ed_global` otherwise.
"""
struct PXPLargerDBenchmarkConfig
    n_values::Vector{Int}
    total_time::Float64
    dt_values::Vector{Float64}
    D_values::Vector{Int}
    cutoff_values::Vector{Float64}
    measure_every::Int
    order::Int
    schedule::Symbol
    initial_state::Symbol
    point_group::Bool
    use_sparse::Bool
    ed_tol::Float64
    ed_m_init::Int
    ed_m_max::Int
    ed_extend_step::Int
    ed_mode::Symbol
    observable_mode::Symbol
    chi_values::Vector{Int}
    ctm_tol::Float64
    ctm_maxiter::Int
    ctm_verbosity::Int
    ctm_seed::Union{Nothing,Int}
    exact_finite_observables::Bool
    exact_finite_max_sites::Int
end

function PXPLargerDBenchmarkConfig(;
    n_values = [3],
    total_time::Real = 0.02,
    dt_values = [0.02],
    D_values = [1, 2, 3, 4],
    cutoff_values = [1e-12],
    measure_every::Integer = 1,
    order::Integer = 1,
    schedule::Symbol = :serial,
    initial_state::Symbol = :down,
    point_group::Bool = true,
    use_sparse::Bool = true,
    ed_tol::Real = 1e-10,
    ed_m_init::Integer = 30,
    ed_m_max::Integer = 60,
    ed_extend_step::Integer = 10,
    ed_mode::Symbol = :symmetric_pbc,
    observable_mode::Symbol = :auto,
    chi_values = Int[],
    ctm_tol::Real = 1e-8,
    ctm_maxiter::Integer = 100,
    ctm_verbosity::Integer = 0,
    ctm_seed = 0,
    exact_finite_observables::Bool = false,
    exact_finite_max_sites::Integer = 12,
)
    n_grid = [_positive_int(n, "n_values entry") for n in n_values]
    dt_grid = [_finite_positive(dt, "dt_values entry") for dt in dt_values]
    D_grid = [_positive_int(D, "D_values entry") for D in D_values]
    cutoff_grid = [_finite_positive(cutoff, "cutoff_values entry") for cutoff in cutoff_values]
    chi_grid = [_positive_int(chi, "chi_values entry") for chi in chi_values]
    isempty(n_grid) && throw(ArgumentError("n_values must be nonempty"))
    isempty(dt_grid) && throw(ArgumentError("dt_values must be nonempty"))
    isempty(D_grid) && throw(ArgumentError("D_values must be nonempty"))
    isempty(cutoff_grid) && throw(ArgumentError("cutoff_values must be nonempty"))
    ed_mode === :symmetric_pbc ||
        throw(ArgumentError("ed_mode must be :symmetric_pbc for M3"))
    observable_mode in (:auto, :exact_finite, :symmetric_pbc_ed_global, :ctm_trusted) ||
        throw(
            ArgumentError(
                "observable_mode must be :auto, :exact_finite, :symmetric_pbc_ed_global, or :ctm_trusted",
            ),
        )
    initial_state in (:down, :all_down) ||
        throw(ArgumentError("initial_state must be :down or :all_down"))

    cadence = _positive_int(measure_every, "measure_every")
    ord = _positive_int(order, "order")
    ord in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    schedule in (:serial, :five_color) ||
        throw(ArgumentError("schedule must be :serial or :five_color"))
    total = _finite_nonnegative(total_time, "total_time")
    tol = _finite_positive(ed_tol, "ed_tol")
    m_init = _positive_int(ed_m_init, "ed_m_init")
    m_max = _positive_int(ed_m_max, "ed_m_max")
    m_max >= m_init || throw(ArgumentError("ed_m_max must be at least ed_m_init"))
    extend = _positive_int(ed_extend_step, "ed_extend_step")
    ctm_tol_f = _finite_positive(ctm_tol, "ctm_tol")
    ctm_maxiter_i = _positive_int(ctm_maxiter, "ctm_maxiter")
    ctm_verbosity_i = _nonnegative_int(ctm_verbosity, "ctm_verbosity")
    seed = ctm_seed === nothing ? nothing : _nonnegative_int(ctm_seed, "ctm_seed")
    exact_limit = _positive_int(exact_finite_max_sites, "exact_finite_max_sites")

    for n in n_grid, dt in dt_grid, D in D_grid, cutoff in cutoff_grid
        PXPValidationConfig(
            n;
            total_time = total,
            dt,
            measure_every = cadence,
            initial_state,
            point_group,
            use_sparse,
            ed_tol = tol,
            ed_m_init = m_init,
            ed_m_max = m_max,
            ed_extend_step = extend,
            order = ord,
            maxdim = D,
            cutoff,
            schedule,
            exact_finite_observables = exact_finite_observables && n^2 <= exact_limit,
            exact_finite_max_sites = exact_limit,
        )
    end

    return PXPLargerDBenchmarkConfig(
        n_grid,
        total,
        dt_grid,
        D_grid,
        cutoff_grid,
        cadence,
        ord,
        schedule,
        initial_state,
        point_group,
        use_sparse,
        tol,
        m_init,
        m_max,
        extend,
        ed_mode,
        observable_mode,
        chi_grid,
        ctm_tol_f,
        ctm_maxiter_i,
        ctm_verbosity_i,
        seed,
        exact_finite_observables,
        exact_limit,
    )
end

"""
    PXPLargerDBenchmarkSummary

Flat per-run M3 benchmark row. ED fields describe the finite PBC symmetric ED
reference; exact finite fields are nullable and are populated only for tiny
iPEPS cells where exact finite contraction was enabled.
"""
struct PXPLargerDBenchmarkSummary
    n::Int
    D::Int
    dt::Float64
    cutoff::Float64
    total_time::Float64
    ed_mode::Symbol
    observable_mode::Symbol
    ed_boundary_condition::Symbol
    ed_symmetry_sector::Symbol
    ed_observable_scope::Symbol
    ed_reference_label::String
    ed_basis_dimension::Int
    ed_constrained_dimension::Int
    ed_group_order::Int
    ed_hamiltonian_nnz::Union{Nothing,Int}
    ed_runtime_seconds::Float64
    ipeps_runtime_seconds::Float64
    reversibility_runtime_seconds::Float64
    density_error_simple::Float64
    density_error_exact_finite::Union{Nothing,Float64}
    density_error_ctm::Union{Nothing,Float64}
    return_probability_error::Union{Nothing,Float64}
    ed_return_probability::Float64
    ed_excitation_density::Float64
    ipeps_simple_density::Float64
    ipeps_exact_finite_density::Union{Nothing,Float64}
    ipeps_ctm_density::Union{Nothing,Float64}
    max_truncerr::Float64
    log_norm_initial::Float64
    log_norm_final::Float64
    log_norm_delta_abs::Float64
    reversibility_density_drift::Float64
    ctm_trust_status::Symbol
    ctm_trust_reason::Symbol
    notes::Vector{String}
    warnings::Vector{String}
end

"""One M3 benchmark run containing the full validation report and flat summary."""
struct PXPLargerDBenchmarkRun
    validation::PXPValidationReport
    reversibility::PXPReversibilityReport
    summary::PXPLargerDBenchmarkSummary
end

"""Complete M3 benchmark campaign report."""
struct PXPLargerDBenchmarkReport
    config::PXPLargerDBenchmarkConfig
    runs::Vector{PXPLargerDBenchmarkRun}
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
        _ctm_blockade(sample.ctm),
        _ctm_trusted(sample.ctm),
        _ctm_reason(sample.ctm),
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
        Tuple(PEPSKitCTMRGParams(chi, 1e-8, 100, 0) for chi in config.chi_values)
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

function _audit_ctm_params(config::PXPAuditConfig)
    isempty(config.chi_values) && return nothing
    if config.ctm_seed === nothing
        return Tuple(
            PEPSKitCTMRGParams(chi, config.ctm_tol, config.ctm_maxiter, config.ctm_verbosity)
            for chi in config.chi_values
        )
    else
        return Tuple(
            PEPSKitCTMRGParams(
                chi,
                config.ctm_tol,
                config.ctm_maxiter,
                config.ctm_verbosity;
                seed = config.ctm_seed,
            ) for chi in config.chi_values
        )
    end
end

function _audit_validation_config(
    config::PXPAuditConfig,
    n::Int,
    dt::Float64,
    D::Int,
    cutoff::Float64,
)
    return PXPValidationConfig(
        n;
        total_time = config.total_time,
        dt,
        measure_every = config.measure_every,
        order = config.order,
        maxdim = D,
        cutoff,
        schedule = config.schedule,
        exact_finite_observables = config.exact_finite_observables,
        exact_finite_max_sites = config.exact_finite_max_sites,
    )
end

function _maximum_or_zero(values)
    isempty(values) && return 0.0
    return maximum(values)
end

function _minimum_or_zero(values)
    isempty(values) && return 0.0
    return minimum(values)
end

function _finite_max_or_nothing(values)
    isempty(values) && return nothing
    return maximum(values)
end

function _finite_chi_max(samples, field::Symbol)
    values = Float64[]
    for sample in samples
        sample.ctm === nothing && continue
        value = getfield(sample.ctm.trust, field)
        value === nothing && continue
        push!(values, value)
    end
    return _finite_max_or_nothing(values)
end

function _audit_trust_status(samples)
    trust_values = [sample.ctm.trust.trusted for sample in samples if sample.ctm !== nothing]
    isempty(trust_values) && return (:not_run, :not_run)
    all(==(true), trust_values) && return (:trusted, :trusted)
    for sample in samples
        sample.ctm === nothing && continue
        sample.ctm.trust.trusted || return (:rejected, sample.ctm.trust.reason)
    end
    return (:rejected, :unknown)
end

function _larger_d_ctm_params(config::PXPLargerDBenchmarkConfig)
    isempty(config.chi_values) && return nothing
    if config.ctm_seed === nothing
        return Tuple(
            PEPSKitCTMRGParams(chi, config.ctm_tol, config.ctm_maxiter, config.ctm_verbosity)
            for chi in config.chi_values
        )
    else
        return Tuple(
            PEPSKitCTMRGParams(
                chi,
                config.ctm_tol,
                config.ctm_maxiter,
                config.ctm_verbosity;
                seed = config.ctm_seed,
            ) for chi in config.chi_values
        )
    end
end

function _larger_d_ed_config(config::PXPLargerDBenchmarkConfig, n::Int, dt::Float64)
    return PXPEEDBenchmarkConfig(
        n;
        total_time = config.total_time,
        dt,
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

function _larger_d_validation_config(
    config::PXPLargerDBenchmarkConfig,
    n::Int,
    dt::Float64,
    D::Int,
    cutoff::Float64,
)
    use_exact = config.exact_finite_observables && n^2 <= config.exact_finite_max_sites
    return PXPValidationConfig(
        n;
        total_time = config.total_time,
        dt,
        measure_every = config.measure_every,
        initial_state = config.initial_state,
        point_group = config.point_group,
        use_sparse = config.use_sparse,
        ed_tol = config.ed_tol,
        ed_m_init = config.ed_m_init,
        ed_m_max = config.ed_m_max,
        ed_extend_step = config.ed_extend_step,
        order = config.order,
        maxdim = D,
        cutoff,
        schedule = config.schedule,
        exact_finite_observables = use_exact,
        exact_finite_max_sites = config.exact_finite_max_sites,
    )
end

function _larger_d_observable_mode(
    config::PXPLargerDBenchmarkConfig,
    run_config::PXPValidationConfig,
)
    config.observable_mode !== :auto && return config.observable_mode
    run_config.exact_finite_observables && return :exact_finite
    !isempty(config.chi_values) && return :ctm_trusted
    return :symmetric_pbc_ed_global
end

function _last_evolution_max_truncerr(samples)
    evolutions = [sample.evolution for sample in samples if sample.evolution !== nothing]
    return _maximum_or_zero([evolution.max_truncerr for evolution in evolutions])
end

function _exact_return_probability_or_nothing(
    sample::PXPIPEPSSample,
    config::PXPValidationConfig,
)
    config.exact_finite_observables || return nothing
    psi = _validation_initial_state(config)
    evolve!(psi, sample.time; params = _validation_trotter(config))
    return exact_all_down_return_probability_finite(psi; max_sites = config.exact_finite_max_sites)
end

function _larger_d_summary(
    config::PXPLargerDBenchmarkConfig,
    run_config::PXPValidationConfig,
    validation::PXPValidationReport,
    reversibility::PXPReversibilityReport,
    ed_seconds::Float64,
    ipeps_seconds::Float64,
    reversibility_seconds::Float64,
)::PXPLargerDBenchmarkSummary
    final_comparison = validation.comparisons[end]
    final_sample = validation.ipeps_samples[end]
    final_ed = validation.ed_result.samples[end]
    log_norms = [sample.log_norm for sample in validation.ipeps_samples]
    trust_status, trust_reason = _audit_trust_status(validation.ipeps_samples)
    mode = _larger_d_observable_mode(config, run_config)
    exact_return = final_sample.exact_finite_density === nothing ?
        nothing : _exact_return_probability_or_nothing(final_sample, run_config)

    warnings = String[]
    run_config.maxdim > 1 && push!(
        warnings,
        "density_simple is a diagnostic for D>1 loopy PEPS, not exact finite truth",
    )
    run_config.exact_finite_observables || push!(
        warnings,
        "exact finite iPEPS observables were not available for this run",
    )
    push!(warnings, "symmetric PBC ED observables are global site averages")

    return PXPLargerDBenchmarkSummary(
        run_config.n,
        run_config.maxdim,
        run_config.dt,
        run_config.cutoff,
        run_config.total_time,
        config.ed_mode,
        mode,
        :periodic,
        validation.ed_result.point_group ? :fully_symmetric_space_group : :translation_symmetric,
        :pbc_global_site_average,
        "finite_pbc_global_density",
        validation.ed_result.basis_dimension,
        validation.ed_result.constrained_dimension,
        validation.ed_result.group_order,
        validation.ed_result.hamiltonian_nnz,
        ed_seconds,
        ipeps_seconds,
        reversibility_seconds,
        final_comparison.density_error_simple,
        final_comparison.density_error_exact_finite,
        final_comparison.density_error_ctm,
        exact_return === nothing ? nothing : exact_return - final_ed.return_probability,
        final_ed.return_probability,
        final_ed.excitation_density,
        final_comparison.ipeps_simple_density,
        final_comparison.ipeps_exact_finite_density,
        final_comparison.ipeps_ctm_density,
        _last_evolution_max_truncerr(validation.ipeps_samples),
        isempty(log_norms) ? 0.0 : first(log_norms),
        isempty(log_norms) ? 0.0 : last(log_norms),
        isempty(log_norms) ? 0.0 : abs(last(log_norms) - first(log_norms)),
        reversibility.density_drift,
        trust_status,
        trust_reason,
        ["M3 larger-D PXP ED benchmark"],
        warnings,
    )
end

function _audit_summary(
    validation::PXPValidationReport,
    reversibility::PXPReversibilityReport,
)::PXPAuditSummary
    config = validation.config
    simple_density_errors = [abs(c.density_error_simple) for c in validation.comparisons]
    ctm_density_errors = [
        abs(c.density_error_ctm) for c in validation.comparisons if c.density_error_ctm !== nothing
    ]
    exact_density_errors = [
        abs(c.density_error_exact_finite) for c in validation.comparisons if
        c.density_error_exact_finite !== nothing
    ]
    simple_blockade = [c.simple_blockade_violation for c in validation.comparisons]
    ctm_blockade = [
        c.ctm_blockade_violation for c in validation.comparisons if c.ctm_blockade_violation !== nothing
    ]
    simple_energy = [sample.simple.pxp_energy_density for sample in validation.ipeps_samples]
    ctm_energy = [
        sample.ctm.measurement.pxp_energy_density for
        sample in validation.ipeps_samples if sample.ctm !== nothing
    ]
    evolutions = [
        sample.evolution for sample in validation.ipeps_samples if sample.evolution !== nothing
    ]
    log_norms = [sample.log_norm for sample in validation.ipeps_samples]
    trust_status, trust_reason = _audit_trust_status(validation.ipeps_samples)
    chi_values = Int[]
    for sample in validation.ipeps_samples
        sample.ctm === nothing && continue
        for point in sample.ctm.points
            point.params.chi in chi_values || push!(chi_values, point.params.chi)
        end
    end

    return PXPAuditSummary(
        config.n,
        config.total_time,
        config.dt,
        config.maxdim,
        config.cutoff,
        config.schedule,
        chi_values,
        _maximum_or_zero(simple_density_errors),
        _finite_max_or_nothing(ctm_density_errors),
        _finite_max_or_nothing(exact_density_errors),
        _maximum_or_zero(simple_blockade),
        _finite_max_or_nothing(ctm_blockade),
        _maximum_or_zero(simple_energy) - _minimum_or_zero(simple_energy),
        isempty(ctm_energy) ? nothing : maximum(ctm_energy) - minimum(ctm_energy),
        trust_status,
        trust_reason,
        _finite_chi_max(validation.ipeps_samples, :finite_chi_density_delta),
        _finite_chi_max(validation.ipeps_samples, :finite_chi_blockade_delta),
        _finite_chi_max(validation.ipeps_samples, :finite_chi_energy_delta),
        _finite_chi_max(validation.ipeps_samples, :observed_max_residual),
        _maximum_or_zero([evolution.max_truncerr for evolution in evolutions]),
        isempty(log_norms) ? 0.0 : first(log_norms),
        isempty(log_norms) ? 0.0 : last(log_norms),
        isempty(log_norms) ? 0.0 : last(log_norms) - first(log_norms),
        isempty(log_norms) ? 0.0 : abs(last(log_norms) - first(log_norms)),
        reversibility.density_drift,
        reversibility.blockade_drift,
        reversibility.energy_drift,
    )
end

"""
    run_pxp_audit_campaign(config = PXPAuditConfig(); trust_policy = CTMTrustPolicy(),
                           ctm_measure = measure_ctm)

Run the small M1 PXP audit grid. Each grid point runs
[`validate_pxp_ed_ipeps`](@ref) for the all-down initial state, optionally
attaches trusted finite-`chi` CTM sweeps when `config.chi_values` is nonempty,
then runs [`validate_pxp_reversibility`](@ref) with matching Trotter controls.
The result is an audit report with full nested reports and one flat summary row
per grid point.
"""
function run_pxp_audit_campaign(
    config::PXPAuditConfig = PXPAuditConfig();
    trust_policy::CTMTrustPolicy = CTMTrustPolicy(),
    ctm_measure = measure_ctm,
)::PXPAuditReport
    ctm_params = _audit_ctm_params(config)
    runs = PXPAuditRun[]

    for n in config.n_values, dt in config.dt_values, D in config.D_values,
        cutoff in config.cutoff_values
        validation_config = _audit_validation_config(config, n, dt, D, cutoff)
        validation = validate_pxp_ed_ipeps(
            validation_config;
            ctm_params,
            trust_policy,
            ctm_measure,
        )
        psi = _validation_initial_state(validation_config)
        reversibility = validate_pxp_reversibility(
            psi,
            validation_config.total_time;
            params = _validation_trotter(validation_config),
        )
        push!(runs, PXPAuditRun(validation, reversibility, _audit_summary(validation, reversibility)))
    end

    return PXPAuditReport(config, runs, _validation_metadata())
end

"""
    run_pxp_larger_d_benchmark(config = PXPLargerDBenchmarkConfig(); kwargs...)

Run the M3 larger-D PXP ED benchmark campaign. ED is run once for each
`(n, dt)` pair and reused for every requested `D` and cutoff at that pair.
"""
function run_pxp_larger_d_benchmark(
    config::PXPLargerDBenchmarkConfig = PXPLargerDBenchmarkConfig();
    trust_policy::CTMTrustPolicy = CTMTrustPolicy(),
    ctm_measure = measure_ctm,
)::PXPLargerDBenchmarkReport
    ctm_params = _larger_d_ctm_params(config)
    runs = PXPLargerDBenchmarkRun[]

    for n in config.n_values, dt in config.dt_values
        ed_result = nothing
        ed_seconds = @elapsed begin
            ed_result = run_pxp_ed_benchmark(_larger_d_ed_config(config, n, dt))
        end
        for D in config.D_values, cutoff in config.cutoff_values
            run_config = _larger_d_validation_config(config, n, dt, D, cutoff)
            validation = nothing
            ipeps_seconds = @elapsed begin
                validation = _validate_pxp_ipeps_against_ed(
                    run_config,
                    ed_result;
                    ctm_params,
                    trust_policy,
                    ctm_measure,
                )
            end
            psi = _validation_initial_state(run_config)
            reversibility = nothing
            reversibility_seconds = @elapsed begin
                reversibility = validate_pxp_reversibility(
                    psi,
                    run_config.total_time;
                    params = _validation_trotter(run_config),
                )
            end
            push!(
                runs,
                PXPLargerDBenchmarkRun(
                    validation,
                    reversibility,
                    _larger_d_summary(
                        config,
                        run_config,
                        validation,
                        reversibility,
                        ed_seconds,
                        ipeps_seconds,
                        reversibility_seconds,
                    ),
                ),
            )
        end
    end

    return PXPLargerDBenchmarkReport(config, runs, _validation_metadata())
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
        exact_finite_observables = config.exact_finite_observables,
        exact_finite_max_sites = config.exact_finite_max_sites,
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
        sublattice_imbalance = summary.sublattice_imbalance,
        checkerboard_structure_factor = summary.checkerboard_structure_factor,
        diagnostics = _ctm_diagnostics_data(summary.diagnostics),
    )
end

function _ctm_point_data(point::CTMValidationPoint)
    return (;
        chi = point.params.chi,
        tol = point.params.tol,
        maxiter = point.params.maxiter,
        verbosity = point.params.verbosity,
        seed = point.params.seed,
        measurement = _ctm_summary_data(point.measurement),
        delta_density = point.delta_density,
        delta_density_even = point.delta_density_even,
        delta_density_odd = point.delta_density_odd,
        delta_blockade_violation = point.delta_blockade_violation,
        delta_pxp_energy_density = point.delta_pxp_energy_density,
    )
end

function _trust_policy_data(policy::CTMTrustPolicy)
    return (;
        min_points = policy.min_points,
        require_accepted_diagnostics = policy.require_accepted_diagnostics,
        max_density_delta = policy.max_density_delta,
        max_blockade_delta = policy.max_blockade_delta,
        max_energy_delta = policy.max_energy_delta,
        max_residual = policy.max_residual,
    )
end

function _trust_data(trust::CTMTrustAssessment, policy::CTMTrustPolicy)
    return (;
        trusted = trust.trusted,
        reason = String(trust.reason),
        message = trust.message,
        compared_points = trust.compared_points,
        finite_chi_density_delta = trust.finite_chi_density_delta,
        finite_chi_blockade_delta = trust.finite_chi_blockade_delta,
        finite_chi_energy_delta = trust.finite_chi_energy_delta,
        observed_max_residual = trust.observed_max_residual,
        policy = _trust_policy_data(policy),
    )
end

function _trusted_ctm_data(ctm::Nothing)
    return nothing
end

function _trusted_ctm_data(ctm::TrustedCTMMeasurement)
    return (;
        measurement = _ctm_summary_data(ctm.measurement),
        points = _ctm_point_data.(ctm.points),
        trust = _trust_data(ctm.trust, ctm.policy),
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
        exact_finite_density = sample.exact_finite_density,
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
        ipeps_exact_finite_density = sample.ipeps_exact_finite_density,
        density_error_simple = sample.density_error_simple,
        density_error_ctm = sample.density_error_ctm,
        density_error_exact_finite = sample.density_error_exact_finite,
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

function _convergence_config_data(config::PXPConvergenceConfig)
    return (;
        base = _config_data(config.base),
        dt_values = config.dt_values,
        D_values = config.D_values,
        chi_values = config.chi_values,
        cutoff_values = config.cutoff_values,
    )
end

function _convergence_report_data(report::PXPConvergenceReport)
    return (;
        config = _convergence_config_data(report.config),
        summary = (;
            max_abs_density_error_simple = report.max_abs_density_error_simple,
            max_abs_density_error_ctm = report.max_abs_density_error_ctm,
            max_abs_density_error_exact_finite = report.max_abs_density_error_exact_finite,
            all_ctm_trusted = report.all_ctm_trusted,
        ),
        runs = _report_data.(report.runs),
    )
end

function _audit_config_data(config::PXPAuditConfig)
    return (;
        n_values = config.n_values,
        total_time = config.total_time,
        dt_values = config.dt_values,
        D_values = config.D_values,
        cutoff_values = config.cutoff_values,
        chi_values = config.chi_values,
        measure_every = config.measure_every,
        order = config.order,
        schedule = String(config.schedule),
        ctm_tol = config.ctm_tol,
        ctm_maxiter = config.ctm_maxiter,
        ctm_verbosity = config.ctm_verbosity,
        ctm_seed = config.ctm_seed,
        exact_finite_observables = config.exact_finite_observables,
        exact_finite_max_sites = config.exact_finite_max_sites,
    )
end

function _reversibility_data(report::PXPReversibilityReport)
    return (;
        before = _simple_data(report.before),
        after_forward = _simple_data(report.after_forward),
        after_reverse = _simple_data(report.after_reverse),
        forward_log = _evolution_data(report.forward_log),
        reverse_log = _evolution_data(report.reverse_log),
        density_drift = report.density_drift,
        blockade_drift = report.blockade_drift,
        energy_drift = report.energy_drift,
    )
end

function _audit_summary_data(summary::PXPAuditSummary)
    return (;
        n = summary.n,
        total_time = summary.total_time,
        dt = summary.dt,
        D = summary.D,
        cutoff = summary.cutoff,
        schedule = String(summary.schedule),
        chi_values = summary.chi_values,
        max_abs_density_error_simple = summary.max_abs_density_error_simple,
        max_abs_density_error_ctm = summary.max_abs_density_error_ctm,
        max_abs_density_error_exact_finite = summary.max_abs_density_error_exact_finite,
        max_blockade_violation_simple = summary.max_blockade_violation_simple,
        max_blockade_violation_ctm = summary.max_blockade_violation_ctm,
        pxp_energy_drift_simple = summary.pxp_energy_drift_simple,
        pxp_energy_drift_ctm = summary.pxp_energy_drift_ctm,
        ctm_trust_status = String(summary.ctm_trust_status),
        ctm_trust_reason = String(summary.ctm_trust_reason),
        finite_chi_density_delta = summary.finite_chi_density_delta,
        finite_chi_blockade_delta = summary.finite_chi_blockade_delta,
        finite_chi_energy_delta = summary.finite_chi_energy_delta,
        finite_chi_max_residual = summary.finite_chi_max_residual,
        max_truncerr = summary.max_truncerr,
        log_norm_initial = summary.log_norm_initial,
        log_norm_final = summary.log_norm_final,
        log_norm_delta = summary.log_norm_delta,
        log_norm_delta_abs = summary.log_norm_delta_abs,
        reversibility_density_drift = summary.reversibility_density_drift,
        reversibility_blockade_drift = summary.reversibility_blockade_drift,
        reversibility_energy_drift = summary.reversibility_energy_drift,
    )
end

function _audit_run_data(run::PXPAuditRun)
    return (;
        summary = _audit_summary_data(run.summary),
        validation = _report_data(run.validation),
        reversibility = _reversibility_data(run.reversibility),
    )
end

function _audit_report_data(report::PXPAuditReport)
    return (;
        config = _audit_config_data(report.config),
        metadata = _metadata_data(report.metadata),
        runs = _audit_run_data.(report.runs),
    )
end

function _larger_d_config_data(config::PXPLargerDBenchmarkConfig)
    return (;
        n_values = config.n_values,
        total_time = config.total_time,
        dt_values = config.dt_values,
        D_values = config.D_values,
        cutoff_values = config.cutoff_values,
        measure_every = config.measure_every,
        order = config.order,
        schedule = String(config.schedule),
        initial_state = String(config.initial_state),
        point_group = config.point_group,
        use_sparse = config.use_sparse,
        ed_tol = config.ed_tol,
        ed_m_init = config.ed_m_init,
        ed_m_max = config.ed_m_max,
        ed_extend_step = config.ed_extend_step,
        ed_mode = String(config.ed_mode),
        observable_mode = String(config.observable_mode),
        chi_values = config.chi_values,
        ctm_tol = config.ctm_tol,
        ctm_maxiter = config.ctm_maxiter,
        ctm_verbosity = config.ctm_verbosity,
        ctm_seed = config.ctm_seed,
        exact_finite_observables = config.exact_finite_observables,
        exact_finite_max_sites = config.exact_finite_max_sites,
    )
end

function _larger_d_summary_data(summary::PXPLargerDBenchmarkSummary)
    return (;
        n = summary.n,
        D = summary.D,
        dt = summary.dt,
        cutoff = summary.cutoff,
        total_time = summary.total_time,
        ed_mode = String(summary.ed_mode),
        observable_mode = String(summary.observable_mode),
        ed_boundary_condition = String(summary.ed_boundary_condition),
        ed_symmetry_sector = String(summary.ed_symmetry_sector),
        ed_observable_scope = String(summary.ed_observable_scope),
        ed_reference_label = summary.ed_reference_label,
        ed_basis_dimension = summary.ed_basis_dimension,
        ed_constrained_dimension = summary.ed_constrained_dimension,
        ed_group_order = summary.ed_group_order,
        ed_hamiltonian_nnz = summary.ed_hamiltonian_nnz,
        ed_runtime_seconds = summary.ed_runtime_seconds,
        ipeps_runtime_seconds = summary.ipeps_runtime_seconds,
        reversibility_runtime_seconds = summary.reversibility_runtime_seconds,
        density_error_simple = summary.density_error_simple,
        density_error_exact_finite = summary.density_error_exact_finite,
        density_error_ctm = summary.density_error_ctm,
        return_probability_error = summary.return_probability_error,
        ed_return_probability = summary.ed_return_probability,
        ed_excitation_density = summary.ed_excitation_density,
        ipeps_simple_density = summary.ipeps_simple_density,
        ipeps_exact_finite_density = summary.ipeps_exact_finite_density,
        ipeps_ctm_density = summary.ipeps_ctm_density,
        max_truncerr = summary.max_truncerr,
        log_norm_initial = summary.log_norm_initial,
        log_norm_final = summary.log_norm_final,
        log_norm_delta_abs = summary.log_norm_delta_abs,
        reversibility_density_drift = summary.reversibility_density_drift,
        ctm_trust_status = String(summary.ctm_trust_status),
        ctm_trust_reason = String(summary.ctm_trust_reason),
        notes = summary.notes,
        warnings = summary.warnings,
    )
end

function _larger_d_run_data(run::PXPLargerDBenchmarkRun)
    return (;
        summary = _larger_d_summary_data(run.summary),
        validation = _report_data(run.validation),
        reversibility = _reversibility_data(run.reversibility),
    )
end

function _larger_d_report_data(report::PXPLargerDBenchmarkReport)
    return (;
        schema_version = 1,
        config = _larger_d_config_data(report.config),
        metadata = _metadata_data(report.metadata),
        runs = _larger_d_run_data.(report.runs),
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

"""
    write_pxp_convergence_json(report, path)

Write a [`PXPConvergenceReport`](@ref) to `path` as JSON containing the swept
configuration, aggregate error-budget summary fields, and full per-run
validation reports.
"""
function write_pxp_convergence_json(report::PXPConvergenceReport, path::AbstractString)
    open(path, "w") do io
        JSON3.write(io, _convergence_report_data(report))
        write(io, '\n')
    end
    return path
end

"""
    write_pxp_audit_json(report, path)

Write a [`PXPAuditReport`](@ref) as nested JSON containing campaign
configuration, per-run flat summaries, full validation reports, and
reversibility reports.
"""
function write_pxp_audit_json(report::PXPAuditReport, path::AbstractString)
    open(path, "w") do io
        JSON3.write(io, _audit_report_data(report))
        write(io, '\n')
    end
    return path
end

"""
    write_pxp_larger_d_benchmark_json(report, path)

Write a nested M3 larger-D PXP ED benchmark report as JSON.
"""
function write_pxp_larger_d_benchmark_json(
    report::PXPLargerDBenchmarkReport,
    path::AbstractString,
)
    open(path, "w") do io
        JSON3.write(io, _larger_d_report_data(report))
        write(io, '\n')
    end
    return path
end

const PXP_AUDIT_CSV_HEADER = [
    "n",
    "total_time",
    "dt",
    "D",
    "cutoff",
    "schedule",
    "chi_values",
    "max_abs_density_error_simple",
    "max_abs_density_error_ctm",
    "max_abs_density_error_exact_finite",
    "max_blockade_violation_simple",
    "max_blockade_violation_ctm",
    "pxp_energy_drift_simple",
    "pxp_energy_drift_ctm",
    "ctm_trust_status",
    "ctm_trust_reason",
    "finite_chi_density_delta",
    "finite_chi_blockade_delta",
    "finite_chi_energy_delta",
    "finite_chi_max_residual",
    "max_truncerr",
    "log_norm_initial",
    "log_norm_final",
    "log_norm_delta",
    "log_norm_delta_abs",
    "reversibility_density_drift",
    "reversibility_blockade_drift",
    "reversibility_energy_drift",
]

function _audit_csv_cell(x::Real)
    isfinite(Float64(x)) || throw(ArgumentError("audit CSV values must be finite"))
    return string(x)
end

_audit_csv_cell(x::Symbol) = String(x)
_audit_csv_cell(::Nothing) = ""

function _audit_csv_cell(xs::Vector{Int})
    return isempty(xs) ? "" : join(xs, ";")
end

function _audit_csv_cell(x::AbstractString)
    if occursin(r"[,\n\"]", x)
        return "\"" * replace(x, "\"" => "\"\"") * "\""
    else
        return x
    end
end

function _audit_csv_row(summary::PXPAuditSummary)
    values = (
        summary.n,
        summary.total_time,
        summary.dt,
        summary.D,
        summary.cutoff,
        summary.schedule,
        summary.chi_values,
        summary.max_abs_density_error_simple,
        summary.max_abs_density_error_ctm,
        summary.max_abs_density_error_exact_finite,
        summary.max_blockade_violation_simple,
        summary.max_blockade_violation_ctm,
        summary.pxp_energy_drift_simple,
        summary.pxp_energy_drift_ctm,
        summary.ctm_trust_status,
        summary.ctm_trust_reason,
        summary.finite_chi_density_delta,
        summary.finite_chi_blockade_delta,
        summary.finite_chi_energy_delta,
        summary.finite_chi_max_residual,
        summary.max_truncerr,
        summary.log_norm_initial,
        summary.log_norm_final,
        summary.log_norm_delta,
        summary.log_norm_delta_abs,
        summary.reversibility_density_drift,
        summary.reversibility_blockade_drift,
        summary.reversibility_energy_drift,
    )
    return join(_audit_csv_cell.(values), ",")
end

"""
    write_pxp_audit_csv(report, path)

Write the flat [`PXPAuditSummary`](@ref) rows from a [`PXPAuditReport`](@ref)
as CSV. Nested validation, CTM, and reversibility details remain available in
the JSON artifact.
"""
function write_pxp_audit_csv(report::PXPAuditReport, path::AbstractString)
    open(path, "w") do io
        println(io, join(PXP_AUDIT_CSV_HEADER, ","))
        for run in report.runs
            println(io, _audit_csv_row(run.summary))
        end
    end
    return path
end

const PXP_LARGER_D_CSV_HEADER = [
    "n",
    "D",
    "dt",
    "cutoff",
    "total_time",
    "ed_mode",
    "observable_mode",
    "ed_boundary_condition",
    "ed_symmetry_sector",
    "ed_observable_scope",
    "ed_reference_label",
    "ed_basis_dimension",
    "ed_constrained_dimension",
    "ed_group_order",
    "ed_hamiltonian_nnz",
    "ed_runtime_seconds",
    "ipeps_runtime_seconds",
    "reversibility_runtime_seconds",
    "density_error_simple",
    "density_error_exact_finite",
    "density_error_ctm",
    "return_probability_error",
    "ed_return_probability",
    "ed_excitation_density",
    "ipeps_simple_density",
    "ipeps_exact_finite_density",
    "ipeps_ctm_density",
    "max_truncerr",
    "log_norm_initial",
    "log_norm_final",
    "log_norm_delta_abs",
    "reversibility_density_drift",
    "ctm_trust_status",
    "ctm_trust_reason",
    "notes",
    "warnings",
]

function _larger_d_csv_cell(xs::Vector{String})
    return _audit_csv_cell(join(xs, ";"))
end

_larger_d_csv_cell(x) = _audit_csv_cell(x)

function _larger_d_csv_row(summary::PXPLargerDBenchmarkSummary)
    values = (
        summary.n,
        summary.D,
        summary.dt,
        summary.cutoff,
        summary.total_time,
        summary.ed_mode,
        summary.observable_mode,
        summary.ed_boundary_condition,
        summary.ed_symmetry_sector,
        summary.ed_observable_scope,
        summary.ed_reference_label,
        summary.ed_basis_dimension,
        summary.ed_constrained_dimension,
        summary.ed_group_order,
        summary.ed_hamiltonian_nnz,
        summary.ed_runtime_seconds,
        summary.ipeps_runtime_seconds,
        summary.reversibility_runtime_seconds,
        summary.density_error_simple,
        summary.density_error_exact_finite,
        summary.density_error_ctm,
        summary.return_probability_error,
        summary.ed_return_probability,
        summary.ed_excitation_density,
        summary.ipeps_simple_density,
        summary.ipeps_exact_finite_density,
        summary.ipeps_ctm_density,
        summary.max_truncerr,
        summary.log_norm_initial,
        summary.log_norm_final,
        summary.log_norm_delta_abs,
        summary.reversibility_density_drift,
        summary.ctm_trust_status,
        summary.ctm_trust_reason,
        summary.notes,
        summary.warnings,
    )
    return join(_larger_d_csv_cell.(values), ",")
end

"""
    write_pxp_larger_d_benchmark_csv(report, path)

Write the flat M3 benchmark summary rows as CSV.
"""
function write_pxp_larger_d_benchmark_csv(
    report::PXPLargerDBenchmarkReport,
    path::AbstractString,
)
    open(path, "w") do io
        println(io, join(PXP_LARGER_D_CSV_HEADER, ","))
        for run in report.runs
            println(io, _larger_d_csv_row(run.summary))
        end
    end
    return path
end

end
