module PXPValidation

import EDKit

using JSON3

using ..SquareIPEPS: SquareIPEPSState, copy_state, log_norm, product_square_ipeps
using ..SquareUnitCells: PeriodicSquareUnitCell
using ..PEPSKitBackend: product_pxp_pepskit_state, PXPIPEPSState
using ..Observables: SimpleObservableSummary, measure_simple
using ..FiniteIPEPSObservables: exact_density_finite, exact_all_down_return_probability_finite
using ..PEPSKitMeasurements:
    CTMObservableSummary,
    CTMRGDiagnostics,
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
using ..Internals:
    _finite_positive,
    _finite_nonnegative,
    _positive_int,
    _nonnegative_int,
    _finite_max_or_nothing

export TrustedCTMMeasurement, measure_ctm_trusted
export PXPValidationConfig, PXPValidationMetadata, PXPIPEPSSample
export PXPEDComparisonSample, PXPValidationReport, validate_pxp_ed_ipeps
export write_pxp_validation_json
export PXPConvergenceConfig, PXPConvergenceReport, validate_pxp_convergence
export write_pxp_convergence_json
export PXPReversibilityReport, validate_pxp_reversibility
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

"""
    PXPValidationConfig(n; kwargs...)

Controls for a short-time finite-ED versus iPEPS PXP validation run. The ED
side uses [`PXPEEDBenchmarkConfig`](@ref); the iPEPS side uses a periodic
`n x n` unit cell, all-down product initialization, real-time PXP evolution,
and a serial star schedule by default so `3 x 3` smoke validation is supported.

`tensor_engine` selects the iPEPS backend: `:legacy_itensors` (default) or
`:pepskit_simple` (the PEPSKit-native simple-update engine, dispatched through
the shared `evolve!`/`measure_simple` generics). The native engine does not yet
support `exact_finite_observables` or CTM trust sweeps; both are rejected.
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
    tensor_engine::Symbol

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
        tensor_engine::Symbol = :legacy_itensors,
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
        tensor_engine in (:legacy_itensors, :pepskit_simple) || throw(
            ArgumentError("tensor_engine must be :legacy_itensors or :pepskit_simple"),
        )
        # The PEPSKit-native engine has no ITensor exact-finite contraction yet;
        # CTM trust sweeps are likewise legacy-only (rejected at the call site).
        tensor_engine === :pepskit_simple && exact_finite_observables && throw(
            ArgumentError(
                "exact_finite_observables is not supported on the :pepskit_simple engine",
            ),
        )

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
            tensor_engine,
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
        tensor_engine = base.tensor_engine,
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
    psi::Union{SquareIPEPSState,PXPIPEPSState},
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
    if config.tensor_engine === :pepskit_simple
        return product_pxp_pepskit_state(cell; state, D = config.maxdim)
    end
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
    psi::Union{SquareIPEPSState,PXPIPEPSState},
    ctm_params,
    trust_policy::CTMTrustPolicy,
    ctm_measure,
)
    # The native engine reaches here only with ctm_params === nothing (CTM sweeps
    # are rejected upstream), so the legacy-typed measure_ctm_trusted is never hit.
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
    config.tensor_engine === :pepskit_simple && ctm_params !== nothing && throw(
        ArgumentError(
            "CTM trust sweeps are not supported on the :pepskit_simple engine; pass ctm_params = nothing",
        ),
    )
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


function _ctm_param_tuple(
    chi_values,
    ctm_tol::Real,
    ctm_maxiter::Integer,
    ctm_verbosity::Integer,
    ctm_seed,
)
    isempty(chi_values) && return nothing
    return Tuple(
        default_ctmrg_params(;
            chi = chi,
            tol = ctm_tol,
            maxiter = ctm_maxiter,
            verbosity = ctm_verbosity,
            seed = ctm_seed,
        ) for chi in chi_values
    )
end

function _maximum_or_zero(values)
    isempty(values) && return 0.0
    return maximum(values)
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

_larger_d_ctm_params(config::PXPLargerDBenchmarkConfig) = _ctm_param_tuple(
    config.chi_values,
    config.ctm_tol,
    config.ctm_maxiter,
    config.ctm_verbosity,
    config.ctm_seed,
)

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

# ---------------------------------------------------------------------------
# Serialization: StructTypes/JSON3/CSV writers for the validation reports.
# Merged from former PXPValidationSerialization.jl.

import StructTypes

# StructTypes declarations let JSON3 walk these types directly without
# hand-written `_*_data` NamedTuple builders. Types that need projection or
# flattening (EvolutionLog drops layer_infos, CTMValidationPoint inlines params,
# TrustedCTMMeasurement combines trust + policy, the report wrappers
# synthesize summary/schema_version fields) get a `StructTypes.lower` instead.

for T in (
    PXPValidationConfig,
    PXPValidationMetadata,
    SimpleObservableSummary,
    CTMObservableSummary,
    CTMRGDiagnostics,
    CTMTrustPolicy,
    CTMTrustAssessment,
    PXPEEDSample,
    PXPEEDBenchmarkResult,
    PXPIPEPSSample,
    PXPEDComparisonSample,
    PXPValidationReport,
    PXPConvergenceConfig,
    PXPReversibilityReport,
    PXPLargerDBenchmarkConfig,
    PXPLargerDBenchmarkSummary,
    PXPLargerDBenchmarkRun,
)
    @eval StructTypes.StructType(::Type{$T}) = StructTypes.Struct()
end

# Project EDKit's KrylovEvolutionDiagnostics down to the seven public fields the
# audit/validation reports care about. Avoids leaking unrelated internal state.
StructTypes.StructType(::Type{<:EDKit.KrylovEvolutionDiagnostics}) = StructTypes.CustomStruct()
StructTypes.lower(d::EDKit.KrylovEvolutionDiagnostics) = (;
    basis_builds = d.basis_builds,
    basis_extensions = d.basis_extensions,
    restarts = d.restarts,
    matvecs = d.matvecs,
    total_times_served = d.total_times_served,
    max_dim_used = d.max_dim_used,
    accepted_intervals = d.accepted_intervals,
)

StructTypes.StructType(::Type{EvolutionLog}) = StructTypes.CustomStruct()
StructTypes.lower(e::EvolutionLog) = (;
    total_time = e.total_time,
    nsteps = e.nsteps,
    dt = e.params.dt,
    order = e.params.order,
    evolution = e.params.evolution,
    maxdim = e.params.maxdim,
    cutoff = e.params.cutoff,
    schedule = e.params.schedule,
    max_truncerr = e.max_truncerr,
    max_bond_entropy = e.max_bond_entropy,
    mean_bond_entropy = e.mean_bond_entropy,
    log_norm_before = e.log_norm_before,
    log_norm_after = e.log_norm_after,
    log_norm_delta = e.log_norm_delta,
    model_metadata = e.model_metadata,
)

StructTypes.StructType(::Type{CTMValidationPoint}) = StructTypes.CustomStruct()
StructTypes.lower(p::CTMValidationPoint) = (;
    chi = p.params.chi,
    tol = p.params.tol,
    maxiter = p.params.maxiter,
    verbosity = p.params.verbosity,
    seed = p.params.seed,
    measurement = p.measurement,
    delta_density = p.delta_density,
    delta_density_even = p.delta_density_even,
    delta_density_odd = p.delta_density_odd,
    delta_blockade_violation = p.delta_blockade_violation,
    delta_pxp_energy_density = p.delta_pxp_energy_density,
)

StructTypes.StructType(::Type{TrustedCTMMeasurement}) = StructTypes.CustomStruct()
StructTypes.lower(t::TrustedCTMMeasurement) = (;
    measurement = t.measurement,
    points = t.points,
    trust = (;
        trusted = t.trust.trusted,
        reason = t.trust.reason,
        message = t.trust.message,
        compared_points = t.trust.compared_points,
        finite_chi_density_delta = t.trust.finite_chi_density_delta,
        finite_chi_blockade_delta = t.trust.finite_chi_blockade_delta,
        finite_chi_energy_delta = t.trust.finite_chi_energy_delta,
        observed_max_residual = t.trust.observed_max_residual,
        policy = t.policy,
    ),
)

StructTypes.StructType(::Type{PXPConvergenceReport}) = StructTypes.CustomStruct()
StructTypes.lower(r::PXPConvergenceReport) = (;
    config = r.config,
    summary = (;
        max_abs_density_error_simple = r.max_abs_density_error_simple,
        max_abs_density_error_ctm = r.max_abs_density_error_ctm,
        max_abs_density_error_exact_finite = r.max_abs_density_error_exact_finite,
        all_ctm_trusted = r.all_ctm_trusted,
    ),
    runs = r.runs,
)

StructTypes.StructType(::Type{PXPLargerDBenchmarkReport}) = StructTypes.CustomStruct()
StructTypes.lower(r::PXPLargerDBenchmarkReport) =
    (; schema_version = 1, config = r.config, metadata = r.metadata, runs = r.runs)

function _write_json(path::AbstractString, value)
    open(path, "w") do io
        JSON3.write(io, value)
        write(io, '\n')
    end
    return path
end

"""
    write_pxp_validation_json(report, path)

Write a [`PXPValidationReport`](@ref) to `path` as a JSON artifact containing
configuration, metadata, ED samples, iPEPS samples, CTM trust data when present,
and matched observable comparisons.
"""
write_pxp_validation_json(report::PXPValidationReport, path::AbstractString) =
    _write_json(path, report)

"""
    write_pxp_convergence_json(report, path)

Write a [`PXPConvergenceReport`](@ref) to `path` as JSON containing the swept
configuration, aggregate error-budget summary fields, and full per-run
validation reports.
"""
write_pxp_convergence_json(report::PXPConvergenceReport, path::AbstractString) =
    _write_json(path, report)

"""
    write_pxp_larger_d_benchmark_json(report, path)

Write a nested M3 larger-D PXP ED benchmark report as JSON.
"""
write_pxp_larger_d_benchmark_json(
    report::PXPLargerDBenchmarkReport,
    path::AbstractString,
) = _write_json(path, report)

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
