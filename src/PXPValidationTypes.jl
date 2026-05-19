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
