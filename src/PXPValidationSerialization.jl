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
