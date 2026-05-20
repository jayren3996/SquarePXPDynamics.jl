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
    PXPAuditConfig,
    PXPAuditSummary,
    PXPAuditRun,
    PXPAuditReport,
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
    write_pxp_audit_json(report, path)

Write a [`PXPAuditReport`](@ref) as nested JSON containing campaign
configuration, per-run flat summaries, full validation reports, and
reversibility reports.
"""
write_pxp_audit_json(report::PXPAuditReport, path::AbstractString) =
    _write_json(path, report)

"""
    write_pxp_larger_d_benchmark_json(report, path)

Write a nested M3 larger-D PXP ED benchmark report as JSON.
"""
write_pxp_larger_d_benchmark_json(
    report::PXPLargerDBenchmarkReport,
    path::AbstractString,
) = _write_json(path, report)

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
