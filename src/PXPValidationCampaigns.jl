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

_audit_ctm_params(config::PXPAuditConfig) = _ctm_param_tuple(
    config.chi_values,
    config.ctm_tol,
    config.ctm_maxiter,
    config.ctm_verbosity,
    config.ctm_seed,
)

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
