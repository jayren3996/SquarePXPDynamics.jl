#!/usr/bin/env julia
#
# ScarFinder audit campaign driver.
#
# Runs `run_scarfinder_audit` on a user-specified grid and emits the full
# JSON+CSV report plus a one-line summary suitable for grepping. The objective
# is chosen via SQUAREPXP_AUDIT_OBJECTIVE; the trust policy is one of
# `default`, `tight`, or `calibrated` (calibrated requires a prior sensitivity
# sweep JSON specified via SQUAREPXP_AUDIT_CALIBRATION_JSON).
#
# Phase 5 of the ScarFinder reliability plan. See
# docs/superpowers/notes/2026-05-26-scarfinder-acceptance.md.

using Pkg

project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using JSON3
using Printf
using SquarePXPDynamics

function _env_value(name, default)
    raw = get(ENV, String(name), "")
    return isempty(strip(raw)) ? String(default) : String(strip(raw))
end

_env_int(name, default) = parse(Int, _env_value(name, string(default)))
_env_float(name, default) = parse(Float64, _env_value(name, string(default)))
_env_symbol(name, default) = Symbol(_env_value(name, String(default)))

function _env_bool(name, default)
    v = lowercase(strip(_env_value(name, string(default))))
    v in ("1", "true", "yes", "on") && return true
    v in ("0", "false", "no", "off") && return false
    throw(ArgumentError("$name must be boolean"))
end

function _env_int_list(name, default)
    raw = strip(_env_value(name, default))
    isempty(raw) && return Int[]
    return parse.(Int, split(raw, ","))
end

function _env_float_list(name, default)
    raw = strip(_env_value(name, default))
    isempty(raw) && return Float64[]
    return parse.(Float64, split(raw, ","))
end

function _select_objective(kind::Symbol)
    if kind === :revival_density
        return RevivalObjective(:density)
    elseif kind === :revival_imbalance
        return RevivalObjective(:sublattice_imbalance)
    elseif kind === :composite
        return CompositeObjective()
    else
        throw(ArgumentError("unknown SQUAREPXP_AUDIT_OBJECTIVE $kind"))
    end
end

function _select_trust_policy(kind::Symbol, calibration_json::String)
    if kind === :default
        return CTMTrustPolicy()
    elseif kind === :tight
        return tight_ctm_trust_policy()
    elseif kind === :calibrated
        isempty(calibration_json) && throw(
            ArgumentError(
                "trust=calibrated requires SQUAREPXP_AUDIT_CALIBRATION_JSON",
            ),
        )
        data = JSON3.read(read(calibration_json, String))
        c = data["calibrated_policy"]
        return CTMTrustPolicy(
            Int(c["min_points"]),
            Bool(c["require_accepted_diagnostics"]),
            Float64(c["max_density_delta"]),
            Float64(c["max_blockade_delta"]),
            Float64(c["max_energy_delta"]),
            c["max_residual"] === nothing ? nothing : Float64(c["max_residual"]),
        )
    else
        throw(ArgumentError("unknown SQUAREPXP_AUDIT_TRUST $kind"))
    end
end

function main()
    configure_ctm_threading_from_env!()

    cell_Lx = _env_int("SQUAREPXP_AUDIT_CELL_LX", 3)
    cell_Ly = _env_int("SQUAREPXP_AUDIT_CELL_LY", 3)
    initial_state = _env_symbol("SQUAREPXP_AUDIT_INITIAL_STATE", :z_up)
    projection_time = _env_float("SQUAREPXP_AUDIT_PROJECTION_TIME", 0.08)
    iterations = _env_int("SQUAREPXP_AUDIT_ITERATIONS", 4)
    dt_values = _env_float_list("SQUAREPXP_AUDIT_DT", "0.02")
    evolve_maxdim = _env_int("SQUAREPXP_AUDIT_EVOLVE_D", 4)
    target_maxdim_values = _env_int_list("SQUAREPXP_AUDIT_TARGET_D", "2")
    compression_cutoff = _env_float("SQUAREPXP_AUDIT_COMPRESSION_CUTOFF", 1e-12)
    chi_values = _env_int_list("SQUAREPXP_AUDIT_CHI", "")
    cutoff_values = _env_float_list("SQUAREPXP_AUDIT_CUTOFF", "1e-12")
    schedule = _env_symbol("SQUAREPXP_AUDIT_SCHEDULE", :serial)
    order = _env_int("SQUAREPXP_AUDIT_ORDER", 1)
    blockade_tol = _env_float("SQUAREPXP_AUDIT_BLOCKADE_TOL", Inf)
    truncation_tol = _env_float("SQUAREPXP_AUDIT_TRUNCATION_TOL", Inf)
    entropy_cap = _env_float("SQUAREPXP_AUDIT_ENTROPY_CAP", Inf)
    ctm_tol = _env_float("SQUAREPXP_AUDIT_CTM_TOL", 1e-8)
    ctm_maxiter = _env_int("SQUAREPXP_AUDIT_CTM_MAXITER", 100)
    ctm_seed = _env_int("SQUAREPXP_AUDIT_CTM_SEED", 0)
    require_trusted_ctm = _env_bool("SQUAREPXP_AUDIT_REQUIRE_TRUSTED_CTM", false)
    top_k = _env_int("SQUAREPXP_AUDIT_TOP_K", 3)

    objective_kind = _env_symbol("SQUAREPXP_AUDIT_OBJECTIVE", :revival_density)
    trust_kind = _env_symbol("SQUAREPXP_AUDIT_TRUST", :default)
    calibration_json = _env_value("SQUAREPXP_AUDIT_CALIBRATION_JSON", "")
    label = _env_value("SQUAREPXP_AUDIT_LABEL", "v1")

    csv_out = _env_value("SQUAREPXP_AUDIT_OUTPUT_CSV",
        joinpath(project_root, "artifacts", "scarfinder_audit_$(label).csv"))
    json_out = _env_value("SQUAREPXP_AUDIT_OUTPUT_JSON",
        joinpath(project_root, "artifacts", "scarfinder_audit_$(label).json"))

    objective = _select_objective(objective_kind)
    trust_policy = _select_trust_policy(trust_kind, calibration_json)

    config = ScarFinderAuditConfig(;
        cell_Lx,
        cell_Ly,
        initial_state,
        projection_time,
        iterations,
        dt_values,
        evolve_maxdim,
        target_maxdim_values,
        cutoff_values,
        compression_cutoff,
        chi_values,
        schedule,
        order,
        blockade_tol,
        truncation_tol,
        entropy_cap,
        ctm_tol,
        ctm_maxiter,
        ctm_seed,
        ctm_trust_policy = trust_policy,
        require_trusted_ctm,
        top_k,
    )

    println("Running ScarFinder audit (label=$label)")
    total_evolution_time = projection_time * iterations
    println("  cell=$(cell_Lx)x$(cell_Ly) initial=$initial_state projection_time=$projection_time iterations=$iterations (total evolution=$total_evolution_time)")
    println("  dt=$dt_values  evolve_maxdim=$evolve_maxdim  target_maxdim=$target_maxdim_values  chi=$chi_values  cutoff=$cutoff_values")
    println("  objective=$objective_kind  trust_policy=$trust_kind")
    println("  julia_threads=$(Threads.nthreads())")

    elapsed = @elapsed report = run_scarfinder_audit(config; objective)

    write_scarfinder_audit_json(report, json_out)
    write_scarfinder_audit_csv(report, csv_out)

    println("--- per-row summary ---")
    for row in report.rows
        chi_str = row.chi === nothing ? "-" : string(row.chi)
        best_iter = row.best_iteration === nothing ? "-" : string(row.best_iteration)
        best_score = row.best_score === nothing ? "-" : @sprintf("%.6e", row.best_score)
        comp_trunc = row.max_compression_truncerr === nothing ? "-" : @sprintf("%.2e", row.max_compression_truncerr)
        @printf("  dt=%g target_D=%d evolve_D=%d chi=%s backend=%s best_iter=%s best_score=%s comp_trunc=%s ok=%s  %.2fs\n",
            row.dt, row.target_maxdim, row.evolve_maxdim, chi_str, row.backend, best_iter,
            best_score, comp_trunc, string(row.ok), row.elapsed_seconds)
    end
    stab = report.stability
    function _fmt(value, fmt)
        return value === nothing ? "-" : Printf.format(Printf.Format(fmt), value)
    end
    @printf("Stability counts: n=%d successful=%d (simple=%d ctm=%d)\n",
        stab.n_rows, stab.n_successful, stab.n_simple, stab.n_ctm)
    @printf("Stability all:    best_iter_agreement=%s score_range=%s score_cv=%s top_k_jaccard_min=%s\n",
        _fmt(stab.best_iteration_agreement_all, "%.3f"),
        _fmt(stab.score_range_all, "%.3e"),
        _fmt(stab.score_cv_all, "%.3e"),
        _fmt(stab.top_k_jaccard_min_all, "%.3f"))
    @printf("Stability simple: best_iter_agreement=%s score_range=%s score_cv=%s top_k_jaccard_min=%s\n",
        _fmt(stab.best_iteration_agreement_simple, "%.3f"),
        _fmt(stab.score_range_simple, "%.3e"),
        _fmt(stab.score_cv_simple, "%.3e"),
        _fmt(stab.top_k_jaccard_min_simple, "%.3f"))
    @printf("Stability ctm:    best_iter_agreement=%s score_range=%s score_cv=%s top_k_jaccard_min=%s\n",
        _fmt(stab.best_iteration_agreement_ctm, "%.3f"),
        _fmt(stab.score_range_ctm, "%.3e"),
        _fmt(stab.score_cv_ctm, "%.3e"),
        _fmt(stab.top_k_jaccard_min_ctm, "%.3f"))
    @printf("Total wall time: %.2fs\n", elapsed)
    println("Wrote $json_out")
    println("Wrote $csv_out")
end

main()
