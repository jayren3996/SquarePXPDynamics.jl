#!/usr/bin/env julia
#
# CTM finite-chi sensitivity sweep.
#
# Build a representative iPEPS state (product :down + optional short PXP
# real-time evolution at given D), call `ctm_chi_sensitivity_sweep` over an
# increasing chi grid, then run `assess_ctm_trust` and also derive a calibrated
# trust policy from the observed adjacent drifts. Writes a CSV (matching
# `write_ctm_trust_csv`) plus a JSON containing the per-chi measurements, the
# assessment, and the calibrated thresholds.
#
# Use this to seed concrete CTMTrustPolicy numerical thresholds (Phase 2 of
# the ScarFinder reliability plan). See
# docs/superpowers/notes/2026-05-26-ctm-trust-calibration.md.

using Pkg

project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using JSON3
using Printf
using SquarePXPDynamics

configure_ctm_threading_from_env!()

function _env_value(name::AbstractString, default::AbstractString)
    raw = get(ENV, String(name), "")
    return isempty(strip(raw)) ? String(default) : String(strip(raw))
end

_env_int(name, default) = parse(Int, _env_value(name, string(default)))
_env_float(name, default) = parse(Float64, _env_value(name, string(default)))

function _env_int_list(name, default)
    raw = strip(_env_value(name, default))
    isempty(raw) && return Int[]
    return parse.(Int, split(raw, ","))
end

function _build_state(Lx, Ly, D, dt, evolve_time)
    cell = PeriodicSquareUnitCell(Lx, Ly)
    psi = product_square_ipeps(cell; state = :down, maxdim = D)
    if evolve_time > 0 && D > 1
        params = TrotterParams(dt, 1, :real, D, 1e-12; schedule = :serial)
        evolve!(psi, evolve_time; params)
    end
    return psi
end

function _measurement_row(point)
    diag = point.diagnostics
    return (;
        chi = point.params.chi,
        tol = point.params.tol,
        maxiter = point.params.maxiter,
        density = point.measurement.density,
        density_even = point.measurement.density_even,
        density_odd = point.measurement.density_odd,
        blockade_violation = point.measurement.blockade_violation,
        pxp_energy_density = point.measurement.pxp_energy_density,
        sublattice_imbalance = point.measurement.sublattice_imbalance,
        delta_density = point.delta_density,
        delta_blockade_violation = point.delta_blockade_violation,
        delta_pxp_energy_density = point.delta_pxp_energy_density,
        iterations = diag === nothing ? nothing : diag.iterations,
        residual = diag === nothing ? nothing : diag.residual,
        converged = diag === nothing ? nothing : diag.converged,
        accepted = diag === nothing ? nothing : diag.accepted,
    )
end

function _policy_namedtuple(policy::CTMTrustPolicy)
    return (;
        min_points = policy.min_points,
        require_accepted_diagnostics = policy.require_accepted_diagnostics,
        max_density_delta = policy.max_density_delta,
        max_blockade_delta = policy.max_blockade_delta,
        max_energy_delta = policy.max_energy_delta,
        max_residual = policy.max_residual,
    )
end

function _assessment_namedtuple(assessment::CTMTrustAssessment)
    return (;
        trusted = assessment.trusted,
        reason = String(assessment.reason),
        message = assessment.message,
        compared_points = assessment.compared_points,
        finite_chi_density_delta = assessment.finite_chi_density_delta,
        finite_chi_blockade_delta = assessment.finite_chi_blockade_delta,
        finite_chi_energy_delta = assessment.finite_chi_energy_delta,
        observed_max_residual = assessment.observed_max_residual,
    )
end

function main()
    Lx = _env_int("SQUAREPXP_SENS_CELL_LX", 3)
    Ly = _env_int("SQUAREPXP_SENS_CELL_LY", 3)
    D = _env_int("SQUAREPXP_SENS_D", 2)
    dt = _env_float("SQUAREPXP_SENS_DT", 0.02)
    evolve_time = _env_float("SQUAREPXP_SENS_EVOLVE_TIME", 0.02)
    chi_values = _env_int_list("SQUAREPXP_SENS_CHI_VALUES", "8,16,24,32")
    tol = _env_float("SQUAREPXP_SENS_TOL", 1e-8)
    maxiter = _env_int("SQUAREPXP_SENS_MAXITER", 100)
    seed = _env_int("SQUAREPXP_SENS_SEED", 0)
    safety = _env_float("SQUAREPXP_SENS_CALIB_SAFETY", 3.0)
    floor_value = _env_float("SQUAREPXP_SENS_CALIB_FLOOR", 1e-8)
    label = _env_value("SQUAREPXP_SENS_LABEL", "")
    csv_out = _env_value("SQUAREPXP_SENS_OUTPUT_CSV",
        joinpath(project_root, "artifacts", "ctm_chi_sensitivity.csv"))
    json_out = _env_value("SQUAREPXP_SENS_OUTPUT_JSON",
        joinpath(project_root, "artifacts", "ctm_chi_sensitivity.json"))

    length(chi_values) >= 2 ||
        throw(ArgumentError("SQUAREPXP_SENS_CHI_VALUES must contain at least 2 chi values"))

    println("Building state at $(Lx)x$(Ly), D=$D, evolve_time=$evolve_time")
    psi = _build_state(Lx, Ly, D, dt, evolve_time)

    println("Running chi-sensitivity sweep over chi = $chi_values (tol=$tol, maxiter=$maxiter)")
    points = ctm_chi_sensitivity_sweep(
        psi;
        chi_values,
        tol,
        maxiter,
        seed,
    )

    default_policy = CTMTrustPolicy()
    tight_policy = tight_ctm_trust_policy()
    calibrated = calibrated_ctm_trust_policy(points; safety, floor = floor_value)

    default_assessment = assess_ctm_trust(points; policy = default_policy)
    tight_assessment = assess_ctm_trust(points; policy = tight_policy)
    calibrated_assessment = assess_ctm_trust(points; policy = calibrated)

    mkpath(dirname(csv_out))
    write_ctm_trust_csv(points, csv_out; policy = calibrated)

    payload = (;
        label,
        cell_Lx = Lx,
        cell_Ly = Ly,
        D,
        dt,
        evolve_time,
        tol,
        maxiter,
        seed,
        chi_values,
        measurements = [_measurement_row(p) for p in points],
        default_policy = _policy_namedtuple(default_policy),
        default_assessment = _assessment_namedtuple(default_assessment),
        tight_policy = _policy_namedtuple(tight_policy),
        tight_assessment = _assessment_namedtuple(tight_assessment),
        calibrated_policy = _policy_namedtuple(calibrated),
        calibrated_assessment = _assessment_namedtuple(calibrated_assessment),
        calibrated_safety = safety,
        calibrated_floor = floor_value,
    )
    mkpath(dirname(json_out))
    open(json_out, "w") do io
        JSON3.pretty(io, payload)
        println(io)
    end

    println("--- per-chi measurements ---")
    for p in points
        @printf("  chi=%3d  density=%.6e  blockade=%.6e  energy=%.6e  iter=%s  resid=%s\n",
            p.params.chi,
            p.measurement.density,
            p.measurement.blockade_violation,
            p.measurement.pxp_energy_density,
            string(p.diagnostics === nothing ? "?" : p.diagnostics.iterations),
            string(p.diagnostics === nothing ? "?" : p.diagnostics.residual),
        )
    end
    println("default policy:    trusted=$(default_assessment.trusted) reason=$(default_assessment.reason)")
    println("tight policy:      trusted=$(tight_assessment.trusted) reason=$(tight_assessment.reason)")
    println("calibrated policy: density=$(calibrated.max_density_delta)  blockade=$(calibrated.max_blockade_delta)  energy=$(calibrated.max_energy_delta)")
    println("                   trusted=$(calibrated_assessment.trusted) reason=$(calibrated_assessment.reason)")

    println("Wrote $csv_out")
    println("Wrote $json_out")
end

main()
