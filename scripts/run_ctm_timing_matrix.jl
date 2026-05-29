#!/usr/bin/env julia
#
# CTM throughput timing matrix.
#
# Builds a representative iPEPS state (product :down + short PXP real-time
# evolution at given D), then times `measure_ctm` across a grid of in-process
# threading knobs: (strided_threads, blas_threads, strided_threaded_mul,
# pepskit_scheduler), at each chi in SQUAREPXP_TIMING_CHI_VALUES.
#
# The JULIA_NUM_THREADS dimension cannot be varied inside one Julia session,
# so this script is meant to be re-invoked once per thread count and have its
# rows accumulated into one CSV by setting SQUAREPXP_TIMING_LABEL and a shared
# SQUAREPXP_TIMING_OUTPUT_CSV. See docs/superpowers/notes/2026-05-26-ctm-throughput-recipe.md.

using Pkg

project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using JSON3
using Printf
using SquarePXPDynamics

function _env_value(name::AbstractString, default::AbstractString)
    raw = get(ENV, String(name), "")
    return isempty(strip(raw)) ? String(default) : String(strip(raw))
end

function _env_int(name::AbstractString, default::Int)
    return parse(Int, _env_value(name, string(default)))
end

function _env_float(name::AbstractString, default::Float64)
    return parse(Float64, _env_value(name, string(default)))
end

function _env_symbol(name::AbstractString, default::Symbol)
    return Symbol(_env_value(name, String(default)))
end

function _env_int_list(name::AbstractString, default::AbstractString)
    raw = strip(_env_value(name, default))
    isempty(raw) && return Int[]
    return parse.(Int, split(raw, ","))
end

function _parse_bool(value::AbstractString)
    v = lowercase(strip(value))
    v in ("1", "true", "yes", "on") && return true
    v in ("0", "false", "no", "off") && return false
    throw(ArgumentError("expected boolean, got \"$value\""))
end

function _env_bool_list(name::AbstractString, default::AbstractString)
    raw = strip(_env_value(name, default))
    isempty(raw) && return Bool[]
    return _parse_bool.(split(raw, ","))
end

function _env_symbol_list(name::AbstractString, default::AbstractString)
    raw = strip(_env_value(name, default))
    isempty(raw) && return Symbol[]
    return Symbol.(split(raw, ","))
end

function _build_state(Lx::Int, Ly::Int, D::Int, dt::Float64, evolve_time::Float64)
    cell = PeriodicSquareUnitCell(Lx, Ly)
    psi = product_square_ipeps(cell; state = :down, maxdim = D)
    if evolve_time > 0 && D > 1
        params = TrotterParams(dt, 1, :real, D, 1e-12; schedule = :serial)
        evolve!(psi, evolve_time; params)
    end
    return psi
end

function _time_once(psi, ctmrg_params)
    elapsed = @elapsed begin
        summary = measure_ctm(psi; params = ctmrg_params)
    end
    return elapsed, summary
end

function _safe_run(psi, ctmrg_params, reps::Int)
    timings = Float64[]
    warm_elapsed = NaN
    warm_summary = nothing
    try
        warm_elapsed, warm_summary = _time_once(psi, ctmrg_params)
        for _ in 1:reps
            t, _ = _time_once(psi, ctmrg_params)
            push!(timings, t)
        end
    catch err
        return (; ok = false, error = sprint(showerror, err), timings, warm_elapsed,
                summary = warm_summary)
    end
    return (; ok = true, error = nothing, timings, warm_elapsed, summary = warm_summary)
end

function _stats(timings)
    isempty(timings) && return (; n = 0, median = NaN, min = NaN, max = NaN, mean = NaN)
    sorted = sort(timings)
    n = length(sorted)
    med = isodd(n) ? sorted[(n + 1) ÷ 2] : 0.5 * (sorted[n ÷ 2] + sorted[n ÷ 2 + 1])
    return (; n, median = med, min = sorted[1], max = sorted[end], mean = sum(sorted) / n)
end

function _diag_field(::Nothing, ::Symbol)
    return nothing
end

function _diag_field(diag, field)
    return getfield(diag, field)
end

function _csv_field(::Nothing)
    return ""
end

function _csv_field(value::Bool)
    return value ? "true" : "false"
end

function _csv_field(value::Real)
    return string(value)
end

function _csv_field(value)
    s = string(value)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function _row_dict(label, julia_threads, blas, strided, threaded_mul, scheduler,
                   chi, reps, run, applied, stats, evolve_time, Lx, Ly, D)
    diag = run.summary === nothing ? nothing : run.summary.diagnostics
    return (;
        label,
        julia_threads,
        blas_threads_requested = blas,
        strided_threads_requested = strided,
        strided_threaded_mul_requested = threaded_mul,
        pepskit_scheduler_requested = scheduler,
        blas_threads_applied = applied.blas_threads,
        strided_threads_applied = applied.strided_threads,
        strided_threaded_mul_applied = applied.strided_threaded_mul,
        pepskit_scheduler_applied = applied.pepskit_scheduler,
        chi,
        reps,
        cell_Lx = Lx,
        cell_Ly = Ly,
        D,
        evolve_time,
        ok = run.ok,
        error = run.error,
        warm_seconds = run.warm_elapsed,
        median_seconds = stats.median,
        min_seconds = stats.min,
        max_seconds = stats.max,
        mean_seconds = stats.mean,
        timings = run.timings,
        ctm_iterations = _diag_field(diag, :iterations),
        ctm_residual = _diag_field(diag, :residual),
        ctm_converged = _diag_field(diag, :converged),
        ctm_accepted = _diag_field(diag, :accepted),
    )
end

const _CSV_COLUMNS = (
    :label, :julia_threads,
    :blas_threads_requested, :strided_threads_requested,
    :strided_threaded_mul_requested, :pepskit_scheduler_requested,
    :blas_threads_applied, :strided_threads_applied,
    :strided_threaded_mul_applied, :pepskit_scheduler_applied,
    :chi, :reps, :cell_Lx, :cell_Ly, :D, :evolve_time,
    :ok, :error, :warm_seconds, :median_seconds, :min_seconds,
    :max_seconds, :mean_seconds,
    :ctm_iterations, :ctm_residual, :ctm_converged, :ctm_accepted,
)

function _write_csv(path::AbstractString, rows; append::Bool)
    mkpath(dirname(path))
    open(path, append ? "a" : "w") do io
        if !append || filesize(path) == 0
            println(io, join(string.(_CSV_COLUMNS), ","))
        end
        for row in rows
            fields = [_csv_field(getproperty(row, c)) for c in _CSV_COLUMNS]
            println(io, join(fields, ","))
        end
    end
    return path
end

function _write_json(path::AbstractString, payload)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.pretty(io, payload)
        println(io)
    end
    return path
end

function main()
    Lx = _env_int("SQUAREPXP_TIMING_CELL_LX", 3)
    Ly = _env_int("SQUAREPXP_TIMING_CELL_LY", 3)
    D = _env_int("SQUAREPXP_TIMING_D", 2)
    dt = _env_float("SQUAREPXP_TIMING_DT", 0.02)
    evolve_time = _env_float("SQUAREPXP_TIMING_EVOLVE_TIME", 0.02)
    chi_values = _env_int_list("SQUAREPXP_TIMING_CHI_VALUES", "16")
    blas_values = _env_int_list("SQUAREPXP_TIMING_BLAS_VALUES", "1")
    strided_default = string(Threads.nthreads())
    strided_values = _env_int_list("SQUAREPXP_TIMING_STRIDED_VALUES", strided_default)
    threaded_mul_values = _env_bool_list("SQUAREPXP_TIMING_THREADED_MUL_VALUES", "false,true")
    scheduler_values = _env_symbol_list("SQUAREPXP_TIMING_SCHEDULER_VALUES", "default,dynamic")
    reps = _env_int("SQUAREPXP_TIMING_REPS", 3)
    label = _env_value("SQUAREPXP_TIMING_LABEL", "")
    csv_out = _env_value("SQUAREPXP_TIMING_OUTPUT_CSV",
                        joinpath(project_root, "artifacts", "ctm_timing_matrix.csv"))
    json_out = _env_value("SQUAREPXP_TIMING_OUTPUT_JSON",
                        joinpath(project_root, "artifacts", "ctm_timing_matrix.json"))
    append = _parse_bool(_env_value("SQUAREPXP_TIMING_APPEND", "false"))

    println("Building reference state at $(Lx)x$(Ly), D=$D, evolve_time=$evolve_time")
    psi = _build_state(Lx, Ly, D, dt, evolve_time)

    rows = NamedTuple[]

    julia_threads = Threads.nthreads()

    println("Grid: blas=$blas_values strided=$strided_values threaded_mul=$threaded_mul_values scheduler=$scheduler_values chi=$chi_values reps=$reps julia_threads=$julia_threads")

    for blas in blas_values, strided in strided_values, threaded_mul in threaded_mul_values, scheduler in scheduler_values
        applied = try
            configure_ctm_threading!(
                blas_threads = blas,
                strided_threads = strided,
                strided_threaded_mul = threaded_mul,
                pepskit_scheduler = scheduler,
            )
        catch err
            println("  config skip (configure failed): blas=$blas strided=$strided threaded_mul=$threaded_mul scheduler=$scheduler: $err")
            continue
        end

        for chi in chi_values
            ctmrg_params = PEPSKitCTMRGParams(chi, 1e-8, 50, 0; seed = 0)
            run = _safe_run(psi, ctmrg_params, reps)
            stats = _stats(run.timings)
            row = _row_dict(label, julia_threads, blas, strided, threaded_mul,
                            scheduler, chi, reps, run, applied, stats,
                            evolve_time, Lx, Ly, D)
            push!(rows, row)
            status = run.ok ? "ok" : "ERR"
            warm = isnan(run.warm_elapsed) ? "n/a" : @sprintf("%.3fs", run.warm_elapsed)
            med = isnan(stats.median) ? "n/a" : @sprintf("%.3fs", stats.median)
            println("  [$(status)] julia=$(julia_threads) blas=$(applied.blas_threads) strided=$(applied.strided_threads) tmul=$(applied.strided_threaded_mul) sched=$(applied.pepskit_scheduler) chi=$chi warm=$warm median=$med n=$(stats.n)")
        end
    end

    payload = (; julia_threads, Lx, Ly, D, dt, evolve_time, chi_values, blas_values,
               strided_values, threaded_mul_values, scheduler_values, reps,
               label, rows)
    _write_csv(csv_out, rows; append)
    _write_json(json_out, payload)
    println("Wrote $(csv_out)")
    println("Wrote $(json_out)")
end

main()
