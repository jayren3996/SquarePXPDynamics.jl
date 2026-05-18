#!/usr/bin/env julia

using Pkg

project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using JSON3
using PEPSKit
using SquarePXPDynamics
using SquarePXPDynamics.PEPSKitMeasurements:
    _pepskit_density_operator, _pepskit_twosite_nn_operator, _pepskit_pxp_energy_operator

function _env_value(name::String, default::AbstractString)
    value = get(ENV, name, "")
    return isempty(value) ? String(default) : value
end

function _env_bool(name::String, default::Bool)
    value = lowercase(strip(_env_value(name, string(default))))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    throw(ArgumentError("$name must be one of 1,true,yes,on,0,false,no,off"))
end

function _env_int(name::String, default::Int)
    return parse(Int, _env_value(name, string(default)))
end

function _env_float(name::String, default::Float64)
    return parse(Float64, _env_value(name, string(default)))
end

function _env_symbol(name::String, default::Symbol)
    return Symbol(_env_value(name, String(default)))
end

function _env_int_list(name::String, default::Vector{Int})
    raw = strip(_env_value(name, join(string.(default), ",")))
    isempty(raw) && return Int[]
    return [parse(Int, strip(tok)) for tok in split(raw, ",") if !isempty(strip(tok))]
end

function _env_float_list(name::String, default::Vector{Float64})
    raw = strip(_env_value(name, join(string.(default), ",")))
    isempty(raw) && return Float64[]
    return [parse(Float64, strip(tok)) for tok in split(raw, ",") if !isempty(strip(tok))]
end

function _is_integer_multiple(total::Float64, dt::Float64)
    nsteps = round(Int, total / dt)
    return isapprox(nsteps * dt, total; atol = 1e-12, rtol = 1e-10)
end

function _csv_value(value::Nothing)
    return ""
end

function _csv_value(value::Bool)
    return string(value)
end

function _csv_value(value)
    return string(value)
end

const CSV_FIELDS = (
    :step,
    :time,
    :density_simple,
    :blockade_violation_simple,
    :pxp_energy_density_simple,
    :exact_finite_density,
    :ctm_density,
    :ctm_density_imag,
    :ctm_density_even,
    :ctm_density_odd,
    :ctm_blockade_violation,
    :ctm_blockade_imag,
    :ctm_pxp_energy_density,
    :ctm_energy_imag,
    :ctm_iterations,
    :ctm_residual,
    :ctm_converged,
    :ctm_accepted,
    :ctm_seconds,
    :max_truncerr,
    :mean_bond_entropy,
    :max_bond_entropy,
    :log_norm,
)

function _write_csv(path::AbstractString, rows)
    open(path, "w") do io
        println(io, join(String.(CSV_FIELDS), ","))
        for row in rows
            println(io, join((_csv_value(row[field]) for field in CSV_FIELDS), ","))
        end
    end
    return path
end

n = _env_int("SQUAREPXP_IPEPS_N", 3)
total_time = _env_float("SQUAREPXP_IPEPS_TOTAL_TIME", 2.0)
dt = _env_float("SQUAREPXP_IPEPS_DT", 0.02)
measure_every = _env_int("SQUAREPXP_IPEPS_MEASURE_EVERY", 5)
D_values = _env_int_list("SQUAREPXP_IPEPS_D", [2])
cutoff_values = _env_float_list("SQUAREPXP_IPEPS_CUTOFF", [1e-10])
order = _env_int("SQUAREPXP_IPEPS_ORDER", 1)
schedule = _env_symbol("SQUAREPXP_IPEPS_SCHEDULE", :serial)
initial_state = _env_symbol("SQUAREPXP_IPEPS_INITIAL_STATE", :down)
exact_finite = _env_bool("SQUAREPXP_IPEPS_EXACT_FINITE", true)
exact_finite_max_sites = _env_int("SQUAREPXP_IPEPS_EXACT_FINITE_MAX_SITES", 12)

ctm_chi_values = _env_int_list("SQUAREPXP_IPEPS_CTM_CHI", Int[8])
ctm_tol = _env_float("SQUAREPXP_IPEPS_CTM_TOL", 1e-8)
ctm_maxiter = _env_int("SQUAREPXP_IPEPS_CTM_MAXITER", 200)
ctm_verbosity = _env_int("SQUAREPXP_IPEPS_CTM_VERBOSITY", 0)
ctm_seed = _env_int("SQUAREPXP_IPEPS_CTM_SEED", 1234)
ctm_skip_initial = _env_bool("SQUAREPXP_IPEPS_CTM_SKIP_INITIAL", true)

_is_integer_multiple(total_time, dt) ||
    throw(ArgumentError("SQUAREPXP_IPEPS_TOTAL_TIME must be an integer multiple of SQUAREPXP_IPEPS_DT"))
measure_every >= 1 || throw(ArgumentError("SQUAREPXP_IPEPS_MEASURE_EVERY must be at least 1"))
isempty(ctm_chi_values) && throw(ArgumentError("SQUAREPXP_IPEPS_CTM_CHI must contain at least one chi"))
isempty(D_values) && throw(ArgumentError("SQUAREPXP_IPEPS_D must contain at least one D"))
isempty(cutoff_values) && throw(ArgumentError("SQUAREPXP_IPEPS_CUTOFF must contain at least one cutoff"))

threading = configure_ctm_threading_from_env!()
@info "CTM threading" threading

ctm_params = [
    PEPSKitCTMRGParams(chi, ctm_tol, ctm_maxiter, ctm_verbosity; seed = ctm_seed)
    for chi in ctm_chi_values
]
primary_idx = argmax([p.chi for p in ctm_params])

artifact_dir = _env_value(
    "SQUAREPXP_IPEPS_ARTIFACT_DIR",
    joinpath(project_root, "artifacts", "m3-systematic-ctm"),
)
mkpath(artifact_dir)
run_label = _env_value("SQUAREPXP_IPEPS_RUN_LABEL", "ctm-traj")

nsteps = round(Int, total_time / dt)

function _try_observable(f)
    try
        return f(), nothing
    catch err
        return nothing, sprint(showerror, err)
    end
end

function _site_density_complex(psi, c, ctx)
    op = _pepskit_density_operator(psi.unitcell, c)
    return ComplexF64(PEPSKit.expectation_value(ctx.peps, op, ctx.env))
end

function _nn_density_complex(psi, c, dir, ctx)
    op = _pepskit_twosite_nn_operator(psi.unitcell, c, dir)
    return ComplexF64(PEPSKit.expectation_value(ctx.peps, op, ctx.env))
end

function _pxp_energy_density_complex(psi, ctx)
    op = _pepskit_pxp_energy_operator(psi.unitcell)
    return ComplexF64(PEPSKit.expectation_value(ctx.peps, op, ctx.env))
end

function _split_re_im(z::Union{Nothing,Complex})
    z === nothing && return (nothing, nothing)
    return (real(z), abs(imag(z)))
end

function _ctm_sample_one(psi, params::PEPSKitCTMRGParams)
    ctx = pepskit_ctmrg_context(psi; params = params)
    reps = psi.unitcell.reps
    raw_density, density_err = _try_observable(
        () -> sum(_site_density_complex(psi, c, ctx) for c in reps) / length(reps),
    )
    raw_density_even, _ = _try_observable(
        () -> begin
            evens = [c for c in reps if iseven(c.x + c.y)]
            isempty(evens) ? nothing : sum(_site_density_complex(psi, c, ctx) for c in evens) / length(evens)
        end,
    )
    raw_density_odd, _ = _try_observable(
        () -> begin
            odds = [c for c in reps if isodd(c.x + c.y)]
            isempty(odds) ? nothing : sum(_site_density_complex(psi, c, ctx) for c in odds) / length(odds)
        end,
    )
    raw_blockade, blockade_err = _try_observable(() -> begin
        total = 0.0 + 0.0im
        count = 0
        for c in reps, dir in (:right, :up)
            total += _nn_density_complex(psi, c, dir, ctx)
            count += 1
        end
        total / count
    end)
    raw_energy, energy_err = _try_observable(() -> _pxp_energy_density_complex(psi, ctx))
    density, density_imag = _split_re_im(raw_density)
    density_even, density_even_imag = _split_re_im(raw_density_even)
    density_odd, density_odd_imag = _split_re_im(raw_density_odd)
    blockade, blockade_imag = _split_re_im(raw_blockade)
    energy, energy_imag = _split_re_im(raw_energy)
    diag = ctx.diagnostics
    return (;
        density,
        density_imag,
        density_err,
        density_even,
        density_even_imag,
        density_odd,
        density_odd_imag,
        blockade,
        blockade_imag,
        blockade_err,
        energy,
        energy_imag,
        energy_err,
        iterations = diag === nothing ? nothing : diag.iterations,
        residual = diag === nothing ? nothing : diag.residual,
        converged = diag === nothing ? nothing : diag.converged,
        accepted = diag === nothing ? nothing : diag.accepted,
    )
end

function _ctm_sample(psi, skip::Bool)
    skip && return (nothing, [], 0.0)
    points = Any[]
    primary_summary = nothing
    elapsed = 0.0
    for (i, params) in enumerate(ctm_params)
        t0 = time()
        summary = _ctm_sample_one(psi, params)
        dt_ctm = time() - t0
        elapsed += dt_ctm
        push!(points, (; chi = params.chi, summary, seconds = dt_ctm))
        if i == primary_idx
            primary_summary = summary
        end
    end
    return primary_summary, points, elapsed
end

function _sample!(rows, step::Int, time::Float64, psi, evolution; skip_ctm = false)
    ctm_summary, ctm_points, ctm_seconds = _ctm_sample(psi, skip_ctm)
    if ctm_summary === nothing
        ctm_density = nothing
        ctm_density_imag = nothing
        ctm_density_even = nothing
        ctm_density_odd = nothing
        ctm_blockade = nothing
        ctm_blockade_imag = nothing
        ctm_energy = nothing
        ctm_energy_imag = nothing
        ctm_iter = nothing
        ctm_residual = nothing
        ctm_converged = nothing
        ctm_accepted = nothing
        ctm_blockade_error = nothing
        ctm_energy_error = nothing
    else
        ctm_density = ctm_summary.density
        ctm_density_imag = ctm_summary.density_imag
        ctm_density_even = ctm_summary.density_even
        ctm_density_odd = ctm_summary.density_odd
        ctm_blockade = ctm_summary.blockade
        ctm_blockade_imag = ctm_summary.blockade_imag
        ctm_energy = ctm_summary.energy
        ctm_energy_imag = ctm_summary.energy_imag
        ctm_iter = ctm_summary.iterations
        ctm_residual = ctm_summary.residual
        ctm_converged = ctm_summary.converged
        ctm_accepted = ctm_summary.accepted
        ctm_blockade_error = ctm_summary.blockade_err
        ctm_energy_error = ctm_summary.energy_err
    end
    push!(
        rows,
        (;
            step,
            time,
            density_simple = density_simple(psi),
            blockade_violation_simple = blockade_violation_simple(psi),
            pxp_energy_density_simple = pxp_energy_density_simple(psi),
            exact_finite_density = exact_finite ?
                exact_density_finite(psi; max_sites = exact_finite_max_sites) : nothing,
            ctm_density,
            ctm_density_imag,
            ctm_density_even,
            ctm_density_odd,
            ctm_blockade_violation = ctm_blockade,
            ctm_blockade_imag,
            ctm_pxp_energy_density = ctm_energy,
            ctm_energy_imag,
            ctm_blockade_error,
            ctm_energy_error,
            ctm_iterations = ctm_iter,
            ctm_residual,
            ctm_converged,
            ctm_accepted,
            ctm_seconds,
            ctm_finite_chi = [
                (;
                    chi = pt.chi,
                    density = pt.summary.density,
                    density_imag = pt.summary.density_imag,
                    density_even = pt.summary.density_even,
                    density_odd = pt.summary.density_odd,
                    blockade_violation = pt.summary.blockade,
                    blockade_imag = pt.summary.blockade_imag,
                    pxp_energy_density = pt.summary.energy,
                    energy_imag = pt.summary.energy_imag,
                    blockade_error = pt.summary.blockade_err,
                    energy_error = pt.summary.energy_err,
                    iterations = pt.summary.iterations,
                    residual = pt.summary.residual,
                    converged = pt.summary.converged,
                    accepted = pt.summary.accepted,
                    seconds = pt.seconds,
                ) for pt in ctm_points
            ],
            max_truncerr = evolution === nothing ? 0.0 : evolution.max_truncerr,
            mean_bond_entropy = evolution === nothing ? 0.0 : evolution.mean_bond_entropy,
            max_bond_entropy = evolution === nothing ? 0.0 : evolution.max_bond_entropy,
            log_norm = log_norm(psi),
        ),
    )
    return rows
end

cell = PeriodicSquareUnitCell(n, n)

case_summaries = []

for D in D_values
    for cutoff in cutoff_values
        case_label = "n$(n)-dt$(dt)-T$(total_time)-D$(D)-cut$(cutoff)"
        println("==> case $(case_label)")
        psi = product_square_ipeps(cell; state = initial_state, maxdim = D)
        trotter = TrotterParams(dt, order, :real, D, cutoff; schedule)

        rows = NamedTuple[]
        case_t0 = time()
        _sample!(rows, 0, 0.0, psi, nothing; skip_ctm = ctm_skip_initial)
        for step = 1:nsteps
            evolution = evolve!(psi, dt; params = trotter)
            if step % measure_every == 0 || step == nsteps
                _sample!(rows, step, step * dt, psi, evolution)
            end
        end
        case_seconds = time() - case_t0

        payload = (;
            config = (;
                n,
                total_time,
                dt,
                measure_every,
                D,
                cutoff,
                order,
                schedule = String(schedule),
                initial_state = String(initial_state),
                exact_finite,
                exact_finite_max_sites,
                ctm_chi_values,
                ctm_tol,
                ctm_maxiter,
                ctm_verbosity,
                ctm_seed,
                ctm_primary_chi = ctm_params[primary_idx].chi,
                ctm_skip_initial,
                run_label,
                case_seconds,
            ),
            threading = (;
                julia_threads = threading.julia_threads,
                blas_threads = threading.blas_threads,
                strided_threads = threading.strided_threads,
                strided_threaded_mul = threading.strided_threaded_mul,
                pepskit_scheduler = String(threading.pepskit_scheduler),
            ),
            samples = rows,
        )

        json_path = joinpath(artifact_dir, "$(run_label)-$(case_label).json")
        csv_path = joinpath(artifact_dir, "$(run_label)-$(case_label).csv")
        open(json_path, "w") do io
            JSON3.write(io, payload)
            write(io, '\n')
        end
        _write_csv(csv_path, rows)
        push!(case_summaries, (; D, cutoff, case_seconds, json = json_path, csv = csv_path))
        println("    wrote $(json_path)")
        println("    wrote $(csv_path)")
        println("    case_seconds = $(round(case_seconds, digits = 2))")
    end
end

println("\n=== campaign summary ===")
for s in case_summaries
    println("D=$(s.D) cutoff=$(s.cutoff) seconds=$(round(s.case_seconds, digits = 2)) json=$(s.json)")
end
