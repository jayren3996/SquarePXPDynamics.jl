#!/usr/bin/env julia
#
# Compare 4x4 PBC PXP dynamics from the Neel/Z2 initial state computed two ways:
#   1. Standalone exact diagonalization in the full constrained basis (no
#      symmetry reduction; the symmetric ED in `FinitePXPEEDBenchmark` lives in
#      the trivial sector and cannot host Neel as a basis vector).
#   2. iPEPS simple-update on a 4x4 unit cell with the existing pipeline.
#
# Convention (matches `SquarePXP.jl` and `SpinOps.projector_up`):
#   bit / digit 0 means |up> = excited (Rydberg-up).
#   bit / digit 1 means |down> = ground.
# Density n_i = projector onto |up> at site i. The reported n(t) is the global
# average (1/N) sum_i <n_i>. PXP blockade forbids adjacent ups (adjacent zeros).
#
# Output: JSON with the ED trajectory plus one iPEPS trajectory per (D, dt) pair.

using Pkg
project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using LinearAlgebra
using SparseArrays
using JSON3
using Printf
using SquarePXPDynamics

# Force unbuffered output so progress is visible when stdout is a file.
say(args...) = (println(args...); flush(stdout))

# ----------- ED: full constrained basis on n x n PBC -----------

site_index(n::Int, x::Int, y::Int) = x + n * y + 1   # 1-based, x,y in 0..n-1
bit_mask(n::Int, x::Int, y::Int) = UInt32(1) << (n * n - site_index(n, x, y))

function is_constrained(state::UInt32, n::Int)
    # `bit == 0` means "up" (excited). Adjacent ups violate the PXP blockade.
    # Equivalent encoding: a site is excited iff its bit is 0. The constraint
    # is "no two adjacent excited sites" with periodic boundary conditions.
    for y = 0:(n - 1), x = 0:(n - 1)
        if (state & bit_mask(n, x, y)) == 0
            xr = mod(x + 1, n)
            yu = mod(y + 1, n)
            if (state & bit_mask(n, xr, y)) == 0
                return false
            end
            if (state & bit_mask(n, x, yu)) == 0
                return false
            end
        end
    end
    return true
end

function enumerate_constrained_basis(n::Int)
    n * n <= 30 || error("constrained basis enumeration limited to 30 sites")
    states = UInt32[]
    sizehint!(states, 2^min(n * n, 16))
    for s::UInt32 = 0:(UInt32(1) << (n * n) - UInt32(1))
        if is_constrained(s, n)
            push!(states, s)
        end
    end
    return states
end

function build_pxp_hamiltonian(states::Vector{UInt32}, n::Int)
    nsites = n * n
    index_of = Dict(state => i for (i, state) in enumerate(states))
    rows = Int[]
    cols = Int[]
    for (i, state) in enumerate(states)
        for site_bit_pos = 0:(nsites - 1)
            flipped = state ⊻ (UInt32(1) << site_bit_pos)
            j = get(index_of, flipped, 0)
            if j != 0
                push!(rows, i)
                push!(cols, j)
            end
        end
    end
    vals = ones(Float64, length(rows))
    return sparse(rows, cols, vals, length(states), length(states))
end

function neel_digit(n::Int; excited_on::Symbol = :even)
    s = UInt32(0)
    for y = 0:(n - 1), x = 0:(n - 1)
        excited = excited_on === :even ? iseven(x + y) : isodd(x + y)
        # bit = 0 means "up" (excited). So set the bit to 1 only when NOT excited.
        if !excited
            s |= bit_mask(n, x, y)
        end
    end
    return s
end

function ed_density(psi::Vector{ComplexF64}, states::Vector{UInt32}, n::Int)
    # n_i counts sites with bit == 0 ("up"). Density operator is diagonal in
    # the computational basis used here, so <n> = sum_s |c_s|^2 popcount_up(s).
    nsites = n * n
    total = 0.0
    for (i, s) in enumerate(states)
        # ups = nsites - popcount(s)
        ups = nsites - count_ones(s)
        total += abs2(psi[i]) * ups
    end
    return total / nsites
end

# ----------- iPEPS: simple-update with checkerboard initial state -----------

function ipeps_neel_trajectory(
    cell_n::Int, dt::Float64, ts::Vector{Float64};
    maxdim::Int = 4, use_ctm::Bool = true, chi::Int = 8,
    on_step::Function = (m -> nothing),
    order::Int = 1,
)
    cell = PeriodicSquareUnitCell(cell_n, cell_n)
    psi = checkerboard_square_ipeps(cell; excited_on = :even, maxdim)
    measurements = NamedTuple{(:t, :density, :density_simple, :max_truncerr),Tuple{Float64,Float64,Float64,Float64}}[]
    measure_ipeps(psi) = begin
        s = density_simple(psi)
        if use_ctm
            ctx = pepskit_ctmrg_context(psi;
                params = default_ctmrg_params(chi = chi, maxiter = 100, tol = 1e-8))
            reps = unitcell_reps(psi)
            c = sum(local_density_ctm(psi, r, ctx) for r in reps) / length(reps)
            (c, s)
        else
            (s, s)
        end
    end
    ctm0, simple0 = measure_ipeps(psi)
    m0 = (t = 0.0, density = ctm0, density_simple = simple0, max_truncerr = 0.0)
    push!(measurements, m0)
    on_step(measurements)
    @printf("    t=%.3f n_ctm=%.6f n_simple=%.6f max_trunc=%.2e\n", 0.0, ctm0, simple0, 0.0)
    flush(stdout)

    params = TrotterParams(dt, order, :real, maxdim, 1e-12; schedule = :serial)
    t_now = 0.0
    target_idx = 2
    max_trunc = 0.0
    while target_idx <= length(ts)
        t_target = ts[target_idx]
        step = t_target - t_now
        step <= 0 && (target_idx += 1; continue)
        log = evolve!(psi, step; params = params)
        max_trunc = max(max_trunc, log.max_truncerr)
        t_now = t_target
        ctm_n, simple_n = measure_ipeps(psi)
        push!(measurements,
            (t = t_now, density = ctm_n, density_simple = simple_n, max_truncerr = max_trunc))
        on_step(measurements)
        @printf("    t=%.3f n_ctm=%.6f n_simple=%.6f max_trunc=%.2e\n",
            t_now, ctm_n, simple_n, max_trunc)
        flush(stdout)
        target_idx += 1
    end
    return measurements
end

# ----------- driver -----------

function main()
    if !isempty(get(ENV, "SQUAREPXP_CTM_BLAS_THREADS", "")) ||
       !isempty(get(ENV, "SQUAREPXP_CTM_STRIDED_THREADS", ""))
        configure_ctm_threading_from_env!()
    end
    n = 4
    total_time = parse(Float64, get(ENV, "NEEL_TOTAL_TIME", "3.0"))
    dt_meas = parse(Float64, get(ENV, "NEEL_DT_MEAS", "0.1"))
    ipeps_dts_str = get(ENV, "NEEL_IPEPS_DTS", "0.01")
    ipeps_dts = parse.(Float64, split(ipeps_dts_str, ","))
    ipeps_maxdims_str = get(ENV, "NEEL_IPEPS_MAXDIMS", "4")
    ipeps_maxdims = parse.(Int, split(ipeps_maxdims_str, ","))
    chi = parse(Int, get(ENV, "NEEL_CTM_CHI", "8"))
    use_ctm = parse(Bool, get(ENV, "NEEL_USE_CTM", "true"))
    order = parse(Int, get(ENV, "NEEL_TROTTER_ORDER", "1"))
    abort_threshold = parse(Float64, get(ENV, "NEEL_ABORT_THRESHOLD", "Inf"))

    ts = collect(0.0:dt_meas:total_time)
    say("Neel vs ED comparison, n=$n, total=$total_time, dt_meas=$dt_meas")

    output_path = get(ENV, "NEEL_OUTPUT_JSON",
        joinpath(project_root, "artifacts", "neel_ed_vs_ipeps_4x4.json"))
    mkpath(dirname(output_path))

    payload = Dict{String,Any}(
        "cell_n" => n,
        "total_time" => total_time,
        "dt_meas" => dt_meas,
        "ts" => ts,
        "ed" => nothing,
        "ipeps" => Dict[],
    )
    save() = open(io -> JSON3.pretty(io, payload), output_path, "w")

    # ED
    say("Enumerating constrained basis...")
    states = enumerate_constrained_basis(n)
    say("  basis size = $(length(states))")
    say("Building Hamiltonian...")
    H = build_pxp_hamiltonian(states, n)
    say("  nnz = $(nnz(H))")
    say("Diagonalizing (dense)...")
    Hd = Matrix(H)
    issymmetric(Hd) || @warn "H not symmetric — check construction"
    F = eigen(Hermitian(Hd))
    say("  spectrum range [$(extrema(F.values))]")

    neel = neel_digit(n; excited_on = :even)
    idx_neel = findfirst(==(neel), states)
    idx_neel === nothing && error("Neel state not in constrained basis")
    psi0 = zeros(ComplexF64, length(states))
    psi0[idx_neel] = 1.0
    coeff = F.vectors' * psi0
    ed_n = Float64[]
    for t in ts
        phase = exp.(-im .* t .* F.values)
        psi_t = F.vectors * (phase .* coeff)
        push!(ed_n, ed_density(psi_t, states, n))
    end

    say("ED done. n(0)=$(ed_n[1])  n(T)=$(ed_n[end])")
    payload["ed"] = Dict("density" => ed_n, "basis_size" => length(states))
    save()

    # iPEPS
    for maxdim in ipeps_maxdims, dt in ipeps_dts
        label = "D=$(maxdim)_dt=$(dt)_chi=$(chi)_order=$(order)"
        say("iPEPS $label …")
        slot = Dict(
            "maxdim" => maxdim, "dt" => dt, "chi" => chi, "order" => order,
            "density" => Float64[], "density_simple" => Float64[],
            "ts" => Float64[], "max_truncerr" => 0.0,
            "elapsed_seconds" => 0.0, "completed" => false,
        )
        push!(payload["ipeps"], slot)
        on_step = (measurements) -> begin
            slot["density"] = [r.density for r in measurements]
            slot["density_simple"] = [r.density_simple for r in measurements]
            slot["ts"] = [r.t for r in measurements]
            slot["max_truncerr"] = measurements[end].max_truncerr
            save()
            # Abort early if |Δn| exceeds threshold against ED reference.
            if !isinf(abort_threshold)
                last_t = measurements[end].t
                idx = findfirst(x -> isapprox(x, last_t; atol = 1e-6), ts)
                if idx !== nothing
                    diff = abs(measurements[end].density - ed_n[idx])
                    if diff > abort_threshold
                        say(@sprintf("  ABORT at t=%.3f: |Δn|=%.2e > %.2e",
                            last_t, diff, abort_threshold))
                        throw(ErrorException("threshold breach |Δn|=$diff at t=$last_t"))
                    end
                end
            end
        end
        t0 = time()
        try
            traj = ipeps_neel_trajectory(n, dt, ts;
                maxdim, use_ctm = use_ctm, chi = chi, on_step = on_step, order = order)
            slot["completed"] = true
            slot["elapsed_seconds"] = time() - t0
            say(@sprintf("  done in %.2fs, max_truncerr=%.2e, n_ctm(T)=%.4f n_simple(T)=%.4f",
                slot["elapsed_seconds"], slot["max_truncerr"],
                isempty(slot["density"]) ? 0.0 : slot["density"][end],
                isempty(slot["density_simple"]) ? 0.0 : slot["density_simple"][end]))
        catch e
            slot["elapsed_seconds"] = time() - t0
            slot["error"] = sprint(showerror, e)
            say("  CRASHED at t=$(isempty(slot["ts"]) ? "?" : last(slot["ts"])): $(typeof(e))")
            say("  partial data through $(length(slot["density"])) points saved")
            # Unwrap TaskFailedException to show the real underlying error.
            e_to_show = e
            while e_to_show isa TaskFailedException
                e_to_show = e_to_show.task.result
            end
            say("  underlying error: ", sprint(showerror, e_to_show))
            for (exc, bt) in current_exceptions()
                println("    stacktrace fragment: ", exc)
                Base.show_backtrace(stdout, bt)
                println()
                break
            end
            flush(stdout)
        end
        save()
    end

    say("Wrote $output_path")

    say("")
    say("--- max |ED - iPEPS| over t in [0, T] ---")
    for r in payload["ipeps"]
        if isempty(r["density"]) || length(r["density"]) > length(ed_n)
            continue
        end
        nlen = length(r["density"])
        max_diff = maximum(abs.(ed_n[1:nlen] .- r["density"]))
        say(@sprintf("  D=%d dt=%.3f  max|Δn|=%.4e (over %d points)",
            r["maxdim"], r["dt"], max_diff, nlen))
    end
end

main()
