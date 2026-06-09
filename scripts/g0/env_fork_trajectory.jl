# env_fork_trajectory.jl — parameterized env-aware fork driver.
#
# Bare-evolve a 4x4 Neel D=2 state to TSTART (rel_floor=0, min bond dim 2), then
# fork and continue with the env-aware star update on the :serial order-2 schedule,
# measuring exact_density_finite every 0.1 vs the canonical same-t ED.
#
# USAGE: julia env_fork_trajectory.jl TSTART TEND DT ENV TRUNC_ALG [TAG]
#   TSTART TEND : fork window (e.g. 1.6 2.8)
#   DT          : Trotter step (0.02 | 0.01 | 0.005)
#   ENV         : :exact_cluster | :bare | :ntu
#   TRUNC_ALG   : :full | :als
#   TAG         : optional label for the header
#
# TWO STANDING RUNS (per review wedw4yosa):
#   dt-convergence (decisive for claim c): 1.6 1.8 0.01 exact_cluster full dtconv
#   effect-A through revival:              1.6 2.8 0.02 exact_cluster full window
# Compare exact_cluster(dt=0.01,t=1.8) vs (dt=0.02,t=1.8)=0.22321 vs true ED(1.8)=0.24370:
# error SHRINKS as dt->0  => env converges to TRUE physics (claim c rescued);
# error tracks D4 / grows => shared-dt Trotter<->truncation coincidence (c refuted).

using SquarePXPDynamics
using TensorKit
using Printf

include(joinpath(@__DIR__, "star_env.jl"))  # -> project_star_env!

# canonical same-t ED (artifacts/neel_to_revival_4x4.json, 0.2 grid), verified.
const ED = Dict(1.0=>0.15549,1.2=>0.10057,1.4=>0.10056,1.6=>0.15558,1.8=>0.24370,
                2.0=>0.33735,2.2=>0.41600,2.4=>0.46696,2.6=>0.48306,2.8=>0.46230)

function evolve_env_serial!(state, total_time; dt, maxdim, cutoff, env, patch, trunc_alg)
    reps = state.unitcell.reps
    h = dt / 2
    subs = vcat([(c, h) for c in reps], [(c, h) for c in Base.reverse(reps)])
    nsteps = round(Int, total_time / dt)
    maxcond = 0.0; nfb = 0; minkept = typemax(Int)
    for _ in 1:nsteps
        for (center, layer_dt) in subs
            info = project_star_env!(state, center; dt = layer_dt, maxdim = maxdim,
                                     cutoff = cutoff, env = env, patch = patch,
                                     trunc_alg = trunc_alg, evolution = :real, projected = true)
            maxcond = max(maxcond, maximum(values(info.conds)))
            nfb += count(values(info.fellbacks))
            minkept = min(minkept, minimum(values(info.keptdims)))
        end
    end
    return (maxcond = maxcond, nfb = nfb, minkept = minkept)
end

function main(tstart, tend, dt, env, trunc_alg, tag, baredt, patch, D)
    @printf("RSS start = %.1f MB ; total = %.1f MB\n", Sys.maxrss()/2^20, Sys.total_memory()/2^20)
    @printf("=== env_fork tag=%s : fork@t=%.2f -> %.2f , env_dt=%.3f (bare_dt=%.3f) , env=%s , patch=%s , alg=%s , D=%d ===\n",
            tag, tstart, tend, dt, baredt, env, patch, trunc_alg, D); flush(stdout)
    chunk = 0.1
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = D)
    _ = exact_density_finite(psi; max_sites = 16)   # warm JIT

    # bare-evolve to the fork point at the FIXED bare_dt (so the fork STATE is
    # identical across env_dt sweeps -> a clean dt-convergence test of the env
    # evolution alone, no bare-Trotter confound). rel_floor=0 -> min bond dim 2.
    t = 0.0
    while t < tstart - 1e-9
        evolve_pepskit!(psi, chunk; dt = baredt, order = 2, schedule = :serial,
                        maxdim = D, cutoff = 1e-12, rel_floor = 0.0)
        t = round(t + chunk, digits = 1)
    end
    @printf("bare-evolved to fork t=%.2f ; n=%.6f\n", tstart, exact_density_finite(psi; max_sites=16))
    flush(stdout)

    # warm the env-build JIT once (absorb cold ncon compile)
    tw0 = time()
    let warm = deepcopy(psi)
        project_star_env!(warm, cell.reps[1]; dt = dt/2, maxdim = D, cutoff = 1e-12,
                          env = env, patch = patch, trunc_alg = trunc_alg,
                          evolution = :real, projected = true)
    end
    @printf("env-JIT warmup = %.1f s\n", time() - tw0); flush(stdout)

    @printf("\n%-5s %-11s %-9s %-9s %-9s %-9s\n", "t", "n_env", "n_ED", "maxcond", "wall_s", "diag")
    @printf("%-5.2f %-11.6f %-9s %-9s %-9s\n",
            tstart, exact_density_finite(psi; max_sites=16),
            haskey(ED, round(tstart,digits=1)) ? @sprintf("%.5f", ED[round(tstart,digits=1)]) : "-", "-", "-")
    flush(stdout)

    t = tstart
    while t < tend - 1e-9
        te0 = time()
        info = evolve_env_serial!(psi, chunk; dt = dt, maxdim = D, cutoff = 1e-12,
                                  env = env, patch = patch, trunc_alg = trunc_alg)
        ewall = time() - te0
        t = round(t + chunk, digits = 1)
        n = exact_density_finite(psi; max_sites = 16)
        edk = round(t, digits = 1)
        ed = haskey(ED, edk) ? @sprintf("%.5f", ED[edk]) : "-"
        @printf("%-5.2f %-11.6f %-9s %-9.2e %-9.1f  [minkept=%d nfb=%d finite=%s]\n",
                t, n, ed, info.maxcond, ewall, info.minkept, info.nfb, isfinite(n)); flush(stdout)
    end
    @printf("RSS peak = %.1f MB\n", Sys.maxrss()/2^20)
    println("DONE"); flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    tstart = parse(Float64, ARGS[1]); tend = parse(Float64, ARGS[2]); dt = parse(Float64, ARGS[3])
    env = Symbol(ARGS[4]); trunc_alg = Symbol(ARGS[5])
    tag = length(ARGS) >= 6 ? ARGS[6] : "run"
    baredt = length(ARGS) >= 7 ? parse(Float64, ARGS[7]) : 0.02
    patch = length(ARGS) >= 8 ? Symbol(ARGS[8]) : :torus
    D = length(ARGS) >= 9 ? parse(Int, ARGS[9]) : 2
    main(tstart, tend, dt, env, trunc_alg, tag, baredt, patch, D)
end
