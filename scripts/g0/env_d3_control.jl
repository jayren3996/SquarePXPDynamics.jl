# env_d3_control.jl — the DECISIVE bond-dimension test for the loop env.
#
# The dt-convergence result (env_dtconv_dt01_to1.8.out) pinned env=:exact_cluster
# at D=2 to ~0.0205 BELOW true ED at t=1.8, dt-converged (gap = D=2 method error,
# NOT Trotter). The review (wbmdbcthf, unanimous) showed the env is ALREADY the
# maximal exact torus env, so the residual is D-BOUNDED, not patch-bounded: the
# accuracy lever is bond dimension D, not dt or env-quality.
#
# This run measures the D=3 point of that lever. From the SAME t=1.6 fork (bare,
# rel_floor=0, but now at maxdim=3), evolve TWO branches to t=2.0 on the 0.1 grid:
#   - native-bare D=3 : project_star_pepskit! (inline bare-SVD, maxdim=3)
#   - env D=3         : project_star_env! env=:exact_cluster, trunc_alg=:full, maxdim=3
# vs true ED. Decision:
#   * env-D3 climbs ABOVE env-D2's 0.2232 toward ED 0.2437 (and >= native-D3)
#       => the loop env's D-lever works; raising D closes the gap to physics.
#   * env-D3 BEATS native-D3 (closer to ED)
#       => effect B alive: a loop env beats the simple-update env at fixed D.
#   * env-D3 ~ native-D3 ~ saturates at 0.223
#       => effect-A-only: the env is a stable integrator, not a physics engine,
#          and the D-lever does not climb (a deeper mean-field ceiling).
#
# NOTE on the native-D3 anchor: native_d34_trajectory used rel_floor=1e-3 (its
# D3==D4 "flat ladder" may be a floor artifact). This run regenerates native-D3
# at rel_floor=0 in-line so the env-vs-native comparison is apples-to-apples.

using SquarePXPDynamics
using TensorKit
using Printf

include(joinpath(@__DIR__, "star_env.jl"))   # -> project_star_env!

# canonical same-t ED (artifacts/neel_to_revival_4x4.json, 0.2 grid), verified.
const ED  = Dict(1.6=>0.15558, 1.8=>0.24370, 2.0=>0.33735)
# env D=2 reference (env_dtconv / onset probe, exact_cluster:full, raw 6-digit)
const ED2_ENV = Dict(1.6=>0.144802, 1.8=>0.223205, 2.0=>0.316522)

function evolve_env_serial!(state, total_time; dt, maxdim, cutoff, env, patch, trunc_alg)
    reps = state.unitcell.reps
    h = dt / 2
    subs = vcat([(c, h) for c in reps], [(c, h) for c in Base.reverse(reps)])
    nsteps = round(Int, total_time / dt)
    minkept = typemax(Int); nfb = 0; maxcond = 0.0
    for _ in 1:nsteps
        for (center, layer_dt) in subs
            info = project_star_env!(state, center; dt = layer_dt, maxdim = maxdim,
                                     cutoff = cutoff, env = env, patch = patch,
                                     trunc_alg = trunc_alg, evolution = :real, projected = true)
            minkept = min(minkept, minimum(values(info.keptdims)))
            nfb += count(values(info.fellbacks))
            maxcond = max(maxcond, maximum(values(info.conds)))
        end
    end
    return (minkept = minkept, nfb = nfb, maxcond = maxcond)
end

bare_chunk!(state, chunk; dt, maxdim) =
    evolve_pepskit!(state, chunk; dt = dt, order = 2, schedule = :serial,
                    maxdim = maxdim, cutoff = 1e-12, rel_floor = 0.0)

function main(tend)
    D = 3; dt = 0.02; chunk = 0.1
    @printf("RSS start = %.1f MB ; total = %.1f MB\n", Sys.maxrss()/2^20, Sys.total_memory()/2^20)
    @printf("=== env_d3_control : bare-D3 vs env-D3(:exact_cluster:full) fork@1.6 -> %.2f, dt=%.3f ===\n",
            tend, dt); flush(stdout)
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = D)
    _ = exact_density_finite(psi; max_sites = 16)        # warm JIT

    t = 0.0
    while t < 1.6 - 1e-9
        bare_chunk!(psi, chunk; dt = dt, maxdim = D); t = round(t + chunk, digits = 1)
    end
    @printf("bare-D3 evolved to fork t=1.6 ; n=%.6f (env-D2 ref %.6f, ED %.5f)\n",
            exact_density_finite(psi; max_sites=16), ED2_ENV[1.6], ED[1.6]); flush(stdout)

    # warm the env-build JIT once at D=3 (absorb cold ncon compile at q-legs 24/12)
    tw0 = time()
    let warm = deepcopy(psi)
        project_star_env!(warm, cell.reps[1]; dt = dt/2, maxdim = D, cutoff = 1e-12,
                          env = :exact_cluster, patch = :torus, trunc_alg = :full,
                          evolution = :real, projected = true)
    end
    @printf("env-D3 JIT warmup = %.1f s\n", time() - tw0); flush(stdout)

    psi_bare = deepcopy(psi)
    psi_env  = deepcopy(psi)

    @printf("\n%-5s %-12s %-12s %-12s %-9s %-9s  %s\n",
            "t", "bare_D3", "env_D3", "env_D2", "ED", "wall_s", "diag")
    @printf("%-5.2f %-12.6f %-12.6f %-12.6f %-9.5f %-9s\n",
            1.6, exact_density_finite(psi_bare;max_sites=16),
            exact_density_finite(psi_env;max_sites=16), ED2_ENV[1.6], ED[1.6], "-"); flush(stdout)

    t = 1.6
    while t < tend - 1e-9
        bare_chunk!(psi_bare, chunk; dt = dt, maxdim = D)
        te0 = time()
        info = evolve_env_serial!(psi_env, chunk; dt = dt, maxdim = D, cutoff = 1e-12,
                                  env = :exact_cluster, patch = :torus, trunc_alg = :full)
        ewall = time() - te0
        t = round(t + chunk, digits = 1)
        nb = exact_density_finite(psi_bare; max_sites = 16)
        ne = exact_density_finite(psi_env;  max_sites = 16)
        e2 = haskey(ED2_ENV, t) ? @sprintf("%.6f", ED2_ENV[t]) : "-"
        ed = haskey(ED, t)      ? @sprintf("%.5f", ED[t])      : "-"
        @printf("%-5.2f %-12.6f %-12.6f %-12s %-9s %-9.1f  [minkept=%d nfb=%d cond=%.1e fin=%s]\n",
                t, nb, ne, e2, ed, ewall, info.minkept, info.nfb, info.maxcond,
                isfinite(ne)); flush(stdout)
    end
    @printf("RSS peak = %.1f MB\n", Sys.maxrss()/2^20)
    println("DONE"); flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    tend = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 2.0
    main(tend)
end
