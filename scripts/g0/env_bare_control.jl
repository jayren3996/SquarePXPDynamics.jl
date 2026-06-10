# env_bare_control.jl — adversarial CAUSAL-ATTRIBUTION control for the onset probe.
#
# The onset probe showed env=:exact_cluster D=2 suppresses the t=1.7-1.9 blowup
# (0.43 -> 0.22, onto the D=3/D=4 curve). But project_star_env! GROWS all bonds
# lossless THEN re-truncates -- so the fix could be the loop-carrying ENVIRONMENT,
# OR merely that lossless-grow-then-retruncate procedure (which differs from the
# native inline-truncate). This control isolates the two:
#
#   - native-bare : project_star_pepskit! (inline truncate, D=2 throughout) -> spikes (known)
#   - env=:bare   : project_star_env! lossless-grow THEN bare-SVD truncate (NO env)
#   - [reference] : env=:exact_cluster from onset_intervention_probe.out (0.223 @ t=1.8)
#
# VERDICT: if env=:bare spikes ~like native-bare while :exact_cluster stays smooth
# (0.22), the loop-carrying ENVIRONMENT is decisively the cause (the procedure is
# not). If env=:bare is itself smooth, the procedure carries (some of) the fix.

using SquarePXPDynamics
using TensorKit
using Printf

include(joinpath(@__DIR__, "star_env.jl"))  # -> project_star_env! (now with env=:bare)

const ED  = Dict(1.6=>0.15558, 1.8=>0.24370, 2.0=>0.33735)
const D34 = Dict(1.6=>0.14615, 1.8=>0.23363, 2.0=>0.32746)
const EXC = Dict(1.6=>0.14480, 1.7=>0.18052, 1.8=>0.22321, 1.9=>0.26977, 2.0=>0.31652)  # :exact_cluster (probe)

function evolve_env_serial!(state, total_time; dt, maxdim, cutoff, env, patch, trunc_alg)
    reps = state.unitcell.reps
    h = dt / 2
    subs = vcat([(c, h) for c in reps], [(c, h) for c in Base.reverse(reps)])
    nsteps = round(Int, total_time / dt)
    minkept = typemax(Int)
    for _ in 1:nsteps
        for (center, layer_dt) in subs
            info = project_star_env!(state, center; dt = layer_dt, maxdim = maxdim,
                                     cutoff = cutoff, env = env, patch = patch,
                                     trunc_alg = trunc_alg, evolution = :real, projected = true)
            minkept = min(minkept, minimum(values(info.keptdims)))
        end
    end
    return minkept
end

bare_chunk!(state, chunk; dt, maxdim) =
    evolve_pepskit!(state, chunk; dt = dt, order = 2, schedule = :serial,
                    maxdim = maxdim, cutoff = 1e-12, rel_floor = 0.0)

function main()
    D = 2; dt = 0.02; chunk = 0.1
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = D)
    _ = exact_density_finite(psi; max_sites = 16)         # warm JIT

    t = 0.0
    while t < 1.6 - 1e-9
        bare_chunk!(psi, chunk; dt = dt, maxdim = D); t = round(t + chunk, digits = 1)
    end
    @printf("evolved bare to t=1.6 ; n=%.6f\n", exact_density_finite(psi; max_sites=16)); flush(stdout)

    psi_native = deepcopy(psi)
    psi_envbar = deepcopy(psi)

    println("\n=== CONTROL: native-bare vs env=:bare (lossless->bareSVD) from t=1.6 ===")
    @printf("%-5s %-12s %-12s %-12s %-9s %-9s  %s\n",
            "t", "native_bare", "env_bare", "exact_clus", "D3D4", "ED", "note")
    @printf("%-5.1f %-12.6f %-12.6f %-12.6f %-9.5f %-9.5f\n",
            1.6, exact_density_finite(psi_native;max_sites=16),
            exact_density_finite(psi_envbar;max_sites=16), EXC[1.6], D34[1.6], ED[1.6]); flush(stdout)

    t = 1.6
    while t < 2.0 - 1e-9
        bare_chunk!(psi_native, chunk; dt = dt, maxdim = D)
        mk = evolve_env_serial!(psi_envbar, chunk; dt = dt, maxdim = D, cutoff = 1e-12,
                                env = :bare, patch = :torus, trunc_alg = :full)
        t = round(t + chunk, digits = 1)
        nn = exact_density_finite(psi_native; max_sites = 16)
        ne = exact_density_finite(psi_envbar; max_sites = 16)
        exc = haskey(EXC, t) ? @sprintf("%.6f", EXC[t]) : "-"
        d3  = haskey(D34, t) ? @sprintf("%.5f", D34[t]) : "-"
        ed  = haskey(ED, t)  ? @sprintf("%.5f", ED[t])  : "-"
        @printf("%-5.1f %-12.6f %-12.6f %-12s %-9s %-9s  [minkept=%d]\n",
                t, nn, ne, exc, d3, ed, mk); flush(stdout)
    end

    println("\n===== VERDICT =====")
    println("env=:bare ~ native-bare (spikes to ~0.43 @ t=1.8)  => the lossless-then-retruncate")
    println("   PROCEDURE does NOT fix it => the loop-carrying ENVIRONMENT is the cause.")
    println("env=:bare ~ exact_cluster (smooth ~0.22 @ t=1.8)   => the PROCEDURE carries the fix")
    println("   (env-quality would then be unproven by this experiment).")
    println("DONE"); flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
