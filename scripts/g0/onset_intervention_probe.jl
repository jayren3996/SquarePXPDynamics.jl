# onset_intervention_probe.jl — the CHEAP decisive-direction gate before any 6-12h
# from-t=0 env trajectory.
#
# PREMISE (re-grounded against the TRUE canonical ED, artifacts/neel_to_revival_4x4.json):
#   native D=2 simple-update TRACKS the true ED to t=1.6 (n=0.145 vs ED 0.156), then
#   has a GENUINE blowup: t=1.7=0.403, t=1.8=0.428 (true ED only 0.244 -> +0.18 overshoot),
#   with min_bd=2 throughout (rel_floor=0, NOT a D=1 drop). D=3==D=4 REMOVE the blowup
#   (smooth, ~1% undershoot of ED). So the D=2 instability is bond-sensitive 2->3 and is a
#   TRUNCATION-ENVIRONMENT effect, not a truncation-magnitude effect.
#
# CONTROLLED EXPERIMENT: evolve bare D=2 to t=1.6 (the disease onset, just before the jump),
# then FORK and continue with the SAME :serial order-2 schedule, swapping ONLY the truncation
# engine:
#   - bare branch  : project_star_pepskit! (bare SVD)        -> should reproduce the spike
#   - env  branch  : project_star_env!(:exact_cluster,:full) -> loop-carrying {2,2} environment
# Holding the t=1.6 state fixed isolates the truncation-ENVIRONMENT effect at the onset
# (a cleaner control than from-t=0, which confounds accumulated history). It ALSO exercises
# fixgauge_benv on the REAL rise (cond(Z'Z)~1e20) -- the regime every frozen-t=1.8 gate left
# unproven.
#
# PASS (direction): at t=1.7/1.8 the env branch stays near the ED/D=3 level (~0.16/0.24)
# instead of the bare 0.40/0.43 spike, AND every env arm is finite (no NaN, fixgauge robust).
# -> then commit to the full from-t=0 trajectory for the publication curve.

using SquarePXPDynamics
using TensorKit
using Printf

include(joinpath(@__DIR__, "star_env.jl"))  # -> project_star_env!, build_star_bond_env

# TRUE canonical ED (0.2 grid, artifacts/neel_to_revival_4x4.json) and the D=3==D=4
# native reference (native_d34_trajectory.out, label-agnostic physics: smooth, no spike).
const ED  = Dict(1.6=>0.15558, 1.8=>0.24370, 2.0=>0.33735)
const D34 = Dict(1.6=>0.14615, 1.8=>0.23363, 2.0=>0.32746)

# --- env-aware analogue of evolve_pepskit!(...; schedule=:serial, order=2) ---
# Substep-for-substep identical to _step_substeps(:serial, order=2):
#   [(rep, dt/2) for rep in reps] ++ [(rep, dt/2) for rep in reverse(reps)]
# only the per-star kernel is swapped (project_star_env! vs project_star_pepskit!).
function evolve_env_serial!(state, total_time; dt, maxdim, cutoff,
                            env, patch, trunc_alg)
    reps = state.unitcell.reps
    h = dt / 2
    subs = vcat([(c, h) for c in reps], [(c, h) for c in Base.reverse(reps)])
    nsteps = round(Int, total_time / dt)
    maxcond = 0.0; nfb = 0; minkept = typemax(Int)
    for _ in 1:nsteps
        for (center, layer_dt) in subs
            info = project_star_env!(state, center; dt = layer_dt, maxdim = maxdim,
                                     cutoff = cutoff, env = env, patch = patch,
                                     trunc_alg = trunc_alg, evolution = :real,
                                     projected = true)
            maxcond = max(maxcond, maximum(values(info.conds)))
            nfb += count(values(info.fellbacks))
            minkept = min(minkept, minimum(values(info.keptdims)))
        end
    end
    return (maxcond = maxcond, nfb = nfb, minkept = minkept)
end

function bare_chunk!(state, chunk; dt, maxdim)
    evolve_pepskit!(state, chunk; dt = dt, order = 2, schedule = :serial,
                    maxdim = maxdim, cutoff = 1e-12, rel_floor = 0.0)
end

function main()
    @printf("RSS start = %.1f MB ; total = %.1f MB\n",
            Sys.maxrss()/2^20, Sys.total_memory()/2^20); flush(stdout)
    D = 2; dt = 0.02; chunk = 0.1
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = D)

    # warm exact_density JIT (separate the ~215s first-call compile from the trajectory)
    _ = exact_density_finite(psi; max_sites = 16)

    # --- evolve bare to the disease onset t=1.6 (rel_floor=0 -> min_bd=2 clean) ---
    println("\n=== bare D=2 to onset t=1.6 (rel_floor=0) ==="); flush(stdout)
    @printf("%-5s %-10s %-9s\n", "t", "n_bare", "n_ED")
    t = 0.0
    while t < 1.6 - 1e-9
        bare_chunk!(psi, chunk; dt = dt, maxdim = D)
        t = round(t + chunk, digits = 1)
        n = exact_density_finite(psi; max_sites = 16)
        ed = haskey(ED, t) ? @sprintf("%.5f", ED[t]) : "-"
        @printf("%-5.1f %-10.6f %-9s\n", t, n, ed); flush(stdout)
    end

    # --- WARM the env-build JIT once on a throwaway copy (absorb ~127s cold ncon) ---
    println("\n--- warming env-build JIT (one throwaway project_star_env!) ---"); flush(stdout)
    tw0 = time()
    let warm = deepcopy(psi)
        project_star_env!(warm, cell.reps[1]; dt = dt/2, maxdim = D, cutoff = 1e-12,
                          env = :exact_cluster, patch = :torus, trunc_alg = :full,
                          evolution = :real, projected = true)
    end
    @printf("warmup wall = %.1f s\n", time() - tw0); flush(stdout)

    # --- FORK at t=1.6 ---
    psi_bare = deepcopy(psi)
    psi_env  = deepcopy(psi)

    println("\n=== ONSET INTERVENTION: bare vs env (:exact_cluster,:full) from t=1.6 ===")
    @printf("%-5s %-11s %-11s %-9s %-9s %-9s %-9s\n",
            "t", "n_bare", "n_env", "n_ED", "n_D3D4", "maxcond", "wall_s")
    @printf("%-5.1f %-11.6f %-11.6f %-9.5f %-9.5f %-9s %-9s\n",
            1.6, exact_density_finite(psi_bare; max_sites=16),
            exact_density_finite(psi_env; max_sites=16), ED[1.6], D34[1.6], "-", "-")
    flush(stdout)

    t = 1.6
    while t < 2.0 - 1e-9
        bare_chunk!(psi_bare, chunk; dt = dt, maxdim = D)
        te0 = time()
        envinfo = evolve_env_serial!(psi_env, chunk; dt = dt, maxdim = D, cutoff = 1e-12,
                                     env = :exact_cluster, patch = :torus, trunc_alg = :full)
        ewall = time() - te0
        t = round(t + chunk, digits = 1)
        nb = exact_density_finite(psi_bare; max_sites = 16)
        ne = exact_density_finite(psi_env;  max_sites = 16)
        ed = haskey(ED, t) ? @sprintf("%.5f", ED[t]) : "-"
        d3 = haskey(D34, t) ? @sprintf("%.5f", D34[t]) : "-"
        finite_ok = isfinite(nb) && isfinite(ne)
        @printf("%-5.1f %-11.6f %-11.6f %-9s %-9s %-9.2e %-9.1f  [minkept=%d nfb=%d finite=%s]\n",
                t, nb, ne, ed, d3, envinfo.maxcond, ewall,
                envinfo.minkept, envinfo.nfb, finite_ok); flush(stdout)
    end

    # --- VERDICT at the spike anchor t=1.8 ---
    println("\n===== VERDICT (spike anchor t=1.8) =====")
    nb18 = exact_density_finite(psi_bare; max_sites = 16)  # note: psi_bare now at t=2.0
    println("Read the t=1.8 row above:")
    println("  bare ~0.43 (the +0.18 overshoot vs true ED 0.244) is the disease.")
    println("  D=3/D=4 = 0.234 (ED-tracking) is the cure-by-bond-dim TARGET.")
    println("  env near 0.234 / 0.244 => loop-carrying environment substitutes for bond dim")
    println("                            => engine works => commit to full from-t=0 trajectory.")
    println("  env still ~0.43        => environment does NOT fix the D=2 instability.")
    @printf("RSS peak = %.1f MB\n", Sys.maxrss()/2^20)
    println("DONE"); flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
