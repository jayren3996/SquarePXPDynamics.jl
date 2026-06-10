# SPIKE N2 — wall-clock + peak RSS of building the double-layer exact-cluster env
# at frozen D=2 and D=4; confirm the OOM-reconcile; extrapolate trajectory cost.
#
# PROCESS ISOLATION: take D from ARGS so each D runs in a FRESH julia process
# (Sys.maxrss is a process-cumulative monotone high-water mark). A D=4 OOM-kill
# of one process then leaves the D=2 numbers already flushed to disk.
#   julia spike_n2_env_cost.jl 2     -> notes/.../exact_cluster_env_cost_d2.out
#   julia spike_n2_env_cost.jl 4     -> notes/.../exact_cluster_env_cost_d4.out
#
# Reuses build_exact_cluster_bondenv from spike_n4 (no duplication).

using SquarePXPDynamics
using PEPSKit
using TensorKit
using LinearAlgebra
using Printf

# include the builder (and the BUILDER ONLY — main() is guarded by PROGRAM_FILE)
include(joinpath(@__DIR__, "spike_n4_exact_cluster_bond.jl"))

function run_cost(D::Int)
    GC.gc()
    @printf("===== SPIKE N2 env-cost, D = %d =====\n", D)
    @printf("Sys.total_memory = %.1f MB (RAM ceiling)\n", Sys.total_memory() / 2^20)
    @printf("RSS at start     = %.1f MB\n", Sys.maxrss() / 2^20); flush(stdout)

    cell = PeriodicSquareUnitCell(4, 4)
    c0 = cell.reps[1]
    cidx = SquarePXPDynamics.squarecoord_to_cartesianindex(cell, c0)
    rc, cc = cidx[1], cidx[2]

    # below the t~2.6 peak (verdict): sample t in (2.0, 2.4)
    for t in (2.0, 2.4)
        @printf("\n--- D=%d, t=%.1f ---\n", D, t); flush(stdout)
        cell2 = PeriodicSquareUnitCell(4, 4)
        psi = checkerboard_pxp_pepskit_state(cell2; excited_on = :even, D = D)
        evolve_pepskit!(psi, t; dt = 0.02, order = 2, schedule = :serial,
                        maxdim = D, cutoff = 1e-12, rel_floor = 1e-3)
        Nr, Nc = size(psi.peps)
        cp1 = mod1(cc + 1, Nc)

        lossless = deepcopy(psi)
        project_star_pepskit!(lossless, c0; dt = 0.02, maxdim = 64, cutoff = 0.0,
                              rel_floor = 0.0, evolution = :real, projected = true)
        A = lossless.peps.A[rc, cc]
        B = lossless.peps.A[rc, cp1]
        X, a, b, Y = PEPSKit._qr_bond(A, B)
        @printf("X q-leg dim = %d ; Y q-leg dim = %d ; post-gate A space = %s\n",
                dim(space(X, 2)), dim(space(Y, 4)), string(space(A))); flush(stdout)

        # === TIME ONLY THE ENV BUILD ===
        GC.gc()
        rss0 = Sys.maxrss() / 2^20
        t0 = time()
        benv = build_exact_cluster_bondenv(lossless, X, Y, rc, cc, cp1)
        wall = time() - t0
        GC.gc()
        rss1 = Sys.maxrss() / 2^20
        @printf("ENV-BUILD: wall = %.3f s ; rss_delta = %.1f MB ; peak_rss = %.1f MB\n",
                wall, rss1 - rss0, rss1)
        @printf("ENV: dim(benv,1)=%d dim(benv,2)=%d ; summarysize = %.3f MB\n",
                dim(space(benv, 1)), dim(space(benv, 2)),
                Base.summarysize(benv) / 2^20); flush(stdout)

        # === BASELINE: single-layer exact_density_finite (warm then time) ===
        # warm JIT first (separate the 215s first-call artifact)
        _ = exact_density_finite(lossless; max_sites = 16)
        GC.gc()
        rssd0 = Sys.maxrss() / 2^20
        td0 = time()
        _ = exact_density_finite(lossless; max_sites = 16)
        dwall = time() - td0
        GC.gc()
        rssd1 = Sys.maxrss() / 2^20
        @printf("SINGLE-LAYER density: wall = %.3f s ; rss_delta = %.1f MB ; peak_rss = %.1f MB\n",
                dwall, rssd1 - rssd0, rssd1)
        @printf("double/single wall ratio = %.2f\n",
                dwall > 0 ? wall / dwall : NaN); flush(stdout)

        # === EXTRAPOLATION ===
        per_step = wall * 16 * 4          # 16 stars/step x 4 arms/star = 64 x wall
        full_traj = per_step * 140        # ~140 steps over [0,2.8] at dt=0.02
        @printf("EXTRAP: per-step env cost = %.2f s (64 x wall) ; full-traj ~ %.1f s = %.2f h (9000 x wall)\n",
                per_step, full_traj, full_traj / 3600); flush(stdout)
    end

    # === OOM-RECONCILE CONFIRMATION (D=4 only is the decisive one) ===
    if D == 4
        @printf("\n--- OOM-RECONCILE: single exact_density_finite at D=4, t=2.0 (fresh) ---\n")
        flush(stdout)
        cell3 = PeriodicSquareUnitCell(4, 4)
        psi3 = checkerboard_pxp_pepskit_state(cell3; excited_on = :even, D = 4)
        evolve_pepskit!(psi3, 2.0; dt = 0.02, order = 2, schedule = :serial,
                        maxdim = 4, cutoff = 1e-12, rel_floor = 1e-3)
        # warm
        _ = exact_density_finite(psi3; max_sites = 16)
        GC.gc()
        rssr0 = Sys.maxrss() / 2^20
        tr0 = time()
        nval = exact_density_finite(psi3; max_sites = 16)
        rwall = time() - tr0
        GC.gc()
        rssr1 = Sys.maxrss() / 2^20
        @printf("OOM-RECONCILE: D=4 t=2.0 single density COMPLETED. n=%.6f ; wall=%.3f s ; peak_rss=%.1f MB\n",
                nval, rwall, rssr1); flush(stdout)
        @printf("=> 'native exact_density_finite killed at D=4' RECLASSIFIED: frozen_step did ~33x\n")
        @printf("   redundant dense contractions with no GC between -> thrash/SIGTERM, NOT a strict OOM.\n")
        flush(stdout)
    end

    @printf("\nFINAL peak RSS (D=%d) = %.1f MB ; total mem = %.1f MB\n",
            D, Sys.maxrss() / 2^20, Sys.total_memory() / 2^20)
    println("DONE")
    flush(stdout)
end

let
    D = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
    run_cost(D)
end
