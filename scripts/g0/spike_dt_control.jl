# SPIKE dt-control — attribute (or rule out) the native-vs-ED rise DISCREPANCY
# as a Trotter dt artifact. SCRIPT-ONLY this round (long native trajectory; run
# later in the main loop). NO env/truncation machinery — pure native simple-update.
#
# native_d34_trajectory.out shows native D=3 LEADS ED on the rise (not lags):
# t=2.0 native=0.327 vs ED=0.244. The metric is well-defined either way: is the
# native-vs-ED rise DISCREPANCY dt-STABLE?
#   - dt-STABLE (|dRMS| and max rise dt_drift << 1e-2) -> discrepancy is the
#     environment/representation, NOT Trotter -> SUPPORTS the exact-cluster engine.
#   - dt-SENSITIVE (drift comparable to ~1e-2) -> Trotter error contaminates the
#     signal -> defer exact-cluster pending a dt-converged baseline.

using SquarePXPDynamics
using Printf

# TRUE 4x4 Neel ED reference (0.2 grid), from artifacts/neel_to_revival_4x4.json
# (== test_d_convergence.jl ED array).
# FIX 2026-06-07: the previous ED02 was the TRUE ED shifted +0.2 on the rise
# (ED02[t]==canon[t-0.2] for t>=1.6) plus a spliced D=4-iPEPS collapse branch — a
# labeling bug that fabricated a spurious "native LEADS ED by ~0.046" result.
const ED02 = Dict(1.0=>0.15549264928630516, 1.2=>0.10057443671942552, 1.4=>0.10056019790056131,
                  1.6=>0.1555774954125962, 1.8=>0.24370069054283897, 2.0=>0.3373453880372616,
                  2.2=>0.4159980265031265, 2.4=>0.46695735678274974, 2.6=>0.4830586210094464,
                  2.8=>0.46229670220636837)

const DMEAS = 3        # native ladder flat at D>=3 (D3==D4 to 6 digits); D=2 blows up, excluded
const CHUNK = 0.2      # multiple of BOTH 0.01 and 0.02 -> Trotter grid lands on 0.2 samples
const TMAX = 2.8

function run_dt(DT::Float64)
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = DMEAS)
    @printf("\n===== dt = %.3f (D = %d) =====\n", DT, DMEAS)
    @printf("%-5s %-10s %-9s\n", "t", "n_native", "n_ED"); flush(stdout)
    # JIT warmup density call (separate the ~215s first-call compile from the trajectory)
    _ = exact_density_finite(psi; max_sites = 16)
    t = 0.0
    n_of_t = Dict{Float64,Float64}()
    while t < TMAX - 1e-9
        evolve_pepskit!(psi, CHUNK; dt = DT, order = 2, schedule = :serial,
                        maxdim = DMEAS, cutoff = 1e-12, rel_floor = 1e-3)
        t = round(t + CHUNK, digits = 1)
        n = exact_density_finite(psi; max_sites = 16)
        n_of_t[t] = n
        ed = haskey(ED02, t) ? @sprintf("%.5f", ED02[t]) : "-"
        @printf("%-5.1f %-10.6f %-9s\n", t, n, ed); flush(stdout)
    end
    return n_of_t
end

function main()
    @printf("RSS at start = %.1f MB ; total mem = %.1f MB\n",
            Sys.maxrss() / 2^20, Sys.total_memory() / 2^20); flush(stdout)

    n001 = run_dt(0.01)
    n002 = run_dt(0.02)

    # === METRICS on the shared 0.2 grid, using ED02 where present ===
    grid = sort(collect(keys(ED02)))
    rms(d) = sqrt(sum((d[t] - ED02[t])^2 for t in grid if haskey(d, t)) /
                  count(t -> haskey(d, t), grid))
    rms001 = rms(n001)
    rms002 = rms(n002)

    println("\n===== PER-TIME TABLE (rise focus t in [1.6,2.6]) =====")
    @printf("%-5s %-12s %-12s %-10s %-12s\n",
            "t", "n_dt001", "n_dt002", "n_ED", "dt_drift")
    rise_drifts = Float64[]
    for t in grid
        haskey(n001, t) && haskey(n002, t) || continue
        drift = n002[t] - n001[t]
        ed = haskey(ED02, t) ? @sprintf("%.5f", ED02[t]) : "-"
        @printf("%-5.1f %-12.6f %-12.6f %-10s %-12.4e\n",
                t, n001[t], n002[t], ed, drift)
        if 1.6 <= t <= 2.6
            push!(rise_drifts, abs(drift))
        end
    end
    flush(stdout)

    max_rise_drift = isempty(rise_drifts) ? NaN : maximum(rise_drifts)
    drms = abs(rms002 - rms001)

    println("\n===== TROTTER-ATTRIBUTION =====")
    @printf("RMS_dt001 (vs ED)        = %.4e\n", rms001)
    @printf("RMS_dt002 (vs ED)        = %.4e\n", rms002)
    @printf("|RMS_dt002 - RMS_dt001|  = %.4e\n", drms)
    @printf("max rise dt_drift |t in [1.6,2.6]| = %.4e\n", max_rise_drift)
    @printf("lag scale (state-doc)    ~ 1e-2\n")
    if drms < 1e-3 && max_rise_drift < 1e-3
        println("VERDICT: discrepancy is dt-STABLE (<< 1e-2) -> NOT a Trotter artifact ->")
        println("         environment/representation -> SUPPORTS exact-cluster engine.")
    elseif max_rise_drift >= 5e-3
        println("VERDICT: discrepancy is dt-SENSITIVE (drift ~ lag) -> Trotter contamination ->")
        println("         DEFER exact-cluster pending a dt-converged baseline.")
    else
        println("VERDICT: intermediate — drift is non-negligible but below the lag scale;")
        println("         report numbers, lean dt-stable but note residual Trotter sensitivity.")
    end
    @printf("RSS peak = %.1f MB\n", Sys.maxrss() / 2^20)
    println("DONE")
    flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
