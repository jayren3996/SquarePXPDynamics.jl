# Disambiguate the native rise divergence: does it persist at higher bond dim?
# Native D=3 and D=4 trajectories vs ED and vs the legacy D=4 probe.
# If native D=4 tracks ED (like legacy D=4, RMS~9.8e-3) -> the D=2 blowup is
# truncation-crudeness/sensitivity (CTM-aware truncation is the candidate fix).
# If native D=4 ALSO blows up on the rise -> a deeper native bug to fix first.

using SquarePXPDynamics
using TensorKit
using Printf

# FIX 2026-06-07: previous ED02 was the TRUE ED shifted +0.2 on the rise (canon[t-0.2])
# + a spliced D=4-iPEPS collapse branch — a labeling bug. Now the true 4x4 Neel ED
# from artifacts/neel_to_revival_4x4.json (== test_d_convergence.jl).
const ED02 = Dict(1.0=>0.15549264928630516,1.2=>0.10057443671942552,1.4=>0.10056019790056131,
                  1.6=>0.1555774954125962,1.8=>0.24370069054283897,2.0=>0.3373453880372616,
                  2.2=>0.4159980265031265,2.4=>0.46695735678274974,2.6=>0.4830586210094464,
                  2.8=>0.46229670220636837)
const LEG4 = Dict(1.6=>0.14286,1.8=>0.23237,2.0=>0.32598,2.2=>0.40511,2.4=>0.46426,
                  2.6=>0.47285,2.8=>0.43770)  # legacy D=4 probe

function run(D::Int; tmax = 2.4, chunk = 0.2)
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = D)
    @printf("\n===== native D=%d =====\n", D)
    @printf("%-5s %-10s %-9s %-9s %-12s\n", "t", "n_native", "n_ED", "n_leg4", "chunk_trnerr")
    t = 0.0
    while t < tmax - 1e-9
        log = evolve_pepskit!(psi, chunk; dt = 0.02, order = 2, schedule = :serial,
                              maxdim = D, cutoff = 1e-12, rel_floor = 1e-3)
        t = round(t + chunk, digits = 1)
        n = exact_density_finite(psi; max_sites = 16)
        ed = haskey(ED02, t) ? @sprintf("%.5f", ED02[t]) : "-"
        l4 = (D == 4 && haskey(LEG4, t)) ? @sprintf("%.5f", LEG4[t]) : "-"
        @printf("%-5.1f %-10.6f %-9s %-9s %-12.4e\n", t, n, ed, l4, log.max_truncerr)
        flush(stdout)
    end
end

run(3; tmax = 2.8, chunk = 0.1)
run(4; tmax = 2.8, chunk = 0.1)
println("\nDONE")
