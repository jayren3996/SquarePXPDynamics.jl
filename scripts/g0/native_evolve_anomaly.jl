# Pinpoint the native evolve_pepskit! divergence on the revival rise.
# One continuous native trajectory per D, sampled every 0.1, printing the exact
# (CTM-free) density, the per-chunk max truncation error, and the smallest bond
# weight (near-degeneracy probe). Compare D=2 (anomalous) vs D=4.
# ED on the 0.2 grid for reference; legacy D=2 from the trajectory probe.

using SquarePXPDynamics
using Printf

# ED reference (0.2 grid) and legacy D=2 trajectory (from scripts/g0/probe_d2)
const ED02 = Dict(
    1.0=>0.15298, 1.2=>0.09521, 1.4=>0.08779, 1.6=>0.10056, 1.8=>0.15558,
    2.0=>0.24370, 2.2=>0.33735, 2.4=>0.41600, 2.6=>0.46696, 2.8=>0.48306)
const LEG2 = Dict(  # legacy native-equivalent simple update, D=2
    1.0=>0.15298, 1.2=>0.09521, 1.4=>0.08779, 1.6=>0.13698, 1.8=>0.20489,
    2.0=>0.28130, 2.2=>0.36063, 2.4=>0.42987, 2.6=>0.47659, 2.8=>0.48685)

function trajectory(D::Int; tmax = 2.2, chunk = 0.1)
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = D)
    @printf("\n===== native D=%d, dt=0.02, order=2, :serial, rel_floor=1e-3 =====\n", D)
    @printf("%-5s %-10s %-10s %-10s %-12s\n", "t", "n_native", "n_legacy2", "n_ED", "chunk_trnerr")
    @printf("%-5.1f %-10.6f %-10s %-10s %-12s\n", 0.0, exact_density_finite(psi; max_sites=16), "-", "-", "-")
    t = 0.0
    while t < tmax - 1e-9
        log = evolve_pepskit!(psi, chunk; dt = 0.02, order = 2, schedule = :serial,
                              maxdim = D, cutoff = 1e-12, rel_floor = 1e-3)
        t += chunk
        n = exact_density_finite(psi; max_sites = 16)
        tk = round(t, digits=1)
        leg = haskey(LEG2, tk) && D == 2 ? @sprintf("%.6f", LEG2[tk]) : "-"
        ed  = haskey(ED02, tk) ? @sprintf("%.6f", ED02[tk]) : "-"
        @printf("%-5.1f %-10.6f %-10s %-10s %-12.4e\n", tk, n, leg, ed, log.max_truncerr)
        flush(stdout)
    end
end

trajectory(2)
trajectory(4)
println("\nDONE")
