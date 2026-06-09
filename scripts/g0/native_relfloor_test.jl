# Decisive test of the native D=2 rise divergence: is it the rel_floor drop-to-D=1?
# Both backends drop singular values below rel_floor*sigma_max (collapsing the bond
# dim). The TensorKit vs ITensors spectra differ, so native may trigger the drop on
# the rise where legacy does not. Compare rel_floor=1e-3 vs rel_floor=0, tracking the
# MINIMUM bond dimension (a drop to 1 = an entanglement channel destroyed) and the
# exact density, on a native D=2 trajectory.

using SquarePXPDynamics
using TensorKit
using Printf

# FIX 2026-06-07: the previous ED02 here was the TRUE 4x4 Neel ED SHIFTED +0.2 on the
# rise (ED02[t] == canon[t-0.2], e.g. old ED02[1.6]=0.10056 == canon[1.4]) — the same
# labeling bug fixed in native_d34_trajectory.jl / spike_dt_control.jl. It made the D=2
# rise look like a "blowup vs a flat ED ~0.10-0.15". Against the TRUE canonical ED
# (artifacts/neel_to_revival_4x4.json, 0.2 grid) the density genuinely RISES, and the
# D=2 EXCURSION is real anyway: D2 tracks ED to t=1.6 (0.145 vs 0.156) then overshoots
# to 0.428 at t=1.8 (true ED 0.244, +0.18) — a genuine D=2 instability that D>=3 removes.
const ED02 = Dict(1.0=>0.15549,1.2=>0.10057,1.4=>0.10056,1.6=>0.15558,1.8=>0.24370,
                  2.0=>0.33735,2.2=>0.41600)

function min_bond_dim(psi)
    m = typemax(Int)
    for T in psi.peps.A, leg in 2:5
        m = min(m, dim(space(T, leg)))
    end
    return m
end

function run(rel_floor; D = 2, tmax = 2.2, chunk = 0.1)
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = D)
    @printf("\n===== native D=%d, rel_floor=%.0e =====\n", D, rel_floor)
    @printf("%-5s %-10s %-9s %-8s %-12s\n", "t", "n_native", "n_ED", "min_bd", "chunk_trnerr")
    t = 0.0
    while t < tmax - 1e-9
        log = evolve_pepskit!(psi, chunk; dt = 0.02, order = 2, schedule = :serial,
                              maxdim = D, cutoff = 1e-12, rel_floor = rel_floor)
        t = round(t + chunk, digits = 1)
        n = exact_density_finite(psi; max_sites = 16)
        ed = haskey(ED02, t) ? @sprintf("%.5f", ED02[t]) : "-"
        @printf("%-5.1f %-10.6f %-9s %-8d %-12.4e\n", t, n, ed, min_bond_dim(psi), log.max_truncerr)
        flush(stdout)
    end
end

run(1e-3)
run(0.0)
println("\nDONE")
