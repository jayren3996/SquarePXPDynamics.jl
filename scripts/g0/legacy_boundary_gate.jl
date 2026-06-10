# Gate the native->legacy adapter + legacy :boundary oracle, then time it at 6x6.
#
# G1 (4x4, value): along a native D=3 rel_floor=0 trajectory, at t=0.5/1.0/2.6
#     n_native(:dense) == n_legacy(:dense) == n_legacy(:boundary) to <=1e-9.
#     This gates the adapter (conventions) AND the boundary sweep (algorithm)
#     against the established native oracle in collapse + entangled regimes.
# G2 (6x6, value+cost): evolve 6x6 D=3 rel_floor=0 to t=0.2, convert, run
#     :boundary at 36 sites. Short-time n(t) is ED-grade checkable: Krylov 6x6
#     gives n(0.2)=0.4802653 (and 4x4/6x6 agree to 1e-10 there). Wall time is
#     the per-point cost of the production 6x6 trajectory observable.
#
# REORDERED 2026-06-10 (gate-runner; tests/params/tolerances/row formats are
# UNCHANGED — see legacy_boundary_gate.jl.orig). The original interleaved the
# legacy :boundary call into each G1 row before anything was flushed; the first
# :boundary call (4x4, D=3, t=0.5, entangled) rank-explodes in the exact
# cutoff=1e-14 ring recompression (~48600^2 zgesdd, ~74 GB RSS, >1 h, no
# output). Now every cheap phase runs and flushes first, and the entangled
# 4x4 :boundary attempts come LAST:
#   phase A: 4x4 trajectory; native-dense vs legacy-dense at t=0.5/1.0/2.6
#            (adapter/convention gate), converted states stashed
#   phase B: original G2, 6x6 t=0.2 convert + 36-site :boundary value/cost
#   phase C: original G1 rows incl. legacy :boundary on the stashed states

include(joinpath(@__DIR__, "native_to_legacy.jl"))
using Printf

println("gate script start, threads=", Threads.nthreads())
flush(stdout)

# --- phase A: 4x4 trajectory + dense adapter check; stash converted states ---
function phase_a()
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = 3)
    t = 0.0
    stash = Tuple{Float64,Float64,Float64,SquareIPEPSState}[]
    println("=== G1 phase A: 4x4 native dense vs legacy dense (adapter gate) ===")
    flush(stdout)
    for tstop in (0.5, 1.0, 2.6)
        while t < tstop - 1e-9
            evolve_pepskit!(psi, 0.1; dt = 0.02, order = 2, schedule = :serial,
                            maxdim = 3, cutoff = 1e-12, rel_floor = 0.0)
            t = round(t + 0.1, digits = 1)
        end
        nnat = exact_density_finite(psi; max_sites = 16)
        leg = native_to_legacy(psi)
        nld = exact_density_finite(leg; max_sites = 16, method = :dense)
        @printf("phase A t=%-4.1f n_nat_dense=%.8f n_leg_dense=%.8f d_dense=%.1e\n",
                t, nnat, nld, nld - nnat)
        flush(stdout)
        push!(stash, (t, nnat, nld, leg))
    end
    return stash
end

# --- phase B: original G2 (unchanged tests and print format) ---
function probe_6x6()
    println("\n=== G2: 6x6 D=3 :boundary value + cost (Krylov n(0.2)=0.4802653) ===")
    flush(stdout)
    cell = PeriodicSquareUnitCell(6, 6)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = 3)
    evolve_pepskit!(psi, 0.2; dt = 0.02, order = 2, schedule = :serial,
                    maxdim = 3, cutoff = 1e-12, rel_floor = 0.0)
    println("6x6 evolved to t=0.2; converting, then 36-site :boundary")
    flush(stdout)
    el_conv = @elapsed leg = native_to_legacy(psi)
    el = @elapsed n = exact_density_finite(leg; max_sites = 36, method = :boundary)
    @printf("t=0.2: n_boundary=%.8f (Krylov 0.48026525)  convert=%.1f s  boundary=%.1f s (%d threads)\n",
            n, el_conv, el, Threads.nthreads())
    flush(stdout)
end

# --- phase C: original G1 rows, legacy :boundary on the stashed states ---
function phase_c(stash)
    println("\n=== G1: 4x4 adapter + :boundary vs native dense ===")
    @printf("%-5s %-12s %-12s %-12s %-9s %-9s\n",
            "t", "n_nat_dense", "n_leg_dense", "n_leg_bnd", "d_dense", "d_bnd")
    flush(stdout)
    for (t, nnat, nld, leg) in stash
        @printf("starting 4x4 :boundary at t=%.1f\n", t)
        flush(stdout)
        el = @elapsed nlb = exact_density_finite(leg; max_sites = 16, method = :boundary)
        @printf("%-5.1f %-12.8f %-12.8f %-12.8f %-9.1e %-9.1e (bnd %.1f s)\n",
                t, nnat, nld, nlb, nld - nnat, nlb - nnat)
        flush(stdout)
    end
end

stash4 = phase_a()
probe_6x6()
phase_c(stash4)
println("DONE")
