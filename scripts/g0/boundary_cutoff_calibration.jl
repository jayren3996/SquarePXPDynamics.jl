# Calibrate the :boundary ring-recompression cutoff against the exact dense
# oracle on 4x4, where both run — at the ENTANGLED times that killed the exact
# (cutoff=1e-14) ring (rank ~5e4, 74 GB, >1h at t>=0.5): which cutoff keeps the
# finite-torus boundary contraction BOTH accurate (vs dense) and fast enough to
# be the 6x6 production observable?
#
# Output rows flush loose-to-tight cutoff and early-to-late t, so the cheap
# verdicts land even if a tight-cutoff/late-t row grinds.

include(joinpath(@__DIR__, "native_to_legacy.jl"))
using Printf

const CUTOFFS = (1e-5, 1e-6, 1e-8, 1e-10)

function main()
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = 3)
    t = 0.0
    println("=== boundary-cutoff calibration: 4x4 D=3 rel_floor=0, threads=$(Threads.nthreads()) ===")
    @printf("%-5s %-11s %-12s %-12s %-10s %-8s\n", "t", "cutoff", "n_dense", "n_boundary", "delta", "wall_s")
    flush(stdout)
    for tstop in (0.5, 1.6, 2.6)
        while t < tstop - 1e-9
            evolve_pepskit!(psi, 0.1; dt = 0.02, order = 2, schedule = :serial,
                            maxdim = 3, cutoff = 1e-12, rel_floor = 0.0)
            t = round(t + 0.1, digits = 1)
        end
        leg = native_to_legacy(psi)
        nd = exact_density_finite(leg; max_sites = 16, method = :dense)
        for co in CUTOFFS
            el = @elapsed nb = exact_density_finite(leg; max_sites = 16, method = :boundary,
                                                    boundary_cutoff = co)
            @printf("%-5.1f %-11.0e %-12.8f %-12.8f %+-10.1e %-8.1f\n", t, co, nd, nb, nb - nd, el)
            flush(stdout)
        end
    end
    println("DONE")
end

main()
