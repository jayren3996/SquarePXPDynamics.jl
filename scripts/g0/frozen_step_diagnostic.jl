# G0 — frozen-state per-step truncation-cost diagnostic (native backend).
#
# Question PR-0/G0 must settle: is the revival-rise lag a per-step TRUNCATION-METRIC
# error (=> a better bond metric, "bare Stage C", could close it) or a structural
# ENVIRONMENT/representation error (=> Stage D / cluster)?
#
# Cheap, zero-new-infrastructure discriminator (reuses the validated native oracle):
# freeze the native D-state at a mid-rise time t, then apply ONE square-star update
# at each center c TWO ways from the SAME frozen state and read the exact CTM-free
# density (exact_density_finite, convention B, ledger-independent):
#   n_lossless(c): project_star_pepskit! at maxdim=BIG, cutoff=0, rel_floor=0
#                  -> NO truncation at all (the exact one-gate ceiling).
#   n_D(c)       : project_star_pepskit! at maxdim=D with the trajectory's own
#                  cutoff/rel_floor -> exactly what the simple update does.
# The per-star truncation cost in density is |n_D(c) - n_lossless(c)|.
#
# DECISION RULE: if max_c |n_D - n_lossless| (the per-step truncation density cost)
# is orders of magnitude below the trajectory lag (~1e-2 on the rise), then per-step
# truncation injects negligible density error and a better per-bond metric has almost
# nothing to redistribute -> bare Stage C cannot close the lag on a per-step basis.
#
# HONEST CAVEAT (do not overclaim): this bounds the PER-STEP truncation headroom,
# not the ACCUMULATED headroom over ~10^3 star updates. A fully decisive test needs
# the exact-environment-optimal truncation run over the trajectory (more machinery,
# its own PR). Reported alongside is the truncerr-vs-lag decoupling across D from the
# trajectory probes, which is independent supporting evidence.

using SquarePXPDynamics
using Printf

const TS   = [1.6, 1.8, 2.0, 2.2]   # mid-rise sample times
const DT   = 0.02
const BIG  = 64                      # lossless maxdim (single star grows 4 bonds by <=2x)
const DS   = [2, 4]

# ED reference density on the 0.2 grid (4x4 Neel), for context only.
const EDGRID = Dict(
    1.6 => 0.10056019790056131, 1.8 => 0.1555774954125962,
    2.0 => 0.24370069054283897, 2.2 => 0.3373453880372616,
)

function run_D(D::Int)
    cell = PeriodicSquareUnitCell(4, 4)
    @printf("\n========== D = %d ==========\n", D)
    @printf("%-5s %-10s %-10s %-12s %-12s %-12s\n",
            "t", "n_frozen", "n_ED", "max|dn_trunc|", "rms|dn_trunc|", "max|dn_gate|")
    for t in TS
        psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = D)
        # freeze: evolve to t with the real trajectory settings
        evolve_pepskit!(psi, t; dt = DT, order = 2, schedule = :serial,
                        maxdim = D, cutoff = 1e-12, rel_floor = 1e-3)
        n_frozen = exact_density_finite(psi; max_sites = 16)

        dn_trunc = Float64[]   # |n_D(c) - n_lossless(c)|  : per-star truncation cost
        dn_gate  = Float64[]   # |n_lossless(c) - n_frozen|: per-star gate effect (context)
        for c in cell.reps
            lossless = deepcopy(psi)
            project_star_pepskit!(lossless, c; dt = DT, maxdim = BIG,
                                  cutoff = 0.0, rel_floor = 0.0,
                                  evolution = :real, projected = true)
            n_lossless = exact_density_finite(lossless; max_sites = 16)

            trunc = deepcopy(psi)
            project_star_pepskit!(trunc, c; dt = DT, maxdim = D,
                                  cutoff = 1e-12, rel_floor = 1e-3,
                                  evolution = :real, projected = true)
            n_D = exact_density_finite(trunc; max_sites = 16)

            push!(dn_trunc, abs(n_D - n_lossless))
            push!(dn_gate, abs(n_lossless - n_frozen))
        end
        @printf("%-5.1f %-10.6f %-10.6f %-12.4e %-12.4e %-12.4e\n",
                t, n_frozen, get(EDGRID, t, NaN),
                maximum(dn_trunc), sqrt(sum(abs2, dn_trunc) / length(dn_trunc)),
                maximum(dn_gate))
        flush(stdout)
    end
end

for D in DS
    run_D(D)
end
println("\nDONE")
