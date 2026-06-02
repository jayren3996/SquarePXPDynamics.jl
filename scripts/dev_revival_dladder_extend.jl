using SquarePXPDynamics
using Printf

# Characterize the mean-field-environment ceiling with the EXACT 16-site oracle:
# does the 4x4 Neel revival error keep improving as D grows past 4, or plateau/
# worsen? If it plateaus/worsens, simple update is fundamentally capped (env
# ceiling) and an environment-aware (cluster/full) update is justified. If it
# keeps decreasing, we just need more D. Production setting rel_floor=1e-3.
# Measured by exact_density_finite(16) (no CTM), vs ED revival peak 0.4825.

const ED = 0.4825
const TARGET_T = 2.6
const DT = 0.02

function revival_error(D)
    psi = checkerboard_square_ipeps(PeriodicSquareUnitCell(4, 4); excited_on = :even, maxdim = D)
    params = TrotterParams(DT, 2, :real, D, 1e-12; schedule = :serial)  # default rel_floor 1e-3
    local maxtrunc = 0.0
    t0 = time()
    for _ = 1:round(Int, TARGET_T / DT)
        log = evolve!(psi, DT; params = params)
        maxtrunc = max(maxtrunc, log.max_truncerr)
    end
    dens = exact_density_finite(psi; max_sites = 16)
    return dens, abs(dens - ED), maxtrunc, time() - t0
end

@printf("%-4s %-12s %-12s %-12s %-8s\n", "D", "exact_dens", "|err vs ED|", "max_truncerr", "secs")
for D in (2, 3, 4, 5, 6)
    try
        d, err, mt, secs = revival_error(D)
        @printf("%-4d %-12.6f %-12.4e %-12.2e %-8.0f\n", D, d, err, mt, secs)
    catch e
        @printf("%-4d CRASH: %s\n", D, first(split(sprint(showerror, e), "\n")))
    end
    flush(stdout)
end
println("Plateau/worsen past D=4 -> mean-field-env ceiling (env-aware update justified).")
println("Keeps decreasing -> just need more D.")
println("DONE")
