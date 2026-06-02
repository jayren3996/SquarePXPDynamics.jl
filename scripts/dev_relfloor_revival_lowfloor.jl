using SquarePXPDynamics
using Printf

# Decisive: is the revival D-ladder D-MONOTONE at a low/zero floor (the value an
# entanglement-adaptive ramp would use at high entanglement)? If D=4 <= D=3 <= D=2
# at rel_floor=0, an adaptive ramp can achieve the hard rule; if D=4 is still worst,
# the revival non-monotonicity is the mean-field-environment limit (needs full update).
# 4x4 Neel, evolve to revival t=2.6, CTM density vs ED peak.

const ED_REVIVAL = 0.4825
const TARGET_T = 2.6
const DT = 0.02

function evolve_and_measure(rel_floor, D; chi = 8)
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_square_ipeps(cell; excited_on = :even, maxdim = D)
    params = TrotterParams(DT, 2, :real, D, 1e-12; schedule = :serial, rel_floor = rel_floor)
    nsteps = round(Int, TARGET_T / DT)
    local maxtrunc = 0.0
    for _ = 1:nsteps
        log = evolve!(psi, DT; params = params)
        maxtrunc = max(maxtrunc, log.max_truncerr)
    end
    ctx = pepskit_ctmrg_context(psi; params = default_ctmrg_params(chi = chi, maxiter = 100, tol = 1e-8))
    reps = unitcell_reps(psi)
    dens = sum(local_density_ctm(psi, r, ctx) for r in reps) / length(reps)
    return dens, maxtrunc
end

@printf("%-10s %-4s %-12s %-12s %-12s\n", "rel_floor", "D", "CTM_density", "|err vs ED|", "max_truncerr")
for rf in (1e-6, 0.0), D in (2, 3, 4)
    try
        d, mt = evolve_and_measure(rf, D)
        @printf("%-10.0e %-4d %-12.6f %-12.4e %-12.2e\n", rf, D, d, abs(d - ED_REVIVAL), mt)
    catch e
        @printf("%-10.0e %-4d CRASH: %s\n", rf, D, sprint(showerror, e) |> x -> first(split(x, "\n")))
    end
end
println("If D=4 <= D=3 <= D=2 at rel_floor=0 -> adaptive ramp can reach the hard rule.")
println("If D=4 still worst -> mean-field-environment limit; only full update fixes it.")
