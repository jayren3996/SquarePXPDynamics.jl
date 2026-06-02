using SquarePXPDynamics
using Printf

# Confirm on the CURRENT code whether rel_floor=1e-3 fixes the revival bad pocket
# that the default 1e-4 has, across the D-ladder. 4x4 Neel, evolve to the first
# n(t) revival t=2.6, CTM density vs ED peak. Short-time sweep already showed
# 1e-3 is safe (D=4 err 8.09e-5) and 1e-4 best (2.5e-5) at t=0.2.

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
for rf in (1e-4, 1e-3), D in (2, 3, 4)
    try
        d, mt = evolve_and_measure(rf, D)
        @printf("%-10.0e %-4d %-12.6f %-12.4e %-12.2e\n", rf, D, d, abs(d - ED_REVIVAL), mt)
    catch e
        @printf("%-10.0e %-4d CRASH: %s\n", rf, D, sprint(showerror, e) |> x -> first(split(x, "\n")))
    end
end
println("Decision: rel_floor=1e-3 is better default iff it is D-monotone AND its")
println("D=4 revival error << the 1e-4 D=4 error (the bad pocket ~0.0215).")
