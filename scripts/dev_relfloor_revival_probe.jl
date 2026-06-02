using SquarePXPDynamics
using Printf

# Disambiguate the D=4 revival evolution error: is rel_floor=1e-4 OVER-truncating
# at high entanglement (fix = looser/adaptive floor, EASY), or is the simple-update
# environment itself the limit (fix = full-update/Vidal regauge, HARD)?
# 4x4 Neel, D=4, evolve to the revival t=2.6, measure CTM density vs ED peak 0.483.

const ED_REVIVAL = 0.4825   # 4x4 Neel ED n(t) first revival peak at t=2.6
const TARGET_T = 2.6
const DT = 0.02

function evolve_and_measure(rel_floor; D = 4, chi = 8)
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_square_ipeps(cell; excited_on = :even, maxdim = D)
    params = TrotterParams(DT, 2, :real, D, 1e-12; schedule = :serial, rel_floor = rel_floor)
    nsteps = round(Int, TARGET_T / DT)
    local maxtrunc = 0.0
    for _ in 1:nsteps
        log = evolve!(psi, DT; params = params)
        maxtrunc = max(maxtrunc, log.max_truncerr)
    end
    ctx = pepskit_ctmrg_context(psi; params = default_ctmrg_params(chi = chi, maxiter = 100, tol = 1e-8))
    reps = unitcell_reps(psi)
    dens = sum(local_density_ctm(psi, r, ctx) for r in reps) / length(reps)
    return dens, maxtrunc
end

println("######## D=4 rel_floor sweep at the revival (t=2.6, ED peak=$(ED_REVIVAL)) ########")
@printf("%-12s %-12s %-12s %-12s\n", "rel_floor", "CTM_density", "|err vs ED|", "max_truncerr")
for rf in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0.0)
    try
        d, mt = evolve_and_measure(rf)
        @printf("%-12.0e %-12.6f %-12.4e %-12.2e\n", rf, d, abs(d - ED_REVIVAL), mt)
    catch e
        @printf("%-12.0e CRASH: %s\n", rf, sprint(showerror, e) |> x -> first(split(x, "\n")))
    end
end
println("######## DONE ########")
println("If a looser rel_floor (1e-5/1e-6) gives lower revival error -> EASY (tune floor).")
println("If error is flat/worse across rel_floor -> HARD (simple-update env limit; need full-update).")
