using SquarePXPDynamics
using Printf

# Stage-2 oracle decontamination. The 4x4 Neel-to-revival benchmark has measured
# density with CTM chi=8, so the D-ladder error mixes evolution error with a CTM
# environment approximation. exact_density_finite(max_sites=16) contracts the
# 16-site cell EXACTLY (no environment). This script:
#   1. verifies the refactored exact_density_finite (build dense state once) equals
#      the old per-site formula on 3x3 to machine precision (correctness), and
#   2. re-runs the revival for (rel_floor in {1e-4,1e-3}) x (D in {2,3,4}) measuring
#      BOTH CTM chi=8 (old) and exact-16 (new), vs ED peak 0.4825 -- so |CTM-EXACT|
#      is the prior CTM contamination and |EXACT-ED| is the true evolution error.

const ED = 0.4825
const TARGET_T = 2.6
const DT = 0.02

function verify_refactor()
    # 3x3 is odd x odd, so checkerboard hits the blockade seam; use :down (the
    # safe quench state) to entangle a tiny cell for the refactor equality check.
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    params = TrotterParams(0.04, 1, :real, 2, 1e-12; schedule = :serial)
    for _ = 1:10
        evolve!(psi, 0.04; params = params)   # entangle (-> t=0.4)
    end
    n = projector_up()
    reps = unitcell_reps(psi)
    old = sum(real(exact_one_site_expectation_finite(psi, c, n; max_sites = 9)) for c in reps) / length(reps)
    new = exact_density_finite(psi; max_sites = 9)
    @printf("refactor check (3x3 entangled): old=%.12f  new=%.12f  |diff|=%.2e\n", old, new, abs(old - new))
    return abs(old - new)
end

function decontaminate(rel_floor, D; chi = 8)
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_square_ipeps(cell; excited_on = :even, maxdim = D)
    params = TrotterParams(DT, 2, :real, D, 1e-12; schedule = :serial, rel_floor = rel_floor)
    for _ = 1:round(Int, TARGET_T / DT)
        evolve!(psi, DT; params = params)
    end
    ctx = pepskit_ctmrg_context(psi; params = default_ctmrg_params(chi = chi, maxiter = 100, tol = 1e-8))
    reps = unitcell_reps(psi)
    ctm = sum(local_density_ctm(psi, r, ctx) for r in reps) / length(reps)
    exact = exact_density_finite(psi; max_sites = 16)
    return ctm, exact
end

verify_refactor()
println()
@printf("%-9s %-3s %-11s %-11s %-12s %-12s %-10s\n",
    "rel_floor", "D", "CTM_chi8", "EXACT_16", "|CTM-ED|", "|EXACT-ED|", "contam")
for rf in (1e-4, 1e-3), D in (2, 3, 4)
    ctm, exact = decontaminate(rf, D)
    @printf("%-9.0e %-3d %-11.6f %-11.6f %-12.4e %-12.4e %-10.2e\n",
        rf, D, ctm, exact, abs(ctm - ED), abs(exact - ED), abs(ctm - exact))
    flush(stdout)
end
println("\nIf |EXACT-ED| is D-monotone where |CTM-ED| was not -> CTM was contaminating.")
println("If EXACT still non-monotone at the peak -> evolution (env-ceiling) error, confirmed exact.")
println("DONE")
