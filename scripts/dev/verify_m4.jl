# Verify Milestone 4: PEPSKit-native star simple update.
# Run:  julia --project=. scripts/dev/verify_m4.jl
using SquarePXPDynamics
using SquarePXPDynamics: PeriodicSquareUnitCell, SquareCoord
using TensorKit
using PEPSKit

const SPD = SquarePXPDynamics
npass = 0; nfail = 0
function check(name, cond)
    global npass, nfail
    cond ? (npass += 1; println("  PASS: ", name)) : (nfail += 1; println("  FAIL: ", name))
end

cell = PeriodicSquareUnitCell(4, 4)
params = SPD.default_ctmrg_params(; chi = 4, tol = 1e-10, maxiter = 80, seed = 0)
center = SquareCoord(2, 2)

function tensors_finite(state)
    for t in state.peps.A
        all(isfinite, convert(Array, t)) || return false
    end
    for w in state.weights.data
        all(isfinite, w.data) || return false
    end
    return true
end

println("== M4: dt=0 identity update leaves density unchanged ==")
st = product_pxp_pepskit_state(cell; state = :down, D = 1)
d_before = measure_ctm_pepskit(st; params).mean_density
info = project_star_pepskit!(st, center; dt = 0.0, maxdim = 4, projected = false)
d_after = measure_ctm_pepskit(st; params).mean_density
println("    density before=$d_before  after=$d_after")
check("dt=0 (identity) preserves density", abs(d_after - d_before) < 1e-8)
check("dt=0 update tensors finite", tensors_finite(st))
check("dt=0 info has 4 directions", length(info.truncerrs) == 4)

println("== M4: dt=0 projected update on blockade-allowed state ==")
st2 = product_pxp_pepskit_state(cell; state = :down, D = 1)
project_star_pepskit!(st2, center; dt = 0.0, maxdim = 4, projected = true)
d2 = measure_ctm_pepskit(st2; params).mean_density
check("dt=0 projected preserves all-down density ~0", abs(d2) < 1e-8)

println("== M4: maxdim cap respected ==")
st3 = product_pxp_pepskit_state(cell; state = :down, D = 1)
info3 = project_star_pepskit!(st3, center; dt = 0.2, maxdim = 2, projected = true)
maxbond = maximum(dim(space(w, 1)) for w in st3.weights.data)
println("    kept dims = ", info3.keptdims, "  max bond dim = ", maxbond)
check("all kept dims <= maxdim", all(<=(2), values(info3.keptdims)))
check("bond dims <= maxdim", maxbond <= 2)
check("nonzero dt update finite", tensors_finite(st3))

println("== M4: D=1 short-time parity vs legacy ITensors backend ==")
# Apply one identical star update at the same center, compare mean density.
dt = 0.1
stp = product_pxp_pepskit_state(cell; state = :down, D = 1)
project_star_pepskit!(stp, center; dt = dt, maxdim = 4, projected = true, rel_floor = 0.0)
dens_native = measure_ctm_pepskit(stp; params).mean_density

psi = product_square_ipeps(cell; state = :down, maxdim = 1)
project_star!(psi, center, dt; projected = true, maxdim = 4, rel_floor = 0.0)
dens_legacy = measure_ctm(psi; params).density
# analytic: center rotates to sin^2(dt) excited; mean over 16 sites
dens_analytic = sin(dt)^2 / length(cell.reps)
println("    native=$dens_native  legacy=$dens_legacy  analytic=$dens_analytic")
check("native vs legacy mean density (1 step)", abs(dens_native - dens_legacy) < 1e-6)
check("native vs analytic mean density (1 step)", abs(dens_native - dens_analytic) < 1e-6)

println("== M4: reversibility (forward dt then -dt) ==")
str = product_pxp_pepskit_state(cell; state = :down, D = 1)
project_star_pepskit!(str, center; dt = 0.1, maxdim = 4, projected = true, rel_floor = 0.0)
project_star_pepskit!(str, center; dt = -0.1, maxdim = 4, projected = true, rel_floor = 0.0)
dens_rev = measure_ctm_pepskit(str; params).mean_density
println("    density after forward+reverse = ", dens_rev)
check("reversible to ~0 density", abs(dens_rev) < 1e-8)

println("\n== SUMMARY: $npass passed, $nfail failed ==")
nfail == 0 || error("M4 verification failed")
