# Cross-check the double-layer oracle against the trusted dense path.
# M1: naive double-layer contraction must equal `exact_density_finite` (dense)
#     on 3x3 to ~1e-10, and reproduce product-state limits.
using SquarePXPDynamics
using SquarePXPDynamics.FiniteIPEPSObservables: exact_density_finite
const FIO = SquarePXPDynamics.FiniteIPEPSObservables
using Printf

# Build a 3x3 D-bond entangled state by serial star projections (mirrors the test
# helper _finite_obs_state_after_serial_stars).
function serial_star_state(Lx, Ly, nstars; D = 2, step = 0.02, cutoff = 1e-12)
    cell = PeriodicSquareUnitCell(Lx, Ly)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    for center in cell.reps[1:nstars]
        project_star!(psi, center, step; evolution = :real, projected = true,
                      maxdim = D, cutoff = cutoff)
    end
    return psi
end

ok = true
check(name, a, b; atol = 1e-10) = begin
    d = abs(a - b)
    pass = d <= atol
    global ok &= pass
    @printf("%-44s dense=% .12f  naive=% .12f  |Δ|=%.2e  %s\n",
            name, a, b, d, pass ? "OK" : "FAIL")
end

println("== product-state limits (3x3) ==")
let cell = PeriodicSquareUnitCell(3, 3)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    up = product_square_ipeps(cell; state = :up, maxdim = 1)
    check("down -> 0", 0.0, FIO._naive_double_layer_density(down); atol = 1e-14)
    check("up   -> 1", 1.0, FIO._naive_double_layer_density(up); atol = 1e-14)
end

println("== naive double-layer vs dense (3x3, serial stars) ==")
for D in (2, 3, 4)
    psi = serial_star_state(3, 3, 3; D = D)
    dense = exact_density_finite(psi; max_sites = 16)
    naive = FIO._naive_double_layer_density(psi; max_sites = 16)
    check("3x3 D=$D", dense, naive)
end

println("== boundary sweep vs naive (3x3 product limits) ==")
let cell = PeriodicSquareUnitCell(3, 3)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    up = product_square_ipeps(cell; state = :up, maxdim = 1)
    check("bnd down -> 0", 0.0, FIO._boundary_density_finite(down); atol = 1e-12)
    check("bnd up   -> 1", 1.0, FIO._boundary_density_finite(up); atol = 1e-12)
end

println("== boundary sweep vs dense (3x3, serial stars) ==")
for D in (2, 3, 4)
    psi = serial_star_state(3, 3, 3; D = D)
    dense = exact_density_finite(psi; max_sites = 16)
    bnd = FIO._boundary_density_finite(psi; max_sites = 16)
    check("3x3 bnd D=$D", dense, bnd)
end

println("== boundary sweep vs dense (4x4, serial stars) ==")
for D in (2, 3, 4)
    psi = serial_star_state(4, 4, 4; D = D)
    t0 = time()
    dense = exact_density_finite(psi; max_sites = 16)
    td = time() - t0
    t0 = time()
    bnd = FIO._boundary_density_finite(psi; max_sites = 16)
    tb = time() - t0
    check("4x4 bnd D=$D", dense, bnd)
    @printf("    (dense %.1fs, boundary %.1fs)\n", td, tb)
end

println(ok ? "\nALL CHECKS PASSED" : "\nCHECKS FAILED")
exit(ok ? 0 : 1)
