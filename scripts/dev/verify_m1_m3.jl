# Verify Milestones 1-3 of the PEPSKit-native backend.
# Run:  julia --project=. scripts/dev/verify_m1_m3.jl
using SquarePXPDynamics
using SquarePXPDynamics: PeriodicSquareUnitCell, SquareCoord
using SquarePXPDynamics.SquarePXP: square_pxp_star_hamiltonian, square_pxp_gate
using TensorKit
using PEPSKit
using LinearAlgebra: I

const SPD = SquarePXPDynamics

npass = 0
nfail = 0
function check(name, cond)
    global npass, nfail
    if cond
        npass += 1
        println("  PASS: ", name)
    else
        nfail += 1
        println("  FAIL: ", name)
    end
end

# Dense matrix from a 5-site gate TensorMap, in (center,right,up,left,down) order.
function tensormap_to_dense(tm, nsites)
    arr = convert(Array, tm)
    dim = 2^nsites
    M = zeros(ComplexF64, dim, dim)
    rng = ntuple(_ -> 1:2, nsites)
    didx(vals) = 1 + sum((vals[i] - 1) * 2^(nsites - i) for i in 1:nsites)
    for outv in Iterators.product(rng...), inv in Iterators.product(rng...)
        M[didx(outv), didx(inv)] = arr[outv..., inv...]
    end
    return M
end

println("== M1: state constructors ==")
cell22 = PeriodicSquareUnitCell(2, 2)
down = product_pxp_pepskit_state(cell22; state = :down, D = 1)
up = product_pxp_pepskit_state(cell22; state = :up, D = 1)
check("InfinitePEPS size == (Ly, Lx)", size(pepskit_peps(down)) == (2, 2))
check("SUWeight size == (2, Ly, Lx)", size(pepskit_weights(down)) == (2, 2, 2))
check("state_version starts 0", state_version(down) == 0)
mark_mutated!(down)
check("mark_mutated! increments", state_version(down) == 1)
add_log_norm!(down, 0.5)
check("add_log_norm! accumulates", log_norm(down) == 0.5)
d2 = copy_state(down)
mark_mutated!(d2)
check("copy_state is independent", state_version(down) == 1 && state_version(d2) == 2)

D3 = product_pxp_pepskit_state(cell22; state = :down, D = 3)
check("D=3 virtual dim on weights", dim(space(pepskit_weights(D3)[1, 1, 1], 1)) == 3)

println("== M1/M2: CTM densities on product states ==")
# 4x4 cell: even -> checkerboard tiles; >=3 -> 5-site stars are distinct.
cell44 = PeriodicSquareUnitCell(4, 4)
params = SPD.default_ctmrg_params(; chi = 4, tol = 1e-10, maxiter = 80, seed = 0)
down4 = product_pxp_pepskit_state(cell44; state = :down, D = 1)
up4 = product_pxp_pepskit_state(cell44; state = :up, D = 1)
sdown = measure_ctm_pepskit(down4; params)
sup = measure_ctm_pepskit(up4; params)
println("    all-down mean density = ", sdown.mean_density)
println("    all-up   mean density = ", sup.mean_density)
check("all-down density ~ 0", abs(sdown.mean_density) < 1e-8)
check("all-up density ~ 1", abs(sup.mean_density - 1) < 1e-8)
check("all-down blockade ~ 0", abs(sdown.blockade_violation) < 1e-8)
check("all-up blockade == nbonds(32)", abs(sup.blockade_violation - 32) < 1e-6)

cb_even = checkerboard_pxp_pepskit_state(cell44; excited_on = :even, D = 1)
ctx = pepskit_native_context(cb_even; params)
d_even = density_ctm_pepskit(ctx, SquareCoord(2, 2))   # x+y=4 even -> excited
d_odd = density_ctm_pepskit(ctx, SquareCoord(1, 2))    # x+y=3 odd  -> vacuum
println("    checkerboard even-site density = ", d_even, ", odd-site = ", d_odd)
check("checkerboard even-site density ~ 1", abs(d_even - 1) < 1e-8)
check("checkerboard odd-site density ~ 0", abs(d_odd) < 1e-8)

println("== M2 parity vs legacy adapter (D=1) ==")
psi_old = checkerboard_square_ipeps(cell44; excited_on = :even, maxdim = 1)
old = measure_ctm(psi_old; params)
println("    legacy mean density   = ", old.density)
cb_summary = measure_ctm_pepskit(cb_even; params)
println("    native mean density   = ", cb_summary.mean_density)
check("native vs legacy mean density", abs(cb_summary.mean_density - old.density) < 1e-6)
check("native vs legacy blockade", abs(cb_summary.blockade_violation - old.blockade_violation) < 1e-6)
check("native vs legacy energy density", abs(cb_summary.pxp_energy_density - old.pxp_energy_density) < 1e-6)

println("== M3: PXP star gate ==")
check("pxp_star_hamiltonian_matrix == legacy", pxp_star_hamiltonian_matrix() == square_pxp_star_hamiltonian())

g0 = pxp_star_gate_tensormap(0.0; projected = false, evolution = :real)
check("dt=0 gate is identity", isapprox(g0, TensorKit.id(domain(g0)); atol = 1e-12))

# roundtrip faithfulness
gr = pxp_star_gate_tensormap(0.123; projected = false, evolution = :real)
Gr = tensormap_to_dense(gr, 5)
Gexpected = square_pxp_gate(0.123; evolution = :real)
check("real gate roundtrip matches dense", isapprox(Gr, Gexpected; atol = 1e-10))

# unitarity of unprojected real-time gate
U = tensormap_to_dense(pxp_star_gate_tensormap(0.37; projected = false, evolution = :real), 5)
check("unprojected real gate is unitary", isapprox(U' * U, Matrix{ComplexF64}(I, 32, 32); atol = 1e-10))

# imaginary-time gate finite + Hermitian-positive
Gi = tensormap_to_dense(pxp_star_gate_tensormap(0.2; projected = false, evolution = :imaginary), 5)
check("imag gate finite", all(isfinite, Gi))
check("imag gate Hermitian", isapprox(Gi, Gi'; atol = 1e-10))

println("\n== SUMMARY: $npass passed, $nfail failed ==")
nfail == 0 || error("M1-M3 verification failed")
