# PEPSKit-native backend tests (Milestones 1-5).
#
# Kept deliberately CTM-cheap: only 1-site (density) and 2-site (blockade) CTM
# expectations are evaluated here. The 5-site star/energy CTM expectation is
# correct but compile-heavy (PEPSKit's `@tensor`/optimaltree for the cross
# network), so its parity is exercised in scripts/dev/verify_m1_m3.jl and
# scripts/dev/verify_m4.jl rather than in the routine suite.
using Test
using SquarePXPDynamics
using SquarePXPDynamics: PeriodicSquareUnitCell, SquareCoord
using SquarePXPDynamics.SquarePXP: square_pxp_star_hamiltonian, square_pxp_gate
using TensorKit
using PEPSKit
using LinearAlgebra: I

const SPD = SquarePXPDynamics

# Dense matrix from a 5-site gate TensorMap, in (center,right,up,left,down) order.
function _tm_to_dense(tm, nsites)
    arr = convert(Array, tm)
    d = 2^nsites
    M = zeros(ComplexF64, d, d)
    rng = ntuple(_ -> 1:2, nsites)
    didx(v) = 1 + sum((v[i] - 1) * 2^(nsites - i) for i in 1:nsites)
    for ov in Iterators.product(rng...), iv in Iterators.product(rng...)
        M[didx(ov), didx(iv)] = arr[ov..., iv...]
    end
    return M
end

ctm = SPD.default_ctmrg_params(; chi = 4, tol = 1e-10, maxiter = 60, seed = 0)

@testset "PEPSKit-native backend" begin
    @testset "M1 constructors" begin
        cell = PeriodicSquareUnitCell(2, 2)
        down = product_pxp_pepskit_state(cell; state = :down, D = 1)
        @test size(pepskit_peps(down)) == (2, 2)
        @test size(pepskit_weights(down)) == (2, 2, 2)
        @test state_version(down) == 0
        mark_mutated!(down)
        @test state_version(down) == 1
        add_log_norm!(down, 0.25)
        @test log_norm(down) == 0.25
        d2 = copy_state(down)
        mark_mutated!(d2)
        @test state_version(down) == 1 && state_version(d2) == 2
        D3 = product_pxp_pepskit_state(cell; state = :down, D = 3)
        @test dim(space(pepskit_weights(D3)[1, 1, 1], 1)) == 3
    end

    @testset "M2 CTM observables (1- and 2-site)" begin
        cell = PeriodicSquareUnitCell(4, 4)
        down = product_pxp_pepskit_state(cell; state = :down, D = 1)
        up = product_pxp_pepskit_state(cell; state = :up, D = 1)
        cxd = pepskit_native_context(down; params = ctm)
        cxu = pepskit_native_context(up; params = ctm)
        @test abs(density_ctm_pepskit(cxd, SquareCoord(1, 1))) < 1e-8
        @test abs(density_ctm_pepskit(cxu, SquareCoord(1, 1)) - 1) < 1e-8
        @test abs(blockade_violation_ctm_pepskit(cxd)) < 1e-8
        @test abs(blockade_violation_ctm_pepskit(cxu) - 32) < 1e-6
        cb = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = 1)
        cx = pepskit_native_context(cb; params = ctm)
        @test abs(density_ctm_pepskit(cx, SquareCoord(2, 2)) - 1) < 1e-8
        @test abs(density_ctm_pepskit(cx, SquareCoord(1, 2))) < 1e-8
        # parity vs legacy adapter on 1-site density (cheap, unambiguous)
        psi_old = checkerboard_square_ipeps(cell; excited_on = :even, maxdim = 1)
        ctx_old = pepskit_ctmrg_context(psi_old; params = ctm)
        @test abs(density_ctm_pepskit(cx, SquareCoord(2, 2)) -
                  local_density_ctm(psi_old, SquareCoord(2, 2), ctx_old)) < 1e-6
        @test abs(density_ctm_pepskit(cx, SquareCoord(1, 2)) -
                  local_density_ctm(psi_old, SquareCoord(1, 2), ctx_old)) < 1e-6
        # operator constructors build without error (5-site evaluated elsewhere)
        @test pxp_star_operator(cell, SquareCoord(2, 2)) isa PEPSKit.LocalOperator
        @test pxp_energy_operator(cell) isa PEPSKit.LocalOperator
    end

    @testset "M3 star gate" begin
        @test pxp_star_hamiltonian_matrix() == square_pxp_star_hamiltonian()
        g0 = pxp_star_gate_tensormap(0.0; projected = false, evolution = :real)
        @test isapprox(g0, TensorKit.id(domain(g0)); atol = 1e-12)
        Gr = _tm_to_dense(pxp_star_gate_tensormap(0.123; projected = false), 5)
        @test isapprox(Gr, square_pxp_gate(0.123; evolution = :real); atol = 1e-10)
        U = _tm_to_dense(pxp_star_gate_tensormap(0.37; projected = false), 5)
        @test isapprox(U' * U, Matrix{ComplexF64}(I, 32, 32); atol = 1e-10)
        Gi = _tm_to_dense(pxp_star_gate_tensormap(0.2; projected = false, evolution = :imaginary), 5)
        @test all(isfinite, Gi) && isapprox(Gi, Gi'; atol = 1e-10)
    end

    @testset "M4 star simple update" begin
        cell = PeriodicSquareUnitCell(4, 4)
        center = SquareCoord(2, 2)
        # dt=0 identity preserves center density
        st = product_pxp_pepskit_state(cell; state = :down, D = 1)
        d0 = density_ctm_pepskit(pepskit_native_context(st; params = ctm), center)
        info = project_star_pepskit!(st, center; dt = 0.0, maxdim = 4, projected = false)
        @test abs(density_ctm_pepskit(pepskit_native_context(st; params = ctm), center) - d0) < 1e-8
        @test length(info.truncerrs) == 4
        # maxdim cap
        st3 = product_pxp_pepskit_state(cell; state = :down, D = 1)
        i3 = project_star_pepskit!(st3, center; dt = 0.2, maxdim = 2, projected = true)
        @test all(<=(2), values(i3.keptdims))
        @test maximum(dim(space(w, 1)) for w in st3.weights.data) <= 2
        # parity vs legacy ITensors backend (1 step, D=1) at the center site
        dt = 0.1
        stp = product_pxp_pepskit_state(cell; state = :down, D = 1)
        project_star_pepskit!(stp, center; dt, maxdim = 4, projected = true, rel_floor = 0.0)
        dn = density_ctm_pepskit(pepskit_native_context(stp; params = ctm), center)
        psi = product_square_ipeps(cell; state = :down, maxdim = 1)
        project_star!(psi, center, dt; projected = true, maxdim = 4, rel_floor = 0.0)
        dl = local_density_ctm(psi, center, pepskit_ctmrg_context(psi; params = ctm))
        @test abs(dn - dl) < 1e-6
        @test abs(dn - sin(dt)^2) < 1e-6           # center rotates to sin^2(dt) excited
        # reversibility
        sr = product_pxp_pepskit_state(cell; state = :down, D = 1)
        project_star_pepskit!(sr, center; dt = 0.1, maxdim = 4, projected = true, rel_floor = 0.0)
        project_star_pepskit!(sr, center; dt = -0.1, maxdim = 4, projected = true, rel_floor = 0.0)
        @test abs(density_ctm_pepskit(pepskit_native_context(sr; params = ctm), center)) < 1e-8
    end

    @testset "M5 evolution schedule" begin
        p1 = TrotterParamsPK(; dt = 0.05, order = 1, maxdim = 2)
        @test trotter_sequence_pk(p1) == [(c, 0.05) for c in 1:5]
        p2 = TrotterParamsPK(; dt = 0.1, order = 2, maxdim = 2)
        seq = trotter_sequence_pk(p2)
        @test first.(seq) == [1, 2, 3, 4, 5, 4, 3, 2, 1]   # palindromic
        @test seq[5] == (5, 0.1)
        @test seq[1] == (1, 0.05)
        # five-colour evolution on a 5-colour-compatible cell, density stays finite
        cell = PeriodicSquareUnitCell(5, 5)
        st = product_pxp_pepskit_state(cell; state = :down, D = 1)
        log = evolve_pepskit!(st, 0.1; dt = 0.05, order = 2, maxdim = 4)
        @test log.nsteps == 2 && isfinite(log.max_truncerr)
        cx = pepskit_native_context(st; params = ctm)
        @test isfinite(density_ctm_pepskit(cx, SquareCoord(2, 2)))
        # serial schedule on a generic (non-5-colour) cell
        cell4 = PeriodicSquareUnitCell(4, 4)
        st4 = product_pxp_pepskit_state(cell4; state = :down, D = 1)
        log4 = evolve_pepskit!(st4, 0.1; dt = 0.05, order = 2, schedule = :serial, maxdim = 4)
        @test log4.nsteps == 2 && isfinite(log4.max_truncerr)
        @test isfinite(density_ctm_pepskit(pepskit_native_context(st4; params = ctm), SquareCoord(1, 1)))
    end
end
