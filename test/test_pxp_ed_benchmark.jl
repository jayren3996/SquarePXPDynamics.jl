using Test
using LinearAlgebra
using SparseArrays
using SquarePXPDynamics

@testset "PXP ED space-group basis counts match periodic hard-square counts" begin
    b3 = pxp_ed_space_group_basis(3)
    b4 = pxp_ed_space_group_basis(4)

    @test size(b3, 1) == 4
    @test size(b4, 1) == 29
    @test pxp_ed_constrained_count(b3) == 34
    @test pxp_ed_constrained_count(b4) == 743
    @test pxp_ed_group_order(b3) == 8 * 3^2
    @test pxp_ed_group_order(b4) == 8 * 4^2
end

@testset "PXP ED basis indexes constrained configurations by orbit" begin
    basis = pxp_ed_space_group_basis(4)
    down = fill(1, 16)
    translated_single = fill(1, 16)
    translated_single[6] = 0

    psi_down = pxp_ed_initial_state(basis; state = :down)
    c_down, pos_down = SquarePXPDynamics.EDKit.index(basis, down)
    c_single, pos_single = SquarePXPDynamics.EDKit.index(basis, translated_single)

    @test norm(psi_down) ≈ 1.0 atol = 1e-12
    @test psi_down[pos_down] ≈ 1.0 atol = 1e-12
    @test c_down > c_single
    @test pos_down != pos_single

    forbidden = fill(1, 16)
    forbidden[1] = 0
    forbidden[2] = 0
    c_forbidden, _ = SquarePXPDynamics.EDKit.index(basis, forbidden)
    @test iszero(c_forbidden)
end

@testset "PXP ED Hamiltonian and Krylov benchmark run through EDKit" begin
    basis = pxp_ed_space_group_basis(3)
    H = pxp_ed_hamiltonian_operator(basis)
    Hs = sparse_pxp_ed_hamiltonian(basis)

    @test size(H) == (4, 4)
    @test size(Hs) == (4, 4)
    @test Matrix(Hs) ≈ Matrix(Hs)' atol = 1e-12

    config = PXPEEDBenchmarkConfig(3; total_time = 0.02, dt = 0.01, measure_every = 1)
    result = run_pxp_ed_benchmark(config)

    @test result.basis_dimension == 4
    @test result.hamiltonian_nnz == nnz(Hs)
    @test [s.step for s in result.samples] == [0, 1, 2]
    @test result.samples[end].time ≈ 0.02 atol = 1e-12
    @test result.samples[end].norm ≈ 1.0 atol = 1e-10
    @test 0.0 <= result.samples[end].return_probability <= 1.0
    @test 0.0 <= result.samples[end].excitation_density <= 1.0
    @test result.diagnostics.matvecs > 0
end
