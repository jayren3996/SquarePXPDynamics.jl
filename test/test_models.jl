using LinearAlgebra

@testset "models" begin
    Hpxp = pxp_star_hamiltonian(projector_down(), pauli_x())
    @test size(Hpxp) == (128, 128)
    @test Hpxp ≈ Hpxp'

    P = blockade_projector()
    @test size(P) == (128, 128)
    @test P * P ≈ P
    @test P ≈ P'
    all_up = zeros(ComplexF64, 128)
    all_up[1] = 1
    @test norm(P * all_up) == 0
    all_down = zeros(ComplexF64, 128)
    all_down[end] = 1
    @test P * all_down ≈ all_down

    Hcluster = cluster_star_hamiltonian()
    @test size(Hcluster) == (128, 128)
    @test Hcluster ≈ Hcluster'
    @test Hcluster * Hcluster ≈ Matrix{ComplexF64}(I, 128, 128)

    Hdiag = diagonal_star_hamiltonian()
    @test Hdiag ≈ Hdiag'
    @test Hdiag * Hdiag ≈ Matrix{ComplexF64}(I, 128, 128)

    Hising = ising_bond_hamiltonian()
    @test size(Hising) == (4, 4)
    @test Hising ≈ Hising'
end
