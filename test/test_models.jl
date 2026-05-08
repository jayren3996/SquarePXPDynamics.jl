using LinearAlgebra

basis_index(bits) = 1 + sum(bits[i] << (7 - i) for i in 1:7)

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

    center_up_neighbor_up = zeros(ComplexF64, 128)
    center_up_neighbor_up[basis_index([0, 0, 1, 1, 1, 1, 1])] = 1
    @test norm(P * center_up_neighbor_up) == 0

    center_up_neighbors_down = zeros(ComplexF64, 128)
    center_up_neighbors_down[basis_index([0, 1, 1, 1, 1, 1, 1])] = 1
    @test P * center_up_neighbors_down ≈ center_up_neighbors_down

    center_down_neighbor_up = zeros(ComplexF64, 128)
    center_down_neighbor_up[basis_index([1, 0, 1, 1, 1, 1, 1])] = 1
    @test P * center_down_neighbor_up ≈ center_down_neighbor_up

    center_down_multiple_neighbor_ups = zeros(ComplexF64, 128)
    center_down_multiple_neighbor_ups[basis_index([1, 0, 0, 1, 1, 1, 1])] = 1
    @test P * center_down_multiple_neighbor_ups ≈ center_down_multiple_neighbor_ups

    @test_throws ArgumentError pxp_star_hamiltonian(ones(ComplexF64, 3, 3), pauli_x())
    @test_throws ArgumentError pxp_star_hamiltonian(projector_down(), ones(ComplexF64, 3, 3))

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
