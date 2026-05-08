using LinearAlgebra

@testset "gates" begin
    H = cluster_star_hamiltonian()
    U = dense_gate(H, 0.2; evolution = :real)
    @test U' * U ≈ Matrix{ComplexF64}(I, 128, 128)

    G = dense_gate(H, 0.2; evolution = :imaginary)
    @test G ≈ G'

    Hpxp = pxp_star_hamiltonian(projector_down(), pauli_x())
    Uproj = projected_gate(Hpxp, 0.05; evolution = :real)
    P = blockade_projector()
    @test size(Uproj) == (128, 128)
    @test Uproj ≈ P * dense_gate(Hpxp, 0.05; evolution = :real)
    @test P * Uproj ≈ Uproj

    bad = zeros(ComplexF64, 128)
    bad[1] = 1
    cleaned = Uproj * bad
    @test norm((I - P) * cleaned) < 1e-12

    @test G ≈ cosh(0.2) * Matrix{ComplexF64}(I, 128, 128) - sinh(0.2) * H
    @test_throws ArgumentError dense_gate(zeros(ComplexF64, 4, 4), 0.1; evolution = :real)
    @test_throws ArgumentError dense_gate(zeros(ComplexF64, 127, 128), 0.1; evolution = :real)
    @test_throws ArgumentError dense_gate(H, 0.1; evolution = :bad)
    @test_throws ArgumentError projected_gate(zeros(ComplexF64, 4, 4), 0.1; evolution = :real)
    @test_throws ArgumentError projected_gate(Hpxp, 0.1; evolution = :real, projector = zeros(ComplexF64, 4, 4))
end
