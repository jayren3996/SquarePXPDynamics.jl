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
    @test P * Uproj ≈ Uproj

    bad = zeros(ComplexF64, 128)
    bad[1] = 1
    cleaned = Uproj * bad
    @test norm((I - P) * cleaned) < 1e-12
end
