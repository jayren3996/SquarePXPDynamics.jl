using LinearAlgebra

@testset "solvable models" begin
    H = cluster_star_hamiltonian()
    for t in (0.0, 0.1, 0.7)
        @test stabilizer_expectation_exact(t; initial = :z_plus) ≈ cos(2t)
    end
    @test_throws ArgumentError stabilizer_expectation_exact(0.1; initial = :minus)

    U1 = dense_gate(H, 0.11; evolution = :real)
    U2 = dense_gate(H, 0.23; evolution = :real)
    U12 = dense_gate(H, 0.34; evolution = :real)
    @test U2 * U1 ≈ U12

    t = 0.19
    U = dense_gate(H, t; evolution = :real)
    @test U ≈ cos(t) * Matrix{ComplexF64}(I, 128, 128) - 1im * sin(t) * H
end
