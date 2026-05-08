@testset "solvable models" begin
    H = cluster_star_hamiltonian()
    for t in (0.0, 0.1, 0.7)
        @test stabilizer_expectation_exact(t; initial = :plus) ≈ cos(2t)
    end

    U1 = dense_gate(H, 0.11; evolution = :real)
    U2 = dense_gate(H, 0.23; evolution = :real)
    U12 = dense_gate(H, 0.34; evolution = :real)
    @test U2 * U1 ≈ U12
end
