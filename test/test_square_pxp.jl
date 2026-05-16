using LinearAlgebra

@testset "square PXP" begin
    H = square_pxp_star_hamiltonian()
    P = square_star_blockade_projector()

    @test size(H) == (32, 32)
    @test size(P) == (32, 32)
    @test P * P ≈ P
    @test P' ≈ P

    @test !square_star_basis_allowed((0, 0, 1, 1, 1))
    @test !square_star_basis_allowed((0, 1, 0, 1, 1))
    @test square_star_basis_allowed((0, 1, 1, 1, 1))
    @test square_star_basis_allowed((1, 0, 0, 0, 0))
    @test_throws ArgumentError square_star_basis_allowed((0, 1))
    @test_throws ArgumentError square_star_basis_allowed((0, 2, 1, 1, 1))
    @test_throws ArgumentError square_star_basis_allowed((0, -1, 1, 1, 1))
    @test_throws ArgumentError square_star_basis_allowed((0, 1, 1, 1, true))

    U1 = square_pxp_gate(0.11; evolution = :real)
    U2 = square_pxp_gate(0.23; evolution = :real)
    U12 = square_pxp_gate(0.34; evolution = :real)
    @test U2 * U1 ≈ U12

    t = 0.19
    U = square_pxp_gate(t; evolution = :real)
    @test U ≈ exp(-1im * t * H)
    @test projected_square_pxp_gate(t; evolution = :real) ≈ P * U
    @test projected_square_pxp_gate(t; evolution = :imaginary) ≈ P * exp(-t * H)
    @test_throws ArgumentError square_pxp_gate(t; evolution = :thermal)
    @test_throws ArgumentError square_pxp_gate(Inf; evolution = :real)
    @test_throws ArgumentError square_pxp_gate(NaN; evolution = :real)
    @test_throws ArgumentError projected_square_pxp_gate(Inf; evolution = :real)
    @test_throws ArgumentError projected_square_pxp_gate(NaN; evolution = :real)
    @test size(square_pxp_gate(-t; evolution = :real)) == (32, 32)
    @test size(projected_square_pxp_gate(-t; evolution = :real)) == (32, 32)
end

@testset "projected PXP gate documents left and sandwich projection" begin
    dt = 0.01
    P = square_star_blockade_projector()
    U = square_pxp_gate(dt; evolution = :real)

    @test projected_square_pxp_gate(dt; evolution = :real) ≈ P * U
    @test projected_square_pxp_gate(dt; evolution = :real, projection = :sandwich) ≈
          P * U * P
    @test projected_square_pxp_gate(dt; evolution = :real) ≈
          projected_square_pxp_gate(dt; evolution = :real, projection = :sandwich)
    @test_throws ArgumentError projected_square_pxp_gate(dt; projection = :bad)
end
