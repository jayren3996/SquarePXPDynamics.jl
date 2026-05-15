using LinearAlgebra

@testset "star model conventions" begin
    @test star_site_order() == (:center, :right, :up, :left, :down)
    @test tfim_pauli_convention() == (:Z_up_is_plus_one, :X_field)
end

@testset "PXP star model reproduces existing gates" begin
    dt = 0.037
    @test star_hamiltonian(PXPStarModel(false)) ≈ square_pxp_star_hamiltonian()
    @test star_gate(PXPStarModel(false), dt; evolution = :real) ≈
          square_pxp_gate(dt; evolution = :real)
    @test star_gate(PXPStarModel(true), dt; evolution = :real) ≈
          projected_square_pxp_gate(dt; evolution = :real)
    @test star_gate(PXPStarModel(false), dt; evolution = :imaginary) ≈
          square_pxp_gate(dt; evolution = :imaginary)
    @test star_gate(PXPStarModel(true), dt; evolution = :imaginary) ≈
          projected_square_pxp_gate(dt; evolution = :imaginary)
end

@testset "TFIM dense star Hamiltonian convention" begin
    model = TFIMStarModel(2.0, 3.0)
    H = star_hamiltonian(model)
    @test size(H) == (32, 32)
    @test H ≈ H'
    @test tfim_product_basis_energy(model, (:up, :up, :up, :up, :up)) ≈ -4.0
    @test tfim_product_basis_energy(model, (:up, :down, :up, :down, :up)) ≈ 0.0
    @test tfim_product_basis_energy(model, (:down, :up, :up, :up, :up)) ≈ 4.0
end

@testset "TFIM dense gates" begin
    model = TFIMStarModel(1.0, 0.7)
    dt = 0.011
    U = star_gate(model, dt; evolution = :real)
    G = star_gate(model, dt; evolution = :imaginary)
    @test U' * U ≈ I atol = 1e-12 rtol = 1e-12
    @test all(isfinite, real.(G))
    @test norm(G' * G - I) > 1e-6
    @test star_gate(model, dt; evolution = :real) *
          star_gate(model, dt; evolution = :real) ≈
          star_gate(model, 2dt; evolution = :real) atol = 1e-12 rtol = 1e-12
    @test_throws ArgumentError star_gate(model, dt; evolution = :bad)
    @test_throws ArgumentError star_gate(model, Inf; evolution = :real)
    @test_throws ArgumentError star_gate(model, NaN; evolution = :imaginary)
    @test_throws ArgumentError TFIMStarModel(Inf, 1.0)
    @test_throws ArgumentError TFIMStarModel(1.0, NaN)
end

@testset "static model protocol" begin
    model = TFIMStarModel(1.0, 2.0)
    protocol = StaticModel(model)
    @test model_at(protocol, 0.0, 1) === model
    @test model_at(protocol, 1.0, 17) === model
end
