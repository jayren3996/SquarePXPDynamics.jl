using LinearAlgebra
using ITensors

@testset "evolution" begin
    @testset "evolve_step! identity preserves D=1 product" begin
        state = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        I128 = Matrix{ComplexF64}(I, 128, 128)
        for order in (:first, :second)
            evolve_step!(state, I128; order = order, update = :simple)
            @test local_expectation(state, Coord(0, 0), pauli_z()) ≈ -1
        end
    end

    @testset "evolve_step! site-symmetric u^⊗7 first-order applies u 7 times" begin
        # u acts on |down⟩ as a small Z rotation: u = exp(-i α Z / 2)
        α = 0.1
        u = cos(α / 2) * Matrix{ComplexF64}(I, 2, 2) - im * sin(α / 2) * pauli_z()
        G = u
        for _ in 2:7
            G = kron(G, u)
        end
        state = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        evolve_step!(state, G; order = :first, update = :simple)
        # Z is preserved exactly (Z rotation about Z axis)
        @test real(local_expectation(state, Coord(0, 0), pauli_z())) ≈ -1 atol = 1e-10
    end

    @testset "invalid update mode" begin
        state = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        I128 = Matrix{ComplexF64}(I, 128, 128)
        @test_throws ArgumentError evolve_step!(state, I128; order = :first, update = :nonexistent)
        @test_throws ArgumentError evolve_step!(state, I128; order = :nonexistent, update = :simple)
    end
end
