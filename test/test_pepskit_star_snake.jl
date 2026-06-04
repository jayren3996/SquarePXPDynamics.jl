# M7 (experimental): star-gate MPO decomposition.
#
# The five-site square-star PXP gate is represented exactly as a matrix product
# operator along the star-site ordering. These tests pin the core correctness
# claim — exact reconstruction of the dense gate — and document the gate's small
# MPO bond structure. The MPO is NOT wired into production evolution.
using Test
using SquarePXPDynamics
using SquarePXPDynamics: pxp_star_gate_dense
using LinearAlgebra: I

@testset "M7 star-gate MPO (experimental)" begin
    @testset "exact reconstruction" begin
        for projected in (false, true), evolution in (:real, :imaginary), dt in (0.0, 0.17, 0.4)
            mpo = star_gate_mpo(dt; projected, evolution)
            G = reconstruct_star_gate(mpo)
            Gref = pxp_star_gate_dense(dt; projected, evolution)
            @test isapprox(G, Gref; atol = 1e-11)
        end
    end

    @testset "MPO structure" begin
        # dt = 0 (real, unprojected) is the identity gate: a product operator,
        # so every internal bond is dimension 1.
        mpo0 = star_gate_mpo(0.0; projected = false)
        @test reconstruct_star_gate(mpo0) ≈ Matrix{ComplexF64}(I, 32, 32)
        @test all(==(1), star_mpo_bond_dims(mpo0))
        @test length(mpo0.tensors) == 5
        @test length(star_mpo_bond_dims(mpo0)) == 4

        # A finite-time gate has bounded MPO rank: H_star = X_c⊗P^4 has
        # (X_c⊗P^4)^2 = I_c⊗P^4, so exp(-iθ H) = I + (cosθ-1)P^4 - i sinθ X_cP^4,
        # i.e. at most two operator terms across any star-site cut → bond ≤ 2.
        mpo = star_gate_mpo(0.37; projected = false)
        @test all(<=(2), star_mpo_bond_dims(mpo))
        @test maximum(star_mpo_bond_dims(mpo)) == 2   # genuinely entangling
    end

    @testset "leg conventions" begin
        mpo = star_gate_mpo(0.25; projected = true)
        @test length(mpo.tensors) == 5
        for W in mpo.tensors
            @test size(W, 1) == 2 && size(W, 2) == 2   # phys_out, phys_in
        end
        @test size(mpo.tensors[1], 3) == 1             # trivial outer-left bond
        @test size(mpo.tensors[end], 4) == 1           # trivial outer-right bond
    end
end
