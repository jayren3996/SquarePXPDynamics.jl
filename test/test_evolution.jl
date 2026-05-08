using LinearAlgebra
using ITensors

@testset "evolution" begin
    @testset "color canonical centers match requested color" begin
        for color in 1:7
            @test star_color(color_canonical_center(color)) == color
        end
    end

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

    @testset "Hamiltonian evolution uses first/second order step weights" begin
        α = 0.03
        Xsum = sum(embed_one_site(pauli_x(), site, 7) for site in 1:7)

        first = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        evolve_step!(first, Xsum, α; order = :first, update = :simple)
        @test real(local_expectation(first, Coord(0, 0), pauli_z())) ≈ -cos(14 * α) atol = 1e-10

        second = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        evolve_step!(second, Xsum, α; order = :second, update = :simple)
        @test real(local_expectation(second, Coord(0, 0), pauli_z())) ≈ -cos(14 * α) atol = 1e-10
    end

    @testset "invalid update mode" begin
        state = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        I128 = Matrix{ComplexF64}(I, 128, 128)
        @test_throws ArgumentError evolve_step!(state, I128; order = :first, update = :nonexistent)
        @test_throws ArgumentError evolve_step!(state, I128; order = :nonexistent, update = :simple)
    end

    @testset "projected PXP evolution smoke test" begin
        state = product_ipeps(ThreeSiteUnitCell(), :down; D = 1)
        H = pxp_star_hamiltonian(projector_down(), pauli_x())
        evolve_step!(state, H, 0.01; order = :first, update = :simple,
                     evolution = :real, projected = true)
        reps = collect(unit_cell_representatives(ThreeSiteUnitCell()))
        @test all(isfinite(real(local_expectation(state, c, pauli_z()))) for c in reps)
        @test mean_blockade_violation(state, reps) < 1e-8
    end

    @testset "projected PXP step helpers report schedule diagnostics" begin
        first = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        d1 = projected_pxp_step!(first, 0.01; order = :first, maxdim = 1, cutoff = 1e-12)
        @test length(d1.layer_diagnostics) == 7
        @test length(d1.discarded_weights) == 7
        @test d1.max_bond_dim == 1
        @test isfinite(d1.blockade_violation)
        @test all(isfinite, values(d1.local_z))
        @test all(isfinite, values(d1.local_x))
        @test all(isfinite, values(d1.local_projector_up))

        second = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        d2 = projected_pxp_step!(second, 0.01; order = :second, maxdim = 1)
        @test length(d2.layer_diagnostics) == 14
        @test length(d2.discarded_weights) == 14
    end

    @testset "imaginary projected PXP step and runs stay finite" begin
        state = product_ipeps(ThreeSiteUnitCell(), :down; D = 1)
        diag = imaginary_projected_pxp_step!(state, 0.01; order = :first, maxdim = 1)
        reps = unit_cell_representatives(ThreeSiteUnitCell())
        @test all(isfinite(tensor_norm(state, c)) for c in reps)
        @test all(isfinite, diag.tensor_norms)

        run_state = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        history = run_projected_pxp!(run_state, 0.005, 3; order = :first, maxdim = 1)
        @test length(history) == 3
        @test all(length(d.layer_diagnostics) == 7 for d in history)
    end

    @testset "projected PXP helpers reject unsupported D>1 non-product updates" begin
        state = random_ipeps(OneSiteUnitCell(), 2; seed = 31)
        @test_throws ArgumentError projected_pxp_step!(state, 0.01; order = :first, maxdim = 2)
    end

    @testset "projected PXP helper invalid options" begin
        state = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        @test_throws ArgumentError projected_pxp_step!(state, 0.01; order = :bad, maxdim = 1)
        @test_throws ArgumentError projected_pxp_step!(state, 0.01; evolution = :bad, maxdim = 1)
        @test_throws ArgumentError projected_pxp_step!(state, 0.01; maxdim = 0)
        @test_throws ArgumentError run_projected_pxp!(state, 0.01, -1; maxdim = 1)
    end
end
