using LinearAlgebra
using ITensors

@testset "observables" begin
    @testset "local Z on product :down (1-site UC)" begin
        state = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        c0 = Coord(0, 0)
        @test local_expectation(state, c0, pauli_z()) ≈ -1
        @test local_expectation(state, c0, pauli_x()) ≈ 0
        @test local_expectation(state, c0, projector_down()) ≈ 1
        @test local_expectation(state, c0, projector_up()) ≈ 0
    end

    @testset "local Z on product :up" begin
        state = product_ipeps(OneSiteUnitCell(), :up; D = 1)
        c0 = Coord(0, 0)
        @test local_expectation(state, c0, pauli_z()) ≈ 1
        @test local_expectation(state, c0, projector_up()) ≈ 1
    end

    @testset "tensor_norm of product :down" begin
        state = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        @test tensor_norm(state, Coord(0, 0)) ≈ 1
    end

    @testset "tensor_norm of 3-site product is 1 per rep" begin
        state = product_ipeps(ThreeSiteUnitCell(), :up; D = 1)
        for c in unit_cell_representatives(ThreeSiteUnitCell())
            @test tensor_norm(state, c) ≈ 1
        end
    end

    @testset "dense_blockade_violations on cluster vectors" begin
        all_down = zeros(ComplexF64, 128)
        all_down[end] = 1
        @test dense_blockade_violations(all_down) ≈ 0

        # |up,up,...> violates every center-neighbor bond
        all_up = zeros(ComplexF64, 128)
        all_up[1] = 1
        @test dense_blockade_violations(all_up) > 0
    end

    @testset "nearest-neighbor blockade diagnostics on product states" begin
        down = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        @test local_blockade_violation(down, Coord(0, 0), 1) ≈ 0
        @test mean_blockade_violation(down, [Coord(0, 0)]) ≈ 0

        up = product_ipeps(OneSiteUnitCell(), :up; D = 1)
        @test local_blockade_violation(up, Coord(0, 0), 1) ≈ 1
        @test mean_blockade_violation(up, [Coord(0, 0)]) ≈ 1
    end
end
