using PEPSKit
using TensorKit

@testset "PEPSKit CTMRG measurement adapter" begin
    @testset "PEPSKit loads" begin
        @test isdefined(SquarePXPDynamics, :PEPSKitCTMRGParams)
    end

    @testset "CTMRG parameter validation" begin
        @test_throws ArgumentError PEPSKitCTMRGParams(0, 1e-8, 10, 0)
        @test_throws ArgumentError PEPSKitCTMRGParams(4, 0.0, 10, 0)
        @test_throws ArgumentError PEPSKitCTMRGParams(4, 1e-8, 0, 0)
        @test_throws ArgumentError PEPSKitCTMRGParams(4, 1e-8, 10, -1)
    end

    @testset "product-state conversion" begin
        cell = PeriodicSquareUnitCell(10, 10)
        psi = product_square_ipeps(cell; state = :down, maxdim = 1)
        before = copy(psi.tensors[SquareCoord(5, 5)])

        peps = to_pepskit_infinitepeps(psi)

        @test peps !== nothing
        @test size(peps) == (10, 10)
        @test PEPSKit.physicalspace(peps, 6, 5) == TensorKit.ComplexSpace(2)
        @test TensorKit.dim(PEPSKit.virtualspace(peps, 6, 5, PEPSKit.NORTH)) == 1
        @test TensorKit.dim(PEPSKit.virtualspace(peps, 6, 5, PEPSKit.EAST)) == 1
        @test TensorKit.dim(PEPSKit.virtualspace(peps, 6, 5, PEPSKit.SOUTH)) == 1
        @test TensorKit.dim(PEPSKit.virtualspace(peps, 6, 5, PEPSKit.WEST)) == 1
        @test psi.tensors[SquareCoord(5, 5)] == before
    end

    @testset "coordinate mapping catches transpose" begin
        cell = PeriodicSquareUnitCell(2, 3)
        @test SquarePXPDynamics.PEPSKitMeasurements._squarecoord_to_cartesianindex(
            cell,
            SquareCoord(1, 1),
        ) == CartesianIndex(3, 1)
        @test SquarePXPDynamics.PEPSKitMeasurements._squarecoord_to_cartesianindex(
            cell,
            SquareCoord(2, 3),
        ) == CartesianIndex(1, 2)
        @test SquarePXPDynamics.PEPSKitMeasurements._squarecoord_to_cartesianindex(
            cell,
            neighbor(cell, SquareCoord(1, 2), :up),
        ) == CartesianIndex(1, 1)
        @test SquarePXPDynamics.PEPSKitMeasurements._squarecoord_to_cartesianindex(
            cell,
            neighbor(cell, SquareCoord(1, 2), :right),
        ) == CartesianIndex(2, 2)
    end

    @testset "local operators" begin
        cell = PeriodicSquareUnitCell(10, 10)
        star = SquarePXPDynamics.PEPSKitMeasurements._star_sites_cartesian(
            cell,
            SquareCoord(5, 5),
        )
        @test star == (
            CartesianIndex(6, 5),
            CartesianIndex(6, 6),
            CartesianIndex(5, 5),
            CartesianIndex(6, 4),
            CartesianIndex(7, 5),
        )
        @test SquarePXPDynamics.PEPSKitMeasurements._pepskit_density_operator(
            cell,
            SquareCoord(5, 5),
        ) !== nothing
        @test SquarePXPDynamics.PEPSKitMeasurements._pepskit_twosite_nn_operator(
            cell,
            SquareCoord(5, 5),
            :right,
        ) !== nothing
        @test SquarePXPDynamics.PEPSKitMeasurements._pepskit_pxp_star_operator(
            cell,
            SquareCoord(5, 5),
        ) !== nothing
    end

    @testset "CTMRG product-state measurements" begin
        cell = PeriodicSquareUnitCell(3, 3)
        params = PEPSKitCTMRGParams(2, 1e-6, 20, 0)
        psi_down = product_square_ipeps(cell; state = :down, maxdim = 1)
        psi_up = product_square_ipeps(cell; state = :up, maxdim = 1)

        ctx_down = pepskit_ctmrg_context(psi_down; params)
        ctx_up = pepskit_ctmrg_context(psi_up; params)

        @test ctx_down isa PEPSKitMeasurementContext
        @test ctx_down.peps !== nothing
        @test ctx_down.env !== nothing
        @test local_density_ctm(psi_down, SquareCoord(2, 2), ctx_down) ≈ 0 atol = 1e-8
        @test local_density_ctm(psi_up, SquareCoord(2, 2), ctx_up) ≈ 1 atol = 1e-8
        @test blockade_violation_ctm(psi_down, ctx_down) ≈ 0 atol = 1e-8
        @test blockade_violation_ctm(psi_up, ctx_up) ≈ 1 atol = 1e-8
        @test pxp_energy_density_ctm(psi_down, ctx_down) ≈ 0 atol = 1e-8
        @test pxp_energy_density_ctm(psi_up, ctx_up) ≈ 0 atol = 1e-8
        @test_throws ArgumentError nearest_neighbor_density_ctm(
            psi_down,
            SquareCoord(1, 1),
            :bad,
            ctx_down,
        )
    end

    @testset "CTM and simple product observables agree" begin
        cell = PeriodicSquareUnitCell(3, 3)
        params = PEPSKitCTMRGParams(2, 1e-6, 20, 0)
        states = (
            product_square_ipeps(cell; state = :down, maxdim = 1),
            product_square_ipeps(cell; state = :up, maxdim = 1),
            checkerboard_square_ipeps(cell; excited_on = :even, maxdim = 1),
        )

        for psi in states
            summary_ctm = measure_ctm(psi; params)
            summary_simple = measure_simple(psi)
            @test summary_ctm.density ≈ summary_simple.density atol = 1e-8
            @test summary_ctm.blockade_violation ≈ summary_simple.blockade_violation atol =
                1e-8
            @test summary_ctm.pxp_energy_density ≈ summary_simple.pxp_energy_density atol =
                1e-8
        end
    end

    @testset "short-evolved CTM smoke test" begin
        cell = PeriodicSquareUnitCell(5, 5)
        psi = product_square_ipeps(cell; state = :down, maxdim = 1)
        evolve!(psi, 0.01; params = TrotterParams(0.01, 1, :real, true, 1, 1e-12))

        params = PEPSKitCTMRGParams(2, 1e-6, 20, 0)
        ctx = pepskit_ctmrg_context(psi; params)
        summary = measure_ctm(psi; params)

        @test ctx isa PEPSKitMeasurementContext
        @test isfinite(summary.density)
        @test isfinite(summary.blockade_violation)
        @test isfinite(summary.pxp_energy_density)
    end
end
