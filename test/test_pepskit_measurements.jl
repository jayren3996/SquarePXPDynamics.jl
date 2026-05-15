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

    @testset "lambda absorption does not mutate source state" begin
        cell = PeriodicSquareUnitCell(3, 3)
        c = SquareCoord(2, 2)
        right = neighbor(cell, c, :right)
        psi = product_square_ipeps(cell; state = :down, maxdim = 2)
        before_center = copy(psi.tensors[c])
        before_right = copy(psi.tensors[right])

        set_link_weight!(psi, c, :right, [9.0, 4.0])

        center_data = SquarePXPDynamics.PEPSKitMeasurements._absorbed_site_array(psi, c)
        right_data =
            SquarePXPDynamics.PEPSKitMeasurements._absorbed_site_array(psi, right)

        @test center_data[2, 1, 1, 1, 1] ≈ 3.0 + 0.0im atol = 1e-12
        @test right_data[2, 1, 1, 1, 1] ≈ 3.0 + 0.0im atol = 1e-12
        @test psi.tensors[c] == before_center
        @test psi.tensors[right] == before_right
        @test link_weight(psi, c, :right) == [9.0, 4.0]
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
        @test_throws ArgumentError SquarePXPDynamics.PEPSKitMeasurements._pepskit_twosite_nn_operator(
            cell,
            SquareCoord(5, 5),
            :bad,
        )
        @test SquarePXPDynamics.PEPSKitMeasurements._pepskit_pxp_star_operator(
            cell,
            SquareCoord(5, 5),
        ) !== nothing
        @test_throws ArgumentError SquarePXPDynamics.PEPSKitMeasurements._pepskit_pxp_star_operator(
            PeriodicSquareUnitCell(2, 2),
            SquareCoord(1, 1),
        )
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
        @test_throws ArgumentError SquarePXPDynamics.PEPSKitMeasurements._local_neighbor_cartesianindex(
            cell,
            SquareCoord(1, 1),
            :bad,
        )
    end

    @testset "CTM and simple product observables agree" begin
        product_cell = PeriodicSquareUnitCell(3, 3)
        checkerboard_cell = PeriodicSquareUnitCell(4, 4)
        params = PEPSKitCTMRGParams(2, 1e-6, 20, 0)
        cases = (
            (
                psi = product_square_ipeps(product_cell; state = :down, maxdim = 1),
                density = 0.0,
                blockade = 0.0,
                energy = 0.0,
            ),
            (
                psi = product_square_ipeps(product_cell; state = :up, maxdim = 1),
                density = 1.0,
                blockade = 1.0,
                energy = 0.0,
            ),
            (
                psi = checkerboard_square_ipeps(
                    checkerboard_cell;
                    excited_on = :even,
                    maxdim = 1,
                ),
                density = 0.5,
                blockade = 0.0,
                energy = 0.0,
            ),
        )

        for case in cases
            psi = case.psi
            summary_ctm = measure_ctm(psi; params)
            summary_simple = measure_simple(psi)
            @test summary_ctm.density ≈ summary_simple.density atol = 1e-8
            @test summary_ctm.blockade_violation ≈ summary_simple.blockade_violation atol =
                1e-8
            @test summary_ctm.pxp_energy_density ≈ summary_simple.pxp_energy_density atol =
                1e-8
            @test summary_ctm.density ≈ case.density atol = 1e-8
            @test summary_ctm.blockade_violation ≈ case.blockade atol = 1e-8
            @test summary_ctm.pxp_energy_density ≈ case.energy atol = 1e-8
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
