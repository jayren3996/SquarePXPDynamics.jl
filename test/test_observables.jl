@testset "PXP scar observables simple product limits" begin
    cell = PeriodicSquareUnitCell(10, 10)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    checker = checkerboard_square_ipeps(cell; maxdim = 1)

    @test sublattice_imbalance_simple(down) ≈ 0.0 atol = 1e-12
    @test abs(sublattice_imbalance_simple(checker)) ≈ 1.0 atol = 1e-12
    @test checkerboard_structure_factor_simple(down) ≈ 0.0 atol = 1e-12
    @test checkerboard_structure_factor_simple(checker) ≈ 1.0 atol = 1e-12
end

@testset "Local on-site observables on simple product states" begin
    cell = PeriodicSquareUnitCell(10, 10)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    up = product_square_ipeps(cell; state = :up, maxdim = 1)
    c = SquareCoord(1, 1)

    # |down>: <n>=0, <z>=-1, <x>=<y>=0 with basis 1=:up, 2=:down so n = (1-z)/2 ? Actually
    # density = <P_up> = 1 on up state, 0 on down state.
    @test local_density_simple(down, c) ≈ 0.0 atol = 1e-12
    @test local_density_simple(up, c) ≈ 1.0 atol = 1e-12

    @test local_z_simple(down, c) ≈ -1.0 atol = 1e-12
    @test local_z_simple(up, c) ≈ 1.0 atol = 1e-12

    @test local_x_simple(down, c) ≈ 0.0 atol = 1e-12
    @test local_x_simple(up, c) ≈ 0.0 atol = 1e-12
    @test local_y_simple(down, c) ≈ 0.0 atol = 1e-12
    @test local_y_simple(up, c) ≈ 0.0 atol = 1e-12
end

@testset "Direction validation rejects bogus symbols" begin
    cell = PeriodicSquareUnitCell(10, 10)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    c = SquareCoord(1, 1)

    @test_throws ArgumentError nearest_neighbor_density_simple(down, c, :diagonal)
    @test_throws ArgumentError nearest_neighbor_zz_simple(down, c, :diagonal)
end

@testset "Nearest-neighbor product observables" begin
    cell = PeriodicSquareUnitCell(10, 10)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    up = product_square_ipeps(cell; state = :up, maxdim = 1)
    c = SquareCoord(1, 1)

    # <n n> on all-down = 0; on all-up = 1.
    @test nearest_neighbor_density_simple(down, c, :right) ≈ 0.0 atol = 1e-12
    @test nearest_neighbor_density_simple(up, c, :right) ≈ 1.0 atol = 1e-12

    # <Z Z> on all-down = (-1)*(-1) = 1; on all-up = (+1)(+1) = 1.
    @test nearest_neighbor_zz_simple(down, c, :right) ≈ 1.0 atol = 1e-12
    @test nearest_neighbor_zz_simple(up, c, :right) ≈ 1.0 atol = 1e-12
end

@testset "PXP energy density on all-down product is zero" begin
    cell = PeriodicSquareUnitCell(10, 10)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)

    # X|down> = |up> but the all-down state has <up|down> = 0, so
    # <down|X|down> = 0, and any star with the projector annihilating
    # leaves us with E_PXP = 0.
    @test pxp_energy_density_simple(down) ≈ 0.0 atol = 1e-12
    @test blockade_violation_simple(down) ≈ 0.0 atol = 1e-12
end

@testset "Bond entropy of D=1 product is zero" begin
    cell = PeriodicSquareUnitCell(10, 10)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)

    @test mean_bond_entropy(down) ≈ 0.0 atol = 1e-12
    @test max_bond_entropy(down) ≈ 0.0 atol = 1e-12
end

@testset "measure_simple aggregates basic fields" begin
    cell = PeriodicSquareUnitCell(10, 10)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    summary = measure_simple(down)

    @test summary isa SimpleObservableSummary
    @test summary.density ≈ 0.0 atol = 1e-12
    @test summary.blockade_violation ≈ 0.0 atol = 1e-12
    @test summary.pxp_energy_density ≈ 0.0 atol = 1e-12
end
