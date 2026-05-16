@testset "PXP scar observables simple product limits" begin
    cell = PeriodicSquareUnitCell(10, 10)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    checker = checkerboard_square_ipeps(cell; maxdim = 1)

    @test sublattice_imbalance_simple(down) ≈ 0.0 atol = 1e-12
    @test abs(sublattice_imbalance_simple(checker)) ≈ 1.0 atol = 1e-12
    @test checkerboard_structure_factor_simple(down) ≈ 0.0 atol = 1e-12
    @test checkerboard_structure_factor_simple(checker) ≈ 1.0 atol = 1e-12
end
