using ITensors

@testset "product square iPEPS constructor" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)

    @test psi.unitcell === cell
    @test length(psi.tensors) == 100
    @test psi.maxdim == 1
    @test psi.gauge === :simple

    @test all(dim(physical_index(psi, c)) == 2 for c in cell.reps)
    @test all(length(lambda) == 1 for lambda in values(psi.link_weights))
    @test all(lambda[1] ≈ 1.0 for lambda in values(psi.link_weights))

    for c in cell.reps
        @test link_index(psi, c, :right) == link_index(psi, neighbor(cell, c, :right), :left)
        @test link_index(psi, c, :up) == link_index(psi, neighbor(cell, c, :up), :down)
    end

    @test density_simple(psi) ≈ 0 atol = 1e-14
    @test blockade_violation_simple(psi) ≈ 0 atol = 1e-14

    up = product_square_ipeps(cell; state = :up, maxdim = 1)
    @test density_simple(up) ≈ 1 atol = 1e-14
    @test blockade_violation_simple(up) ≈ 1 atol = 1e-14
end

@testset "checkerboard square iPEPS constructor" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = checkerboard_square_ipeps(cell; excited_on = :even, maxdim = 1)

    dens = sublattice_densities(psi)
    @test dens.even ≈ 1 atol = 1e-14
    @test dens.odd ≈ 0 atol = 1e-14
    @test density_simple(psi) ≈ 0.5 atol = 1e-14
    @test blockade_violation_simple(psi) ≈ 0 atol = 1e-14
end

@testset "iPEPS constructor validation and D>1 simple observables" begin
    cell = PeriodicSquareUnitCell(10, 10)

    @test_throws ArgumentError product_square_ipeps(cell; state = :bad, maxdim = 1)
    @test_throws ArgumentError product_square_ipeps(cell; state = :down, maxdim = 0)
    @test_throws ArgumentError checkerboard_square_ipeps(cell; excited_on = :bad, maxdim = 1)
    @test_throws ArgumentError product_square_ipeps(PeriodicSquareUnitCell(1, 10); state = :down, maxdim = 1)
    @test_throws ArgumentError product_square_ipeps(PeriodicSquareUnitCell(10, 1); state = :down, maxdim = 1)

    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    @test density_simple(psi) ≈ 0 atol = 1e-14
    @test blockade_violation_simple(psi) ≈ 0 atol = 1e-14
end
