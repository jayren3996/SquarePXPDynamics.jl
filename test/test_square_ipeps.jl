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
        @test link_index(psi, c, :right) ==
              link_index(psi, neighbor(cell, c, :right), :left)
        @test link_index(psi, c, :up) == link_index(psi, neighbor(cell, c, :up), :down)
    end

    @test density_simple(psi) ≈ 0 atol = 1e-14
    @test blockade_violation_simple(psi) ≈ 0 atol = 1e-14

    up = product_square_ipeps(cell; state = :up, maxdim = 1)
    @test density_simple(up) ≈ 1 atol = 1e-14
    @test blockade_violation_simple(up) ≈ 1 atol = 1e-14
end

@testset "square iPEPS public helper APIs" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)

    @test unitcell_reps(psi) == cell.reps
    @test unitcell_reps(psi) !== cell.reps
    @test physical_dim(psi, SquareCoord(1, 1)) == 2
    @test physical_dim(psi, SquareCoord(5, 5)) == 2
    @test simple_weight_dim(psi, SquareCoord(1, 1), :right) == 2
    @test simple_weight_dim(psi, SquareCoord(1, 1), :up) == 2
    @test simple_weight_dim(psi, SquareCoord(1, 1), :left) == 2
    @test simple_weight_dim(psi, SquareCoord(1, 1), :down) == 2
    @test_throws ArgumentError simple_weight_dim(psi, SquareCoord(1, 1), :diagonal)

    copied = copy_state(psi)
    @test copied !== psi
    @test copied.unitcell == psi.unitcell
    @test copied.maxdim == psi.maxdim
    @test copied.gauge == psi.gauge
    @test state_version(copied) == state_version(psi)
    @test log_norm(copied) == log_norm(psi)

    set_link_weight!(copied, SquareCoord(1, 1), :right, [0.6, 0.8])
    @test link_weight(copied, SquareCoord(1, 1), :right) == [0.6, 0.8]
    @test link_weight(psi, SquareCoord(1, 1), :right) == [1.0, 0.0]
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
    @test_throws ArgumentError checkerboard_square_ipeps(
        cell;
        excited_on = :bad,
        maxdim = 1,
    )
    @test_throws ArgumentError product_square_ipeps(
        PeriodicSquareUnitCell(1, 10);
        state = :down,
        maxdim = 1,
    )
    @test_throws ArgumentError product_square_ipeps(
        PeriodicSquareUnitCell(10, 1);
        state = :down,
        maxdim = 1,
    )

    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    @test density_simple(psi) ≈ 0 atol = 1e-14
    @test blockade_violation_simple(psi) ≈ 0 atol = 1e-14
end

@testset "benchmark product-state aliases" begin
    cell = PeriodicSquareUnitCell(2, 2)
    z_up = product_square_ipeps(cell; state = :z_up, maxdim = 1)
    z_down = product_square_ipeps(cell; state = :z_down, maxdim = 1)
    x_plus = product_square_ipeps(cell; state = :x_plus, maxdim = 1)

    @test local_density_simple(z_up, SquareCoord(1, 1)) ≈ 1.0 atol = 1e-12
    @test local_density_simple(z_down, SquareCoord(1, 1)) ≈ 0.0 atol = 1e-12
    @test local_density_simple(x_plus, SquareCoord(1, 1)) ≈ 0.5 atol = 1e-12
end
