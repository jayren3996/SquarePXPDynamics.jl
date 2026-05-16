using LinearAlgebra

function _dense_square_star_index_obs(values)
    idx = 1
    for (site, value) in enumerate(values)
        idx += (value - 1) * 2^(SQUARE_STAR_SITES - site)
    end
    return idx
end

@testset "simple observables on product states" begin
    cell = PeriodicSquareUnitCell(10, 10)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    @test density_simple(down) ≈ 0 atol = 1e-14
    @test blockade_violation_simple(down) ≈ 0 atol = 1e-14
    @test pxp_energy_density_simple(down) ≈ 0 atol = 1e-14
    @test mean_bond_entropy(down) ≈ 0 atol = 1e-14
    @test max_bond_entropy(down) ≈ 0 atol = 1e-14

    up = product_square_ipeps(cell; state = :up, maxdim = 1)
    @test density_simple(up) ≈ 1 atol = 1e-14
    @test blockade_violation_simple(up) ≈ 1 atol = 1e-14
    @test pxp_energy_density_simple(up) ≈ 0 atol = 1e-14

    cb = checkerboard_square_ipeps(cell; excited_on = :even, maxdim = 1)
    dens = sublattice_densities(cb)
    @test dens.even ≈ 1 atol = 1e-14
    @test dens.odd ≈ 0 atol = 1e-14
    @test density_simple(cb) ≈ 0.5 atol = 1e-14
    @test blockade_violation_simple(cb) ≈ 0 atol = 1e-14
    @test pxp_energy_density_simple(cb) ≈ 0 atol = 1e-14
end

@testset "PXP x-plus star expectation matches dense product reference" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :x_plus, maxdim = 1)
    c = SquareCoord(5, 5)
    Hstar = square_pxp_star_hamiltonian()

    plus = fill(inv(sqrt(2)), 2)
    dense = zeros(ComplexF64, 2^SQUARE_STAR_SITES)
    for values in Iterators.product((1:2 for _ = 1:SQUARE_STAR_SITES)...)
        dense[_dense_square_star_index_obs(values)] = prod(plus[value] for value in values)
    end
    expected = dot(dense, Hstar * dense) / dot(dense, dense)

    @test star_expectation_simple(psi, c, Hstar) ≈ expected atol = 1e-12
end

@testset "local density after one star update" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    c = SquareCoord(5, 5)
    dt = 0.05
    project_star!(psi, c, dt; evolution = :real, projected = true, maxdim = 1)

    @test local_density_simple(psi, c) ≈ sin(dt)^2 atol = 1e-10 rtol = 1e-8
    for leaf in (
        neighbor(cell, c, :right),
        neighbor(cell, c, :up),
        neighbor(cell, c, :left),
        neighbor(cell, c, :down),
    )
        @test local_density_simple(psi, leaf) ≈ 0 atol = 1e-12
    end
end

@testset "star expectation simple matches dense D=1 reference" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    c = SquareCoord(5, 5)
    dt = 0.05
    Hstar = square_pxp_star_hamiltonian()

    before = zeros(ComplexF64, 2^SQUARE_STAR_SITES)
    before[_dense_square_star_index_obs((2, 2, 2, 2, 2))] = 1
    after = projected_square_pxp_gate(dt; evolution = :real) * before
    expected = dot(after, Hstar * after) / dot(after, after)

    project_star!(psi, c, dt; evolution = :real, projected = true, maxdim = 1)
    @test star_expectation_simple(psi, c, Hstar) ≈ expected atol = 1e-10 rtol = 1e-8
end

@testset "two-site blockade after local star update" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    c = SquareCoord(5, 5)
    project_star!(psi, c, 0.05; evolution = :real, projected = true, maxdim = 1)

    @test nearest_neighbor_density_simple(psi, c, :right) ≈ 0 atol = 1e-12
end

@testset "measure simple after short evolution" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    evolve!(psi, 0.01; params = params)

    summary = measure_simple(psi)
    @test summary isa SimpleObservableSummary
    @test all(
        isfinite,
        (
            summary.density,
            summary.density_even,
            summary.density_odd,
            summary.blockade_violation,
            summary.pxp_energy_density,
            summary.mean_bond_entropy,
            summary.max_bond_entropy,
        ),
    )
    @test 0 <= summary.density <= 1
    @test summary.blockade_violation >= 0
end

@testset "simple observables finite for D=2 short evolution" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    params = TrotterParams(0.01, 1, :real, true, 2, 1e-12)
    evolve!(psi, 0.01; params = params)

    summary = measure_simple(psi)
    @test isfinite(summary.density)
    @test isfinite(summary.blockade_violation)
    @test isfinite(summary.pxp_energy_density)
    @test isfinite(summary.mean_bond_entropy)
    @test isfinite(summary.max_bond_entropy)
end

@testset "simple observable validation" begin
    psi = product_square_ipeps(PeriodicSquareUnitCell(10, 10); state = :down, maxdim = 1)

    @test_throws ArgumentError density_simple(psi; sublattice = :bad)
    @test_throws ArgumentError nearest_neighbor_density_simple(psi, SquareCoord(1, 1), :bad)
    @test_throws ArgumentError SquarePXPDynamics.Observables._real_expectation(1 + 1e-3im)
end
