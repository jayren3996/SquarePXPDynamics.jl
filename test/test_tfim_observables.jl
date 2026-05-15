@testset "TFIM product-state observables" begin
    cell = PeriodicSquareUnitCell(10, 10)
    c = SquareCoord(5, 5)
    model = TFIMStarModel(1.0, 0.5)

    up = product_square_ipeps(cell; state = :z_up, maxdim = 1)
    down = product_square_ipeps(cell; state = :z_down, maxdim = 1)
    xplus = product_square_ipeps(cell; state = :x_plus, maxdim = 1)

    @test local_z_simple(up, c) ≈ 1.0 atol = 1e-12
    @test local_z_simple(down, c) ≈ -1.0 atol = 1e-12
    @test local_x_simple(xplus, c) ≈ 1.0 atol = 1e-12
    @test local_y_simple(xplus, c) ≈ 0.0 atol = 1e-12
    @test nearest_neighbor_zz_simple(up, c, :right) ≈ 1.0 atol = 1e-12
    @test nearest_neighbor_zz_simple(up, c, :up) ≈ 1.0 atol = 1e-12
    @test_throws ArgumentError nearest_neighbor_zz_simple(up, c, :left)
    @test tfim_energy_density_decomposed_simple(up, model) ≈ -2.0 atol = 1e-12
    @test tfim_energy_density_star_simple(up, model) ≈ -2.0 atol = 1e-12
    @test tfim_energy_density_decomposed_simple(xplus, model) ≈ -0.5 atol = 1e-12

    summary = measure_tfim_simple(up, model)
    @test summary.mean_x ≈ 0.0 atol = 1e-12
    @test summary.mean_y ≈ 0.0 atol = 1e-12
    @test summary.mean_z ≈ 1.0 atol = 1e-12
    @test summary.z_even ≈ 1.0 atol = 1e-12
    @test summary.z_odd ≈ 1.0 atol = 1e-12
    @test summary.zz_right ≈ 1.0 atol = 1e-12
    @test summary.zz_up ≈ 1.0 atol = 1e-12
    @test summary.energy_density_discrepancy ≈ 0.0 atol = 1e-12
    @test summary.zz_imag_abs ≤ 1e-12
    @test summary.max_imag_abs ≤ 1e-10
    @test summary.mean_bond_entropy ≈ 0.0 atol = 1e-12
    @test summary.max_bond_entropy ≈ 0.0 atol = 1e-12
end

@testset "TFIM h=0 stationary product observables" begin
    cell = PeriodicSquareUnitCell(10, 10)
    model = TFIMStarModel(1.0, 0.0)
    psi = product_square_ipeps(cell; state = :z_up, maxdim = 1)
    before = measure_tfim_simple(psi, model)
    evolve!(
        psi,
        0.01;
        params = TrotterParams(0.01, 1, :real, 1, 1e-12),
        protocol = StaticModel(model),
    )
    after = measure_tfim_simple(psi, model)
    @test after.mean_z ≈ before.mean_z atol = 1e-10
    @test after.zz_right ≈ before.zz_right atol = 1e-10
    @test after.zz_up ≈ before.zz_up atol = 1e-10
    @test after.energy_density_star ≈ before.energy_density_star atol = 1e-10
end

@testset "TFIM J=0 independent spin dynamics" begin
    cell = PeriodicSquareUnitCell(10, 10)
    h = 1.0
    model = TFIMStarModel(0.0, h)
    psi = product_square_ipeps(cell; state = :z_up, maxdim = 1)
    total_time = 0.02
    evolve!(
        psi,
        total_time;
        params = TrotterParams(0.01, 1, :real, 1, 1e-12),
        protocol = StaticModel(model),
    )
    summary = measure_tfim_simple(psi, model)
    @test summary.mean_z ≈ cos(2h * total_time) atol = 1e-6 rtol = 1e-6
    @test abs(abs(summary.mean_y) - abs(sin(2h * total_time))) ≤ 1e-6
end
