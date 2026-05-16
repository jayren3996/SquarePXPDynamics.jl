using Test
using LinearAlgebra
using SquarePXPDynamics

@testset "finite TFIM reference validates inputs" begin
    model = TFIMStarModel(1.0, 0.5)

    @test_throws ArgumentError finite_tfim_hamiltonian(PeriodicSquareUnitCell(4, 4), model)
    @test_throws ArgumentError finite_tfim_product_state(PeriodicSquareUnitCell(3, 3); state = :bad)
    @test_throws ArgumentError run_finite_tfim_reference(
        PeriodicSquareUnitCell(3, 3),
        model;
        initial_state = :z_up,
        total_time = 0.015,
        dt = 0.01,
    )
end

@testset "finite TFIM reference product observables match infinite product checks" begin
    cell = PeriodicSquareUnitCell(3, 3)
    model = TFIMStarModel(1.0, 0.5)

    up = finite_tfim_product_state(cell; state = :z_up)
    xplus = finite_tfim_product_state(cell; state = :x_plus)
    up_obs = measure_finite_tfim(up, cell, model)
    xplus_obs = measure_finite_tfim(xplus, cell, model)

    @test length(up) == 2^9
    @test norm(up) ≈ 1.0 atol = 1e-12
    @test up_obs.mean_x ≈ 0.0 atol = 1e-12
    @test up_obs.mean_z ≈ 1.0 atol = 1e-12
    @test up_obs.zz_right ≈ 1.0 atol = 1e-12
    @test up_obs.zz_up ≈ 1.0 atol = 1e-12
    @test up_obs.energy_density_decomposed ≈ -2.0 atol = 1e-12
    @test up_obs.energy_density_star ≈ up_obs.energy_density_decomposed atol = 1e-12
    @test up_obs.max_imag_abs ≤ 1e-12
    @test xplus_obs.mean_x ≈ 1.0 atol = 1e-12
    @test xplus_obs.mean_z ≈ 0.0 atol = 1e-12
    @test xplus_obs.energy_density_decomposed ≈ -0.5 atol = 1e-12
end

@testset "finite TFIM exact reference reproduces independent spin dynamics" begin
    cell = PeriodicSquareUnitCell(3, 3)
    h = 1.0
    dt = 0.01
    total_time = 0.02
    samples = run_finite_tfim_reference(
        cell,
        TFIMStarModel(0.0, h);
        initial_state = :z_up,
        total_time,
        dt,
        measure_every = 1,
    )

    @test [s.step for s in samples] == [0, 1, 2]
    @test samples[end].time ≈ total_time atol = 1e-12
    @test samples[end].observables.mean_z ≈ cos(2h * total_time) atol = 1e-12
    @test abs(abs(samples[end].observables.mean_y) - abs(sin(2h * total_time))) ≤ 1e-12
end
