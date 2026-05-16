using Test
using SquarePXPDynamics

@testset "finite MPS square lattice bonds use open boundaries" begin
    bonds = finite_mps_square_lattice_bonds(3, 2)

    @test finite_mps_site_index(1, 1, 3, 2) == 1
    @test finite_mps_site_index(3, 1, 3, 2) == 3
    @test finite_mps_site_index(1, 2, 3, 2) == 6
    @test finite_mps_site_index(3, 2, 3, 2) == 4
    @test length(bonds) == 7
    @test (1, 2, :right) in bonds
    @test (2, 3, :right) in bonds
    @test (1, 6, :up) in bonds
    @test (3, 4, :up) in bonds
    @test_throws ArgumentError finite_mps_site_index(0, 1, 3, 2)
    @test_throws ArgumentError finite_mps_square_lattice_bonds(1, 2)
end

@testset "finite MPS TFIM benchmark validates inputs" begin
    @test_throws ArgumentError run_finite_mps_tfim_reference(1, 6; total_time = 0.1)
    @test_throws ArgumentError run_finite_mps_tfim_reference(2, 2; initial_state = :bad)
    @test_throws ArgumentError run_finite_mps_tfim_reference(2, 2; total_time = 0.015, dt = 0.01)
    @test_throws ArgumentError run_finite_mps_tfim_reference(2, 2; maxdim = 0)
end

@testset "finite MPS TFIM reference matches exact independent spin dynamics" begin
    result = run_finite_mps_tfim_reference(
        2,
        2;
        J = 0.0,
        h = 1.0,
        initial_state = :z_up,
        total_time = 0.02,
        dt = 0.01,
        measure_every = 1,
        maxdim = 8,
        cutoff = 1e-12,
    )

    @test result.metadata.Lx == 2
    @test result.metadata.Ly == 2
    @test result.metadata.boundary == :open
    @test result.metadata.method == :tdvp
    @test [sample.step for sample in result.samples] == [0, 1, 2]
    @test result.samples[end].time ≈ 0.02 atol = 1e-12
    @test result.samples[end].mean_z ≈ cos(0.04) atol = 1e-8 rtol = 1e-8
    @test abs(abs(result.samples[end].mean_y) - abs(sin(0.04))) ≤ 1e-8
    @test result.samples[end].maxlinkdim <= 8
    @test isfinite(result.samples[end].energy_density)
end
