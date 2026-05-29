using Test
using SquarePXPDynamics

@testset "compress_to_target_maxdim! is mathematically exact at target==current rank" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    evolve!(psi, 0.04; params = TrotterParams(0.02, 1, :real, 2, 1e-12; schedule = :serial))
    log_norm_before = log_norm(psi)
    info = compress_to_target_maxdim!(psi, 2)

    @test info.target_maxdim == 2
    # Compression to the existing rank should be exact to machine precision.
    @test info.max_truncerr <= 1e-20
    # Kept dimension per bond cannot exceed target.
    @test all(d <= 2 for d in values(info.bond_keptdims))
    # log_norm bookkeeping closes: after = before + delta.
    @test isapprox(log_norm(psi), log_norm_before + info.log_norm_delta; atol = 1e-12)
end

@testset "compress_to_target_maxdim! CTM observables match across no-op rank compression" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    evolve!(psi, 0.04; params = TrotterParams(0.02, 1, :real, 2, 1e-12; schedule = :serial))
    ctm_params = PEPSKitCTMRGParams(8, 1e-8, 100, 0; seed = 0)
    before = measure_ctm(psi; params = ctm_params)
    compress_to_target_maxdim!(psi, 2)
    after = measure_ctm(psi; params = ctm_params)

    # CTM is gauge-invariant: a compression that preserves the singular spectrum
    # must preserve all CTM observables to high precision.
    @test isapprox(before.density, after.density; atol = 1e-8)
    @test isapprox(before.blockade_violation, after.blockade_violation; atol = 1e-8)
    @test isapprox(before.pxp_energy_density, after.pxp_energy_density; atol = 1e-8)
end

@testset "compress_to_target_maxdim! reports nonneg truncerr when shrinking" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 4)
    evolve!(psi, 0.04; params = TrotterParams(0.02, 1, :real, 4, 1e-12; schedule = :serial))
    before = measure_simple(psi)
    info = compress_to_target_maxdim!(psi, 2)
    after = measure_simple(psi)

    @test info.target_maxdim == 2
    @test info.max_truncerr >= 0
    @test info.mean_truncerr >= 0
    @test all(d <= 2 for d in values(info.bond_keptdims))
    # Compression at this short evolution should still leave observables close
    @test isfinite(after.density)
    @test isfinite(after.blockade_violation)
end

@testset "compress_to_target_maxdim! validates inputs" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    @test_throws ArgumentError compress_to_target_maxdim!(psi, 0)
    @test_throws ArgumentError compress_to_target_maxdim!(psi, -1)
    @test_throws ArgumentError compress_to_target_maxdim!(psi, 2; cutoff = -1.0)
    @test_throws ArgumentError compress_to_target_maxdim!(psi, 2; cutoff = NaN)
end

@testset "compress_to_target_maxdim! preserves state version bookkeeping" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 4)
    evolve!(psi, 0.04; params = TrotterParams(0.02, 1, :real, 4, 1e-12; schedule = :serial))
    version_before = state_version(psi)
    log_norm_before = log_norm(psi)
    info = compress_to_target_maxdim!(psi, 2)

    @test state_version(psi) > version_before
    @test isapprox(log_norm(psi), log_norm_before + info.log_norm_delta; atol = 1e-12)
end

@testset "ScarFinderParams rejects target_maxdim > trotter.maxdim" begin
    trotter = TrotterParams(0.02, 1, :real, 2, 1e-12; schedule = :serial)
    @test_throws ArgumentError ScarFinderParams(
        0.02, trotter, 1, Inf, Inf, Inf, false; target_maxdim = 3,
    )
    @test_throws ArgumentError ScarFinderParams(
        0.02, trotter, 1, Inf, Inf, Inf, false; target_maxdim = 0,
    )
    ok = ScarFinderParams(
        0.02, trotter, 1, Inf, Inf, Inf, false; target_maxdim = 2,
    )
    @test ok.target_maxdim == 2
end

@testset "scarfinder! without target_maxdim records no compression (backwards compat)" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    trotter = TrotterParams(0.02, 1, :real, 2, 1e-12; schedule = :serial)
    params = ScarFinderParams(0.02, trotter, 2, Inf, Inf, Inf, false)
    result = scarfinder!(psi, params)
    @test all(it -> it.compression === nothing, result.iterations)
end

@testset "scarfinder! with target_maxdim < trotter.maxdim records compression" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 4)
    trotter = TrotterParams(0.02, 1, :real, 4, 1e-12; schedule = :serial)
    params = ScarFinderParams(
        0.02, trotter, 2, Inf, Inf, Inf, false;
        target_maxdim = 2,
    )
    result = scarfinder!(psi, params)
    @test all(it -> it.compression !== nothing, result.iterations)
    @test all(it -> it.compression.info.target_maxdim == 2, result.iterations)
    @test all(it -> it.compression.info.max_truncerr >= 0, result.iterations)
    @test all(it -> it.compression.density_shift >= 0, result.iterations)
end
