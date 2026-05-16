using Test
using SquarePXPDynamics

function _validation_fake_ctm_summary(params; density, blockade, energy, accepted = true)
    return CTMObservableSummary(
        density,
        density,
        density,
        blockade,
        energy,
        CTMRGDiagnostics(
            params.chi,
            params.tol,
            params.maxiter,
            params.maxiter,
            params.tol / 10,
            true,
            accepted,
        ),
    )
end

@testset "trusted CTM measurement composes finite chi sweep and trust policy" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0),
    )

    trusted = measure_ctm_trusted(
        psi;
        params,
        measure = (state; params) -> _validation_fake_ctm_summary(
            params;
            density = 0.2 + params.chi * 1e-5,
            blockade = params.chi * 1e-6,
            energy = -0.1 - params.chi * 1e-5,
        ),
    )

    @test trusted isa TrustedCTMMeasurement
    @test length(trusted.points) == 2
    @test trusted.measurement === trusted.points[end].measurement
    @test trusted.trust.trusted === true
    @test trusted.trust.reason === :trusted
    @test trusted.measurement.diagnostics.chi == 4
end

@testset "trusted CTM measurement records rejected finite chi drift" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0),
    )

    rejected = measure_ctm_trusted(
        psi;
        params,
        measure = (state; params) -> _validation_fake_ctm_summary(
            params;
            density = 0.2,
            blockade = 0.0,
            energy = params.chi == 2 ? -0.1 : -0.2,
        ),
    )

    @test rejected.trust.trusted === false
    @test rejected.trust.reason === :energy_delta_too_large
    @test rejected.trust.finite_chi_energy_delta > CTMTrustPolicy().max_energy_delta
end
