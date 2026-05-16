using Test
using JSON3
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
    policy = CTMTrustPolicy(2, true, 1e-2, 1e-3, 1e-2, 1e-4)

    trusted = measure_ctm_trusted(
        psi;
        params,
        policy,
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
    @test trusted.policy.min_points == policy.min_points
    @test trusted.policy.require_accepted_diagnostics == policy.require_accepted_diagnostics
    @test trusted.policy.max_density_delta == policy.max_density_delta
    @test trusted.policy.max_blockade_delta == policy.max_blockade_delta
    @test trusted.policy.max_energy_delta == policy.max_energy_delta
    @test trusted.policy.max_residual == policy.max_residual
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

@testset "PXP validation config validates ED and iPEPS controls" begin
    config = PXPValidationConfig(3; total_time = 0.02, dt = 0.01, measure_every = 1)

    @test config.n == 3
    @test config.total_time ≈ 0.02
    @test config.dt ≈ 0.01
    @test config.measure_every == 1
    @test config.initial_state === :down
    @test config.order == 1
    @test config.maxdim == 1
    @test config.schedule === :serial

    @test_throws ArgumentError PXPValidationConfig(2; total_time = 0.02, dt = 0.01)
    @test_throws ArgumentError PXPValidationConfig(3; total_time = 0.025, dt = 0.01)
    @test_throws ArgumentError PXPValidationConfig(3; total_time = 0.02, dt = 0.0)
    @test_throws ArgumentError PXPValidationConfig(3; total_time = 0.02, dt = 0.01, order = 3)
    @test_throws ArgumentError PXPValidationConfig(3; total_time = 0.02, dt = 0.01, schedule = :five_color)
end

@testset "ED-vs-iPEPS validation report samples matched times" begin
    config = PXPValidationConfig(3; total_time = 0.02, dt = 0.01, measure_every = 1)
    report = validate_pxp_ed_ipeps(config; ctm_params = nothing)

    @test report isa PXPValidationReport
    @test report.config === config
    @test report.ed_result.lattice_size == (3, 3)
    @test length(report.ed_result.samples) == 3
    @test length(report.ipeps_samples) == 3
    @test length(report.comparisons) == 3
    @test [s.step for s in report.ipeps_samples] == [0, 1, 2]
    @test [s.time for s in report.ipeps_samples] ≈ [0.0, 0.01, 0.02] atol = 1e-12
    @test report.ipeps_samples[1].evolution === nothing
    @test report.ipeps_samples[2].evolution isa EvolutionLog
    @test report.ipeps_samples[1].ctm === nothing
    @test report.comparisons[1].density_error_simple ≈ 0.0 atol = 1e-12
    @test all(c -> isfinite(c.density_error_simple), report.comparisons)
    @test all(c -> c.ipeps_ctm_density === nothing, report.comparisons)
    @test report.metadata.julia_version == string(VERSION)
end

@testset "PXP validation metadata uses package checkout git commit" begin
    config = PXPValidationConfig(3; total_time = 0.0, dt = 0.01, measure_every = 1)
    package_root = abspath(joinpath(@__DIR__, ".."))
    expected_commit = chomp(read(`git -C $package_root rev-parse HEAD`, String))

    cd(mktempdir()) do
        report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
        @test report.metadata.git_commit == expected_commit
    end
end

@testset "ED-vs-iPEPS validation can attach trusted fake CTM" begin
    config = PXPValidationConfig(3; total_time = 0.01, dt = 0.01, measure_every = 1)
    params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0),
    )

    report = validate_pxp_ed_ipeps(
        config;
        ctm_params = params,
        ctm_measure = (state; params) -> begin
            simple = measure_simple(state)
            return _validation_fake_ctm_summary(
                params;
                density = simple.density + params.chi * 1e-5,
                blockade = simple.blockade_violation + params.chi * 1e-6,
                energy = simple.pxp_energy_density - params.chi * 1e-5,
            )
        end,
    )

    @test all(sample -> sample.ctm isa TrustedCTMMeasurement, report.ipeps_samples)
    @test all(comparison -> comparison.ctm_trusted === true, report.comparisons)
    @test all(comparison -> comparison.ctm_reason === :trusted, report.comparisons)
    @test all(comparison -> comparison.ipeps_ctm_density !== nothing, report.comparisons)
    for (sample, comparison) in zip(report.ipeps_samples, report.comparisons)
        expected_density = sample.simple.density + 4e-5
        expected_blockade = sample.simple.blockade_violation + 4e-6
        @test comparison.ipeps_ctm_density ≈ expected_density atol = 1e-12
        @test comparison.density_error_ctm ≈
            expected_density - comparison.ed_excitation_density atol = 1e-12
        @test comparison.ctm_blockade_violation ≈ expected_blockade atol = 1e-12
    end
end

@testset "PXP reversibility report measures forward and reverse drift" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = TrotterParams(0.01, 1, :real, 1, 1e-12; schedule = :serial)

    report = validate_pxp_reversibility(psi, 0.01; params)

    @test report isa PXPReversibilityReport
    @test report.before isa SimpleObservableSummary
    @test report.after_forward isa SimpleObservableSummary
    @test report.after_reverse isa SimpleObservableSummary
    @test report.forward_log isa EvolutionLog
    @test report.reverse_log isa EvolutionLog
    @test isfinite(report.density_drift)
    @test isfinite(report.blockade_drift)
    @test isfinite(report.energy_drift)
    @test report.density_drift >= 0
    @test report.blockade_drift >= 0
    @test report.energy_drift >= 0
    @test report.density_drift <= 1e-10
    @test report.blockade_drift <= 1e-10
    @test report.energy_drift <= 1e-10
end

@testset "PXP validation report writes JSON artifact" begin
    config = PXPValidationConfig(3; total_time = 0.01, dt = 0.01, measure_every = 1)
    report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
    path = tempname() * ".json"

    written = write_pxp_validation_json(report, path)
    data = read(path, String)
    parsed = JSON3.read(data)

    @test written == path
    @test endswith(data, '\n')
    @test parsed.config.initial_state == "down"
    @test parsed.config.schedule == "serial"
    @test parsed.ed_result.diagnostics.matvecs > 0
    @test length(parsed.ed_result.diagnostics.accepted_intervals) >= 0
    @test parsed.ipeps_samples[1].ctm === nothing
    @test parsed.comparisons[1].ctm_reason === nothing
end

@testset "PXP validation report writes CTM trust policy JSON" begin
    config = PXPValidationConfig(3; total_time = 0.0, dt = 0.01, measure_every = 1)
    params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0; seed = 1234),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0; seed = 1234),
    )
    policy = CTMTrustPolicy(2, true, 1e-2, 1e-3, 1e-2, 1e-4)
    report = validate_pxp_ed_ipeps(
        config;
        ctm_params = params,
        trust_policy = policy,
        ctm_measure = (state; params) -> _validation_fake_ctm_summary(
            params;
            density = 0.2 + params.chi * 1e-5,
            blockade = params.chi * 1e-6,
            energy = -0.1 - params.chi * 1e-5,
        ),
    )
    path = tempname() * ".json"

    write_pxp_validation_json(report, path)
    parsed = JSON3.read(read(path, String))
    trust = parsed.ipeps_samples[1].ctm.trust

    @test trust.reason == "trusted"
    @test length(parsed.ipeps_samples[1].ctm.points) == 2
    @test parsed.ipeps_samples[1].ctm.points[1].seed == 1234
    @test parsed.ipeps_samples[1].ctm.points[2].seed == 1234
    @test parsed.ipeps_samples[1].ctm.measurement.diagnostics.chi == 4
    @test parsed.ipeps_samples[1].ctm.measurement.sublattice_imbalance == 0
    @test parsed.ipeps_samples[1].ctm.measurement.checkerboard_structure_factor == 0
    @test trust.policy.min_points == policy.min_points
    @test trust.policy.max_density_delta == policy.max_density_delta
    @test trust.policy.max_residual == policy.max_residual
end

@testset "PXP convergence report aggregates validation grid" begin
    base = PXPValidationConfig(3; total_time = 0.01, dt = 0.01, measure_every = 1)
    @test_throws ArgumentError PXPConvergenceConfig(
        base;
        dt_values = [0.01],
        D_values = [1],
        chi_values = Int[],
        cutoff_values = [0.0],
    )
    @test_throws ArgumentError PXPConvergenceConfig(
        base;
        dt_values = [0.01],
        D_values = [1],
        chi_values = [0],
        cutoff_values = [1e-12],
    )
    sweep = PXPConvergenceConfig(
        base;
        dt_values = [0.01, 0.005],
        D_values = [1],
        chi_values = Int[],
        cutoff_values = [1e-12],
    )

    report = validate_pxp_convergence(
        sweep;
        ctm_measure = (state; params) -> CTMObservableSummary(0.0, 0.0, 0.0, 0.0, 0.0),
    )

    @test length(report.runs) == 2
    @test report.runs[1].config.dt == 0.01
    @test report.runs[2].config.dt == 0.005
    @test isfinite(report.max_abs_density_error_simple)
end

@testset "PXP convergence report aggregates CTM chi sweep" begin
    base = PXPValidationConfig(3; total_time = 0.01, dt = 0.01, measure_every = 1)
    sweep = PXPConvergenceConfig(
        base;
        dt_values = [0.01],
        D_values = [1],
        chi_values = [2, 4],
        cutoff_values = [1e-12],
    )

    report = validate_pxp_convergence(
        sweep;
        ctm_measure = (state; params) -> _validation_fake_ctm_summary(
            params;
            density = measure_simple(state).density + (params.chi == 2 ? 0.5 : 0.0),
            blockade = 0.0,
            energy = -0.1 - params.chi * 1e-5,
        ),
    )

    @test report.max_abs_density_error_ctm !== nothing
    @test report.all_ctm_trusted !== nothing
    @test all(comparison -> comparison.density_error_ctm !== nothing, report.runs[1].comparisons)
    final_chi_max = maximum(abs(c.density_error_ctm) for c in report.runs[1].comparisons)
    @test report.max_abs_density_error_ctm > final_chi_max + 0.1
end

@testset "PXP convergence report writes JSON artifact" begin
    base = PXPValidationConfig(3; total_time = 0.01, dt = 0.01, measure_every = 1)
    sweep = PXPConvergenceConfig(
        base;
        dt_values = [0.01, 0.005],
        D_values = [1],
        chi_values = Int[],
        cutoff_values = [1e-12],
    )
    report = validate_pxp_convergence(sweep)
    path = tempname() * ".json"

    written = write_pxp_convergence_json(report, path)
    parsed = JSON3.read(read(path, String))

    @test written == path
    @test parsed.config.base.dt == 0.01
    @test collect(parsed.config.dt_values) == [0.01, 0.005]
    @test collect(parsed.config.D_values) == [1]
    @test collect(parsed.config.chi_values) == Int[]
    @test collect(parsed.config.cutoff_values) == [1e-12]
    @test haskey(parsed.summary, :max_abs_density_error_simple)
    @test parsed.summary.max_abs_density_error_ctm === nothing
    @test parsed.summary.all_ctm_trusted === nothing
    @test length(parsed.runs) == 2
    @test parsed.runs[1].config.dt == 0.01
end
