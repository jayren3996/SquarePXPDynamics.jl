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

function _validation_csv_cell(csv, name)
    header = split(csv[1], ','; keepempty = true)
    row = split(csv[2], ','; keepempty = true)
    index = findfirst(==(name), header)
    index === nothing && error("missing CSV column $name")
    return row[index]
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

@testset "PXP validation can attach exact finite tiny-cell density" begin
    config = PXPValidationConfig(
        3;
        total_time = 0.02,
        dt = 0.02,
        measure_every = 1,
        maxdim = 2,
        cutoff = 1e-12,
        exact_finite_observables = true,
        exact_finite_max_sites = 12,
    )
    report = validate_pxp_ed_ipeps(config; ctm_params = nothing)

    @test report.config.exact_finite_observables === true
    @test all(sample -> sample.exact_finite_density !== nothing, report.ipeps_samples)
    @test all(comparison -> comparison.ipeps_exact_finite_density !== nothing, report.comparisons)
    @test all(comparison -> comparison.density_error_exact_finite !== nothing, report.comparisons)
    @test report.comparisons[end].ipeps_exact_finite_density ≈
          0.0003996269892620211 atol = 1e-12
    @test abs(report.comparisons[end].density_error_exact_finite) < 1e-6
    @test abs(report.comparisons[end].density_error_simple) > 1e-4
end

@testset "PXP validation rejects exact finite contraction above configured max sites" begin
    @test_throws ArgumentError PXPValidationConfig(
        3;
        total_time = 0.02,
        dt = 0.02,
        maxdim = 1,
        exact_finite_observables = true,
        exact_finite_max_sites = 8,
    )
    @test_throws ArgumentError PXPValidationConfig(
        4;
        total_time = 0.02,
        dt = 0.02,
        maxdim = 1,
        exact_finite_observables = true,
        exact_finite_max_sites = 12,
    )
end

@testset "PXP validation metadata uses package checkout git commit" begin
    config = PXPValidationConfig(3; total_time = 0.0, dt = 0.01, measure_every = 1)
    package_root = abspath(joinpath(@__DIR__, ".."))
    git_dir = joinpath(package_root, ".git")
    expected_commit = chomp(
        read(`git --git-dir $git_dir --work-tree $package_root rev-parse HEAD`, String),
    )

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
    @test haskey(parsed.config, :exact_finite_observables)
    @test parsed.config.exact_finite_observables === false
    @test parsed.config.exact_finite_max_sites == 12
    @test parsed.ipeps_samples[1].exact_finite_density === nothing
    @test parsed.comparisons[1].ipeps_exact_finite_density === nothing
    @test parsed.comparisons[1].density_error_exact_finite === nothing
end

@testset "PXP validation JSON preserves opt-in exact finite density" begin
    config = PXPValidationConfig(
        3;
        total_time = 0.02,
        dt = 0.02,
        measure_every = 1,
        maxdim = 2,
        cutoff = 1e-12,
        exact_finite_observables = true,
        exact_finite_max_sites = 9,
    )
    report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
    path = tempname() * ".json"

    write_pxp_validation_json(report, path)
    parsed = JSON3.read(read(path, String))

    @test parsed.config.exact_finite_observables === true
    @test parsed.config.exact_finite_max_sites == 9
    @test parsed.ipeps_samples[end].exact_finite_density ≈ 0.0003996269892620211 atol =
        1e-12
    @test abs(parsed.comparisons[end].density_error_exact_finite) < 1e-6
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

@testset "PXP convergence propagates exact finite density opt-in" begin
    base = PXPValidationConfig(
        3;
        total_time = 0.02,
        dt = 0.02,
        measure_every = 1,
        maxdim = 2,
        cutoff = 1e-12,
        exact_finite_observables = true,
        exact_finite_max_sites = 9,
    )
    sweep = PXPConvergenceConfig(
        base;
        dt_values = [0.02],
        D_values = [2],
        chi_values = Int[],
        cutoff_values = [1e-12],
    )

    report = validate_pxp_convergence(sweep)

    @test report.runs[1].config.exact_finite_observables === true
    @test report.runs[1].config.exact_finite_max_sites == 9
    @test report.runs[1].ipeps_samples[end].exact_finite_density !== nothing
    @test report.max_abs_density_error_exact_finite !== nothing
    @test report.max_abs_density_error_exact_finite < 1e-6
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
    @test haskey(parsed.summary, :max_abs_density_error_exact_finite)
    @test parsed.summary.max_abs_density_error_ctm === nothing
    @test parsed.summary.all_ctm_trusted === nothing
    @test length(parsed.runs) == 2
    @test parsed.runs[1].config.dt == 0.01
end

@testset "PXP audit campaign produces machine-readable summaries" begin
    config = PXPAuditConfig(;
        n_values = [3],
        total_time = 0.01,
        dt_values = [0.01],
        D_values = [1],
        cutoff_values = [1e-12],
        chi_values = Int[],
    )

    report = run_pxp_audit_campaign(config)

    @test report isa PXPAuditReport
    @test length(report.runs) == 1
    @test report.runs[1].validation isa PXPValidationReport
    @test report.runs[1].reversibility isa PXPReversibilityReport
    summary = report.runs[1].summary
    @test summary.n == 3
    @test summary.dt == 0.01
    @test summary.D == 1
    @test summary.cutoff == 1e-12
    @test summary.ctm_trust_status === :not_run
    @test summary.ctm_trust_reason === :not_run
    @test summary.max_abs_density_error_simple >= 0
    @test summary.max_abs_density_error_ctm === nothing
    @test summary.max_abs_density_error_exact_finite === nothing
    @test summary.max_blockade_violation_simple >= 0
    @test summary.pxp_energy_drift_simple >= 0
    @test summary.max_truncerr >= 0
    @test summary.log_norm_delta_abs >= 0
    @test summary.reversibility_density_drift >= 0
    @test summary.reversibility_blockade_drift >= 0
    @test summary.reversibility_energy_drift >= 0
end

@testset "PXP audit campaign can opt into exact finite density summaries" begin
    config = PXPAuditConfig(;
        n_values = [3],
        total_time = 0.02,
        dt_values = [0.02],
        D_values = [2],
        cutoff_values = [1e-12],
        chi_values = Int[],
        exact_finite_observables = true,
        exact_finite_max_sites = 9,
    )
    report = run_pxp_audit_campaign(config; ctm_measure = _validation_fake_ctm_summary)
    summary = report.runs[1].summary

    @test report.config.exact_finite_observables === true
    @test report.config.exact_finite_max_sites == 9
    @test report.runs[1].validation.config.exact_finite_max_sites == 9
    @test summary.max_abs_density_error_exact_finite !== nothing
    @test summary.max_abs_density_error_exact_finite < 1e-6
    @test summary.max_abs_density_error_simple > 1e-4
end

@testset "PXP audit campaign records CTM trust summaries with fake CTM" begin
    config = PXPAuditConfig(;
        n_values = [3],
        total_time = 0.0,
        dt_values = [0.01],
        D_values = [1],
        cutoff_values = [1e-12],
        chi_values = [2, 4],
    )

    report = run_pxp_audit_campaign(
        config;
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

    summary = only(report.runs).summary
    @test summary.max_abs_density_error_ctm !== nothing
    @test summary.max_blockade_violation_ctm !== nothing
    @test summary.ctm_trust_status === :trusted
    @test summary.ctm_trust_reason === :trusted
    @test summary.finite_chi_density_delta !== nothing
    @test summary.finite_chi_energy_delta !== nothing
end

@testset "PXP audit campaign writes JSON and CSV artifacts" begin
    config = PXPAuditConfig(;
        n_values = [3],
        total_time = 0.01,
        dt_values = [0.01],
        D_values = [1],
        cutoff_values = [1e-12],
        chi_values = Int[],
    )
    report = run_pxp_audit_campaign(config)
    json_path = tempname() * ".json"
    csv_path = tempname() * ".csv"

    written_json = write_pxp_audit_json(report, json_path)
    written_csv = write_pxp_audit_csv(report, csv_path)
    parsed = JSON3.read(read(json_path, String))
    csv = split(chomp(read(csv_path, String)), '\n')

    @test written_json == json_path
    @test written_csv == csv_path
    @test parsed.config.total_time == 0.01
    @test length(parsed.runs) == 1
    @test haskey(parsed.runs[1].summary, :max_abs_density_error_simple)
    @test haskey(parsed.runs[1].summary, :reversibility_energy_drift)
    @test startswith(csv[1], "n,total_time,dt,D,cutoff")
    @test occursin("max_abs_density_error_exact_finite", csv[1])
    @test _validation_csv_cell(csv, "max_abs_density_error_exact_finite") == ""
    @test length(csv) == 2
    @test occursin(",not_run,not_run,", csv[2])
end

@testset "PXP audit JSON and CSV preserve opt-in exact finite summaries" begin
    config = PXPAuditConfig(;
        n_values = [3],
        total_time = 0.02,
        dt_values = [0.02],
        D_values = [2],
        cutoff_values = [1e-12],
        chi_values = Int[],
        exact_finite_observables = true,
        exact_finite_max_sites = 9,
    )
    report = run_pxp_audit_campaign(config)
    json_path = tempname() * ".json"
    csv_path = tempname() * ".csv"

    write_pxp_audit_json(report, json_path)
    write_pxp_audit_csv(report, csv_path)
    parsed = JSON3.read(read(json_path, String))
    csv = split(chomp(read(csv_path, String)), '\n')

    @test parsed.config.exact_finite_observables === true
    @test parsed.config.exact_finite_max_sites == 9
    @test parsed.runs[1].validation.config.exact_finite_max_sites == 9
    @test parsed.runs[1].summary.max_abs_density_error_exact_finite < 1e-6
    @test occursin("max_abs_density_error_exact_finite", csv[1])
    @test parse(Float64, _validation_csv_cell(csv, "max_abs_density_error_exact_finite")) < 1e-6
    @test occursin("not_run,not_run", csv[2])
end
