@testset "ScarFinder parameter validation" begin
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    params = ScarFinderParams(0.01, trotter, 1, Inf, Inf, Inf, false)
    @test params.trotter == TrotterParams(0.01, 1, :real, 1, 1e-12)
    @test params.target_energy === nothing
    @test params.correction_time == 0.0
    @test params.correction_attempts == 0

    corrected = ScarFinderParams(
        0.01,
        trotter,
        1,
        Inf,
        Inf,
        Inf,
        false;
        target_energy = -0.1,
        correction_time = 0.01,
        correction_attempts = 2,
    )
    @test corrected.target_energy == -0.1
    @test corrected.correction_time == 0.01
    @test corrected.correction_attempts == 2

    @test_throws ArgumentError ScarFinderParams(-0.1, trotter, 1, Inf, Inf, Inf, false)
    @test_throws ArgumentError ScarFinderParams(Inf, trotter, 1, Inf, Inf, Inf, false)
    @test_throws ArgumentError ScarFinderParams(NaN, trotter, 1, Inf, Inf, Inf, false)
    @test_throws ArgumentError ScarFinderParams(0.1, trotter, -1, Inf, Inf, Inf, false)
    @test_throws ArgumentError ScarFinderParams(0.1, trotter, 1, -1.0, Inf, Inf, false)
    @test_throws ArgumentError ScarFinderParams(0.1, trotter, 1, Inf, -1.0, Inf, false)
    @test_throws ArgumentError ScarFinderParams(0.1, trotter, 1, Inf, Inf, -1.0, false)
    @test_throws ArgumentError ScarFinderParams(0.01, trotter, 1, NaN, Inf, Inf, false)
    @test_throws ArgumentError ScarFinderParams(0.01, trotter, 1, Inf, NaN, Inf, false)
    @test_throws ArgumentError ScarFinderParams(0.01, trotter, 1, Inf, Inf, NaN, false)
    @test_throws ArgumentError ScarFinderParams(
        0.01,
        trotter,
        1,
        Inf,
        Inf,
        Inf,
        false;
        target_energy = NaN,
    )
    @test_throws ArgumentError ScarFinderParams(
        0.01,
        trotter,
        1,
        Inf,
        Inf,
        Inf,
        false;
        target_energy = 0.0,
        correction_time = -0.01,
    )
    @test_throws ArgumentError ScarFinderParams(
        0.01,
        trotter,
        1,
        Inf,
        Inf,
        Inf,
        false;
        target_energy = 0.0,
        correction_attempts = -1,
    )
    @test_throws ArgumentError ScarFinderParams(0.01, "bad", 1, Inf, Inf, Inf, false)
end

@testset "ScarFinder zero iterations do not mutate state" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    weights_before = deepcopy(psi.link_weights)
    obs_before = measure_simple(psi)
    params = ScarFinderParams(0.01, trotter, 0, Inf, Inf, Inf, false)

    result = scarfinder!(psi, params)

    @test result.state === psi
    @test isempty(result.iterations)
    @test result.accepted_iterations == 0
    @test result.rejected_iterations == 0
    @test psi.link_weights == weights_before
    @test measure_simple(psi).density ≈ obs_before.density
end

@testset "ScarFinder one accepted iteration" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = ScarFinderParams(0.01, trotter, 1, Inf, Inf, Inf, false)

    result = scarfinder!(psi, params)

    @test length(result.iterations) == 1
    @test result.accepted_iterations == 1
    @test result.rejected_iterations == 0
    @test result.iterations[1].accepted
    @test result.iterations[1].reject_reason === nothing
    @test result.iterations[1].evolution isa EvolutionLog
    @test result.iterations[1].observables isa SimpleObservableSummary
    @test result.iterations[1].correction_accepted === nothing
    @test result.iterations[1].correction_energy_before === nothing
    @test result.iterations[1].correction_energy_after === nothing
end

@testset "ScarFinder rejects non-improving simple energy correction" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = ScarFinderParams(
        0.0,
        trotter,
        1,
        Inf,
        Inf,
        Inf,
        false;
        target_energy = measure_simple(psi).pxp_energy_density,
        correction_time = 0.01,
        correction_attempts = 1,
    )

    result = scarfinder!(psi, params)
    iteration = only(result.iterations)

    @test iteration.correction_accepted === false
    @test iteration.correction_energy_before ≈ 0.0 atol = 1e-12
    @test iteration.correction_energy_after ≈ iteration.correction_energy_before atol = 1e-12
    @test iteration.simple_score.correction_accepted === false
end

@testset "ScarFinder simple energy correction never worsens recorded objective" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    target = -0.1
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = ScarFinderParams(
        0.0,
        trotter,
        1,
        Inf,
        Inf,
        Inf,
        false;
        target_energy = target,
        correction_time = 0.01,
        correction_attempts = 2,
    )

    result = scarfinder!(psi, params)
    iteration = only(result.iterations)

    @test iteration.correction_accepted isa Bool
    @test isfinite(iteration.correction_energy_before)
    @test isfinite(iteration.correction_energy_after)
    @test abs(iteration.correction_energy_after - target) <=
          abs(iteration.correction_energy_before - target) + 1e-12
    @test iteration.simple_score.correction_energy_after == iteration.correction_energy_after
end

@testset "ScarFinder multiple accepted iterations" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = ScarFinderParams(0.01, trotter, 3, Inf, Inf, Inf, false)

    result = scarfinder!(psi, params)

    @test length(result.iterations) == 3
    @test result.accepted_iterations == 3
    @test all(it.iteration == i for (i, it) in enumerate(result.iterations))
    @test all(it.accepted for it in result.iterations)
    @test all(isfinite(it.observables.density) for it in result.iterations)
end

@testset "ScarFinder deterministic rejection by blockade threshold" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :up, maxdim = 1)
    params = ScarFinderParams(0.0, trotter, 1, Inf, 0.5, Inf, false)

    result = scarfinder!(psi, params)

    @test length(result.iterations) == 1
    @test result.accepted_iterations == 0
    @test result.rejected_iterations == 1
    @test !result.iterations[1].accepted
    @test occursin("blockade", result.iterations[1].reject_reason)
end

@testset "ScarFinder stops on rejection" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :up, maxdim = 1)
    params = ScarFinderParams(0.0, trotter, 5, Inf, 0.5, Inf, true)

    result = scarfinder!(psi, params)

    @test length(result.iterations) == 1
    @test result.accepted_iterations == 0
    @test result.rejected_iterations == 1
end

@testset "ScarFinder continues after rejection when requested" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :up, maxdim = 1)
    params = ScarFinderParams(0.0, trotter, 3, Inf, 0.5, Inf, false)

    result = scarfinder!(psi, params)

    @test length(result.iterations) == 3
    @test result.accepted_iterations == 0
    @test result.rejected_iterations == 3
    @test all(!it.accepted for it in result.iterations)
end

@testset "ScarFinder keyword convenience path" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)

    result = scarfinder!(psi; projection_time = 0.01, trotter = trotter, iterations = 2)

    @test length(result.iterations) == 2
    @test result.accepted_iterations == 2
end

@testset "ScarFinder preserves legacy unprojected TrotterParams" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, false, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :up, maxdim = 1)
    params = ScarFinderParams(0.01, trotter, 1, Inf, Inf, Inf, false)

    result = scarfinder!(psi, params)

    @test params.trotter == TrotterParams(0.01, 1, :real, 1, 1e-12)
    @test length(result.iterations) == 1
    @test result.iterations[1].evolution.params == params.trotter
    @test isfinite(result.iterations[1].evolution.max_truncerr)
end

@testset "ScarFinder simple candidate scores rank without CTM" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = ScarFinderParams(0.01, trotter, 2, Inf, Inf, Inf, false)

    result = scarfinder!(psi, params)
    ranked = rank_scarfinder_candidates(result; diagnostics = :simple)

    @test length(ranked) == 2
    @test all(score -> score isa ScarFinderCandidateScore, ranked)
    @test all(score -> score.diagnostics === :simple, ranked)
    @test all(score -> isfinite(score.score), ranked)
    @test all(score -> isfinite(score.log_norm_after), ranked)
    @test all(score -> score.log_norm_delta ≈ score.log_norm_after - score.log_norm_before, ranked)
    @test all(iteration -> iteration.simple_score.diagnostics === :simple, result.iterations)
    @test all(iteration -> iteration.ctm_score === nothing, result.iterations)
end

@testset "ScarFinder physics objectives score candidates" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = ScarFinderParams(0.0, trotter, 1, Inf, Inf, Inf, false)
    objective = CompositeObjective(;
        revival = RevivalObjective(:sublattice_imbalance, 1.0),
        target_energy = TargetEnergyObjective(-0.25, 4.0),
        low_variance = LowVarianceObjective(0.75),
        blockade_weight = 10.0,
        truncation_weight = 2.0,
        finite_chi_weight = 3.0,
        entropy_weight = 0.5,
    )

    result = scarfinder!(psi, params; objective)
    score = only(rank_scarfinder_candidates(result; diagnostics = :simple))

    @test score.objective_name == "CompositeObjective"
    @test isfinite(score.score)
    @test score.revival_strength !== nothing
    @test score.finite_chi_drift === nothing
    @test score.energy_variance_proxy === nothing
    @test score.objective_parameters ==
          "revival_observable=sublattice_imbalance;revival_weight=1.0;target_energy=-0.25;target_energy_weight=4.0;low_variance_weight=0.75;blockade_weight=10.0;truncation_weight=2.0;finite_chi_weight=3.0;entropy_weight=0.5"

    expected_score =
        score.blockade_violation * 10.0 +
        score.max_truncerr * 2.0 +
        score.max_bond_entropy * 0.5 -
        score.revival_strength * 1.0 +
        abs(score.pxp_energy_density - (-0.25)) * 4.0
    @test score.score ≈ expected_score atol = 1e-12

    csv_path = tempname() * ".csv"
    json_path = tempname() * ".json"
    write_scarfinder_log(result, csv_path; format = :csv)
    write_scarfinder_log(result, json_path; format = :json)
    csv = read(csv_path, String)
    json = read(json_path, String)
    @test occursin("objective_parameters", csv)
    @test occursin("revival_observable=sublattice_imbalance", csv)
    @test occursin("target_energy=-0.25", csv)
    @test occursin("\"objective_parameters\"", json)
    @test occursin("low_variance_weight=0.75", json)
end

@testset "ScarFinder CTM callback is optional and scheduled" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = ScarFinderParams(0.0, trotter, 3, Inf, Inf, Inf, false)
    calls = Int[]

    result = scarfinder!(
        psi,
        params;
        ctm_every = 2,
        ctm_at_end = true,
        ctm_callback = (state, iteration, simple_score) -> begin
            push!(calls, iteration)
            return measure_simple(state)
        end,
    )
    ranked_ctm = rank_scarfinder_candidates(result; diagnostics = :ctm)

    @test calls == [2, 3]
    @test result.iterations[1].ctm_score === nothing
    @test result.iterations[2].ctm_score isa ScarFinderCandidateScore
    @test result.iterations[3].ctm_score isa ScarFinderCandidateScore
    @test length(ranked_ctm) == 2
    @test all(score -> score.diagnostics === :ctm, ranked_ctm)
    @test result.iterations[1].simple_score.ctm_accepted === nothing
end

@testset "ScarFinder measurement backends construct and validate" begin
    simple = SimpleBackend()
    @test simple isa MeasurementBackend

    params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0),
    )
    policy = CTMTrustPolicy(2, true, 1e-2, 1e-3, 1e-2, 1e-4)
    trusted = TrustedCTMBackend(params, policy)

    @test trusted isa MeasurementBackend
    @test trusted.params === params
    @test trusted.policy === policy
    @test_throws ArgumentError TrustedCTMBackend((), policy)
end

@testset "ScarFinder measurement backends measure and integrate with CTM callbacks" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)

    simple_measurement = measure_scarfinder(psi, SimpleBackend())
    @test simple_measurement isa SimpleObservableSummary
    @test simple_measurement == measure_simple(psi)

    params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0),
    )
    policy = CTMTrustPolicy(2, true, 1e-2, 1e-3, 1e-2, 1e-4)
    calls = PEPSKitCTMRGParams[]
    fake_measure = (state; params) -> begin
        push!(calls, params)
        return CTMObservableSummary(
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            CTMRGDiagnostics(
                params.chi,
                params.tol,
                params.maxiter,
                params.maxiter,
                params.tol / 10,
                true,
                true,
            ),
        )
    end
    backend = TrustedCTMBackend(params, policy; measure = fake_measure)

    trusted_measurement = measure_scarfinder(psi, backend)
    @test trusted_measurement isa TrustedCTMMeasurement
    @test calls == collect(params)
    @test trusted_measurement.measurement isa CTMObservableSummary
    @test trusted_measurement.measurement.diagnostics.chi == 4
    @test trusted_measurement.trust.trusted

    empty!(calls)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    scar_params = ScarFinderParams(0.0, trotter, 1, Inf, Inf, Inf, false)
    result = scarfinder!(
        psi,
        scar_params;
        ctm_every = 1,
        ctm_callback = (state, iteration, simple_score) -> measure_scarfinder(state, backend),
    )
    ctm_score = only(rank_scarfinder_candidates(result; diagnostics = :ctm))

    @test calls == collect(params)
    @test ctm_score.ctm_chi == 4
    @test ctm_score.ctm_residual == 1e-7
    @test ctm_score.ctm_accepted === true
end

@testset "ScarFinder trusted CTM backend gates ranking" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    params = ScarFinderParams(0.0, trotter, 1, Inf, Inf, Inf, false)
    ctm_params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0),
    )
    backend = TrustedCTMBackend(
        ctm_params,
        CTMTrustPolicy(2, true, 1e-2, 1e-3, 1e-2, 1e-4);
        measure = (state; params) -> CTMObservableSummary(
            0.1 + params.chi * 1e-5,
            0.12,
            0.08,
            params.chi * 1e-6,
            -0.01,
            CTMRGDiagnostics(params.chi, params.tol, params.maxiter, 3, params.tol / 10, true, true),
        ),
    )

    result = scarfinder!(
        product_square_ipeps(cell; state = :down, maxdim = 1),
        params;
        measurement = backend,
        ctm_every = 1,
        require_trusted_ctm = true,
    )

    ranked = rank_scarfinder_candidates(result; diagnostics = :ctm, require_ctm_trusted = true)
    @test length(ranked) == 1
    @test ranked[1].ctm_trusted === true
    @test ranked[1].ctm_trust_reason == "trusted"
    @test ranked[1].finite_chi_drift !== nothing

    csv_path = tempname() * ".csv"
    json_path = tempname() * ".json"
    write_scarfinder_log(result, csv_path; format = :csv)
    write_scarfinder_log(result, json_path; format = :json)
    csv = read(csv_path, String)
    json = read(json_path, String)
    @test occursin("ctm_trusted", csv)
    @test occursin("ctm_trust_reason", csv)
    @test occursin(",ctm,", csv)
    @test occursin(",true,trusted,", csv)
    @test occursin("\"ctm_trusted\":true", json)
    @test occursin("\"ctm_trust_reason\":\"trusted\"", json)
end

@testset "ScarFinder require trusted CTM validates scheduled observations" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    params = ScarFinderParams(0.0, trotter, 1, Inf, Inf, Inf, false)

    @test_throws ArgumentError scarfinder!(
        product_square_ipeps(cell; state = :down, maxdim = 1),
        params;
        require_trusted_ctm = true,
    )

    @test_throws ArgumentError scarfinder!(
        product_square_ipeps(cell; state = :down, maxdim = 1),
        params;
        ctm_every = 1,
        require_trusted_ctm = true,
        ctm_callback = (state, iteration, simple_score) -> nothing,
    )
end

@testset "ScarFinder untrusted CTM policy rejects iteration" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    params = ScarFinderParams(0.0, trotter, 1, Inf, Inf, Inf, false)
    ctm_params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0),
    )
    backend = TrustedCTMBackend(
        ctm_params,
        CTMTrustPolicy(2, true, 1e-6, 1e-6, 1e-6, 1e-4);
        measure = (state; params) -> CTMObservableSummary(
            0.1 + params.chi * 1e-3,
            0.12,
            0.08,
            params.chi * 1e-6,
            -0.01,
            CTMRGDiagnostics(params.chi, params.tol, params.maxiter, 3, params.tol / 10, true, true),
        ),
    )

    result = scarfinder!(
        product_square_ipeps(cell; state = :down, maxdim = 1),
        params;
        measurement = backend,
        ctm_every = 1,
        require_trusted_ctm = true,
    )
    ctm_score = only(rank_scarfinder_candidates(result; diagnostics = :ctm))

    @test result.accepted_iterations == 0
    @test result.rejected_iterations == 1
    @test result.iterations[1].reject_reason == "trusted CTM policy rejected iteration"
    @test ctm_score.accepted === false
    @test ctm_score.ctm_trusted === false
    @test ctm_score.ctm_trust_reason == "density_delta_too_large"
    @test isempty(rank_scarfinder_candidates(result; diagnostics = :ctm, require_ctm_trusted = true))
end

@testset "ScarFinder CTM diagnostics are logged and can require trust" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = ScarFinderParams(0.0, trotter, 1, Inf, Inf, Inf, false)
    diag = CTMRGDiagnostics(4, 1e-8, 10, 10, 1e-9, true, true)
    ctm_summary = CTMObservableSummary(0.0, 0.0, 0.0, 0.0, 0.0, diag)

    result = scarfinder!(
        psi,
        params;
        ctm_every = 1,
        ctm_callback = (state, iteration, simple_score) -> ctm_summary,
    )
    ctm_score = only(rank_scarfinder_candidates(result; diagnostics = :ctm))

    @test ctm_score.ctm_chi == 4
    @test ctm_score.ctm_tol == 1e-8
    @test ctm_score.ctm_maxiter == 10
    @test ctm_score.ctm_iterations == 10
    @test ctm_score.ctm_residual == 1e-9
    @test ctm_score.ctm_converged === true
    @test ctm_score.ctm_accepted === true
    @test isfinite(ctm_score.log_norm_after)
    @test only(rank_scarfinder_candidates(result; diagnostics = :ctm, require_ctm_accepted = true)) ===
          ctm_score

    csv_path = tempname() * ".csv"
    json_path = tempname() * ".json"
    write_scarfinder_log(result, csv_path; format = :csv)
    write_scarfinder_log(result, json_path; format = :json)
    csv = read(csv_path, String)
    json = read(json_path, String)
    @test occursin("ctm_chi", csv)
    @test occursin("log_norm_delta", csv)
    @test occursin("correction_accepted", csv)
    ctm_row = split(split(chomp(csv), '\n')[3], ',')
    @test ctm_row[4] == "ctm"
    @test ctm_row[22:28] == ["4", "1.0e-8", "10", "10", "1.0e-9", "true", "true"]
    @test occursin("\"log_norm_delta\"", json)
    @test occursin("\"ctm_accepted\":true", json)
    @test occursin("\"correction_accepted\"", json)

    untrusted = CTMObservableSummary(
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        CTMRGDiagnostics(4, 1e-8, 10, 10, 1e-2, false, false),
    )
    result_untrusted = scarfinder!(
        product_square_ipeps(cell; state = :down, maxdim = 1),
        params;
        ctm_every = 1,
        ctm_callback = (state, iteration, simple_score) -> untrusted,
    )
    @test isempty(
        rank_scarfinder_candidates(
            result_untrusted;
            diagnostics = :ctm,
            require_ctm_accepted = true,
        ),
    )
end

@testset "ScarFinder nullable CTM sort keys put missing values last" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = ScarFinderParams(0.0, trotter, 2, Inf, Inf, Inf, false)
    summaries = Dict(
        1 => CTMObservableSummary(
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            CTMRGDiagnostics(4, 1e-8, 10, 10, nothing, true, true),
        ),
        2 => CTMObservableSummary(
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            CTMRGDiagnostics(4, 1e-8, 10, 7, 1e-9, true, true),
        ),
    )

    result = scarfinder!(
        psi,
        params;
        ctm_every = 1,
        ctm_callback = (state, iteration, simple_score) -> summaries[iteration],
    )

    ranked_residual = rank_scarfinder_candidates(result; diagnostics = :ctm, by = :ctm_residual)
    @test [score.iteration for score in ranked_residual] == [2, 1]

    ranked_iterations = rank_scarfinder_candidates(result; diagnostics = :ctm, by = :ctm_iterations, rev = true)
    @test [score.iteration for score in ranked_iterations] == [1, 2]
end

@testset "ScarFinder iteration logs write CSV and JSON" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    params = ScarFinderParams(0.0, trotter, 1, Inf, Inf, Inf, false)

    csv_path = tempname() * ".csv"
    json_path = tempname() * ".json"

    csv_result = scarfinder!(
        product_square_ipeps(cell; state = :down, maxdim = 1),
        params;
        log_path = csv_path,
        log_format = :csv,
    )
    json_result = scarfinder!(
        product_square_ipeps(cell; state = :down, maxdim = 1),
        params;
        log_path = json_path,
        log_format = :json,
    )

    csv = read(csv_path, String)
    json = read(json_path, String)
    @test occursin("iteration,accepted,reject_reason", csv)
    @test occursin("log_norm_before,log_norm_after,log_norm_delta", csv)
    @test occursin("correction_energy_after", csv)
    @test occursin(",simple,", csv)
    @test occursin("\"iterations\"", json)
    @test occursin("\"diagnostics\":\"simple\"", json)
    @test occursin("\"log_norm_after\"", json)
    @test occursin("\"correction_energy_after\"", json)
    @test length(csv_result.iterations) == 1
    @test length(json_result.iterations) == 1
end

@testset "ScarFinder source stays above star projection layer" begin
    source = read(joinpath(@__DIR__, "..", "src", "ScarFinder.jl"), String)

    @test occursin("evolve!", source)
    @test !occursin("project_star!", source)
end
