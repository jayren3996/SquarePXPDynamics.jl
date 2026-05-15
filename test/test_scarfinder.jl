@testset "ScarFinder parameter validation" begin
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)

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
    @test all(iteration -> iteration.simple_score.diagnostics === :simple, result.iterations)
    @test all(iteration -> iteration.ctm_score === nothing, result.iterations)
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
    @test occursin(",simple,", csv)
    @test occursin("\"iterations\"", json)
    @test occursin("\"diagnostics\":\"simple\"", json)
    @test length(csv_result.iterations) == 1
    @test length(json_result.iterations) == 1
end

@testset "ScarFinder source stays above star projection layer" begin
    source = read(joinpath(@__DIR__, "..", "src", "ScarFinder.jl"), String)

    @test occursin("evolve!", source)
    @test !occursin("project_star!", source)
end
