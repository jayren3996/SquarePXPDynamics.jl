@testset "ScarFinder parameter validation" begin
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)

    @test_throws ArgumentError ScarFinderParams(-0.1, trotter, 1, Inf, Inf, Inf, false)
    @test_throws ArgumentError ScarFinderParams(0.1, trotter, -1, Inf, Inf, Inf, false)
    @test_throws ArgumentError ScarFinderParams(0.1, trotter, 1, -1.0, Inf, Inf, false)
    @test_throws ArgumentError ScarFinderParams(0.1, trotter, 1, Inf, -1.0, Inf, false)
    @test_throws ArgumentError ScarFinderParams(0.1, trotter, 1, Inf, Inf, -1.0, false)
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

    result = scarfinder!(
        psi;
        projection_time = 0.01,
        trotter = trotter,
        iterations = 2,
    )

    @test length(result.iterations) == 2
    @test result.accepted_iterations == 2
end

@testset "ScarFinder source stays above star projection layer" begin
    source = read(joinpath(@__DIR__, "..", "src", "ScarFinder.jl"), String)

    @test occursin("evolve!", source)
    @test !occursin("project_star!", source)
end
