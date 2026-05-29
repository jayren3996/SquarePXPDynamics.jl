using Test
using JSON3
using SquarePXPDynamics

function _audit_config_simple_only(; cell_Lx = 3, cell_Ly = 3, iterations = 2)
    return ScarFinderAuditConfig(;
        cell_Lx,
        cell_Ly,
        initial_state = :down,
        projection_time = 0.02,
        iterations,
        dt_values = [0.02, 0.01],
        target_maxdim_values = [2], evolve_maxdim = 2,
        cutoff_values = [1e-12],
    )
end

@testset "ScarFinderAuditConfig validates inputs" begin
    base = _audit_config_simple_only()
    @test base.iterations == 2
    @test length(base.dt_values) == 2

    @test_throws ArgumentError ScarFinderAuditConfig(;
        cell_Lx = 1,
        cell_Ly = 3,
        projection_time = 0.02,
        iterations = 1,
        dt_values = [0.02],
        target_maxdim_values = [1], evolve_maxdim = 1,
    )
    @test_throws ArgumentError ScarFinderAuditConfig(;
        cell_Lx = 3,
        cell_Ly = 3,
        projection_time = 0.02,
        iterations = 0,
        dt_values = [0.02],
        target_maxdim_values = [1], evolve_maxdim = 1,
    )
    @test_throws ArgumentError ScarFinderAuditConfig(;
        cell_Lx = 3,
        cell_Ly = 3,
        projection_time = 0.02,
        iterations = 1,
        dt_values = Float64[],
        target_maxdim_values = [1], evolve_maxdim = 1,
    )
    @test_throws ArgumentError ScarFinderAuditConfig(;
        cell_Lx = 3,
        cell_Ly = 3,
        projection_time = 0.02,
        iterations = 1,
        dt_values = [-0.01],
        target_maxdim_values = [1], evolve_maxdim = 1,
    )
    @test_throws ArgumentError ScarFinderAuditConfig(;
        cell_Lx = 3,
        cell_Ly = 3,
        projection_time = 0.02,
        iterations = 1,
        dt_values = [0.02],
        target_maxdim_values = [1], evolve_maxdim = 1,
        require_trusted_ctm = true,
    )
end

@testset "run_scarfinder_audit produces simple-backend rows with stability" begin
    config = _audit_config_simple_only(; iterations = 2)
    report = run_scarfinder_audit(config)

    @test length(report.rows) == 2
    @test all(row -> row.ok, report.rows)
    @test all(row -> row.backend === :simple, report.rows)
    @test all(row -> row.chi === nothing, report.rows)
    @test all(row -> row.iterations_run == config.iterations, report.rows)
    @test all(row -> length(row.per_iteration_scores) == config.iterations, report.rows)
    @test all(row -> row.best_iteration !== nothing, report.rows)
    @test report.stability.n_rows == 2
    @test report.stability.n_successful == 2
    @test report.stability.best_iteration_agreement_all !== nothing
    @test report.stability.top_k_jaccard_min_all !== nothing
    @test report.stability.top_k_jaccard_min_simple !== nothing
    @test report.stability.n_simple == length(report.rows)
    @test report.stability.n_ctm == 0
    @test all(row -> row.max_truncerr_used !== nothing, report.rows)
    @test all(row -> row.reversibility_density_drift !== nothing, report.rows)
    @test all(row -> row.reversibility_density_drift >= 0, report.rows)
    @test all(row -> row.reversibility_blockade_drift !== nothing, report.rows)
    @test all(row -> row.reversibility_energy_drift !== nothing, report.rows)
end

@testset "projection_time is per-iteration, not aggregate" begin
    # Regression: ScarFinderAuditConfig.projection_time must be the per-iteration
    # evolution time, matching ScarFinderParams's first positional argument.
    # Treating it as total time across iterations historically caused over-evolution
    # and non-Hermitian CTM observables.
    config = _audit_config_simple_only(; iterations = 3)
    report = run_scarfinder_audit(config)
    successful = filter(row -> row.ok, report.rows)
    @test !isempty(successful)
    row = successful[1]
    # iter index runs 1..iterations and at least one iteration should produce
    # a finite per-iteration score, proving the loop ran the intended count.
    finite_count = count(s -> isfinite(s), row.per_iteration_scores)
    @test finite_count == config.iterations
end

@testset "ScarFinderAuditReport CSV and JSON round-trip cleanly" begin
    config = _audit_config_simple_only(; iterations = 2)
    report = run_scarfinder_audit(config)

    csv_path = tempname() * ".csv"
    json_path = tempname() * ".json"
    try
        write_scarfinder_audit_csv(report, csv_path)
        write_scarfinder_audit_json(report, json_path)

        csv_lines = readlines(csv_path)
        @test length(csv_lines) == length(report.rows) + 1
        @test occursin("dt,target_maxdim,evolve_maxdim,cutoff,chi,backend", csv_lines[1])

        loaded = JSON3.read(read(json_path, String))
        @test loaded.config.cell_Lx == config.cell_Lx
        @test length(loaded.rows) == length(report.rows)
        @test loaded.stability.n_successful == report.stability.n_successful
    finally
        rm(csv_path; force = true)
        rm(json_path; force = true)
    end
end
