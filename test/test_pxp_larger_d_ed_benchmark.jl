using Test
using JSON3
using SquarePXPDynamics

const RUN_EXTENDED_PXP_ED_TESTS =
    get(ENV, "SQUAREPXP_EXTENDED_PXP_ED_TESTS", "") == "1" ||
    get(ENV, "SQUAREPXP_EXTENDED_TESTS", "") == "1"

function _csv_cell(lines, name)
    header = split(lines[1], ','; keepempty = true)
    row = split(lines[2], ','; keepempty = true)
    index = findfirst(==(name), header)
    index === nothing && error("missing CSV column $name")
    return row[index]
end

@testset "PXP ED observable provenance rejects central-region claims" begin
    basis = pxp_ed_space_group_basis(5)

    @test pxp_ed_boundary_condition(basis) === :periodic
    @test pxp_ed_symmetry_sector(basis) === :fully_symmetric_space_group
    @test pxp_ed_observable_scope(basis) === :pbc_global_site_average
    @test pxp_ed_reference_label(basis) == "finite_pbc_global_density"
    @test pxp_ed_group_order(basis) == 8 * 5^2

    @test_throws ArgumentError pxp_ed_site_density_operator(basis, 13)
    @test_throws ArgumentError pxp_ed_region_density_operator(basis, 1:9)
end

@testset "exact finite return probability is available only through tiny-cell contraction" begin
    cell = PeriodicSquareUnitCell(3, 3)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    up = product_square_ipeps(cell; state = :up, maxdim = 1)

    @test exact_all_down_return_probability_finite(down; max_sites = 9) ≈ 1.0 atol = 1e-15
    @test exact_all_down_return_probability_finite(up; max_sites = 9) ≈ 0.0 atol = 1e-15

    large = product_square_ipeps(PeriodicSquareUnitCell(4, 4); state = :down, maxdim = 1)
    @test_throws ArgumentError exact_all_down_return_probability_finite(large; max_sites = 9)
end

@testset "larger-D PXP benchmark config validates controls" begin
    config = PXPLargerDBenchmarkConfig(;
        n_values = [3],
        total_time = 0.02,
        dt_values = [0.02],
        D_values = [1, 2, 3],
        cutoff_values = [1e-12],
        exact_finite_observables = true,
        exact_finite_max_sites = 9,
    )

    @test config.n_values == [3]
    @test config.D_values == [1, 2, 3]
    @test config.observable_mode === :auto
    @test config.ed_mode === :symmetric_pbc
    @test_throws ArgumentError PXPLargerDBenchmarkConfig(; n_values = Int[])
    @test_throws ArgumentError PXPLargerDBenchmarkConfig(; D_values = Int[])
    @test_throws ArgumentError PXPLargerDBenchmarkConfig(; observable_mode = :central_region)
    @test_throws ArgumentError PXPLargerDBenchmarkConfig(; ed_mode = :open_boundary)
end

@testset "larger-D PXP benchmark separates exact finite and simple diagnostics" begin
    config = PXPLargerDBenchmarkConfig(;
        n_values = [3],
        total_time = 0.02,
        dt_values = [0.02],
        D_values = [1, 2, 3],
        cutoff_values = [1e-12],
        exact_finite_observables = true,
        exact_finite_max_sites = 9,
    )

    report = run_pxp_larger_d_benchmark(config)

    @test length(report.runs) == 3
    @test all(run -> run.summary.observable_mode === :exact_finite, report.runs)
    @test all(run -> run.summary.ed_observable_scope === :pbc_global_site_average, report.runs)
    @test all(run -> run.summary.density_error_exact_finite !== nothing, report.runs)
    @test all(run -> run.summary.return_probability_error !== nothing, report.runs)
    @test all(run -> run.summary.density_error_simple !== nothing, report.runs)
    @test all(run -> run.summary.ed_runtime_seconds >= 0, report.runs)
    @test all(run -> run.summary.ipeps_runtime_seconds >= 0, report.runs)
    @test all(run -> run.summary.max_truncerr >= 0, report.runs)
    @test all(run -> run.summary.log_norm_delta_abs >= 0, report.runs)
    @test all(run -> run.summary.reversibility_density_drift >= 0, report.runs)

    d2 = only(run for run in report.runs if run.summary.D == 2)
    @test abs(d2.summary.density_error_exact_finite) < 1e-6
    @test abs(d2.summary.density_error_simple) > 1e-4
end

@testset "larger-D PXP benchmark JSON and CSV preserve required schema" begin
    config = PXPLargerDBenchmarkConfig(;
        n_values = [3],
        total_time = 0.02,
        dt_values = [0.02],
        D_values = [1, 2],
        cutoff_values = [1e-12],
        exact_finite_observables = true,
        exact_finite_max_sites = 9,
    )
    report = run_pxp_larger_d_benchmark(config)
    json_path = tempname() * ".json"
    csv_path = tempname() * ".csv"

    @test write_pxp_larger_d_benchmark_json(report, json_path) == json_path
    @test write_pxp_larger_d_benchmark_csv(report, csv_path) == csv_path

    parsed = JSON3.read(read(json_path, String))
    @test parsed.schema_version == 1
    @test parsed.config.ed_mode == "symmetric_pbc"
    @test parsed.config.observable_mode == "auto"
    @test length(parsed.runs) == 2
    @test parsed.runs[1].summary.ed_observable_scope == "pbc_global_site_average"
    @test parsed.runs[1].summary.observable_mode == "exact_finite"
    @test parsed.runs[2].summary.D == 2
    @test !any(k -> occursin("central", lowercase(String(k))), keys(parsed.runs[1].summary))

    csv = split(chomp(read(csv_path, String)), '\n')
    header = split(csv[1], ','; keepempty = true)
    required = [
        "n",
        "D",
        "dt",
        "cutoff",
        "total_time",
        "ed_basis_dimension",
        "ed_constrained_dimension",
        "ed_group_order",
        "ed_runtime_seconds",
        "ipeps_runtime_seconds",
        "observable_mode",
        "density_error_simple",
        "density_error_exact_finite",
        "density_error_ctm",
        "return_probability_error",
        "max_truncerr",
        "log_norm_initial",
        "log_norm_final",
        "log_norm_delta_abs",
        "reversibility_density_drift",
        "ctm_trust_status",
        "ctm_trust_reason",
        "notes",
        "warnings",
    ]
    @test all(name -> name in header, required)
    @test !any(name -> occursin("central", lowercase(name)), header)
    @test length(csv) == 3
    @test parse(Float64, _csv_cell(csv, "density_error_exact_finite")) < 1e-5
end

@testset "larger-D PXP benchmark script exists" begin
    script = joinpath(dirname(@__DIR__), "scripts", "pxp_larger_d_ed_benchmark.jl")
    @test isfile(script)
    text = read(script, String)
    @test occursin("SQUAREPXP_LARGERD_N", text)
    @test occursin("SQUAREPXP_LARGERD_EXACT_FINITE", text)
    @test occursin("write_pxp_larger_d_benchmark_json", text)
    @test occursin("write_pxp_larger_d_benchmark_csv", text)
end
