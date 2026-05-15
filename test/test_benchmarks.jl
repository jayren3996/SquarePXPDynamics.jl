using Test
using JSON3
using SquarePXPDynamics

struct BenchmarkRampProtocol <: AbstractModelProtocol end

SquarePXPDynamics.model_at(::BenchmarkRampProtocol, time, step) =
    TFIMStarModel(0.0, 1.0 + time + step)

@testset "benchmark result records time series" begin
    spec = BenchmarkSpec(
        "tfim-j0",
        StaticModel(TFIMStarModel(0.0, 1.0)),
        PeriodicSquareUnitCell(10, 10),
        :z_up,
        0.02,
        TrotterParams(0.01, 1, :real, 1, 1e-12),
        1,
    )
    result = run_benchmark(spec; run_label = "unit-test")

    @test result.name == "tfim-j0"
    @test result.run_label == "unit-test"
    @test result.metadata.package_version == string(Base.pkgversion(SquarePXPDynamics))
    @test length(result.samples) == 3
    @test [s.step for s in result.samples] == [0, 1, 2]
    @test result.samples[end].time ≈ 0.02 atol = 1e-12
    @test result.final_state_summary === result.samples[end]
    @test result.samples[end].observables.mean_z ≈ cos(0.04) atol = 1e-6 rtol = 1e-6
end

@testset "benchmark JSON and CSV writers are deterministic" begin
    spec = BenchmarkSpec(
        "tfim-static",
        StaticModel(TFIMStarModel(1.0, 0.0)),
        PeriodicSquareUnitCell(10, 10),
        :z_up,
        0.01,
        TrotterParams(0.01, 1, :real, 1, 1e-12, (:up, :right, :down, :left)),
        1,
    )
    result = run_benchmark(spec; run_label = "serialize-test")
    dir = mktempdir()
    json_path = joinpath(dir, "result.json")
    csv_path = joinpath(dir, "result.csv")

    write_benchmark_json(result, json_path)
    write_benchmark_csv([result], csv_path)

    parsed = JSON3.read(read(json_path, String))
    @test parsed[:name] == "tfim-static"
    @test parsed[:run_label] == "serialize-test"
    @test parsed[:metadata][:observable_source] == "simple"
    @test parsed[:metadata][:package_version] == string(Base.pkgversion(SquarePXPDynamics))
    @test collect(parsed[:metadata][:split_order]) == ["up", "right", "down", "left"]
    @test length(parsed[:samples]) == 2

    csv = read(csv_path, String)
    expected_header = [
        "name",
        "run_label",
        "observable_source",
        "step",
        "time",
        "J",
        "h",
        "D",
        "dt",
        "order",
        "mean_x",
        "mean_y",
        "mean_z",
        "zz_right",
        "zz_up",
        "energy_star",
        "energy_decomposed",
        "energy_discrepancy",
        "max_truncerr",
        "mean_bond_entropy",
        "max_bond_entropy",
        "lognorm_delta",
    ]
    @test split(split(csv, '\n')[1], ',') == expected_header
    rows = split(chomp(csv), '\n')
    @test length(rows) == 3
    @test length.(split.(rows, ',')) == fill(length(expected_header), 3)
    @test split(rows[2], ',')[1:10] ==
          ["tfim-static", "serialize-test", "simple", "0", "0.0", "1.0", "0.0", "1", "0.01", "1"]
    @test split(rows[3], ',')[1:10] ==
          ["tfim-static", "serialize-test", "simple", "1", "0.01", "1.0", "0.0", "1", "0.01", "1"]
end

@testset "benchmark cadence records final sample" begin
    spec = BenchmarkSpec(
        "tfim-final",
        StaticModel(TFIMStarModel(0.0, 1.0)),
        PeriodicSquareUnitCell(10, 10),
        :z_up,
        0.03,
        TrotterParams(0.01, 1, :real, 1, 1e-12),
        2,
    )
    result = run_benchmark(spec)

    @test [s.step for s in result.samples] == [0, 2, 3]
    @test result.samples[end].time ≈ 0.03 atol = 1e-12
end

@testset "benchmark writers reject nonfinite records" begin
    spec = BenchmarkSpec(
        "tfim-finite",
        StaticModel(TFIMStarModel(1.0, 0.0)),
        PeriodicSquareUnitCell(10, 10),
        :z_up,
        0.0,
        TrotterParams(0.01, 1, :real, 1, 1e-12),
        1,
    )
    result = run_benchmark(spec)
    bad_diag = EvolutionDiagnostics(NaN, 0.0, 0.0, 0.0, 0.0, 0.0)
    bad_sample = BenchmarkSample(0, 0.0, result.samples[1].observables, bad_diag)
    obs = result.samples[1].observables
    bad_obs = TFIMObservableSummary(
        obs.mean_x,
        obs.mean_y,
        obs.mean_z,
        NaN,
        obs.z_odd,
        obs.zz_right,
        obs.zz_up,
        obs.energy_density_star,
        obs.energy_density_decomposed,
        obs.energy_density_discrepancy,
        obs.x_imag_abs,
        obs.y_imag_abs,
        obs.z_imag_abs,
        obs.zz_imag_abs,
        obs.energy_imag_abs,
        obs.max_imag_abs,
        obs.mean_bond_entropy,
        obs.max_bond_entropy,
    )
    hidden_bad_sample = BenchmarkSample(0, 0.0, bad_obs, result.samples[1].diagnostics)
    bad = BenchmarkResult(
        result.name,
        result.run_label,
        result.metadata,
        [bad_sample],
        bad_sample,
    )
    hidden_bad = BenchmarkResult(
        result.name,
        result.run_label,
        result.metadata,
        [hidden_bad_sample],
        hidden_bad_sample,
    )
    dir = mktempdir()

    @test_throws ArgumentError write_benchmark_json(bad, joinpath(dir, "bad.json"))
    @test_throws ArgumentError write_benchmark_csv([bad], joinpath(dir, "bad.csv"))
    @test_throws ArgumentError write_benchmark_json(hidden_bad, joinpath(dir, "hidden-bad.json"))
    @test_throws ArgumentError write_benchmark_csv([hidden_bad], joinpath(dir, "hidden-bad.csv"))
end

@testset "benchmark validation rejects unsupported inputs" begin
    @test_throws ArgumentError BenchmarkSpec(
        "bad-measure",
        StaticModel(TFIMStarModel(1.0, 0.0)),
        PeriodicSquareUnitCell(10, 10),
        :z_up,
        0.01,
        TrotterParams(0.01, 1, :real, 1, 1e-12),
        0,
    )
    @test_throws ArgumentError BenchmarkSpec(
        "bad-state",
        StaticModel(TFIMStarModel(1.0, 0.0)),
        PeriodicSquareUnitCell(10, 10),
        :bad,
        0.01,
        TrotterParams(0.01, 1, :real, 1, 1e-12),
        1,
    )
    @test_throws ArgumentError BenchmarkSpec(
        "bad-time",
        StaticModel(TFIMStarModel(1.0, 0.0)),
        PeriodicSquareUnitCell(10, 10),
        :z_up,
        0.015,
        TrotterParams(0.01, 1, :real, 1, 1e-12),
        1,
    )
    @test_throws ArgumentError BenchmarkSpec(
        "bad-time",
        StaticModel(TFIMStarModel(1.0, 0.0)),
        PeriodicSquareUnitCell(10, 10),
        :z_up,
        NaN,
        TrotterParams(0.01, 1, :real, 1, 1e-12),
        1,
    )
    @test_throws ArgumentError BenchmarkSpec(
        "dynamic-protocol",
        BenchmarkRampProtocol(),
        PeriodicSquareUnitCell(10, 10),
        :z_up,
        0.01,
        TrotterParams(0.01, 1, :real, 1, 1e-12),
        1,
    )
end
