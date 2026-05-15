module Benchmarks

using JSON3

using ..IPEPSEvolution: EvolutionLog, TrotterParams, evolve!
using ..Observables: TFIMObservableSummary, measure_tfim_simple
using ..SquareIPEPS: product_square_ipeps
using ..SquareUnitCells: PeriodicSquareUnitCell
using ..StarModels: AbstractModelProtocol, StaticModel, TFIMStarModel, model_at

export BenchmarkSpec,
    BenchmarkMetadata,
    EvolutionDiagnostics,
    BenchmarkSample,
    BenchmarkResult,
    run_benchmark,
    write_benchmark_json,
    write_benchmark_csv

const BENCHMARK_INITIAL_STATES = (:z_up, :z_down, :x_plus)

"""
    BenchmarkSpec(name, protocol, cell, initial_state, total_time, trotter, measure_every)

Configuration for a v1 TFIM benchmark run. `protocol` must be a static
[`TFIMStarModel`](@ref) protocol, `initial_state` must be one of `:z_up`,
`:z_down`, or `:x_plus`, `total_time` must be a finite nonnegative integer
multiple of `trotter.dt`, and `measure_every` is the positive full-step
cadence for recorded samples.
"""
struct BenchmarkSpec{P<:AbstractModelProtocol,C}
    name::String
    protocol::P
    cell::C
    initial_state::Symbol
    total_time::Float64
    trotter::TrotterParams
    measure_every::Int

    function BenchmarkSpec(
        name::AbstractString,
        protocol::P,
        cell::C,
        initial_state::Symbol,
        total_time::Real,
        trotter::TrotterParams,
        measure_every::Integer,
    ) where {P<:AbstractModelProtocol,C<:PeriodicSquareUnitCell}
        isempty(name) && throw(ArgumentError("benchmark name must be nonempty"))
        initial_state in BENCHMARK_INITIAL_STATES ||
            throw(ArgumentError("initial_state must be :z_up, :z_down, or :x_plus"))
        total = Float64(total_time)
        isfinite(total) && total >= 0 ||
            throw(ArgumentError("total_time must be finite and nonnegative"))
        _is_integer_multiple(total, trotter.dt) ||
            throw(ArgumentError("total_time must be an integer multiple of trotter.dt"))
        cadence = Int(measure_every)
        cadence >= 1 || throw(ArgumentError("measure_every must be at least 1"))
        protocol isa StaticModel ||
            throw(ArgumentError("benchmark protocol must be StaticModel(TFIMStarModel)"))
        model = protocol.model
        model isa TFIMStarModel ||
            throw(ArgumentError("benchmark protocol must be StaticModel(TFIMStarModel)"))
        return new{P,C}(String(name), protocol, cell, initial_state, total, trotter, cadence)
    end
end

"""
    BenchmarkMetadata

Reproducibility metadata attached to a benchmark result, including TFIM model
parameters, package version, unit-cell dimensions, initial state, evolution
parameters, and the simple observable source tag.
"""
struct BenchmarkMetadata
    observable_source::String
    package_version::Union{Nothing,String}
    protocol_type::String
    model_type::String
    J::Union{Nothing,Float64}
    h::Union{Nothing,Float64}
    cell_Lx::Int
    cell_Ly::Int
    initial_state::Symbol
    total_time::Float64
    dt::Float64
    order::Int
    evolution::Symbol
    maxdim::Int
    cutoff::Float64
    split_order::NTuple{4,Symbol}
    measure_every::Int
end

"""
    EvolutionDiagnostics

Per-sample evolution diagnostics distilled from an [`EvolutionLog`](@ref):
truncation error, link entropy summaries, and log-normalization before/after
the measured full step.
"""
struct EvolutionDiagnostics
    max_truncerr::Float64
    mean_bond_entropy::Float64
    max_bond_entropy::Float64
    log_norm_before::Float64
    log_norm_after::Float64
    log_norm_delta::Float64
end

"""
    BenchmarkSample

A measured benchmark time point containing the full-step index, physical time,
TFIM simple observables, and evolution diagnostics for that sample.
"""
struct BenchmarkSample
    step::Int
    time::Float64
    observables::TFIMObservableSummary
    diagnostics::EvolutionDiagnostics
end

"""
    BenchmarkResult

Complete benchmark trajectory returned by [`run_benchmark`](@ref), including
run metadata, all recorded samples, and the final sample repeated as
`final_state_summary` for convenient table construction.
"""
struct BenchmarkResult
    name::String
    run_label::String
    metadata::BenchmarkMetadata
    samples::Vector{BenchmarkSample}
    final_state_summary::BenchmarkSample
end

function _is_integer_multiple(total::Float64, dt::Float64)
    total == 0.0 && return true
    nsteps = round(Int, total / dt)
    return isapprox(nsteps * dt, total; atol = 1e-12, rtol = 1e-10)
end

function _finite_float(x, name)
    value = Float64(x)
    isfinite(value) || throw(ArgumentError("$name must be finite"))
    return value
end

function _diagnostics_from_log(log::EvolutionLog)
    return EvolutionDiagnostics(
        _finite_float(log.max_truncerr, "max_truncerr"),
        _finite_float(log.mean_bond_entropy, "mean_bond_entropy"),
        _finite_float(log.max_bond_entropy, "max_bond_entropy"),
        _finite_float(log.log_norm_before, "log_norm_before"),
        _finite_float(log.log_norm_after, "log_norm_after"),
        _finite_float(log.log_norm_delta, "log_norm_delta"),
    )
end

_zero_diagnostics() = EvolutionDiagnostics(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

function _package_version()
    version = Base.pkgversion(@__MODULE__)
    return version === nothing ? nothing : string(version)
end

function _tfim_model(spec::BenchmarkSpec, time::Real, step::Integer)
    model = model_at(spec.protocol, time, step)
    model isa TFIMStarModel ||
        throw(ArgumentError("benchmark protocol must select a TFIMStarModel"))
    return model
end

function _metadata(spec::BenchmarkSpec)
    model = _tfim_model(spec, 0.0, 0)
    return BenchmarkMetadata(
        "simple",
        _package_version(),
        string(typeof(spec.protocol)),
        string(typeof(model)),
        _finite_float(model.J, "J"),
        _finite_float(model.h, "h"),
        spec.cell.Lx,
        spec.cell.Ly,
        spec.initial_state,
        spec.total_time,
        spec.trotter.dt,
        spec.trotter.order,
        spec.trotter.evolution,
        spec.trotter.maxdim,
        spec.trotter.cutoff,
        spec.trotter.split_order,
        spec.measure_every,
    )
end

"""
    run_benchmark(spec; run_label = "manual")

Run a TFIM benchmark from a product iPEPS initial state, returning samples at
step 0, every `spec.measure_every` full Trotter steps, and the final step.
Evolution is performed one full step at a time through [`evolve!`](@ref).
"""
function run_benchmark(spec::BenchmarkSpec; run_label::AbstractString = "manual")
    psi = product_square_ipeps(
        spec.cell;
        state = spec.initial_state,
        maxdim = spec.trotter.maxdim,
    )
    samples = BenchmarkSample[]
    model0 = _tfim_model(spec, 0.0, 0)
    push!(
        samples,
        BenchmarkSample(0, 0.0, measure_tfim_simple(psi, model0), _zero_diagnostics()),
    )

    nsteps = round(Int, spec.total_time / spec.trotter.dt)
    for step in 1:nsteps
        log = evolve!(psi, spec.trotter.dt; params = spec.trotter, protocol = spec.protocol)
        if step % spec.measure_every == 0 || step == nsteps
            time = step * spec.trotter.dt
            model = _tfim_model(spec, time, step)
            push!(
                samples,
                BenchmarkSample(
                    step,
                    time,
                    measure_tfim_simple(psi, model),
                    _diagnostics_from_log(log),
                ),
            )
        end
    end

    return BenchmarkResult(spec.name, String(run_label), _metadata(spec), samples, samples[end])
end

_json_scalar(x::Symbol) = String(x)
_json_scalar(x::Tuple) = map(_json_scalar, x)
_json_scalar(x) = x

function _summary_nt(obs::TFIMObservableSummary)
    return NamedTuple{fieldnames(TFIMObservableSummary)}(
        Tuple(_json_scalar(getfield(obs, f)) for f in fieldnames(TFIMObservableSummary)),
    )
end

function _diagnostics_nt(diag::EvolutionDiagnostics)
    return NamedTuple{fieldnames(EvolutionDiagnostics)}(
        Tuple(getfield(diag, f) for f in fieldnames(EvolutionDiagnostics)),
    )
end

function _metadata_nt(meta::BenchmarkMetadata)
    return NamedTuple{fieldnames(BenchmarkMetadata)}(
        Tuple(_json_scalar(getfield(meta, f)) for f in fieldnames(BenchmarkMetadata)),
    )
end

function _sample_nt(sample::BenchmarkSample)
    return (
        step = sample.step,
        time = sample.time,
        observables = _summary_nt(sample.observables),
        diagnostics = _diagnostics_nt(sample.diagnostics),
    )
end

function _result_nt(result::BenchmarkResult)
    return (
        name = result.name,
        run_label = result.run_label,
        metadata = _metadata_nt(result.metadata),
        samples = [_sample_nt(s) for s in result.samples],
        final_state_summary = _sample_nt(result.final_state_summary),
    )
end

function _assert_serializable_finite(x, path::String)
    return nothing
end

function _assert_serializable_finite(x::Real, path::String)
    isfinite(Float64(x)) || throw(ArgumentError("$path must be finite"))
    return nothing
end

function _assert_serializable_finite(x::Union{AbstractString,Symbol,Nothing}, path::String)
    return nothing
end

function _assert_serializable_finite(x::NamedTuple, path::String)
    for name in keys(x)
        _assert_serializable_finite(getfield(x, name), "$path.$name")
    end
    return nothing
end

function _assert_serializable_finite(xs::Tuple, path::String)
    for (idx, x) in pairs(xs)
        _assert_serializable_finite(x, "$path[$idx]")
    end
    return nothing
end

function _assert_serializable_finite(xs::AbstractVector, path::String)
    for (idx, x) in pairs(xs)
        _assert_serializable_finite(x, "$path[$idx]")
    end
    return nothing
end

"""
    write_benchmark_json(result, path)

Write a benchmark result as deterministic JSON3 output. Nonfinite numeric
fields are rejected before serialization.
"""
function write_benchmark_json(result::BenchmarkResult, path::AbstractString)
    data = _result_nt(result)
    _assert_serializable_finite(data, "result")
    open(path, "w") do io
        JSON3.write(io, data)
        println(io)
    end
    return path
end

const CSV_HEADER = [
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

function _csv_cell(x::Real)
    _assert_serializable_finite(x, "CSV value")
    return string(x)
end

_csv_cell(x::Symbol) = String(x)

function _csv_cell(x::AbstractString)
    if occursin(r"[,\n\"]", x)
        return "\"" * replace(x, "\"" => "\"\"") * "\""
    else
        return x
    end
end

_csv_cell(::Nothing) = ""

function _csv_row(result::BenchmarkResult, sample::BenchmarkSample)
    meta = result.metadata
    obs = sample.observables
    diag = sample.diagnostics
    values = (
        result.name,
        result.run_label,
        meta.observable_source,
        sample.step,
        sample.time,
        meta.J,
        meta.h,
        meta.maxdim,
        meta.dt,
        meta.order,
        obs.mean_x,
        obs.mean_y,
        obs.mean_z,
        obs.zz_right,
        obs.zz_up,
        obs.energy_density_star,
        obs.energy_density_decomposed,
        obs.energy_density_discrepancy,
        diag.max_truncerr,
        diag.mean_bond_entropy,
        diag.max_bond_entropy,
        diag.log_norm_delta,
    )
    return join(_csv_cell.(values), ",")
end

"""
    write_benchmark_csv(results, path)

Write one flattened deterministic CSV time-series row per benchmark sample.
The header starts with
`name,run_label,observable_source,step,time,J,h,D,dt,order`.
"""
function write_benchmark_csv(results, path::AbstractString)
    checked_results = collect(results)
    for result in checked_results
        _assert_serializable_finite(_result_nt(result), "result")
    end

    open(path, "w") do io
        println(io, join(CSV_HEADER, ","))
        for result in checked_results
            for sample in result.samples
                println(io, _csv_row(result, sample))
            end
        end
    end
    return path
end

end
