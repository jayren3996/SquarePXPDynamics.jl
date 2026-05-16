#!/usr/bin/env julia

using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using SquarePXPDynamics

function _env_float(name::String, default::Float64)
    value = get(ENV, name, "")
    isempty(value) && return default
    return parse(Float64, value)
end

function _env_int(name::String, default::Int)
    value = get(ENV, name, "")
    isempty(value) && return default
    return parse(Int, value)
end

function _env_bool(name::String, default::Bool)
    value = lowercase(get(ENV, name, ""))
    isempty(value) && return default
    value in ("1", "true", "yes", "y") && return true
    value in ("0", "false", "no", "n") && return false
    error("environment variable $name must be a boolean")
end

config = PXPEEDBenchmarkConfig(
    7;
    total_time = _env_float("PXP_ED_TOTAL_TIME", 0.1),
    dt = _env_float("PXP_ED_DT", 0.01),
    measure_every = _env_int("PXP_ED_MEASURE_EVERY", 1),
    initial_state = :down,
    point_group = _env_bool("PXP_ED_POINT_GROUP", true),
    use_sparse = _env_bool("PXP_ED_USE_SPARSE", true),
    tol = _env_float("PXP_ED_TOL", 1e-10),
    m_init = _env_int("PXP_ED_M_INIT", 30),
    m_max = _env_int("PXP_ED_M_MAX", 60),
    extend_step = _env_int("PXP_ED_EXTEND_STEP", 10),
)

output_path = get(ENV, "PXP_ED_OUTPUT", joinpath(@__DIR__, "pxp-ed-7x7.json"))

println("Running 7x7 PBC PXP ED benchmark")
println("  point group: ", config.point_group)
println("  sparse Hamiltonian: ", config.use_sparse)
println("  total_time: ", config.total_time)
println("  dt: ", config.dt)
println("  measure_every: ", config.measure_every)
println("  output: ", output_path)

result = run_pxp_ed_benchmark(config)
write_pxp_ed_benchmark_json(result, output_path)

println("Finished")
println("  basis dimension: ", result.basis_dimension)
println("  constrained dimension: ", result.constrained_dimension)
println("  group order: ", result.group_order)
println("  Hamiltonian nnz: ", something(result.hamiltonian_nnz, -1))
println("  Krylov matvecs: ", result.diagnostics.matvecs)
println("  final norm: ", result.samples[end].norm)
println("  final return probability: ", result.samples[end].return_probability)
println("  final excitation density: ", result.samples[end].excitation_density)
