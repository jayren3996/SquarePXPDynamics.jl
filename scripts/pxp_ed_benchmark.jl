#!/usr/bin/env julia

using Pkg

project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using SquarePXPDynamics

function _env_value(name::String, default::AbstractString)
    value = get(ENV, name, "")
    return isempty(value) ? String(default) : value
end

function _env_bool(name::String, default::Bool)
    value = lowercase(strip(_env_value(name, string(default))))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    throw(ArgumentError("$name must be one of 1,true,yes,on,0,false,no,off"))
end

function _env_int(name::String, default::Int)
    return parse(Int, _env_value(name, string(default)))
end

function _env_float(name::String, default::Float64)
    return parse(Float64, _env_value(name, string(default)))
end

function _env_symbol(name::String, default::Symbol)
    return Symbol(_env_value(name, String(default)))
end

function _csv_value(value)
    return string(value)
end

function _write_csv(path::AbstractString, result)
    fields = (:step, :time, :norm, :return_probability, :excitation_density)
    open(path, "w") do io
        println(io, join(String.(fields), ","))
        for sample in result.samples
            println(io, join((_csv_value(getfield(sample, field)) for field in fields), ","))
        end
    end
    return path
end

n = _env_int("SQUAREPXP_ED_N", 3)
total_time = _env_float("SQUAREPXP_ED_TOTAL_TIME", 2.0)
dt = _env_float("SQUAREPXP_ED_DT", 0.02)
measure_every = _env_int("SQUAREPXP_ED_MEASURE_EVERY", 1)
initial_state = _env_symbol("SQUAREPXP_ED_INITIAL_STATE", :down)
point_group = _env_bool("SQUAREPXP_ED_POINT_GROUP", true)
use_sparse = _env_bool("SQUAREPXP_ED_USE_SPARSE", true)
tol = _env_float("SQUAREPXP_ED_TOL", 1e-10)
m_init = _env_int("SQUAREPXP_ED_M_INIT", 30)
m_max = _env_int("SQUAREPXP_ED_M_MAX", 60)
extend_step = _env_int("SQUAREPXP_ED_EXTEND_STEP", 10)

json_out = _env_value(
    "SQUAREPXP_ED_JSON",
    joinpath(project_root, "artifacts", "pxp_ed_benchmark.json"),
)
csv_out = _env_value(
    "SQUAREPXP_ED_CSV",
    joinpath(project_root, "artifacts", "pxp_ed_benchmark.csv"),
)
mkpath(dirname(json_out))
mkpath(dirname(csv_out))

config = PXPEEDBenchmarkConfig(
    n;
    total_time,
    dt,
    measure_every,
    initial_state,
    point_group,
    use_sparse,
    tol,
    m_init,
    m_max,
    extend_step,
)

result = run_pxp_ed_benchmark(config)
write_pxp_ed_benchmark_json(result, json_out)
_write_csv(csv_out, result)

println(json_out)
println(csv_out)
