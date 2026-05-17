#!/usr/bin/env julia

using Pkg

project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using SquarePXPDynamics

ctm_threading = configure_ctm_threading_from_env!()

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

function _env_list(::Type{T}, name::String, default::String) where {T}
    raw = strip(_env_value(name, default))
    isempty(raw) && return T[]
    return parse.(T, split(raw, ","))
end

function _env_optional_int(name::String, default::String)
    raw = strip(_env_value(name, default))
    return isempty(raw) || lowercase(raw) == "nothing" ? nothing : parse(Int, raw)
end

config = PXPLargerDBenchmarkConfig(;
    n_values = _env_list(Int, "SQUAREPXP_LARGERD_N", "3"),
    total_time = _env_float("SQUAREPXP_LARGERD_TOTAL_TIME", 0.02),
    dt_values = _env_list(Float64, "SQUAREPXP_LARGERD_DT", "0.02"),
    D_values = _env_list(Int, "SQUAREPXP_LARGERD_D", "1,2,3,4"),
    cutoff_values = _env_list(Float64, "SQUAREPXP_LARGERD_CUTOFF", "1e-12"),
    measure_every = _env_int("SQUAREPXP_LARGERD_MEASURE_EVERY", 1),
    order = _env_int("SQUAREPXP_LARGERD_ORDER", 1),
    schedule = _env_symbol("SQUAREPXP_LARGERD_SCHEDULE", :serial),
    initial_state = _env_symbol("SQUAREPXP_LARGERD_INITIAL_STATE", :down),
    point_group = _env_bool("SQUAREPXP_LARGERD_POINT_GROUP", true),
    use_sparse = _env_bool("SQUAREPXP_LARGERD_USE_SPARSE", true),
    ed_tol = _env_float("SQUAREPXP_LARGERD_ED_TOL", 1e-10),
    ed_m_init = _env_int("SQUAREPXP_LARGERD_ED_M_INIT", 30),
    ed_m_max = _env_int("SQUAREPXP_LARGERD_ED_M_MAX", 60),
    ed_extend_step = _env_int("SQUAREPXP_LARGERD_ED_EXTEND_STEP", 10),
    ed_mode = _env_symbol("SQUAREPXP_LARGERD_ED_MODE", :symmetric_pbc),
    observable_mode = _env_symbol("SQUAREPXP_LARGERD_OBSERVABLE_MODE", :auto),
    chi_values = _env_list(Int, "SQUAREPXP_LARGERD_CHI", ""),
    ctm_tol = _env_float("SQUAREPXP_LARGERD_CTM_TOL", 1e-8),
    ctm_maxiter = _env_int("SQUAREPXP_LARGERD_CTM_MAXITER", 100),
    ctm_verbosity = _env_int("SQUAREPXP_LARGERD_CTM_VERBOSITY", 0),
    ctm_seed = _env_optional_int("SQUAREPXP_LARGERD_CTM_SEED", "0"),
    exact_finite_observables = _env_bool("SQUAREPXP_LARGERD_EXACT_FINITE", false),
    exact_finite_max_sites = _env_int("SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES", 12),
)

json_out = _env_value(
    "SQUAREPXP_LARGERD_JSON",
    joinpath(project_root, "artifacts", "pxp_larger_d_ed_benchmark.json"),
)
csv_out = _env_value(
    "SQUAREPXP_LARGERD_CSV",
    joinpath(project_root, "artifacts", "pxp_larger_d_ed_benchmark.csv"),
)

mkpath(dirname(json_out))
mkpath(dirname(csv_out))

report = run_pxp_larger_d_benchmark(config)
write_pxp_larger_d_benchmark_json(report, json_out)
write_pxp_larger_d_benchmark_csv(report, csv_out)

println("CTM threading: ", ctm_threading)
println(json_out)
println(csv_out)
