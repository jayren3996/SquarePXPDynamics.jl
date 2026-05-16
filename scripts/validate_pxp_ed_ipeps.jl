using Pkg

project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using SquarePXPDynamics

function _env_value(name::String, default::AbstractString)
    value = get(ENV, name, "")
    return isempty(value) ? String(default) : value
end

function _env_int(name::String, default::Int)
    return parse(Int, _env_value(name, string(default)))
end

function _env_float(name::String, default::Float64)
    return parse(Float64, _env_value(name, string(default)))
end

function _env_bool(name::String, default::Bool)
    value = lowercase(strip(_env_value(name, string(default))))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    throw(ArgumentError("$name must be one of 1,true,yes,on,0,false,no,off"))
end

function _env_symbol(name::String, default::Symbol)
    return Symbol(_env_value(name, String(default)))
end

out = get(
    ENV,
    "SQUAREPXP_PXP_VALIDATION_OUT",
    joinpath(project_root, "artifacts", "pxp_validation_report.json"),
)

config = PXPValidationConfig(
    _env_int("SQUAREPXP_PXP_VALIDATION_N", 3);
    total_time = _env_float("SQUAREPXP_PXP_VALIDATION_TOTAL_TIME", 0.02),
    dt = _env_float("SQUAREPXP_PXP_VALIDATION_DT", 0.01),
    measure_every = _env_int("SQUAREPXP_PXP_VALIDATION_MEASURE_EVERY", 1),
    initial_state = _env_symbol("SQUAREPXP_PXP_VALIDATION_INITIAL_STATE", :down),
    point_group = _env_bool("SQUAREPXP_PXP_VALIDATION_POINT_GROUP", true),
    use_sparse = _env_bool("SQUAREPXP_PXP_VALIDATION_USE_SPARSE", true),
    ed_tol = _env_float("SQUAREPXP_PXP_VALIDATION_ED_TOL", 1e-10),
    ed_m_init = _env_int("SQUAREPXP_PXP_VALIDATION_ED_M_INIT", 30),
    ed_m_max = _env_int("SQUAREPXP_PXP_VALIDATION_ED_M_MAX", 60),
    ed_extend_step = _env_int("SQUAREPXP_PXP_VALIDATION_ED_EXTEND_STEP", 10),
    order = _env_int("SQUAREPXP_PXP_VALIDATION_ORDER", 1),
    maxdim = _env_int("SQUAREPXP_PXP_VALIDATION_MAXDIM", 1),
    cutoff = _env_float("SQUAREPXP_PXP_VALIDATION_CUTOFF", 1e-12),
    schedule = _env_symbol("SQUAREPXP_PXP_VALIDATION_SCHEDULE", :serial),
)

report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
mkpath(dirname(out))
write_pxp_validation_json(report, out)
println(out)
