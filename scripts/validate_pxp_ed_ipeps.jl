using Pkg

project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using SquarePXPDynamics

function _env_int(name::String, default::Int)
    return parse(Int, get(ENV, name, string(default)))
end

function _env_float(name::String, default::Float64)
    return parse(Float64, get(ENV, name, string(default)))
end

function _env_symbol(name::String, default::Symbol)
    return Symbol(get(ENV, name, String(default)))
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
    order = _env_int("SQUAREPXP_PXP_VALIDATION_ORDER", 1),
    maxdim = _env_int("SQUAREPXP_PXP_VALIDATION_MAXDIM", 1),
    cutoff = _env_float("SQUAREPXP_PXP_VALIDATION_CUTOFF", 1e-12),
    schedule = _env_symbol("SQUAREPXP_PXP_VALIDATION_SCHEDULE", :serial),
)

report = validate_pxp_ed_ipeps(config; ctm_params = nothing)
mkpath(dirname(out))
write_pxp_validation_json(report, out)
println(out)
