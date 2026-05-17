using SquarePXPDynamics

function _parse_list(::Type{T}, name::String, default::String) where {T}
    raw = strip(get(ENV, name, default))
    isempty(raw) && return T[]
    return parse.(T, split(raw, ","))
end

function _parse_optional_int(name::String, default::String)
    raw = strip(get(ENV, name, default))
    return isempty(raw) || lowercase(raw) == "nothing" ? nothing : parse(Int, raw)
end

config = PXPAuditConfig(;
    n_values = _parse_list(Int, "SQUAREPXP_AUDIT_N", "3"),
    total_time = parse(Float64, get(ENV, "SQUAREPXP_AUDIT_TOTAL_TIME", "0.02")),
    dt_values = _parse_list(Float64, "SQUAREPXP_AUDIT_DT", "0.02,0.01"),
    D_values = _parse_list(Int, "SQUAREPXP_AUDIT_D", "1,2"),
    cutoff_values = _parse_list(Float64, "SQUAREPXP_AUDIT_CUTOFF", "1e-12"),
    chi_values = _parse_list(Int, "SQUAREPXP_AUDIT_CHI", ""),
    measure_every = parse(Int, get(ENV, "SQUAREPXP_AUDIT_MEASURE_EVERY", "1")),
    order = parse(Int, get(ENV, "SQUAREPXP_AUDIT_ORDER", "1")),
    schedule = Symbol(get(ENV, "SQUAREPXP_AUDIT_SCHEDULE", "serial")),
    ctm_tol = parse(Float64, get(ENV, "SQUAREPXP_AUDIT_CTM_TOL", "1e-8")),
    ctm_maxiter = parse(Int, get(ENV, "SQUAREPXP_AUDIT_CTM_MAXITER", "100")),
    ctm_verbosity = parse(Int, get(ENV, "SQUAREPXP_AUDIT_CTM_VERBOSITY", "0")),
    ctm_seed = _parse_optional_int("SQUAREPXP_AUDIT_CTM_SEED", "0"),
)

json_out = get(ENV, "SQUAREPXP_AUDIT_JSON", "artifacts/pxp_audit_report.json")
csv_out = get(ENV, "SQUAREPXP_AUDIT_CSV", "artifacts/pxp_audit_summary.csv")
mkpath(dirname(json_out))
mkpath(dirname(csv_out))

report = run_pxp_audit_campaign(config)
write_pxp_audit_json(report, json_out)
write_pxp_audit_csv(report, csv_out)

println(json_out)
println(csv_out)
