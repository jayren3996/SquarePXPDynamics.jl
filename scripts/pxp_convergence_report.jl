using SquarePXPDynamics

base = PXPValidationConfig(
    parse(Int, get(ENV, "SQUAREPXP_CONVERGENCE_N", "3"));
    total_time = parse(Float64, get(ENV, "SQUAREPXP_CONVERGENCE_TOTAL_TIME", "0.02")),
    dt = parse(Float64, get(ENV, "SQUAREPXP_CONVERGENCE_BASE_DT", "0.01")),
)
sweep = PXPConvergenceConfig(
    base;
    dt_values = parse.(Float64, split(get(ENV, "SQUAREPXP_CONVERGENCE_DT", "0.01,0.005"), ",")),
    D_values = parse.(Int, split(get(ENV, "SQUAREPXP_CONVERGENCE_D", "1"), ",")),
    chi_values = Int[],
    cutoff_values = parse.(Float64, split(get(ENV, "SQUAREPXP_CONVERGENCE_CUTOFF", "1e-12"), ",")),
)
report = validate_pxp_convergence(sweep)
haskey(ENV, "SQUAREPXP_CONVERGENCE_OUT") ||
    error("set SQUAREPXP_CONVERGENCE_OUT to the output JSON path")
out = ENV["SQUAREPXP_CONVERGENCE_OUT"]
write_pxp_convergence_json(report, out)
println(out)
