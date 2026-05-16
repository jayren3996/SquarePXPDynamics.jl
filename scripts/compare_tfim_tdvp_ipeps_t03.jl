using Printf
using SquarePXPDynamics

println("="^86)
println("TFIM TDVP-vs-iPEPS Comparison to t = 0.30")
println("="^86)

J = 1.0
h = 1.0
initial_state = :x_plus

tdvp_Lx, tdvp_Ly = 6, 6
tdvp_dt = 0.05
tdvp_total_time = 0.30
tdvp_maxdim = 64
tdvp_cutoff = 1e-9

ipeps_Lx, ipeps_Ly = 3, 3
ipeps_dt = 0.01
ipeps_total_time = tdvp_total_time
ipeps_measure_every = round(Int, tdvp_dt / ipeps_dt)
ipeps_maxdim = 4
ipeps_cutoff = 1e-12
ipeps_order = 2
ipeps_schedule = :serial

println("\nShared physical setup:")
println("  TFIM convention: H = -h * sum_i X_i - J * sum_<ij> Z_i Z_j")
println("  J = $J, h = $h, initial_state = $initial_state")
println("  comparison grid: t = 0:$(tdvp_dt):$(tdvp_total_time)")

println("\nTDVP reference:")
println("  finite MPS = $(tdvp_Lx)x$(tdvp_Ly), boundary = open")
println("  dt = $tdvp_dt, maxdim = $tdvp_maxdim, cutoff = $tdvp_cutoff")

println("\niPEPS run:")
println("  periodic unit cell = $(ipeps_Lx)x$(ipeps_Ly)")
println("  dt = $ipeps_dt, order = $ipeps_order, schedule = $ipeps_schedule")
println("  maxdim = $ipeps_maxdim, cutoff = $ipeps_cutoff")

tdvp = run_finite_mps_tfim_reference(
    tdvp_Lx,
    tdvp_Ly;
    J = J,
    h = h,
    initial_state = initial_state,
    total_time = tdvp_total_time,
    dt = tdvp_dt,
    measure_every = 1,
    maxdim = tdvp_maxdim,
    cutoff = tdvp_cutoff,
)

ipeps_spec = BenchmarkSpec(
    "tfim-serial-3x3-t03",
    StaticModel(TFIMStarModel(J, h)),
    PeriodicSquareUnitCell(ipeps_Lx, ipeps_Ly),
    initial_state,
    ipeps_total_time,
    TrotterParams(
        ipeps_dt,
        ipeps_order,
        :real,
        ipeps_maxdim,
        ipeps_cutoff;
        schedule = ipeps_schedule,
    ),
    ipeps_measure_every,
)
ipeps = run_benchmark(ipeps_spec; run_label = "tdvp-comparison")

tdvp_by_time = Dict(round(sample.time; digits = 12) => sample for sample in tdvp.samples)
ipeps_by_time = Dict(round(sample.time; digits = 12) => sample for sample in ipeps.samples)
times = sort(collect(intersect(keys(tdvp_by_time), keys(ipeps_by_time))))

isempty(times) && error("no shared sample times between TDVP and iPEPS trajectories")

println("\nComparison table:")
println(
    "  time    TDVP <X>    iPEPS <X>   |dX|        TDVP <Z>    iPEPS <Z>   |dZ|        TDVP E/N    iPEPS E/N   |dE|        D_MPS  D_iPEPS",
)
for time in times
    mps_sample = tdvp_by_time[time]
    ipeps_sample = ipeps_by_time[time]
    obs = ipeps_sample.observables
    println(@sprintf(
        "  %.2f    % .6f   % .6f   %.3e   % .6f   % .6f   %.3e   % .6f   % .6f   %.3e   %5d  %7d",
        time,
        mps_sample.mean_x,
        obs.mean_x,
        abs(mps_sample.mean_x - obs.mean_x),
        mps_sample.mean_z,
        obs.mean_z,
        abs(mps_sample.mean_z - obs.mean_z),
        mps_sample.energy_density,
        obs.energy_density_star,
        abs(mps_sample.energy_density - obs.energy_density_star),
        mps_sample.maxlinkdim,
        ipeps.metadata.maxdim,
    ))
end

dx = Float64[]
dz = Float64[]
de = Float64[]
for time in times
    mps_sample = tdvp_by_time[time]
    obs = ipeps_by_time[time].observables
    push!(dx, abs(mps_sample.mean_x - obs.mean_x))
    push!(dz, abs(mps_sample.mean_z - obs.mean_z))
    push!(de, abs(mps_sample.energy_density - obs.energy_density_star))
end

tdvp_energies = [sample.energy_density for sample in tdvp.samples]
tdvp_norms = [sample.norm for sample in tdvp.samples]
ipeps_energies = [sample.observables.energy_density_star for sample in ipeps.samples]
ipeps_truncerrs = [sample.diagnostics.max_truncerr for sample in ipeps.samples]

println("\nSummary:")
println(@sprintf("  max |TDVP <X> - iPEPS <X>|:  %.6e", maximum(dx)))
println(@sprintf("  max |TDVP <Z> - iPEPS <Z>|:  %.6e", maximum(dz)))
println(@sprintf("  max |TDVP E/N - iPEPS E/N|:  %.6e", maximum(de)))
println(@sprintf("  final |TDVP <X> - iPEPS <X>|: %.6e", dx[end]))
println(@sprintf("  final |TDVP E/N - iPEPS E/N|: %.6e", de[end]))
println(@sprintf("  TDVP energy max-min:          %.6e", maximum(tdvp_energies) - minimum(tdvp_energies)))
println(@sprintf("  TDVP norm drift:              %.6e", maximum(abs.(tdvp_norms .- 1.0))))
println(@sprintf("  TDVP hit maxdim:              %s", any(sample.maxlinkdim >= tdvp_maxdim for sample in tdvp.samples)))
println(@sprintf("  iPEPS energy max-min:         %.6e", maximum(ipeps_energies) - minimum(ipeps_energies)))
println(@sprintf("  iPEPS max truncerr:           %.6e", maximum(ipeps_truncerrs)))

all_finite = all(isfinite, dx) &&
             all(isfinite, dz) &&
             all(isfinite, de) &&
             all(isfinite, tdvp_energies) &&
             all(isfinite, tdvp_norms) &&
             all(isfinite, ipeps_energies) &&
             all(isfinite, ipeps_truncerrs)

println("\nInterpretation:")
println("  TDVP is a finite open-boundary 6x6 reference; iPEPS is infinite periodic 3x3.")
println("  The comparison is meaningful as a short-time local-observable benchmark,")
println("  not as an exact thermodynamic-limit equality test.")

println("\nOVERALL: $(all_finite ? "PASS" : "FAIL")")

(; all_finite, tdvp, ipeps, times, dx, dz, de)
