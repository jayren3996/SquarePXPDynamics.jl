using Printf
using SquarePXPDynamics

println("="^72)
println("Finite MPS TFIM 6x6 Benchmark")
println("="^72)

J = 1.0
h = 1.0
Lx, Ly = 6, 6
initial_state = :x_plus
total_time = 0.5
dt = 0.05
measure_every = 1
maxdim = 64
cutoff = 1e-9

println("\nParameters:")
println("  J = $J, h = $h")
println("  Cell = $(Lx)x$(Ly), boundary = open")
println("  initial_state = $initial_state")
println("  total_time = $total_time, dt = $dt")
println("  measure_every = $measure_every")
println("  TDVP maxdim = $maxdim, cutoff = $cutoff")

result = run_finite_mps_tfim_reference(
    Lx,
    Ly;
    J = J,
    h = h,
    initial_state = initial_state,
    total_time = total_time,
    dt = dt,
    measure_every = measure_every,
    maxdim = maxdim,
    cutoff = cutoff,
)

println("\nSamples:")
for sample in result.samples
    println(@sprintf(
        "  t = %.2f step = %3d  <X> = %.6f  <Z> = %.6f  E/N = %.6f  Dmax = %d  norm = %.8f",
        sample.time,
        sample.step,
        sample.mean_x,
        sample.mean_z,
        sample.energy_density,
        sample.maxlinkdim,
        sample.norm,
    ))
end

initial = first(result.samples)
final = last(result.samples)
energies = [sample.energy_density for sample in result.samples]
norms = [sample.norm for sample in result.samples]

println("\n" * "="^72)
println("Finite MPS Summary")
println("="^72)
println(@sprintf("  Initial <X>:       %.6f", initial.mean_x))
println(@sprintf("  Final <X>:         %.6f", final.mean_x))
println(@sprintf("  Initial E/N:       %.6f", initial.energy_density))
println(@sprintf("  Final E/N:         %.6f", final.energy_density))
println(@sprintf("  Energy max-min:    %.6e", maximum(energies) - minimum(energies)))
println(@sprintf("  Final maxlinkdim:  %d", final.maxlinkdim))
println("  Hit maxdim:        $(any(sample.maxlinkdim >= maxdim for sample in result.samples))")
println(@sprintf("  Norm drift:        %.6e", maximum(abs.(norms .- 1.0))))

pass = all(isfinite, energies) &&
       all(isfinite, norms) &&
       maximum(abs.(norms .- 1.0)) < 1e-4 &&
       final.maxlinkdim <= maxdim

println("\nOVERALL: $(pass ? "PASS" : "FAIL")")

(; pass, result)
