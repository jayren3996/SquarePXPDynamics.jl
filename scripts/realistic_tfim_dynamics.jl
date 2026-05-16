using SquarePXPDynamics
using Printf

println("="^70)
println("Realistic TFIM Dynamics Test")
println("="^70)

# --- Parameters: short-time interacting TFIM ---
J = 1.0
h = 1.0
Lx, Ly = 3, 3
maxdim = 4
dt = 0.01
order = 2
schedule = :serial
total_time = 0.5
measure_every = 5

println("\nParameters:")
println("  J = $J, h = $h")
println("  Cell = $(Lx)x$(Ly)")
println("  maxdim = $maxdim")
println("  dt = $dt, order = $order, schedule = $schedule")
println("  total_time = $total_time")
println("  measure_every = $measure_every")

# --- Setup ---
cell = PeriodicSquareUnitCell(Lx, Ly)
psi = product_square_ipeps(cell; state = :x_plus, maxdim = maxdim)
model = TFIMStarModel(J, h)
protocol = StaticModel(model)
params = TrotterParams(dt, order, :real, maxdim, 1e-12; schedule = schedule)
reference_samples = run_finite_tfim_reference(
    cell,
    model;
    initial_state = :x_plus,
    total_time = total_time,
    dt = dt,
    measure_every = measure_every,
)
reference_by_step = Dict(sample.step => sample.observables for sample in reference_samples)

println("\nInitial state: |x+> product state")
println("Initial maxdim = $(psi.maxdim)")

# --- Measure initial state ---
obs0 = measure_tfim_simple(psi, model)
ref0 = reference_by_step[0]
println("\n--- t = 0.0 ---")
println(@sprintf("  mean_x  = %.6f", obs0.mean_x))
println(@sprintf("  mean_z  = %.6f", obs0.mean_z))
println(@sprintf("  zz_right = %.6f", obs0.zz_right))
println(@sprintf("  zz_up    = %.6f", obs0.zz_up))
println(@sprintf("  energy (star)        = %.6f", obs0.energy_density_star))
println(@sprintf("  energy (decomposed)  = %.6f", obs0.energy_density_decomposed))
println(@sprintf("  energy discrepancy   = %.6f", obs0.energy_density_discrepancy))
println(@sprintf("  |mean_x - ED|        = %.6e", abs(obs0.mean_x - ref0.mean_x)))
println(@sprintf("  |energy - ED|        = %.6e", abs(obs0.energy_density_star - ref0.energy_density_star)))
println(@sprintf("  max_bond_entropy     = %.6f", obs0.max_bond_entropy))
println(@sprintf("  mean_bond_entropy    = %.6f", obs0.mean_bond_entropy))

# --- Time evolution with periodic measurement ---
nsteps = round(Int, total_time / dt)
samples = [(0.0, obs0, ref0, 0.0)]  # (time, observables, ED observables, max_truncerr)

println("\n--- Evolving $(nsteps) steps ---")

for step in 1:nsteps
    log = evolve!(psi, dt; params = params, protocol = protocol)

    if step % measure_every == 0 || step == nsteps
        time = step * dt
        obs = measure_tfim_simple(psi, model)
        ref = reference_by_step[step]
        push!(samples, (time, obs, ref, log.max_truncerr))

        println(@sprintf("\n--- t = %.2f (step %d) ---", time, step))
        println(@sprintf("  mean_x  = %.6f", obs.mean_x))
        println(@sprintf("  mean_z  = %.6f", obs.mean_z))
        println(@sprintf("  ED mean_x             = %.6f", ref.mean_x))
        println(@sprintf("  ED mean_z             = %.6f", ref.mean_z))
        println(@sprintf("  |mean_x - ED|        = %.6e", abs(obs.mean_x - ref.mean_x)))
        println(@sprintf("  |mean_z - ED|        = %.6e", abs(obs.mean_z - ref.mean_z)))
        println(@sprintf("  energy (star)        = %.6f", obs.energy_density_star))
        println(@sprintf("  ED energy            = %.6f", ref.energy_density_star))
        println(@sprintf("  |energy - ED|        = %.6e", abs(obs.energy_density_star - ref.energy_density_star)))
        println(@sprintf("  energy discrepancy   = %.6f", obs.energy_density_discrepancy))
        println(@sprintf("  max_truncerr         = %.2e", log.max_truncerr))
        println(@sprintf("  max_bond_entropy     = %.6f", obs.max_bond_entropy))
        println(@sprintf("  mean_bond_entropy    = %.6f", obs.mean_bond_entropy))
        println(@sprintf("  log_norm_delta       = %.6f", log.log_norm_delta))
    end
end

# --- Physical sanity checks ---
println("\n" * "="^70)
println("Physical Sanity Checks")
println("="^70)

energies = [obs.energy_density_star for (_, obs, _, _) in samples]
energy_drift = maximum(energies) - minimum(energies)
energy_mean = sum(energies) / length(energies)
ed_mean_x_deltas = [abs(obs.mean_x - ref.mean_x) for (_, obs, ref, _) in samples]
ed_mean_z_deltas = [abs(obs.mean_z - ref.mean_z) for (_, obs, ref, _) in samples]
ed_energy_deltas = [
    abs(obs.energy_density_star - ref.energy_density_star) for (_, obs, ref, _) in samples
]

println(@sprintf("\nEnergy conservation:"))
println(@sprintf("  Initial energy:  %.6f", energies[1]))
println(@sprintf("  Final energy:    %.6f", energies[end]))
println(@sprintf("  Max - Min:       %.6f", energy_drift))
println(@sprintf("  Relative drift:  %.2e", energy_drift / abs(energy_mean)))

println(@sprintf("\nFinite ED reference deltas:"))
println(@sprintf("  Final |mean_x - ED|:   %.6e", ed_mean_x_deltas[end]))
println(@sprintf("  Final |mean_z - ED|:   %.6e", ed_mean_z_deltas[end]))
println(@sprintf("  Final |energy - ED|:   %.6e", ed_energy_deltas[end]))
println(@sprintf("  Max |mean_x - ED|:     %.6e", maximum(ed_mean_x_deltas)))
println(@sprintf("  Max |energy - ED|:     %.6e", maximum(ed_energy_deltas)))

max_entropies = [obs.max_bond_entropy for (_, obs, _, _) in samples]
mean_entropies = [obs.mean_bond_entropy for (_, obs, _, _) in samples]
max_truncerrs = [tr for (_, _, _, tr) in samples]

println(@sprintf("\nBond entropy growth:"))
println(@sprintf("  Initial max entropy:  %.6f", max_entropies[1]))
println(@sprintf("  Final max entropy:    %.6f", max_entropies[end]))
println(@sprintf("  Initial mean entropy: %.6f", mean_entropies[1]))
println(@sprintf("  Final mean entropy:   %.6f", mean_entropies[end]))

println(@sprintf("\nTruncation errors:"))
println(@sprintf("  Max over trajectory:  %.2e", maximum(max_truncerrs)))

# Check all observables finite
all_finite = all(s -> all(isfinite, (
    s.mean_x, s.mean_y, s.mean_z,
    s.zz_right, s.zz_up,
    s.energy_density_star, s.energy_density_decomposed,
    s.mean_bond_entropy, s.max_bond_entropy
)), [obs for (_, obs, _, _) in samples])

println("\nAll observables finite: $all_finite")

# Approximate pass/fail criteria for simple-update dynamics
# Energy drift < 10% is a rough sanity check for short-time simple update
# Bond entropy < 2.0 for maxdim=4 is reasonable (log(4) ≈ 1.39)
pass = all_finite &&
       energy_drift / abs(energy_mean) < 0.1 &&
       max_entropies[end] < 2.0 &&
       maximum(max_truncerrs) < 1.0

println("\n" * "="^70)
println("OVERALL: $(pass ? "PASS" : "FAIL")")
println("="^70)

# Return results
(; pass, samples, energy_drift, ed_mean_x_deltas, ed_mean_z_deltas, ed_energy_deltas, max_entropies, mean_entropies, max_truncerrs)
