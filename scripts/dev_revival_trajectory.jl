using SquarePXPDynamics
using Printf

# Trajectory comparison (NOT just the t=2.6 endpoint): iPEPS n(t) measured with the
# EXACT 16-site oracle vs the ED n(t) curve, over the full collapse-and-revival
# [0, 2.8]. Tests whether "D=5 best at t=2.6" is a real trajectory win or a
# curve-crossing coincidence (a frequency/phase error can match one instant while
# the whole curve is worse). ED reference: artifacts/neel_to_revival_4x4.json
# (ed.density, basis_size 743), sampled every 0.2.

const TS = collect(0.0:0.2:2.8)   # 15 points
const ED = [0.5, 0.48026525133620573, 0.4241792325244558, 0.34070816732757286,
            0.2442464623167832, 0.15549264928630516, 0.10057443671942552,
            0.10056019790056131, 0.1555774954125962, 0.24370069054283897,
            0.3373453880372616, 0.4159980265031265, 0.46695735678274974,
            0.4830586210094464, 0.46229670220636837]
const DT = 0.02
const STEP_PER_SAMPLE = 10   # 0.2 / 0.02

function ipeps_trajectory(D)
    psi = checkerboard_square_ipeps(PeriodicSquareUnitCell(4, 4); excited_on = :even, maxdim = D)
    params = TrotterParams(DT, 2, :real, D, 1e-12; schedule = :serial)  # default rel_floor 1e-3
    measure() = (GC.gc(); exact_density_finite(psi; max_sites = 16))   # GC: 2^16 contraction is memory-heavy at D=5
    dens = Float64[measure()]   # t=0
    for _ = 2:length(TS)
        for _ = 1:STEP_PER_SAMPLE
            evolve!(psi, DT; params = params)
        end
        push!(dens, measure())
    end
    return dens
end

results = Dict{Int,Vector{Float64}}()
const I26 = findfirst(==(2.6), TS)
function save_results()
    open("/tmp/traj_data.json", "w") do io
        print(io, "{\"ts\":", TS, ",\"ed\":", ED)
        for D in sort(collect(keys(results)))
            print(io, ",\"D", D, "\":", results[D])
        end
        print(io, "}")
    end
end

@printf("%-4s %-12s %-12s %-12s\n", "D", "max|err|_traj", "rms|err|_traj", "err@2.6")
for D in (2, 3, 4, 5)
    try
        traj = ipeps_trajectory(D)
        results[D] = traj
        errs = abs.(traj .- ED)
        @printf("%-4d %-12.4e %-12.4e %-12.4e\n", D, maximum(errs),
            sqrt(sum(errs .^ 2) / length(errs)), errs[I26])
        save_results()   # incremental: a later-D OOM must not lose earlier D
    catch e
        @printf("%-4d FAILED: %s\n", D, first(split(sprint(showerror, e), "\n")))
    end
    GC.gc()
    flush(stdout)
end
save_results()
println("SAVED /tmp/traj_data.json")
println("DONE")
