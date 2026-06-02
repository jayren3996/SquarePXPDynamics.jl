using SquarePXPDynamics
using Printf

# Re-evaluate the rel_floor default by the n(t) TRAJECTORY (not the misleading
# t=2.6 endpoint). The committed 1e-4 -> 1e-3 change was decided at t=2.6, a
# revival-peak crossing point. Does 1e-3 actually beat 1e-4 (and 0) by trajectory
# RMS error? 4x4 Neel, exact oracle, every 0.2 over [0,2.8], D=2,3,4 (D=5 OOMs).

const TS = collect(0.0:0.2:2.8)
const ED = [0.5, 0.48026525133620573, 0.4241792325244558, 0.34070816732757286,
            0.2442464623167832, 0.15549264928630516, 0.10057443671942552,
            0.10056019790056131, 0.1555774954125962, 0.24370069054283897,
            0.3373453880372616, 0.4159980265031265, 0.46695735678274974,
            0.4830586210094464, 0.46229670220636837]
const DT = 0.02
const SPS = 10
const I26 = findfirst(==(2.6), TS)

function traj_errs(rel_floor, D)
    psi = checkerboard_square_ipeps(PeriodicSquareUnitCell(4, 4); excited_on = :even, maxdim = D)
    params = TrotterParams(DT, 2, :real, D, 1e-12; schedule = :serial, rel_floor = rel_floor)
    dens = Float64[(GC.gc(); exact_density_finite(psi; max_sites = 16))]
    for _ = 2:length(TS)
        for _ = 1:SPS
            evolve!(psi, DT; params = params)
        end
        push!(dens, (GC.gc(); exact_density_finite(psi; max_sites = 16)))
    end
    e = abs.(dens .- ED)
    return sqrt(sum(e .^ 2) / length(e)), maximum(e), e[I26]
end

@printf("%-9s %-3s %-12s %-12s %-12s\n", "rel_floor", "D", "rms_traj", "max_traj", "err@2.6")
for rf in (1e-4, 1e-3, 0.0), D in (2, 3, 4)
    try
        rms, mx, e26 = traj_errs(rf, D)
        @printf("%-9.0e %-3d %-12.4e %-12.4e %-12.4e\n", rf, D, rms, mx, e26)
    catch e
        @printf("%-9.0e %-3d FAILED: %s\n", rf, D, first(split(sprint(showerror, e), "\n")))
    end
    flush(stdout)
end
println("Pick the floor with the lowest TRAJECTORY rms (not err@2.6). If 1e-4 ~ 1e-3")
println("by trajectory, the committed default change was an endpoint artifact.")
println("DONE")
