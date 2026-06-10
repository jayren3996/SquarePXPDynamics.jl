# The simplest-thing-that-works candidate: plain native simple-update at D>=3 with
# rel_floor=0 (no floor artifact), full 4x4 Neel trajectory t=0 -> 3.0 vs exact ED.
#
# Context: the rel_floor=1e-3 grid run (native_d34_0p1grid_to_2.8.out) made D=3==D=4
# to 6 digits (floor caps effective rank) and lagged ED by ~1e-2 on the rise; the
# rel_floor=0 D=3 control reached only t=1.6 (env_d3_control_to2.0.out) where it was
# essentially exact (0.155540 vs ED 0.15558). This run closes the gap: full-trajectory
# bare-D3/D4 at rel_floor=0, plus a dt=0.01 D=3 pass as the Trotter control.
# Pass order puts D=4 last: rel_floor=0 conditioning at D=4 is the one untested risk
# (the legacy backend crashed there; native is unproven beyond D=3).

using SquarePXPDynamics
using TensorKit
using Printf

# True 4x4 Neel ED from artifacts/neel_to_revival_4x4.json (0.2 grid, t=0->3.0).
const ED = Dict(
    0.2=>0.48026525133620573, 0.4=>0.4241792325244558, 0.6=>0.34070816732757286,
    0.8=>0.2442464623167832,  1.0=>0.15549264928630516, 1.2=>0.10057443671942552,
    1.4=>0.10056019790056131, 1.6=>0.1555774954125962,  1.8=>0.24370069054283897,
    2.0=>0.3373453880372616,  2.2=>0.4159980265031265,  2.4=>0.46695735678274974,
    2.6=>0.4830586210094464,  2.8=>0.46229670220636837, 3.0=>0.40822211833244476)

function run(D::Int, dt::Float64; tmax = 3.0, chunk = 0.1)
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = D)
    @printf("\n===== native D=%d, dt=%.3g, rel_floor=0 =====\n", D, dt)
    @printf("%-5s %-10s %-9s %-9s %-12s %-8s\n", "t", "n_native", "n_ED", "delta", "chunk_trnerr", "wall_s")
    t = 0.0
    devs = Float64[]
    while t < tmax - 1e-9
        el = @elapsed log = evolve_pepskit!(psi, chunk; dt = dt, order = 2, schedule = :serial,
                                            maxdim = D, cutoff = 1e-12, rel_floor = 0.0)
        t = round(t + chunk, digits = 1)
        n = exact_density_finite(psi; max_sites = 16)
        if haskey(ED, t)
            d = n - ED[t]
            push!(devs, d)
            @printf("%-5.1f %-10.6f %-9.5f %+-9.5f %-12.4e %-8.1f\n", t, n, ED[t], d, log.max_truncerr, el)
        else
            @printf("%-5.1f %-10.6f %-9s %-9s %-12.4e %-8.1f\n", t, n, "-", "-", log.max_truncerr, el)
        end
        flush(stdout)
    end
    rms = sqrt(sum(abs2, devs) / length(devs))
    @printf("SUMMARY D=%d dt=%.3g: trajectory RMS(0.2-grid)=%.3e  max|delta|=%.3e\n",
            D, dt, rms, maximum(abs, devs))
    flush(stdout)
end

run(3, 0.02)   # the headline candidate
run(3, 0.01)   # Trotter control: if n(t) moves, dt matters; if not, D is the only lever
run(4, 0.02)   # D-ladder: is D=4 closer than D=3? (conditioning risk pass, hence last)
println("\nDONE")
