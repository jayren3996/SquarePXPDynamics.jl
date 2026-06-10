# 6x6-torus Neel quench, plain D rel_floor=0 simple update, vs the exact 6x6
# Krylov ED oracle (artifacts/neel_ed_krylov_6.json; m=40 converged Lanczos).
#
# Observable: sector_density (scripts/g0/sector_density.jl) — the EXACT
# constrained-sector meet-in-the-middle contraction, gated at machine precision
# against the dense oracle on 4x4 and to 1.8e-9 against the 6x6 Krylov value at
# t=0.2 (sector_density_gate.out; adversarially reviewed). ~80 s / ~41 GB per
# 6x6 D=3 eval. Also reports the seam constraint-leakage fraction per point.
#
# Env knobs (parallel staggered-grid lanes / D-controls):
#   PXP_D       bond dimension (default 3)
#   MEAS_EVERY  measurement interval (default 0.1); evolution dt=0.02 throughout
#   MEAS_OFFSET first measurement time (default = MEAS_EVERY); a lane with
#               MEAS_EVERY=0.2, MEAS_OFFSET=0.1 covers the odd 0.1-grid points
#               so two lanes fill the full grid in half the wall-clock.
#   BLAS_T      BLAS threads for the observable (default 16)

using SquarePXPDynamics
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "sector_density.jl"))

# 6x6 Krylov ED oracle, 0.1 grid.
const ED6 = Dict(
    0.1=>0.49501664446317023, 0.2=>0.4802652513362058,  0.3=>0.45633405877827876,
    0.4=>0.42417923263430746, 0.5=>0.38509718401185467, 0.6=>0.3407082117509682,
    0.7=>0.29297303898284155, 0.8=>0.24424863202324137, 0.9=>0.19734863197300032,
    1.0=>0.15552170701479606, 1.1=>0.12223806948455485, 1.2=>0.10072585850681608,
    1.3=>0.09333270965600968, 1.4=>0.10093633578255576, 1.5=>0.12267381112668196,
    1.6=>0.15613395062461333, 1.7=>0.1979272820254735,  1.8=>0.24437552753721997,
    1.9=>0.29205964682135716, 2.0=>0.3381018959507266,  2.1=>0.3802121259426343,
    2.2=>0.41660816030552783, 2.3=>0.4459142560506063,  2.4=>0.4670937134447127,
    2.5=>0.47942591924386085, 2.6=>0.4825137728423614,  2.7=>0.47630227858494406,
    2.8=>0.46109325160775944, 2.9=>0.43754763063008517, 3.0=>0.40667424323431356)

function main(; tmax = 3.0, dt = 0.02)
    D = parse(Int, get(ENV, "PXP_D", "3"))
    meas = parse(Float64, get(ENV, "MEAS_EVERY", "0.1"))
    offset = parse(Float64, get(ENV, "MEAS_OFFSET", string(meas)))
    blas_t = parse(Int, get(ENV, "BLAS_T", "16"))
    cell = PeriodicSquareUnitCell(6, 6)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = D)
    @printf("=== 6x6 native D=%d rel_floor=0 dt=%.3g meas=%.2g offset=%.2g vs Krylov ED (sector observable, blas=%d) ===\n",
            D, dt, meas, offset, blas_t)
    @printf("%-5s %-10s %-9s %+-10s %-9s %-10s %-8s %-8s\n",
            "t", "n_sector", "n_ED6", "delta", "leak", "n_simple", "ev_s", "obs_s")
    t = 0.0
    devs = Float64[]
    while true
        step = t == 0.0 ? offset : meas
        tnext = round(t + step, digits = 1)
        tnext > tmax + 1e-9 && break
        ev = @elapsed evolve_pepskit!(psi, step; dt = dt, order = 2, schedule = :serial,
                                      maxdim = D, cutoff = 1e-12, rel_floor = 0.0)
        t = tnext
        nsim = measure_simple(psi).density
        r = sector_density(psi; blas_threads = blas_t, verbose = false)
        d = r.n - ED6[t]
        push!(devs, d)
        @printf("%-5.1f %-10.6f %-9.5f %+-10.5f %-9.2e %-10.6f %-8.1f %-8.1f\n",
                t, r.n, ED6[t], d, r.leak, nsim, ev, r.secs)
        flush(stdout)
    end
    rms = sqrt(sum(abs2, devs) / length(devs))
    @printf("SUMMARY 6x6 D=%d (offset %.2g, every %.2g): RMS=%.3e  max|delta|=%.3e over %d points\n",
            D, offset, meas, rms, maximum(abs, devs), length(devs))
    println("DONE")
end

main()
