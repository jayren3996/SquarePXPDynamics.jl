# Gate for scripts/g0/sector_density.jl (exact constrained-sector density).
#
# 4x4 D=3 rel_floor=0 trajectory; at t = 0.5, 1.6, 2.6 compare
#   (i)   n_full = exact_density_finite (full dense oracle, all configs),
#   (ii)  an INDEPENDENT sector reference from dense_state_finite: filter the
#         2^16 configs by the hard-square constraint ON THE TORUS (bond pairs
#         from PeriodicSquareUnitCell neighbor(:right/:up), dense bit positions
#         per _dense_index/_local_positions conventions) -> n_dense_sector,
#         plus the seam-comparable and full-torus leakage references,
#   (iii) (n_sector, leak) = sector_density(psi).
# PASS: |n_sector - n_dense_sector| <= 1e-9 AND |leak - leak_dense| <=
# max(1e-9, 0.05*leak_dense) at all three times, where leak_dense is the
# seam-comparable dense reference (identical config universe as the sector
# method: every PEPS row ring-valid, vertical bonds valid except at the two
# seams rows h|h+1 and Nr|1). The full-torus leakage (all configs / all bonds)
# is printed alongside for information. n_full - n_sector printed, no threshold.
#
# Extra structural checks (free at 4x4): universe/valid config COUNTS must
# match between the two paths, the dense-path density over ALL configs must
# reproduce exact_density_finite (bit-convention guard), and the per-config
# amplitudes Psi[p,s]*exp(logscale) must match dense_state_finite entries.
#
# Then 6x6: evolve to t=0.2 (expect |n - 0.48026525| <= 2e-4, Krylov ED ref;
# SU truncation only), continue to t=1.6 (Krylov n=0.15613395, ~1e-3 expected,
# no hard threshold), report wall seconds and peak RAM.

using SquarePXPDynamics
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "sector_density.jl"))

say(args...) = (println(args...); flush(stdout))

const KRYLOV_6X6 = Dict(0.2 => 0.48026525, 1.6 => 0.15613395)

# ---------------------------------------------------------------------------
# Independent dense-path sector reference (4x4 only).
# ---------------------------------------------------------------------------
function dense_sector_reference(psi::PXPIPEPSState, h::Int)
    cell = psi.unitcell
    Nr, Nc = cell.Ly, cell.Lx
    nsites = Nr * Nc
    sv = dense_state_finite(psi; max_sites = 16)

    pos(coord) = findfirst(==(coord), cell.reps)
    physpos(r, c) = (Nr - r) * Nc + c              # _native_dense_state convention
    # cross-check the (r,c) -> reps-position map against the cell geometry
    for r = 1:Nr, c = 1:Nc
        coord = SquareCoord(c, Nr - r + 1)
        pos(coord) == physpos(r, c) || error("physpos mismatch at ($r,$c)")
        squarecoord_to_cartesianindex(cell, coord) == CartesianIndex(r, c) ||
            error("coordinate map mismatch at ($r,$c)")
    end
    npair(t) = minmax(t[1], t[2])
    geo_up = Set(npair((pos(c), pos(neighbor(cell, c, :up)))) for c in cell.reps)
    geo_right = Set(npair((pos(c), pos(neighbor(cell, c, :right)))) for c in cell.reps)
    vert(r, c) = (physpos(r, c), physpos(mod1(r + 1, Nr), c))   # S(r,c) ~ N(r+1,c)
    horiz(r, c) = (physpos(r, c), physpos(r, mod1(c + 1, Nc)))  # E(r,c) ~ W(r,c+1)
    Set(npair(vert(r, c)) for r = 1:Nr, c = 1:Nc) == geo_up ||
        error("vertical bond set disagrees with neighbor(:up) geometry")
    Set(npair(horiz(r, c)) for r = 1:Nr, c = 1:Nc) == geo_right ||
        error("horizontal bond set disagrees with neighbor(:right) geometry")

    seam = [vert(r, c) for r in (h, Nr) for c = 1:Nc]
    nonseam = vcat(
        [horiz(r, c) for r = 1:Nr for c = 1:Nc],
        [vert(r, c) for r = 1:Nr for c = 1:Nc if r != h && r != Nr],
    )
    # dense bits (_dense_index): idx-1 = sum_p (v_p - 1) * 2^(nsites - p), and
    # v_p = 1 means EXCITED -> excited at reps position p iff bit (nsites-p) is 0.
    bondmask(b) = (UInt64(1) << (nsites - b[1])) | (UInt64(1) << (nsites - b[2]))
    seam_masks = bondmask.(seam)
    nonseam_masks = bondmask.(nonseam)
    fullmask = (UInt64(1) << nsites) - 1

    S_tot = 0.0; S_tot_n = 0.0
    S_univ = 0.0
    S_valid = 0.0; S_valid_n = 0.0
    cnt_univ = 0; cnt_valid = 0
    @inbounds for idx = 1:(1<<nsites)
        w = abs2(sv[idx])
        exc = (~UInt64(idx - 1)) & fullmask
        nex = count_ones(exc)
        S_tot += w
        S_tot_n += w * nex
        violn = false
        for mk in nonseam_masks
            if exc & mk == mk
                violn = true
                break
            end
        end
        violn && continue
        cnt_univ += 1
        S_univ += w
        viols = false
        for mk in seam_masks
            if exc & mk == mk
                viols = true
                break
            end
        end
        viols && continue
        cnt_valid += 1
        S_valid += w
        S_valid_n += w * nex
    end
    return (
        n_full_check = S_tot_n / S_tot / nsites,
        n_sector = S_valid_n / S_valid / nsites,
        leak_seam = (S_univ - S_valid) / S_univ,
        leak_full = (S_tot - S_valid) / S_tot,
        cnt_univ = cnt_univ,
        cnt_valid = cnt_valid,
        sv = sv,
    )
end

# Per-config amplitude cross-check: sector amplitude (rescaled to physical)
# vs the dense statevector entry of the SAME configuration.
function amp_crosscheck(res, sv, Nr::Int, Nc::Int, h::Int)
    nsites = Nr * Nc
    physpos(r, c) = (Nr - r) * Nc + c
    scale = exp(res.logscale)
    maxdiff = 0.0
    maxamp = 0.0
    for s = 1:res.nbot, p = 1:res.ntop
        idxm1 = UInt64(0)
        for (half, off) in ((res.top_stacks[p], 0), (res.bot_stacks[s], h))
            for (k, mi) in enumerate(half)
                m = res.masks[mi]
                r = off + k
                for c = 1:Nc
                    excited = ((m >> (c - 1)) & 0x1) == 0x1
                    if !excited   # v=2 -> bit (nsites - p) set
                        idxm1 |= UInt64(1) << (nsites - physpos(r, c))
                    end
                end
            end
        end
        a_sector = res.amps[p, s] * scale
        a_dense = sv[Int(idxm1)+1]
        maxdiff = max(maxdiff, abs(a_sector - a_dense))
        maxamp = max(maxamp, abs(a_dense))
    end
    return maxdiff, maxamp
end

function checkpoint_4x4(psi, t, h; gating::Bool)
    n_full = exact_density_finite(psi; max_sites = 16)
    ref = dense_sector_reference(psi, h)
    res = sector_density(psi; keep_amplitudes = true, verbose = false)
    # structural identities between the two enumerations
    res.npairs == ref.cnt_univ ||
        error("universe count mismatch: $(res.npairs) vs $(ref.cnt_univ)")
    res.nvalid_pairs == ref.cnt_valid ||
        error("valid count mismatch: $(res.nvalid_pairs) vs $(ref.cnt_valid)")
    abs(ref.n_full_check - n_full) < 1e-10 ||
        error("dense-path full density check failed: $(ref.n_full_check) vs $n_full")
    ampdiff, maxamp = amp_crosscheck(res, ref.sv, 4, 4, h)

    dn = abs(res.n - ref.n_sector)
    dl = abs(res.leak - ref.leak_seam)
    tol_l = max(1e-9, 0.05 * ref.leak_seam)
    ok = dn <= 1e-9 && dl <= tol_l
    say(@sprintf(
        "t=%.1f  n_full=%.12f  n_dense_sector=%.12f  n_sector=%.12f  |dn|=%.3e",
        t, n_full, ref.n_sector, res.n, dn))
    say(@sprintf(
        "t=%.1f  leak_sector=%.6e  leak_dense_seam=%.6e  |dleak|=%.3e (tol %.3e)  leak_dense_full=%.6e",
        t, res.leak, ref.leak_seam, dl, tol_l, ref.leak_full))
    say(@sprintf(
        "t=%.1f  n_full-n_sector=%+.6e  ampdiff=%.3e (maxamp %.3e)  pairs=%d/%d valid=%d  secs=%.2f  %s",
        t, n_full - res.n, ampdiff, maxamp, res.npairs, 2^16, res.nvalid_pairs,
        res.secs, gating ? (ok ? "PASS" : "FAIL") : "(info)"))
    return ok, dn, res.leak, n_full - res.n
end

function gate_4x4()
    say("=== 4x4 gate: D=3 rel_floor=0, dt=0.02 order=2 serial, maxdim=3 cutoff=1e-12 ===")
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = 3)
    h = cell.Ly ÷ 2
    say(@sprintf("masks: Nc=4 -> %d (expect 7)", length(sector_row_masks(4))))
    checkpoint_4x4(psi, 0.0, h; gating = false)   # info-only sanity at t=0
    allpass = true
    maxdn = 0.0
    leak26 = NaN
    gaps = Float64[]
    t = 0.0
    for (tgt, nchunks) in ((0.5, 5), (1.6, 11), (2.6, 10))
        for _ = 1:nchunks
            evolve_pepskit!(
                psi, 0.1;
                dt = 0.02, order = 2, schedule = :serial,
                maxdim = 3, cutoff = 1e-12, rel_floor = 0.0,
            )
        end
        t = tgt
        ok, dn, leak, gap = checkpoint_4x4(psi, t, h; gating = true)
        allpass &= ok
        maxdn = max(maxdn, dn)
        push!(gaps, gap)
        t == 2.6 && (leak26 = leak)
    end
    say(@sprintf("4x4 gate: %s  (max |n_sector - n_dense_sector| = %.3e)",
        allpass ? "PASS" : "FAIL", maxdn))
    return allpass, maxdn, leak26, gaps
end

function run_6x6()
    say("")
    say("=== 6x6: D=3 rel_floor=0, dt=0.02 order=2 serial, maxdim=3 cutoff=1e-12 ===")
    say(@sprintf("masks: Nc=6 -> %d (expect 18)", length(sector_row_masks(6))))
    cell = PeriodicSquareUnitCell(6, 6)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = 3)
    t = 0.0
    evolve_chunk!() = begin
        tw = @elapsed log = evolve_pepskit!(
            psi, 0.1;
            dt = 0.02, order = 2, schedule = :serial,
            maxdim = 3, cutoff = 1e-12, rel_floor = 0.0,
        )
        t = round(t + 0.1; digits = 10)
        say(@sprintf("  evolved to t=%.1f  (%.1f s, chunk trunc=%.2e)", t, tw, log.max_truncerr))
    end
    for _ = 1:2
        evolve_chunk!()
    end
    say("-- sector_density at t=0.2 --")
    res02 = sector_density(psi; verbose = true)
    d02 = abs(res02.n - KRYLOV_6X6[0.2])
    say(@sprintf(
        "6x6 t=0.2: n=%.10f  leak=%.6e  secs=%.1f  |n - Krylov(%.8f)| = %.3e  (expect <= 2e-4) %s",
        res02.n, res02.leak, res02.secs, KRYLOV_6X6[0.2], d02, d02 <= 2e-4 ? "OK" : "EXCEEDED"))
    say(@sprintf("maxrss after t=0.2 eval: %.2f GB", Sys.maxrss() / 2^30))
    GC.gc()
    for _ = 1:14
        evolve_chunk!()
    end
    say("-- sector_density at t=1.6 --")
    res16 = sector_density(psi; verbose = true)
    d16 = abs(res16.n - KRYLOV_6X6[1.6])
    say(@sprintf(
        "6x6 t=1.6: n=%.10f  leak=%.6e  secs=%.1f  |n - Krylov(%.8f)| = %.3e  (~1e-3 expected, no hard threshold)",
        res16.n, res16.leak, res16.secs, KRYLOV_6X6[1.6], d16))
    ram_gb = Sys.maxrss() / 2^30
    say(@sprintf("maxrss after t=1.6 eval: %.2f GB", ram_gb))
    return res02, res16, d02, d16, ram_gb
end

function main()
    say("julia $(VERSION), threads=$(Threads.nthreads()), BLAS threads set to 16")
    LinearAlgebra.BLAS.set_num_threads(16)
    t_start = time()
    pass4, maxdn, leak26, gaps = gate_4x4()
    if !pass4
        say("4x4 GATE FAILED - aborting before 6x6")
        say(@sprintf("RESULT pass=false gate4x4_max_diff=%.6e leak4x4_t26=%.6e", maxdn, leak26))
        exit(1)
    end
    res02, res16, d02, d16, ram_gb = run_6x6()
    pass = pass4 && d02 <= 2e-4
    say("")
    say(@sprintf("n_full - n_sector gaps at 4x4 t=0.5/1.6/2.6: %+.3e %+.3e %+.3e",
        gaps[1], gaps[2], gaps[3]))
    say(@sprintf(
        "RESULT pass=%s gate4x4_max_diff=%.6e leak4x4_t26=%.6e n6_t02=%.10f n6_diff_krylov_t02=%.6e n6_t16=%.10f n6_diff_krylov_t16=%.6e leak6_t02=%.6e leak6_t16=%.6e cost6_s=%.1f ram6_gb=%.2f total_s=%.0f",
        pass, maxdn, leak26, res02.n, d02, res16.n, d16,
        res02.leak, res16.leak, res16.secs, ram_gb, time() - t_start))
end

main()
