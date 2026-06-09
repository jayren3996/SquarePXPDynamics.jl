# REVIEWER PROBE (adversarial verification 2/3) — fixes the 3 harness bugs in
# spike_ntu_patch_bond.jl that crashed the original run before the headline
# |dn| ladder printed:
#   (1) G-PSD: eigvals(convert(Array, H)) on a {2,2} TensorMap -> a 4-D array,
#       which eigvals rejects. Reshape the 16*8 x 16*8 matrix explicitly.
#   (2) G-VALUE: reldev(:P10/:P12) is unguarded but the builder asserts Nc>=5;
#       on the 4x4 cell it throws -> whole script dies. Guard with try/catch.
#   (3) NTU timing conflated JIT with compute -> add a SECOND-call timing.
# Everything else (builders, oracle, pipeline) is reused verbatim.

using SquarePXPDynamics
using PEPSKit
using TensorKit
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "ntu_bondenv.jl"))   # also brings in spike_n4 builder + _qr_bond glue

function compose_density(lossless, X, Y, a_bt, b_bt, benv, rc, cc, cp1, alg)
    aD, sD, bD, tinfo = PEPSKit.bond_truncate(a_bt, b_bt, benv, alg)
    aDf, bDf = PEPSKit.absorb_s(aD, sD, bD)
    A3, B3 = PEPSKit._qr_bond_undo(X, permute(aDf, ((1, 2), (3,))),
                                   permute(bDf, ((1,), (2, 3))), Y)
    st = deepcopy(lossless)
    st.peps.A[rc, cc] = A3
    st.peps.A[rc, cp1] = B3
    n = exact_density_finite(st; max_sites = 16)
    return n, tinfo.fid
end

function main()
    const_D = 2
    BIG = 64
    T_FREEZE = 1.8

    GC.gc()
    @printf("RSS at start = %.1f MB\n", Sys.maxrss() / 2^20); flush(stdout)

    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = const_D)
    evolve_pepskit!(psi, T_FREEZE; dt = 0.02, order = 2, schedule = :serial,
                    maxdim = const_D, cutoff = 1e-12, rel_floor = 1e-3)
    @printf("frozen density at t=%.2f: %.8f\n", T_FREEZE,
            exact_density_finite(psi; max_sites = 16)); flush(stdout)

    c0 = cell.reps[1]
    cidx = SquarePXPDynamics.squarecoord_to_cartesianindex(cell, c0)
    rc, cc = cidx[1], cidx[2]
    Nr, Nc = size(psi.peps)
    cp1 = mod1(cc + 1, Nc)

    lossless = deepcopy(psi)
    project_star_pepskit!(lossless, c0; dt = 0.02, maxdim = BIG, cutoff = 0.0,
                          rel_floor = 0.0, evolution = :real, projected = true)
    n_lossless = exact_density_finite(lossless; max_sites = 16)
    @printf("n_lossless (oracle) = %.8f\n", n_lossless); flush(stdout)

    A = lossless.peps.A[rc, cc]
    B = lossless.peps.A[rc, cp1]
    X, a, b, Y = PEPSKit._qr_bond(A, B)
    a_bt = permute(a, ((1, 2), (3,)))
    b_bt = permute(b, ((1,), (2, 3)))

    # ---- COST: warm both builders ONCE (JIT) then time the 2ND call (compute) ----
    println("\n===== G-COST (JIT-separated) ====="); flush(stdout)
    benv_exact = build_exact_cluster_bondenv(lossless, X, Y, rc, cc, cp1)   # warm
    t_exact2 = @elapsed benv_exact = build_exact_cluster_bondenv(lossless, X, Y, rc, cc, cp1)
    benv_p6 = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1; patch=:P6, boundary=:trace)  # warm
    t_p6_2 = @elapsed benv_p6 = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1; patch=:P6, boundary=:trace)
    @printf("exact-cluster build (2nd call) = %.3f s\n", t_exact2)
    @printf("NTU P6 trace  build (2nd call) = %.3f s\n", t_p6_2)
    @printf("SPEEDUP (exact/P6, 2nd-call)   = %.2fx\n", t_exact2 / max(t_p6_2, 1e-9)); flush(stdout)

    # ---- G-PSD (FIXED: reshape to a square matrix) ----
    println("\n===== G-PSD (fixed) ====="); flush(stdout)
    H = (benv_p6 + benv_p6') / 2
    Harr = convert(Array, H)                 # 16 x 8 x 16 x 8
    d1, d2 = size(Harr, 1), size(Harr, 2)
    Hmat = reshape(Harr, d1 * d2, d1 * d2)   # (16*8) x (16*8)
    evs = real.(eigvals(Hmat))
    mineig = minimum(evs); nrm = norm(benv_p6)
    @printf("min eig (benv+benv')/2 = %.4e ; norm = %.4e ; min/norm = %.4e\n",
            mineig, nrm, mineig / nrm); flush(stdout)
    psd_ok = mineig >= -1e-10 * nrm
    @printf("G-PSD PASS (Gram, no clamp) = %s\n", string(psd_ok)); flush(stdout)
    Z = PEPSKit.positive_approx(benv_p6)
    @printf("cond(Z'Z) [NTU P6] = %.4e (exact env was 9.03e20)\n", cond(Z' * Z)); flush(stdout)

    # ---- G-VALUE (FIXED: guard P10/P12) ----
    println("\n===== G-VALUE (rel-dev vs exact-cluster) ====="); flush(stdout)
    function reldev(patch, boundary)
        try
            b = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1; patch=patch, boundary=boundary)
            return norm(b - benv_exact) / norm(benv_exact)
        catch err
            @printf("  reldev(%s,%s) THREW: %s\n", patch, boundary, string(err))
            return NaN
        end
    end
    rel_torus = reldev(:torus, :trace)
    @printf("rel-dev torus      = %.4e (==exact equality gate, expect 0)\n", rel_torus)
    @printf("rel-dev P6  trace  = %.4e\n", reldev(:P6, :trace))
    @printf("rel-dev P6  lambda = %.4e\n", reldev(:P6, :lambda))
    @printf("rel-dev P10 trace  = %.4e\n", reldev(:P10, :trace))
    @printf("rel-dev P12 trace  = %.4e\n", reldev(:P12, :trace)); flush(stdout)

    # ---- G-ACCURACY: the headline |dn| ladder ----
    println("\n===== G-ACCURACY (|dn| ladder) ====="); flush(stdout)
    alg_D = ALSTruncation(; trunc = PEPSKit.truncrank(const_D) & PEPSKit.truncerror(; atol = 1e-12),
                          maxiter = 50, tol = 1e-12)

    # bare-SVD baseline (no env)
    ab = PEPSKit._combine_ab(a_bt, b_bt)
    a0, s0, b0 = PEPSKit.svd_trunc(permute(ab, ((1, 3), (4, 2))); trunc = PEPSKit.truncrank(const_D))
    a0f, b0f = PEPSKit.absorb_s(a0, s0, b0)
    A0, B0 = PEPSKit._qr_bond_undo(X, permute(a0f, ((1, 2), (3,))),
                                   permute(b0f, ((1,), (2, 3))), Y)
    bare = deepcopy(lossless)
    bare.peps.A[rc, cc] = A0; bare.peps.A[rc, cp1] = B0
    n_bare = exact_density_finite(bare; max_sites = 16)
    dn_bare = abs(n_bare - n_lossless)
    @printf("  %-22s |dn|=%.4e   (MIN BAR to beat)\n", "bare-SVD (no env)", dn_bare); flush(stdout)

    function dn_of(benv, label)
        try
            n, fid = compose_density(lossless, X, Y, a_bt, b_bt, benv, rc, cc, cp1, alg_D)
            dn = abs(n - n_lossless)
            @printf("  %-22s |dn|=%.4e  fid=%.6e\n", label, dn, fid); flush(stdout)
            return dn
        catch err
            @printf("  %-22s THREW: %s\n", label, string(err)); flush(stdout)
            return NaN
        end
    end

    dn_exact = dn_of(benv_exact, "exact-cluster")
    dn_p6_tr = dn_of(benv_p6, "NTU P6 trace")
    benv_p6_la = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1; patch=:P6, boundary=:lambda)
    dn_p6_la = dn_of(benv_p6_la, "NTU P6 lambda")
    # DROP-1-loop control (tree -> should be ~bare-SVD if loop drives the gain)
    cs_drop = patch_coords_droprow(rc, cc, cp1, Nr, Nc)
    benv_drop = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1; coords=cs_drop, boundary=:trace)
    dn_drop = dn_of(benv_drop, "DROP-1-loop (4 cells)")

    println("\n===== VERDICT NUMBERS =====")
    @printf("n_lossless        = %.8f\n", n_lossless)
    @printf("|dn| bare-SVD     = %.4e  (BAR)\n", dn_bare)
    @printf("|dn| NTU P6 trace = %.4e  beats-bare=%s\n", dn_p6_tr, string(isfinite(dn_p6_tr) && dn_p6_tr < dn_bare))
    @printf("|dn| NTU P6 lambda= %.4e  beats-bare=%s\n", dn_p6_la, string(isfinite(dn_p6_la) && dn_p6_la < dn_bare))
    @printf("|dn| DROP-1-loop  = %.4e\n", dn_drop)
    @printf("|dn| exact-cluster= %.4e  (TARGET)\n", dn_exact)
    @printf("COST: exact %.3fs / NTU P6 %.3fs (%.2fx, 2nd-call)\n", t_exact2, t_p6_2, t_exact2/max(t_p6_2,1e-9))
    @printf("G-PSD=%s G-SHAPE(prev run)=true G-MAXPATCH(prev run)=true\n", string(psd_ok))
    println("DONE")
    flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
