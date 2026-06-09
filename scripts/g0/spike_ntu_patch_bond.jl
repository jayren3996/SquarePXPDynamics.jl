# SPIKE NTU PATCH — validate build_ntu_bondenv against the in-process
# build_exact_cluster_bondenv GOLD reference on a FROZEN D=2 state at a RISE time
# (t≈1.8, where the blowup lives), on ONE horizontal center->right bond.
#
# GATES (in order):
#  G-MAXPATCH : patch=:torus MUST equal build_exact_cluster_bondenv to ~1e-12.
#  G-SHAPE    : space == (C^16 (x) C^8)<-(C^16 (x) C^8), flags [0,0,1,1].
#  G-PSD      : min eig of (benv+benv')/2 >= -tol; report cond(Z'Z).
#  G-VALUE    : rel-dev vs exact for P6/P10/P12 (trace & lambda); monotone toward torus.
#  G-COMPOSE  : feed each env through the EXACT spike pipeline; no-op-at-BIG <1e-9.
#  G-ACCURACY : |dn| ladder: bare-SVD < NTU(P6) < ... < exact-cluster.
#  G-COST     : @time + Sys.maxrss of ONE NTU build vs ONE exact build.

using SquarePXPDynamics
using PEPSKit
using TensorKit
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "ntu_bondenv.jl"))   # also brings in spike_n4 builder

# ---- compose a truncated density from an env (the EXACT spike pipeline) ----
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
    T_FREEZE = 1.8     # RISE time where the blowup lives

    GC.gc()
    @printf("RSS at start = %.1f MB ; total mem = %.1f MB\n",
            Sys.maxrss() / 2^20, Sys.total_memory() / 2^20); flush(stdout)

    # === FROZEN-STATE SETUP at the RISE time ===
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
    @printf("center (r,c)=(%d,%d), right leaf (r,c)=(%d,%d) on %dx%d\n",
            rc, cc, rc, cp1, Nr, Nc); flush(stdout)

    # lossless projected state -> the oracle density
    lossless = deepcopy(psi)
    project_star_pepskit!(lossless, c0; dt = 0.02, maxdim = BIG, cutoff = 0.0,
                          rel_floor = 0.0, evolution = :real, projected = true)
    n_lossless = exact_density_finite(lossless; max_sites = 16)
    @printf("n_lossless (oracle) = %.8f\n", n_lossless); flush(stdout)

    A = lossless.peps.A[rc, cc]
    B = lossless.peps.A[rc, cp1]
    X, a, b, Y = PEPSKit._qr_bond(A, B)
    @printf("X q-leg space(X,2) = %s isdual=%s ; Y q-leg space(Y,4) = %s isdual=%s\n",
            string(space(X, 2)), isdual(space(X, 2)),
            string(space(Y, 4)), isdual(space(Y, 4))); flush(stdout)

    a_bt = permute(a, ((1, 2), (3,)))
    b_bt = permute(b, ((1,), (2, 3)))

    # ============================================================
    # G-COST : time ONE exact build vs ONE NTU(P6) build
    # ============================================================
    println("\n===== G-COST ====="); flush(stdout)
    GC.gc(); rss0 = Sys.maxrss() / 2^20
    t_exact = @elapsed benv_exact = build_exact_cluster_bondenv(lossless, X, Y, rc, cc, cp1)
    rss_exact = Sys.maxrss() / 2^20
    @printf("exact-cluster build: %.2f s ; peak RSS now %.1f MB (Δ %.1f MB)\n",
            t_exact, rss_exact, rss_exact - rss0); flush(stdout)

    GC.gc(); rss1 = Sys.maxrss() / 2^20
    t_p6 = @elapsed benv_p6 = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1;
                                                patch = :P6, boundary = :trace)
    rss_p6 = Sys.maxrss() / 2^20
    @printf("NTU P6 (trace) build: %.4f s ; peak RSS now %.1f MB (Δ %.1f MB)\n",
            t_p6, rss_p6, rss_p6 - rss1); flush(stdout)
    @printf("SPEEDUP (exact/P6) = %.1fx wall\n", t_exact / max(t_p6, 1e-9)); flush(stdout)

    # ============================================================
    # G-MAXPATCH : :torus must equal the exact builder to ~1e-12
    # ============================================================
    println("\n===== G-MAXPATCH ====="); flush(stdout)
    benv_torus = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1; patch = :torus, boundary = :trace)
    rel_torus = norm(benv_torus - benv_exact) / norm(benv_exact)
    @printf("rel-dev(:torus vs exact) = %.4e  (gate < 1e-12)\n", rel_torus); flush(stdout)
    gate_maxpatch = rel_torus < 1e-12
    @printf("G-MAXPATCH PASS = %s\n", string(gate_maxpatch)); flush(stdout)

    # ============================================================
    # G-SHAPE
    # ============================================================
    println("\n===== G-SHAPE ====="); flush(stdout)
    target = (ComplexSpace(16) ⊗ ComplexSpace(8)) ← (ComplexSpace(16) ⊗ ComplexSpace(8))
    flags_p6 = [isdual(space(benv_p6, ax)) for ax in 1:numind(benv_p6)]
    shape_ok = space(benv_p6) == target && flags_p6 == [false, false, true, true] &&
               codomain(benv_p6) == domain(benv_p6)
    @printf("space(benv_p6) = %s  flags %s\n", string(space(benv_p6)), string(flags_p6))
    @printf("G-SHAPE PASS = %s (matches exact target [0,0,1,1])\n", string(shape_ok)); flush(stdout)

    # ============================================================
    # G-PSD
    # ============================================================
    println("\n===== G-PSD ====="); flush(stdout)
    H = (benv_p6 + benv_p6') / 2
    evs = real.(eigvals(convert(Array, H)))
    mineig = minimum(evs)
    nrm = norm(benv_p6)
    @printf("min eig (benv+benv')/2 = %.4e ; norm = %.4e ; min/norm = %.4e\n",
            mineig, nrm, mineig / nrm); flush(stdout)
    psd_ok = mineig >= -1e-10 * nrm
    Z = PEPSKit.positive_approx(benv_p6)
    @printf("cond(Z'Z) [NTU P6] = %.4e (exact env was 9.03e20)\n", cond(Z' * Z)); flush(stdout)
    @printf("G-PSD PASS (Gram, no clamp) = %s\n", string(psd_ok)); flush(stdout)

    # ============================================================
    # G-VALUE : rel-dev vs exact for the patch ladder (trace & lambda)
    # ============================================================
    println("\n===== G-VALUE (rel-dev vs exact-cluster) ====="); flush(stdout)
    function reldev(patch, boundary)
        b = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1; patch = patch, boundary = boundary)
        return norm(b - benv_exact) / norm(benv_exact)
    end
    rd_p6_tr = reldev(:P6, :trace)
    rd_p6_la = reldev(:P6, :lambda)
    rd_p10_tr = reldev(:P10, :trace)
    rd_p10_la = reldev(:P10, :lambda)
    rd_p12_tr = reldev(:P12, :trace)
    rd_p12_la = reldev(:P12, :lambda)
    @printf("rel-dev P6  trace = %.4e ; lambda = %.4e\n", rd_p6_tr, rd_p6_la)
    @printf("rel-dev P10 trace = %.4e ; lambda = %.4e\n", rd_p10_tr, rd_p10_la)
    @printf("rel-dev P12 trace = %.4e ; lambda = %.4e\n", rd_p12_tr, rd_p12_la)
    @printf("rel-dev torus     = %.4e\n", rel_torus); flush(stdout)
    monotone = (rd_p12_tr <= rd_p10_tr <= rd_p6_tr) || (rd_p12_tr < rd_p6_tr)
    @printf("G-VALUE monotone-toward-exact (P12<=P10<=P6 trace) = %s\n", string(monotone)); flush(stdout)

    # ============================================================
    # G-COMPOSE + G-ACCURACY : the |dn| ladder
    # ============================================================
    println("\n===== G-COMPOSE + G-ACCURACY ====="); flush(stdout)
    alg_big = ALSTruncation(; trunc = PEPSKit.truncrank(BIG) & PEPSKit.truncerror(; atol = 1e-12),
                            maxiter = 50, tol = 1e-12)
    alg_D = ALSTruncation(; trunc = PEPSKit.truncrank(const_D) & PEPSKit.truncerror(; atol = 1e-12),
                          maxiter = 50, tol = 1e-12)

    # no-op-at-BIG (plumbing) with the NTU P6 env
    n_noop, _ = compose_density(lossless, X, Y, a_bt, b_bt, benv_p6, rc, cc, cp1, alg_big)
    noop_dev = abs(n_noop - n_lossless)
    @printf("no-op-at-BIG (NTU P6): n=%.8f ; dev=%.4e (plumbing <1e-9)\n",
            n_noop, noop_dev); flush(stdout)

    function dn_of(benv; label = "")
        try
            n, fid = compose_density(lossless, X, Y, a_bt, b_bt, benv, rc, cc, cp1, alg_D)
            dn = abs(n - n_lossless)
            @printf("  %-22s n_D=%.8f  |dn|=%.4e  fid=%.6e\n", label, n, dn, fid); flush(stdout)
            return dn
        catch err
            @printf("  %-22s THREW: %s\n", label, string(err)); flush(stdout)
            return NaN
        end
    end

    # bare-SVD (no env)
    ab = PEPSKit._combine_ab(a_bt, b_bt)
    a0, s0, b0 = PEPSKit.svd_trunc(permute(ab, ((1, 3), (4, 2))); trunc = PEPSKit.truncrank(const_D))
    a0f, b0f = PEPSKit.absorb_s(a0, s0, b0)
    A0, B0 = PEPSKit._qr_bond_undo(X, permute(a0f, ((1, 2), (3,))),
                                   permute(b0f, ((1,), (2, 3))), Y)
    bare = deepcopy(lossless)
    bare.peps.A[rc, cc] = A0
    bare.peps.A[rc, cp1] = B0
    n_bare = exact_density_finite(bare; max_sites = 16)
    dn_bare = abs(n_bare - n_lossless)
    @printf("  %-22s n_D=%.8f  |dn|=%.4e\n", "bare-SVD (no env)", n_bare, dn_bare); flush(stdout)

    dn_exact = dn_of(benv_exact; label = "exact-cluster")
    dn_p6_tr = dn_of(benv_p6; label = "NTU P6 trace")
    benv_p6_la = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1; patch = :P6, boundary = :lambda)
    dn_p6_la = dn_of(benv_p6_la; label = "NTU P6 lambda")
    benv_p10_tr = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1; patch = :P10, boundary = :trace)
    dn_p10_tr = dn_of(benv_p10_tr; label = "NTU P10 trace")
    benv_p10_la = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1; patch = :P10, boundary = :lambda)
    dn_p10_la = dn_of(benv_p10_la; label = "NTU P10 lambda")
    benv_p12_tr = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1; patch = :P12, boundary = :trace)
    dn_p12_tr = dn_of(benv_p12_tr; label = "NTU P12 trace")
    dn_torus = dn_of(benv_torus; label = "NTU torus (==exact)")

    # LOOP-CONTENT control: drop one plaquette row -> tree -> ~bare-SVD
    cs_drop = patch_coords_droprow(rc, cc, cp1, Nr, Nc)
    benv_drop = build_ntu_bondenv(lossless, X, Y, rc, cc, cp1; coords = cs_drop, boundary = :trace)
    dn_drop = dn_of(benv_drop; label = "DROP-1-loop (4 cells)")

    # ============================================================
    # SUMMARY
    # ============================================================
    println("\n===== NTU SUMMARY (t=1.8 rise, bond (4,1)->(4,2)) =====")
    @printf("n_lossless              = %.8f\n", n_lossless)
    @printf("|dn| bare-SVD           = %.4e   (MIN BAR to beat)\n", dn_bare)
    @printf("|dn| NTU P6 trace       = %.4e\n", dn_p6_tr)
    @printf("|dn| NTU P6 lambda      = %.4e\n", dn_p6_la)
    @printf("|dn| NTU P10 trace      = %.4e\n", dn_p10_tr)
    @printf("|dn| NTU P10 lambda     = %.4e\n", dn_p10_la)
    @printf("|dn| NTU P12 trace      = %.4e\n", dn_p12_tr)
    @printf("|dn| NTU torus          = %.4e\n", dn_torus)
    @printf("|dn| exact-cluster      = %.4e   (TARGET)\n", dn_exact)
    @printf("|dn| DROP-1-loop control= %.4e   (~bare-SVD if loop drives gain)\n", dn_drop)
    @printf("\nGATES: MAXPATCH=%s SHAPE=%s PSD=%s VALUE-monotone=%s\n",
            string(gate_maxpatch), string(shape_ok), string(psd_ok), string(monotone))
    beats_bare = isfinite(dn_p6_tr) && dn_p6_tr < dn_bare
    @printf("NTU P6 trace beats bare-SVD = %s\n", string(beats_bare))
    @printf("COST: exact %.2fs / NTU P6 %.4fs (%.0fx)\n", t_exact, t_p6, t_exact / max(t_p6, 1e-9))
    @printf("RSS peak = %.1f MB\n", Sys.maxrss() / 2^20)
    println("DONE")
    flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
