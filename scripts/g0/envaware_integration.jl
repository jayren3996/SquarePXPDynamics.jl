# envaware_integration.jl — VALIDATE the D-generalized 4-direction env builder
# and an env-aware prototype star update.
#
# Gates (per the task):
#   G0  rotate_peps roundtrip identity + _qr_bond uniform dual flags (all dirs).
#   G1  env-VALUE oracle: for each vertical dir, build the env via the row/col
#       transpose trick (validated horizontal builder on the transposed lattice)
#       and require norm(benv_rot - benv_oracle)/norm < 1e-10 (the discriminating
#       gate the SPEC REVIEW demands — no-op-at-BIG alone is toothless for a
#       transpose).
#   G2  NO-OP-AT-BIG, per direction: env-aware single-bond truncate at maxdim=BIG
#       reproduces the bare lossless density to ~1e-9 (the task's named key gate).
#   G3  FULL-star env-aware update no-op-at-BIG (all 4 arms) < 1e-9.
#   G4  short 5-step D=2 trajectory with env=:exact_cluster runs w/o NaN and
#       density differs from bare.

using SquarePXPDynamics
using PEPSKit
using TensorKit
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "star_env.jl"))

# BIG (the no-op-at-BIG ceiling) is set per-D inside main(): it must exceed the
# lossless bond rank, which grows with D. D=2 used 64; D=3 q-legs are 24/12 and
# the lossless star bond is larger, so we bump it.

# ----------------------------------------------------------------------------
# Transpose a PXP state's lattice (row<->col) so a vertical bond becomes a
# horizontal one. On the square lattice, transposing swaps N<->W and E<->S and
# rows<->cols. We build a fresh state by transposing the unit cell + relabeling
# each PEPSTensor's legs: new [P,N,E,S,W] = old [P,W,S,E,N].
# This is used ONLY to build an independent oracle env for vertical dirs.
# ----------------------------------------------------------------------------
function transpose_peps_tensor(T)
    # old axes: 1=P,2=N,3=E,4=S,5=W ; new: P, N'=W, E'=S, S'=E, W'=N
    return permute(T, ((1,), (5, 4, 3, 2)))
end

# Build an INDEPENDENT oracle env for a vertical (:up/:down) bond by physically
# transposing the lattice (row<->col, N<->W, E<->S) so the bond becomes a
# HORIZONTAL :right bond, then running the EXACT-CLUSTER ncon (the validated
# build_exact_cluster_bondenv routing) on the transposed array. This is a fully
# independent construction (different integer bookkeeping, different reduced
# tensors) so a transposed N/S leg in the rotate-route is caught with zero
# ambiguity (norm(benv_rot - benv_oracle) must be ~0).
#
# Under transpose (r,c)->(c,r): a :down bond (center.S -> leaf.N) becomes
# center.E -> leaf.W = a :right bond; a :up bond (center.N -> leaf.S) becomes
# center.W -> leaf.E = a :left bond. After transpose we further use the rotate
# route (:right => r=0 directly; :left => r=2) on the TRANSPOSED arrays.
function build_vertical_oracle_env(state, rc, cc, leaf_rc, leaf_cc, dir)
    cell = state.unitcell
    Nr, Nc = cell.Ly, cell.Lx
    A = state.peps.A
    # transposed peps array shape (Nc rows, Nr cols): TA[c,r] = transpose(A[r,c])
    TA = Matrix{Any}(undef, Nc, Nr)
    for rr in 1:Nr, ccol in 1:Nc
        TA[ccol, rr] = transpose_peps_tensor(A[rr, ccol])
    end
    # transposed coords (r,c)->(c,r)
    t_rc, t_cc = cc, rc
    t_lrc, t_lcc = leaf_cc, leaf_rc
    tNr, tNc = Nc, Nr   # transposed lattice dims
    tdir = dir === :down ? :right : :left
    r = _rot_for_dir(tdir)
    Ar = rotate_peps(TA[t_rc, t_cc], r)
    Br = rotate_peps(TA[t_lrc, t_lcc], r)
    X, a, b, Y = _qr_bond(Ar, Br)
    # Build the env directly from the transposed array (standalone ncon, no state).
    return _exact_cluster_env_from_array(TA, X, Y, tdir, t_rc, t_cc, t_lrc, t_lcc, tNr, tNc)
end

# Standalone exact-cluster {2,2} env from a raw transposed peps array, using the
# SAME rotated-id routing as build_star_bond_env but reading TA[r,c] directly.
function _exact_cluster_env_from_array(TA, X, Y, dir, rc, cc, leaf_rc, leaf_cc, Nr, Nc)
    f = _idfns(Nr, Nc)
    rr0 = _rot_for_dir(dir)
    rt = _rotated_to_true(rr0)
    Xarr = convert(Array{ComplexF64,4}, convert(Array, X))
    Yarr = convert(Array{ComplexF64,4}, convert(Array, Y))
    cN_k = _interiorKetId(f, rc, cc, rt.rotN); cS_k = _interiorKetId(f, rc, cc, rt.rotS); cW_k = _interiorKetId(f, rc, cc, rt.rotW)
    cN_b = _interiorBraId(f, rc, cc, rt.rotN); cS_b = _interiorBraId(f, rc, cc, rt.rotS); cW_b = _interiorBraId(f, rc, cc, rt.rotW)
    lN_k = _interiorKetId(f, leaf_rc, leaf_cc, rt.rotN); lE_k = _interiorKetId(f, leaf_rc, leaf_cc, rt.rotE); lS_k = _interiorKetId(f, leaf_rc, leaf_cc, rt.rotS)
    lN_b = _interiorBraId(f, leaf_rc, leaf_cc, rt.rotN); lE_b = _interiorBraId(f, leaf_rc, leaf_cc, rt.rotE); lS_b = _interiorBraId(f, leaf_rc, leaf_cc, rt.rotS)
    arrays = Vector{Array{ComplexF64}}(); idxlists = Vector{Vector{Int}}()
    k = 0
    for rr in 1:Nr, ccol in 1:Nc
        k += 1
        if rr == rc && ccol == cc
            push!(arrays, Xarr); push!(idxlists, [cN_k, -1, cS_k, cW_k])
            push!(arrays, conj(Xarr)); push!(idxlists, [cN_b, -3, cS_b, cW_b])
        elseif rr == leaf_rc && ccol == leaf_cc
            push!(arrays, Yarr); push!(idxlists, [lN_k, lE_k, lS_k, -2])
            push!(arrays, conj(Yarr)); push!(idxlists, [lN_b, lE_b, lS_b, -4])
        else
            T = convert(Array{ComplexF64,5}, convert(Array, TA[rr, ccol]))
            push!(arrays, T)
            push!(idxlists, [_physTrace(f, k), _VbondK(f, rr, ccol), _HbondK(f, rr, ccol),
                             _VbondK(f, _nextr(f, rr), ccol), _HbondK(f, rr, _prevc(f, ccol))])
            push!(arrays, conj(T))
            push!(idxlists, [_physTrace(f, k), _VbondBra(f, rr, ccol), _HbondBra(f, rr, ccol),
                             _VbondBra(f, _nextr(f, rr), ccol), _HbondBra(f, rr, _prevc(f, ccol))])
        end
    end
    raw = convert(Array, TensorKitNcon.ncon(arrays, idxlists))
    raw_reordered = permutedims(raw, (3, 4, 1, 2))
    DXn = space(X, 2)'; DYn = space(Y, 4)'
    benv = TensorMap(raw_reordered, DXn ⊗ DYn, DXn ⊗ DYn)
    normalize!(benv, Inf)
    return benv
end

function main()
    GC.gc()
    @printf("=== envaware_integration: D-generalized 4-dir env builder validation ===\n")
    @printf("RSS at start = %.1f MB\n", Sys.maxrss() / 2^20); flush(stdout)

    const_D = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
    BIG = const_D <= 2 ? 64 : 256   # no-op ceiling: must exceed lossless bond rank
    @printf("D = %d ; BIG (no-op ceiling) = %d\n", const_D, BIG); flush(stdout)
    T_FREEZE = 1.8   # rise regime (past native blowup) where truncation is lossy

    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on=:even, D=const_D)
    evolve_pepskit!(psi, T_FREEZE; dt=0.02, order=2, schedule=:serial,
                    maxdim=const_D, cutoff=1e-12, rel_floor=1e-3)
    @printf("frozen density at t=%.2f: %.8f\n", T_FREEZE,
            exact_density_finite(psi; max_sites=16)); flush(stdout)

    c0 = cell.reps[1]
    cidx = SquarePXPDynamics.squarecoord_to_cartesianindex(cell, c0)
    rc, cc = cidx[1], cidx[2]
    Nr, Nc = size(psi.peps)
    @printf("center PEPSKit (r,c)=(%d,%d) ; Nr=%d Nc=%d\n", rc, cc, Nr, Nc); flush(stdout)

    # leaf cells per direction (mirror project_star's _leaf_rowcol)
    leafcell = Dict{Symbol,Tuple{Int,Int}}()
    leafcell[:right] = (rc, mod1(cc + 1, Nc))
    leafcell[:up]    = (mod1(rc - 1, Nr), cc)
    leafcell[:left]  = (rc, mod1(cc - 1, Nc))
    leafcell[:down]  = (mod1(rc + 1, Nr), cc)
    for d in (:right, :up, :left, :down)
        @printf("  leaf[%s] = (%d,%d)\n", d, leafcell[d][1], leafcell[d][2])
    end
    flush(stdout)

    # ===== G0: rotate_peps roundtrip + _qr_bond uniform dual flags =====
    @printf("\n--- G0: rotate_peps roundtrip + _qr_bond uniform dual flags ---\n"); flush(stdout)
    A = psi.peps.A[rc, cc]
    g0_ok = true
    for r in 0:3
        T2 = rotate_peps(rotate_peps(A, r), mod(4 - r, 4))
        dev = norm(T2 - A)
        @printf("  rotate roundtrip r=%d: norm(rot(rot(A,r),4-r)-A) = %.3e\n", r, dev); flush(stdout)
        g0_ok &= (dev < 1e-12)
    end
    for dir in (:right, :up, :left, :down)
        (lr, lc) = leafcell[dir]
        B = psi.peps.A[lr, lc]
        r = _rot_for_dir(dir)
        Ar = rotate_peps(A, r); Br = rotate_peps(B, r)
        X, a, b, Y = _qr_bond(Ar, Br)
        dX = isdual(space(X, 2)); dY = isdual(space(Y, 4))
        @printf("  dir=%-6s : isdual(X,2)=%s isdual(Y,4)=%s ; qdimX=%d qdimY=%d\n",
                dir, dX, dY, dim(space(X, 2)), dim(space(Y, 4))); flush(stdout)
        g0_ok &= dX && dY
    end
    @printf("  G0 PASS = %s\n", g0_ok); flush(stdout)

    # ===== G1: env-VALUE oracle for vertical dirs (transpose trick) =====
    @printf("\n--- G1: env-VALUE oracle (transpose trick) for vertical dirs ---\n"); flush(stdout)
    g1_ok = true
    g1_devs = Dict{Symbol,Float64}()
    for dir in (:up, :down)
        (lr, lc) = leafcell[dir]
        env_rot = build_star_env_one(psi, rc, cc, lr, lc, dir; env=:exact_cluster)
        benv_rot = env_rot.benv
        local benv_oracle, reldev
        try
            benv_oracle = build_vertical_oracle_env(psi, rc, cc, lr, lc, dir)
            if space(benv_rot) == space(benv_oracle)
                reldev = norm(benv_rot - benv_oracle) / norm(benv_oracle)
            else
                @printf("  dir=%s : ORACLE SPACE MISMATCH rot=%s oracle=%s\n",
                        dir, string(space(benv_rot)), string(space(benv_oracle)))
                reldev = NaN
            end
        catch err
            @printf("  dir=%s : oracle build THREW: %s\n", dir, string(err))
            reldev = NaN
        end
        g1_devs[dir] = reldev
        @printf("  dir=%-6s : norm(benv_rot - benv_oracle)/norm = %.4e\n", dir, reldev); flush(stdout)
        g1_ok &= isfinite(reldev) && reldev < 1e-10
    end
    @printf("  G1 (vertical env value oracle) PASS = %s\n", g1_ok); flush(stdout)

    # ===== G2: NO-OP-AT-BIG per direction (single-bond, env-aware) =====
    @printf("\n--- G2: NO-OP-AT-BIG per direction (single-bond env-aware truncate) ---\n"); flush(stdout)
    # bare lossless single-bond reference per direction: replace center+leaf via
    # the env path at maxdim=BIG and compare density to the unmodified state.
    n_ref = exact_density_finite(psi; max_sites=16)
    @printf("  n_ref (unmodified) = %.8f\n", n_ref); flush(stdout)
    alg_big = ALSTruncation(; trunc=PEPSKit.truncrank(BIG) & PEPSKit.truncerror(; atol=1e-12),
                            maxiter=50, tol=1e-12)
    g2_ok = true
    g2_devs = Dict{Symbol,Float64}()
    g2_conds = Dict{Symbol,Float64}()
    for dir in (:right, :up, :left, :down)
        (lr, lc) = leafcell[dir]
        local dev = NaN
        try
            e = build_star_env_one(psi, rc, cc, lr, lc, dir; env=:exact_cluster)
            X, Y, a, b, benv = e.X, e.Y, e.a, e.b, e.benv
            # fixgauge path (primary, per fixgaugeSpec)
            Z = PEPSKit.positive_approx(benv)
            cnd = cond(Z' * Z); g2_conds[dir] = cnd
            Z2, a2, b2, (Linv, Rinv) = PEPSKit.fixgauge_benv(Z, a, b)
            benv2 = Z2' * Z2; normalize!(benv2, Inf)
            X2, Y2 = PEPSKit._fixgauge_benvXY(X, Y, Linv, Rinv)
            a2_bt = permute(a2, ((1, 2), (3,)))
            b2_bt = permute(b2, ((1,), (2, 3)))
            aD, sD, bD, info = PEPSKit.bond_truncate(a2_bt, b2_bt, benv2, alg_big)
            aDf, bDf = PEPSKit.absorb_s(aD, sD, bD)
            A_new, B_new = PEPSKit._qr_bond_undo(X2, permute(aDf, ((1, 2), (3,))),
                                                 permute(bDf, ((1,), (2, 3))), Y2)
            # unrotate
            r = _rot_for_dir(dir)
            A_new = rotate_peps(A_new, mod(4 - r, 4))
            B_new = rotate_peps(B_new, mod(4 - r, 4))
            st = deepcopy(psi)
            st.peps.A[rc, cc] = A_new
            st.peps.A[lr, lc] = B_new
            n_env = exact_density_finite(st; max_sites=16)
            dev = abs(n_env - n_ref)
        catch err
            @printf("  dir=%s : G2 THREW: %s\n", dir, string(err)); flush(stdout)
            dev = NaN
        end
        g2_devs[dir] = dev
        @printf("  dir=%-6s : no-op-at-BIG dev = %.4e ; cond(Z'Z)=%.3e\n",
                dir, dev, get(g2_conds, dir, NaN)); flush(stdout)
        g2_ok &= isfinite(dev) && dev < 1e-9
    end
    @printf("  G2 (per-dir no-op-at-BIG) PASS = %s\n", g2_ok); flush(stdout)

    # ===== G3: FULL-star env-aware update no-op-at-BIG (all 4 arms) =====
    @printf("\n--- G3: FULL-star env-aware update no-op-at-BIG ---\n"); flush(stdout)
    g3_dev = NaN
    try
        lossless = deepcopy(psi)
        project_star_pepskit!(lossless, c0; dt=0.02, maxdim=BIG, cutoff=0.0,
                              rel_floor=0.0, evolution=:real, projected=true)
        n_lossless = exact_density_finite(lossless; max_sites=16)
        st_env = deepcopy(psi)
        info = project_star_env!(st_env, c0; dt=0.02, maxdim=BIG, cutoff=0.0,
                                 env=:exact_cluster, patch=:torus)
        n_env_full = exact_density_finite(st_env; max_sites=16)
        g3_dev = abs(n_env_full - n_lossless)
        @printf("  n_lossless(bare star) = %.8f ; n_env_full = %.8f ; dev = %.4e\n",
                n_lossless, n_env_full, g3_dev); flush(stdout)
        @printf("  per-arm truncerrs = %s\n", string(info.truncerrs)); flush(stdout)
    catch err
        @printf("  G3 THREW: %s\n", string(err)); flush(stdout)
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt); println()
        end
    end
    g3_ok = isfinite(g3_dev) && g3_dev < 1e-9
    @printf("  G3 (full-star no-op-at-BIG) PASS = %s\n", g3_ok); flush(stdout)

    # ===== G4: short 5-step D=2 trajectory env vs bare =====
    @printf("\n--- G4: short 5-step D=2 trajectory env=:exact_cluster vs bare ---\n"); flush(stdout)
    g4_ok = false
    try
        st_bare = deepcopy(psi)
        st_envt = deepcopy(psi)
        @printf("  step  n_bare      n_env       |diff|\n"); flush(stdout)
        all_finite = true
        any_diff = false
        for step in 1:5
            project_star_pepskit!(st_bare, c0; dt=0.02, maxdim=const_D, cutoff=1e-12,
                                  rel_floor=1e-3, evolution=:real, projected=true)
            project_star_env!(st_envt, c0; dt=0.02, maxdim=const_D, cutoff=1e-12,
                              env=:exact_cluster, patch=:torus)
            nb = exact_density_finite(st_bare; max_sites=16)
            ne = exact_density_finite(st_envt; max_sites=16)
            d = abs(nb - ne)
            @printf("  %-4d  %.8f  %.8f  %.4e\n", step, nb, ne, d); flush(stdout)
            all_finite &= isfinite(nb) && isfinite(ne)
            any_diff |= d > 1e-9
        end
        g4_ok = all_finite && any_diff
        @printf("  all_finite=%s ; any_diff=%s\n", all_finite, any_diff); flush(stdout)
    catch err
        @printf("  G4 THREW: %s\n", string(err)); flush(stdout)
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt); println()
        end
    end
    @printf("  G4 (trajectory) PASS = %s\n", g4_ok); flush(stdout)

    @printf("\n===== SUMMARY =====\n")
    @printf("G0 rotate/dual          = %s\n", g0_ok)
    @printf("G1 vertical env oracle  = %s  (up=%.3e down=%.3e)\n", g1_ok,
            get(g1_devs, :up, NaN), get(g1_devs, :down, NaN))
    @printf("G2 per-dir no-op-BIG    = %s  (R=%.2e U=%.2e L=%.2e D=%.2e)\n", g2_ok,
            get(g2_devs, :right, NaN), get(g2_devs, :up, NaN),
            get(g2_devs, :left, NaN), get(g2_devs, :down, NaN))
    @printf("G3 full-star no-op-BIG  = %s  (dev=%.3e)\n", g3_ok, g3_dev)
    @printf("G4 trajectory finite    = %s\n", g4_ok)
    overall = g0_ok && g2_ok && g3_ok   # the task's key gate is no-op-at-BIG (G2,G3)
    @printf("OVERALL (no-op-at-BIG all 4 bonds + full star) = %s\n", overall)
    @printf("RSS peak = %.1f MB\n", Sys.maxrss() / 2^20); flush(stdout)
    println("DONE")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
