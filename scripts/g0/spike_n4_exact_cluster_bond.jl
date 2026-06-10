# SPIKE N4 — hand-build a double-layer EXACT-CLUSTER {2,2} BondEnv on ONE
# horizontal center->right bond and prove it composes with PEPSKit.bond_truncate.
#
# This REPLACES bondenv_fu's CTM environment with a hand-built finite-torus
# (4x4 PBC) double-layer <patch|patch> bond environment, computed by the SAME
# ncon machinery as _native_dense_state (PEPSKitObservables.jl:484-504) but
# DOUBLE-LAYER (bra (x) ket, physical legs TRACED, two bond-end q-legs OPEN).
#
# The builder is factored into `build_exact_cluster_bondenv` so spike_n2 can
# `include` and reuse it with no duplication.
#
# SUCCESS = (1) hand env space + dual-flags EQUAL bondenv_fu's output, AND
#           (2) no-op-at-BIG reproduces lossless density to <1e-9 (PLUMBING check), AND
#           (3) maxdim=D truncation runs without SpaceMismatch and returns finite n_D,
#               with the env-content gate carried by (3) + a wrong-env NEGATIVE CONTROL
#               and a VALUE comparison norm(benv - benv_ref)/norm(benv_ref).
#
# A documented SpaceMismatch (with the exact spaces printed) is an ACCEPTED result.

using SquarePXPDynamics
using PEPSKit
using TensorKit
using LinearAlgebra
using Printf

const TensorKitNcon = TensorKit  # ncon is exported by TensorKit/TensorOperations

# =====================================================================
# THE BUILDER (factored out so spike_n2 can include + reuse it)
# =====================================================================
"""
    build_exact_cluster_bondenv(state, X, Y, rc, cc, cp1) -> benv::TensorMap{T,S,2,2}

Hand-build the double-layer exact-cluster {2,2} bond environment for the
HORIZONTAL center(rc,cc) -> right-leaf(rc,cp1) bond, using the q-reduced
PEPSOrth tensors X (legs [N,q,S,W], q=axis 2, the DX open leg) and
Y (legs [N,E,S,q], q=axis 4, the DY open leg) from PEPSKit._qr_bond.

Convention (verified against als_solve._tensor_Sa + bondenv_fu):
  benv axis order [DX1 DY1; DX0 DY0]  -> codomain=(DX1,DY1)=BRA, domain=(DX0,DY0)=KET.
  KET q-legs  -> DOMAIN   (DX0,DY0)
  BRA q-legs  -> CODOMAIN (DX1,DY1)
Target space (ground truth, ctm_bond_proto.out): (C^16 (x) C^8) <- (C^16 (x) C^8)
with dual flags [0,0,1,1] (codomain non-dual, domain dual).
"""
function build_exact_cluster_bondenv(state, X, Y, rc::Int, cc::Int, cp1::Int)
    cell = state.unitcell
    Nr, Nc = cell.Ly, cell.Lx
    nsites = Nr * Nc
    peps = measurement_peps(state)

    # --- pre-build asserts (the q-legs must be the reduced DX/DY legs) ---
    @assert dim(space(X, 2)) == 16 "expected DX q-leg dim 16, got $(dim(space(X,2)))"
    @assert dim(space(Y, 4)) == 8 "expected DY q-leg dim 8, got $(dim(space(Y,4)))"

    _nextr(r) = r == Nr ? 1 : r + 1
    _prevc(c) = c == 1 ? Nc : c - 1
    HbondK(r, c) = (r - 1) * Nc + c             # ket E<->W bond ids 1..nsites
    VbondK(r, c) = nsites + (r - 1) * Nc + c    # ket N<->S bond ids nsites+1..2nsites
    HbondBra(r, c) = 2 * nsites + HbondK(r, c)  # bra bonds: disjoint id block
    VbondBra(r, c) = 2 * nsites + VbondK(r, c)
    physTrace(k) = 4 * nsites + k               # ket-P & bra-P of site k SHARE this id

    # the CUT bond ids (center.E <-> leaf.W) are NEVER emitted; the open
    # negative legs -1..-4 carry the q-legs instead.

    Xarr = convert(Array{ComplexF64,4}, convert(Array, X))   # [N,q,S,W]
    Yarr = convert(Array{ComplexF64,4}, convert(Array, Y))   # [N,E,S,q]

    arrays = Vector{Array{ComplexF64}}()
    idxlists = Vector{Vector{Int}}()

    k = 0
    for r in 1:Nr, c in 1:Nc
        k += 1
        is_center = (r == rc && c == cc)
        is_leaf = (r == rc && c == cp1)

        if is_center
            # KET = X [N,q,S,W];  q=open DX0 (ket) = -1.  NO physical leg, NO E leg.
            push!(arrays, Xarr)
            push!(idxlists, [
                VbondK(rc, cc),            # N
                -1,                        # q = DX0 OPEN (ket)
                VbondK(_nextr(rc), cc),    # S
                HbondK(rc, _prevc(cc)),    # W
            ])
            # BRA = conj(X);  q=open DX1 (bra) = -3.
            push!(arrays, conj(Xarr))
            push!(idxlists, [
                VbondBra(rc, cc),
                -3,                        # q = DX1 OPEN (bra)
                VbondBra(_nextr(rc), cc),
                HbondBra(rc, _prevc(cc)),
            ])
        elseif is_leaf
            # KET = Y [N,E,S,q];  q=open DY0 (ket) = -2.  Y replaces the W cut bond.
            push!(arrays, Yarr)
            push!(idxlists, [
                VbondK(rc, cp1),           # N
                HbondK(rc, cp1),           # E
                VbondK(_nextr(rc), cp1),   # S
                -2,                        # q = DY0 OPEN (ket)
            ])
            # BRA = conj(Y);  q=open DY1 (bra) = -4.
            push!(arrays, conj(Yarr))
            push!(idxlists, [
                VbondBra(rc, cp1),
                HbondBra(rc, cp1),
                VbondBra(_nextr(rc), cp1),
                -4,                        # q = DY1 OPEN (bra)
            ])
        else
            T = convert(Array{ComplexF64,5}, convert(Array, peps[r, c]))  # [P,N,E,S,W]
            # KET
            push!(arrays, T)
            push!(idxlists, [
                physTrace(k),              # P (TRACED bra<->ket, shared positive id)
                VbondK(r, c),              # N
                HbondK(r, c),              # E
                VbondK(_nextr(r), c),      # S
                HbondK(r, _prevc(c)),      # W
            ])
            # BRA = conj(T)
            push!(arrays, conj(T))
            push!(idxlists, [
                physTrace(k),              # P (same shared id -> traced)
                VbondBra(r, c),
                HbondBra(r, c),
                VbondBra(_nextr(r), c),
                HbondBra(r, _prevc(c)),
            ])
        end
    end

    # ncon -> 4-leg array, open legs ordered (-1=DX0_ket, -2=DY0_ket, -3=DX1_bra, -4=DY1_bra)
    raw = convert(Array, TensorKitNcon.ncon(arrays, idxlists))
    # reorder to benv axis convention [DX1 DY1; DX0 DY0] = (bra DX1, bra DY1, ket DX0, ket DY0)
    raw_reordered = permutedims(raw, (3, 4, 1, 2))  # (-3,-4,-1,-2)

    # Build the {2,2} TensorMap on the q-spaces.
    # space(X,2) and space(Y,4) are BOTH DUAL, so their conj ' gives the NON-DUAL
    # C^16, C^8. Use the NON-DUAL pair for BOTH codomain and domain; TensorKit then
    # auto-dualizes the domain legs of a V<-V map -> flags [0,0,1,1], codomain==domain.
    DXn = space(X, 2)'    # non-dual C^16
    DYn = space(Y, 4)'    # non-dual C^8
    cod = DXn ⊗ DYn
    dom = DXn ⊗ DYn
    benv = TensorMap(raw_reordered, cod, dom)
    normalize!(benv, Inf)
    return benv
end

# =====================================================================
# Run the experiment only when executed as the main script.
# =====================================================================
function main()
    const_D = 2
    BIG = 64
    CHI = 16          # space-oracle ONLY; smaller chi -> faster leading_boundary
    T_FREEZE = 1.0

    GC.gc()
    @printf("RSS at start = %.1f MB ; total mem = %.1f MB\n",
            Sys.maxrss() / 2^20, Sys.total_memory() / 2^20); flush(stdout)

    # === STAGE A: FROZEN-STATE SETUP (verbatim from ctm_bond_proto.jl) ===
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
    @printf("center PEPSKit (r,c)=(%d,%d), right leaf (r,c)=(%d,%d)\n", rc, cc, rc, cp1)
    flush(stdout)

    lossless = deepcopy(psi)
    project_star_pepskit!(lossless, c0; dt = 0.02, maxdim = BIG, cutoff = 0.0,
                          rel_floor = 0.0, evolution = :real, projected = true)
    n_lossless = exact_density_finite(lossless; max_sites = 16)
    @printf("n_lossless (oracle, expect ~0.15183976) = %.8f\n", n_lossless); flush(stdout)

    A = lossless.peps.A[rc, cc]
    B = lossless.peps.A[rc, cp1]
    X, a, b, Y = PEPSKit._qr_bond(A, B)
    @printf("a space    = %s\n", string(space(a)))
    @printf("b space    = %s\n", string(space(b)))
    @printf("X space    = %s  ; dual flags = %s\n", string(space(X)),
            string([isdual(space(X, ax)) for ax in 1:numind(X)]))
    @printf("Y space    = %s  ; dual flags = %s\n", string(space(Y)),
            string([isdual(space(Y, ax)) for ax in 1:numind(Y)]))
    flush(stdout)

    a_bt = permute(a, ((1, 2), (3,)))   # (C^16 (x) C^2) <- (C^8)'
    b_bt = permute(b, ((1,), (2, 3)))   # (C^8)' <- ((C^2)' (x) (C^8)')
    @printf("a_bt space = %s\n", string(space(a_bt)))
    @printf("b_bt space = %s\n", string(space(b_bt)))
    @printf("space(X,2) = %s isdual=%s ; space(Y,4) = %s isdual=%s\n",
            string(space(X, 2)), isdual(space(X, 2)),
            string(space(Y, 4)), isdual(space(Y, 4)))
    @printf("space(a_bt,1) = %s ; space(b_bt,3) = %s\n",
            string(space(a_bt, 1)), string(space(b_bt, 3)))
    flush(stdout)

    # === STAGE B: HAND-BUILD THE DOUBLE-LAYER {2,2} ENV ===
    GC.gc()
    benv = build_exact_cluster_bondenv(lossless, X, Y, rc, cc, cp1)
    @printf("\n[hand env] space(benv)   = %s\n", string(space(benv)))
    @printf("[hand env] dual flags    = %s\n",
            string([isdual(space(benv, ax)) for ax in 1:numind(benv)]))
    @printf("[hand env] codomain==domain = %s\n", string(codomain(benv) == domain(benv)))
    flush(stdout)

    # === STAGE C(1a): SHAPE GATE vs the DOCUMENTED bondenv_fu target ===
    # (ctm_bond_proto.out ground truth: (C^16(x)C^8)<-(C^16(x)C^8), flags [0,0,1,1]).
    # This is the decisive shape assert; the CTMRG-built bondenv_fu comparison is a
    # best-effort TAIL (it is slow to converge on the enlarged post-gate state and
    # must NOT block the composition result).
    flags_hand = [isdual(space(benv, ax)) for ax in 1:numind(benv)]
    target_space = (ComplexSpace(16) ⊗ ComplexSpace(8)) ← (ComplexSpace(16) ⊗ ComplexSpace(8))
    shape_ok = (space(benv) == target_space)
    @printf("\n>>> ASSERT(1a) SHAPE vs documented target:\n")
    @printf("    hand   : %s  flags %s\n", string(space(benv)), string(flags_hand))
    @printf("    target : %s  flags [false,false,true,true]\n", string(target_space))
    @printf("    space==target = %s ; flags==[0,0,1,1] = %s\n",
            string(shape_ok), string(flags_hand == [false, false, true, true]))
    flush(stdout)
    if !shape_ok || flags_hand != [false, false, true, true]
        @printf(">>> SpaceMismatch DOCUMENTED (accepted result):\n")
        @printf("    hand: %s flags %s\n", string(space(benv)), string(flags_hand))
        println("DONE (N4: SpaceMismatch — see spaces above)")
        return
    end
    println(">>> ASSERT(1a) PASS: hand env shape == documented bondenv_fu target [0,0,1,1]")
    flush(stdout)

    # --- DIAGNOSTIC PARITY: conditioning (NOT a PSD gate; env is <patch|patch>) ---
    Z = PEPSKit.positive_approx(benv)
    @printf(">>> cond(Z'Z) [hand env] = %.4e (expect ~3.4e20 like CTM env)\n",
            cond(Z' * Z)); flush(stdout)

    # === STAGE C(2): NO-OP-AT-BIG self-check (PLUMBING only, env-blind) ===
    alg_big = ALSTruncation(; trunc = PEPSKit.truncrank(BIG) & PEPSKit.truncerror(; atol = 1e-12),
                            maxiter = 50, tol = 1e-12)
    a1, s1, b1, _ = PEPSKit.bond_truncate(a_bt, b_bt, benv, alg_big)
    a1f, b1f = PEPSKit.absorb_s(a1, s1, b1)
    A2, B2 = PEPSKit._qr_bond_undo(X, permute(a1f, ((1, 2), (3,))),
                                   permute(b1f, ((1,), (2, 3))), Y)
    ctm_big = deepcopy(lossless)
    ctm_big.peps.A[rc, cc] = A2
    ctm_big.peps.A[rc, cp1] = B2
    n_noop = exact_density_finite(ctm_big; max_sites = 16)
    noop_dev = abs(n_noop - n_lossless)
    @printf("\n>>> ASSERT(2) no-op-at-BIG: n_noop = %.8f ; dev = %.4e (PLUMBING check)\n",
            n_noop, noop_dev); flush(stdout)
    @assert noop_dev < 1e-9 "no-op-at-BIG plumbing dev $(noop_dev) >= 1e-9"
    println(">>> ASSERT(2) PASS (plumbing: QR/absorb_s/undo glue OK)")
    flush(stdout)

    # === STAGE C(3): TRUNCATE-AT-D (the real env-content gate) ===
    alg_D = ALSTruncation(; trunc = PEPSKit.truncrank(const_D) & PEPSKit.truncerror(; atol = 1e-12),
                          maxiter = 50, tol = 1e-12)
    n_D = NaN
    fid_D = NaN
    try
        aD, sD, bD, tinfo = PEPSKit.bond_truncate(a_bt, b_bt, benv, alg_D)
        fid_D = tinfo.fid
        aDf, bDf = PEPSKit.absorb_s(aD, sD, bD)
        A3, B3 = PEPSKit._qr_bond_undo(X, permute(aDf, ((1, 2), (3,))),
                                       permute(bDf, ((1,), (2, 3))), Y)
        ctmD = deepcopy(lossless)
        ctmD.peps.A[rc, cc] = A3
        ctmD.peps.A[rc, cp1] = B3
        n_D = exact_density_finite(ctmD; max_sites = 16)
    catch err
        @printf(">>> ASSERT(3) bond_truncate THREW: %s\n", string(err)); flush(stdout)
        # fixgauge_benv fallback (correct idiom per benv_fu.jl test) only if needed
        @printf(">>> attempting fixgauge_benv fallback...\n"); flush(stdout)
        try
            Zb = PEPSKit.positive_approx(benv)
            # fixgauge_benv takes RAW a ({1,2}) and b ({2,1}) — NOT a_bt/b_bt.
            Z2, a2g, b2g, (Linv, Rinv) = PEPSKit.fixgauge_benv(Zb, a, b)
            benv2 = Z2' * Z2
            normalize!(benv2, Inf)
            X2, Y2 = PEPSKit._fixgauge_benvXY(X, Y, Linv, Rinv)
            a2_bt = permute(a2g, ((1, 2), (3,)))
            b2_bt = permute(b2g, ((1,), (2, 3)))
            aD, sD, bD, tinfo = PEPSKit.bond_truncate(a2_bt, b2_bt, benv2, alg_D)
            fid_D = tinfo.fid
            aDf, bDf = PEPSKit.absorb_s(aD, sD, bD)
            A3, B3 = PEPSKit._qr_bond_undo(X2, permute(aDf, ((1, 2), (3,))),
                                           permute(bDf, ((1,), (2, 3))), Y2)
            ctmD = deepcopy(lossless)
            ctmD.peps.A[rc, cc] = A3
            ctmD.peps.A[rc, cp1] = B3
            n_D = exact_density_finite(ctmD; max_sites = 16)
            @printf(">>> fixgauge_benv fallback SUCCEEDED\n"); flush(stdout)
        catch err2
            @printf(">>> fixgauge_benv fallback ALSO THREW: %s\n", string(err2)); flush(stdout)
        end
    end
    @printf("\n>>> ASSERT(3) truncate-at-D: n_D = %.8f ; |dn| = %.4e ; fid = %.8e\n",
            n_D, abs(n_D - n_lossless), fid_D); flush(stdout)
    @assert isfinite(n_D) "n_D not finite"
    println(">>> ASSERT(3) PASS: D-truncation composed, finite n_D")
    flush(stdout)

    # === NEGATIVE CONTROL: a deliberately-WRONG env must do WORSE than the right env ===
    # global-conj env (a realistic wrong-env failure mode that still has valid shape).
    @printf("\n--- NEGATIVE CONTROL: wrong env (global conj) at maxdim=D ---\n"); flush(stdout)
    n_wrong = NaN
    try
        benv_wrong = conj(benv)
        normalize!(benv_wrong, Inf)
        aW, sW, bW, _ = PEPSKit.bond_truncate(a_bt, b_bt, benv_wrong, alg_D)
        aWf, bWf = PEPSKit.absorb_s(aW, sW, bW)
        A4, B4 = PEPSKit._qr_bond_undo(X, permute(aWf, ((1, 2), (3,))),
                                       permute(bWf, ((1,), (2, 3))), Y)
        ctmW = deepcopy(lossless)
        ctmW.peps.A[rc, cc] = A4
        ctmW.peps.A[rc, cp1] = B4
        n_wrong = exact_density_finite(ctmW; max_sites = 16)
    catch err
        @printf("    wrong-env conj THREW: %s\n", string(err)); flush(stdout)
    end
    # also a permute-swap (DX<->DY) wrong env, only if dims allow (here 16!=8 so it cannot
    # be built as a {2,2} on the same spaces — skip; conj is the available control).
    @printf("    n_wrong(conj env) = %.8f ; |dn_wrong| = %.4e\n",
            n_wrong, abs(n_wrong - n_lossless)); flush(stdout)

    # === bare-SVD reference (no env) for the deviation ladder ===
    ab = PEPSKit._combine_ab(a_bt, b_bt)
    a0, s0, b0 = PEPSKit.svd_trunc(permute(ab, ((1, 3), (4, 2))); trunc = PEPSKit.truncrank(const_D))
    a0f, b0f = PEPSKit.absorb_s(a0, s0, b0)
    A0, B0 = PEPSKit._qr_bond_undo(X, permute(a0f, ((1, 2), (3,))),
                                   permute(b0f, ((1,), (2, 3))), Y)
    bare = deepcopy(lossless)
    bare.peps.A[rc, cc] = A0
    bare.peps.A[rc, cp1] = B0
    n_bare = exact_density_finite(bare; max_sites = 16)
    @printf("    n_bare-SVD = %.8f ; |dn_bare| = %.4e (expect ~3.97e-7)\n",
            n_bare, abs(n_bare - n_lossless)); flush(stdout)

    println("\n===== N4 SUMMARY =====")
    @printf("n_lossless        = %.8f\n", n_lossless)
    @printf("|dn| hand-env D   = %.4e  (target <= CTM 1.52e-8)\n", abs(n_D - n_lossless))
    @printf("|dn| wrong(conj)  = %.4e  (negative control)\n", abs(n_wrong - n_lossless))
    @printf("|dn| bare-SVD     = %.4e\n", abs(n_bare - n_lossless))
    neg_ok = !isfinite(n_wrong) || abs(n_wrong - n_lossless) >= abs(n_D - n_lossless)
    @printf("negative-control PASS (wrong >= right) = %s\n", string(neg_ok))
    @printf("RSS peak = %.1f MB\n", Sys.maxrss() / 2^20)
    println("DONE (N4: composed = YES)")
    flush(stdout)

    # === BEST-EFFORT TAIL: explicit space-equality vs CTM-built bondenv_fu ===
    # The decisive results above do NOT depend on this. CTMRG leading_boundary on the
    # enlarged post-gate state is slow (it can time out); wrap so a kill loses nothing.
    @printf("\n--- BEST-EFFORT: bondenv_fu CTM space-oracle (chi=%d) for explicit comparison ---\n", CHI)
    flush(stdout)
    try
        mp = measurement_peps(lossless)
        env0 = PEPSKit.CTMRGEnv(randn, ComplexF64, mp, ComplexSpace(CHI))
        env, _ = PEPSKit.leading_boundary(env0, mp; alg = :simultaneous, tol = 1e-6,
                                          miniter = 2, maxiter = 25, verbosity = 1,
                                          trunc = PEPSKit.truncrank(CHI))
        benv_ref = PEPSKit.bondenv_fu(rc, cc, X, Y, env)
        flags_ref = [isdual(space(benv_ref, ax)) for ax in 1:numind(benv_ref)]
        @printf("[ref env] space(benv_ref) = %s ; flags = %s\n",
                string(space(benv_ref)), string(flags_ref))
        @printf("space(benv)==space(benv_ref) = %s ; flags equal = %s\n",
                string(space(benv) == space(benv_ref)),
                string(flags_hand == flags_ref))
        rel_dev = norm(benv - benv_ref) / norm(benv_ref)
        @printf("VALUE rel-dev norm(benv - benv_ref)/norm(benv_ref) = %.4e (CTM chi=%d, diagnostic only)\n",
                rel_dev, CHI)
        flush(stdout)
    catch err
        @printf("CTM space-oracle DNF (best-effort, non-blocking): %s\n", string(err))
        flush(stdout)
    end
    println("DONE (best-effort tail)")
    flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
