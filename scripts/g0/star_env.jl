# star_env.jl — D-generalized, 4-direction {2,2} bond-environment builder for the
# square-star simple update, plus an env-aware prototype star update.
#
# Builds on spike_n4_exact_cluster_bond.jl (build_exact_cluster_bondenv, the
# validated gold-standard horizontal-right env) and ntu_bondenv.jl
# (build_ntu_bondenv patch route). It generalizes BOTH to:
#   (1) arbitrary D — the 16/8 dim asserts are dropped; q-dims come from spaces.
#   (2) ALL FOUR star bonds — right/left (horizontal) AND up/down (VERTICAL).
#
# ROUTE (verticalSpec): rotate the center & leaf PEPSTensors so the CUT bond is
# in the EAST(center)/WEST(leaf) slot, then call PEPSKit._qr_bond (which is
# geometry-blind and ALWAYS splits A's east / B's west). After rotation the q-leg
# sits at X-axis-2 / Y-axis-4 in ALL directions (uniform dual flags). The builder
# then emits the TRUE-geometry interior bond ids for the three non-cut center legs
# and three non-cut leaf legs, so NO explicit un-rotate of X/Y is needed — the
# unrotation is absorbed by emitting the right ids.
#
# Per the SPEC REVIEW, the env is built from the PRE-STAR (plain convention-B)
# state — peps.A[center]/peps.A[leaf] directly, mirroring spike:189-191 — NOT from
# the loop's double-weighted Qfactors/reduced. The env-aware star update here uses
# the "all-4-envs-once-from-pre-star-state" STALE variant (well-posed, matches the
# validated object). Per-arm mid-loop rebuild is intentionally NOT attempted.

include(joinpath(@__DIR__, "ntu_bondenv.jl"))  # -> build_exact_cluster_bondenv, build_ntu_bondenv

using PEPSKit: _qr_bond, _qr_bond_undo, positive_approx, fixgauge_benv,
               _fixgauge_benvXY, bond_truncate, ALSTruncation, truncrank,
               truncerror, absorb_s
using TensorKit: permute, isdual, space, numind, TensorMap, normalize!
import TensorKit

# ----------------------------------------------------------------------------
# rotate_peps — cyclic CCW rotation of the 4 virtual legs of a PEPSTensor.
# PEPS leg order is [P, N, E, S, W] = axes [1,2,3,4,5].
# One CCW step maps (N,E,S,W) -> (W,N,E,S): new-N = old-W, new-E = old-N,
# new-S = old-E, new-W = old-S. (P fixed at axis 1.)
# We need r CCW steps where r maps the CUT direction onto EAST:
#   :right => 0, :up => 1, :left => 2, :down => 3.
# After r CCW steps, the leg that was at direction d sits at d shifted... we build
# the permutation directly. With axes (N,E,S,W) = (2,3,4,5):
#   1 CCW: new[N,E,S,W] = old[W,N,E,S] = old axes (5,2,3,4)
# so rotate_peps(T,1) permutes virtual axes to (5,2,3,4); compose for r>1.
# ----------------------------------------------------------------------------
const _CCW1 = (5, 2, 3, 4)  # new (N,E,S,W) <- old axes

function _ccw_axes(r::Int)
    # returns the old-axis indices (among 2..5) for new (N,E,S,W) after r CCW steps
    axes = collect(2:5)  # identity: new[N,E,S,W] = old (2,3,4,5)
    for _ in 1:mod(r, 4)
        # apply one CCW: new(N,E,S,W) = old(W,N,E,S) => index by _CCW1
        axes = [axes[a - 1] for a in _CCW1]  # _CCW1 are absolute axes 5,2,3,4 -> local idx a-1
    end
    return axes
end

"""
    rotate_peps(T, r) -> PEPSTensor

Rotate a PEPSTensor `T` (P<-N,E,S,W) by `r` CCW quarter-turns of its 4 virtual
legs (physical leg fixed). r=0 is identity. Implemented as a permute that keeps
codomain={P} and domain={N,E,S,W}, so the result is again a valid PEPSTensor.
"""
function rotate_peps(T, r::Int)
    r = mod(r, 4)
    r == 0 && return T
    va = _ccw_axes(r)  # old axes (2..5) for the new N,E,S,W slots
    return permute(T, ((1,), (va[1], va[2], va[3], va[4])))
end

# CCW rotation count that maps the cut direction onto EAST.
_rot_for_dir(dir::Symbol) = dir === :right ? 0 :
                            dir === :up    ? 1 :
                            dir === :left  ? 2 :
                            dir === :down  ? 3 :
                            error("bad dir $dir")

# ----------------------------------------------------------------------------
# Lattice id helpers (mirror build_exact_cluster_bondenv exactly, parameterized).
# ----------------------------------------------------------------------------
struct _IdFns
    Nr::Int
    Nc::Int
    nsites::Int
end
_idfns(Nr, Nc) = _IdFns(Nr, Nc, Nr * Nc)
_nextr(f::_IdFns, r) = r == f.Nr ? 1 : r + 1
_prevr(f::_IdFns, r) = r == 1 ? f.Nr : r - 1
_nextc(f::_IdFns, c) = c == f.Nc ? 1 : c + 1
_prevc(f::_IdFns, c) = c == 1 ? f.Nc : c - 1
_HbondK(f::_IdFns, r, c) = (r - 1) * f.Nc + c
_VbondK(f::_IdFns, r, c) = f.nsites + (r - 1) * f.Nc + c
_HbondBra(f::_IdFns, r, c) = 2 * f.nsites + _HbondK(f, r, c)
_VbondBra(f::_IdFns, r, c) = 2 * f.nsites + _VbondK(f, r, c)
_physTrace(f::_IdFns, k) = 4 * f.nsites + k

# interior id for a TRUE-geometry virtual leg in [N,E,S,W] direction.
_interiorKetId(f::_IdFns, r, c, dir) =
    dir === :N ? _VbondK(f, r, c) :
    dir === :E ? _HbondK(f, r, c) :
    dir === :S ? _VbondK(f, _nextr(f, r), c) :
                 _HbondK(f, r, _prevc(f, c))   # :W
_interiorBraId(f::_IdFns, r, c, dir) =
    dir === :N ? _VbondBra(f, r, c) :
    dir === :E ? _HbondBra(f, r, c) :
    dir === :S ? _VbondBra(f, _nextr(f, r), c) :
                 _HbondBra(f, r, _prevc(f, c))  # :W

# ----------------------------------------------------------------------------
# Per-direction routing table.
# Returns, for the center and the leaf, the (cut_dir_center, cut_dir_leaf,
# leaf cell, the three non-cut directions of center & leaf in ROTATED-frame
# axis order). After rotation the X/Y axes correspond to rotated (N,E,S,W);
# the cut sits at center.E / leaf.W in the rotated frame, with q at X-axis2,
# Y-axis4. We map each ROTATED axis of X (N,E,S,W after rotation) back to the
# TRUE direction of the unrotated tensor, and emit its true interior id.
#
# For the CENTER, the rotation is r CCW. After r CCW steps, the rotated leg at
# rotated-direction d came from the TRUE direction that is r steps CW of d.
# One CW step (inverse of CCW): (N,E,S,W) -> (E,S,W,N), i.e. rotated-N<-true-?
# We computed _ccw_axes(r): new (N,E,S,W) <- old axes. So rotated-N comes from
# the true direction given by which of (2..5) _ccw_axes(r)[1] is, etc.
# ----------------------------------------------------------------------------
const _AX2DIR = Dict(2 => :N, 3 => :E, 4 => :S, 5 => :W)

# For rotation r, return mapping rotated-direction (in N,E,S,W order) -> true dir.
function _rotated_to_true(r::Int)
    va = _ccw_axes(r)  # old axes for rotated N,E,S,W
    return (rotN = _AX2DIR[va[1]], rotE = _AX2DIR[va[2]],
            rotS = _AX2DIR[va[3]], rotW = _AX2DIR[va[4]])
end

"""
    build_star_bond_env(state, X, Y, dir, rc, cc, leaf_rc, leaf_cc; env=:exact_cluster, patch=:torus)
        -> benv::TensorMap{T,S,2,2}

D-generalized, 4-direction {2,2} bond environment for the center(rc,cc) ->
leaf(leaf_rc,leaf_cc) bond in direction `dir` ∈ (:right,:up,:left,:down).

`X`, `Y` are the q-reduced PEPSOrth tensors from `_qr_bond(rotate_peps(center,r),
rotate_peps(leaf,r))` — i.e. in the ROTATED frame, q at X-axis-2 / Y-axis-4 in
all directions. The builder routes each rotated X/Y axis to its TRUE interior id.

`env=:exact_cluster` uses the full torus (gold standard); `:ntu` uses a local
patch (passes `patch` through). For vertical dirs only `:torus` patch is
validated.
"""
function build_star_bond_env(state, X, Y, dir::Symbol, rc::Int, cc::Int,
                             leaf_rc::Int, leaf_cc::Int;
                             env::Symbol = :exact_cluster, patch::Symbol = :torus)
    cell = state.unitcell
    Nr, Nc = cell.Ly, cell.Lx
    f = _idfns(Nr, Nc)
    peps = measurement_peps(state)

    # D-agnostic dual-flag asserts (replace the dropped ==16/==8 dim asserts).
    @assert isdual(space(X, 2)) "X q-leg (axis 2) must be dual"
    @assert isdual(space(Y, 4)) "Y q-leg (axis 4) must be dual"
    # vertical self-aliasing guard (review strengthened to Nr>=4 && Nc>=4)
    @assert Nr >= 4 && Nc >= 4 "vertical/horizontal star env needs Nr>=4 && Nc>=4 (self-aliasing guard)"

    r = _rot_for_dir(dir)
    rt = _rotated_to_true(r)  # rotated (N,E,S,W) -> true dir for the CENTER

    # The LEAF is rotated by the SAME r. In the rotated frame the leaf's cut sits
    # at WEST (its q is at axis 4). Its three non-cut rotated dirs (N,E,S) map to
    # true leaf dirs by the SAME _rotated_to_true(r).

    # Build the ncon id lists, iterating cells in (r,c) order with the two open
    # q-legs carrying the cut bond. The center array X has axes [N,E,S,W]_rotated
    # = [rotN, q(=rotE→cut), rotS, rotW]? No: _qr_bond gives X legs [N,q,S,W] with
    # q at axis 2 (it split A's EAST). So X axes are: ax1=rotN, ax2=q (cut/east),
    # ax3=rotS, ax4=rotW. The non-cut X axes are ax1(rotN),ax3(rotS),ax4(rotW).
    # Y axes [N,E,S,q]: ax1=rotN, ax2=rotE, ax3=rotS, ax4=q (cut/west).
    # Non-cut Y axes: ax1(rotN),ax2(rotE),ax3(rotS).

    Xarr = convert(Array{ComplexF64,4}, convert(Array, X))   # [rotN, q, rotS, rotW]
    Yarr = convert(Array{ComplexF64,4}, convert(Array, Y))   # [rotN, rotE, rotS, q]

    # true center ids for the three non-cut rotated axes (rotN, rotS, rotW)
    cN_id_k = _interiorKetId(f, rc, cc, rt.rotN)
    cS_id_k = _interiorKetId(f, rc, cc, rt.rotS)
    cW_id_k = _interiorKetId(f, rc, cc, rt.rotW)
    cN_id_b = _interiorBraId(f, rc, cc, rt.rotN)
    cS_id_b = _interiorBraId(f, rc, cc, rt.rotS)
    cW_id_b = _interiorBraId(f, rc, cc, rt.rotW)

    # true leaf ids for the three non-cut rotated axes (rotN, rotE, rotS)
    lN_id_k = _interiorKetId(f, leaf_rc, leaf_cc, rt.rotN)
    lE_id_k = _interiorKetId(f, leaf_rc, leaf_cc, rt.rotE)
    lS_id_k = _interiorKetId(f, leaf_rc, leaf_cc, rt.rotS)
    lN_id_b = _interiorBraId(f, leaf_rc, leaf_cc, rt.rotN)
    lE_id_b = _interiorBraId(f, leaf_rc, leaf_cc, rt.rotE)
    lS_id_b = _interiorBraId(f, leaf_rc, leaf_cc, rt.rotS)

    if env === :ntu && patch !== :torus
        # NTU patch route: only validated for horizontal; for vertical use the
        # rotated-id machinery with a patch set. We restrict to :torus elsewhere.
        return _build_star_bond_env_ntu(state, X, Y, dir, rc, cc, leaf_rc, leaf_cc,
                                        f, peps, Xarr, Yarr, rt, patch)
    end

    # ---- full-torus (exact-cluster) assembly ----
    arrays = Vector{Array{ComplexF64}}()
    idxlists = Vector{Vector{Int}}()

    k = 0
    for rr in 1:Nr, ccol in 1:Nc
        k += 1
        is_center = (rr == rc && ccol == cc)
        is_leaf = (rr == leaf_rc && ccol == leaf_cc)

        if is_center
            # KET center = X [rotN, q(=-1), rotS, rotW]
            push!(arrays, Xarr)
            push!(idxlists, [cN_id_k, -1, cS_id_k, cW_id_k])
            push!(arrays, conj(Xarr))
            push!(idxlists, [cN_id_b, -3, cS_id_b, cW_id_b])
        elseif is_leaf
            # KET leaf = Y [rotN, rotE, rotS, q(=-2)]
            push!(arrays, Yarr)
            push!(idxlists, [lN_id_k, lE_id_k, lS_id_k, -2])
            push!(arrays, conj(Yarr))
            push!(idxlists, [lN_id_b, lE_id_b, lS_id_b, -4])
        else
            T = convert(Array{ComplexF64,5}, convert(Array, peps[rr, ccol]))  # [P,N,E,S,W]
            push!(arrays, T)
            push!(idxlists, [
                _physTrace(f, k),
                _VbondK(f, rr, ccol),
                _HbondK(f, rr, ccol),
                _VbondK(f, _nextr(f, rr), ccol),
                _HbondK(f, rr, _prevc(f, ccol)),
            ])
            push!(arrays, conj(T))
            push!(idxlists, [
                _physTrace(f, k),
                _VbondBra(f, rr, ccol),
                _HbondBra(f, rr, ccol),
                _VbondBra(f, _nextr(f, rr), ccol),
                _HbondBra(f, rr, _prevc(f, ccol)),
            ])
        end
    end

    raw = convert(Array, TensorKitNcon.ncon(arrays, idxlists))
    raw_reordered = permutedims(raw, (3, 4, 1, 2))  # [DX1 DY1; DX0 DY0]
    DXn = space(X, 2)'   # non-dual
    DYn = space(Y, 4)'   # non-dual
    benv = TensorMap(raw_reordered, DXn ⊗ DYn, DXn ⊗ DYn)
    normalize!(benv, Inf)
    return benv
end

# NTU patch variant — reuses the rotated-id routing with a patch set + boundary
# trace. Only used when env=:ntu && patch != :torus.
function _build_star_bond_env_ntu(state, X, Y, dir, rc, cc, leaf_rc, leaf_cc,
                                  f::_IdFns, peps, Xarr, Yarr, rt, patch::Symbol)
    Nr, Nc = f.Nr, f.Nc
    nsites = f.nsites
    boundary_start = 5 * nsites + 1
    boundary_counter = Ref(boundary_start - 1)
    next_boundary() = (boundary_counter[] += 1; boundary_counter[])

    # P6 patch cells (the loop-carrying 6-cell band straddling the cut bond),
    # reached for ALL non-torus dirs (horizontal AND vertical). The band is the
    # 90-deg analog of the validated horizontal P6: horizontal = rows{rm1,rc,rp1}
    # x cols{cc,leaf_cc}; vertical = cols{cm1,cc,cp1} x rows{rc,leaf_rc}. Both
    # contain the TWO length-4 plaquettes that thread the bond (audit-verified:
    # geometry correct, boundary traces clean, PSD min-eig>0, torus==exact bitwise
    # for up+down). Non-patch neighbours get a fresh ket==bra boundary id (delta
    # trace = open-boundary approx). Validated at VALUE level by the vertical-P6
    # transpose oracle (G5, envaware_integration.jl).
    rm1 = _prevr(f, rc); rp1 = _nextr(f, rc)
    cm1 = _prevc(f, cc); cp1c = _nextc(f, cc)
    # Build a 6-cell P6-equivalent around the cut: center, leaf, and the 4 cells
    # adjacent to the cut perpendicular to the bond.
    if dir in (:right, :left)
        cs = [(rm1, cc), (rm1, leaf_cc), (rc, cc), (rc, leaf_cc), (rp1, cc), (rp1, leaf_cc)]
    else # vertical: perpendicular is columns
        cs = [(rc, cm1), (leaf_rc, cm1), (rc, cc), (leaf_rc, cc), (rc, cp1c), (leaf_rc, cp1c)]
    end
    patch_set = Set(cs)

    neighbour(rr, ccol, d) = d === :N ? (_prevr(f, rr), ccol) :
                             d === :E ? (rr, _nextc(f, ccol)) :
                             d === :S ? (_nextr(f, rr), ccol) :
                                        (rr, _prevc(f, ccol))
    is_internal(rr, ccol, d) = neighbour(rr, ccol, d) in patch_set
    boundary_ids = Dict{Tuple{Int,Int,Symbol},Int}()
    boundary_id(rr, ccol, d) = get!(boundary_ids, (rr, ccol, d)) do; next_boundary() end
    legid(rr, ccol, d, layer) = is_internal(rr, ccol, d) ?
        (layer === :ket ? _interiorKetId(f, rr, ccol, d) : _interiorBraId(f, rr, ccol, d)) :
        boundary_id(rr, ccol, d)

    arrays = Vector{Array{ComplexF64}}()
    idxlists = Vector{Vector{Int}}()
    k = 0
    for (rr, ccol) in cs
        k += 1
        is_center = (rr == rc && ccol == cc)
        is_leaf = (rr == leaf_rc && ccol == leaf_cc)
        if is_center
            push!(arrays, Xarr)
            push!(idxlists, [legid(rc, cc, rt.rotN, :ket), -1,
                             legid(rc, cc, rt.rotS, :ket), legid(rc, cc, rt.rotW, :ket)])
            push!(arrays, conj(Xarr))
            push!(idxlists, [legid(rc, cc, rt.rotN, :bra), -3,
                             legid(rc, cc, rt.rotS, :bra), legid(rc, cc, rt.rotW, :bra)])
        elseif is_leaf
            push!(arrays, Yarr)
            push!(idxlists, [legid(leaf_rc, leaf_cc, rt.rotN, :ket),
                             legid(leaf_rc, leaf_cc, rt.rotE, :ket),
                             legid(leaf_rc, leaf_cc, rt.rotS, :ket), -2])
            push!(arrays, conj(Yarr))
            push!(idxlists, [legid(leaf_rc, leaf_cc, rt.rotN, :bra),
                             legid(leaf_rc, leaf_cc, rt.rotE, :bra),
                             legid(leaf_rc, leaf_cc, rt.rotS, :bra), -4])
        else
            T = convert(Array{ComplexF64,5}, convert(Array, peps[rr, ccol]))
            push!(arrays, T)
            push!(idxlists, [_physTrace(f, k), legid(rr, ccol, :N, :ket),
                             legid(rr, ccol, :E, :ket), legid(rr, ccol, :S, :ket),
                             legid(rr, ccol, :W, :ket)])
            push!(arrays, conj(T))
            push!(idxlists, [_physTrace(f, k), legid(rr, ccol, :N, :bra),
                             legid(rr, ccol, :E, :bra), legid(rr, ccol, :S, :bra),
                             legid(rr, ccol, :W, :bra)])
        end
    end
    # structural safety (copied from build_ntu_bondenv:262-275): every boundary id
    # must be used EXACTLY twice (this site's ket + bra -> delta trace); a mis-paired
    # boundary leg (e.g. a future widened band that wraps onto itself) is caught here
    # before the contraction, on exactly this previously-unguarded path.
    idcount = Dict{Int,Int}()
    for idl in idxlists, id in idl
        id >= boundary_start && (idcount[id] = get(idcount, id, 0) + 1)
    end
    for (id, n) in idcount
        @assert n == 2 "boundary id $id used $n times (expected 2: ket+bra)"
    end
    @assert boundary_start > 4 * nsites + nsites "boundary ids collide with physTrace block"
    raw = convert(Array, TensorKitNcon.ncon(arrays, idxlists))
    raw_reordered = permutedims(raw, (3, 4, 1, 2))
    DXn = space(X, 2)'
    DYn = space(Y, 4)'
    benv = TensorMap(raw_reordered, DXn ⊗ DYn, DXn ⊗ DYn)
    normalize!(benv, Inf)
    return benv
end

# ----------------------------------------------------------------------------
# build_star_envs_from_state — build all 4 (or a chosen subset) bond envs from
# the PRE-STAR (plain convention-B) state, returning a Dict dir -> (benv, X, Y).
# This mirrors spike:189-191 per direction via rotate_peps + _qr_bond.
# ----------------------------------------------------------------------------
function build_star_env_one(state, rc, cc, leaf_rc, leaf_cc, dir;
                            env=:exact_cluster, patch=:torus)
    A = state.peps.A[rc, cc]
    B = state.peps.A[leaf_rc, leaf_cc]
    r = _rot_for_dir(dir)
    Ar = rotate_peps(A, r)
    Br = rotate_peps(B, r)
    X, a, b, Y = _qr_bond(Ar, Br)
    benv = build_star_bond_env(state, X, Y, dir, rc, cc, leaf_rc, leaf_cc; env=env, patch=patch)
    return (benv=benv, X=X, Y=Y, a=a, b=b)
end

# ----------------------------------------------------------------------------
# env_truncate_bond! — env-aware truncation of ONE center<->leaf bond, in place.
# Mirrors the spike's fixgauge+bond_truncate+absorb_s+_qr_bond_undo pipeline,
# generalized to direction `dir` via the rotate route. The env is REBUILT from
# the CURRENT state (so per-arm rebuild captures the previous arm's truncation —
# peel-staleness handled). Writes the truncated center+leaf back into state.peps.
# Returns (kept_dim, trunc_err, cond_ZZ, fell_back::Bool).
# ----------------------------------------------------------------------------
function env_truncate_bond!(state, rc, cc, leaf_rc, leaf_cc, dir, maxdim, cutoff;
                            env=:exact_cluster, patch=:torus, trunc_alg=:als)
    A = state.peps.A[rc, cc]
    B = state.peps.A[leaf_rc, leaf_cc]
    Nr, Nc = size(state.peps)
    # FIX (onset probe): keep the SUWeight ledger in lockstep with the truncated
    # bond. The lossless star (project_star_pepskit! maxdim=BIG inside
    # project_star_env!) GROWS this bond's weight to the full untruncated size;
    # the env truncation below shrinks peps.A but, without this, leaves
    # weights.data[slot] stale -> the NEXT center's absorb_weight throws
    # SpaceMismatch(ℂ^kept ≠ ℂ^grown). Mirror project_star_pepskit!:367-368:
    # store the non-dual truncated singular values S at this bond's ledger slot.
    bslot = SquarePXPDynamics.PEPSKitStarUpdate._bond_slot(rc, cc, Nr, Nc, dir)
    r = _rot_for_dir(dir)
    Ar = rotate_peps(A, r); Br = rotate_peps(B, r)
    X, a, b, Y = _qr_bond(Ar, Br)
    if env === :bare
        # CONTROL (causal attribution): lossless-star-then-bare-SVD-truncate, NO
        # environment. Isolates whether the loop-carrying ENVIRONMENT — vs merely
        # the lossless-grow-then-retruncate PROCEDURE that project_star_env! shares
        # with :exact_cluster — is what suppresses the D=2 blowup. If :bare spikes
        # like native-bare while :exact_cluster stays smooth, the env is the cause.
        ab0 = PEPSKit._combine_ab(permute(a, ((1, 2), (3,))), permute(b, ((1,), (2, 3))))
        a0, s0, b0 = PEPSKit.svd_trunc(permute(ab0, ((1, 3), (4, 2))); trunc=truncrank(maxdim))
        if isdual(space(s0, 1))
            a0, s0, b0 = PEPSKit.flip_svd(a0, s0, b0)
        end
        a0f, b0f = absorb_s(a0, s0, b0)
        A_new, B_new = _qr_bond_undo(X, permute(a0f, ((1, 2), (3,))),
                                     permute(b0f, ((1,), (2, 3))), Y)
        A_new = rotate_peps(A_new, mod(4 - r, 4))
        B_new = rotate_peps(B_new, mod(4 - r, 4))
        state.peps.A[rc, cc] = A_new
        state.peps.A[leaf_rc, leaf_cc] = B_new
        state.weights.data[bslot[1], bslot[2], bslot[3]] = s0
        return (dim(space(s0, 1)), NaN, NaN, false)
    end
    atol = cutoff > 0 ? Float64(cutoff) : 1e-12
    # BOTH engines are iterative. FullEnvTruncation alternates inner linear solves
    # (KrylovKit GMRES, fullenv_truncation.jl:245/255 -> the "GMRES stopped" warnings
    # at residual ~1e-10 are benign) but converges to a MUCH tighter no-op floor
    # (~1e-10 full-star vs ALS ~5e-9 at D=2; D=3: G3 als 6.8e-9). ALS is cheaper but
    # leaves a higher residual floor, so any trajectory science result MUST use :full.
    alg = trunc_alg === :full ?
        PEPSKit.FullEnvTruncation(; trunc=truncrank(maxdim) & truncerror(; atol=atol)) :
        ALSTruncation(; trunc=truncrank(maxdim) & truncerror(; atol=atol), maxiter=50, tol=1e-12)
    fell_back = false
    local aD, sD, bD, info, X2, Y2
    cnd = NaN
    benv = build_star_bond_env(state, X, Y, dir, rc, cc, leaf_rc, leaf_cc; env=env, patch=patch)
    try
        Z = positive_approx(benv)
        cnd = cond(Z' * Z)
        Z2, a2, b2, (Linv, Rinv) = fixgauge_benv(Z, a, b)
        benv2 = Z2' * Z2; normalize!(benv2, Inf)
        X2, Y2 = _fixgauge_benvXY(X, Y, Linv, Rinv)
        a2_bt = permute(a2, ((1, 2), (3,)))
        b2_bt = permute(b2, ((1,), (2, 3)))
        aD, sD, bD, info = bond_truncate(a2_bt, b2_bt, benv2, alg)
    catch err
        # fallback: bare SVD split of the same bond (no env)
        fell_back = true
        ab = PEPSKit._combine_ab(permute(a, ((1, 2), (3,))), permute(b, ((1,), (2, 3))))
        a0, s0, b0 = PEPSKit.svd_trunc(permute(ab, ((1, 3), (4, 2))); trunc=truncrank(maxdim))
        if isdual(space(s0, 1))
            a0, s0, b0 = PEPSKit.flip_svd(a0, s0, b0)
        end
        a0f, b0f = absorb_s(a0, s0, b0)
        A_new, B_new = _qr_bond_undo(X, permute(a0f, ((1, 2), (3,))),
                                     permute(b0f, ((1,), (2, 3))), Y)
        A_new = rotate_peps(A_new, mod(4 - r, 4))
        B_new = rotate_peps(B_new, mod(4 - r, 4))
        state.peps.A[rc, cc] = A_new
        state.peps.A[leaf_rc, leaf_cc] = B_new
        state.weights.data[bslot[1], bslot[2], bslot[3]] = s0   # ledger lockstep
        return (dim(space(s0, 1)), NaN, cnd, fell_back)
    end
    if isdual(space(sD, 1))          # SUWeight requires a non-dual stored weight
        aD, sD, bD = PEPSKit.flip_svd(aD, sD, bD)
    end
    aDf, bDf = absorb_s(aD, sD, bD)
    A_new, B_new = _qr_bond_undo(X2, permute(aDf, ((1, 2), (3,))),
                                 permute(bDf, ((1,), (2, 3))), Y2)
    A_new = rotate_peps(A_new, mod(4 - r, 4))
    B_new = rotate_peps(B_new, mod(4 - r, 4))
    state.peps.A[rc, cc] = A_new
    state.peps.A[leaf_rc, leaf_cc] = B_new
    state.weights.data[bslot[1], bslot[2], bslot[3]] = sD   # ledger lockstep
    kept = dim(space(sD, 1))
    terr = sqrt(max(0.0, 1.0 - info.fid))
    return (kept, terr, cnd, fell_back)
end

# ----------------------------------------------------------------------------
# project_star_env! — env-aware prototype star update.
#
# Well-posed design (per SPEC REVIEW): run the LOSSLESS bare star update first
# (maxdim=BIG, no truncation) so the center+4 leaves are updated by the gate
# without any bond loss; THEN apply env-aware truncation to each of the 4
# center<->leaf bonds, REBUILDING the env per arm from the current (partially
# truncated) state so arm k+1's env sees arm k's freshly-truncated bond. The
# env is built from the plain convention-B post-gate tensors (peps.A directly),
# matching the validated object — NOT from double-weighted reduced cores.
#
# Returns a NamedTuple with per-dir truncerrs, keptdims, conds, fellbacks.
# ----------------------------------------------------------------------------
function project_star_env!(state, center; dt, maxdim, cutoff=1e-12,
                           env=:exact_cluster, patch=:torus, trunc_alg=:als,
                           evolution=:real, projected=true,
                           order=(:right, :up, :left, :down))
    cell = state.unitcell
    cidx = SquarePXPDynamics.squarecoord_to_cartesianindex(cell, center)
    rc, cc = cidx[1], cidx[2]
    Nr, Nc = size(state.peps)
    # 1. lossless bare star update (no truncation): the gate is applied, all bonds
    #    grow to their exact (untruncated) size.
    project_star_pepskit!(state, center; dt=dt, maxdim=10^6, cutoff=0.0,
                          rel_floor=0.0, evolution=evolution, projected=projected)
    # 2. env-aware truncate each arm, rebuilding env per arm from current state.
    leafcell = Dict(
        :right => (rc, mod1(cc + 1, Nc)),
        :up    => (mod1(rc - 1, Nr), cc),
        :left  => (rc, mod1(cc - 1, Nc)),
        :down  => (mod1(rc + 1, Nr), cc),
    )
    truncerrs = Dict{Symbol,Float64}(); keptdims = Dict{Symbol,Int}()
    conds = Dict{Symbol,Float64}(); fellbacks = Dict{Symbol,Bool}()
    for dir in order
        (lr, lc) = leafcell[dir]
        kept, terr, cnd, fb = env_truncate_bond!(state, rc, cc, lr, lc, dir, maxdim, cutoff;
                                                 env=env, patch=patch, trunc_alg=trunc_alg)
        truncerrs[dir] = terr; keptdims[dir] = kept; conds[dir] = cnd; fellbacks[dir] = fb
    end
    return (truncerrs=truncerrs, keptdims=keptdims, conds=conds, fellbacks=fellbacks)
end
