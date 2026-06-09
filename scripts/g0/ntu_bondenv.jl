# NTU FINITE-PATCH bond environment — the cheap, trajectory-feasible production
# engine for CTM-aware bond truncation.
#
# build_ntu_bondenv builds the {2,2} BondEnv for ONE horizontal center->right
# bond as the bra-ket overlap <patch|patch> of a LOCAL PATCH (the cross + a
# one-plaquette ring), with the active bond's two virtual q-legs OPEN, all
# physical legs TRACED, and the patch-BOUNDARY virtual bonds closed LOCALLY
# (traced, or SUWeight-lambda mean-field-dressed). PSD by construction (it is a
# literal Gram form <patch|patch>, no positive_approx).
#
# It REUSES build_exact_cluster_bondenv's id-allocation + assembly skeleton
# VERBATIM (only the iteration set + per-leg internal-vs-boundary routing change).
# patch=:torus reproduces build_exact_cluster_bondenv EXACTLY (the equality
# oracle G-MAXPATCH).
#
# IMPORTANT: the contraction is over RAW Arrays (convert(Array,...) + ncon), which
# is plain TensorOperations delta-contraction and IGNORES TensorKit dual flags
# entirely. So a boundary leg traced ket-to-bra via a shared positive integer id
# is byte-identical in mechanism to the validated physical trace (physTrace). The
# only [0,0,1,1] dual flags come from the final TensorMap rewrap of the OPEN q-legs
# (the unchanged, validated assembly tail). DO NOT reintroduce a TensorMap-based
# contraction here — it WOULD be dual-sensitive.

include(joinpath(@__DIR__, "spike_n4_exact_cluster_bond.jl"))

using PEPSKit: NORTH, EAST, SOUTH, WEST, absorb_weight

"""
    patch_coords(rc, cc, cp1, Nr, Nc, patch) -> Vector{Tuple{Int,Int}}

The list of (r,c) cells in the patch for the horizontal center(rc,cc)->leaf(rc,cp1)
bond. patch ∈ (:P6, :P10, :P12, :torus).

  :P6  = the 2-col x 3-row block straddling the cut bond (the two length-4
         plaquettes that thread the bond — UPPER + LOWER — tile exactly this set).
  :P10 = P6 + the two cross arms re-coupling center.W and leaf.E (8 cells).
  :P12 = P10 + the 4 corner cells giving those arms a double-layer environment.
  :torus = all (r,c) (reproduces build_exact_cluster_bondenv exactly).
"""
function patch_coords(rc::Int, cc::Int, cp1::Int, Nr::Int, Nc::Int, patch::Symbol)
    _prevr(r) = r == 1 ? Nr : r - 1
    _nextr(r) = r == Nr ? 1 : r + 1
    _prevc(c) = c == 1 ? Nc : c - 1
    _nextc(c) = c == Nc ? 1 : c + 1
    rm1 = _prevr(rc)
    rp1 = _nextr(rc)
    if patch === :torus
        return [(r, c) for r in 1:Nr for c in 1:Nc]
    elseif patch === :P6
        return [(rm1, cc), (rm1, cp1), (rc, cc), (rc, cp1), (rp1, cc), (rp1, cp1)]
    elseif patch === :P10
        cm1 = _prevc(cc)
        cp2 = _nextc(cp1)
        return [(rm1, cc), (rm1, cp1), (rc, cc), (rc, cp1), (rp1, cc), (rp1, cp1),
                (rc, cm1), (rc, cp2)]
    elseif patch === :P12
        cm1 = _prevc(cc)
        cp2 = _nextc(cp1)
        return [(rm1, cc), (rm1, cp1), (rc, cc), (rc, cp1), (rp1, cc), (rp1, cp1),
                (rc, cm1), (rc, cp2),
                (rm1, cm1), (rp1, cm1), (rm1, cp2), (rp1, cp2)]
    else
        error("unknown patch $patch (use :P6, :P10, :P12, :torus)")
    end
end

# A deliberately-WRONG patch dropping the LOWER plaquette row (only the UPPER
# loop -> still a tree through the bond? no: it keeps UPPER as a real loop but the
# control we want is "drop one plaquette" -> 4 cells, ONE loop). Used by the
# LOOP-CONTENT control.
function patch_coords_droprow(rc::Int, cc::Int, cp1::Int, Nr::Int, Nc::Int)
    _prevr(r) = r == 1 ? Nr : r - 1
    rm1 = _prevr(rc)
    # UPPER plaquette only: {center, leaf, (rm1,cc), (rm1,cp1)} — drops the rp1 row,
    # so the LOWER loop is severed.
    return [(rm1, cc), (rm1, cp1), (rc, cc), (rc, cp1)]
end

"""
    build_ntu_bondenv(state, X, Y, rc, cc, cp1; patch=:P6, boundary=:trace, coords=nothing)
        -> benv::TensorMap{T,S,2,2}

NTU finite-patch {2,2} bond environment for the HORIZONTAL center(rc,cc) ->
leaf(rc,cp1) bond. See file header. `coords` overrides the patch list (used by
the loop-content control).
"""
function build_ntu_bondenv(state, X, Y, rc::Int, cc::Int, cp1::Int;
                           patch::Symbol = :P6, boundary::Symbol = :trace,
                           coords = nothing)
    cell = state.unitcell
    Nr, Nc = cell.Ly, cell.Lx
    nsites = Nr * Nc
    peps = measurement_peps(state)

    # --- pre-build asserts: the q-legs must be the reduced DX/DY legs (D-agnostic;
    #     replaces the ==16/==8 hardcode so D=3 q-dims 24/12 pass) ---
    @assert isdual(space(X, 2)) "expected X q-leg (axis 2) to be dual"
    @assert isdual(space(Y, 4)) "expected Y q-leg (axis 4) to be dual"
    # --- self-aliasing guard (finite-patch analog of the E+ aliasing hazard) ---
    # On Nr>=4 rows {rm1,rc,rp1} are 3 distinct rows leaving >=1 row outside, so
    # N-of-rm1 / S-of-rp1 boundary legs point to genuinely external bonds. Nc>=3 so
    # cols {cc,cp1} plus the wrap do not self-alias. P10/P12 widen the columns to
    # {cm1,cc,cp1,cp2}; on a 4x4 those are ALL 4 cols (full-width band) — guard.
    @assert Nr >= 4 && Nc >= 3 "patch needs Nr>=4 && Nc>=3 (self-aliasing guard)"
    if patch in (:P10, :P12) && coords === nothing
        @assert Nc >= 5 "P10/P12 need Nc>=5 so the added arm columns do not self-alias on the wrap (got Nc=$Nc)"
    end

    _prevr(r) = r == 1 ? Nr : r - 1
    _nextr(r) = r == Nr ? 1 : r + 1
    _prevc(c) = c == 1 ? Nc : c - 1
    _nextc(c) = c == Nc ? 1 : c + 1

    # id-allocation skeleton — VERBATIM from build_exact_cluster_bondenv.
    HbondK(r, c) = (r - 1) * Nc + c
    VbondK(r, c) = nsites + (r - 1) * Nc + c
    HbondBra(r, c) = 2 * nsites + HbondK(r, c)
    VbondBra(r, c) = 2 * nsites + VbondK(r, c)
    physTrace(k) = 4 * nsites + k
    # boundary ids start STRICTLY above the physTrace block (4*nsites + k, k<=nsites)
    boundary_start = 5 * nsites + 1
    boundary_counter = Ref(boundary_start - 1)
    next_boundary() = (boundary_counter[] += 1; boundary_counter[])

    cs = coords === nothing ? patch_coords(rc, cc, cp1, Nr, Nc, patch) : coords
    patch_set = Set(cs)

    # neighbour cell in PEPSKit [N,E,S,W] direction (matches the Hbond/Vbond geometry)
    function neighbour(r, c, dir)
        if dir === :N
            return (_prevr(r), c)
        elseif dir === :E
            return (r, _nextc(c))
        elseif dir === :S
            return (_nextr(r), c)
        else # :W
            return (r, _prevc(c))
        end
    end

    # interior id for a virtual leg in [N,E,S,W] direction (the EXACT expressions
    # build_exact_cluster_bondenv uses).
    interiorKetId(r, c, dir) = dir === :N ? VbondK(r, c) :
                               dir === :E ? HbondK(r, c) :
                               dir === :S ? VbondK(_nextr(r), c) :
                               HbondK(r, _prevc(c))                    # :W
    interiorBraId(r, c, dir) = dir === :N ? VbondBra(r, c) :
                               dir === :E ? HbondBra(r, c) :
                               dir === :S ? VbondBra(_nextr(r), c) :
                               HbondBra(r, _prevc(c))                  # :W

    # PEPSKit absorb_weight direction constant for a [N,E,S,W] leg direction.
    dir_const(dir) = dir === :N ? NORTH : dir === :E ? EAST : dir === :S ? SOUTH : WEST

    # boundary-leg id cache: a boundary leg (neighbour outside patch) gets ONE fresh
    # id SHARED between this site's ket leg and its OWN bra leg -> delta trace.
    boundary_ids = Dict{Tuple{Int,Int,Symbol},Int}()
    function boundary_id(r, c, dir)
        get!(boundary_ids, (r, c, dir)) do
            next_boundary()
        end
    end

    # is this virtual leg internal (both endpoints in patch) or boundary?
    is_internal(r, c, dir) = neighbour(r, c, dir) in patch_set

    # the ket/bra id for a NON-q virtual leg: interior id if internal, else a fresh
    # per-site boundary id (same id for ket & bra -> trace).
    function legid(r, c, dir, layer)
        if is_internal(r, c, dir)
            return layer === :ket ? interiorKetId(r, c, dir) : interiorBraId(r, c, dir)
        else
            return boundary_id(r, c, dir)   # SAME id for ket and bra -> delta trace
        end
    end

    Xarr = convert(Array{ComplexF64,4}, convert(Array, X))   # [N,q,S,W]
    Yarr = convert(Array{ComplexF64,4}, convert(Array, Y))   # [N,E,S,q]

    arrays = Vector{Array{ComplexF64}}()
    idxlists = Vector{Vector{Int}}()

    # optional lambda dressing of a BULK boundary leg (bulk PEPSTensor only — NOT the
    # q-reduced X/Y, which lack a PEPSTensor leg structure).
    function maybe_dress(T_tm, r, c, boundary_dirs)
        if boundary === :lambda && !isempty(boundary_dirs)
            Td = T_tm
            for d in boundary_dirs
                Td = absorb_weight(Td, state.weights, r, c, dir_const(d); inv = false)
            end
            return Td
        else
            return T_tm
        end
    end

    k = 0
    for (r, c) in cs
        k += 1
        is_center = (r == rc && c == cc)
        is_leaf = (r == rc && c == cp1)

        if is_center
            # KET = X [N,q,S,W]; q=-1 (DX0 ket). The cut E bond is carried by q. The W
            # leg routes through legid (boundary in P6, interior in P10/P12/torus).
            push!(arrays, Xarr)
            push!(idxlists, [
                legid(rc, cc, :N, :ket),   # N
                -1,                        # q = DX0 OPEN (ket)
                legid(rc, cc, :S, :ket),   # S
                legid(rc, cc, :W, :ket),   # W
            ])
            push!(arrays, conj(Xarr))
            push!(idxlists, [
                legid(rc, cc, :N, :bra),
                -3,                        # q = DX1 OPEN (bra)
                legid(rc, cc, :S, :bra),
                legid(rc, cc, :W, :bra),
            ])
        elseif is_leaf
            # KET = Y [N,E,S,q]; q=-2 (DY0 ket). The cut W bond is carried by q. The E
            # leg routes through legid (boundary in P6, interior in P10/P12/torus).
            push!(arrays, Yarr)
            push!(idxlists, [
                legid(rc, cp1, :N, :ket),  # N
                legid(rc, cp1, :E, :ket),  # E
                legid(rc, cp1, :S, :ket),  # S
                -2,                        # q = DY0 OPEN (ket)
            ])
            push!(arrays, conj(Yarr))
            push!(idxlists, [
                legid(rc, cp1, :N, :bra),
                legid(rc, cp1, :E, :bra),
                legid(rc, cp1, :S, :bra),
                -4,                        # q = DY1 OPEN (bra)
            ])
        else
            # BULK cell: full peps[r,c] [P,N,E,S,W], physical leg TRACED bra<->ket.
            T_tm = peps[r, c]
            bdirs = [d for d in (:N, :E, :S, :W) if !is_internal(r, c, d)]
            T_ket_tm = maybe_dress(T_tm, r, c, bdirs)
            T = convert(Array{ComplexF64,5}, convert(Array, T_ket_tm))  # [P,N,E,S,W]
            push!(arrays, T)
            push!(idxlists, [
                physTrace(k),              # P (TRACED, shared positive id)
                legid(r, c, :N, :ket),
                legid(r, c, :E, :ket),
                legid(r, c, :S, :ket),
                legid(r, c, :W, :ket),
            ])
            push!(arrays, conj(T))         # same lambda on bra (real sqrt(lambda))
            push!(idxlists, [
                physTrace(k),              # P (same shared id -> traced)
                legid(r, c, :N, :bra),
                legid(r, c, :E, :bra),
                legid(r, c, :S, :bra),
                legid(r, c, :W, :bra),
            ])
        end
    end

    # --- structural safety check: every boundary id is used EXACTLY TWICE
    #     (once ket, once the SAME site's bra). A mis-paired boundary trace is
    #     caught here, before the contraction. ---
    idcount = Dict{Int,Int}()
    for idl in idxlists, id in idl
        if id >= boundary_start
            idcount[id] = get(idcount, id, 0) + 1
        end
    end
    for (id, n) in idcount
        @assert n == 2 "boundary id $id used $n times (expected exactly 2: ket+bra)"
    end
    # phys-trace / boundary-id non-collision (the gap is exactly 0; assert it holds)
    @assert boundary_start > 4 * nsites + nsites "boundary ids collide with physTrace block"

    # === ASSEMBLY TAIL — BYTE-FOR-BYTE from build_exact_cluster_bondenv ===
    raw = convert(Array, TensorKitNcon.ncon(arrays, idxlists))
    raw_reordered = permutedims(raw, (3, 4, 1, 2))   # (-3,-4,-1,-2) = [DX1 DY1; DX0 DY0]
    DXn = space(X, 2)'    # non-dual
    DYn = space(Y, 4)'    # non-dual
    benv = TensorMap(raw_reordered, DXn ⊗ DYn, DXn ⊗ DYn)
    normalize!(benv, Inf)
    return benv
end
