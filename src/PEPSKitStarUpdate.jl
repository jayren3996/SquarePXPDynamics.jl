module PEPSKitStarUpdate

# PEPSKit-native five-site square-star PXP gate and simple update. The star term
# is a genuine cross  P_right P_up P_left P_down X_center  in site order
# (center, right, up, left, down), basis 1 = excited/up, 2 = down/vacuum.

using TensorKit
using PEPSKit
using LinearAlgebra

using ..SquareGeometry: SquareCoord
using ..SquareUnitCells: PeriodicSquareUnitCell, BondKey, wrap, neighbor, bondkey
using ..SquarePXP:
    SQUARE_STAR_SITES,
    square_pxp_star_hamiltonian,
    square_pxp_gate,
    projected_square_pxp_gate
using ..PEPSKitBackend:
    PXPIPEPSState,
    squarecoord_to_cartesianindex,
    mark_mutated!,
    add_log_norm!,
    _dense_operator_tensor_map

export pxp_star_hamiltonian_matrix, pxp_star_gate_tensormap
export pxp_star_gate_dense
export StarUpdateInfoPK, project_star_pepskit!

"""
    StarUpdateInfoPK

Diagnostics from one PEPSKit-native five-site square-star simple update. The
dicts are keyed by center-to-leaf directions `:right`, `:up`, `:left`, `:down`.
"""
struct StarUpdateInfoPK
    center::SquareCoord
    max_truncerr::Float64
    truncerrs::Dict{Symbol, Float64}
    keptdims::Dict{Symbol, Int}
    min_lambda::Dict{Symbol, Float64}
    norm_factors::Dict{Symbol, Float64}
end

function project_star_pepskit! end

# ---------------------------------------------------------------------------
# Milestone 3: gate as a TensorMap
# ---------------------------------------------------------------------------
# `_dense_operator_tensor_map` is shared from PEPSKitBackend.

"""
    pxp_star_hamiltonian_matrix()

Return the dense 32x32 five-site square-star PXP Hamiltonian
`X_center * P_down_right * P_down_up * P_down_left * P_down_down` in site order
`(center, right, up, left, down)`. Identical to the legacy
`square_pxp_star_hamiltonian`.
"""
pxp_star_hamiltonian_matrix() = square_pxp_star_hamiltonian()

"""
    pxp_star_gate_dense(dt; projected = true, evolution = :real)

Return the dense 32x32 square-star PXP evolution gate. `evolution = :real` gives
`exp(-im*dt*H)`, `:imaginary` gives `exp(-dt*H)`. When `projected`, the local
blockade projector is applied (`P*U`).
"""
function pxp_star_gate_dense(dt::Real; projected::Bool = true, evolution::Symbol = :real)
    return projected ? projected_square_pxp_gate(dt; evolution) :
           square_pxp_gate(dt; evolution)
end

"""
    pxp_star_gate_tensormap(dt; projected = true, evolution = :real)

Return the five-site square-star PXP evolution gate as a TensorMap on
`(P)^5 ← (P)^5` with physical legs in order `(center, right, up, left, down)`,
`P = ComplexSpace(2)`. Real time uses `exp(-im*dt*H)`, imaginary time uses
`exp(-dt*H)`; `projected` applies the local blockade projector.
"""
function pxp_star_gate_tensormap(
    dt::Real;
    projected::Bool = true,
    evolution::Symbol = :real,
)
    gate = pxp_star_gate_dense(dt; projected, evolution)
    return _dense_operator_tensor_map(gate, SQUARE_STAR_SITES)
end

# ---------------------------------------------------------------------------
# Milestone 4: PEPSKit-native star simple update
# ---------------------------------------------------------------------------
# The five-site square-star simple update operates on InfinitePEPS + SUWeight in
# PEPSKit's native simple-update gauge (stored Γ tensors carry sqrt(λ) on each
# bond).
#
# Algorithm (mirrors the legacy ITensors `project_star!` but on PEPSKit objects
# and PEPSKit's gauge / absorb_weight conventions):
#   1. Locate center + right/up/left/down leaves in PEPSKit (row, col) coords.
#   2. Absorb sqrt of each leaf's three *environment* bond weights (everything
#      except the center-facing bond), normalise.
#   3. QR-reduce each leaf: split the 3 environment legs (X) from the reduced
#      core (physical + center-facing bond + internal q).
#   4. Contract center with the 4 reduced cores over the active bonds and apply
#      the 5-site star gate on the physical legs.
#   5. Sequentially SVD-split the reduced object into center + 4 leaf cores,
#      truncating each new center-leaf bond; absorb sqrt(s) symmetrically.
#   6. Reattach the QR isometries, remove the environment weights, write the new
#      InfinitePEPS tensors and SUWeight bonds.
#
# PEPSKit axis convention: 1,2,3,4 = N,E,S,W. TensorKit PEPS leg order:
# 1 = P (physical), 2 = N, 3 = E, 4 = S, 5 = W.

# StarUpdateInfoPK is defined in the enclosing module before this include.

const _STAR_DIRS = (:right, :up, :left, :down)

_nextidx(i, n) = mod1(i + 1, n)
_previdx(i, n) = mod1(i - 1, n)

# Per-direction static geometry, from the center's point of view.
#   leaf_offset    : (drow, dcol) from center to leaf in PEPSKit coords
#   cface_leg      : leaf TensorKit leg index facing the center (2=N,3=E,4=S,5=W)
#   env_legs       : leaf's three TensorKit env legs (the other PEPS legs)
#   env_ax         : PEPSKit absorb axes (1..4) for those env bonds
#   center_leg     : center TensorKit leg index toward this leaf
function _dir_spec(dir::Symbol)
    if dir === :right          # leaf east of center; faces via its WEST leg
        return (offset = (0, 1), cface_leg = 5, env_legs = (2, 3, 4),
            env_ax = (1, 2, 3), center_leg = 3)
    elseif dir === :up         # leaf north of center; faces via its SOUTH leg
        return (offset = (-1, 0), cface_leg = 4, env_legs = (2, 3, 5),
            env_ax = (1, 2, 4), center_leg = 2)
    elseif dir === :left       # leaf west of center; faces via its EAST leg
        return (offset = (0, -1), cface_leg = 3, env_legs = (2, 4, 5),
            env_ax = (1, 3, 4), center_leg = 5)
    elseif dir === :down       # leaf south of center; faces via its NORTH leg
        return (offset = (1, 0), cface_leg = 2, env_legs = (3, 4, 5),
            env_ax = (2, 3, 4), center_leg = 4)
    else
        throw(ArgumentError("direction must be :right, :up, :left, or :down"))
    end
end

function _leaf_rowcol(rc::Int, cc::Int, Nr::Int, Nc::Int, dir::Symbol)
    o = _dir_spec(dir).offset
    return (mod1(rc + o[1], Nr), mod1(cc + o[2], Nc))
end

# SUWeight slot (axis d in 1:2, row, col) holding the center-leaf bond for `dir`.
function _bond_slot(rc::Int, cc::Int, Nr::Int, Nc::Int, dir::Symbol)
    if dir === :right
        return (1, rc, cc)                    # east bond of center
    elseif dir === :up
        return (2, rc, cc)                    # north bond of center
    elseif dir === :left
        return (1, rc, _previdx(cc, Nc))      # west bond of center == east of left leaf
    else # :down
        return (2, _nextidx(rc, Nr), cc)      # south bond of center == north of down leaf
    end
end

function _validate_finite_step(step::Real)
    isfinite(step) || throw(ArgumentError("dt must be finite"))
    return step
end

function _validate_maxdim(maxdim::Integer)
    maxdim >= 1 || throw(ArgumentError("maxdim must be at least 1"))
    return Int(maxdim)
end

function _validate_split_order(split_order)
    order = Tuple(split_order)
    (length(order) == 4 && Set(order) == Set(_STAR_DIRS)) ||
        throw(ArgumentError("split_order must be a permutation of $(_STAR_DIRS)"))
    return order
end

function _assert_distinct_star(rc, cc, Nr, Nc)
    coords = Set{Tuple{Int, Int}}()
    push!(coords, (rc, cc))
    for dir in _STAR_DIRS
        push!(coords, _leaf_rowcol(rc, cc, Nr, Nc, dir))
    end
    length(coords) == 5 || throw(
        ArgumentError(
            "wrapped square star must contain five distinct unit-cell sites; " *
            "increase the unit cell (need Lx, Ly >= 3)",
        ),
    )
    return nothing
end

# Truncation strategy for one bond split.
function _bond_trunc(maxdim::Int, cutoff::Float64)
    if cutoff > 0
        return truncdim(maxdim) & truncbelow(cutoff)
    else
        return truncdim(maxdim)
    end
end

# tsvd of `t` (codomain = the two legs to split off), returning the leaf-active
# tensor (carrying sqrt(s) on the new bond), the new core (sqrt(s) on the new
# bond too), the non-dual weight s, and the truncation error.
function _split_bond(t, maxdim::Int, cutoff::Float64, rel_floor::Float64)
    # full spectrum (untruncated) for truncation-error bookkeeping
    _, Sfull, _ = tsvd(t)
    full_sq = sum(abs2, Sfull.data)
    U, S, V = tsvd(t; trunc = _bond_trunc(maxdim, cutoff))
    if rel_floor > 0
        svals = S.data
        if !isempty(svals)
            σmax = maximum(svals)
            if σmax > 0
                keep = max(1, count(>=(rel_floor * σmax), svals))
                if keep < length(svals)
                    U, S, V = tsvd(t; trunc = truncdim(keep))
                end
            end
        end
    end
    kept_sq = sum(abs2, S.data)
    ϵ = sqrt(max(0.0, full_sq - kept_sq)) / sqrt(max(full_sq, eps(Float64)))
    # Ensure the stored weight space is non-dual (SUWeight requirement).
    if isdual(space(S, 1))
        U, S, V = PEPSKit.flip_svd(U, S, V)
    end
    leaf_active, new_core = PEPSKit.absorb_s(U, S, V)
    return leaf_active, new_core, S, ϵ
end

"""
    project_star_pepskit!(state, center; dt, maxdim, cutoff,
                          rel_floor = 1e-3, evolution = :real, projected = true,
                          split_order = (:right, :up, :left, :down))

Apply one PEPSKit-native five-site square-star simple update at `center`
(a `SquareCoord`). Mutates `state.peps` and `state.weights` in place and returns
a `StarUpdateInfoPK` with per-direction truncation errors, kept dimensions,
minimum weights, and norm factors. Only this single star is updated.
"""
function project_star_pepskit!(
    state::PXPIPEPSState,
    center::SquareCoord;
    dt::Real,
    maxdim::Integer,
    cutoff::Real = 1e-12,
    rel_floor::Real = 1.0e-3,
    evolution::Symbol = :real,
    projected::Bool = true,
    split_order = _STAR_DIRS,
)::StarUpdateInfoPK
    step = _validate_finite_step(dt)
    kept_maxdim = _validate_maxdim(maxdim)
    trunc_cutoff = Float64(cutoff)
    trunc_rel_floor = Float64(rel_floor)
    order = _validate_split_order(split_order)

    cell = state.unitcell
    peps = state.peps
    weights = state.weights
    Nr, Nc = size(peps)
    cidx = squarecoord_to_cartesianindex(cell, center)
    rc, cc = cidx[1], cidx[2]
    _assert_distinct_star(rc, cc, Nr, Nc)

    gate = pxp_star_gate_tensormap(step; projected, evolution)

    # --- absorb env weights + QR-reduce each leaf ---
    Qfactors = Dict{Symbol, Any}()
    reduced = Dict{Symbol, Any}()
    leaf_rowcol = Dict{Symbol, Tuple{Int, Int}}()
    for dir in _STAR_DIRS
        spec = _dir_spec(dir)
        (lr, lc) = _leaf_rowcol(rc, cc, Nr, Nc, dir)
        leaf_rowcol[dir] = (lr, lc)
        T = peps.A[lr, lc]
        T = PEPSKit.absorb_weight(T, weights, lr, lc, spec.env_ax; inv = false)
        normalize!(T, Inf)
        # codomain = env legs, domain = (physical, center-facing bond)
        Tp = permute(T, (spec.env_legs, (1, spec.cface_leg)))
        Q, R = leftorth(Tp)
        Qfactors[dir] = Q          # env legs <- q
        reduced[dir] = R           # q <- (P, cbond)
    end

    # --- contract center with reduced cores, apply gate ---
    C = peps.A[rc, cc]
    aR, aU, aL, aD = reduced[:right], reduced[:up], reduced[:left], reduced[:down]
    # C legs: P(-1), N=cbond_up(2), E=cbond_right(1), S=cbond_down(4), W=cbond_left(3)
    Cidx = [-1, 2, 1, 4, 3]
    θpre = ncon(
        (C, aR, aU, aL, aD),
        (Cidx, [-6, -2, 1], [-7, -3, 2], [-8, -4, 3], [-9, -5, 4]),
    )
    # apply gate on the 5 physical legs (center, right, up, left, down)
    θ = ncon((gate, θpre), ([-1, -2, -3, -4, -5, 1, 2, 3, 4, 5], [1, 2, 3, 4, 5, -6, -7, -8, -9]))

    # θ legs (linear order): 1=Pc 2=Pr 3=Pu 4=Pl 5=Pd 6=qr 7=qu 8=ql 9=qd
    legs = [:Pc, :Pr, :Pu, :Pl, :Pd, :qr, :qu, :ql, :qd]
    phys_of = Dict(:right => :Pr, :up => :Pu, :left => :Pl, :down => :Pd)
    q_of = Dict(:right => :qr, :up => :qu, :left => :ql, :down => :qd)
    bond_of = Dict(:right => :bR, :up => :bU, :left => :bL, :down => :bD)

    core = θ
    leaf_active = Dict{Symbol, Any}()
    new_weights = Dict{Symbol, Any}()
    truncerrs = Dict{Symbol, Float64}()
    keptdims = Dict{Symbol, Int}()
    min_lambda = Dict{Symbol, Float64}()
    norm_factors = Dict{Symbol, Float64}()

    for dir in order
        # move (phys_dir, q_dir) into codomain, everything else into domain
        pname, qname = phys_of[dir], q_of[dir]
        pi = findfirst(==(pname), legs)
        qi = findfirst(==(qname), legs)
        rest_positions = [i for i in eachindex(legs) if i != pi && i != qi]
        core_perm = permute(core, ((pi, qi), Tuple(rest_positions)))
        la, new_core, S, ϵ = _split_bond(core_perm, kept_maxdim, trunc_cutoff, trunc_rel_floor)
        leaf_active[dir] = la
        new_weights[dir] = S
        svals = S.data
        scale = norm(svals)
        truncerrs[dir] = ϵ
        keptdims[dir] = length(svals)
        min_lambda[dir] = isempty(svals) ? 0.0 : minimum(svals)
        norm_factors[dir] = scale
        # new_core legs: [bond_dir, rest...] (bond in codomain from sqrt(S)*V)
        core = new_core
        legs = vcat([bond_of[dir]], legs[rest_positions])
    end

    # Now `core` is the new center. Its legs are the 4 new bonds + Pc.
    # Assemble center as P <- N(bU) E(bR) S(bD) W(bL).
    posC = Dict(name => findfirst(==(name), legs) for name in legs)
    new_center = permute(
        core,
        ((posC[:Pc],), (posC[:bU], posC[:bR], posC[:bD], posC[:bL])),
    )

    # --- reattach QR isometries, remove env weights, write back ---
    add_log = 0.0
    peps.A[rc, cc] = new_center
    for dir in _STAR_DIRS
        spec = _dir_spec(dir)
        (lr, lc) = leaf_rowcol[dir]
        Q = Qfactors[dir]
        la = leaf_active[dir]
        # la legs: (P_dir, q_dir, bond_dir); Q legs: (env1, env2, env3, q_dir)
        # contract over q -> full leaf with legs (env1,env2,env3, P_dir, bond_dir)
        full = ncon((Q, la), ([-1, -2, -3, 1], [-4, 1, -5]))
        # full linear legs: env_legs[1],env_legs[2],env_legs[3], P, cbond
        # Build a PEPS tensor P <- N E S W: place cbond at cface_leg, env at others.
        # Map: we know spec.env_legs (target PEPS leg ids) and spec.cface_leg.
        leg_targets = (spec.env_legs[1], spec.env_legs[2], spec.env_legs[3], 1, spec.cface_leg)
        # current linear order of `full`: positions 1..5 correspond to leg_targets
        # Find permutation so output order is PEPS legs (1=P,2=N,3=E,4=S,5=W).
        order_pos = [findfirst(==(target), leg_targets) for target in (1, 2, 3, 4, 5)]
        leaf_tensor = permute(full, ((order_pos[1],), Tuple(order_pos[2:5])))
        # remove env weights (inv) that we absorbed earlier
        leaf_tensor = PEPSKit.absorb_weight(leaf_tensor, weights, lr, lc, spec.env_ax; inv = true)
        peps.A[lr, lc] = leaf_tensor
        # write new bond weight
        slot = _bond_slot(rc, cc, Nr, Nc, dir)
        weights.data[slot[1], slot[2], slot[3]] = new_weights[dir]
        add_log += log(max(norm_factors[dir], eps(Float64)))
    end

    add_log_norm!(state, add_log)
    mark_mutated!(state)

    max_truncerr = maximum(values(truncerrs))
    return StarUpdateInfoPK(center, max_truncerr, truncerrs, keptdims, min_lambda, norm_factors)
end

end
