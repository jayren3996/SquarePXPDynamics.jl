module PEPSKitBackend

# PEPSKit-native iPEPS state for square PXP dynamics. PEPSKit owns the tensor
# network (InfinitePEPS + SUWeight); this wrapper adds only the bookkeeping the
# repo needs (unit cell, mutation version, log-norm ledger) and the PXP-specific
# product/checkerboard constructors. No site tensors, link indices, or link
# weights are duplicated outside PEPSKit.

using TensorKit
using PEPSKit

using ..SquareGeometry: SquareCoord
using ..SquareUnitCells: PeriodicSquareUnitCell, wrap, neighbor
# Extend the existing generics so dispatch on PXPIPEPSState coexists with the
# legacy SquareIPEPSState methods (single generic, no name clash).
import ..SquareIPEPS: copy_state, state_version, log_norm

export PXPIPEPSState
export product_pxp_pepskit_state, checkerboard_pxp_pepskit_state
export copy_state, state_version, log_norm, mark_mutated!, add_log_norm!
export pepskit_peps, pepskit_weights
export squarecoord_to_cartesianindex, cartesianindex_to_squarecoord
export measurement_peps

# PEPSKit axis convention for absorb_weight and PEPS domain order: 1..4 = N,E,S,W.

"""
    PXPIPEPSState

Lightweight wrapper bundling a PEPSKit `InfinitePEPS` (bare Γ tensors) and its
`SUWeight` (the λ ledger) for square PXP dynamics, in Vidal/Γ-λ simple-update
form. `unitcell` is the repo's `PeriodicSquareUnitCell`. `mutation_version`
counts state-changing operations so cached CTMRG contexts can detect staleness;
`log_norm_value` accumulates the diagnostic log-normalization bookkeeping from
simple-update splits (not a physical observable).

This wrapper is for bookkeeping only. The site tensors and link weights live
inside `peps` and `weights`; nothing is duplicated here.
"""
struct PXPIPEPSState
    peps::PEPSKit.InfinitePEPS
    weights::PEPSKit.SUWeight
    unitcell::PeriodicSquareUnitCell
    mutation_version::Base.RefValue{Int}
    log_norm_value::Base.RefValue{Float64}
end

"""
    pepskit_peps(state)

Return the underlying PEPSKit `InfinitePEPS` of `state`.
"""
pepskit_peps(state::PXPIPEPSState) = state.peps

"""
    pepskit_weights(state)

Return the underlying PEPSKit `SUWeight` (the λ ledger) of `state`.
"""
pepskit_weights(state::PXPIPEPSState) = state.weights

"""
    state_version(state)

Return the mutation counter for `state`. Operations that change tensors or
stored link weights increment this counter so cached measurement contexts can
detect stale states.
"""
state_version(state::PXPIPEPSState)::Int = state.mutation_version[]

"""
    log_norm(state)

Return the accumulated logarithmic normalization ledger recorded during
simple-update star splits. Diagnostic bookkeeping, not a physical observable.
"""
log_norm(state::PXPIPEPSState)::Float64 = state.log_norm_value[]

"""
    mark_mutated!(state)

Increment the mutation counter of `state`.
"""
function mark_mutated!(state::PXPIPEPSState)
    state.mutation_version[] += 1
    return state
end

"""
    add_log_norm!(state, delta)

Add the finite increment `delta` to the log-normalization ledger of `state`.
"""
function add_log_norm!(state::PXPIPEPSState, delta::Real)
    value = Float64(delta)
    isfinite(value) || throw(ArgumentError("log-norm increment must be finite"))
    state.log_norm_value[] += value
    return state
end

"""
    copy_state(state)

Return a deep mutable copy of `state` with independent PEPSKit tensors, weights,
mutation counter, and log-normalization ledger.
"""
function copy_state(state::PXPIPEPSState)::PXPIPEPSState
    return PXPIPEPSState(
        deepcopy(state.peps),
        deepcopy(state.weights),
        state.unitcell,
        Ref(state_version(state)),
        Ref(log_norm(state)),
    )
end

# --- coordinate mapping (identical to the legacy to_pepskit adapter) ---

"""
    squarecoord_to_cartesianindex(cell, c)

Map `SquareCoord(x, y)` to PEPSKit matrix coordinate `CartesianIndex(row, col)`
with `row = Ly - y + 1`, `col = x` (after periodic wrapping). The custom `:up`
neighbour maps to PEPSKit north (`row - 1`); `:right` maps to PEPSKit east
(`col + 1`).
"""
function squarecoord_to_cartesianindex(cell::PeriodicSquareUnitCell, c::SquareCoord)
    site = wrap(cell, c)
    return CartesianIndex(cell.Ly - site.y + 1, site.x)
end

"""
    cartesianindex_to_squarecoord(cell, idx)

Inverse of [`squarecoord_to_cartesianindex`](@ref) for in-cell indices.
"""
function cartesianindex_to_squarecoord(cell::PeriodicSquareUnitCell, idx::CartesianIndex{2})
    row, col = idx[1], idx[2]
    return SquareCoord(col, cell.Ly - row + 1)
end

# --- dense PXP operator -> TensorMap (shared by the native modules) ---

# Big-endian dense computational-basis index for a tuple of 1/2 values (site 1
# most significant). Matches kron_all / square_pxp_star_hamiltonian ordering.
function _dense_index(values)
    idx = 1
    n = length(values)
    for (site, value) in enumerate(values)
        1 <= value <= 2 || throw(ArgumentError("basis values must be 1 or 2"))
        idx += (value - 1) * 2^(n - site)
    end
    return idx
end

# Convert a dense (2^n x 2^n) operator into a TensorMap on (P)^n <- (P)^n with
# physical legs in the same order as the dense kron ordering. Shared by the
# native star gate (PEPSKitStarUpdate) and the native observables
# (PEPSKitPXPObservables).
function _dense_operator_tensor_map(O::AbstractMatrix, nsites::Int)
    size(O) == (2^nsites, 2^nsites) ||
        throw(ArgumentError("dense operator must be $(2^nsites)x$(2^nsites)"))
    p = ComplexSpace(2)
    data = zeros(ComplexF64, ntuple(Returns(2), 2nsites))
    rng = ntuple(_ -> 1:2, nsites)
    for out_values in Iterators.product(rng...)
        out_idx = _dense_index(out_values)
        for in_values in Iterators.product(rng...)
            in_idx = _dense_index(in_values)
            data[out_values..., in_values...] = O[out_idx, in_idx]
        end
    end
    spaces = ntuple(_ -> p, nsites)
    product_space = reduce(⊗, spaces)
    return TensorMap(data, product_space ← product_space)
end

# --- construction ---

function _validate_cell(cell::PeriodicSquareUnitCell)
    cell.Lx >= 2 || throw(ArgumentError("iPEPS construction requires Lx >= 2"))
    cell.Ly >= 2 || throw(ArgumentError("iPEPS construction requires Ly >= 2"))
    return cell
end

function _validate_bonddim(D::Integer)
    D >= 1 || throw(ArgumentError("D must be at least 1"))
    return Int(D)
end

function _product_amplitudes(state::Symbol)
    if state in (:up, :z_up)
        return ComplexF64[1, 0]      # basis 1 = excited/up
    elseif state in (:down, :z_down)
        return ComplexF64[0, 1]      # basis 2 = down/vacuum
    elseif state === :x_plus
        return ComplexF64[inv(sqrt(2)), inv(sqrt(2))]
    else
        throw(ArgumentError("state must be :up, :down, :z_up, :z_down, or :x_plus"))
    end
end

# Build a single PEPS Γ tensor for a product amplitude vector at bond dim D.
# Domain order N,E,S,W with N,E non-dual and S,W dual (PEPSKit convention).
# Only the leading virtual index (1) carries the product state; higher virtual
# components are zero, matching the leading-1 bond weight.
function _product_peps_tensor(amplitudes::AbstractVector, D::Int)
    P = ComplexSpace(2)
    V = ComplexSpace(D)
    data = zeros(ComplexF64, 2, D, D, D, D)
    for s in 1:2
        data[s, 1, 1, 1, 1] = amplitudes[s]
    end
    return TensorMap(data, P ← V ⊗ V ⊗ V' ⊗ V')
end

# Leading-1 Schmidt weight of length D, as a non-dual DiagonalTensorMap.
function _leading_weight(D::Int)
    vals = zeros(Float64, D)
    vals[1] = 1.0
    return DiagonalTensorMap(vals, ComplexSpace(D))
end

# Build the SUWeight whose [1,r,c] (east) and [2,r,c] (north) entries are the
# leading-1 product weights. With D=1 this is exactly the identity weight.
function _product_weights(Nr::Int, Nc::Int, D::Int)
    data = Array{DiagonalTensorMap{Float64, ComplexSpace, Vector{Float64}}, 3}(undef, 2, Nr, Nc)
    for d in 1:2, r in 1:Nr, c in 1:Nc
        data[d, r, c] = _leading_weight(D)
    end
    return SUWeight(data)
end

function _build_product_state(cell::PeriodicSquareUnitCell, amplitude_at, D::Int)
    Nr, Nc = cell.Ly, cell.Lx
    tensors = Matrix{TensorMap{ComplexF64, ComplexSpace, 1, 4, Vector{ComplexF64}}}(undef, Nr, Nc)
    for row in 1:Nr, col in 1:Nc
        c = cartesianindex_to_squarecoord(cell, CartesianIndex(row, col))
        tensors[row, col] = _product_peps_tensor(amplitude_at(c), D)
    end
    peps = InfinitePEPS(tensors)
    weights = _product_weights(Nr, Nc, D)
    return PXPIPEPSState(peps, weights, cell, Ref(0), Ref(0.0))
end

"""
    product_pxp_pepskit_state(cell; state = :down, D = 1)

Construct a PEPSKit-native square iPEPS product state on `cell`. Physical basis
`:up`/`:z_up` is index `1` (excited), `:down`/`:z_down` is index `2` (vacuum),
`:x_plus` is equal amplitudes. Bond dimension `D` defaults to `1`.
"""
function product_pxp_pepskit_state(
    cell::PeriodicSquareUnitCell;
    state::Symbol = :down,
    D::Integer = 1,
)
    _validate_cell(cell)
    d = _validate_bonddim(D)
    amplitudes = _product_amplitudes(state)
    return _build_product_state(cell, Returns(amplitudes), d)
end

"""
    checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = 1)

Construct a PEPSKit-native square iPEPS checkerboard product state. Even parity
is `iseven(c.x + c.y)`; excited sites use physical basis index `1`.
"""
function checkerboard_pxp_pepskit_state(
    cell::PeriodicSquareUnitCell;
    excited_on::Symbol = :even,
    D::Integer = 1,
)
    _validate_cell(cell)
    d = _validate_bonddim(D)
    excited_on === :even ||
        excited_on === :odd ||
        throw(ArgumentError("excited_on must be :even or :odd"))
    up = ComplexF64[1, 0]
    down = ComplexF64[0, 1]
    function amplitude_at(c::SquareCoord)
        even = iseven(c.x + c.y)
        excited = excited_on === :even ? even : !even
        return excited ? up : down
    end
    return _build_product_state(cell, amplitude_at, d)
end

# --- measurement helpers ---

"""
    measurement_peps(state)

Return the `InfinitePEPS` ready for CTMRG. The backend uses PEPSKit's native
simple-update gauge (convention "B"): every stored Γ tensor already carries
`sqrt(λ)` on each bond, so neighbouring tensors reproduce the full bond weight λ
and the `InfinitePEPS` is the measurement state directly. The SUWeight ledger is
used only by the simple update, not for measurement. `state` is not mutated.
"""
measurement_peps(state::PXPIPEPSState)::PEPSKit.InfinitePEPS = state.peps

end
