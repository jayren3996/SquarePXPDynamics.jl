module GaugeDiagnostics

using ITensors
using LinearAlgebra

using ..SquareGeometry: SquareCoord
using ..SquareUnitCells: BondKey, bondkey, neighbor
using ..SquareIPEPS:
    SquareIPEPSState, physical_index, link_index, absorb_link_weight

export SimpleGaugeDiagnostic
export gauge_diagnostic_simple, gauge_deviation_simple, all_gauge_deviations_simple

const _DIRECTIONS = (:right, :up, :left, :down)

"""
    SimpleGaugeDiagnostic

Read-only diagnostic summary for one canonical simple-gauge bond. `deviation`
is the relative Frobenius norm of the off-diagonal part of the local bond norm
matrix, while the diagonal fields summarize its real diagonal entries.
"""
struct SimpleGaugeDiagnostic
    bond::BondKey
    deviation::Float64
    frobenius_norm::Float64
    diagonal_min::Float64
    diagonal_max::Float64
    diagonal_condition::Float64
end

function _require_simple_gauge(psi::SquareIPEPSState)
    psi.gauge === :simple ||
        throw(ArgumentError("simple-gauge diagnostics require psi.gauge === :simple"))
    return psi
end

function _validate_direction(dir::Symbol)
    dir in _DIRECTIONS ||
        throw(ArgumentError("direction must be :right, :up, :left, or :down"))
    return dir
end

function _opposite_direction(dir::Symbol)
    if dir === :right
        return :left
    elseif dir === :left
        return :right
    elseif dir === :up
        return :down
    elseif dir === :down
        return :up
    else
        throw(ArgumentError("direction must be :right, :up, :left, or :down"))
    end
end

function _external_absorbed_tensor(
    psi::SquareIPEPSState,
    c::SquareCoord,
    target_dir::Symbol,
)
    T = psi.tensors[c]
    for dir in _DIRECTIONS
        if dir !== target_dir
            T = absorb_link_weight(T, psi, c, dir)
        end
    end
    return T
end

function _endpoint_gram(psi::SquareIPEPSState, c::SquareCoord, target_dir::Symbol)
    T = _external_absorbed_tensor(psi, c, target_dir)
    target = link_index(psi, c, target_dir)
    physical = physical_index(psi, c)
    hasind(T, target) || throw(ArgumentError("endpoint tensor is missing target link index"))
    hasind(T, physical) || throw(ArgumentError("endpoint tensor is missing physical index"))

    gram = dag(prime(T, target)) * T
    return Matrix(Array(gram, prime(target), target))
end

"""
    gauge_diagnostic_simple(psi, c, dir)::SimpleGaugeDiagnostic

Compute a read-only simple-gauge diagnostic for the nearest-neighbor bond from
coordinate `c` in direction `dir`. The requested bond is canonicalized with
`bondkey(psi.unitcell, c, dir)`, the target bond's link weight is not absorbed,
and all three external link weights at each endpoint are absorbed into local
endpoint tensors before forming the bond norm matrix.

Only `psi.gauge === :simple` states are accepted. `dir` must be one of
`:right`, `:up`, `:left`, or `:down`.
"""
function gauge_diagnostic_simple(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
)::SimpleGaugeDiagnostic
    _require_simple_gauge(psi)
    _validate_direction(dir)

    key = bondkey(psi.unitcell, c, dir)
    canonical_site = key.site
    canonical_dir = key.dir
    neighbor_site = neighbor(psi.unitcell, canonical_site, canonical_dir)
    neighbor_dir = _opposite_direction(canonical_dir)

    GA = _endpoint_gram(psi, canonical_site, canonical_dir)
    GB = _endpoint_gram(psi, neighbor_site, neighbor_dir)
    size(GA) == size(GB) ||
        throw(ArgumentError("endpoint Gram matrices have incompatible dimensions"))

    N = GA .* GB
    N = (N + N') / 2
    frobenius_norm = Float64(norm(N))
    isfinite(frobenius_norm) && frobenius_norm > 0 ||
        throw(ArgumentError("local bond norm matrix must have finite nonzero norm"))

    diagonal = diag(N)
    diagonal_values = Float64.(real.(diagonal))
    diagonal_min = minimum(diagonal_values)
    diagonal_max = maximum(diagonal_values)
    diagonal_condition =
        diagonal_min > 0 ? diagonal_max / diagonal_min : Inf
    deviation = Float64(norm(N - Diagonal(diagonal)) / frobenius_norm)

    return SimpleGaugeDiagnostic(
        key,
        deviation,
        frobenius_norm,
        diagonal_min,
        diagonal_max,
        diagonal_condition,
    )
end

"""
    gauge_deviation_simple(psi, c, dir)::Float64

Return only the off-diagonal relative deviation from
[`gauge_diagnostic_simple`](@ref) for the nearest-neighbor bond from `c` in
direction `dir`.
"""
function gauge_deviation_simple(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
)::Float64
    return gauge_diagnostic_simple(psi, c, dir).deviation
end

"""
    all_gauge_deviations_simple(psi)::Dict{BondKey,Float64}

Return a dictionary containing one simple-gauge deviation for every canonical
bond stored in `psi.link_weights`. The state must satisfy `psi.gauge ===
:simple`.
"""
function all_gauge_deviations_simple(psi::SquareIPEPSState)::Dict{BondKey,Float64}
    _require_simple_gauge(psi)
    return Dict{BondKey,Float64}(
        key => gauge_deviation_simple(psi, key.site, key.dir) for
        key in keys(psi.link_weights)
    )
end

end
