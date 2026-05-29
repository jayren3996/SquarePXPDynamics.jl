module IPEPSCompression

using ITensors
using LinearAlgebra: norm

using ..SquareGeometry: SquareCoord
using ..SquareUnitCells: PeriodicSquareUnitCell, BondKey, neighbor, bondkey
using ..SquareIPEPS:
    SquareIPEPSState,
    physical_index,
    link_index,
    link_weight,
    absorb_link_weight,
    deabsorb_link_weight,
    _mark_mutated!,
    _add_log_norm!

export IPEPSCompressionInfo, compress_to_target_maxdim!

const _PRIMARY_BOND_DIRECTIONS = (:right, :up)
const _OPPOSITE_DIRECTION = Dict(
    :right => :left,
    :left => :right,
    :up => :down,
    :down => :up,
)
const _ALL_DIRECTIONS = (:right, :up, :left, :down)

"""
    IPEPSCompressionInfo

Diagnostics from a [`compress_to_target_maxdim!`](@ref) call. `target_maxdim`
is the requested per-bond cap. `bond_truncerrs` maps `BondKey` to the SVD
truncation error returned by the per-bond compression. `bond_keptdims` maps
`BondKey` to the kept singular-value count (`<= target_maxdim`).
`max_truncerr` / `mean_truncerr` summarize over all compressed bonds.
`log_norm_delta` is the sum of log normalization scales accumulated by the
per-bond SVDs (added to `psi.log_norm_value` by the call).
"""
struct IPEPSCompressionInfo
    target_maxdim::Int
    max_truncerr::Float64
    mean_truncerr::Float64
    bond_truncerrs::Dict{BondKey,Float64}
    bond_keptdims::Dict{BondKey,Int}
    log_norm_delta::Float64
end

function _singular_values_from_S(S::ITensor)
    indices = inds(S)
    length(indices) == 2 ||
        throw(ArgumentError("singular-value tensor must have two indices"))
    a, b = indices
    n = min(dim(a), dim(b))
    values = Vector{Float64}(undef, n)
    for i in 1:n
        values[i] = abs(S[a => i, b => i])
    end
    return values
end

function _normalize_singular_values(svals::Vector{Float64})
    !isempty(svals) ||
        throw(ArgumentError("empty singular spectrum encountered during compression"))
    all(isfinite, svals) ||
        throw(ArgumentError("singular values must be finite during compression"))
    all(>=(0), svals) ||
        throw(ArgumentError("singular values must be nonnegative during compression"))
    scale = norm(svals)
    isfinite(scale) && scale > 0 || throw(
        ArgumentError(
            "zero or non-finite singular spectrum encountered during compression",
        ),
    )
    lambda = svals ./ scale
    return lambda, scale
end

function _site_order_indices(psi::SquareIPEPSState, c::SquareCoord, link_lookup)
    p = physical_index(psi, c)
    left = link_lookup(:left)
    right = link_lookup(:right)
    up = link_lookup(:up)
    down = link_lookup(:down)
    return (p, left, right, up, down)
end

function _compress_bond!(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
    target_maxdim::Int,
    cutoff::Float64,
)
    c_other = neighbor(psi.unitcell, c, dir)
    opp_dir = _OPPOSITE_DIRECTION[dir]
    other_dirs_c = filter(d -> d !== dir, _ALL_DIRECTIONS)
    other_dirs_n = filter(d -> d !== opp_dir, _ALL_DIRECTIONS)

    T_a = psi.tensors[c]
    T_b = psi.tensors[c_other]
    p_a = physical_index(psi, c)

    # Absorb the diagonal weights on all non-shared bonds into the respective
    # endpoint tensors. The shared bond's weight is absorbed into T_a so the
    # SVD's singular values become the new bond weight (no further rescaling
    # needed).
    for d in other_dirs_c
        T_a = absorb_link_weight(T_a, psi, c, d)
    end
    for d in other_dirs_n
        T_b = absorb_link_weight(T_b, psi, c_other, d)
    end
    T_a = absorb_link_weight(T_a, psi, c, dir)

    theta = T_a * T_b

    left_links = [link_index(psi, c, d) for d in other_dirs_c]
    U, S, V, spec, u, v =
        svd(theta, p_a, left_links...; maxdim = target_maxdim, cutoff = cutoff)

    svals = _singular_values_from_S(S)
    lambda_new, scale = _normalize_singular_values(svals)
    keptdim = length(lambda_new)
    truncerr = Float64(spec.truncerr)

    new_link = Index(keptdim, "link,$(c.x),$(c.y),$(dir)")
    U_relabeled = replaceind(U, u, new_link)
    V_relabeled = replaceind(V, v, new_link)

    # Restore canonical (Γ) form: deabsorb the non-shared neighbor weights and
    # absorb the global scale on the U side so log-norm bookkeeping is exact.
    T_a_new = U_relabeled * scale
    for d in other_dirs_c
        T_a_new = deabsorb_link_weight(T_a_new, psi, c, d)
    end
    T_b_new = V_relabeled
    for d in other_dirs_n
        T_b_new = deabsorb_link_weight(T_b_new, psi, c_other, d)
    end

    # Replace the shared link Index in psi consistently.
    a_link_lookup = function (side_dir)
        if side_dir === dir
            return new_link
        else
            return link_index(psi, c, side_dir)
        end
    end
    b_link_lookup = function (side_dir)
        if side_dir === opp_dir
            return new_link
        else
            return link_index(psi, c_other, side_dir)
        end
    end
    psi.tensors[c] = permute(T_a_new, _site_order_indices(psi, c, a_link_lookup)...)
    psi.tensors[c_other] =
        permute(T_b_new, _site_order_indices(psi, c_other, b_link_lookup)...)
    psi.link_indices[(c, dir)] = new_link
    psi.link_indices[(c_other, opp_dir)] = new_link
    psi.link_weights[bondkey(psi.unitcell, c, dir)] = lambda_new

    return (truncerr = truncerr, keptdim = keptdim, scale = scale)
end

"""
    compress_to_target_maxdim!(psi, target_maxdim; cutoff = 1e-12)

Compress every canonical (`:right` and `:up`) bond of `psi` to at most
`target_maxdim` singular values via a per-bond SVD truncation in Γ-λ
simple-update form. Endpoint tensors and the shared link weight are updated
in place; the discarded weight is accumulated into per-bond truncation
errors. Returns an [`IPEPSCompressionInfo`](@ref).

`target_maxdim` must be a positive integer; `cutoff` is forwarded to the
ITensors SVD as the absolute singular-value floor. The state's
`log_norm_value` ledger is incremented by `sum(log, per-bond scales)` so the
compressed state remains normalization-tracked.

This primitive is intended for ScarFinder's evolve-then-compress loop:
evolve at a faithful `D_evol = trotter.maxdim`, then compress back to
`D_target = target_maxdim`. The accumulated truncation error is the natural
scar-quality signal — small means the trajectory genuinely lived in the
`D_target` manifold; large means we forced it.
"""
function compress_to_target_maxdim!(
    psi::SquareIPEPSState,
    target_maxdim::Integer;
    cutoff::Real = 1e-12,
)::IPEPSCompressionInfo
    target = Int(target_maxdim)
    target >= 1 || throw(ArgumentError("target_maxdim must be at least 1"))
    cutoff_value = Float64(cutoff)
    isfinite(cutoff_value) && cutoff_value >= 0 ||
        throw(ArgumentError("cutoff must be finite and nonnegative"))

    bond_truncerrs = Dict{BondKey,Float64}()
    bond_keptdims = Dict{BondKey,Int}()
    log_norm_delta = 0.0

    for c in psi.unitcell.reps
        for dir in _PRIMARY_BOND_DIRECTIONS
            result = _compress_bond!(psi, c, dir, target, cutoff_value)
            key = bondkey(psi.unitcell, c, dir)
            bond_truncerrs[key] = result.truncerr
            bond_keptdims[key] = result.keptdim
            log_norm_delta += log(result.scale)
        end
    end

    if log_norm_delta != 0.0
        _add_log_norm!(psi, log_norm_delta)
    end
    _mark_mutated!(psi)

    truncerrs = collect(values(bond_truncerrs))
    max_truncerr = isempty(truncerrs) ? 0.0 : maximum(truncerrs)
    mean_truncerr = isempty(truncerrs) ? 0.0 : sum(truncerrs) / length(truncerrs)

    return IPEPSCompressionInfo(
        target,
        max_truncerr,
        mean_truncerr,
        bond_truncerrs,
        bond_keptdims,
        log_norm_delta,
    )
end

end # module
