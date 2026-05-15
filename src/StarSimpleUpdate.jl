module StarSimpleUpdate

using ITensors
using LinearAlgebra

using ..SquareGeometry: SquareCoord
using ..SquarePXP: SQUARE_STAR_SITES
using ..SquareUnitCells: PeriodicSquareUnitCell, wrap, neighbor, bondkey
using ..SquareIPEPS:
    SquareIPEPSState,
    physical_index,
    link_index,
    link_weight,
    absorb_link_weight,
    deabsorb_link_weight,
    square_pxp_gate_itensor,
    projected_square_pxp_gate_itensor

export StarUpdateInfo, project_star!

const _STAR_DIRECTIONS = (:right, :up, :left, :down)

"""
    StarUpdateInfo

Diagnostics returned by [`project_star!`](@ref) for one QR-reduced five-site
square-star simple update. `truncerrs`, `keptdims`, `min_lambda`, and
`norm_factors` are keyed by the center-to-leaf directions `:right`, `:up`,
`:left`, and `:down`.
"""
struct StarUpdateInfo
    center::SquareCoord
    max_truncerr::Float64
    truncerrs::Dict{Symbol,Float64}
    keptdims::Dict{Symbol,Int}
    min_lambda::Dict{Symbol,Float64}
    norm_factors::Dict{Symbol,Float64}
end

function _validate_maxdim(maxdim::Integer)
    maxdim >= 1 || throw(ArgumentError("maxdim must be at least 1"))
    return Int(maxdim)
end

function _validate_cutoff(cutoff::Real)
    value = Float64(cutoff)
    isfinite(value) || throw(ArgumentError("cutoff must be finite"))
    value >= 0 || throw(ArgumentError("cutoff must be nonnegative"))
    return value
end

function _validate_split_order(split_order)
    order_values = try
        collect(split_order)
    catch
        throw(ArgumentError("split_order must be a permutation of (:right, :up, :left, :down)"))
    end
    all(dir -> dir isa Symbol, order_values) ||
        throw(ArgumentError("split_order must be a permutation of (:right, :up, :left, :down)"))
    order = Tuple(order_values)
    length(order) == length(_STAR_DIRECTIONS) ||
        throw(ArgumentError("split_order must contain four directions"))
    all(dir -> dir in _STAR_DIRECTIONS, order) ||
        throw(ArgumentError("split_order must be a permutation of (:right, :up, :left, :down)"))
    Set(order) == Set(_STAR_DIRECTIONS) ||
        throw(ArgumentError("split_order must be a permutation of (:right, :up, :left, :down)"))
    return order
end

function _star_coords(psi::SquareIPEPSState, center::SquareCoord)
    cell = psi.unitcell
    c = wrap(cell, center)
    return (
        center = c,
        right = neighbor(cell, c, :right),
        up = neighbor(cell, c, :up),
        left = neighbor(cell, c, :left),
        down = neighbor(cell, c, :down),
    )
end

function _validate_distinct_star!(psi::SquareIPEPSState, center::SquareCoord)
    coords = _star_coords(psi, center)
    distinct = Set((coords.center, coords.right, coords.up, coords.left, coords.down))
    length(distinct) == SQUARE_STAR_SITES ||
        throw(ArgumentError("wrapped square star must contain five distinct unit-cell representatives"))
    return coords
end

function _opposite_dir(dir::Symbol)
    dir === :right && return :left
    dir === :up && return :down
    dir === :left && return :right
    dir === :down && return :up
    throw(ArgumentError("direction must be :right, :up, :left, or :down"))
end

function _center_dir_for_leaf(dir::Symbol)
    dir === :right && return :left
    dir === :up && return :down
    dir === :left && return :right
    dir === :down && return :up
    throw(ArgumentError("leaf direction must be :right, :up, :left, or :down"))
end

function _external_dirs_for_leaf(dir::Symbol)
    dir === :right && return (:right, :up, :down)
    dir === :up && return (:left, :right, :up)
    dir === :left && return (:left, :up, :down)
    dir === :down && return (:left, :right, :down)
    throw(ArgumentError("leaf direction must be :right, :up, :left, or :down"))
end

function _star_physical_indices(psi::SquareIPEPSState, coords)
    return (
        physical_index(psi, coords.center),
        physical_index(psi, coords.right),
        physical_index(psi, coords.up),
        physical_index(psi, coords.left),
        physical_index(psi, coords.down),
    )
end

function _site_order_indices(psi::SquareIPEPSState, c::SquareCoord, links_by_dir)
    p = physical_index(psi, c)
    return (p, links_by_dir[:left], links_by_dir[:right], links_by_dir[:up], links_by_dir[:down])
end

function _singular_values_from_S(S)::Vector{Float64}
    sinds = inds(S)
    length(sinds) == 2 || throw(ArgumentError("S tensor must have exactly two indices"))
    left, right = sinds
    n = min(dim(left), dim(right))
    values = Vector{Float64}(undef, n)
    for i in 1:n
        values[i] = abs(S[left => i, right => i])
    end
    return values
end

function _new_link_index(center::SquareCoord, dir::Symbol, keptdim::Int)
    return Index(keptdim, "link,$(center.x),$(center.y),$dir,star")
end

function _absorb_star_weights(psi::SquareIPEPSState, coords)
    leaves = Dict{Symbol,ITensor}()
    external_absorbed = Dict{Symbol,NTuple{3,Symbol}}()
    for dir in _STAR_DIRECTIONS
        leaf = getproperty(coords, dir)
        T = copy(psi.tensors[leaf])
        to_center = _center_dir_for_leaf(dir)
        T = absorb_link_weight(T, psi, leaf, to_center)
        external_dirs = _external_dirs_for_leaf(dir)
        for ext in external_dirs
            T = absorb_link_weight(T, psi, leaf, ext)
        end
        leaves[dir] = T
        external_absorbed[dir] = external_dirs
    end
    return copy(psi.tensors[coords.center]), leaves, external_absorbed
end

function _qr_reduce_leaves(psi::SquareIPEPSState, coords, absorbed_leaves)
    qfactors = Dict{Symbol,ITensor}()
    rfactors = Dict{Symbol,ITensor}()
    qinds = Dict{Symbol,Index}()
    for dir in _STAR_DIRECTIONS
        leaf = getproperty(coords, dir)
        external = _external_dirs_for_leaf(dir)
        external_inds = map(ext -> link_index(psi, leaf, ext), external)
        Q, R, _, q = factorize(
            absorbed_leaves[dir],
            external_inds...;
            ortho = "left",
            which_decomp = "qr",
            tags = "Link,qr,$dir",
        )
        qfactors[dir] = Q
        rfactors[dir] = R
        qinds[dir] = q
    end
    return qfactors, rfactors, qinds
end

function _assert_reduced_theta(theta::ITensor, psi::SquareIPEPSState, coords, phys, qinds)
    for dir in _STAR_DIRECTIONS
        leaf = getproperty(coords, dir)
        for ext in _external_dirs_for_leaf(dir)
            hasind(theta, link_index(psi, leaf, ext)) &&
                error("reduced star core unexpectedly contains external virtual leg $ext for $dir leaf")
        end
        hasind(theta, qinds[dir]) ||
            error("reduced star core is missing QR index for $dir leaf")
    end
    for p in phys
        hasind(theta, p) || error("reduced star core is missing a physical index")
    end
    return nothing
end

function _split_reduced_theta(theta::ITensor, psi::SquareIPEPSState, coords, qinds, split_order, maxdim, cutoff)
    leaf_active = Dict{Symbol,ITensor}()
    new_links = Dict{Symbol,Index}()
    new_weights = Dict{Symbol,Vector{Float64}}()
    truncerrs = Dict{Symbol,Float64}()
    keptdims = Dict{Symbol,Int}()
    min_lambda = Dict{Symbol,Float64}()
    norm_factors = Dict{Symbol,Float64}()

    core = theta
    for dir in split_order
        leaf = getproperty(coords, dir)
        p_leaf = physical_index(psi, leaf)
        q_leaf = qinds[dir]
        U, S, V, spec, u, v = svd(core, p_leaf, q_leaf; maxdim = maxdim, cutoff = cutoff)
        svals = _singular_values_from_S(S)
        scale = norm(svals)
        scale > 0 || throw(ArgumentError("zero singular spectrum encountered during star split"))

        lambda_new = svals ./ scale
        keptdim = length(lambda_new)
        new_link = _new_link_index(coords.center, dir, keptdim)

        leaf_active[dir] = replaceind(U, u, new_link)
        core = replaceind(V, v, new_link) * scale
        new_links[dir] = new_link
        new_weights[dir] = lambda_new
        truncerrs[dir] = Float64(spec.truncerr)
        keptdims[dir] = keptdim
        min_lambda[dir] = minimum(lambda_new)
        norm_factors[dir] = scale
    end

    max_truncerr = maximum(values(truncerrs))
    info = StarUpdateInfo(
        coords.center,
        max_truncerr,
        truncerrs,
        keptdims,
        min_lambda,
        norm_factors,
    )
    return core, leaf_active, new_links, new_weights, info
end

function _reconstruct_leaf(psi::SquareIPEPSState, coords, dir::Symbol, Q, leaf_active, new_link)
    leaf = getproperty(coords, dir)
    T = Q * leaf_active
    for ext in _external_dirs_for_leaf(dir)
        T = deabsorb_link_weight(T, psi, leaf, ext)
    end

    links_by_dir = Dict{Symbol,Index}()
    for site_dir in _STAR_DIRECTIONS
        links_by_dir[site_dir] =
            site_dir === _center_dir_for_leaf(dir) ? new_link : link_index(psi, leaf, site_dir)
    end
    return permute(T, _site_order_indices(psi, leaf, links_by_dir)...)
end

function _reconstruct_center(psi::SquareIPEPSState, coords, core, new_links)
    links_by_dir = Dict(dir => new_links[dir] for dir in _STAR_DIRECTIONS)
    return permute(core, _site_order_indices(psi, coords.center, links_by_dir)...)
end

function _commit_star_update!(psi::SquareIPEPSState, coords, new_tensors, new_links, new_weights)
    for (site, tensor) in new_tensors
        psi.tensors[site] = tensor
    end
    for dir in _STAR_DIRECTIONS
        leaf = getproperty(coords, dir)
        to_center = _center_dir_for_leaf(dir)
        link = new_links[dir]
        psi.link_indices[(coords.center, dir)] = link
        psi.link_indices[(leaf, to_center)] = link
        psi.link_weights[bondkey(psi.unitcell, coords.center, dir)] = new_weights[dir]
    end
    return psi
end

"""
    project_star!(
        psi::SquareIPEPSState,
        center::SquareCoord,
        step::Real;
        evolution::Symbol = :real,
        projected::Bool = true,
        maxdim::Integer = psi.maxdim,
        cutoff::Real = 1e-12,
        split_order = (:right, :up, :left, :down),
    )::StarUpdateInfo

Apply one QR-reduced five-site square-star simple update to the custom
ITensors iPEPS backend. The star order and gate physical order are
`(center, right, up, left, down)`, and only this one local star is updated.
The update is transactional: `psi` is mutated only after local weight
absorption, QR reduction, gate application, SVD splitting, and reconstruction
all succeed.
"""
function project_star!(
    psi::SquareIPEPSState,
    center::SquareCoord,
    step::Real;
    evolution::Symbol = :real,
    projected::Bool = true,
    maxdim::Integer = psi.maxdim,
    cutoff::Real = 1e-12,
    split_order = _STAR_DIRECTIONS,
)::StarUpdateInfo
    kept_maxdim = _validate_maxdim(maxdim)
    trunc_cutoff = _validate_cutoff(cutoff)
    order = _validate_split_order(split_order)
    coords = _validate_distinct_star!(psi, center)
    phys = _star_physical_indices(psi, coords)

    gate = projected ?
        projected_square_pxp_gate_itensor(phys, step; evolution = evolution) :
        square_pxp_gate_itensor(phys, step; evolution = evolution)

    center_tensor, absorbed_leaves, _ = _absorb_star_weights(psi, coords)
    qfactors, rfactors, qinds = _qr_reduce_leaves(psi, coords, absorbed_leaves)

    theta = gate * center_tensor * rfactors[:right] * rfactors[:up] * rfactors[:left] * rfactors[:down]
    theta = noprime(theta)
    _assert_reduced_theta(theta, psi, coords, phys, qinds)

    center_core, leaf_active, new_links, new_weights, info =
        _split_reduced_theta(theta, psi, coords, qinds, order, kept_maxdim, trunc_cutoff)

    new_tensors = Dict{SquareCoord,ITensor}()
    new_tensors[coords.center] = _reconstruct_center(psi, coords, center_core, new_links)
    for dir in _STAR_DIRECTIONS
        new_tensors[getproperty(coords, dir)] =
            _reconstruct_leaf(psi, coords, dir, qfactors[dir], leaf_active[dir], new_links[dir])
    end

    _commit_star_update!(psi, coords, new_tensors, new_links, new_weights)
    return info
end

end
