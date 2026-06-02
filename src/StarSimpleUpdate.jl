module StarSimpleUpdate

using ITensors
using LinearAlgebra

using ..SquareGeometry: SquareCoord
using ..SquarePXP: SQUARE_STAR_SITES
using ..SquareUnitCells: PeriodicSquareUnitCell, BondKey, wrap, neighbor, bondkey
using ..SquareIPEPS:
    SquareIPEPSState,
    physical_index,
    link_index,
    link_weight,
    log_norm,
    _link_weight_view,
    absorb_link_weight,
    deabsorb_link_weight,
    _mark_mutated!,
    _add_log_norm!
using ..StarModels: AbstractStarModel, PXPStarModel, star_gate_itensor

export StarUpdateInfo, project_star!, canonicalize_simple!

const _STAR_DIRECTIONS = (:right, :up, :left, :down)

# State-preserving regauge cutoff: small enough that dropping the discarded
# weight leaves the wavefunction unchanged to ~1e-10, large enough to discard
# the genuinely-null directions whose near-zero weights would otherwise have to
# be inverted during later de-absorption.
const _GAUGE_CUTOFF = 1e-20

"""
    StarUpdateInfo

Diagnostics returned by [`project_star!`](@ref) for one QR-reduced five-site
square-star simple update. `truncerrs`, `keptdims`, `min_lambda`, and
`norm_factors` are keyed by the center-to-leaf directions `:right`, `:up`,
`:left`, and `:down`. `touched_min_lambda` records the pre-update minimum
stored link weight for every canonical bond touched by the star patch.
"""
struct StarUpdateInfo
    center::SquareCoord
    max_truncerr::Float64
    truncerrs::Dict{Symbol,Float64}
    keptdims::Dict{Symbol,Int}
    min_lambda::Dict{Symbol,Float64}
    touched_min_lambda::Dict{BondKey,Float64}
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

function _validate_step(step::Real)
    isfinite(step) || throw(ArgumentError("step must be finite"))
    return step
end

function _validate_split_order(split_order)
    order_values = try
        collect(split_order)
    catch
        throw(
            ArgumentError(
                "split_order must be a permutation of (:right, :up, :left, :down)",
            ),
        )
    end
    all(dir -> dir isa Symbol, order_values) || throw(
        ArgumentError("split_order must be a permutation of (:right, :up, :left, :down)"),
    )
    order = Tuple(order_values)
    length(order) == length(_STAR_DIRECTIONS) ||
        throw(ArgumentError("split_order must contain four directions"))
    all(dir -> dir in _STAR_DIRECTIONS, order) || throw(
        ArgumentError("split_order must be a permutation of (:right, :up, :left, :down)"),
    )
    Set(order) == Set(_STAR_DIRECTIONS) || throw(
        ArgumentError("split_order must be a permutation of (:right, :up, :left, :down)"),
    )
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
    length(distinct) == SQUARE_STAR_SITES || throw(
        ArgumentError(
            "wrapped square star must contain five distinct unit-cell representatives",
        ),
    )
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

function _pre_update_touched_min_lambda(psi::SquareIPEPSState, coords)
    keys = Set{BondKey}()
    for dir in _STAR_DIRECTIONS
        push!(keys, bondkey(psi.unitcell, coords.center, dir))
        leaf = getproperty(coords, dir)
        for ext in _external_dirs_for_leaf(dir)
            push!(keys, bondkey(psi.unitcell, leaf, ext))
        end
    end
    values = Dict{BondKey,Float64}()
    for key in keys
        lambda = _link_weight_view(psi, key.site, key.dir)
        minimum_lambda = minimum(lambda)
        isfinite(minimum_lambda) && minimum_lambda >= 0 ||
            throw(ArgumentError("touched link weights must be finite and nonnegative"))
        values[key] = minimum_lambda
    end
    return values
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
    return (
        p,
        links_by_dir[:left],
        links_by_dir[:right],
        links_by_dir[:up],
        links_by_dir[:down],
    )
end

function _singular_values_from_S(S)::Vector{Float64}
    sinds = inds(S)
    length(sinds) == 2 || throw(ArgumentError("S tensor must have exactly two indices"))
    left, right = sinds
    n = min(dim(left), dim(right))
    values = Vector{Float64}(undef, n)
    for i = 1:n
        values[i] = abs(S[left=>i, right=>i])
    end
    return values
end

function _normalized_link_weight_from_singular_values(svals::Vector{Float64})
    !isempty(svals) ||
        throw(ArgumentError("empty singular spectrum encountered during star split"))
    all(isfinite, svals) ||
        throw(ArgumentError("singular values must be finite during star split"))
    all(>=(0), svals) ||
        throw(ArgumentError("singular values must be nonnegative during star split"))

    scale = norm(svals)
    isfinite(scale) && scale > 0 || throw(
        ArgumentError("zero or non-finite singular spectrum encountered during star split"),
    )

    lambda = svals ./ scale
    all(isfinite, lambda) ||
        throw(ArgumentError("normalized link weights must be finite during star split"))
    isapprox(norm(lambda), 1; atol = 1e-12, rtol = 1e-12) ||
        throw(ArgumentError("normalized link weights must have unit norm"))
    return lambda, scale
end

function _validate_split_diagnostic(truncerr::Real, min_lambda::Real, norm_factor::Real)
    isfinite(truncerr) && truncerr >= 0 ||
        throw(ArgumentError("truncation error must be finite and nonnegative"))
    isfinite(min_lambda) && min_lambda >= 0 ||
        throw(ArgumentError("minimum link weight must be finite and nonnegative"))
    isfinite(norm_factor) && norm_factor > 0 ||
        throw(ArgumentError("star split norm factor must be finite and positive"))
    return nothing
end

function _new_link_index(center::SquareCoord, dir::Symbol, keptdim::Int)
    return Index(keptdim, "link,$(center.x),$(center.y),$dir,star")
end

function _validate_rel_floor(rel_floor::Real)
    value = Float64(rel_floor)
    isfinite(value) || throw(ArgumentError("rel_floor must be finite"))
    0 <= value < 1 || throw(ArgumentError("rel_floor must be in [0, 1)"))
    return value
end

# SVD with a relative singular-value condition floor. ITensors' `cutoff` bounds
# the *discarded* relative weight, which still admits retained singular values
# as small as ~sqrt(cutoff) of the leading one; deabsorbing the full lambda then
# divides by those near-zero weights and amplifies noise on later star updates
# (the documented larger-D-worsens-at-tight-cutoff instability). `rel_floor`
# additionally drops any retained singular value below `rel_floor * sigma_max`,
# capping the bond condition number at 1/rel_floor regardless of `cutoff`.
function _svd_with_rel_floor(core, left, right, maxdim::Int, cutoff::Float64, rel_floor::Float64)
    U, S, V, spec, u, v = svd(core, left, right; maxdim = maxdim, cutoff = cutoff)
    if rel_floor > 0
        svals = _singular_values_from_S(S)
        if !isempty(svals) && isfinite(svals[1]) && svals[1] > 0
            keep = max(1, count(>=(rel_floor * svals[1]), svals))
            if keep < length(svals)
                U, S, V, spec, u, v =
                    svd(core, left, right; maxdim = keep, cutoff = cutoff)
            end
        end
    end
    return U, S, V, spec, u, v
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
            hasind(theta, link_index(psi, leaf, ext)) && error(
                "reduced star core unexpectedly contains external virtual leg $ext for $dir leaf",
            )
        end
        hasind(theta, qinds[dir]) ||
            error("reduced star core is missing QR index for $dir leaf")
    end
    for p in phys
        hasind(theta, p) || error("reduced star core is missing a physical index")
    end
    return nothing
end

function _split_reduced_theta(
    theta::ITensor,
    psi::SquareIPEPSState,
    coords,
    qinds,
    split_order,
    maxdim,
    cutoff,
    rel_floor,
    touched_min_lambda,
)
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
        U, S, V, spec, u, v =
            _svd_with_rel_floor(core, p_leaf, q_leaf, maxdim, cutoff, rel_floor)
        svals = _singular_values_from_S(S)
        lambda_new, scale = _normalized_link_weight_from_singular_values(svals)
        keptdim = length(lambda_new)
        new_link = _new_link_index(coords.center, dir, keptdim)
        truncerr = Float64(spec.truncerr)
        lambda_min = minimum(lambda_new)
        _validate_split_diagnostic(truncerr, lambda_min, scale)

        leaf_active[dir] = replaceind(U, u, new_link)
        core = replaceind(V, v, new_link) * scale
        new_links[dir] = new_link
        new_weights[dir] = lambda_new
        truncerrs[dir] = truncerr
        keptdims[dir] = keptdim
        min_lambda[dir] = lambda_min
        norm_factors[dir] = scale
    end

    max_truncerr = maximum(values(truncerrs))
    info = StarUpdateInfo(
        coords.center,
        max_truncerr,
        truncerrs,
        keptdims,
        min_lambda,
        touched_min_lambda,
        norm_factors,
    )
    return core, leaf_active, new_links, new_weights, info
end

function _reconstruct_leaf(
    psi::SquareIPEPSState,
    coords,
    dir::Symbol,
    Q,
    leaf_active,
    new_link,
)
    leaf = getproperty(coords, dir)
    T = Q * leaf_active
    for ext in _external_dirs_for_leaf(dir)
        T = deabsorb_link_weight(T, psi, leaf, ext)
    end

    links_by_dir = Dict{Symbol,Index}()
    for site_dir in _STAR_DIRECTIONS
        links_by_dir[site_dir] =
            site_dir === _center_dir_for_leaf(dir) ? new_link :
            link_index(psi, leaf, site_dir)
    end
    return permute(T, _site_order_indices(psi, leaf, links_by_dir)...)
end

function _reconstruct_center(psi::SquareIPEPSState, coords, core, new_links)
    links_by_dir = Dict(dir => new_links[dir] for dir in _STAR_DIRECTIONS)
    return permute(core, _site_order_indices(psi, coords.center, links_by_dir)...)
end

function _commit_star_update!(
    psi::SquareIPEPSState,
    coords,
    new_tensors,
    new_links,
    new_weights,
    info::StarUpdateInfo,
)
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
    _add_log_norm!(psi, sum(log, values(info.norm_factors)))
    _mark_mutated!(psi)
    return psi
end

"""
    project_star!(
        psi::SquareIPEPSState,
        center::SquareCoord,
        step::Real;
        model::Union{AbstractStarModel,Nothing} = nothing,
        evolution::Symbol = :real,
        projected::Bool = true,
        maxdim::Integer = psi.maxdim,
        cutoff::Real = 1e-12,
        rel_floor::Real = 0.0,
        split_order = (:right, :up, :left, :down),
    )::StarUpdateInfo

Apply one QR-reduced five-site square-star simple update to the custom
ITensors iPEPS backend. The star order and gate physical order are
`(center, right, up, left, down)`, and only this one local star is updated.
The update is transactional: `psi` is mutated only after local weight
absorption, QR reduction, gate application, SVD splitting, and reconstruction
all succeed. `maxdim` is the cap for this update and may be larger than
`psi.maxdim`; in that case affected links may grow while `psi.maxdim` remains
the state's construction/default cap. `rel_floor` (in `[0, 1)`) drops any
retained singular value below `rel_floor * sigma_max` on each bond split, capping
the bond condition number at `1/rel_floor` so the de-absorb step never inverts a
near-zero link weight. `rel_floor = 0` disables the floor (legacy behavior).

# Example

```julia
using SquarePXPDynamics
cell = PeriodicSquareUnitCell(10, 10)
psi = product_square_ipeps(cell; state = :down, maxdim = 1)
info = project_star!(psi, SquareCoord(1, 1), 0.01)
# info.max_truncerr summarizes the largest singular-value truncation
# error across the four bond splits of this star.
```
"""
function project_star!(
    psi::SquareIPEPSState,
    center::SquareCoord,
    step::Real;
    model::Union{AbstractStarModel,Nothing} = nothing,
    evolution::Symbol = :real,
    projected::Bool = true,
    maxdim::Integer = psi.maxdim,
    cutoff::Real = 1e-12,
    rel_floor::Real = 0.0,
    split_order = _STAR_DIRECTIONS,
)::StarUpdateInfo
    finite_step = _validate_step(step)
    kept_maxdim = _validate_maxdim(maxdim)
    trunc_cutoff = _validate_cutoff(cutoff)
    trunc_rel_floor = _validate_rel_floor(rel_floor)
    order = _validate_split_order(split_order)
    coords = _validate_distinct_star!(psi, center)
    touched_min_lambda = _pre_update_touched_min_lambda(psi, coords)
    phys = _star_physical_indices(psi, coords)

    actual_model = model === nothing ? PXPStarModel(projected) : model
    gate = star_gate_itensor(actual_model, phys, finite_step; evolution = evolution)

    center_tensor, absorbed_leaves, _ = _absorb_star_weights(psi, coords)
    qfactors, rfactors, qinds = _qr_reduce_leaves(psi, coords, absorbed_leaves)

    theta =
        gate *
        center_tensor *
        rfactors[:right] *
        rfactors[:up] *
        rfactors[:left] *
        rfactors[:down]
    theta = noprime(theta)
    _assert_reduced_theta(theta, psi, coords, phys, qinds)

    # A fully blockade-forbidden star (excited center with an excited neighbor)
    # is annihilated by the projected PXP gate, leaving theta == 0. Detect it here
    # and report the offending center clearly, instead of throwing an opaque
    # "zero singular spectrum" deep inside the SVD split. This is the failure mode
    # for :z_up (all-up) and for checkerboard seams on odd-period cells.
    theta_norm = norm(theta)
    isfinite(theta_norm) && theta_norm > sqrt(eps(Float64)) || throw(
        ArgumentError(
            "projected square-star update annihilated the state at center " *
            "$(coords.center) (||theta|| = $theta_norm): this star is " *
            "blockade-forbidden under the projected PXP gate (an excited center " *
            "adjacent to an excited site). Initial states such as :z_up (all-up) " *
            "and checkerboard seams on odd-period cells trigger this; use :down " *
            "or an even-tileable checkerboard cell with schedule = :serial.",
        ),
    )

    center_core, leaf_active, new_links, new_weights, info =
        _split_reduced_theta(
            theta,
            psi,
            coords,
            qinds,
            order,
            kept_maxdim,
            trunc_cutoff,
            trunc_rel_floor,
            touched_min_lambda,
        )

    new_tensors = Dict{SquareCoord,ITensor}()
    new_tensors[coords.center] = _reconstruct_center(psi, coords, center_core, new_links)
    for dir in _STAR_DIRECTIONS
        new_tensors[getproperty(coords, dir)] = _reconstruct_leaf(
            psi,
            coords,
            dir,
            qfactors[dir],
            leaf_active[dir],
            new_links[dir],
        )
    end

    _commit_star_update!(psi, coords, new_tensors, new_links, new_weights, info)
    return info
end

function _max_link_weight_change(before, after)
    change = 0.0
    for (key, lambda_after) in after
        lambda_before = get(before, key, nothing)
        (lambda_before === nothing || length(lambda_before) != length(lambda_after)) &&
            return Inf
        change = max(change, maximum(abs.(lambda_after .- lambda_before)))
    end
    return change
end

"""
    canonicalize_simple!(psi::SquareIPEPSState; max_sweeps = 8, tol = 1e-10)::Int

Regauge `psi` toward the mean-field super-orthogonal (Vidal canonical) form so
each stored link weight approximates the Schmidt spectrum of its bond. This is
done with idle square-star updates — [`project_star!`](@ref) at `step = 0` with
the unprojected identity gate — which are exact local gauge transformations:
the represented wavefunction is left unchanged while the link weights relax
toward canonical values. Serial sweeps over all unit-cell centers repeat until
the largest change in any stored weight is below `tol` or `max_sweeps` sweeps
have run; the return value is the number of sweeps performed. The log-norm
ledger is restored afterward because a gauge transformation does not change the
physical norm.

# Example

```julia
using SquarePXPDynamics
cell = PeriodicSquareUnitCell(4, 4)
psi = checkerboard_square_ipeps(cell; maxdim = 2)
evolve!(psi, 0.1; dt = 0.05, order = 1, maxdim = 2, schedule = :serial)
canonicalize_simple!(psi)  # link weights now approximate Schmidt spectra
```
"""
function canonicalize_simple!(
    psi::SquareIPEPSState;
    max_sweeps::Integer = 8,
    tol::Real = 1e-10,
)::Int
    max_sweeps >= 1 || throw(ArgumentError("max_sweeps must be at least 1"))
    tol > 0 || throw(ArgumentError("tol must be positive"))
    saved_log_norm = log_norm(psi)
    gauge_maxdim = maximum(dim(link) for link in values(psi.link_indices))
    sweeps = 0
    for _ = 1:max_sweeps
        before = Dict(key => copy(lambda) for (key, lambda) in psi.link_weights)
        for center in psi.unitcell.reps
            project_star!(
                psi,
                center,
                0.0;
                projected = false,
                maxdim = gauge_maxdim,
                cutoff = _GAUGE_CUTOFF,
                rel_floor = 0.0,
            )
        end
        sweeps += 1
        _max_link_weight_change(before, psi.link_weights) < tol && break
    end
    _add_log_norm!(psi, saved_log_norm - log_norm(psi))
    return sweeps
end

end
