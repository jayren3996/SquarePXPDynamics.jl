module SquareIPEPS

using ITensors
using ..SquareGeometry
using ..SquarePXP
using ..SquareUnitCells
import ..SquarePEPS: physical_index, link_index

export SquareIPEPSState
export product_square_ipeps, checkerboard_square_ipeps
export unitcell_reps, physical_dim, simple_weight_dim, copy_state
export physical_index, link_index
export link_weight, set_link_weight!, link_weight_tensor
export state_version, log_norm
export absorb_link_weight, deabsorb_link_weight
export weight_entropy, bond_entropy, all_bond_entropies
export square_pxp_gate_itensor, projected_square_pxp_gate_itensor

"""
    SquareIPEPSState

Periodic square-lattice iPEPS state in a simple-update Γ-λ representation.
The site tensors are bare Γ tensors, `link_weights` store explicit λ spectra
for canonical bonds, and `gauge` is currently `:simple`. `maxdim` records the
construction/default bond-dimension cap; individual updates may request a
larger `maxdim` and grow affected links without rewriting this default.
"""
struct SquareIPEPSState
    unitcell::PeriodicSquareUnitCell
    tensors::Dict{SquareCoord,ITensor}
    physical_indices::Dict{SquareCoord,Index}
    link_indices::Dict{Tuple{SquareCoord,Symbol},Index}
    link_weights::Dict{BondKey,Vector{Float64}}
    maxdim::Int
    gauge::Symbol
    mutation_version::Base.RefValue{Int}
    log_norm_value::Base.RefValue{Float64}
end

function _validate_ipeps_cell(cell::PeriodicSquareUnitCell)
    cell.Lx >= 2 || throw(ArgumentError("iPEPS construction requires Lx >= 2"))
    cell.Ly >= 2 || throw(ArgumentError("iPEPS construction requires Ly >= 2"))
    return cell
end

function _validate_maxdim(maxdim::Integer)
    maxdim >= 1 || throw(ArgumentError("maxdim must be at least 1"))
    return Int(maxdim)
end

function _product_state_vector(state::Symbol)
    if state in (:up, :z_up)
        return ComplexF64[1, 0]
    elseif state in (:down, :z_down)
        return ComplexF64[0, 1]
    elseif state === :x_plus
        return ComplexF64[inv(sqrt(2)), inv(sqrt(2))]
    else
        throw(ArgumentError("state must be :up, :down, :z_up, :z_down, or :x_plus"))
    end
end

function _lambda_vector(maxdim::Int)
    lambda = zeros(Float64, maxdim)
    lambda[1] = 1.0
    return lambda
end

function _validate_link_direction(dir::Symbol)
    dir in (:right, :up, :left, :down) ||
        throw(ArgumentError("direction must be :right, :up, :left, or :down"))
    return dir
end

function _build_periodic_links(cell::PeriodicSquareUnitCell, maxdim::Int)
    links = Dict{Tuple{SquareCoord,Symbol},Index}()
    weights = Dict{BondKey,Vector{Float64}}()

    for c in cell.reps
        right_key = BondKey(c, :right)
        right_index = Index(maxdim, "link,$(c.x),$(c.y),right")
        links[(c, :right)] = right_index
        links[(neighbor(cell, c, :right), :left)] = right_index
        weights[right_key] = _lambda_vector(maxdim)

        up_key = BondKey(c, :up)
        up_index = Index(maxdim, "link,$(c.x),$(c.y),up")
        links[(c, :up)] = up_index
        links[(neighbor(cell, c, :up), :down)] = up_index
        weights[up_key] = _lambda_vector(maxdim)
    end

    return links, weights
end

function _product_square_ipeps(cell::PeriodicSquareUnitCell, state_at; maxdim::Integer)
    _validate_ipeps_cell(cell)
    dim = _validate_maxdim(maxdim)
    links, weights = _build_periodic_links(cell, dim)
    physical = Dict(c => Index(2, "phys,$(c.x),$(c.y)") for c in cell.reps)
    tensors = Dict{SquareCoord,ITensor}()

    for c in cell.reps
        p = physical[c]
        left = links[(c, :left)]
        right = links[(c, :right)]
        up = links[(c, :up)]
        down = links[(c, :down)]
        tensor = ITensor(ComplexF64, p, left, right, up, down)
        state = state_at(c)
        if state isa Integer
            tensor[p=>state, left=>1, right=>1, up=>1, down=>1] = 1.0 + 0.0im
        else
            length(state) == 2 ||
                throw(ArgumentError("product state vector must have length 2"))
            for s = 1:2
                tensor[p=>s, left=>1, right=>1, up=>1, down=>1] = state[s]
            end
        end
        tensors[c] = tensor
    end

    return SquareIPEPSState(cell, tensors, physical, links, weights, dim, :simple, Ref(0), Ref(0.0))
end

"""
    state_version(psi)

Return the mutation counter for `psi`. Operations that change tensors or
stored link weights increment this counter so cached measurement contexts can
detect stale states.
"""
state_version(psi::SquareIPEPSState)::Int = psi.mutation_version[]

"""
    log_norm(psi)

Return the accumulated logarithmic normalization ledger recorded during
simple-update star splits. This is diagnostic bookkeeping for normalized link
spectra; it is not a physical observable.
"""
log_norm(psi::SquareIPEPSState)::Float64 = psi.log_norm_value[]

"""
    unitcell_reps(psi)

Return a copy of the periodic unit-cell representatives used by `psi`.
"""
unitcell_reps(psi::SquareIPEPSState)::Vector{SquareCoord} = copy(psi.unitcell.reps)

"""
    physical_dim(psi, c)

Return the physical Hilbert-space dimension at coordinate `c`, after periodic
wrapping into the unit cell.
"""
physical_dim(psi::SquareIPEPSState, c::SquareCoord)::Int = dim(physical_index(psi, c))

"""
    simple_weight_dim(psi, c, dir)

Return the length of the simple-update link-weight vector on the nearest
neighbor bond from `c` in direction `dir`.
"""
function simple_weight_dim(psi::SquareIPEPSState, c::SquareCoord, dir::Symbol)::Int
    _validate_link_direction(dir)
    return length(_validated_link_weight(psi, c, dir))
end

"""
    copy_state(psi)

Return a deep mutable copy of `psi` with independent tensors, link-weight
vectors, mutation counter, and log-normalization ledger.
"""
function copy_state(psi::SquareIPEPSState)::SquareIPEPSState
    return SquareIPEPSState(
        psi.unitcell,
        Dict(c => copy(T) for (c, T) in psi.tensors),
        copy(psi.physical_indices),
        copy(psi.link_indices),
        Dict(key => copy(lambda) for (key, lambda) in psi.link_weights),
        psi.maxdim,
        psi.gauge,
        Ref(state_version(psi)),
        Ref(log_norm(psi)),
    )
end

function _mark_mutated!(psi::SquareIPEPSState)
    psi.mutation_version[] += 1
    return psi
end

function _add_log_norm!(psi::SquareIPEPSState, delta::Real)
    value = Float64(delta)
    isfinite(value) || throw(ArgumentError("log-norm increment must be finite"))
    psi.log_norm_value[] += value
    return psi
end

"""
    product_square_ipeps(cell; state = :down, maxdim = 1)

Construct a periodic square iPEPS product state on `cell` using physical basis
`:up`/`:z_up` as index `1`, `:down`/`:z_down` as index `2`, or `:x_plus`
as equal amplitudes on both basis states. Neighboring representatives share
virtual `Index` objects, so `cell.Lx` and `cell.Ly` must both be at least `2`.
"""
function product_square_ipeps(
    cell::PeriodicSquareUnitCell;
    state::Symbol = :down,
    maxdim::Integer = 1,
)
    amplitudes = _product_state_vector(state)
    return _product_square_ipeps(cell, Returns(amplitudes); maxdim)
end

"""
    checkerboard_square_ipeps(cell; excited_on = :even, maxdim = 1)

Construct a periodic square iPEPS checkerboard product state. Even parity is
defined by `iseven(c.x + c.y)`; excited sites use physical basis index `1`.
"""
function checkerboard_square_ipeps(
    cell::PeriodicSquareUnitCell;
    excited_on::Symbol = :even,
    maxdim::Integer = 1,
)
    excited_on === :even ||
        excited_on === :odd ||
        throw(ArgumentError("excited_on must be :even or :odd"))

    function state_at(c::SquareCoord)
        even = iseven(c.x + c.y)
        excited = excited_on === :even ? even : !even
        return excited ? 1 : 2
    end

    return _product_square_ipeps(cell, state_at; maxdim)
end

"""
    physical_index(psi::SquareIPEPSState, c)

Return the physical ITensors index for coordinate `c`, after wrapping `c` into
the iPEPS unit cell.
"""
physical_index(psi::SquareIPEPSState, c::SquareCoord) =
    psi.physical_indices[wrap(psi.unitcell, c)]

"""
    link_index(psi::SquareIPEPSState, c, dir)

Return the virtual link index at coordinate `c` in `dir`, after wrapping `c`
into the iPEPS unit cell. All four directions are supported.
"""
function link_index(psi::SquareIPEPSState, c::SquareCoord, dir::Symbol)
    _validate_link_direction(dir)
    return psi.link_indices[(wrap(psi.unitcell, c), dir)]
end

"""
    link_weight(psi, c, dir)

Return a copy of the simple-update link-weight vector on the periodic nearest
neighbor bond from coordinate `c` in direction `dir`.
"""
function link_weight(psi::SquareIPEPSState, c::SquareCoord, dir::Symbol)
    _validate_link_direction(dir)
    return copy(_validated_link_weight(psi, c, dir))
end

function _validate_link_weight_values(link::Index, values)
    length(values) == dim(link) || throw(
        ArgumentError(
            "link weight length must match the corresponding link index dimension",
        ),
    )
    all(isfinite, values) || throw(ArgumentError("link weights must be finite"))
    all(x -> x >= 0, values) || throw(ArgumentError("link weights must be nonnegative"))
    any(!iszero, values) || throw(ArgumentError("link weights must not all be zero"))
    return Float64.(values)
end

function _validated_link_weight(psi::SquareIPEPSState, c::SquareCoord, dir::Symbol)
    link = link_index(psi, c, dir)
    return _validate_link_weight_values(
        link,
        psi.link_weights[bondkey(psi.unitcell, c, dir)],
    )
end

"""
    set_link_weight!(psi, c, dir, values)

Replace the simple-update link-weight vector on the periodic nearest-neighbor
bond from `c` in `dir`. Values must match the link dimension, be finite,
nonnegative, and not all zero. The stored vector is copied.
"""
function set_link_weight!(
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
    values::AbstractVector{<:Real},
)
    _validate_link_direction(dir)
    psi.link_weights[bondkey(psi.unitcell, c, dir)] =
        _validate_link_weight_values(link_index(psi, c, dir), values)
    _mark_mutated!(psi)
    return psi
end

"""
    link_weight_tensor(psi, c, dir)

Return the diagonal ITensor representation of `link_weight(psi, c, dir)` on the
link index for `(c, dir)` and its primed copy.
"""
function link_weight_tensor(psi::SquareIPEPSState, c::SquareCoord, dir::Symbol)
    link = link_index(psi, c, dir)
    lambda = link_weight(psi, c, dir)
    length(lambda) == dim(link) || throw(
        ArgumentError(
            "link weight length must match the corresponding link index dimension",
        ),
    )
    tensor = ITensor(ComplexF64, link, prime(link))
    for i in eachindex(lambda)
        tensor[link=>i, prime(link)=>i] = lambda[i]
    end
    return tensor
end

function _link_weight_tensor_from_values(link::Index, values::AbstractVector{<:Real})
    length(values) == dim(link) || throw(
        ArgumentError(
            "link weight length must match the corresponding link index dimension",
        ),
    )
    tensor = ITensor(ComplexF64, link, prime(link))
    for i in eachindex(values)
        tensor[link=>i, prime(link)=>i] = values[i]
    end
    return tensor
end

function _apply_link_weight_values(
    T::ITensor,
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
    values,
)
    link = link_index(psi, c, dir)
    hasind(T, link) || throw(
        ArgumentError(
            "tensor does not contain the virtual link index for direction $dir at coordinate $c",
        ),
    )
    scaled = T * _link_weight_tensor_from_values(link, values)
    return permute(replaceind(scaled, prime(link), link), inds(T))
end

"""
    absorb_link_weight(T::ITensor, psi::SquareIPEPSState, c::SquareCoord, dir::Symbol)::ITensor

Multiply tensor `T` by the diagonal simple-update link weight on the virtual
leg `(c, dir)`. The returned tensor preserves the original indices and index
ordering of `T`.
"""
function absorb_link_weight(
    T::ITensor,
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol,
)::ITensor
    return _apply_link_weight_values(T, psi, c, dir, link_weight(psi, c, dir))
end

"""
    deabsorb_link_weight(T::ITensor, psi::SquareIPEPSState, c::SquareCoord, dir::Symbol; atol = 1e-14)::ITensor

Divide tensor `T` by the diagonal simple-update link weight on the virtual leg
`(c, dir)`. Entries with `abs(lambda[i]) <= atol` use inverse value `0.0`, so
the result never contains `Inf` from zero or near-zero link weights.
"""
function deabsorb_link_weight(
    T::ITensor,
    psi::SquareIPEPSState,
    c::SquareCoord,
    dir::Symbol;
    atol = 1e-14,
)::ITensor
    isfinite(atol) && atol >= 0 ||
        throw(ArgumentError("atol must be finite and nonnegative"))
    inverse =
        map(lambda -> abs(lambda) <= atol ? 0.0 : inv(lambda), link_weight(psi, c, dir))
    return _apply_link_weight_values(T, psi, c, dir, inverse)
end

"""
    weight_entropy(lambda)

Return the normalized von Neumann entropy `-sum(p_i * log(p_i))` of a
nonnegative link-weight vector, with `p_i = lambda_i^2 / sum(lambda_j^2)`.
Zero probabilities do not contribute.
"""
function weight_entropy(lambda)
    values = Float64.(lambda)
    !isempty(values) || throw(ArgumentError("link weight vector must not be empty"))
    all(isfinite, values) || throw(ArgumentError("link weights must be finite"))
    all(x -> x >= 0, values) || throw(ArgumentError("link weights must be nonnegative"))
    normsq = sum(abs2, values)
    normsq > 0 || throw(ArgumentError("link weights must not all be zero"))

    entropy = 0.0
    for value in values
        probability = abs2(value) / normsq
        if probability > 0
            entropy -= probability * log(probability)
        end
    end
    return entropy
end

"""
    bond_entropy(psi, c, dir)

Return `weight_entropy(link_weight(psi, c, dir))` for the periodic
nearest-neighbor bond from `c` in direction `dir`.
"""
bond_entropy(psi::SquareIPEPSState, c::SquareCoord, dir::Symbol) =
    weight_entropy(link_weight(psi, c, dir))

"""
    all_bond_entropies(psi)

Return a dictionary mapping each canonical `BondKey` in `psi` to the entropy
of its stored simple-update link-weight vector.
"""
function all_bond_entropies(psi::SquareIPEPSState)
    return Dict(key => weight_entropy(lambda) for (key, lambda) in psi.link_weights)
end

function _square_star_dense_index(values)
    length(values) == SQUARE_STAR_SITES ||
        throw(ArgumentError("square-star basis value tuple must have 5 sites"))
    idx = 1
    for (site, value) in enumerate(values)
        1 <= value <= 2 || throw(ArgumentError("square-star basis values must be 1 or 2"))
        idx += (value - 1) * 2^(SQUARE_STAR_SITES - site)
    end
    return idx
end

function _square_star_indices(site_indices)
    sites = collect(site_indices)
    length(sites) == SQUARE_STAR_SITES || throw(
        ArgumentError(
            "square-star gate requires 5 physical indices in (center, right, up, left, down) order",
        ),
    )
    all(i -> dim(i) == 2, sites) ||
        throw(ArgumentError("square-star physical indices must all have dimension 2"))
    return sites
end

function _square_star_gate_itensor(dense_gate::AbstractMatrix, site_indices)
    size(dense_gate) == (2^SQUARE_STAR_SITES, 2^SQUARE_STAR_SITES) ||
        throw(ArgumentError("dense square-star gate must be 32x32"))
    sites = _square_star_indices(site_indices)
    out = prime.(sites)
    data = zeros(ComplexF64, ntuple(Returns(2), 2 * SQUARE_STAR_SITES))

    basis_values = Iterators.product((1:2 for _ = 1:SQUARE_STAR_SITES)...)
    for out_values in basis_values
        out_idx = _square_star_dense_index(out_values)
        for in_values in Iterators.product((1:2 for _ = 1:SQUARE_STAR_SITES)...)
            in_idx = _square_star_dense_index(in_values)
            data[out_values..., in_values...] = dense_gate[out_idx, in_idx]
        end
    end

    return ITensor(data, out..., sites...)
end

"""
    square_pxp_gate_itensor(site_indices, step; evolution = :real)

Wrap `square_pxp_gate(step; evolution)` as an ITensor with primed output
indices followed by unprimed input indices. `site_indices` must contain five
dimension-2 indices in dense square-star order `(center, right, up, left, down)`.

Example:

```julia
sites = (center, right, up, left, down)
gate = square_pxp_gate_itensor(sites, 0.01; evolution = :real)
```
"""
function square_pxp_gate_itensor(site_indices, step::Real; evolution::Symbol = :real)
    return _square_star_gate_itensor(square_pxp_gate(step; evolution), site_indices)
end

function square_pxp_gate_itensor(
    step::Real,
    site_indices::Union{Tuple,AbstractVector};
    evolution::Symbol = :real,
)
    return square_pxp_gate_itensor(site_indices, step; evolution)
end

"""
    projected_square_pxp_gate_itensor(site_indices, step; evolution = :real)

Wrap `projected_square_pxp_gate(step; evolution)` as an ITensor with primed
output indices followed by unprimed input indices. `site_indices` must be in
dense square-star order `(center, right, up, left, down)`.

Example:

```julia
sites = (center, right, up, left, down)
gate = projected_square_pxp_gate_itensor(sites, 0.01; evolution = :real)
```
"""
function projected_square_pxp_gate_itensor(
    site_indices,
    step::Real;
    evolution::Symbol = :real,
)
    return _square_star_gate_itensor(
        projected_square_pxp_gate(step; evolution),
        site_indices,
    )
end

function projected_square_pxp_gate_itensor(
    step::Real,
    site_indices::Union{Tuple,AbstractVector};
    evolution::Symbol = :real,
)
    return projected_square_pxp_gate_itensor(site_indices, step; evolution)
end

end
