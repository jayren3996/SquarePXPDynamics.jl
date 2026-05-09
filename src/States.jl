module States

using ITensors
using Random
using ..Geometry: Coord, neighbor, star_color

export AbstractUnitCell, OneSiteUnitCell, ThreeSiteUnitCell, SevenSiteUnitCell
export TriangularIPEPS
export unit_cell_representatives, wrap_coord
export product_ipeps, random_ipeps
export site_tensor, phys_index, bond_index, bond_indices, bond_lambda
export opposite_direction
export StateTruncationDiagnostics, truncate_state!

"""
    opposite_direction(d) -> Int

Return the index of the direction antiparallel to `d` in `TRIANGULAR_DIRECTIONS`.
Directions 1↔4, 2↔5, 3↔6.
"""
opposite_direction(d::Integer) = mod(d - 1 + 3, 6) + 1

abstract type AbstractUnitCell end

struct OneSiteUnitCell <: AbstractUnitCell end

unit_cell_representatives(::OneSiteUnitCell) = (Coord(0, 0),)
wrap_coord(::OneSiteUnitCell, ::Coord) = Coord(0, 0)

"""
    ThreeSiteUnitCell()

Three-sublattice partition of the triangular lattice using `(q - r) mod 3`.
Representatives: `(0,0)`, `(1,0)`, `(2,0)`.
"""
struct ThreeSiteUnitCell <: AbstractUnitCell end

const _THREE_REPS = (Coord(0, 0), Coord(1, 0), Coord(2, 0))

unit_cell_representatives(::ThreeSiteUnitCell) = _THREE_REPS

function wrap_coord(::ThreeSiteUnitCell, c::Coord)
    s = mod(c.q - c.r, 3)
    return _THREE_REPS[s + 1]
end

"""
    SevenSiteUnitCell()

Seven-sublattice partition aligned with the 7-color triangular star schedule.
Representatives are `(0,0)` through `(6,0)`, matching `star_color(c)`.
"""
struct SevenSiteUnitCell <: AbstractUnitCell end

const _SEVEN_REPS = (Coord(0, 0), Coord(1, 0), Coord(2, 0), Coord(3, 0),
                     Coord(4, 0), Coord(5, 0), Coord(6, 0))

unit_cell_representatives(::SevenSiteUnitCell) = _SEVEN_REPS
wrap_coord(::SevenSiteUnitCell, c::Coord) = _SEVEN_REPS[star_color(c)]

"""
    TriangularIPEPS

Native triangular iPEPS state container. Each representative tensor has
indices `(phys, n1, ..., n6)` where `n_d` is the bond index in direction `d`.

Fields:
- `unitcell`: `AbstractUnitCell` describing translational structure.
- `phys_inds`: physical `Index` per representative `Coord`.
- `bond_inds`: bond `Index` per `(rep, dir)`. Opposite-direction bonds across
  representatives share the same `Index` object.
- `tensors`: site `ITensor` per representative.
- `lambdas`: diagonal bond spectrum (Vector{Float64}) per `(rep, dir)`.
  Opposite-direction entries reference the same vector.
"""
struct TriangularIPEPS{UC<:AbstractUnitCell}
    unitcell::UC
    phys_inds::Dict{Coord,Index}
    bond_inds::Dict{Tuple{Coord,Int},Index}
    tensors::Dict{Coord,ITensor}
    lambdas::Dict{Tuple{Coord,Int},Vector{Float64}}
end

struct StateTruncationDiagnostics
    target_maxdim::Int
    discarded_weight_by_bond::Dict{Tuple{Coord,Int},Float64}
    output_bond_dims::Dict{Tuple{Coord,Int},Int}
    max_output_bond_dim::Int
end

site_tensor(state::TriangularIPEPS, c::Coord) = state.tensors[wrap_coord(state.unitcell, c)]
phys_index(state::TriangularIPEPS, c::Coord) = state.phys_inds[wrap_coord(state.unitcell, c)]
bond_index(state::TriangularIPEPS, c::Coord, d::Integer) =
    state.bond_inds[(wrap_coord(state.unitcell, c), Int(d))]
bond_lambda(state::TriangularIPEPS, c::Coord, d::Integer) =
    state.lambdas[(wrap_coord(state.unitcell, c), Int(d))]

function bond_indices(state::TriangularIPEPS, c::Coord)
    rep = wrap_coord(state.unitcell, c)
    return ntuple(d -> state.bond_inds[(rep, d)], 6)
end

function _build_indices(uc::AbstractUnitCell, D::Integer)
    phys_inds = Dict{Coord,Index}()
    bond_inds = Dict{Tuple{Coord,Int},Index}()
    lambdas = Dict{Tuple{Coord,Int},Vector{Float64}}()

    for c in unit_cell_representatives(uc)
        phys_inds[c] = Index(2, "phys,$(c.q),$(c.r)")
    end

    # Allocate one Index per (rep, dir). When the opposite-direction bond on a
    # *distinct* representative is the same logical bond, share the Index; for
    # self-loops (e.g. 1-site UC where every neighbor wraps back to the same
    # rep), give each direction its own unique Index so a tensor never repeats
    # an Index.
    for c in unit_cell_representatives(uc)
        for d in 1:6
            haskey(bond_inds, (c, d)) && continue
            opp_c = wrap_coord(uc, neighbor(c, d))
            opp_d = opposite_direction(d)
            idx = Index(D, "bond,d=$d,from=$(c.q),$(c.r)")
            bond_inds[(c, d)] = idx
            if opp_c != c && !haskey(bond_inds, (opp_c, opp_d))
                bond_inds[(opp_c, opp_d)] = idx
            end
            lambda = get(lambdas, (c, d), nothing)
            if lambda === nothing
                lambda = ones(Float64, D)
                lambdas[(c, d)] = lambda
            end
            if !haskey(lambdas, (opp_c, opp_d))
                lambdas[(opp_c, opp_d)] = lambda
            end
        end
    end

    return phys_inds, bond_inds, lambdas
end

"""
    product_ipeps(uc, state_symbol; D=1)

Build a product-state iPEPS where each representative is in the local
state given by `state_symbol`. Currently supports `:up` (|0>) and `:down` (|1>).
"""
function product_ipeps(uc::AbstractUnitCell, state_symbol::Symbol; D::Integer = 1)
    D >= 1 || throw(ArgumentError("D must be >= 1"))
    phys_inds, bond_inds, lambdas = _build_indices(uc, D)

    local_vec = if state_symbol === :up
        ComplexF64[1, 0]
    elseif state_symbol === :down
        ComplexF64[0, 1]
    else
        throw(ArgumentError("state_symbol must be :up or :down"))
    end

    tensors = Dict{Coord,ITensor}()
    for c in unit_cell_representatives(uc)
        ph = phys_inds[c]
        binds = ntuple(d -> bond_inds[(c, d)], 6)
        T = ITensor(ComplexF64, ph, binds...)
        # Place |state> at the (1,1,1,1,1,1) bond corner.
        for k in 1:2
            T[ph => k, (binds[d] => 1 for d in 1:6)...] = local_vec[k]
        end
        tensors[c] = T
    end
    return TriangularIPEPS(uc, phys_inds, bond_inds, tensors, lambdas)
end

"""
    random_ipeps(uc, D; seed=nothing)

Build an iPEPS with random Gaussian (complex) entries per representative.
"""
function random_ipeps(uc::AbstractUnitCell, D::Integer; seed::Union{Nothing,Integer} = nothing)
    D >= 1 || throw(ArgumentError("D must be >= 1"))
    phys_inds, bond_inds, lambdas = _build_indices(uc, D)
    rng = seed === nothing ? Random.default_rng() : Random.MersenneTwister(seed)
    tensors = Dict{Coord,ITensor}()
    for c in unit_cell_representatives(uc)
        ph = phys_inds[c]
        binds = ntuple(d -> bond_inds[(c, d)], 6)
        data = randn(rng, ComplexF64, dim(ph), ntuple(d -> dim(binds[d]), 6)...)
        tensors[c] = ITensor(data, ph, binds...)
    end
    return TriangularIPEPS(uc, phys_inds, bond_inds, tensors, lambdas)
end

function truncate_state!(state::TriangularIPEPS, target_maxdim::Integer)
    target_maxdim >= 1 || throw(ArgumentError("target_maxdim must be >= 1"))

    processed = Set{Tuple{Coord,Int}}()
    discarded = Dict{Tuple{Coord,Int},Float64}()
    output_dims = Dict{Tuple{Coord,Int},Int}()

    for rep in unit_cell_representatives(state.unitcell)
        for d in 1:6
            key = (rep, d)
            key in processed && continue

            opp_rep = wrap_coord(state.unitcell, neighbor(rep, d))
            opp_d = opposite_direction(d)
            opp_key = (opp_rep, opp_d)

            λ_old = state.lambdas[key]
            current_dim = length(λ_old)
            if target_maxdim > current_dim
                throw(ArgumentError(
                    "target_maxdim $target_maxdim exceeds current bond dimension $current_dim"))
            end

            order = sort(collect(eachindex(λ_old)); by = i -> abs(λ_old[i]), rev = true)
            keep = order[1:target_maxdim]
            drop = target_maxdim == current_dim ? Int[] : order[(target_maxdim + 1):end]
            total_weight = sum(abs2, λ_old)
            discarded_weight = total_weight == 0 ? 0.0 : sum(abs2, λ_old[drop]) / total_weight
            λ_new = Float64.(λ_old[keep])

            old_idx = state.bond_inds[key]
            new_idx = Index(target_maxdim, "bond,d=$d,from=$(rep.q),$(rep.r),truncated")
            state.tensors[rep] = _replace_tensor_index_slice(state.tensors[rep], old_idx, new_idx, keep)
            state.bond_inds[key] = new_idx

            if opp_key != key
                old_opp_idx = state.bond_inds[opp_key]
                if old_opp_idx === old_idx
                    state.tensors[opp_rep] = _replace_tensor_index_slice(
                        state.tensors[opp_rep], old_idx, new_idx, keep)
                    state.bond_inds[opp_key] = new_idx
                else
                    new_opp_idx = Index(target_maxdim,
                                        "bond,d=$opp_d,from=$(opp_rep.q),$(opp_rep.r),truncated")
                    state.tensors[opp_rep] = _replace_tensor_index_slice(
                        state.tensors[opp_rep], old_opp_idx, new_opp_idx, keep)
                    state.bond_inds[opp_key] = new_opp_idx
                end
            end

            state.lambdas[key] = λ_new
            if opp_key != key && state.lambdas[opp_key] === λ_old
                state.lambdas[opp_key] = λ_new
            end

            discarded[key] = Float64(discarded_weight)
            output_dims[key] = target_maxdim
            push!(processed, key)
            push!(processed, opp_key)
        end
    end

    maxdim = isempty(output_dims) ? 0 : maximum(values(output_dims))
    return StateTruncationDiagnostics(Int(target_maxdim), discarded, output_dims, maxdim)
end

function _replace_tensor_index_slice(T::ITensor, old_idx::Index, new_idx::Index,
                                     keep::Vector{Int})
    idxs = collect(inds(T))
    pos = findfirst(i -> i == old_idx, idxs)
    pos === nothing && error("bond index to truncate was not found in tensor")
    selectors = ntuple(i -> i == pos ? keep : Colon(), length(idxs))
    data = array(T)[selectors...]
    idxs[pos] = new_idx
    return ITensor(data, idxs...)
end

end
