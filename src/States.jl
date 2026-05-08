module States

using ITensors
using Random
using ..Geometry: Coord, neighbor

export AbstractUnitCell, OneSiteUnitCell, ThreeSiteUnitCell
export TriangularIPEPS
export unit_cell_representatives, wrap_coord
export product_ipeps, random_ipeps
export site_tensor, phys_index, bond_index, bond_indices, bond_lambda
export opposite_direction

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
        T = ITensor(ComplexF64, ph, binds...)
        data = randn(rng, ComplexF64, dim(ph), ntuple(d -> dim(binds[d]), 6)...)
        T .= ITensor(data, ph, binds...)
        tensors[c] = T
    end
    return TriangularIPEPS(uc, phys_inds, bond_inds, tensors, lambdas)
end

end
