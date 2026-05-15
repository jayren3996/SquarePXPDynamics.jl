module SquarePEPS

using ITensors
using ..SquareGeometry

export SquarePEPSState, product_square_peps, site_tensor, physical_index, link_index

"""
    SquarePEPSState

Finite square-lattice PEPS state with one ITensor per site, physical indices,
nearest-neighbor link indices, and a target maximum link dimension.
"""
struct SquarePEPSState
    sites::Vector{SquareCoord}
    tensors::Dict{SquareCoord,ITensor}
    physical_indices::Dict{SquareCoord,Index}
    link_indices::Dict{Tuple{SquareCoord,Symbol},Index}
    maxdim::Int
end

function _state_vector(state::Symbol)
    if state === :up
        return ComplexF64[1, 0]
    elseif state === :down
        return ComplexF64[0, 1]
    else
        throw(ArgumentError("state must be :up or :down"))
    end
end

function _site_grid(width::Integer, height::Integer)
    width >= 1 || throw(ArgumentError("width must be positive"))
    height >= 1 || throw(ArgumentError("height must be positive"))
    return [SquareCoord(x, y) for y = 1:Int(height) for x = 1:Int(width)]
end

function _has_site(site_set, c::SquareCoord)
    return c in site_set
end

function _build_link_indices(sites::Vector{SquareCoord}, maxdim::Integer)
    maxdim >= 1 || throw(ArgumentError("maxdim must be positive"))
    site_set = Set(sites)
    links = Dict{Tuple{SquareCoord,Symbol},Index}()
    for c in sites
        right = square_neighbor(c, :right)
        up = square_neighbor(c, :up)
        if _has_site(site_set, right)
            idx = Index(Int(maxdim), "link,$(c.x),$(c.y),right")
            links[(c, :right)] = idx
            links[(right, :left)] = idx
        end
        if _has_site(site_set, up)
            idx = Index(Int(maxdim), "link,$(c.x),$(c.y),up")
            links[(c, :up)] = idx
            links[(up, :down)] = idx
        end
    end
    return links
end

function _link_or_boundary(links, c::SquareCoord, direction::Symbol)
    return get(links, (c, direction), Index(1, "boundary,$(c.x),$(c.y),$direction"))
end

"""
    product_square_peps(width, height; state = :down, maxdim = 1)

Construct a finite `width` by `height` square PEPS product state. `state` may be
`:up` or `:down`, using the convention `|up> = |0>` and `|down> = |1>`.
Interior links have dimension `maxdim`; open-boundary links have dimension `1`.
"""
function product_square_peps(
    width::Integer,
    height::Integer;
    state::Symbol = :down,
    maxdim::Integer = 1,
)
    sites = _site_grid(width, height)
    links = _build_link_indices(sites, maxdim)
    physical = Dict(c => Index(2, "phys,$(c.x),$(c.y)") for c in sites)
    tensors = Dict{SquareCoord,ITensor}()
    amplitudes = _state_vector(state)

    for c in sites
        p = physical[c]
        left = _link_or_boundary(links, c, :left)
        right = _link_or_boundary(links, c, :right)
        up = _link_or_boundary(links, c, :up)
        down = _link_or_boundary(links, c, :down)
        tensor = ITensor(ComplexF64, p, left, right, up, down)
        for s = 1:2
            tensor[p=>s, left=>1, right=>1, up=>1, down=>1] = amplitudes[s]
        end
        tensors[c] = tensor
    end

    return SquarePEPSState(sites, tensors, physical, links, Int(maxdim))
end

"""
    site_tensor(psi, c)

Return the ITensor stored at square-lattice coordinate `c`.
"""
site_tensor(psi::SquarePEPSState, c::SquareCoord) = psi.tensors[c]

"""
    physical_index(psi, c)

Return the physical ITensors index for square-lattice coordinate `c`.
"""
physical_index(psi::SquarePEPSState, c::SquareCoord) = psi.physical_indices[c]

"""
    link_index(psi, c, direction)

Return the nearest-neighbor link index at site `c` in `direction`.
"""
link_index(psi::SquarePEPSState, c::SquareCoord, direction::Symbol) =
    psi.link_indices[(c, direction)]

end
