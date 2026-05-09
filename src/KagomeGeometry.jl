module KagomeGeometry

export KagomeCoord, KagomeTriangleCoord
export KagomeUnitCell, NineSiteKagomeUC
export kagome_neighbor, kagome_star_sites, kagome_star_color, disjoint_kagome_stars
export up_triangle_of, down_triangle_of, triangle_sites
export unit_cell_representatives, wrap_kagome_coord

const _SUBLATTICES = (:A, :B, :C)

struct KagomeCoord
    n1::Int
    n2::Int
    sublat::Symbol

    function KagomeCoord(n1::Integer, n2::Integer, sublat::Symbol)
        sublat in _SUBLATTICES || throw(ArgumentError("sublat must be :A, :B, or :C"))
        return new(Int(n1), Int(n2), sublat)
    end
end

Base.:(==)(a::KagomeCoord, b::KagomeCoord) =
    a.n1 == b.n1 && a.n2 == b.n2 && a.sublat === b.sublat
Base.hash(c::KagomeCoord, h::UInt) = hash((c.n1, c.n2, c.sublat), h)

struct KagomeTriangleCoord
    n1::Int
    n2::Int
    orientation::Symbol

    function KagomeTriangleCoord(n1::Integer, n2::Integer, orientation::Symbol)
        orientation in (:up, :down) || throw(ArgumentError("orientation must be :up or :down"))
        return new(Int(n1), Int(n2), orientation)
    end
end

Base.:(==)(a::KagomeTriangleCoord, b::KagomeTriangleCoord) =
    a.n1 == b.n1 && a.n2 == b.n2 && a.orientation === b.orientation
Base.hash(t::KagomeTriangleCoord, h::UInt) = hash((t.n1, t.n2, t.orientation), h)

up_triangle_of(c::KagomeCoord) = KagomeTriangleCoord(c.n1, c.n2, :up)

function down_triangle_of(c::KagomeCoord)
    if c.sublat === :A
        return KagomeTriangleCoord(c.n1, c.n2, :down)
    elseif c.sublat === :B
        return KagomeTriangleCoord(c.n1 + 1, c.n2, :down)
    else
        return KagomeTriangleCoord(c.n1, c.n2 + 1, :down)
    end
end

function triangle_sites(t::KagomeTriangleCoord)
    if t.orientation === :up
        return (
            KagomeCoord(t.n1, t.n2, :A),
            KagomeCoord(t.n1, t.n2, :B),
            KagomeCoord(t.n1, t.n2, :C),
        )
    else
        return (
            KagomeCoord(t.n1, t.n2, :A),
            KagomeCoord(t.n1 - 1, t.n2, :B),
            KagomeCoord(t.n1, t.n2 - 1, :C),
        )
    end
end

function kagome_neighbor(c::KagomeCoord, d::Integer)
    1 <= d <= 4 || throw(ArgumentError("direction must be in 1:4"))
    up_sites = filter(!=(c), triangle_sites(up_triangle_of(c)))
    down_sites = filter(!=(c), triangle_sites(down_triangle_of(c)))
    neighbors = (up_sites..., down_sites...)
    length(neighbors) == 4 || error("expected 4 kagome neighbors, got $(length(neighbors))")
    return neighbors[Int(d)]
end

kagome_star_sites(c::KagomeCoord) = (c, (kagome_neighbor(c, d) for d in 1:4)...)

kagome_star_color(c::KagomeCoord) = findfirst(==(c.sublat), _SUBLATTICES)

function disjoint_kagome_stars(a::KagomeCoord, b::KagomeCoord)
    a == b && return false
    kagome_star_color(a) == kagome_star_color(b) && return true
    return isempty(intersect(Set(kagome_star_sites(a)), Set(kagome_star_sites(b))))
end

abstract type KagomeUnitCell end

struct NineSiteKagomeUC <: KagomeUnitCell end

function unit_cell_representatives(::NineSiteKagomeUC)
    return [KagomeCoord(k, 0, s) for k in 0:2 for s in _SUBLATTICES]
end

function wrap_kagome_coord(::NineSiteKagomeUC, c::KagomeCoord)
    cell = mod(c.n1 + 2 * c.n2, 3)
    return KagomeCoord(cell, 0, c.sublat)
end

end
