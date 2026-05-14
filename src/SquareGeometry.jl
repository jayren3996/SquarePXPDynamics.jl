module SquareGeometry

export SquareCoord, SquareUnitCell, OneSiteSquareUC, FiveSiteSquareUC
export square_neighbor, square_star_sites, square_star_color, disjoint_square_stars
export unit_cell_representatives, wrap_square_coord

const SQUARE_DIRECTIONS = (:right, :up, :left, :down)

"""
    SquareCoord(x, y)

Integer coordinate for a square-lattice site.
"""
struct SquareCoord
    x::Int
    y::Int
end

SquareCoord(x::Integer, y::Integer) = SquareCoord(Int(x), Int(y))

"""
    square_neighbor(c, direction)

Return the nearest neighbor of `c` in `:right`, `:up`, `:left`, or `:down`.
Integer directions `1:4` use that same order.
"""
function square_neighbor(c::SquareCoord, direction::Symbol)
    if direction === :right
        return SquareCoord(c.x + 1, c.y)
    elseif direction === :up
        return SquareCoord(c.x, c.y + 1)
    elseif direction === :left
        return SquareCoord(c.x - 1, c.y)
    elseif direction === :down
        return SquareCoord(c.x, c.y - 1)
    else
        throw(ArgumentError("direction must be :right, :up, :left, or :down"))
    end
end

function square_neighbor(c::SquareCoord, direction::Integer)
    1 <= direction <= 4 || throw(ArgumentError("direction must be in 1:4"))
    return square_neighbor(c, SQUARE_DIRECTIONS[Int(direction)])
end

"""
    square_star_sites(c)

Return the dense square-star sites in the convention `(center, right, up, left, down)`.
"""
square_star_sites(c::SquareCoord) = (c, (square_neighbor(c, d) for d in SQUARE_DIRECTIONS)...)

"""
    square_star_color(c)

Return the color `1:5` used to schedule non-overlapping square-star updates.
"""
square_star_color(c::SquareCoord) = mod(c.x + 2c.y, 5) + 1

"""
    disjoint_square_stars(a, b)

Return whether the square stars centered at `a` and `b` have no shared sites.
Stars with the same 5-color schedule color are treated as disjoint, except when
the centers are identical.
"""
function disjoint_square_stars(a::SquareCoord, b::SquareCoord)
    a == b && return false
    square_star_color(a) == square_star_color(b) && return true
    return isempty(intersect(Set(square_star_sites(a)), Set(square_star_sites(b))))
end

"""
    SquareUnitCell

Abstract marker type for finite square-lattice unit-cell conventions.
"""
abstract type SquareUnitCell end

"""
    OneSiteSquareUC()

One-site square unit cell where every coordinate wraps to one representative.
"""
struct OneSiteSquareUC <: SquareUnitCell end

"""
    FiveSiteSquareUC()

Five-site square unit cell aligned with the square-star 5-color schedule.
"""
struct FiveSiteSquareUC <: SquareUnitCell end

"""
    unit_cell_representatives(unit_cell)

Return representative coordinates for a square unit-cell convention.
"""
unit_cell_representatives(::OneSiteSquareUC) = [SquareCoord(0, 0)]
unit_cell_representatives(::FiveSiteSquareUC) = [SquareCoord(k, 0) for k in 0:4]

"""
    wrap_square_coord(unit_cell, c)

Map a square-lattice coordinate to its unit-cell representative.
"""
wrap_square_coord(::OneSiteSquareUC, c::SquareCoord) = SquareCoord(0, 0)
wrap_square_coord(::FiveSiteSquareUC, c::SquareCoord) = SquareCoord(square_star_color(c) - 1, 0)

end
