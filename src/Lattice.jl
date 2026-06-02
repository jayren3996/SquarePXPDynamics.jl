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
square_star_sites(c::SquareCoord) =
    (c, (square_neighbor(c, d) for d in SQUARE_DIRECTIONS)...)

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
unit_cell_representatives(::FiveSiteSquareUC) = [SquareCoord(k, 0) for k = 0:4]

"""
    wrap_square_coord(unit_cell, c)

Map a square-lattice coordinate to its unit-cell representative.
"""
wrap_square_coord(::OneSiteSquareUC, c::SquareCoord) = SquareCoord(0, 0)
wrap_square_coord(::FiveSiteSquareUC, c::SquareCoord) =
    SquareCoord(square_star_color(c) - 1, 0)

end


module SquareUnitCells

using ..SquareGeometry

export PeriodicSquareUnitCell
export wrap, neighbor, update_centers, assert_five_color_compatible
export stars_are_disjoint_mod_unitcell
export BondKey, bondkey

"""
    PeriodicSquareUnitCell(Lx, Ly)

Periodic rectangular square-lattice unit cell with one-based representatives
ordered as `[SquareCoord(x, y) for y in 1:Ly for x in 1:Lx]`.
"""
struct PeriodicSquareUnitCell <: SquareUnitCell
    Lx::Int
    Ly::Int
    reps::Vector{SquareCoord}

    function PeriodicSquareUnitCell(Lx::Integer, Ly::Integer)
        Lx >= 1 || throw(ArgumentError("Lx must be at least 1"))
        Ly >= 1 || throw(ArgumentError("Ly must be at least 1"))
        nx = Int(Lx)
        ny = Int(Ly)
        reps = [SquareCoord(x, y) for y = 1:ny for x = 1:nx]
        return new(nx, ny, reps)
    end
end

"""
    wrap(cell, c)

Wrap square-lattice coordinate `c` into the one-based representatives of
periodic unit cell `cell`.
"""
wrap(cell::PeriodicSquareUnitCell, c::SquareCoord) =
    SquareCoord(mod1(c.x, cell.Lx), mod1(c.y, cell.Ly))

"""
    neighbor(cell, c, dir)

Return the periodic nearest neighbor of `c` in `:right`, `:up`, `:left`, or
`:down`, wrapping the result into `cell`.
"""
neighbor(cell::PeriodicSquareUnitCell, c::SquareCoord, dir::Symbol) =
    wrap(cell, square_neighbor(c, dir))

"""
    update_centers(cell, color)

Return representatives whose square-star schedule color is `color`, where
`color` must be in `1:5`.
"""
function update_centers(cell::PeriodicSquareUnitCell, color::Integer)
    1 <= color <= 5 || throw(ArgumentError("color must be in 1:5"))
    return [c for c in cell.reps if square_star_color(c) == Int(color)]
end

"""
    stars_are_disjoint_mod_unitcell(cell, centers)

Return whether the wrapped five-site square stars centered at `centers` share
no representative sites inside periodic unit cell `cell`.
"""
function stars_are_disjoint_mod_unitcell(cell::PeriodicSquareUnitCell, centers)
    seen = Set{SquareCoord}()
    for center in centers
        for site in square_star_sites(center)
            wrapped = wrap(cell, site)
            wrapped in seen && return false
            push!(seen, wrapped)
        end
    end
    return true
end

"""
    assert_five_color_compatible(cell)

Validate that `cell` is compatible with the five-color square-star update
schedule. Both dimensions must be multiples of five, and each color layer must
have disjoint wrapped stars. Returns `cell` on success.
"""
function assert_five_color_compatible(cell::PeriodicSquareUnitCell)
    cell.Lx % 5 == 0 || throw(ArgumentError("Lx must be a multiple of 5"))
    cell.Ly % 5 == 0 || throw(ArgumentError("Ly must be a multiple of 5"))
    for color = 1:5
        stars_are_disjoint_mod_unitcell(cell, update_centers(cell, color)) ||
            throw(ArgumentError("color layer $color has overlapping wrapped stars"))
    end
    return cell
end

"""
    BondKey(site, dir)

Canonical key for an undirected periodic nearest-neighbor bond. `dir` must be
`:right` or `:up`; left and down bonds are represented by the neighboring
right or up key.
"""
struct BondKey
    site::SquareCoord
    dir::Symbol

    function BondKey(site::SquareCoord, dir::Symbol)
        dir === :right ||
            dir === :up ||
            throw(ArgumentError("BondKey direction must be :right or :up"))
        return new(site, dir)
    end
end

"""
    bondkey(cell, c, dir)

Return the canonical `BondKey` for the periodic nearest-neighbor bond from
`c` in `dir`. Right and up bonds are stored at `c`; left and down bonds are
stored as the corresponding neighbor's right or up bond.
"""
function bondkey(cell::PeriodicSquareUnitCell, c::SquareCoord, dir::Symbol)
    if dir === :right
        return BondKey(wrap(cell, c), :right)
    elseif dir === :up
        return BondKey(wrap(cell, c), :up)
    elseif dir === :left
        return BondKey(neighbor(cell, c, :left), :right)
    elseif dir === :down
        return BondKey(neighbor(cell, c, :down), :up)
    else
        throw(ArgumentError("direction must be :right, :up, :left, or :down"))
    end
end

end
