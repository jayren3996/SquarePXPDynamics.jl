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
