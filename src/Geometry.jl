module Geometry

export Coord, TRIANGULAR_DIRECTIONS, triangular_distance, neighbor, star_sites
export star_color, disjoint_stars

struct Coord
    q::Int
    r::Int
end

const TRIANGULAR_DIRECTIONS = (
    Coord(1, 0),
    Coord(0, 1),
    Coord(-1, 1),
    Coord(-1, 0),
    Coord(0, -1),
    Coord(1, -1),
)

Base.:+(a::Coord, b::Coord) = Coord(a.q + b.q, a.r + b.r)
Base.:-(a::Coord, b::Coord) = Coord(a.q - b.q, a.r - b.r)
Base.:(==)(a::Coord, b::Coord) = a.q == b.q && a.r == b.r
Base.hash(c::Coord, h::UInt) = hash((c.q, c.r), h)

function triangular_distance(a::Coord, b::Coord)
    d = b - a
    return (abs(d.q) + abs(d.r) + abs(d.q + d.r)) ÷ 2
end

function neighbor(c::Coord, direction::Integer)
    1 <= direction <= 6 || throw(ArgumentError("direction must be in 1:6"))
    return c + TRIANGULAR_DIRECTIONS[direction]
end

function star_sites(center::Coord)
    return [center; [neighbor(center, dir) for dir in 1:6]]
end

function star_color(c::Coord)
    return mod(c.q + 3c.r, 7) + 1
end

function disjoint_stars(a::Coord, b::Coord)
    return isempty(intersect(Set(star_sites(a)), Set(star_sites(b))))
end

end
