@testset "geometry" begin
    c = Coord(2, -1)
    @test neighbor(c, 1) == Coord(3, -1)
    @test neighbor(c, 4) == Coord(1, -1)
    @test triangular_distance(Coord(0, 0), Coord(2, -1)) == 2
    @test all(triangular_distance(Coord(0, 0), dir) == 1 for dir in TRIANGULAR_DIRECTIONS)
    @test triangular_distance(Coord(0, 0), Coord(1, 1)) == 2
    @test triangular_distance(Coord(-2, 3), Coord(1, -1)) ==
        triangular_distance(Coord(1, -1), Coord(-2, 3))

    radius_2_distances = Dict(
        Coord(0, 0) => 0,
        Coord(1, 0) => 1,
        Coord(0, 1) => 1,
        Coord(-1, 1) => 1,
        Coord(-1, 0) => 1,
        Coord(0, -1) => 1,
        Coord(1, -1) => 1,
        Coord(2, 0) => 2,
        Coord(1, 1) => 2,
        Coord(0, 2) => 2,
        Coord(-1, 2) => 2,
        Coord(-2, 2) => 2,
        Coord(-2, 1) => 2,
        Coord(-2, 0) => 2,
        Coord(-1, -1) => 2,
        Coord(0, -2) => 2,
        Coord(1, -2) => 2,
        Coord(2, -2) => 2,
        Coord(2, -1) => 2,
    )
    for (site, distance) in radius_2_distances
        @test triangular_distance(Coord(0, 0), site) == distance
        @test triangular_distance(site, Coord(0, 0)) == distance
    end

    s = star_sites(Coord(0, 0))
    expected_star = [Coord(0, 0); collect(TRIANGULAR_DIRECTIONS)]
    @test length(s) == 7
    @test s == expected_star
    @test length(Set(s)) == 7
    @test Set(s) == Set(expected_star)
    @test Coord(0, 0) in s
    @test Coord(1, 0) in s
    @test Coord(0, 1) in s

    centers = [Coord(q, r) for q in -4:4 for r in -4:4]
    for a in centers, b in centers
        if a != b && star_color(a) == star_color(b)
            @test disjoint_stars(a, b)
        end
    end
end
