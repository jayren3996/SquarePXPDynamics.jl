@testset "geometry" begin
    c = Coord(2, -1)
    @test neighbor(c, 1) == Coord(3, -1)
    @test neighbor(c, 4) == Coord(1, -1)
    @test triangular_distance(Coord(0, 0), Coord(2, -1)) == 2

    s = star_sites(Coord(0, 0))
    @test length(s) == 7
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
