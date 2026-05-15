@testset "square geometry" begin
    c = SquareCoord(0, 0)

    @test square_neighbor(c, :right) == SquareCoord(1, 0)
    @test square_neighbor(c, :up) == SquareCoord(0, 1)
    @test square_neighbor(c, :left) == SquareCoord(-1, 0)
    @test square_neighbor(c, :down) == SquareCoord(0, -1)
    @test square_neighbor(c, 1) == SquareCoord(1, 0)
    @test_throws ArgumentError square_neighbor(c, :north)
    @test_throws ArgumentError square_neighbor(c, 5)

    star = square_star_sites(c)
    @test star == (
        SquareCoord(0, 0),
        SquareCoord(1, 0),
        SquareCoord(0, 1),
        SquareCoord(-1, 0),
        SquareCoord(0, -1),
    )
    @test length(unique(star)) == 5

    centers = [SquareCoord(x, y) for x = -4:4 for y = -4:4]
    for color = 1:5
        same = filter(site -> square_star_color(site) == color, centers)
        for i in eachindex(same), j = (i+1):lastindex(same)
            @test disjoint_square_stars(same[i], same[j])
        end
    end

    uc = FiveSiteSquareUC()
    @test length(unit_cell_representatives(uc)) == 5
    @test wrap_square_coord(uc, SquareCoord(3, 7)) ==
          SquareCoord(square_star_color(SquareCoord(3, 7)) - 1, 0)
    @test wrap_square_coord(OneSiteSquareUC(), SquareCoord(3, 7)) == SquareCoord(0, 0)
end
