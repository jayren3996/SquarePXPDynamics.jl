@testset "schedules" begin
    @test first_order_colors() == collect(1:7)
    @test second_order_colors() == [1, 2, 3, 4, 5, 6, 7, 7, 6, 5, 4, 3, 2, 1]
    @test schedule_layers(:first) == [(color = c, scale = 1.0) for c in 1:7]
    @test schedule_layers(:second) == vcat(
        [(color = c, scale = 0.5) for c in 1:7],
        [(color = c, scale = 0.5) for c in 7:-1:1],
    )
    @test_throws ArgumentError schedule_layers(:bad)

    centers = [Coord(q, r) for q in -3:3 for r in -3:3]
    for color in first_order_colors()
        layer = [c for c in centers if star_color(c) == color]
        for i in eachindex(layer), j in (i + 1):lastindex(layer)
            @test disjoint_stars(layer[i], layer[j])
        end
    end
end
