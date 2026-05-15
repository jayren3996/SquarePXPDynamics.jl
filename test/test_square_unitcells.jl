@testset "periodic square unit cell wrapping" begin
    cell = PeriodicSquareUnitCell(10, 10)

    @test wrap(cell, SquareCoord(11, 1)) == SquareCoord(1, 1)
    @test wrap(cell, SquareCoord(0, 1)) == SquareCoord(10, 1)
    @test wrap(cell, SquareCoord(1, 11)) == SquareCoord(1, 1)
    @test wrap(cell, SquareCoord(1, 0)) == SquareCoord(1, 10)

    @test neighbor(cell, SquareCoord(10, 1), :right) == SquareCoord(1, 1)
    @test neighbor(cell, SquareCoord(1, 1), :left) == SquareCoord(10, 1)
    @test neighbor(cell, SquareCoord(1, 10), :up) == SquareCoord(1, 1)
    @test neighbor(cell, SquareCoord(1, 1), :down) == SquareCoord(1, 10)

    @test length(cell.reps) == 100
end

@testset "periodic square unit cell validation" begin
    @test_throws ArgumentError PeriodicSquareUnitCell(0, 10)
    @test_throws ArgumentError PeriodicSquareUnitCell(10, 0)
end

@testset "five-color schedule compatibility" begin
    cell = PeriodicSquareUnitCell(10, 10)
    @test assert_five_color_compatible(cell) === cell

    for color in 1:5
        centers = update_centers(cell, color)
        @test !isempty(centers)
        @test stars_are_disjoint_mod_unitcell(cell, centers)
    end

    @test_throws ArgumentError update_centers(cell, 0)
    @test_throws ArgumentError update_centers(cell, 6)
    @test_throws ArgumentError assert_five_color_compatible(PeriodicSquareUnitCell(4, 4))
end

@testset "canonical bond keys" begin
    cell = PeriodicSquareUnitCell(10, 10)
    c = SquareCoord(1, 1)

    @test bondkey(cell, c, :right) == BondKey(SquareCoord(1, 1), :right)
    @test bondkey(cell, c, :up) == BondKey(SquareCoord(1, 1), :up)
    @test bondkey(cell, c, :left) == BondKey(SquareCoord(10, 1), :right)
    @test bondkey(cell, c, :down) == BondKey(SquareCoord(1, 10), :up)

    @test_throws ArgumentError BondKey(c, :left)
    @test_throws ArgumentError bondkey(cell, c, :diagonal)
end
