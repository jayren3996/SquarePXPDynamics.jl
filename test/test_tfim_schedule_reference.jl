using Test
using LinearAlgebra
using SquarePXPDynamics

function _finite_site_index(x, y, L)
    return mod1(x, L) + (mod1(y, L) - 1) * L
end

function _finite_star_sites(center, L)
    x, y = center
    return (
        _finite_site_index(x, y, L),
        _finite_site_index(x + 1, y, L),
        _finite_site_index(x, y + 1, L),
        _finite_site_index(x - 1, y, L),
        _finite_site_index(x, y - 1, L),
    )
end

function _production_finite_star_sites(center, L)
    cell = PeriodicSquareUnitCell(L, L)
    c = SquareCoord(center...)
    return Tuple(
        _finite_site_index(wrap(cell, s).x, wrap(cell, s).y, L) for
        s in square_star_sites(c)
    )
end

@testset "trotter coefficient sums" begin
    p1 = TrotterParams(0.1, 1, :real, 2, 1e-12)
    p2 = TrotterParams(0.1, 2, :real, 2, 1e-12)
    seq1 = trotter_sequence(p1)
    seq2 = trotter_sequence(p2)

    @test sum(step for (_, step) in seq1) ≈ 5 * p1.dt
    @test sum(step for (_, step) in seq2) ≈ 5 * p2.dt
    @test [color for (color, _) in seq1] == collect(1:5)
    @test [color for (color, _) in seq2] == [1, 2, 3, 4, 5, 4, 3, 2, 1]
    @test [sum(step for (color, step) in seq1 if color == c) for c = 1:5] ≈
          fill(p1.dt, 5)
    @test [sum(step for (color, step) in seq2 if color == c) for c = 1:5] ≈
          fill(p2.dt, 5)
end

@testset "five-color schedule has disjoint finite periodic stars" begin
    cell = PeriodicSquareUnitCell(10, 10)
    covered_centers = Set{Tuple{Int,Int}}()

    for color in 1:5
        centers = update_centers(cell, color)
        @test stars_are_disjoint_mod_unitcell(cell, centers)
        occupied = Set{Tuple{Int,Int}}()

        for c in centers
            push!(covered_centers, (wrap(cell, c).x, wrap(cell, c).y))
            for s in square_star_sites(c)
                wrapped = wrap(cell, s)
                key = (wrapped.x, wrapped.y)
                @test !(key in occupied)
                push!(occupied, key)
            end
        end
    end

    @test length(covered_centers) == length(cell.reps)
    @test covered_centers == Set((c.x, c.y) for c in cell.reps)
end

@testset "finite reference center mapping" begin
    L = 5
    cell = PeriodicSquareUnitCell(L, L)
    expected_centers = Dict(
        1 => [3, 6, 14, 17, 25],
        2 => [4, 7, 15, 18, 21],
        3 => [5, 8, 11, 19, 22],
        4 => [1, 9, 12, 20, 23],
        5 => [2, 10, 13, 16, 24],
    )

    for color = 1:5
        centers = update_centers(cell, color)
        @test [
            _finite_site_index(wrap(cell, c).x, wrap(cell, c).y, L) for
            c in centers
        ] == expected_centers[color]
    end

    center = (3, 3)
    @test _finite_star_sites(center, L) == (13, 14, 18, 12, 8)
    @test _finite_star_sites((1, 1), L) == (1, 2, 6, 5, 21)
    @test _production_finite_star_sites(center, L) == _finite_star_sites(center, L)
    @test _production_finite_star_sites((1, 1), L) == _finite_star_sites((1, 1), L)
end
