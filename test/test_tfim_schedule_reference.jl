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

@testset "trotter coefficient sums" begin
    p1 = TrotterParams(0.1, 1, :real, 2, 1e-12)
    p2 = TrotterParams(0.1, 2, :real, 2, 1e-12)

    @test sum(step for (_, step) in trotter_sequence(p1)) ≈ 5 * p1.dt
    @test [color for (color, _) in trotter_sequence(p1)] == collect(1:5)
    @test [color for (color, _) in trotter_sequence(p2)] == [1, 2, 3, 4, 5, 4, 3, 2, 1]
    @test sum(step for (color, step) in trotter_sequence(p2) if color == 1) ≈ p2.dt
    @test sum(step for (color, step) in trotter_sequence(p2) if color == 5) ≈ p2.dt
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
    center = (3, 3)

    @test _finite_star_sites(center, L) == (13, 14, 18, 12, 8)
    @test _finite_star_sites((1, 1), L) == (1, 2, 6, 5, 21)
end
