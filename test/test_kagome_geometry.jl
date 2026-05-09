@testset "kagome geometry" begin
    @testset "coord and sublattice" begin
        c = KagomeCoord(0, 0, :A)
        @test c.n1 == 0 && c.n2 == 0 && c.sublat === :A
        @test_throws ArgumentError KagomeCoord(0, 0, :X)  # invalid sublattice
    end

    @testset "coordination 4: every site has exactly 4 neighbors" begin
        for sublat in (:A, :B, :C)
            c = KagomeCoord(0, 0, sublat)
            nbrs = [kagome_neighbor(c, d) for d in 1:4]
            @test length(unique(nbrs)) == 4
        end
    end

    @testset "two-triangles-per-site invariant" begin
        # Each kagome site belongs to one up-triangle and one down-triangle.
        for sublat in (:A, :B, :C)
            c = KagomeCoord(1, 2, sublat)
            up = up_triangle_of(c)
            down = down_triangle_of(c)
            @test up !== down
            @test c in triangle_sites(up)
            @test c in triangle_sites(down)
        end
    end

    @testset "5-site star occupies 5 distinct positions" begin
        c = KagomeCoord(0, 0, :A)
        star = kagome_star_sites(c)
        @test length(star) == 5
        @test length(unique(star)) == 5
        @test c == star[1]   # center first
    end

    @testset "9-site UC: 5-star occupies 5 distinct reps" begin
        uc = NineSiteKagomeUC()
        reps = unit_cell_representatives(uc)
        @test length(reps) == 9
        for c in reps
            star = kagome_star_sites(c)
            wrapped = [wrap_kagome_coord(uc, sc) for sc in star]
            @test length(unique(wrapped)) == 5
        end
    end

    @testset "3-coloring: same-color stars are vertex-disjoint" begin
        # Stars centered on the same sublattice (A, B, or C) form one color.
        # Disjointness is checked over a small window of centers.
        centers = [KagomeCoord(n1, n2, s) for n1 in -2:2 for n2 in -2:2 for s in (:A, :B, :C)]
        for color in 1:3
            same = filter(c -> kagome_star_color(c) == color, centers)
            for i in eachindex(same), j in (i+1):lastindex(same)
                @test disjoint_kagome_stars(same[i], same[j])
            end
        end
    end
end
