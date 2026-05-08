using ITensors

@testset "states" begin
    uc1 = OneSiteUnitCell()
    @test unit_cell_representatives(uc1) == (Coord(0, 0),)
    @test wrap_coord(uc1, Coord(5, -3)) == Coord(0, 0)

    uc3 = ThreeSiteUnitCell()
    @test length(unit_cell_representatives(uc3)) == 3
    reps3 = unit_cell_representatives(uc3)
    @test wrap_coord(uc3, Coord(0, 0)) in reps3
    # Sublattice partitioning: any two reps must lie on different sublattices.
    @test length(unique(wrap_coord(uc3, c) for c in reps3)) == 3
    # Going around the 6 neighbors of any rep visits the other two sublattices in alternating fashion.
    for c in reps3
        wrapped_neighbors = [wrap_coord(uc3, neighbor(c, d)) for d in 1:6]
        @test length(unique(wrapped_neighbors)) == 2
        @test !(c in wrapped_neighbors)
    end

    @testset "product :down at D=1 (1-site unit cell)" begin
        state = product_ipeps(uc1, :down; D = 1)
        @test state isa TriangularIPEPS

        c0 = Coord(0, 0)
        T = site_tensor(state, c0)
        ph = phys_index(state, c0)
        @test dim(ph) == 2
        # All bond indices have dim 1
        for d in 1:6
            @test dim(bond_index(state, c0, d)) == 1
        end
        # Local physical vector: |down> at index 2
        v_down = ITensor(ComplexF64[0, 1], ph)
        # Project onto v_down by contracting; remaining is bond-only (all dim 1)
        proj = T * dag(v_down)
        @test abs(scalar(proj) - 1) < 1e-12

        v_up = ITensor(ComplexF64[1, 0], ph)
        proj_up = T * dag(v_up)
        @test abs(scalar(proj_up)) < 1e-12
    end

    @testset "product :up at D=1" begin
        state = product_ipeps(uc1, :up; D = 1)
        c0 = Coord(0, 0)
        T = site_tensor(state, c0)
        ph = phys_index(state, c0)
        v_up = ITensor(ComplexF64[1, 0], ph)
        @test abs(scalar(T * dag(v_up)) - 1) < 1e-12
    end

    @testset "opposite virtual bonds match" begin
        for uc in (OneSiteUnitCell(), ThreeSiteUnitCell())
            state = product_ipeps(uc, :down; D = 1)
            for c in unit_cell_representatives(uc)
                for d in 1:6
                    opp_d = opposite_direction(d)
                    nbr_rep = wrap_coord(uc, neighbor(c, d))
                    b_here = bond_index(state, c, d)
                    b_there = bond_index(state, nbr_rep, opp_d)
                    # Dimensions of opposite bonds always match.
                    @test dim(b_here) == dim(b_there)
                    # Lambda spectra are shared across the bond.
                    @test bond_lambda(state, c, d) === bond_lambda(state, nbr_rep, opp_d)
                    if nbr_rep != c
                        # Distinct-rep neighbors share the same Index object.
                        @test b_here === b_there
                    end
                end
            end
        end
    end

    @testset "random_ipeps: shape and reproducibility" begin
        s1 = random_ipeps(OneSiteUnitCell(), 2; seed = 42)
        s2 = random_ipeps(OneSiteUnitCell(), 2; seed = 42)
        c0 = Coord(0, 0)
        T1 = site_tensor(s1, c0)
        T2 = site_tensor(s2, c0)
        @test dim(phys_index(s1, c0)) == 2
        @test dim(bond_index(s1, c0, 1)) == 2
        # Same seed → same tensor data
        @test array(T1) ≈ array(T2)

        s3 = random_ipeps(OneSiteUnitCell(), 2; seed = 43)
        T3 = site_tensor(s3, c0)
        @test !(array(T1) ≈ array(T3))
    end

    @testset "3-site product state has all reps in :down" begin
        state = product_ipeps(ThreeSiteUnitCell(), :down; D = 1)
        for c in unit_cell_representatives(ThreeSiteUnitCell())
            ph = phys_index(state, c)
            v_down = ITensor(ComplexF64[0, 1], ph)
            T = site_tensor(state, c)
            @test abs(scalar(T * dag(v_down)) - 1) < 1e-12
        end
    end
end
