using ITensors

@testset "scar finder" begin
    @testset "config validation" begin
        cfg = ScarFinderConfig(0.01, 1, 2, 2, 1e-12, OneSiteUnitCell(), 3, 1.0)
        @test cfg.seed_count == 3
        @test cfg.dynamics_maxdim == 2
        @test cfg.scar_maxdim == 2
        cfg_split = ScarFinderConfig(0.01, 1, 2, 3, 1e-12, OneSiteUnitCell(), 3, 1.0; scar_maxdim = 2)
        @test cfg_split.dynamics_maxdim == 3
        @test cfg_split.scar_maxdim == 2
        @test_throws ArgumentError ScarFinderConfig(0.0, 1, 2, 2, 1e-12, OneSiteUnitCell(), 3, 1.0)
        @test_throws ArgumentError ScarFinderConfig(0.01, 0, 2, 2, 1e-12, OneSiteUnitCell(), 3, 1.0)
        @test_throws ArgumentError ScarFinderConfig(0.01, 1, -1, 2, 1e-12, OneSiteUnitCell(), 3, 1.0)
        @test_throws ArgumentError ScarFinderConfig(0.01, 1, 2, 0, 1e-12, OneSiteUnitCell(), 3, 1.0)
        @test_throws ArgumentError ScarFinderConfig(0.01, 1, 2, 2, 1e-12, OneSiteUnitCell(), 3, 1.0; scar_maxdim = 3)
        @test_throws ArgumentError ScarFinderConfig(0.01, 1, 2, 2, -1e-12, OneSiteUnitCell(), 3, 1.0)
        @test_throws ArgumentError ScarFinderConfig(0.01, 1, 2, 2, 1e-12, OneSiteUnitCell(), 0, 1.0)
        @test_throws ArgumentError ScarFinderConfig(0.01, 1, 2, 2, 1e-12, OneSiteUnitCell(), 3, -0.1)
    end

    @testset "search is deterministic and returns diagnostics" begin
        cfg = ScarFinderConfig(0.005, 1, 2, 1, 1e-12, OneSiteUnitCell(), 1, 1.0)
        a = scar_search(cfg; seed = 123)
        b = scar_search(cfg; seed = 123)

        @test length(a) == cfg.seed_count
        @test [c.score for c in a] ≈ [c.score for c in b]
        @test all(!isempty(c.diagnostics) for c in a)
        @test all(isfinite(c.score) for c in a)
        @test all(c.entanglement_proxy >= 0 for c in a)
    end

    @testset "search runs D=2 candidates across supported unit cells" begin
        for uc in (OneSiteUnitCell(), ThreeSiteUnitCell(), SevenSiteUnitCell())
            cfg = ScarFinderConfig(0.005, 1, 2, 2, 1e-12, uc, 4, 1.0)
            candidates = scar_search(cfg; seed = 321)
            @test length(candidates) == 4
            @test all(!isempty(c.diagnostics) for c in candidates)
            @test all(c.diagnostics[end].max_bond_dim <= 2 for c in candidates)
            @test all(isfinite(c.blockade_violation) for c in candidates)
            @test issorted([c.score for c in candidates])
        end
    end

    @testset "search truncates from dynamics dimension to scar dimension" begin
        cfg = ScarFinderConfig(0.005, 1, 1, 3, 1e-12, ThreeSiteUnitCell(), 2, 1.0; scar_maxdim = 2)
        candidates = scar_search(cfg; seed = 11)
        @test length(candidates) == 2
        @test maximum([
            dim(bond_index(c.state, rep, d))
            for c in candidates
            for rep in unit_cell_representatives(c.state.unitcell)
            for d in 1:6
        ]) <= 2
    end

    @testset "blockade tolerance flags candidates and ranking is stable" begin
        cfg = ScarFinderConfig(0.005, 1, 1, 1, 1e-12, OneSiteUnitCell(), 3, 0.0)
        candidates = scar_search(cfg; seed = 9)
        @test any(!c.accepted for c in candidates)

        ranked1 = rank_candidates(candidates)
        ranked2 = rank_candidates(reverse(candidates))
        @test [c.seed_index for c in ranked1] == [c.seed_index for c in ranked2]
        @test issorted([c.score for c in ranked1])
    end
end
