@testset "scar finder" begin
    @testset "config validation" begin
        cfg = ScarFinderConfig(0.01, 1, 2, 2, 1e-12, OneSiteUnitCell(), 3, 1.0)
        @test cfg.seed_count == 3
        @test_throws ArgumentError ScarFinderConfig(0.0, 1, 2, 2, 1e-12, OneSiteUnitCell(), 3, 1.0)
        @test_throws ArgumentError ScarFinderConfig(0.01, 0, 2, 2, 1e-12, OneSiteUnitCell(), 3, 1.0)
        @test_throws ArgumentError ScarFinderConfig(0.01, 1, -1, 2, 1e-12, OneSiteUnitCell(), 3, 1.0)
        @test_throws ArgumentError ScarFinderConfig(0.01, 1, 2, 0, 1e-12, OneSiteUnitCell(), 3, 1.0)
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

    @testset "search rejects unsupported fixed-D greater than one" begin
        cfg = ScarFinderConfig(0.005, 1, 1, 2, 1e-12, OneSiteUnitCell(), 2, 1.0)
        @test_throws ArgumentError scar_search(cfg; seed = 123)
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
