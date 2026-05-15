using Aqua

@testset "package quality" begin
    Aqua.test_all(
        SquarePXPDynamics;
        ambiguities = false,
        # JSON3 is added before the benchmark writer that imports it.
        stale_deps = (ignore = [:JSON3],),
    )
end
