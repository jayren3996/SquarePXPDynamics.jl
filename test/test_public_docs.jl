@testset "public docstrings" begin
    undocumented = [
        name for name in names(SquarePXPDynamics; all = false, imported = false) if
        !Docs.hasdoc(SquarePXPDynamics, name)
    ]

    @test isempty(undocumented)
    @test Docs.hasdoc(SquarePXPDynamics, :StarUpdateInfo)
    @test Docs.hasdoc(SquarePXPDynamics, :project_star!)
end

@testset "load-bearing API docstrings include an Example block" begin
    # A curated set of public functions whose docstrings should contain a
    # runnable example. Adding a name here is a deliberate decision that the
    # function is load-bearing enough to deserve worked usage.
    @test occursin("# Example", string(@doc SquarePXPDynamics.evolve!))
    @test occursin("# Example", string(@doc SquarePXPDynamics.project_star!))
    @test occursin("# Example", string(@doc SquarePXPDynamics.measure_simple))
    @test occursin("# Example", string(@doc SquarePXPDynamics.measure_ctm))
    @test occursin("# Example", string(@doc SquarePXPDynamics.scarfinder!))
    @test occursin(
        "# Example",
        string(@doc SquarePXPDynamics.validate_pxp_ed_ipeps),
    )
end

@testset "public exports do not expose internal star helpers" begin
    public_names = names(SquarePXPDynamics; all = false)
    @test all(name -> !startswith(String(name), "_"), public_names)
    @test !(:_validate_split_order in public_names)
    @test !(:_absorb_star_weights in public_names)
    @test !(:_split_reduced_theta in public_names)
    @test !(:_commit_star_update! in public_names)
end
