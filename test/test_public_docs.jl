@testset "public docstrings" begin
    undocumented = [
        name for name in names(SquarePXPDynamics; all = false, imported = false) if
        !Docs.hasdoc(SquarePXPDynamics, name)
    ]

    @test isempty(undocumented)
    @test Docs.hasdoc(SquarePXPDynamics, :StarUpdateInfo)
    @test Docs.hasdoc(SquarePXPDynamics, :project_star!)
end

@testset "public exports do not expose internal star helpers" begin
    public_names = names(SquarePXPDynamics; all = false)
    @test all(name -> !startswith(String(name), "_"), public_names)
    @test !(:_validate_split_order in public_names)
    @test !(:_absorb_star_weights in public_names)
    @test !(:_split_reduced_theta in public_names)
    @test !(:_commit_star_update! in public_names)
end
