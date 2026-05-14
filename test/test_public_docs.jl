@testset "public docstrings" begin
    undocumented = [
        name for name in names(SquarePXPDynamics; all = false, imported = false)
        if !Docs.hasdoc(SquarePXPDynamics, name)
    ]

    @test isempty(undocumented)
end
