using Aqua

@testset "package quality" begin
    Aqua.test_all(SquarePXPDynamics; ambiguities = false)
end
