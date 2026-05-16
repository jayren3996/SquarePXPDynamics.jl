using Aqua

@testset "package quality" begin
    Aqua.test_all(SquarePXPDynamics)
end
