using Test
using KagomePXPDynamics

@testset "KagomePXPDynamics" begin
    include("test_spinops.jl")
    include("test_solvable_models.jl")
    include("test_kagome_geometry.jl")
end
