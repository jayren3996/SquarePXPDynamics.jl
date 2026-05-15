using Test
using SquarePXPDynamics

@testset "SquarePXPDynamics" begin
    include("test_spinops.jl")
    include("test_square_geometry.jl")
    include("test_square_pxp.jl")
    include("test_square_peps.jl")
    include("test_public_docs.jl")
    include("test_aqua.jl")
end
