using Test
using SquarePXPDynamics

@testset "SquarePXPDynamics" begin
    include("test_spinops.jl")
    include("test_square_geometry.jl")
    include("test_square_pxp.jl")
    include("test_square_peps.jl")
    include("test_square_unitcells.jl")
    include("test_square_ipeps.jl")
    include("test_square_ipeps_s2.jl")
    include("test_star_simple_update.jl")
    include("test_ipeps_evolution.jl")
    include("test_observables_evolved.jl")
    include("test_pepskit_measurements.jl")
    include("test_scarfinder.jl")
    include("test_public_docs.jl")
    include("test_aqua.jl")
end
