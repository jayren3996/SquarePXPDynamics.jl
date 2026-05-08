using Test
using TriangularPEPSDynamics

@testset "TriangularPEPSDynamics" begin
    include("test_geometry.jl")
    include("test_spinops.jl")
    include("test_models.jl")
    include("test_gates.jl")
    include("test_schedules.jl")
    include("test_solvable_models.jl")
end
