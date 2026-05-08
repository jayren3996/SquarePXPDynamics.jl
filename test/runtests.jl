using Test
using TriangularPEPSDynamics

@testset "TriangularPEPSDynamics" begin
    include("test_geometry.jl")
    include("test_spinops.jl")
    include("test_models.jl")
    include("test_gates.jl")
    include("test_schedules.jl")
    include("test_solvable_models.jl")
    include("test_states.jl")
    include("test_observables.jl")
    include("test_simple_update.jl")
    include("test_evolution.jl")
    include("test_scar_finder.jl")
end
