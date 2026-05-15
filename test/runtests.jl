using Test
using SquarePXPDynamics

const TEST_FILES = [
    "test_spinops.jl",
    "test_square_geometry.jl",
    "test_square_pxp.jl",
    "test_star_models.jl",
    "test_square_peps.jl",
    "test_square_unitcells.jl",
    "test_square_ipeps.jl",
    "test_square_ipeps_s2.jl",
    "test_star_simple_update.jl",
    "test_ipeps_evolution.jl",
    "test_observables_evolved.jl",
    "test_pepskit_measurements.jl",
    "test_scarfinder.jl",
    "test_public_docs.jl",
    "test_aqua.jl",
]

@testset "SquarePXPDynamics" begin
    files = isempty(ARGS) ? TEST_FILES : ARGS
    for file in files
        file in TEST_FILES || throw(ArgumentError("unknown test file: $file"))
        include(file)
    end
end
