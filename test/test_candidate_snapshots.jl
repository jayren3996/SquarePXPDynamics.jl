using Test
using JLD2
using SquarePXPDynamics

function _approx_equal_observables(a, b; atol = 1e-10)
    return isapprox(a.density, b.density; atol) &&
           isapprox(a.density_even, b.density_even; atol) &&
           isapprox(a.density_odd, b.density_odd; atol) &&
           isapprox(a.blockade_violation, b.blockade_violation; atol) &&
           isapprox(a.pxp_energy_density, b.pxp_energy_density; atol)
end

@testset "SquareIPEPSState snapshot round-trips a product state" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    before = measure_simple(psi)

    path = tempname() * ".jld2"
    try
        write_square_ipeps_snapshot(path, psi)
        @test isfile(path)

        restored = load_square_ipeps_snapshot(path)
        @test restored.unitcell.Lx == cell.Lx
        @test restored.unitcell.Ly == cell.Ly
        @test restored.maxdim == psi.maxdim
        @test restored.gauge === psi.gauge

        after = measure_simple(restored)
        @test _approx_equal_observables(before, after)
    finally
        rm(path; force = true)
    end
end

@testset "SquareIPEPSState snapshot round-trips an evolved D=2 state" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    evolve!(psi, 0.04; params = TrotterParams(0.02, 1, :real, 2, 1e-12; schedule = :serial))
    before = measure_simple(psi)
    log_norm_before = log_norm(psi)
    version_before = state_version(psi)

    path = tempname() * ".jld2"
    try
        write_square_ipeps_snapshot(path, psi)
        restored = load_square_ipeps_snapshot(path)

        @test restored.unitcell.Lx == cell.Lx
        @test restored.unitcell.Ly == cell.Ly
        @test restored.maxdim == psi.maxdim
        @test log_norm(restored) ≈ log_norm_before atol = 1e-12
        @test state_version(restored) == version_before

        after = measure_simple(restored)
        @test _approx_equal_observables(before, after; atol = 1e-10)
    finally
        rm(path; force = true)
    end
end

@testset "SquareIPEPSState snapshot rejects unknown format_version" begin
    path = tempname() * ".jld2"
    try
        JLD2.jldsave(
            path;
            format_version = 0,
            Lx = 3,
            Ly = 3,
            maxdim = 1,
            gauge = :simple,
            log_norm_value = 0.0,
            mutation_version = 0,
            site_arrays = Dict{Tuple{Int,Int},Array{ComplexF64,5}}(),
            link_weights = Dict{Tuple{Int,Int,Symbol},Vector{Float64}}(),
            link_dims = Dict{Tuple{Int,Int,Symbol},Int}(),
        )
        @test_throws ArgumentError load_square_ipeps_snapshot(path)
    finally
        rm(path; force = true)
    end
end

@testset "JLD2CandidateStore writes both metadata and snapshot" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    evolve!(psi, 0.02; params = TrotterParams(0.02, 1, :real, 2, 1e-12; schedule = :serial))

    directory = mktempdir()
    try
        store = JLD2CandidateStore(directory)
        trotter = TrotterParams(0.02, 1, :real, 2, 1e-12; schedule = :serial)
        params = ScarFinderParams(0.02, trotter, 1, Inf, Inf, Inf, false)
        result = scarfinder!(psi, params; candidate_store = store)
        @test result.accepted_iterations + result.rejected_iterations == 1

        metadata_files = filter(endswith(".json"), readdir(directory))
        snapshot_files = filter(endswith(".jld2"), readdir(directory))
        @test length(metadata_files) == 1
        @test length(snapshot_files) == 1

        snapshot_path = joinpath(directory, snapshot_files[1])
        restored = load_square_ipeps_snapshot(snapshot_path)
        @test restored.unitcell.Lx == cell.Lx
        @test restored.unitcell.Ly == cell.Ly

        original = measure_simple(psi)
        restored_obs = measure_simple(restored)
        @test _approx_equal_observables(original, restored_obs)
    finally
        rm(directory; recursive = true, force = true)
    end
end

