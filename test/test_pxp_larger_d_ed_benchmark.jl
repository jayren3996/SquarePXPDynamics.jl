using Test
using JSON3
using SquarePXPDynamics

const RUN_EXTENDED_PXP_ED_TESTS =
    get(ENV, "SQUAREPXP_EXTENDED_PXP_ED_TESTS", "") == "1" ||
    get(ENV, "SQUAREPXP_EXTENDED_TESTS", "") == "1"

function _csv_cell(lines, name)
    header = split(lines[1], ','; keepempty = true)
    row = split(lines[2], ','; keepempty = true)
    index = findfirst(==(name), header)
    index === nothing && error("missing CSV column $name")
    return row[index]
end

@testset "PXP ED observable provenance rejects central-region claims" begin
    basis = pxp_ed_space_group_basis(5)

    @test pxp_ed_boundary_condition(basis) === :periodic
    @test pxp_ed_symmetry_sector(basis) === :fully_symmetric_space_group
    @test pxp_ed_observable_scope(basis) === :pbc_global_site_average
    @test pxp_ed_reference_label(basis) == "finite_pbc_global_density"
    @test pxp_ed_group_order(basis) == 8 * 5^2

    @test_throws ArgumentError pxp_ed_site_density_operator(basis, 13)
    @test_throws ArgumentError pxp_ed_region_density_operator(basis, 1:9)
end

@testset "exact finite return probability is available only through tiny-cell contraction" begin
    cell = PeriodicSquareUnitCell(3, 3)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    up = product_square_ipeps(cell; state = :up, maxdim = 1)

    @test exact_all_down_return_probability_finite(down; max_sites = 9) ≈ 1.0 atol = 1e-15
    @test exact_all_down_return_probability_finite(up; max_sites = 9) ≈ 0.0 atol = 1e-15

    large = product_square_ipeps(PeriodicSquareUnitCell(4, 4); state = :down, maxdim = 1)
    @test_throws ArgumentError exact_all_down_return_probability_finite(large; max_sites = 9)
end
