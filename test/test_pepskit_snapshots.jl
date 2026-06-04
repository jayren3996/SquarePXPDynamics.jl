# PEPSKit-native candidate snapshots: write/load round trip for PXPIPEPSState.
#
# Native increment-2 feature: persist the InfinitePEPS + SUWeight state so
# ScarFinder candidates can be reloaded on the :pepskit_simple engine, the
# native counterpart of write_square_ipeps_snapshot. The round trip must
# reproduce the state exactly so all observables and subsequent evolution match.
using Test
using JLD2
using SquarePXPDynamics
using SquarePXPDynamics: pepskit_peps, pepskit_weights

# Compare two native states tensor-for-tensor and weight-for-weight.
function _states_bitmatch(a, b; atol = 1e-12)
    pa, pb = pepskit_peps(a), pepskit_peps(b)
    wa, wb = pepskit_weights(a), pepskit_weights(b)
    cell = a.unitcell
    Nr, Nc = cell.Ly, cell.Lx
    for row in 1:Nr, col in 1:Nc
        isapprox(convert(Array, pa[row, col]), convert(Array, pb[row, col]); atol) ||
            return false
    end
    for slot in 1:2, row in 1:Nr, col in 1:Nc
        isapprox(wa.data[slot, row, col].data, wb.data[slot, row, col].data; atol) ||
            return false
    end
    return true
end

_max_bond(state) = maximum(length(w.data) for w in pepskit_weights(state).data)

@testset "PEPSKit-native snapshots (PXPIPEPSState)" begin
    @testset "D=1 round trip" begin
        cell = PeriodicSquareUnitCell(4, 4)
        st = product_pxp_pepskit_state(cell; state = :x_plus, D = 1)
        path = tempname() * ".jld2"
        @test write_pxp_pepskit_snapshot(path, st) == path

        loaded = load_pxp_pepskit_snapshot(path)
        @test loaded isa PXPIPEPSState
        @test loaded.unitcell.Lx == 4 && loaded.unitcell.Ly == 4
        @test state_version(loaded) == state_version(st)
        @test log_norm(loaded) ≈ log_norm(st) atol = 1e-12
        @test _states_bitmatch(st, loaded)

        a = measure_simple(st)
        b = measure_simple(loaded)
        @test b.density ≈ a.density atol = 1e-12
        @test b.blockade_violation ≈ a.blockade_violation atol = 1e-12
        @test b.pxp_energy_density ≈ a.pxp_energy_density atol = 1e-12
    end

    @testset "D>1 round trip + evolution parity" begin
        cell = PeriodicSquareUnitCell(4, 4)
        st = product_pxp_pepskit_state(cell; state = :down, D = 1)
        # Evolve to genuine D>1 with small maxdim (measure_simple's dense energy
        # patch is only feasible at D<=2-3).
        evolve_pepskit!(st, 0.02; dt = 0.01, order = 1, schedule = :serial,
                        maxdim = 2, rel_floor = 0.0)
        @test _max_bond(st) >= 2          # genuinely entangled state on disk

        path = tempname() * ".jld2"
        write_pxp_pepskit_snapshot(path, st)
        loaded = load_pxp_pepskit_snapshot(path)
        @test _states_bitmatch(st, loaded; atol = 1e-12)
        @test _max_bond(loaded) == _max_bond(st)
        @test state_version(loaded) == state_version(st)
        @test log_norm(loaded) ≈ log_norm(st) atol = 1e-12

        # Observables agree on the reloaded state.
        a = measure_simple(st)
        b = measure_simple(loaded)
        @test b.density ≈ a.density atol = 1e-10
        @test b.max_bond_entropy ≈ a.max_bond_entropy atol = 1e-10

        # Evolving the original and the reloaded state one more step gives the
        # same state: the snapshot preserves all dynamical content.
        st2 = load_pxp_pepskit_snapshot(path)
        evolve_pepskit!(st, 0.01; dt = 0.01, order = 1, schedule = :serial,
                        maxdim = 2, rel_floor = 0.0)
        evolve_pepskit!(st2, 0.01; dt = 0.01, order = 1, schedule = :serial,
                        maxdim = 2, rel_floor = 0.0)
        @test _states_bitmatch(st, st2; atol = 1e-10)
    end

    @testset "format version guard" begin
        # The version check runs before any other key is read, so a stub file
        # with only an unsupported format_version triggers the guard.
        bad = tempname() * ".jld2"
        JLD2.jldsave(bad; format_version = PXP_PEPSKIT_SNAPSHOT_FORMAT_VERSION + 1)
        @test_throws ArgumentError load_pxp_pepskit_snapshot(bad)
    end
end
