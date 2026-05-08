using LinearAlgebra
using ITensors

@testset "simple update" begin
    @testset "identity gate: no-op on D=1 product" begin
        state = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        T_before = copy(site_tensor(state, Coord(0, 0)))

        I128 = Matrix{ComplexF64}(I, 128, 128)
        diag = apply_star_gate_simple_update!(state, I128, Coord(0, 0))

        T_after = site_tensor(state, Coord(0, 0))
        # Local Z is unchanged
        @test local_expectation(state, Coord(0, 0), pauli_z()) ≈ -1
        @test diag isa SimpleUpdateDiagnostics
        @test diag.discarded_weight == 0
        @test array(T_before) ≈ array(T_after)
    end

    @testset "site-symmetric single-site gate (u^⊗7) on D=1" begin
        # u = exp(-i theta X / 2): rotation about X by 2*theta.
        theta = 0.37
        u = cos(theta) * Matrix{ComplexF64}(I, 2, 2) - im * sin(theta) * pauli_x()
        # G = u ⊗ u ⊗ ... ⊗ u (7 factors)
        G = u
        for _ in 2:7
            G = kron(G, u)
        end

        state = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        diag = apply_star_gate_simple_update!(state, G, Coord(0, 0))
        @test diag isa SimpleUpdateDiagnostics

        # Dense expectation: <down| u' Z u |down> = -cos(2 theta)
        @test local_expectation(state, Coord(0, 0), pauli_z()) ≈ -cos(2 * theta) atol = 1e-10
        # X expectation: <down| u' X u |down>
        v_after = u * ComplexF64[0, 1]
        expected_X = real(v_after' * pauli_x() * v_after)
        @test real(local_expectation(state, Coord(0, 0), pauli_x())) ≈ expected_X atol = 1e-10
    end

    @testset "non-product gate raises informative error at D=1" begin
        # PXP gate at finite t: not a tensor product across the 7 sites
        H = pxp_star_hamiltonian(projector_down(), pauli_x())
        Upxp = dense_gate(H, 0.1; evolution = :real)
        state = product_ipeps(OneSiteUnitCell(), :down; D = 1)
        @test_throws ErrorException apply_star_gate_simple_update!(state, Upxp, Coord(0, 0))
    end

    @testset "D=2 random state: identity gate is a no-op" begin
        state = random_ipeps(OneSiteUnitCell(), 2; seed = 7)
        T_before = copy(site_tensor(state, Coord(0, 0)))
        I128 = Matrix{ComplexF64}(I, 128, 128)
        diag = apply_star_gate_simple_update!(state, I128, Coord(0, 0))
        T_after = site_tensor(state, Coord(0, 0))
        @test array(T_before) ≈ array(T_after)
        @test diag.discarded_weight == 0
    end

    @testset "PXP-on-allowed-product preserves blockade (dense level)" begin
        # All-down product is blockade-allowed. The projected PXP gate must keep it allowed.
        all_down = zeros(ComplexF64, 128)
        all_down[end] = 1
        @test dense_blockade_violations(all_down) ≈ 0

        H = pxp_star_hamiltonian(projector_down(), pauli_x())
        for t in (0.0, 0.05, 0.3, 1.1)
            Uproj = projected_gate(H, t; evolution = :real)
            evolved = Uproj * all_down
            @test dense_blockade_violations(evolved) < 1e-12
        end
    end
end
