using LinearAlgebra
using ITensors

@testset "simple update" begin
    @testset "product-gate factorization preserves kron order for asymmetric factors" begin
        factors = [
            ComplexF64[1.0 0.2; 0.3 1.1],
            ComplexF64[0.7 -0.1im; 0.4 0.9],
            ComplexF64[1.2 0.0; 0.1im 0.8],
            ComplexF64[0.6 0.3; -0.2 1.0],
            ComplexF64[1.0 0.5im; 0.25 0.75],
            ComplexF64[0.9 -0.15; 0.2im 1.05],
            ComplexF64[1.1 0.05; -0.35im 0.65],
        ]
        G = factors[1]
        for f in factors[2:end]
            G = kron(G, f)
        end

        recovered = TriangularPEPSDynamics.SimpleUpdate._try_factorize_product_gate(G, 7)

        @test recovered !== nothing
        if recovered !== nothing
            reconstructed = recovered[1]
            for f in recovered[2:end]
                reconstructed = kron(reconstructed, f)
            end
            @test reconstructed ≈ G
        end
    end

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

    @testset "global-phase identity gate is a no-op" begin
        state = random_ipeps(OneSiteUnitCell(), 2; seed = 17)
        T_before = copy(site_tensor(state, Coord(0, 0)))

        G = cis(0.37) * Matrix{ComplexF64}(I, 128, 128)
        diag = apply_star_gate_simple_update!(state, G, Coord(0, 0); maxdim = 2)

        @test diag isa SimpleUpdateDiagnostics
        @test diag.discarded_weight == 0
        @test array(T_before) ≈ array(site_tensor(state, Coord(0, 0)))
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

    @testset "projected PXP gate updates D=1 three-site product exactly" begin
        t = 0.17
        H = pxp_star_hamiltonian(projector_down(), pauli_x())
        Uproj = projected_gate(H, t; evolution = :real)
        state = product_ipeps(ThreeSiteUnitCell(), :down; D = 1)

        diag = apply_star_gate_simple_update!(state, Uproj, Coord(0, 0); maxdim = 1)

        @test diag isa SimpleUpdateDiagnostics
        @test diag.discarded_weight ≈ 0 atol = 1e-12
        @test real(local_expectation(state, Coord(0, 0), pauli_z())) ≈ -cos(2t) atol = 1e-10
        @test real(local_expectation(state, Coord(1, 0), pauli_z())) ≈ -1 atol = 1e-10
        @test real(local_expectation(state, Coord(2, 0), pauli_z())) ≈ -1 atol = 1e-10
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

    @testset "simple update diagnostics report affected bonds and dimensions" begin
        state = random_ipeps(OneSiteUnitCell(), 2; seed = 7)
        I128 = Matrix{ComplexF64}(I, 128, 128)
        diag = apply_star_gate_simple_update!(state, I128, Coord(0, 0); maxdim = 2)
        @test Set(diag.affected_bonds) == Set((Coord(0, 0), d) for d in 1:6)
        @test diag.output_bond_dims == fill(2, 6)
        for d in 1:6
            λ = bond_lambda(state, Coord(0, 0), d)
            @test all(λ .>= 0)
            @test norm(λ) ≈ sqrt(length(λ))
        end
    end

    @testset "general non-product star updates at D>1 use local product projection" begin
        state = random_ipeps(OneSiteUnitCell(), 2; seed = 11)
        T_before = copy(site_tensor(state, Coord(0, 0)))
        H = pxp_star_hamiltonian(projector_down(), pauli_x())
        Uproj = projected_gate(H, 0.02; evolution = :real)

        diag = apply_star_gate_simple_update!(
            state, Uproj, Coord(0, 0); maxdim = 2, cutoff = 1e-12,
        )
        @test Set(diag.affected_bonds) == Set((Coord(0, 0), d) for d in 1:6)
        @test diag.output_bond_dims == fill(2, 6)
        @test isfinite(diag.discarded_weight)
        @test 0 <= diag.discarded_weight <= 1
        @test !(array(T_before) ≈ array(site_tensor(state, Coord(0, 0))))
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

    @testset "general star kernel: random unitary with maxdim large enough is exact" begin
        using Random
        rng = MersenneTwister(2026_05_09)

        state = random_ipeps(ThreeSiteUnitCell(), 2; seed = 11)
        psi_before = cluster_vector_from_state(state, Coord(0, 0))
        @test isapprox(norm(psi_before), norm(psi_before); atol = 0)  # sanity

        # Haar-random unitary on the 128-dim cluster Hilbert space.
        A = randn(rng, ComplexF64, 128, 128)
        Q, _ = qr(A)
        G = Matrix{ComplexF64}(Q)
        @test G' * G ≈ Matrix{ComplexF64}(I, 128, 128) atol = 1e-10

        # Apply via the general path. maxdim large enough that bond growth from
        # 2 -> up to 2 * 2 = 4 (per peel) is uncapped.
        diag = apply_star_gate_simple_update!(state, G, Coord(0, 0); maxdim = 64, cutoff = 0.0)

        psi_after = cluster_vector_from_state(state, Coord(0, 0))
        expected = G * psi_before

        # Compare up to a global phase / norm because the iPEPS may absorb a
        # state-wide scalar that cancels in any expectation value.
        overlap = abs(dot(psi_after, expected))
        @test overlap ≈ norm(psi_after) * norm(expected) atol = 1e-8
        @test diag isa SimpleUpdateDiagnostics
        @test diag.discarded_weight ≈ 0 atol = 1e-10
    end
end
