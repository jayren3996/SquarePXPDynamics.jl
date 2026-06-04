# M6: tensor-engine selector through ScarFinder + audit.
#
# Validates the backend-agnostic simple/evolve interface (native measure_simple
# and evolve! adapters), end-to-end ScarFinder parity between the legacy
# ITensors engine and the PEPSKit-native engine, and the audit tensor_engine
# metadata. Deliberately CTM-free: every check uses the cheap simple-update /
# local observables, so the suite stays fast.
using Test
using SquarePXPDynamics
using SquarePXPDynamics: PeriodicSquareUnitCell, SquareCoord
using SquarePXPDynamics.Observables: measure_simple
using SquarePXPDynamics.IPEPSEvolution: TrotterParams, evolve!
using SquarePXPDynamics.ScarFinder: scarfinder!, ScarFinderParams, SimpleBackend, CompositeObjective

const SPD = SquarePXPDynamics

_obs_fields(s) = (
    s.density, s.density_even, s.density_odd, s.blockade_violation,
    s.pxp_energy_density, s.mean_bond_entropy, s.max_bond_entropy,
)

function _summaries_match(a, b; atol = 1e-8)
    return all(abs(x - y) <= atol for (x, y) in zip(_obs_fields(a), _obs_fields(b)))
end

@testset "M6 tensor-engine selector" begin
    @testset "native measure_simple parity (D=1)" begin
        cell = PeriodicSquareUnitCell(4, 4)
        for (nat, leg) in (
            (product_pxp_pepskit_state(cell; state = :down, D = 1),
             product_square_ipeps(cell; state = :down, maxdim = 1)),
            (product_pxp_pepskit_state(cell; state = :up, D = 1),
             product_square_ipeps(cell; state = :up, maxdim = 1)),
            (product_pxp_pepskit_state(cell; state = :x_plus, D = 1),
             product_square_ipeps(cell; state = :x_plus, maxdim = 1)),
            (checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = 1),
             checkerboard_square_ipeps(cell; excited_on = :even, maxdim = 1)),
        )
            @test _summaries_match(measure_simple(nat), measure_simple(leg); atol = 1e-12)
        end
        # checkerboard density convention: even sublattice excited
        cb = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = 1)
        s = measure_simple(cb)
        @test isapprox(s.density, 0.5; atol = 1e-12)
        @test isapprox(s.density_even, 1.0; atol = 1e-12)
        @test isapprox(s.density_odd, 0.0; atol = 1e-12)
        # :x_plus pins the two- and five-site patch WIRING with nonzero analytic
        # values: ⟨n⟩ = 1/2, ⟨n_i n_j⟩ = 1/4 (blockade), star energy
        # ⟨+|X|+⟩·⟨+|P|+⟩^4 = (1/2)^4 = 1/16.
        sx = measure_simple(product_pxp_pepskit_state(cell; state = :x_plus, D = 1))
        @test isapprox(sx.density, 0.5; atol = 1e-12)
        @test isapprox(sx.blockade_violation, 0.25; atol = 1e-12)
        @test isapprox(sx.pxp_energy_density, 1 / 16; atol = 1e-12)
    end

    @testset "native measure_simple parity (D>1 via converter)" begin
        # Grow a genuinely D>1 legacy state, convert it to a native state in the
        # SAME Vidal gauge, and require the simple observables to agree. This
        # validates the squared-Schmidt cut-leg environment on the D>1 single-,
        # two-, and five-site patches (the patch *wiring* with nonzero analytic
        # values is pinned at D=1 by the :x_plus checks above). The PXP density
        # is nonzero; blockade/energy stay small because PXP suppresses adjacent
        # excitations — the point here is native==legacy, not their magnitude.
        cell = PeriodicSquareUnitCell(4, 4)
        psi = product_square_ipeps(cell; state = :down, maxdim = 2)
        evolve!(psi, 0.5; params = TrotterParams(0.1, 1, :real, 2, 1e-12; schedule = :serial))
        sl = measure_simple(psi)
        @test sl.max_bond_entropy > 1e-6          # genuinely D>1
        @test sl.density > 1e-3                    # nonzero single-site signal
        conv = pxpipeps_from_square_ipeps(psi)
        @test _summaries_match(measure_simple(conv), sl; atol = 1e-8)
    end

    @testset "native evolve! -> EvolutionLog + density parity (D=1)" begin
        # Cross-engine simple-observable parity holds in the gauge-trivial D = 1
        # regime: with maxdim = 1 each star truncates back to a product state, so
        # (per the single-star M4 parity) the native and legacy trajectories stay
        # bit-for-bit identical. At D > 1 the simple observables are gauge-
        # dependent and the two engines pick different gauges; gauge-invariant
        # D > 1 parity (CTM) lives in scripts/dev/verify_m4.jl.
        cell = PeriodicSquareUnitCell(5, 5)
        tp = TrotterParams(0.05, 2, :real, 1, 1e-12; schedule = :five_color)
        nat = product_pxp_pepskit_state(cell; state = :down, D = 1)
        log = evolve!(nat, 0.1; params = tp)
        @test log.nsteps == 2
        @test isfinite(log.max_truncerr) && isfinite(log.max_bond_entropy)
        @test isempty(log.layer_infos)
        leg = product_square_ipeps(cell; state = :down, maxdim = 1)
        evolve!(leg, 0.1; params = tp)
        @test _summaries_match(measure_simple(nat), measure_simple(leg); atol = 1e-8)
    end

    @testset "end-to-end scarfinder! parity (legacy vs native)" begin
        # D = 1 (maxdim = 1) keeps both engines bit-for-bit identical, so the
        # whole orchestration — evolve, measure_simple, accept/reject, score —
        # must agree per iteration. See the note above on the D > 1 gauge.
        cell = PeriodicSquareUnitCell(5, 5)
        trotter = TrotterParams(0.05, 2, :real, 1, 1e-12; schedule = :five_color)
        params = ScarFinderParams(0.1, trotter, 3, Inf, Inf, Inf, false)
        # Score on physical observables only: the raw `max_truncerr` term has an
        # engine-specific accumulation convention, so a truncerr-free objective
        # keeps the score comparison a pure-physics parity check.
        objective = CompositeObjective(; truncation_weight = 0.0, entropy_weight = 0.0,
                                       finite_chi_weight = 0.0)

        nat = product_pxp_pepskit_state(cell; state = :down, D = 1)
        leg = product_square_ipeps(cell; state = :down, maxdim = 1)
        rn = scarfinder!(nat, params; objective, measurement = SimpleBackend())
        rl = scarfinder!(leg, params; objective, measurement = SimpleBackend())
        @test length(rn.iterations) == length(rl.iterations) == 3
        for (itn, itl) in zip(rn.iterations, rl.iterations)
            @test itn.accepted == itl.accepted
            @test _summaries_match(itn.observables, itl.observables; atol = 1e-8)
            @test isapprox(itn.simple_score.score, itl.simple_score.score; atol = 1e-8)
        end
        @test rn.state isa SPD.PXPIPEPSState   # native state threaded through
    end

    @testset "native scarfinder! rejects unsupported features" begin
        cell = PeriodicSquareUnitCell(5, 5)
        trotter = TrotterParams(0.05, 2, :real, 4, 1e-12; schedule = :five_color)
        nat = product_pxp_pepskit_state(cell; state = :down, D = 1)
        # target_energy (imaginary-time correction) unsupported on native
        params_e = ScarFinderParams(0.05, trotter, 1, Inf, Inf, Inf, false;
                                    target_energy = -1.0, correction_time = 0.05,
                                    correction_attempts = 1)
        @test_throws ArgumentError scarfinder!(nat, params_e; measurement = SimpleBackend())
    end

    @testset "audit tensor_engine metadata" begin
        cfg = ScarFinderAuditConfig(;
            cell_Lx = 5, cell_Ly = 5, initial_state = :down,
            projection_time = 0.05, iterations = 2,
            dt_values = [0.05], evolve_maxdim = 1, target_maxdim_values = [1],
            schedule = :five_color, order = 2, tensor_engine = :pepskit_simple,
        )
        report = run_scarfinder_audit(cfg)
        @test !isempty(report.rows)
        @test all(row -> row.tensor_engine === :pepskit_simple, report.rows)
        @test all(row -> row.backend === :simple, report.rows)   # no CTM rows

        dir = mktempdir()
        csv = joinpath(dir, "audit.csv")
        json = joinpath(dir, "audit.json")
        write_scarfinder_audit_csv(report, csv)
        write_scarfinder_audit_json(report, json)
        csv_text = read(csv, String)
        json_text = read(json, String)
        @test occursin("tensor_engine", csv_text)
        @test occursin("pepskit_simple", csv_text)
        @test occursin("tensor_engine", json_text)
        @test occursin("pepskit_simple", json_text)

        # legacy default still records its engine
        cfg_legacy = ScarFinderAuditConfig(;
            cell_Lx = 5, cell_Ly = 5, projection_time = 0.05, iterations = 1,
            dt_values = [0.05], evolve_maxdim = 1, target_maxdim_values = [1],
            schedule = :five_color, order = 2,
        )
        @test cfg_legacy.tensor_engine === :legacy_itensors
    end

    @testset "audit native-engine restrictions" begin
        base = (; cell_Lx = 5, cell_Ly = 5, projection_time = 0.05, iterations = 1,
                dt_values = [0.05], schedule = :five_color, order = 2)
        # CTM (chi) not supported on native
        @test_throws ArgumentError ScarFinderAuditConfig(;
            base..., evolve_maxdim = 1, target_maxdim_values = [1],
            chi_values = [4], tensor_engine = :pepskit_simple)
        # per-iteration compression (target < evolve) not supported on native
        @test_throws ArgumentError ScarFinderAuditConfig(;
            base..., evolve_maxdim = 2, target_maxdim_values = [1],
            tensor_engine = :pepskit_simple)
        # unknown engine rejected
        @test_throws ArgumentError ScarFinderAuditConfig(;
            base..., evolve_maxdim = 1, target_maxdim_values = [1],
            tensor_engine = :nonsense)
    end
end
