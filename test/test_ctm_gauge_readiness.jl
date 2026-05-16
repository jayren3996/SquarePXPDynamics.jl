using LinearAlgebra

function _trusted_ctm_assessment()
    return CTMTrustAssessment(true, :trusted, "trusted", 2, 0.0, 0.0, 0.0, 0.0)
end

function _untrusted_ctm_assessment()
    return CTMTrustAssessment(
        false,
        :residual_too_large,
        "residual too large",
        2,
        0.0,
        0.0,
        0.0,
        1e-2,
    )
end

@testset "CTM gauge readiness public API exports" begin
    @test isdefined(SquarePXPDynamics, :CTMGaugePolicy)
    @test isdefined(SquarePXPDynamics, :CTMBondNormDiagnostic)
    @test isdefined(SquarePXPDynamics, :CTMGaugeReadiness)
    @test isdefined(SquarePXPDynamics, :BondGaugeFixInfo)
    @test isdefined(SquarePXPDynamics, :ctm_bond_norm_matrix)
    @test isdefined(SquarePXPDynamics, :ctm_bond_norm_diagnostic)
    @test isdefined(SquarePXPDynamics, :all_ctm_bond_norm_diagnostics)
    @test isdefined(SquarePXPDynamics, :ctm_ready_for_gauge_updates)
    @test isdefined(SquarePXPDynamics, :fix_bond_gauge!)
    @test isdefined(SquarePXPDynamics, :assert_fresh_pepskit_context)
end

@testset "CTM bond norm diagnostic validates synthetic matrices" begin
    bond = BondKey(SquareCoord(1, 1), :right)
    policy = CTMGaugePolicy()

    good = CTMBondNormDiagnostic(bond, ComplexF64[1.0 0.0; 0.0 0.25]; policy)
    @test good.accepted === true
    @test good.reject_reason === nothing
    @test good.hermiticity_residual ≈ 0.0 atol = 1e-14
    @test good.eigen_min ≈ 0.25 atol = 1e-14
    @test good.eigen_max ≈ 1.0 atol = 1e-14
    @test good.rcond ≈ 0.25 atol = 1e-14
    @test good.matrix ≈ adjoint(good.matrix) atol = 1e-14

    nonfinite = CTMBondNormDiagnostic(bond, ComplexF64[1.0 NaN; 0.0 1.0]; policy)
    @test nonfinite.accepted === false
    @test nonfinite.reject_reason === :nonfinite_entries

    nonhermitian = CTMBondNormDiagnostic(bond, ComplexF64[1.0 1.0; 0.0 1.0]; policy)
    @test nonhermitian.accepted === false
    @test nonhermitian.reject_reason === :nonhermitian
    @test nonhermitian.hermiticity_residual > policy.max_hermiticity_residual

    indefinite = CTMBondNormDiagnostic(bond, ComplexF64[1.0 0.0; 0.0 -0.1]; policy)
    @test indefinite.accepted === false
    @test indefinite.reject_reason === :indefinite

    strict = CTMGaugePolicy(; min_rcond = 0.5)
    ill_conditioned =
        CTMBondNormDiagnostic(bond, ComplexF64[1.0 0.0; 0.0 0.1]; policy = strict)
    @test ill_conditioned.accepted === false
    @test ill_conditioned.reject_reason === :ill_conditioned

    @test_throws ArgumentError CTMBondNormDiagnostic(bond, ComplexF64[1.0 0.0])
end

@testset "CTM bond norm diagnostics use fresh PEPSKit contexts" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = PEPSKitCTMRGParams(1, 1e-4, 1, 0)
    ctx = pepskit_ctmrg_context(psi; params)
    c = SquareCoord(2, 2)

    assert_fresh_pepskit_context(psi, ctx)
    right = ctm_bond_norm_diagnostic(psi, c, :right, ctx)
    left_alias = ctm_bond_norm_diagnostic(psi, neighbor(cell, c, :right), :left, ctx)
    up = ctm_bond_norm_diagnostic(psi, c, :up, ctx)
    down_alias = ctm_bond_norm_diagnostic(psi, neighbor(cell, c, :up), :down, ctx)

    for diag in (right, left_alias, up, down_alias)
        @test diag isa CTMBondNormDiagnostic
        @test size(diag.matrix) == (1, 1)
        @test diag.matrix ≈ adjoint(diag.matrix) atol = 1e-10
        @test isfinite(diag.frobenius_norm)
        @test diag.frobenius_norm > 0
        @test diag.accepted === true
        @test diag.reject_reason === nothing
        @test diag.eigen_min >= -1e-10
        @test diag.rcond >= 0
    end
    @test left_alias.bond == right.bond == bondkey(cell, c, :right)
    @test down_alias.bond == up.bond == bondkey(cell, c, :up)

    all_diags = all_ctm_bond_norm_diagnostics(psi, ctx)
    @test all_diags isa Dict{BondKey,CTMBondNormDiagnostic}
    @test Set(keys(all_diags)) == Set(keys(psi.link_weights))
    @test all(diag -> diag.accepted, values(all_diags))

    set_link_weight!(psi, c, :right, [1.0])
    @test_throws ArgumentError assert_fresh_pepskit_context(psi, ctx)
    @test_throws ArgumentError ctm_bond_norm_diagnostic(psi, c, :right, ctx)
end

@testset "CTM readiness combines trust freshness and norm quality" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    ctx = pepskit_ctmrg_context(psi; params = PEPSKitCTMRGParams(1, 1e-4, 1, 0))
    trusted = _trusted_ctm_assessment()
    untrusted = _untrusted_ctm_assessment()
    diags = all_ctm_bond_norm_diagnostics(psi, ctx)

    ready = ctm_ready_for_gauge_updates(psi, ctx, trusted; diagnostics = diags)
    @test ready isa CTMGaugeReadiness
    @test ready.ready === true
    @test ready.reason === :ready
    @test ready.trust === trusted

    untrusted_ready = ctm_ready_for_gauge_updates(psi, ctx, untrusted; diagnostics = diags)
    @test untrusted_ready.ready === false
    @test untrusted_ready.reason === :untrusted_ctm

    bad_bond = first(keys(diags))
    bad_diag = Dict(
        bad_bond => CTMBondNormDiagnostic(
            bad_bond,
            ComplexF64[1.0 1.0; 0.0 1.0],
        ),
    )
    bad_ready = ctm_ready_for_gauge_updates(
        psi,
        ctx,
        trusted;
        diagnostics = bad_diag,
        policy = CTMGaugePolicy(; require_all_bonds = false),
    )
    @test bad_ready.ready === false
    @test bad_ready.reason === :bad_bond_norm

    missing_ready = ctm_ready_for_gauge_updates(
        psi,
        ctx,
        trusted;
        diagnostics = Dict(bad_bond => diags[bad_bond]),
    )
    @test missing_ready.ready === false
    @test missing_ready.reason === :missing_bond_norm

    stale_psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    stale_ready = ctm_ready_for_gauge_updates(stale_psi, ctx, trusted; diagnostics = diags)
    @test stale_ready.ready === false
    @test stale_ready.reason === :stale_context
end

@testset "fix_bond_gauge product path is transactional no-op" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    ctx = pepskit_ctmrg_context(psi; params = PEPSKitCTMRGParams(1, 1e-4, 1, 0))
    diags = all_ctm_bond_norm_diagnostics(psi, ctx)
    before = measure_simple(psi)
    before_version = state_version(psi)

    info = fix_bond_gauge!(
        psi,
        SquareCoord(2, 2),
        :right,
        ctx,
        _trusted_ctm_assessment();
        diagnostics = diags,
    )
    after = measure_simple(psi)

    @test info isa BondGaugeFixInfo
    @test info.bond == bondkey(cell, SquareCoord(2, 2), :right)
    @test info.mutated === false
    @test info.reason === :product_noop
    @test info.readiness.ready === true
    @test state_version(psi) == before_version
    @test after.density ≈ before.density atol = 1e-14
    @test after.blockade_violation ≈ before.blockade_violation atol = 1e-14
    @test after.pxp_energy_density ≈ before.pxp_energy_density atol = 1e-14
    assert_fresh_pepskit_context(psi, ctx)
end

@testset "fix_bond_gauge D greater than one is explicit nonmutation in Slice 4" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    ctx = pepskit_ctmrg_context(psi; params = PEPSKitCTMRGParams(1, 1e-4, 1, 0))
    bond = bondkey(cell, SquareCoord(2, 2), :right)
    diag = CTMBondNormDiagnostic(bond, Matrix{ComplexF64}(I, 2, 2))
    before_version = state_version(psi)

    info = fix_bond_gauge!(
        psi,
        SquareCoord(2, 2),
        :right,
        ctx,
        _trusted_ctm_assessment();
        diagnostics = Dict(bond => diag),
        policy = CTMGaugePolicy(; require_all_bonds = false),
    )

    @test info.mutated === false
    @test info.reason === :d_greater_than_one_not_implemented
    @test info.readiness.ready === true
    @test state_version(psi) == before_version
end

