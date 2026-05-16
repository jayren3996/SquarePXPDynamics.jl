const TRUST_CSV_HEADER = join(
    [
        "chi",
        "tol",
        "maxiter",
        "verbosity",
        "ctm_density",
        "ctm_blockade_violation",
        "ctm_pxp_energy_density",
        "ctm_iterations",
        "ctm_residual",
        "ctm_converged",
        "ctm_accepted",
        "trust_policy_min_points",
        "trust_policy_require_accepted_diagnostics",
        "trust_policy_max_density_delta",
        "trust_policy_max_blockade_delta",
        "trust_policy_max_energy_delta",
        "trust_policy_max_residual",
        "trust_trusted",
        "trust_reason",
        "trust_compared_points",
        "trust_finite_chi_density_delta",
        "trust_finite_chi_blockade_delta",
        "trust_finite_chi_energy_delta",
        "trust_observed_max_residual",
    ],
    ",",
)

const TRUST_REASONS = (
    :trusted,
    :too_few_points,
    :nonmonotonic_sweep,
    :missing_diagnostics,
    :unaccepted_diagnostics,
    :missing_residual,
    :residual_too_large,
    :density_delta_too_large,
    :blockade_delta_too_large,
    :energy_delta_too_large,
)

function _trust_diag_test(chi; tol = 1e-8, residual = tol / 10, accepted = true)
    return CTMRGDiagnostics(chi, tol, 100, 12, residual, true, accepted)
end

function _trust_summary_test(; density, blockade, energy, diag = _trust_diag_test(4))
    return CTMObservableSummary(density, density, density, blockade, energy, diag)
end

function _trust_point_test(;
    chi,
    tol,
    density,
    blockade,
    energy,
    reference_density = 0.0,
    reference_blockade = 0.0,
    reference_energy = 0.0,
    diag = _trust_diag_test(chi; tol),
)
    reference = CTMObservableSummary(
        reference_density,
        reference_density,
        reference_density,
        reference_blockade,
        reference_energy,
        CTMRGDiagnostics(chi, tol, 100, 12, tol / 10, true, true),
    )
    measurement = _trust_summary_test(
        density = density,
        blockade = blockade,
        energy = energy,
        diag = diag,
    )
    return CTMValidationPoint(PEPSKitCTMRGParams(chi, tol, 100, 0), reference, measurement)
end

@testset "CTM trust policy validation" begin
    policy = CTMTrustPolicy()

    @test policy.min_points == 2
    @test policy.require_accepted_diagnostics === true
    @test policy.max_density_delta ≈ 1e-3
    @test policy.max_blockade_delta ≈ 1e-4
    @test policy.max_energy_delta ≈ 1e-3
    @test policy.max_residual === nothing

    @test_throws ArgumentError CTMTrustPolicy(1, true, 1e-3, 1e-4, 1e-3, nothing)
    @test_throws ArgumentError CTMTrustPolicy(2, true, -1.0, 1e-4, 1e-3, nothing)
    @test_throws ArgumentError CTMTrustPolicy(2, true, 1e-3, NaN, 1e-3, nothing)
    @test_throws ArgumentError CTMTrustPolicy(2, true, 1e-3, 1e-4, Inf, nothing)
    @test_throws ArgumentError CTMTrustPolicy(2, true, 1e-3, 1e-4, 1e-3, -1e-8)
end

@testset "CTM trust assessment validation" begin
    assessment = CTMTrustAssessment(true, :trusted, "trusted", 2, 1e-5, 1e-6, 1e-5, 1e-9)

    @test assessment.trusted === true
    @test assessment.reason === :trusted
    @test assessment.compared_points == 2
    @test assessment.finite_chi_density_delta ≈ 1e-5
    @test assessment.observed_max_residual ≈ 1e-9

    for reason in TRUST_REASONS
        trusted = reason === :trusted
        stable = CTMTrustAssessment(trusted, reason, String(reason), 0, nothing, nothing, nothing, nothing)
        @test stable.reason === reason
    end

    @test_throws ArgumentError CTMTrustAssessment(
        true,
        :density_delta_too_large,
        "bad",
        2,
        1e-5,
        1e-6,
        1e-5,
        1e-9,
    )
    @test_throws ArgumentError CTMTrustAssessment(
        false,
        :bad_reason,
        "bad",
        2,
        nothing,
        nothing,
        nothing,
        nothing,
    )
    @test_throws ArgumentError CTMTrustAssessment(
        false,
        :too_few_points,
        "bad",
        -1,
        nothing,
        nothing,
        nothing,
        nothing,
    )
    @test_throws ArgumentError CTMTrustAssessment(
        false,
        :too_few_points,
        "bad",
        0,
        Inf,
        nothing,
        nothing,
        nothing,
    )
end

@testset "CTM trust accepts stable final chi window" begin
    points = [
        _trust_point_test(chi = 2, tol = 1e-5, density = 9.0, blockade = 3.0, energy = -4.0),
        _trust_point_test(
            chi = 4,
            tol = 1e-6,
            density = 10.0,
            blockade = 5.0,
            energy = -7.0,
            reference_density = -100.0,
            reference_blockade = -100.0,
            reference_energy = 100.0,
        ),
        _trust_point_test(
            chi = 8,
            tol = 1e-8,
            density = 10.0002,
            blockade = 5.00002,
            energy = -7.0002,
            reference_density = -100.0,
            reference_blockade = -100.0,
            reference_energy = 100.0,
        ),
    ]

    assessment = assess_ctm_trust(points)

    @test assessment.trusted === true
    @test assessment.reason === :trusted
    @test assessment.compared_points == 2
    @test assessment.finite_chi_density_delta ≈ 2e-4 atol = 1e-12
    @test assessment.finite_chi_blockade_delta ≈ 2e-5 atol = 1e-12
    @test assessment.finite_chi_energy_delta ≈ 2e-4 atol = 1e-12
    @test assessment.observed_max_residual !== nothing
    @test assessment.observed_max_residual <= 1e-6
    @test abs(points[end].delta_density) > 100
    @test abs(points[end].delta_pxp_energy_density) > 100
end

@testset "CTM trust uses maximum adjacent drift across custom final chi window" begin
    points = [
        _trust_point_test(chi = 2, tol = 1e-5, density = 9.0, blockade = 3.0, energy = -4.0),
        _trust_point_test(chi = 4, tol = 1e-6, density = 1.0, blockade = 0.01, energy = -0.3),
        _trust_point_test(chi = 8, tol = 1e-7, density = 1.0009, blockade = 0.01007, energy = -0.3008),
        _trust_point_test(chi = 16, tol = 1e-8, density = 1.0004, blockade = 0.01002, energy = -0.3003),
    ]
    policy = CTMTrustPolicy(3, true, 1e-3, 1e-4, 1e-3, nothing)

    assessment = assess_ctm_trust(points; policy)

    @test assessment.trusted === true
    @test assessment.reason === :trusted
    @test assessment.compared_points == 3
    @test assessment.finite_chi_density_delta ≈ 9e-4 atol = 1e-12
    @test assessment.finite_chi_blockade_delta ≈ 7e-5 atol = 1e-12
    @test assessment.finite_chi_energy_delta ≈ 8e-4 atol = 1e-12
    @test assessment.finite_chi_density_delta != abs(points[end].measurement.density - points[end - 2].measurement.density)
    @test assessment.finite_chi_density_delta != abs(points[end].measurement.density - points[end - 1].measurement.density)
end

@testset "CTM trust rejects expected trust failures" begin
    empty = assess_ctm_trust(CTMValidationPoint[])
    @test empty.trusted === false
    @test empty.reason === :too_few_points
    @test empty.compared_points == 0

    stable4 = _trust_point_test(chi = 4, tol = 1e-6, density = 0.1, blockade = 0.01, energy = -0.2)
    stable8 =
        _trust_point_test(chi = 8, tol = 1e-8, density = 0.1001, blockade = 0.01001, energy = -0.2001)

    too_few = assess_ctm_trust([stable4])
    @test too_few.trusted === false
    @test too_few.reason === :too_few_points
    @test too_few.compared_points == 1

    nonmonotonic = assess_ctm_trust([stable8, stable4])
    @test nonmonotonic.trusted === false
    @test nonmonotonic.reason === :nonmonotonic_sweep

    loose_tol = assess_ctm_trust([
        _trust_point_test(chi = 4, tol = 1e-8, density = 0.1, blockade = 0.01, energy = -0.2),
        _trust_point_test(chi = 8, tol = 1e-6, density = 0.1001, blockade = 0.01001, energy = -0.2001),
    ])
    @test loose_tol.trusted === false
    @test loose_tol.reason === :nonmonotonic_sweep

    missing_diag = assess_ctm_trust([
        _trust_point_test(chi = 4, tol = 1e-6, density = 0.1, blockade = 0.01, energy = -0.2),
        _trust_point_test(
            chi = 8,
            tol = 1e-8,
            density = 0.1001,
            blockade = 0.01001,
            energy = -0.2001,
            diag = nothing,
        ),
    ])
    @test missing_diag.trusted === false
    @test missing_diag.reason === :missing_diagnostics

    unaccepted = assess_ctm_trust([
        stable4,
        _trust_point_test(
            chi = 8,
            tol = 1e-8,
            density = 0.1001,
            blockade = 0.01001,
            energy = -0.2001,
            diag = _trust_diag_test(8; tol = 1e-8, accepted = false),
        ),
    ])
    @test unaccepted.trusted === false
    @test unaccepted.reason === :unaccepted_diagnostics

    missing_residual = assess_ctm_trust(
        [
            stable4,
            _trust_point_test(
                chi = 8,
                tol = 1e-8,
                density = 0.1001,
                blockade = 0.01001,
                energy = -0.2001,
                diag = CTMRGDiagnostics(8, 1e-8, 100, 12, nothing, true, true),
            ),
        ];
        policy = CTMTrustPolicy(2, true, 1e-3, 1e-4, 1e-3, 1e-6),
    )
    @test missing_residual.trusted === false
    @test missing_residual.reason === :missing_residual

    residual_large = assess_ctm_trust(
        [
            stable4,
            _trust_point_test(
                chi = 8,
                tol = 1e-8,
                density = 0.1001,
                blockade = 0.01001,
                energy = -0.2001,
                diag = _trust_diag_test(8; tol = 1e-8, residual = 1e-4),
            ),
        ];
        policy = CTMTrustPolicy(2, true, 1e-3, 1e-4, 1e-3, 1e-6),
    )
    @test residual_large.trusted === false
    @test residual_large.reason === :residual_too_large
end

@testset "CTM trust rejects finite-chi observable drift" begin
    base = _trust_point_test(chi = 4, tol = 1e-6, density = 0.2, blockade = 0.01, energy = -0.3)

    density_drift = assess_ctm_trust([
        base,
        _trust_point_test(chi = 8, tol = 1e-8, density = 0.202, blockade = 0.01001, energy = -0.3001),
    ])
    @test density_drift.trusted === false
    @test density_drift.reason === :density_delta_too_large

    blockade_drift = assess_ctm_trust([
        base,
        _trust_point_test(chi = 8, tol = 1e-8, density = 0.2001, blockade = 0.0102, energy = -0.3001),
    ])
    @test blockade_drift.trusted === false
    @test blockade_drift.reason === :blockade_delta_too_large

    energy_drift = assess_ctm_trust([
        base,
        _trust_point_test(chi = 8, tol = 1e-8, density = 0.2001, blockade = 0.01001, energy = -0.302),
    ])
    @test energy_drift.trusted === false
    @test energy_drift.reason === :energy_delta_too_large
end

@testset "CTM trust CSV audit output" begin
    points = [
        _trust_point_test(chi = 4, tol = 1e-6, density = 0.2, blockade = 0.01, energy = -0.3),
        _trust_point_test(chi = 8, tol = 1e-8, density = 0.2001, blockade = 0.01001, energy = -0.3001),
    ]
    policy = CTMTrustPolicy(2, true, 5e-5, 5e-6, 5e-5, nothing)
    path = tempname() * ".csv"

    write_ctm_trust_csv(points, path; policy)
    csv = read(path, String)
    lines = split(chomp(csv), '\n')
    rows = [split(line, ',') for line in lines[2:end]]
    header = split(lines[1], ',')

    policy_min_points_idx = findfirst(==("trust_policy_min_points"), header)
    policy_density_idx = findfirst(==("trust_policy_max_density_delta"), header)
    policy_blockade_idx = findfirst(==("trust_policy_max_blockade_delta"), header)
    policy_energy_idx = findfirst(==("trust_policy_max_energy_delta"), header)
    trusted_idx = findfirst(==("trust_trusted"), header)
    reason_idx = findfirst(==("trust_reason"), header)
    compared_points_idx = findfirst(==("trust_compared_points"), header)
    density_delta_idx = findfirst(==("trust_finite_chi_density_delta"), header)
    blockade_delta_idx = findfirst(==("trust_finite_chi_blockade_delta"), header)
    energy_delta_idx = findfirst(==("trust_finite_chi_energy_delta"), header)

    @test lines[1] == TRUST_CSV_HEADER
    @test length(lines) == 3
    @test occursin("8,1.0e-8,100", csv)
    @test !occursin("message", lowercase(csv))
    @test all(length(row) == length(header) for row in rows)
    @test all(row[policy_min_points_idx] == "2" for row in rows)
    @test all(row[policy_density_idx] == "5.0e-5" for row in rows)
    @test all(row[policy_blockade_idx] == "5.0e-6" for row in rows)
    @test all(row[policy_energy_idx] == "5.0e-5" for row in rows)
    @test all(row[trusted_idx] == "false" for row in rows)
    @test all(row[reason_idx] == "density_delta_too_large" for row in rows)
    @test all(row[compared_points_idx] == "2" for row in rows)
    @test all(!isempty(row[density_delta_idx]) for row in rows)
    @test all(!isempty(row[blockade_delta_idx]) for row in rows)
    @test all(!isempty(row[energy_delta_idx]) for row in rows)
end

@testset "CTM trust malformed inputs throw" begin
    @test_throws ArgumentError assess_ctm_trust(Any["not a point"])
end
