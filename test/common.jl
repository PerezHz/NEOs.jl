using NEOs
using Dates
using Preferences
using TaylorSeries
using LinearAlgebra
using PlanetaryEphemeris
using Test

@testset "Common" begin

    @testset "Solar system ephemerides" begin
        using NEOs: SSEPH_SOURCE, SSEPH_ARTIFACT_PATH, print_sseph_summary, set_sseph_source

        print_sseph_summary()

        @test SSEPH_SOURCE == SSEPH_ARTIFACT_PATH

        set_sseph_source("sseph343ast016_p100y_et.jld2")

        @test has_preference(NEOs, "SSEPH_SOURCE")
        @test load_preference(NEOs, "SSEPH_SOURCE") == "sseph343ast016_p100y_et.jld2"

        set_sseph_source(SSEPH_ARTIFACT_PATH)

        @test has_preference(NEOs, "SSEPH_SOURCE")
        @test load_preference(NEOs, "SSEPH_SOURCE") == SSEPH_ARTIFACT_PATH
    end

    @testset "Time conversions" begin
        t0 = now()
        @test et2dtutc(dtutc2et(t0)) == t0
        @test jdtdb2dtutc(dtutc2jdtdb(t0)) == t0
    end

    @testset "JPL ephemerides" begin
        using NEOs: sseph, ttmtdb

        loadjpleph()

        @test kmsec2auday(sunposvel(0.0)) ≈ sseph(su, 0.0)
        @test kmsec2auday(earthposvel(0.0)) ≈ sseph(ea, 0.0)
        @test kmsec2auday(moonposvel(0.0)) ≈ sseph(mo, 0.0)
        @test tt_tdb(0.0) ≈ ttmtdb(0.0)
        # Note: By Sep 26, 2025 the time derivative of TT-TDB differs
        # between JPL and PlanetaryEphemeris
        @test !(dtt_tdb(0.0) ≈ TS.differentiate(1, ttmtdb.x[1]))

        q_Apophis_JPL_220 = [-1.045062875223473E+00, -1.294565996082367E-01, -7.496573820257184E-02,
                              4.232752764404431E-03, -1.412783025595556E-02, -5.148374014688117E-03]
        s_Apophis_JPL_220 = [2.94843651E-09, 7.02530459E-09, 3.17546972E-09,
                             1.27304960E-10, 3.15390825E-11, 3.99604517E-11]

        @test norm((kmsec2auday(apophisposvel197(0.0)) - q_Apophis_JPL_220) ./ s_Apophis_JPL_220) < 250
        @test norm((kmsec2auday(apophisposvel199(0.0)) - q_Apophis_JPL_220) ./ s_Apophis_JPL_220) < 290
    end

    @testset "fasttaylors" begin
        using NEOs: AbstractBuffer, OpticalBuffer, auxzero, scalingfactor, auday2kmsec!,
              dot3D, euclid3D, evaleph, evaleph!, cte, sseph

        buffer = OpticalBuffer(0.0)
        @test isa(string(buffer), String)
        @test isa(buffer, AbstractBuffer)
        @test isa(buffer, OpticalBuffer{Float64})

        q00 = rand(6)
        scalings = fill(1E-8, 6)
        vars = set_variables(Float64, "dx"; order = 2, numvars = 6)
        q0T1 = q00 + scalings * Taylor1(2)
        q0TN = q00 + scalings .* vars
        q0T11 = [Taylor1(Taylor1([q00[i], rand()], 2), 2) for i in 1:6]
        q0T1N =  [Taylor1(q00[i] + dot(rand(6), vars), 2) for i in 1:6]

        @test all(iszero, auxzero.(q0T1))
        @test all(iszero, auxzero.(q0TN))
        @test all(iszero, auxzero.(q0T11))
        @test all(iszero, auxzero.(q0T1N))
        @test all(@. scalingfactor(q0TN) == scalings)
        @test all(@. q00 == cte(q0T1) == cte(q0TN) == cte(cte(q0T11)) == cte(cte(q0TN)))

        @test q00[1]*q0T1[1] + q00[2]*q0T1[2] + q00[3]*q0T1[3] == dot3D(q00, q0T1) ==
              dot3D(q0T1, q00)
        @test q00[1]*q0TN[1] + q00[2]*q0TN[2] + q00[3]*q0TN[3] == dot3D(q00, q0TN) ==
              dot3D(q0TN, q00)
        @test sqrt(q00[1]^2 + q00[2]^2 + q00[3]^2) == sqrt(dot3D(q00, q00)) ==
              euclid3D(q00)
        @test sqrt(q0T1[1]^2 + q0T1[2]^2 + q0T1[3]^2) == sqrt(dot3D(q0T1, q0T1)) ==
              euclid3D(q0T1)
        @test sqrt(q0TN[1]^2 + q0TN[2]^2 + q0TN[3]^2) == sqrt(dot3D(q0TN, q0TN)) ==
              euclid3D(q0TN)

        @test evaleph(sseph, Taylor1(2), Taylor1(2)) ≈ evaleph(sseph, Taylor1(2))

        evaleph!.(q0T1, q0T11, 0.0)
        evaleph!.(q0TN, q0T1N, 0.0)

        @test q0T1 == cte(q0T11)
        @test q0TN == cte(q0T1N)

        Q00 = deepcopy(q00)
        Q0T1 = deepcopy(q0T1)
        Q0TN = deepcopy(q0TN)

        auday2kmsec!(q00)
        auday2kmsec!(q0T1)
        auday2kmsec!(q0TN)

        @test auday2kmsec(Q00) ≈ q00
        @test auday2kmsec(Q0T1) ≈ q0T1
        @test auday2kmsec(Q0TN) ≈ q0TN

    end

end