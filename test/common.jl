using NEOs
using Dates
using TaylorSeries
using LinearAlgebra
using PlanetaryEphemeris
using Test

@testset "Common" begin

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
end