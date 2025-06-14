# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using LinearAlgebra
using Test

using NEOs: src_path

@testset "Observations" begin

    @testset "Time conversions" begin
        t0 = now()
        @test et2dtutc(dtutc2et(t0)) == t0
        @test jdtdb2dtutc(dtutc2jdtdb(t0)) == t0
    end

    @testset "Topocentric" begin
        using NEOs: TimeOfDay, sunriseset, obsposECEF, obsposvelECI

        # Ground observation
        radec_1 = read_radec_mpc("""
        99942        |C2012 12 12.33230011 28 40.300-26 29 32.10         17.70Vu~0mfl807
        99942        |C2012 12 12.33730011 28 38.970-26 29 34.80         17.60Vu~0mfl807
        99942        |C2012 12 12.34221011 28 37.640-26 29 37.50         17.50Vu~0mfl807
        99942        |C2012 12 12.34712011 28 36.330-26 29 40.00         17.50Vu~0mfl807
        99942        |C2012 12 12.35054011 28 35.400-26 29 41.80         17.50Vu~0mfl807
        """)
        # Sattellite observation
        radec_2 = read_radec_mpc("""
        99942         S2020 12 18.97667011 30 15.530-10 46 20.20         19.00RL~4ROFC51
        99942         s2020 12 18.9766701 - 5634.1734 - 2466.2657 - 3038.3924   ~4ROFC51
        99942         S2020 12 19.10732011 30 22.510-10 48 20.00               L~4ROFC51
        99942         s2020 12 19.1073201 - 5654.1816 - 2501.9465 - 2971.1902   ~4ROFC51
        99942         S2020 12 19.23810011 30 29.500-10 50 19.60               L~4ROFC51
        99942         s2020 12 19.2381001 - 5645.7831 - 2512.1036 - 2978.6411   ~4ROFC51
        99942         S2020 12 19.23822011 30 29.570-10 50 19.20               L~4ROFC51
        99942         s2020 12 19.2382201 - 5617.3465 - 2486.4031 - 3053.2209   ~4ROFC51
        """)

        # Check parsing
        @test length(radec_1) == 5
        @test all( map(x -> x.observatory.code, radec_1) .== "807")
        @test length(radec_2) == 4
        @test all( map(x -> x.observatory.code, radec_2) .== "C51")

        # TimeOfDay
        tod_1 = TimeOfDay.(radec_1)
        tod_2 = TimeOfDay.(radec_2)
        # Check
        @test allequal(tod_1)
        @test tod_1[1].light == :night
        @test tod_1[1].start == Date(2012, 12, 11)
        @test tod_1[1].stop == Date(2012, 12, 12)
        @test tod_1[1].utc == -5
        @test allunique(tod_2)
        @test all( getfield.(tod_2, :light) .== :space )
        @test all( date.(radec_2) .== getfield.(tod_2, :start) .== getfield.(tod_2, :start) )
        @test all( getfield.(tod_2, :utc) .== 0 )

        # Sunrise and sunset
        radec = read_radec_mpc("99942        8C2020 12 08.15001011 20 07.510-08 02 54.20         18.50GV~4ROF094")
        sun = sunriseset(radec[1])
        @test dtutc2jdtdb(sun[1]) ≈ dtutc2jdtdb(DateTime("2020-12-08T05:05:59.384"))
        @test dtutc2jdtdb(sun[2]) ≈ dtutc2jdtdb(DateTime("2020-12-08T14:05:49.386"))

        # obsposECEF
        ecef_2 = obsposECEF.(radec_2)
        @test ecef_2[1] ≈ [-3462.643557087632, 5076.197661798687, -3049.6756672719907]
        @test ecef_2[2] ≈ [1351.315736765706, 6027.937408384214, -2982.5146167937583]
        @test ecef_2[3] ≈ [5332.067839021762, 3112.403799578623, -2989.9547254809945]
        @test ecef_2[4] ≈ [5308.786202404402, 3079.725220466387, -3064.4773721684687]

        # obsposvelECI
        eci_2 = obsposvelECI.(radec_2)
        @test eci_2[1] == [-5634.1734, -2466.2657, -3038.3924, 0.0, 0.0, 0.0]
        @test eci_2[2] == [-5654.1816, -2501.9465, -2971.1902, 0.0, 0.0, 0.0]
        @test eci_2[3] == [-5645.7831, -2512.1036, -2978.6411, 0.0, 0.0, 0.0]
        @test eci_2[4] == [-5617.3465, -2486.4031, -3053.2209, 0.0, 0.0, 0.0]

    end

    @testset "Tracklet" begin
        using NEOs: reduce_tracklets

        # Choose this example because of the discontinuity in α

        # Fetch optical astrometry
        radec = fetch_radec_mpc("2020 TJ6")
        # Reduce tracklets
        tracklets = reduce_tracklets(radec)

        # Values by December 19, 2023
        @test length(tracklets) == 5
        @test getfield.(tracklets, :nobs) == [4, 4, 3, 4, 3]
        @test tracklets[3].α ≈ 6.2831 atol = 1e-4
        @test tracklets[3].observatory == search_obs_code("J04")
        @test tracklets[3].night.utc == -1
    end

    @testset "AbstractErrorModel" begin
        # Download optical astrometry
        radec = fetch_radec_mpc("2014 AA")

        # Weighting schemes
        w1 = UniformWeights(radec)
        w2 = Veres17(radec)
        w3 = ADESWeights(radec)
        w4 = NEOCCWeights(radec)
        w5 = NEODyS2Weights(radec)

        @test isa(w1, UniformWeights{Float64})
        @test isa(w2, Veres17{Float64})
        @test isa(w3, ADESWeights{Float64})
        @test isa(w4, NEOCCWeights{Float64})
        @test isa(w5, NEODyS2Weights{Float64})

        @test length(radec) == length(w1.w8s) == length(w2.w8s) == length(w3.w8s) ==
            length(w4.w8s) == length(w5.w8s)
        @test all(==((1.0, 1.0)), w1.w8s)
        @test all(x -> isapprox(x[1], 2.285, atol = 1e-3) &&
            isapprox(x[2], 2.285, atol = 1e-3), w2.w8s)
        @test all(==((4.0, 4.0)), w3.w8s)
        @test all(x -> isapprox(x[1], 2.288, atol = 1e-3) &&
            isapprox(x[2], 2.288, atol = 1e-3), w4.w8s)
        @test all(x -> isapprox(x[1], 2.288, atol = 1e-3) &&
            isapprox(x[2], 2.288, atol = 1e-3), w5.w8s)

        # Debiasing schemes
        d1 = Farnocchia15(radec)
        d2 = Eggl20(radec)
        d3 = ZeroDebiasing(radec)
        d4 = NEOCCDebiasing(radec)
        d5 = NEODyS2Debiasing(radec)

        @test isa(d1, Farnocchia15{Float64})
        @test isa(d2, Eggl20{Float64})
        @test isa(d3, ZeroDebiasing{Float64})
        @test isa(d4, NEOCCDebiasing{Float64})
        @test isa(d5, NEODyS2Debiasing{Float64})

        @test length(radec) == length(d1.bias) == length(d2.bias) == length(d3.bias) ==
            length(d4.bias) == length(d5.bias)
        @test all(x -> isapprox(x[1], 0.017, atol = 1e-3) &&
            isapprox(x[2], 0.027, atol = 1e-3), d1.bias)
        @test all(x -> isapprox(x[1], 0.021, atol = 1e-3) &&
            isapprox(x[2], -0.019, atol = 1e-3), d2.bias)
        @test all(==((0.0, 0.0)), d3.bias)
        @test all(==((0.017, 0.028)), d4.bias)
        @test all(==((0.017, 0.028)), d5.bias)

    end
end
