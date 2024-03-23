# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using PlanetaryEphemeris
using LinearAlgebra
using Test

using NEOs: NEOSolution, numberofdays

@testset "Orbit Determination" begin
    @testset "Gauss Method" begin
        # Load observations
        radec = read_radec_mpc(joinpath("data", "RADEC_2023_DW.dat"))
        # Parameters
        params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007, parse_eqs = false)
        params = NEOParameters(params, parse_eqs=true)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by February 24, 2024

        # Vector of observations
        @test length(radec) == 123
        @test numberofdays(radec) < 21.0
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 44
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test datetime2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 123
        @test iszero(count(outlier.(sol.res)))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 4e-7 )
        @test nrms(sol) < 0.36
        # Scalig factors
        @test all(sol.scalings .== 1e-6)
        # Compatibility with JPL
        JPL = [-1.1003339484439327, 0.20772201506095814, 0.04202338912370205,
               -0.004735538686138557, -0.010626685053348663, -0.006016258344003866]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 6)

        # Refine orbit
        sol1 = orbitdetermination(radec, sol, params)

        # Orbit solution
        @test isa(sol1, NEOSolution{Float64, Float64})
        # Tracklets
        @test sol1.tracklets == sol.tracklets
        # Backward integration
        @test datetime2days(date(radec[1])) > sol1.bwd.t0 + sol1.bwd.t[end]
        @test all( norm.(sol1.bwd.x, Inf) .< 1.2 )
        @test isempty(sol1.t_bwd) && isempty(sol1.x_bwd) && isempty(sol1.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[end])) < sol1.fwd.t0 + sol1.fwd.t[end]
        @test all( norm.(sol1.fwd.x, Inf) .< 1.2 )
        @test isempty(sol1.t_fwd) && isempty(sol1.x_fwd) && isempty(sol1.g_fwd)
        # Vector of residuals
        @test length(sol1.res) == 123
        @test iszero(count(outlier.(sol1.res)))
        # Least squares fit
        @test sol1.fit.success
        @test all( sigmas(sol1) .< 5e-5 )
        @test nrms(sol1) < 0.36
        # Scalig factors
        @test all(sol1.scalings .< 2e-6)
        # Compatibility with JPL
        @test all(abs.(sol1() - JPL) ./ sigmas(sol1) .< 0.3)
    end

    @testset "Too Short Arc" begin
        # Fetch optical astrometry
        radec = fetch_radec_mpc("designation" => "2008 EK68")

        # Parameters
        params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by February 24, 2024

        # Vector of observations
        @test length(radec) == 10
        @test numberofdays(radec) < 0.05
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 1
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test datetime2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 10
        @test iszero(count(outlier.(sol.res)))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 6e-3 )
        @test nrms(sol) < 0.85
        # Scalig factors
        @test all(sol.scalings .< 1e-5)
        # Compatibility with JPL
        JPL = [-0.9698333701500199, 0.24036461256880043, 0.10288887522619743,
               -0.009512521373861719, -0.015325432152904881, -0.008094623534198382]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.1)
    end

    @testset "Outlier Rejection" begin
        # Fetch optical astrometry
        radec = fetch_radec_mpc("designation" => "2007 VV7")

        # Parameters
        params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007)

        # Orbit Determination
        sol = orbitdetermination(radec, params)
        # Outlier rejection
        sol = outlier_rejection(radec, sol, params)

        # Values by February 24, 2024

        # Vector of observations
        @test length(radec) == 21
        @test numberofdays(radec) < 3.03
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 4
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test datetime2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 21
        @test count(outlier.(sol.res)) == 2
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 5e-4 )
        @test nrms(sol) < 0.25
        # Scalig factors
        @test all(sol.scalings .< 8e-7)
        # Compatibility with JPL
        JPL = [0.7673358221902306, 0.6484904294813807, 0.2932331617634889,
               -0.011023358761553661, 0.015392684491034429, 0.006528836324700449]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.1)
    end

    @testset "Interesting NEOs" begin

        # 2014 AA hit the Earth around January 2, 2014, 02:49 UTC

        # Fetch optical astrometry
        radec = fetch_radec_mpc("designation" => "2014 AA")

        # Parameters
        params = NEOParameters(coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by February 24, 2024

        # Vector of observations
        @test length(radec) == 7
        @test numberofdays(radec) < 0.05
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 1
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test datetime2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 1e9 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 7
        @test iszero(count(outlier.(sol.res)))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 1e-3 )
        @test nrms(sol) < 0.13
        # Scalig factors
        @test all(sol.scalings .< 1e-5)
        # Compatibility with JPL
        JPL = [-0.17932853771087842, 0.8874166708545763, 0.38414497114153867,
               -0.01755788350351527, -0.005781328974619869, -0.0020073946363600814]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.3)

        # 2008 TC3 entered the Earth's atmosphere around October 7, 2008, 02:46 UTC

        # Fetch optical astrometry
        radec = fetch_radec_mpc("designation" => "2008 TC3")

        # Parameters
        params = NEOParameters(coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007)

        # Observations with <1" weight
        idxs = findall(x -> x < 1, w8sveres17.(radec))
        # Restricted Orbit Determination
        sol = orbitdetermination(radec[idxs], params)

        # Values by February 24, 2024

        # Vector of observations
        @test length(radec) == 883
        @test numberofdays(radec) < 0.80
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 2
        @test sol.tracklets[1].radec[1] == radec[idxs[1]]
        @test sol.tracklets[end].radec[end] == radec[idxs[end]]
        @test issorted(sol.tracklets)
        # Backward integration
        @test datetime2days(date(radec[idxs[1]])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[idxs[end]])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 1e4 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 20
        @test iszero(count(outlier.(sol.res)))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 2e-5 )
        @test nrms(sol) < 0.30
        # Scalig factors
        @test all(sol.scalings .< 1e-5)
        # Compatibility with JPL
        JPL = [0.9741070119227359, 0.21515061351517384, 0.09390897837680391,
               -0.007890445009307178, 0.016062726197198392, 0.006136042043681892]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.3)

        # Update orbit with more observations
        sol1 = orbitdetermination(radec[1:30], sol, params)

        # Orbit solution
        @test isa(sol1, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol1.tracklets) == 5
        @test sol1.tracklets[1].radec[1] == radec[idxs[1]]
        @test sol1.tracklets[end].radec[end] == radec[30]
        @test issorted(sol1.tracklets)
        # Backward integration
        @test datetime2days(date(radec[idxs[1]])) > sol1.bwd.t0 + sol1.bwd.t[end]
        @test all( norm.(sol1.bwd.x, Inf) .< 1 )
        @test isempty(sol1.t_bwd) && isempty(sol1.x_bwd) && isempty(sol1.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[idxs[end]])) < sol1.fwd.t0 + sol1.fwd.t[end]
        @test all( norm.(sol1.fwd.x, Inf) .< 1e4 )
        @test isempty(sol1.t_fwd) && isempty(sol1.x_fwd) && isempty(sol1.g_fwd)
        # Vector of residuals
        @test length(sol1.res) == 30
        @test iszero(count(outlier.(sol1.res)))
        # Least squares fit
        @test sol1.fit.success
        @test all( sigmas(sol1) .< 2e-6 )
        @test nrms(sol1) < 0.37
        # Scalig factors
        @test all(sol1.scalings .< 1e-6)
        # Compatibility with JPL
        @test all(abs.(sol1() - JPL) ./ sigmas(sol1) .< 1.6)
    end

end