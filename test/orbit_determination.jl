# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using PlanetaryEphemeris
using LinearAlgebra
using Test

using NEOs: NEOSolution, numberofdays

@testset "Orbit Determination" begin
    @testset "Gauss Method" begin
        # Load observations 
        radec = read_radec_mpc(joinpath("data", "RADEC_2023_DW.dat"))
        # Parameters
        params = Parameters(abstol = 1e-20, order = 25, parse_eqs = true)
    
        # Orbit Determination
        sol = orbitdetermination(radec, params)
    
        # Values by December 21, 2023
        
        # Vector of observations
        @test length(radec) == 123
        @test numberofdays(radec) < 21.0
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Observation nights
        @test length(sol.nights) == 44
        @test sol.nights[1].radec[1] == radec[1]
        @test sol.nights[end].radec[end] == radec[end]
        @test issorted(sol.nights)
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
        @test all( sigmas(sol) .< 1e-4 )
        @test nrms(sol) < 0.38
        # Scalig factors
        @test all(sol.scalings .== 1e-6)
    end

    @testset "Too Short Arc" begin
        # Optical astrometry file
        filename = joinpath("data", "2008_EK68.txt")
        # Download observations
        get_radec_mpc("designation" => "2008 EK68", filename)
        # Load observations 
        radec = read_radec_mpc(filename)
        # Delete astrometry file
        rm(filename)
        # Parameters
        params = Parameters(abstol = 1e-20, order = 25, parse_eqs = true)
    
        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by December 22, 2023
        
        # Vector of observations
        @test length(radec) == 10
        @test numberofdays(radec) < 0.05
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Observation nights
        @test length(sol.nights) == 1
        @test sol.nights[1].radec[1] == radec[1]
        @test sol.nights[end].radec[end] == radec[end]
        @test issorted(sol.nights)
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
    end

    @testset "Outlier Rejection" begin
        # Optical astrometry file
        filename = joinpath("data", "2007_VV7.txt")
        # Download observations
        get_radec_mpc("designation" => "2007 VV7", filename)
        # Load observations 
        radec = read_radec_mpc(filename)
        # Delete astrometry file
        rm(filename)
        # Parameters
        params = Parameters(abstol = 1e-20, order = 25, parse_eqs = true)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by December 23, 2023

        # Vector of observations
        @test length(radec) == 21
        @test numberofdays(radec) < 3.03
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Observation nights
        @test length(sol.nights) == 4
        @test sol.nights[1].radec[1] == radec[1]
        @test sol.nights[end].radec[end] == radec[end]
        @test issorted(sol.nights)
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
    end

end