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
        params = NEOParameters(abstol = 1e-20, order = 25, parse_eqs = true)
    
        # Orbit Determination
        sol = orbitdetermination(radec, params)
    
        # Values by December 21, 2023
        
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
        params = NEOParameters(abstol = 1e-20, order = 25, parse_eqs = true)
    
        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by December 22, 2023
        
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
        params = NEOParameters(abstol = 1e-20, order = 25, parse_eqs = true)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by December 23, 2023

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
    end

    @testset "Interesting NEOs" begin
        using NEOs: reduce_tracklets, selecteph, sseph, propres, evalfit

        # 2014 AA hit the Earth around January 2, 2014, 02:49 UTC

        # Optical astrometry file
        filename = joinpath("data", "2014_AA.txt")
        # Download observations
        get_radec_mpc("designation" => "2014 AA", filename)
        # Load observations 
        radec = read_radec_mpc(filename)
        # Delete astrometry file
        rm(filename)
        # Parameters
        params = NEOParameters(abstol = 1e-20, order = 25, parse_eqs = true)
    
        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by January 26, 2024
        
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

        # 2008 TC3 entered the Earth's atmosphere around October 7, 2008, 02:46 UTC

        # Optical astrometry file
        filename = joinpath("data", "2008_TC3.txt")
        # Download observations
        get_radec_mpc("designation" => "2008 TC3", filename)
        # Load observations 
        radec = read_radec_mpc(filename)
        # Delete astrometry file
        rm(filename)
        # Parameters
        params = NEOParameters(abstol = 1e-20, order = 25, parse_eqs = true,
                               coeffstol = Inf)

        # Observations before  October 6, 2008, 10:00 UTC
        idxs = findall(x -> date(x) < DateTime(2008, 10, 06, 10), radec)
        # Restricted Orbit Determination
        sol = orbitdetermination(radec[idxs], params)

        # Tracklets 
        tracklets = reduce_tracklets(radec)
        # Time of first (last) observation [julian days]
        t0, tf = datetime2julian(date(radec[1])), datetime2julian(date(radec[end]))
        # Sun (earth's) ephemeris
        eph_su = selecteph(NEOs.sseph, su)
        eph_ea = selecteph(sseph, ea)
        
        # Initial time of propagation [julian days]
        jd0 = sol.bwd.t0 + J2000
        # Years in backward (forward) propagation
        nyears_bwd = -(jd0 - t0 + 0.5) / yr
        nyears_fwd = (tf - jd0 + 0.5) / yr
        # Initial conditions
        q0 = sol()
        scalings = abs.(q0) ./ 10^6
        dq = NEOs.scaled_variables("dx", scalings, order = 5)
        q = q0 + dq
        # Propagation and residuals
        bwd = NEOs.propagate(RNp1BP_pN_A_J23E_J2S_eph_threads!, jd0, nyears_bwd, q, params)
        fwd = NEOs.propagate(RNp1BP_pN_A_J23E_J2S_eph_threads!, jd0, nyears_fwd, q, params)
        res = residuals(radec, params;
                        xvs = et -> auday2kmsec(eph_su(et/daysec)),
                        xve = et -> auday2kmsec(eph_ea(et/daysec)),
                        xva = et -> bwdfwdeph(et, bwd, fwd))
        # Least squares fit
        fit = tryls(res, zeros(6), params.niter)
        # Orbit with all the observations
        sol = evalfit(NEOSolution(tracklets, bwd, fwd, res, fit, scalings))

        # Initial time of propagation [julian days]
        jd0 = sol.bwd.t0 + J2000
        # Initial conditions
        q0 = sol()
        scalings = abs.(q0) ./ 10^6
        dq = NEOs.scaled_variables("dx", scalings, order = 5)
        q = q0 + dq
        # Propagation and residuals
        bwd = NEOs.propagate(RNp1BP_pN_A_J23E_J2S_eph_threads!, jd0, nyears_bwd, q, params)
        fwd = NEOs.propagate(RNp1BP_pN_A_J23E_J2S_eph_threads!, jd0, nyears_fwd, q, params)
        res = residuals(radec, params;
                        xvs = et -> auday2kmsec(eph_su(et/daysec)),
                        xve = et -> auday2kmsec(eph_ea(et/daysec)),
                        xva = et -> bwdfwdeph(et, bwd, fwd))
        # Least squares fit
        fit = tryls(res, zeros(6), params.niter)
        # Orbit refinement
        sol = evalfit(NEOSolution(tracklets, bwd, fwd, res, fit, scalings))

        # Values by January 26, 2024
        
        # Vector of observations
        @test length(radec) == 883
        @test numberofdays(radec) < 0.80
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 29
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test datetime2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 1e55 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 883
        @test iszero(count(outlier.(sol.res)))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 1e-7 )
        @test nrms(sol) < 0.61
        # Scalig factors
        @test all(sol.scalings .< 1e-6)
    end

end