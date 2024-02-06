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
        params = NEOParameters(abstol = 1e-20, order = 25, parse_eqs = true,
                               bwdoffset = 0.007, fwdoffset = 0.007)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by February 4, 2024

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
        # Compatibility with JPL
        JPL = [-1.100331943890894E+00, 2.077277940539948E-01, 4.202679172087372E-02,
               -4.735673360137306E-03, -1.062665933519790E-02, -6.016253165957075E-03]
        @test all(abs.(sol() - JPL) .< 27)
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
        params = NEOParameters(abstol = 1e-20, order = 25, parse_eqs = true,
                               bwdoffset = 0.007, fwdoffset = 0.007)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by February 4, 2024

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
        JPL = [-9.698293924982635E-01, 2.403718659163320E-01, 1.028928812918412E-01,
               -9.512665381191419E-03, -1.532539714362475E-02,  -8.094608965098163E-03]
        @test all(abs.(sol() - JPL) .< 1)
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
        params = NEOParameters(abstol = 1e-20, order = 25, parse_eqs = true,
                               bwdoffset = 0.007, fwdoffset = 0.007)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by February 4, 2024

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
        JPL = [7.673391189462907E-01, 6.484845064680100E-01, 2.932307930882956E-01,
               -1.102328578432002E-02, 1.539274577937339E-02, 6.528864069204219E-03]
        @test all(abs.(sol() - JPL) .< 1)
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
        params = NEOParameters(abstol = 1e-20, order = 25, parse_eqs = true,
                               bwdoffset = 0.007, fwdoffset = 0.007)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by February 4, 2024

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
        JPL = [-1.793234664811157E-01, 8.874182631964184E-01, 3.841456198855505E-01,
               -1.755789723380711E-02, -5.781199629059458E-03, -2.007345498817972E-03]
        @test all(abs.(sol() - JPL) .< 1)

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
                               coeffstol = Inf, fwdoffset = 0.007)

        # Observations with <1" weight
        idxs = findall(x -> x < 1, w8sveres17.(radec))
        # Restricted Orbit Determination
        sol = orbitdetermination(radec[idxs], params)

        # Values by February 4, 2024

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
        JPL = [9.741084137931751E-01, 2.151459813682736E-01, 9.390733559030655E-02,
               -7.890343235423515E-03, 1.606273839327320E-02, 6.136052979678407E-03]
        @test all(abs.(sol() - JPL) .< 9)
    end

end