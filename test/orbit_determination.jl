# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using PlanetaryEphemeris
using LinearAlgebra
using Test

using NEOs: NEOSolution, numberofdays
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

    # Gauss  refinement (with outlier rejection)
    #=
    sol = gauss_refinement(radec, sol, params)

    @test isa(sol, NEOSolution{Float64, Float64})
    @test datetime2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
    @test datetime2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
    @test all( norm.(sol.bwd.x, Inf) .< 2 )
    @test all( norm.(sol.fwd.x, Inf) .< 2 )
    @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
    @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
    @test length(sol.res) == 123
    @test sol.fit.success
    @test all( sqrt.(diag(sol.fit.Γ)) .< 100 )
    @test nrms(sol) < 0.4
    @test count(outlier.(sol.res)) == 0
    =#
end