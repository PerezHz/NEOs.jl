# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using PlanetaryEphemeris
using LinearAlgebra
using Test

using NEOs: NEOSolution, adaptative_maxsteps, scaled_variables
@testset "Gauss initial conditions" begin
    # Load observations 
    radec = read_radec_mpc(joinpath("data", "RADEC_2023_DW.dat"))
    # Abstolute tolerance
    abstol = 1e-20
    # Order (wrt time)
    order = 25
    # Use @taylorize
    parse_eqs = true
    
    # Gauss initial conditions
    sol = gaussinitcond(radec; order = order, abstol = abstol)

    @test isa(sol, NEOSolution{Float64, Float64})
    @test datetime2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
    @test datetime2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
    @test all( norm.(sol.bwd.x, Inf) .< 2 )
    @test all( norm.(sol.fwd.x, Inf) .< 2 )
    @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
    @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
    @test length(sol.res) == 123
    @test sol.fit.success
    @test all( sqrt.(diag(sol.fit.Γ)) .< 10 )
    @test nrms(sol) < 0.4

    # Orbit determination (with outlier rejection)
    sol = orbitdetermination(radec, sol; order = order, abstol = abstol, parse_eqs = parse_eqs)

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
end