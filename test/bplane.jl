# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using PlanetaryEphemeris
using LinearAlgebra
using Roots
using Test

using NEOs: rvelea

@testset "B-plane" begin

    # Fetch optical astrometry
    radec = fetch_optical_mpc80("2018 LA", MPC)
    # Parameters
    params = Parameters(coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007)
    # Orbit determination problem
    od = ODProblem(newtonian!, radec)

    # Initial Orbit Determination
    orbit = initialorbitdetermination(od, params)

    # Values by Jul 4, 2025

    # Time of close approach
    params = Parameters(params; fwdoffset = 0.3)
    bwd, fwd, res = propres(od, orbit(), epoch(orbit) + J2000, params)
    t_CA = find_zeros(t -> rvelea(fwd, params, t), fwd.t0, fwd.t0 + fwd.t[end-1])[1]
    # Asteroid's geocentric state vector
    xae = fwd(t_CA) - params.eph_ea(t_CA)
    # Earth's heliocentric state vector
    xes = params.eph_ea(t_CA) - params.eph_su(t_CA)

    # Öpik's coordinates
    B = bopik(xae, xes)
    # Modified Target Plane
    X, Y = mtp(xae)

    # Impact parameter
    @test B.b >= 1.0
    # Impact condition
    @test hypot(B.ξ, B.ζ) <= B.b
    @test hypot(X, Y) <= 1

end