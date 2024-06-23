# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using PlanetaryEphemeris
using LinearAlgebra
using Roots
using Test

using NEOs: propres

@testset "2018 LA" begin

    # Fetch optical astrometry
    radec = fetch_radec_mpc("designation" => "2018 LA")
    # Parameters
    params = NEOParameters(coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007)

    # Orbit Determination
    sol = orbitdetermination(radec[1:8], params)
    params = NEOParameters(params; jtlsorder = 6)
    sol = orbitdetermination(radec, sol, params)

    # Radial velocity with respect to the Earth.
    function rvelea(t, fwd, params)
        # Geocentric state vector
        rv = fwd(t) - params.eph_ea(t)
        # Derivative of geocentric distance
        return dot(rv[1:3], rv[4:6])
    end

    # Values by June 23, 2024

    # Time of close approach
    params = NEOParameters(params; fwdoffset = 0.3)
    bwd, fwd, res = propres(radec, epoch(sol) + J2000, sol(), params)
    t_CA = find_zeros(t -> rvelea(t, fwd, params), fwd.t0, fwd.t0 + fwd.t[end-1])[1]
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