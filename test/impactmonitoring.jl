# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using PlanetaryEphemeris
using LinearAlgebra
using Roots
using Test

@testset "Impact monitoring" begin

    # Fetch optical astrometry
    optical = fetch_optical_mpc80("2018 LA", MPC)
    # Parameters
    params = Parameters(
        coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
        tsaorder = 2, adamiter = 500, adamQtol = 1e-5,
        jtlsorder = 2, jtlsiter = 200, lsiter = 1,
        significance = 0.99, outrej = false
    )
    # Orbit determination problem
    od = ODProblem(newtonian!, optical)

    # Initial Orbit Determination
    orbit = initialorbitdetermination(od, params)

    # Values by Sep 16, 2025

    # Line of variations
    order, σmax = 12, 5.0
    lov = lineofvariations(od, orbit, params; order, σmax)
    @test lov.dynamics == od.dynamics
    @test epoch(lov) == epoch(orbit) + PE.J2000
    @test lov.domain == (-σmax, σmax)
    @test -σmax in lov
    @test 0.0 in lov
    @test σmax in lov
    @test lov(0.0) == orbit()

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
    BP = targetplane(B)
    # Modified Target Plane
    M = mtp(xae)
    MP = targetplane(M)

    # Impact parameter
    @test BP[3] >= MP[3] == 1.0
    # Impact condition
    @test hypot(BP[1], BP[2]) <= BP[3]
    @test hypot(MP[1], MP[2]) <= MP[3]

end