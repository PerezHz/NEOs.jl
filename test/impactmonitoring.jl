# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using TaylorSeries
using PlanetaryEphemeris
using Test

using NEOs: center, lbound, ubound, nominaltime, nominalstate, domain_radius,
          convergence_radius, isconvergent, timeofca

@testset "Impact monitoring" begin

    @testset "BPlane" begin
        # Fetch optical astrometry
        optical = fetch_optical_mpc80("2018 LA", MPC)
        # Parameters
        params = Parameters(
            coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
            gaussorder = 2, tsaorder = 2, adamiter = 500, adamQtol = 1e-5,
            jtlsorder = 2, jtlsiter = 200, lsiter = 1,
            significance = 0.99, outrej = false
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, optical)

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by Sep 19, 2025

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
        @test get_order(lov) == order

        # Close approaches
        σ, domain = 0.0, (-1.0, 1.0)
        nyears = 0.4 / yr
        ctol = 0.01
        CAs = closeapproaches(lov, σ, domain, nyears, params; ctol)
        @test length(CAs) == 1
        CA = CAs[1]

        @test σ in CA
        @test center(CA) == σ
        @test lbound(CA) == -1.0
        @test ubound(CA) == 1.0
        @test nominaltime(CA) == timeofca(CA, σ)
        @test nominalstate(CA) == targetplane(CA, σ)
        @test convergence_radius(CA, ctol) > domain_radius(CA) == 1.0
        @test isconvergent(CA, ctol)
        @test CA.tp == BPlane{Taylor1{Float64}}
        @test hypot(CA.x(), CA.y()) < CA.z()
    end

    @testset "Modified target plane" begin
        # Fetch optical astrometry
        optical = fetch_optical_mpc80("2024 PT5", MPC)
        filter!(x -> Date(2024, 10) < date(x) < Date(2024, 11), optical)
        # Parameters
        params = Parameters(
            coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
            gaussorder = 2, tsaorder = 2, adamiter = 500, adamQtol = 1e-5,
            jtlsorder = 2, jtlsiter = 200, lsiter = 1,
            significance = 0.99, outrej = false
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, optical[1:9])

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by Sep 19, 2025

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
        @test get_order(lov) == order

        # Close approaches
        σ, domain = 0.0, (-1.0, 1.0)
        nyears = 26 / yr
        ctol = 0.01
        CAs = closeapproaches(lov, σ, domain, nyears, params; ctol)
        @test length(CAs) == 1
        CA = CAs[1]

        @test σ in CA
        @test center(CA) == σ
        @test lbound(CA) == -1.0
        @test ubound(CA) == 1.0
        @test nominaltime(CA) == timeofca(CA, σ)
        @test nominalstate(CA) == targetplane(CA, σ)
        @test convergence_radius(CA, ctol) > domain_radius(CA) == 1.0
        @test isconvergent(CA, ctol)
        @test CA.tp == MTP{Taylor1{Float64}}
        @test hypot(CA.x(), CA.y()) > CA.z()
    end

end