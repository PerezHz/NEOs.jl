# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using TaylorSeries
using PlanetaryEphemeris
using Test

using NEOs: nominaltime, nominalstate, domain_radius, convergence_radius,
      convergence_domain, isconvergent, timeofca, distance, rvelea,
      concavity, width, ismarginal, isoutlov

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

        # Values by Sep 25, 2025

        # Line of variations
        order, σmax = 12, 5.0
        lov = lineofvariations(od, orbit, params; order, σmax)

        @test round(date(lov), Minute) == DateTime(2018, 06, 02, 09, 34)
        @test sigma(lov) == 0.0
        @test lbound(lov) == -σmax
        @test ubound(lov) == σmax
        @test width(lov) == 2σmax

        @test lov.dynamics == od.dynamics
        @test nominaltime(lov) == epoch(lov) == epoch(orbit)
        @test lov.domain == (-σmax, σmax)
        @test -σmax in lov
        @test 0.0 in lov
        @test σmax in lov
        @test lov(0.0) == orbit()
        @test get_order(lov) == order

        # Close approaches
        σ, domain = 0.0, (-σmax, σmax)
        nyears = 0.4 / yr
        lovorder = 6
        ctol = 0.01
        CAs = closeapproaches(lov, σ, domain, nyears, params; lovorder, ctol)
        @test length(CAs) == 1
        CA = CAs[1]

        # Virtual asteroids
        VAs = virtualasteroids(CAs)
        @test length(VAs) == 1
        VA = VAs[1]

        @test round(date(CA), Minute) == round(date(VA), Minute) == DateTime(2018, 06, 02, 16, 49)
        @test sigma(CA) == sigma(VA) == σ
        @test lbound(CA) == lbound(VA) == -σmax
        @test ubound(CA) == ubound(VA) == σmax

        @test σ in CA && σ in VA
        @test get_order(CA) == get_order(VA) == lovorder
        @test nominaltime(CA) == nominaltime(VA) == timeofca(CA, σ) == timeofca(VA, σ, ctol)
        @test nominalstate(CA) == nominalstate(VA, ctol) == targetplane(CA, σ) == targetplane(VA, σ, ctol)
        @test domain_radius(CA) == σmax
        @test convergence_radius(CA, ctol) > 1
        @test convergence_domain(CA, ctol) == convergence_domain(CA, ctol)
        @test isconvergent(CA, ctol) && isconvergent(CA, ctol)
        @test CA.tp == BPlane{Taylor1{Float64}}
        @test distance(CA, σ) == distance(VA, σ, ctol) < 0
        @test rvelea(CA, σ) == rvelea(VA, σ, ctol)
        @test concavity(CA, σ) == concavity(VA, σ, ctol) > 0

        # Virtual impactors
        VIs = virtualimpactors(VAs, ctol, lov, od, orbit, params)
        @test isempty(VIs)

        VI1 = VirtualImpactor(lov, od, orbit, params, σ, nominaltime(VA), domain)
        VI2 = VirtualImpactor(lov, od, orbit, params, σ, nominaltime(VA), (σ, σ))

        @test date(VI1) == date(VI2) == date(CA) == date(VA)
        @test sigma(VI1) == sigma(VI2) == σ
        @test impact_probability(VI1) ≈ impact_probability(VI2) atol = 0.01
        @test width(VI1) == 2σmax
        @test width(VI2) == 0.0
        @test !ismarginal(VI1) && !ismarginal(VI2)
        @test !isoutlov(VI1) && isoutlov(VI2)

        VIs = [VI1, VI2]
        impactor_table(VIs)
        @test isa(summary(VIs), String)
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

        # Values by Sep 25, 2025

        # Line of variations
        order, σmax = 12, 5.0
        lov = lineofvariations(od, orbit, params; order, σmax)

        @test round(date(lov), Minute) == DateTime(2024, 10, 04, 20, 08)
        @test sigma(lov) == 0.0
        @test lbound(lov) == -σmax
        @test ubound(lov) == σmax
        @test width(lov) == 2σmax

        @test lov.dynamics == od.dynamics
        @test nominaltime(lov) == epoch(lov) == epoch(orbit)
        @test lov.domain == (-σmax, σmax)
        @test -σmax in lov
        @test 0.0 in lov
        @test σmax in lov
        @test lov(0.0) == orbit()
        @test get_order(lov) == order

        # Close approaches
        σ, domain = 0.0, (-σmax, σmax)
        nyears = 26 / yr
        lovorder = 6
        ctol = 0.01
        CAs = closeapproaches(lov, σ, domain, nyears, params; lovorder, ctol)
        @test length(CAs) == 1
        CA = CAs[1]

        # Virtual asteroids
        VAs = virtualasteroids(CAs)
        @test length(VAs) == 1
        VA = VAs[1]

        @test round(date(CA), Minute) == round(date(VA), Minute) == DateTime(2024, 10, 28, 23, 41)
        @test sigma(CA) == sigma(VA) == σ
        @test lbound(CA) == lbound(VA) == -σmax
        @test ubound(CA) == ubound(VA) == σmax

        @test σ in CA && σ in VA
        @test get_order(CA) == get_order(VA) == lovorder
        @test nominaltime(CA) == nominaltime(VA) == timeofca(CA, σ) == timeofca(VA, σ, ctol)
        @test nominalstate(CA) == nominalstate(VA, ctol) == targetplane(CA, σ) == targetplane(VA, σ, ctol)
        @test domain_radius(CA) == σmax
        @test convergence_radius(CA, ctol) > 1
        @test convergence_domain(CA, ctol) == convergence_domain(CA, ctol)
        @test isconvergent(CA, ctol) && isconvergent(CA, ctol)
        @test CA.tp == MTP{Taylor1{Float64}}
        @test distance(CA, σ) == distance(VA, σ, ctol) > 0
        @test rvelea(CA, σ) == rvelea(VA, σ, ctol)
        @test concavity(CA, σ) == concavity(VA, σ, ctol) < 0

        # Virtual impactors
        VIs = virtualimpactors(VAs, ctol, lov, od, orbit, params)
        @test isempty(VIs)

        VI1 = VirtualImpactor(lov, od, orbit, params, σ, nominaltime(VA), domain)
        VI2 = VirtualImpactor(lov, od, orbit, params, σ, nominaltime(VA), (σ, σ))

        @test date(VI1) == date(VI2) == date(CA) == date(VA)
        @test sigma(VI1) == sigma(VI2) == σ
        @test impact_probability(VI1) > impact_probability(VI2) == 0.0
        @test width(VI1) == 2σmax
        @test width(VI2) == 0.0
        @test !ismarginal(VI1) && !ismarginal(VI2)
        @test !isoutlov(VI1) && isoutlov(VI2)

        VIs = [VI1, VI2]
        impactor_table(VIs)
        @test isa(summary(VIs), String)
    end

end