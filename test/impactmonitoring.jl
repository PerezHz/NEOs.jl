# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using Distributions
using TaylorSeries
using PlanetaryEphemeris
using Test

using NEOs: dynamicalmodel, opticalindices, numtypes, nominaltime, radius,
      initialcondition, nominalstate, difft, domain_radius, convergence_radius,
      convergence_domain, isconvergent, timeofca, distance, radialvelocity,
      concavity, width, isoutlov, vinf, ismarginal, closeapproaches, issamereturn

@testset "Impact monitoring" begin

    @testset "Common" begin
        using NEOs: PLANET_NAMES_TO_INDEX, PLANET_RADII, escapevelocity, sseph, numtype

        # Impact monitoring scales
        Es  = [1E2,  1E4,  1E5,  1E1,  1E3,  1E4,  1E7,  1E7,  1E1,   1E4,   1E7] # Mt
        IPs = [1E-6, 1E-4, 1E-3, 1E-1, 1E-1, 1E-1, 1E-3, 1E-1, 0.999, 0.999, 0.999]
        ΔT = 50.0

        TS = torinoscale.(Es, IPs)
        @test issorted(TS)
        @test TS == collect(0:10)

        PS = palermoscale.(Es, IPs, Ref(ΔT))
        perm = sortperm(@. IPs / Es^(-0.8))
        @test issorted(PS[perm])

        # Impact target
        t = dtutc2days(now())
        planets = Symbol.(first.(sort(collect(PLANET_NAMES_TO_INDEX), by = last)))
        for (i, planet) in enumerate(planets)
            target = ImpactTarget(:($planet))
            @test gm(target) == PE.μ[i]
            @test radius(target) == PLANET_RADII[i]
            @test escapevelocity(target) == sqrt(2 * gm(target) / radius(target)) * (au/daysec)
            @test target(t) == sseph(i, t)
        end

        # Target plane
        TP1 = BPlane(rand(5)...)
        TP2 = MTP(rand(2)...)

        @test numtype(TP1) == numtype(TP2) == Float64
        @test isa(string(TP1), String)
        @test isa(string(TP2), String)

        R, D = valsecchi_circle(1.0, 0.5, 0.0, 1, 1)
        @test isinf(R) && isinf(D)
    end

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

        # Values by Jan 15, 2026

        # Impact target
        target = ImpactTarget(:earth)
        # Impact monitoring problem
        IM = IMProblem(orbit, target)

        @test dynamicalmodel(IM) == newtonian!
        @test NEOs.dof(IM) == 6
        @test epoch(IM) == epoch(orbit)
        @test noptical(IM) == length(optical)
        @test opticalindices(IM) == eachindex(optical)
        @test gm(IM) == gm(target)
        @test mass(IM, params) == mass(orbit, params)
        @test escapevelocity(IM) == escapevelocity(target)
        @test minmaxdates(IM) == (date(optical[1]), date(optical[end]))

        # Line of variations
        order, σmax = 12, 5.0
        lov = lineofvariations(IM, params; order, σmax)

        @test numtypes(lov) == (typeof(newtonian!), Float64)
        @test dynamicalmodel(lov) == newtonian!
        @test epoch(lov) == nominaltime(lov) == epoch(orbit)
        @test round(date(lov), Minute) == round(days2dtutc(epoch(orbit)), Minute)
        @test sigma(lov) == 0.0
        @test get_order(lov) == order
        @test lov(0.0) == orbit()
        @test constant_term(lov(0.0, (-σmax, σmax))) == orbit()

        @test lov.domain == (-σmax, σmax)
        @test lbound(lov) == -σmax
        @test ubound(lov) == σmax
        @test width(lov) == 2σmax
        @test all(Base.Fix2(in, lov), [-σmax, 0.0, σmax])

        # Virtual asteroids
        N = 11
        R_TP, R_P, IP, Δσmax = 0.2, radius(target), 1E-5, 0.1
        VAs1 = virtualasteroids(lov, :uniform; N)
        VAs2 = virtualasteroids(lov, :normal; N)
        VAs3 = virtualasteroids(lov, :DelVigna19; R_TP, R_P, IP, Δσmax)

        @test sigma(VAs1[(N ÷ 2) + 1]) == 0.0
        @test lbound(VAs1[1]) == -σmax
        @test ubound(VAs1[end]) == σmax
        @test all(ubound(VAs1[i]) == lbound(VAs1[i+1]) for i in 1:N-1)
        d1 = Uniform(-σmax, σmax)
        ps = @. cdf(d1, lbound(VAs1)) - cdf(d1, ubound(VAs1))
        @test all(Base.Fix1(isapprox, ps[1]), ps)

        @test sigma(VAs2[(N ÷ 2) + 1]) == 0.0
        @test lbound(VAs2[1]) == -σmax
        @test ubound(VAs2[end]) == σmax
        @test all(ubound(VAs2[i]) == lbound(VAs2[i+1]) for i in 1:N-1)
        d2 = Normal(0.0, 1.0)
        ps = @. cdf(d2, lbound(VAs2)) - cdf(d2, ubound(VAs2))
        @test all(Base.Fix1(isapprox, ps[1]), ps)

        @test lbound(VAs3[1]) == -σmax
        @test ubound(VAs3[end]) == σmax
        @test all(ubound(VAs3[i]) == lbound(VAs3[i+1]) for i in 1:length(VAs3)-1)
        @test maximum(width, VAs3) < Δσmax

        N = 1
        VAs = virtualasteroids(lov, :uniform; N)
        VA = VAs[1]
        @test epoch(VA) == epoch(lov)
        @test nominaltime(VA) == nominaltime(lov)
        @test sigma(VA) == sigma(lov)
        @test initialcondition(VA) == lov(0.0, (-σmax, σmax))
        @test get_order(VA) == get_order(lov)

        # Close approaches
        nyears = 0.4 / yr
        vaorder = 6
        ctol = 0.01
        CAs = closeapproaches(IM, VA, nyears, params; vaorder, ctol)
        @test length(CAs) == 1
        CA = CAs[1]

        @test nominaltime(CA) == constant_term(CA.t) == timeofca(CA, 0.0)
        @test sigma(CA) == sigma(VA)
        @test get_order(CA) == vaorder < get_order(VA)
        @test nominalstate(CA) == [CA.x[0], CA.y[0], CA.z[0]] == targetplane(CA, 0.0)
        @test difft(CA, CA) == 0.0
        @test domain_radius(CA) == σmax

        @test 0.0 in CA
        @test lbound(CA) == lbound(VA)
        @test ubound(CA) == ubound(VA)
        @test convergence_radius(CA, ctol) > 1
        a, b = convergence_domain(CA, ctol)
        @test a < -σmax && σmax < b
        @test isconvergent(CA, ctol)
        @test CA.tp == BPlane{Taylor1{Float64}}

        # Returns
        RTs = showersnreturns(CAs)
        @test length(RTs) == 1
        RT = RTs[1]

        @test isa(string(RT), String)
        @test RT[1] == CA
        @test closeapproaches(RT) == RT[1:end] == [CA]
        @test firstindex(RT) == lastindex(RT) == 1
        @test !issamereturn(CA, CA, 45.0)
        @test showersnreturns(CAs) == showersnreturns([CAs])

        @test round(date(CA), Minute) == round(date(RT), Minute) == DateTime(2018, 06, 02, 16, 49)
        @test sigma(RT) == sigma(CA)
        @test lbound(RT) == lbound(CA)
        @test ubound(RT) == ubound(CA)

        @test 0.0 in RT
        @test get_order(RT) == get_order(CA)
        @test nominaltime(RT) == timeofca(RT, 0.0, ctol) == nominaltime(CA)
        @test nominalstate(RT, ctol) == targetplane(RT, 0.0, ctol) == nominalstate(CA)
        @test convergence_domain(RT, ctol) == convergence_domain(CA, ctol)
        @test isconvergent(RT, ctol)

        @test distance(CA, 0.0) == distance(RT, 0.0, ctol) < 0
        @test radialvelocity(CA, 0.0) == radialvelocity(RT, 0.0, ctol)
        @test concavity(CA, 0.0) ≈ concavity(RT, 0.0, ctol) > 0

        # Virtual impactors
        VIs = virtualimpactors(IM, lov, RTs, ctol, params)
        @test length(VIs) == 1
        VI1 = VIs[1]

        VI2 = VirtualImpactor(IM, lov, params, 0.0, nominaltime(RT), (0.0, 0.0))

        @test date(VI1) == date(VI2) == date(CA) == date(RT)
        @test sigma(VI1) == sigma(VI2) == 0.0
        @test impact_probability(VI1) ≈ impact_probability(VI2) atol = 0.01
        @test width(VI1) == 2σmax
        @test width(VI2) == 0.0
        @test !ismarginal(VI1) && !ismarginal(VI2)
        @test !isoutlov(VI1) && isoutlov(VI2)

        @test semiwidth(VI1) == semiwidth(VI2) > 0
        @test stretching(VI1) == stretching(VI2) > 0
        @test semimajoraxis(VI1) == semimajoraxis(VI2) < 0
        @test ishyperbolic(VI1) && ishyperbolic(VI2)

        @test vinf(IM, VI1) == vinf(IM, VI2) > 0
        @test impactenergy(IM, VI1, params) == impactenergy(IM, VI2, params) > 0
        @test palermoscale(IM, VI1, params) ≈ palermoscale(IM, VI2, params) atol=0.01
        @test torinoscale(IM, VI1, params) == torinoscale(IM, VI1, params) == 0

        VIs = [VI1, VI2]
        impactor_table(VIs)
        @test isa(summary(VIs), String)
        @test isa(string(VI1), String) && isa(string(VI2), String)
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

        # Values by Jan 15, 2026

        # Impact target
        target = ImpactTarget(:earth)
        # Impact monitoring problem
        IM = IMProblem(orbit, target)

        @test dynamicalmodel(IM) == newtonian!
        @test NEOs.dof(IM) == 6
        @test epoch(IM) == epoch(orbit)
        @test noptical(IM) == 9
        @test opticalindices(IM) == 1:9
        @test gm(IM) == gm(target)
        @test mass(IM, params) == mass(orbit, params)
        @test escapevelocity(IM) == escapevelocity(target)
        @test minmaxdates(IM) == (date(optical[1]), date(optical[9]))

        # Line of variations
        order, σmax = 12, 5.0
        lov = lineofvariations(IM, params; order, σmax)

        @test numtypes(lov) == (typeof(newtonian!), Float64)
        @test dynamicalmodel(lov) == newtonian!
        @test epoch(lov) == nominaltime(lov) == epoch(orbit)
        @test round(date(lov), Minute) == round(days2dtutc(epoch(orbit)), Minute)
        @test sigma(lov) == 0.0
        @test get_order(lov) == order
        @test lov(0.0) == orbit()
        @test constant_term(lov(0.0, (-σmax, σmax))) == orbit()

        @test lov.domain == (-σmax, σmax)
        @test lbound(lov) == -σmax
        @test ubound(lov) == σmax
        @test width(lov) == 2σmax
        @test all(Base.Fix2(in, lov), [-σmax, 0.0, σmax])

        # Virtual asteroids
        N = 11
        R_TP, R_P, IP, Δσmax = 0.2, radius(target), 1E-5, 0.1
        VAs1 = virtualasteroids(lov, :uniform; N)
        VAs2 = virtualasteroids(lov, :normal; N)
        VAs3 = virtualasteroids(lov, :DelVigna19; R_TP, R_P, IP, Δσmax)

        @test sigma(VAs1[(N ÷ 2) + 1]) == 0.0
        @test lbound(VAs1[1]) == -σmax
        @test ubound(VAs1[end]) == σmax
        @test all(ubound(VAs1[i]) == lbound(VAs1[i+1]) for i in 1:N-1)
        d1 = Uniform(-σmax, σmax)
        ps = @. cdf(d1, lbound(VAs1)) - cdf(d1, ubound(VAs1))
        @test all(Base.Fix1(isapprox, ps[1]), ps)

        @test sigma(VAs2[(N ÷ 2) + 1]) == 0.0
        @test lbound(VAs2[1]) == -σmax
        @test ubound(VAs2[end]) == σmax
        @test all(ubound(VAs2[i]) == lbound(VAs2[i+1]) for i in 1:N-1)
        d2 = Normal(0.0, 1.0)
        ps = @. cdf(d2, lbound(VAs2)) - cdf(d2, ubound(VAs2))
        @test all(Base.Fix1(isapprox, ps[1]), ps)

        @test lbound(VAs3[1]) == -σmax
        @test ubound(VAs3[end]) == σmax
        @test all(ubound(VAs3[i]) == lbound(VAs3[i+1]) for i in 1:length(VAs3)-1)
        @test maximum(width, VAs3) < Δσmax

        N = 1
        VAs = virtualasteroids(lov, :uniform; N)
        VA = VAs[1]
        @test epoch(VA) == epoch(lov)
        @test nominaltime(VA) == nominaltime(lov)
        @test sigma(VA) == sigma(lov)
        @test initialcondition(VA) == lov(0.0, (-σmax, σmax))
        @test get_order(VA) == get_order(lov)

        # Close approaches
        nyears = 26 / yr
        vaorder = 6
        ctol = 0.01
        CAs = closeapproaches(IM, VA, nyears, params; vaorder, ctol)
        @test length(CAs) == 1
        CA = CAs[1]

        @test nominaltime(CA) == constant_term(CA.t) == timeofca(CA, 0.0)
        @test sigma(CA) == sigma(VA)
        @test get_order(CA) == vaorder < get_order(VA)
        @test nominalstate(CA) == [CA.x[0], CA.y[0], CA.z[0]] == targetplane(CA, 0.0)
        @test difft(CA, CA) == 0.0
        @test domain_radius(CA) == σmax

        @test 0.0 in CA
        @test lbound(CA) == lbound(VA)
        @test ubound(CA) == ubound(VA)
        @test convergence_radius(CA, ctol) > 1
        a, b = convergence_domain(CA, ctol)
        @test a < -σmax && σmax < b
        @test isconvergent(CA, ctol)
        @test CA.tp == MTP{Taylor1{Float64}}

        # Returns
        RTs = showersnreturns(CAs)
        @test length(RTs) == 1
        RT = RTs[1]

        @test isa(string(RT), String)
        @test RT[1] == CA
        @test closeapproaches(RT) == RT[1:end] == [CA]
        @test firstindex(RT) == lastindex(RT) == 1
        @test !issamereturn(CA, CA, 45.0)
        @test showersnreturns(CAs) == showersnreturns([CAs])

        @test round(date(CA), Minute) == round(date(RT), Minute) == DateTime(2024, 10, 28, 23, 41)
        @test sigma(RT) == sigma(CA)
        @test lbound(RT) == lbound(CA)
        @test ubound(RT) == ubound(CA)

        @test 0.0 in RT
        @test get_order(RT) == get_order(CA)
        @test nominaltime(RT) == timeofca(RT, 0.0, ctol) == nominaltime(CA)
        @test nominalstate(RT, ctol) == targetplane(RT, 0.0, ctol) == nominalstate(CA)
        @test convergence_domain(RT, ctol) == convergence_domain(CA, ctol)
        @test isconvergent(RT, ctol)

        @test distance(CA, 0.0) == distance(RT, 0.0, ctol) > 0
        @test radialvelocity(CA, 0.0) == radialvelocity(RT, 0.0, ctol)
        @test concavity(CA, 0.0) ≈ concavity(RT, 0.0, ctol) < 0

        # Virtual impactors
        VIs = virtualimpactors(IM, lov, RTs, ctol, params)
        @test isempty(VIs)

        VI1 = VirtualImpactor(IM, lov, params, 0.0, nominaltime(RT), (-σmax, σmax))
        VI2 = VirtualImpactor(IM, lov, params, 0.0, nominaltime(RT), (0.0, 0.0))

        @test date(VI1) == date(VI2) == date(CA) == date(RT)
        @test sigma(VI1) == sigma(VI2) == 0.0
        @test impact_probability(VI1) > impact_probability(VI2) == 0.0
        @test width(VI1) == 2σmax
        @test width(VI2) == 0.0
        @test !ismarginal(VI1) && !ismarginal(VI2)
        @test !isoutlov(VI1) && isoutlov(VI2)

        @test semiwidth(VI1) == semiwidth(VI2) > 0
        @test stretching(VI1) == stretching(VI2) > 0
        @test semimajoraxis(VI1) == semimajoraxis(VI2) > 0
        @test !ishyperbolic(VI1) && !ishyperbolic(VI2)

        @test vinf(IM, VI1) == vinf(IM, VI2) == 0.0
        @test impactenergy(IM, VI1, params) == impactenergy(IM, VI2, params) > 0
        @test palermoscale(IM, VI1, params) > palermoscale(IM, VI2, params) == -Inf
        @test torinoscale(IM, VI1, params) == torinoscale(IM, VI1, params) == 0

        VIs = [VI1, VI2]
        impactor_table(VIs)
        @test isa(summary(VIs), String)
        @test isa(string(VI1), String) && isa(string(VI2), String)
    end

end