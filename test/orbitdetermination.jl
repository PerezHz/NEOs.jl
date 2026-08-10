# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using PlanetaryEphemeris
using LinearAlgebra
using StaticArraysCore
using Test

using NEOs: OpticalMPC80, RadarJPL, AbstractOpticalVector, RadarResidual, KeplerianElements,
      EquinoctialElements, AttributableElements, μ_S, indices, equatorial2ecliptic,
      ecliptic2equatorial, numtypes, sseph, covariance, scalartype, opticaltype, radartype,
      dof, hasradar

using Statistics: mean

const TEST_DATA = joinpath(pkgdir(NEOs), "test", "data")
const OpticalOrbit{T} = LeastSquaresOrbit{typeof(newtonian!), T, T, Vector{OpticalMPC80{T}}}
const RadarOrbit{T} = LeastSquaresOrbit{typeof(newtonian!), T, T, Vector{OpticalMPC80{T}},
    Vector{RadarJPL{T}}, Vector{RadarResidual{T, T}}}

function iodsuboptical(optical::AbstractOpticalVector, N::Int = 3)
    tracklets = reduce_tracklets(optical)
    idxs = indices(tracklets[1:N])
    suboptical = optical[idxs]
    return suboptical
end

# Maximum relative error
mre(x, y, z) = maximum(@. abs(x - y) / z)

# Mahalanobis distance
mahalanobis(x::AbstractVector, μ::AbstractVector, Σ::AbstractMatrix) =
        sqrt((x - μ)' * inv(Diagonal(Σ)) * (x - μ))

function jpl_compatibility_tests(orbit, params, bounds, JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
    @testset "Compatibility with JPL" begin
        # Reference epoch
        t0 = epoch(orbit)
        mjd0 = t0 + MJD2000
        # Barycentric cartesian sate vector
        q0, σ0 = orbit(), sigmas(orbit)
        @test mre(q0, JPL_CAR, σ0) < bounds[1]
        # Osculating elements
        kep = keplerian(orbit, params)
        eqn = equinoctial(orbit, params)
        attr = attributable(orbit, params)
        @test conicsection(kep) == conicsection(eqn) == :elliptic
        @test conicsection(attr) == :hyperbolic
        @test gm(kep) == gm(eqn) == μ_S
        @test gm(attr) == PE.μ[ea]
        @test epoch(kep) == epoch(eqn) == epoch(attr) == mjd0
        @test frame(kep) == frame(eqn) == :ecliptic
        @test frame(attr) == :equatorial
        @test isa(string(kep), String) && isa(string(eqn), String) &&
              isa(string(attr), String)
        kep0, σkep0 = elements(kep), sigmas(kep)
        eqn0, σeqn0 = elements(eqn), sigmas(eqn)
        attr0, σattr0 = elements(attr), sigmas(attr)
        @test mre(kep0, JPL_KEP, σkep0) < bounds[2]
        @test mre(eqn0, JPL_EQN, σeqn0) < bounds[2]
        @test mre(attr0, JPL_ATTR, σattr0) < bounds[2]
        q1 = equatorial2ecliptic(q0 - params.eph_su(t0))
        q2 = q0 - params.eph_ea(t0)
        @test mre(q1, keplerian2cartesian(kep0, mjd0; μ = μ_S), σ0) < bounds[3]
        @test mre(q1, equinoctial2cartesian(eqn0; μ = μ_S), σ0) < bounds[3]
        @test mre(q2, attributable2cartesian(attr0), σ0) < bounds[3]
        @test mre(kep(mjd0 + yr), eqn(mjd0 + yr), σ0) < bounds[3]
        @test mre(cartesian2keplerian(q1, mjd0; μ = μ_S), kep0, σkep0) < bounds[4]
        @test mre(cartesian2equinoctial(q1; μ = μ_S), eqn0, σeqn0) < bounds[4]
        @test mre(cartesian2attributable(q2), attr0, σattr0) < bounds[4]
        @test mre(keplerian2equinoctial(kep0, mjd0; μ = μ_S), eqn0, σeqn0) < bounds[5]
        @test mre(equinoctial2keplerian(eqn0, mjd0; μ = μ_S), kep0, σkep0) < bounds[5]
        _kep_ = keplerian(orbit, params, t0 - 30)
        _eqn_ = equinoctial(orbit, params, t0 - 30)
        _attr_ = attributable(orbit, params, t0 - 30)
        @test mahalanobis(elements(kep), elements(_kep_), covariance(kep)) > 1
        @test mahalanobis(elements(kep), elements(_kep_), covariance(_kep_)) > 1
        @test mahalanobis(elements(eqn), elements(_eqn_), covariance(eqn)) > 1
        @test mahalanobis(elements(eqn), elements(_eqn_), covariance(_eqn_)) > 1
        @test mahalanobis(elements(attr), elements(_attr_), covariance(attr)) > 1
        @test mahalanobis(elements(attr), elements(_attr_), covariance(_attr_)) > 1
    end
end

@testset "AbstractOsculatingElements" begin

    @testset "Rotation between the equatorial plane and the ecliptic" begin
        using NEOs: eq2ecl, ecl2eq, ϵ0_deg

        rv = [−0.18034828526, 0.94069105951, 0.34573599029,
              −0.0162659397882, 4.39154800E−5, −0.000395204013]
        drv = [7.12E−9, 1.94E−9, 5.41E−9, 4.79E−11, 6.72E−11, 1.38E−10]

        @test eq2ecl() * ecl2eq() ≈ I
        @test eq2ecl(-ϵ0_deg) * ecl2eq(-ϵ0_deg) ≈ I

        _rv_ = ecliptic2equatorial(equatorial2ecliptic(rv))
        @test mre(rv, _rv_, drv) < 6E-8
    end

    @testset "Pérez-Hernández & Benet (2022) Apophis OR7 orbit" begin
        # Pérez-Hernández & Benet (2022) Apophis OR7 orbit
        # See Tables 2 and 3 of the Supplementary Information in
        # https://doi.org/10.1038/s43247-021-00337-x

        # Reference epoch [TDB]
        jd0 = 2459200.5                      # JD
        mjd0 = jd0 + (MJD2000 - PE.J2000)       # MJD

        # Cartesian state vector [au, au/day]
        rv = [−0.18034828526, 0.94069105951, 0.34573599029,
            −0.0162659397882, 4.39154800E−5, −0.000395204013]
        drv = [7.12E−9, 1.94E−9, 5.41E−9, 4.79E−11, 6.72E−11, 1.38E−10]

        # Keplerian elements
        e, de = 0.19150886716, 1.60E-9
        q, dq = 0.74585305033, 1.54E−9                           # au
        tp, dtp = 2459101.04092537 + (MJD2000 - PE.J2000), 1.17E−6  # MJDTDB
        Ω, dΩ = 204.04199116, 8.81E−6                            # deg
        ω, dω = 126.65396094, 9.37E−6                            # deg
        i, di = 3.336773201, 1.74E−7                             # deg
        a, da = q / (1 - e), hypot(dq/(1-e), q*de/(1-e)^2)       # au
        M = rad2deg(sqrt(μ_S / a^3)) * (mjd0 - tp)               # deg
        dM = hypot(-3*M*da/(2a), -M*dtp/(jd0 - tp))
        A2, dA2 = -2.8988E-14, 2.48E-16                          # au/day^2
        adot, dadot = −199.0, 1.5                                # m/yr
        kep = KeplerianElements(μ_S, mjd0, :ecliptic,
            SVector{6}(a, e, i, ω, Ω, M),
            SMatrix{6, 6}(diagm([da^2, de^2, di^2, dω^2, dΩ^2, dM^2]))
        )

        @test numtypes(kep) == (Float64, Float64)
        @test isa(string(kep), String)
        @test !iscircular(kep) && iselliptic(kep) &&
              !isparabolic(kep) && !ishyperbolic(kep)
        @test conicsection(kep) == :elliptic
        @test gm(kep) == μ_S
        @test epoch(kep) == mjd0
        @test date(kep) == DateTime(2020, 12, 17, 00, 00)
        @test frame(kep) == :ecliptic
        @test sigmas(kep) ≈ [da, de, di, dω, dΩ, dM]

        @test semimajoraxis(kep) ≈ a
        @test pericenter(kep) ≈ q
        @test eccentricity(kep) ≈ e
        @test inclination(kep) ≈ i
        @test argperi(kep) ≈ ω
        @test longascnode(kep) ≈ Ω
        @test meananomaly(kep) ≈ M
        @test timeperipass(kep) ≈ tp

        # Average semimajor axis drift [m/yr]
        _adot_ = 1E3au * yr * yarkp2adot(A2, a, e)
        _dadot_ = hypot(_adot_*dA2/A2, -_adot_*da/(2a), 2*_adot_*e*de/(1 - e^2))
        @test abs((adot - _adot_) / adot) < 1.6E-4
        @test abs((dadot - _dadot_) / dadot) < 1.4E-1
    end

    @testset "JPL 3I/ATLAS #44 Orbit" begin
        # JPL 3I/ATLAS #44 orbit
        # See https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=3I

        # Reference epoch [TDB]
        jd0 = 2460894.5                      # JD
        mjd0 = jd0 + (MJD2000 - PE.J2000)       # MJD

        # Keplerian elements
        a, da = -0.2639243367163182, 5.5149E-7                          # au
        q, dq = 1.356418761995381, 1.6452E-6                            # au
        e, de = 6.139422831829797, 1.4616E-5
        i, di = 175.1130917268881, 1.3912E-5	                        # deg
        ω, dω = 128.0096924001076, 0.00021962	                        # deg
        Ω, dΩ =	322.1566239181344, 0.00021152	                        # deg
        M, dM = -606.8465596139827, 0.0016883                           # deg
        tp, dtp = 2460977.982217865578 + (MJD2000 - PE.J2000), 3.133E-5    # MJDTDB
        kep = KeplerianElements(μ_S, mjd0, :ecliptic,
            SVector{6}(q, e, i, ω, Ω, tp),
            SMatrix{6, 6}(diagm([dq^2, de^2, di^2, dω^2, dΩ^2, dtp^2]))
        )

        @test numtypes(kep) == (Float64, Float64)
        @test isa(string(kep), String)
        @test !iscircular(kep) && !iselliptic(kep) &&
              !isparabolic(kep) && ishyperbolic(kep)
        @test conicsection(kep) == :hyperbolic
        @test gm(kep) == μ_S
        @test epoch(kep) == mjd0
        @test date(kep) == DateTime(2025, 8, 7, 00, 00)
        @test frame(kep) == :ecliptic
        @test sigmas(kep) ≈ [dq, de, di, dω, dΩ, dtp]

        @test semimajoraxis(kep) ≈ a
        @test pericenter(kep) ≈ q
        @test eccentricity(kep) ≈ e
        @test inclination(kep) ≈ i
        @test argperi(kep) ≈ ω
        @test longascnode(kep) ≈ Ω
        @test meananomaly(kep) ≈ M
        @test timeperipass(kep) ≈ tp

        # Conversions to other sets of elements
        kep0, σkep0 = elements(kep), sigmas(kep)
        kep1 = cartesian2keplerian(keplerian2cartesian(kep0, mjd0; μ = μ_S), mjd0; μ = μ_S)
        kep2 = equinoctial2keplerian(keplerian2equinoctial(kep0, mjd0; μ = μ_S), mjd0; μ = μ_S)
        @test mre(kep0, kep1, σkep0) < 3.1E-10
        @test mre(kep0, kep2, σkep0) < 1.3E-10
    end

end

@testset "Orbit Determination" begin

    @testset "Straight Gauss Method" begin
        # Load observations
        optical = read_optical_mpc80(joinpath(TEST_DATA, "2023DW.txt"))
        # Subset of optical for IOD
        suboptical = iodsuboptical(optical, 3)

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           gaussorder = 2, jtlsorder = 2, jtlsiter = 20, lsiter = 20,
           significance = 0.99, outrej = false, parse_eqs = false
        )
        params = Parameters(params, parse_eqs = true)

        # Orbit determination problem
        od = ODProblem(newtonian!, suboptical)

        @test isa(string(od), String)
        @test scalartype(od) == Float64
        @test opticaltype(od) == OpticalMPC80{Float64}
        @test radartype(od) == Nothing && !hasradar(od)
        @test dof(od) == 6
        @test NEOs.optical(od) == suboptical && isnothing(NEOs.radar(od))

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by August 9, 2026

        # Check type
        @test isa(orbit, OpticalOrbit{Float64})
        # Tracklets
        @test designation(orbit) == "2023 DW"
        @test length(suboptical) == nobs(od) == nobs(orbit) == 9
        @test numberofdays(suboptical) == numberofdays(orbit) < 0.18
        @test minmaxdates(orbit) == (date(suboptical[1]), date(suboptical[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 3
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].indices[1] == 1
        @test orbit.tracklets[end].indices[end] == length(suboptical)
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), dtutc2days(date(od.tracklets[2])), atol = 4e-4)
        @test dtutc2days(date(suboptical[1])) > firsttime(orbit)
        @test all( norm.(orbit.bwd.p, Inf) .< 2 )
        @test dtutc2days(date(suboptical[end])) < lasttime(orbit)
        @test all( norm.(orbit.fwd.p, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 9
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test isa(string(orbit.fit), String)
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 9e-4 )
        @test all( snr(orbit) .> 14.5)
        @test chi2(orbit) < 0.53
        @test nrms(orbit) < 0.18
        # Covariance matrix
        Γ = covariance(orbit)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL_CAR = [-9.867822372310442E-01, 3.781927033166271E-01, 1.409429196233261E-01,
            -8.773480716697608E-03, -9.470221311888708E-03, -5.654113212253916E-03]
        JPL_GEO = [-6.640076709353979E-02, 2.398619117376862E-02, -1.283857914434059E-02,
            -1.782970529237847E-03, 5.133750583261109E-03, 6.771648773014439E-04]
        JPL_KEP = [8.198694885290753E-01, 3.963590496724497E-01, 5.808386294295408E+00,
            4.045222425970146E+01, 3.261410820553623E+02, 1.216619661909267E+02]
        JPL_EQN = keplerian2equinoctial(JPL_KEP, epoch(orbit) + MJD2000; μ = μ_S)
        JPL_ATTR = cartesian2attributable(JPL_GEO)
        jpl_compatibility_tests(orbit, params, (7.6E-01, 8.3E-01, 7.0E-12, 4.6E-14, 1.1E-13),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 24.3 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit, params)
        @test 34.2 < Da < Dc < Db < 77.0
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit, params)
        @test 5.4E7 < Ma < Mc < Mb < 6.3E8
        # MPEC
        @test isnothing(print_mpec(orbit, params))
        println()

        # Add observations
        suboptical = iodsuboptical(optical, 15)

        # Refine orbit
        NEOs.update!(od, suboptical)
        orbit1 = orbitdetermination(od, orbit, params)

        # Check type
        @test isa(orbit1, OpticalOrbit{Float64})
        # Tracklets
        @test designation(orbit1) == "2023 DW"
        @test length(suboptical) == nobs(od) == nobs(orbit1) == 43
        @test numberofdays(suboptical) == numberofdays(orbit1) < 2.76
        @test minmaxdates(orbit1) == (date(suboptical[1]), date(suboptical[end]))
        @test length(od.tracklets) == length(orbit1.tracklets) == 15
        @test od.tracklets == orbit1.tracklets
        @test orbit1.tracklets[1].indices[1] == 1
        @test orbit1.tracklets[end].indices[end] == length(suboptical)
        @test issorted(orbit1.tracklets)
        # Backward (forward) integration
        @test epoch(orbit1) == epoch(orbit)
        @test dtutc2days(date(suboptical[1])) > firsttime(orbit1)
        @test all( norm.(orbit1.bwd.p, Inf) .< 1.2 )
        @test dtutc2days(date(suboptical[end])) < lasttime(orbit1)
        @test all( norm.(orbit1.fwd.p, Inf) .< 1.2 )
        # Vector of residuals
        @test notout(orbit1.ores) == 43
        @test nout(orbit1.ores) == 0
        # Least squares fit
        @test isa(string(orbit1.fit), String)
        @test orbit1.fit.success
        @test all( sigmas(orbit1) .< 2e-4 )
        @test all( snr(orbit1) .> 866 )
        @test chi2(orbit1) < 11.64
        @test nrms(orbit1) < 0.37
        # Covariance matrix
        Γ = covariance(orbit1)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit1.qs, 1) == 6
        @test size(orbit1.qs, 2) == length(orbit1.Qs) <= 5
        @test issorted(orbit1.Qs, rev = true)
        @test orbit1.Qs[end] == nrms(orbit1)
        # Compatibility with JPL
        jpl_compatibility_tests(orbit1, params, (3.1E-01, 4.5E-01, 5.7E-11, 8.2E-12, 6.3E-12),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit1, params)
        @test H - dH ≤ 24.3 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit1, params) == 7
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit1, params)
        @test 36.2 < Da < Dc < Db < 81.4
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit1, params)
        @test 6.4E7 < Ma < Mc < Mb < 7.4E8
        # MPEC
        @test isnothing(print_mpec(orbit1, params))
        println()
    end

    @testset "Unsafe Gauss Method" begin
        # Load observations
        optical = read_optical_mpc80(joinpath(TEST_DATA, "2005TM173.txt"))
        # Subset of optical for IOD
        # suboptical = iodsuboptical(optical, 3)
        # In this case: optical == suboptical

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           gaussorder = 2, jtlsorder = 2, jtlsiter = 20, lsiter = 10,
           significance = 0.99, outrej = false, safegauss = false
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, optical)

        @test isa(string(od), String)
        @test scalartype(od) == Float64
        @test opticaltype(od) == OpticalMPC80{Float64}
        @test radartype(od) == Nothing && !hasradar(od)
        @test dof(od) == 6
        @test NEOs.optical(od) == optical && isnothing(NEOs.radar(od))

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by August 9, 2026

        # Check type
        @test isa(orbit, OpticalOrbit{Float64})
        # Tracklets
        @test designation(orbit) == "2005 TM173"
        @test length(optical) == nobs(od) == nobs(orbit) == 6
        @test numberofdays(optical) == numberofdays(orbit) < 1.95
        @test minmaxdates(orbit) == (date(optical[1]), date(optical[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 2
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].indices[1] == 1
        @test orbit.tracklets[end].indices[end] == length(optical)
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), dtutc2days(date(optical[4])), atol = 3e-4)
        @test dtutc2days(date(optical[1])) > firsttime(orbit)
        @test all( norm.(orbit.bwd.p, Inf) .< 2 )
        @test dtutc2days(date(optical[end])) < lasttime(orbit)
        @test all( norm.(orbit.fwd.p, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 6
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test isa(string(orbit.fit), String)
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 5e-3 )
        @test all( snr(orbit) .> 21.4)
        @test chi2(orbit) < 2.53
        @test nrms(orbit) < 0.46
        # Covariance matrix
        Γ = covariance(orbit)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL_CAR = [1.004256184830538E+00, 2.231635895303758E-01, 1.151385598814907E-01,
            -1.082421281819635E-02, 1.742879823300478E-02, 7.104678045458838E-03]
        JPL_GEO = [3.977197755123790E-02, -3.014438907068757E-02, 5.444636306828039E-03,
            -5.809136993789294E-03, 2.303264464230069E-03, 5.470042538831260E-04]
        JPL_KEP = [2.872424697642789E+00, 6.749395051551541E-01, 1.282355986214476E+00,
            1.725712172730245E+02, 2.413793589329106E+02, 3.538743602962668E+02]
        JPL_EQN = keplerian2equinoctial(JPL_KEP, epoch(orbit) + MJD2000; μ = μ_S)
        JPL_ATTR = cartesian2attributable(JPL_GEO)
        jpl_compatibility_tests(orbit, params, (6.0E-03, 6.0E-03, 2.4E-11, 1.4E-14, 3.8E-12),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 24.0 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit, params)
        @test 41.1 < Da < Dc < Db < 92.3
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit, params)
        @test 9.4E7 < Ma < Mc < Mb < 1.2E9
        # MPEC
        @test isnothing(print_mpec(orbit, params))
        println()
    end

    @testset "Gauss Method with ADAM refinement" begin
        # Load observations
        optical = read_optical_mpc80(joinpath(TEST_DATA, "2024MK.txt"))
        # Subset of optical for IOD
        suboptical = optical[10:21]

        # Parameters
        params = Parameters(
            coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
            gaussorder = 2, safegauss = true,
            tsaorder = 2, adamiter = 500, adamQtol = 1e-5, jtlsorder = 4,
            jtlsmask = false, jtlsiter = 20, lsiter = 10, significance = 0.99,
            outrej = true, χ2_rec = sqrt(9.21), χ2_rej = sqrt(10),
            fudge = 100.0, max_per = 34.0,
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, suboptical)

        @test isa(string(od), String)
        @test scalartype(od) == Float64
        @test opticaltype(od) == OpticalMPC80{Float64}
        @test radartype(od) == Nothing && !hasradar(od)
        @test dof(od) == 6
        @test NEOs.optical(od) == suboptical && isnothing(NEOs.radar(od))

        # Initial Orbit Determination
        orbit = gaussiod(od, params)

        # Values by August 9, 2026

        # Check type
        @test isa(orbit, OpticalOrbit{Float64})
        # Tracklets
        @test designation(orbit) == "2024 MK"
        @test length(suboptical) == nobs(od) == nobs(orbit) == 12
        @test numberofdays(suboptical) == numberofdays(orbit) < 42.8
        @test minmaxdates(orbit) == (date(suboptical[1]), date(suboptical[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 3
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].indices[1] == 1
        @test orbit.tracklets[end].indices[end] == length(suboptical)
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), dtutc2days(date(od.tracklets[2])), atol = 5e-4)
        @test dtutc2days(date(suboptical[1])) > firsttime(orbit)
        @test all( norm.(orbit.bwd.p, Inf) .< 2 )
        @test dtutc2days(date(suboptical[end])) < lasttime(orbit)
        @test all( norm.(orbit.fwd.p, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 12
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test isa(string(orbit.fit), String)
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 6.6e-4 )
        @test all( snr(orbit) .> 38.8)
        @test chi2(orbit) < 2.43
        @test nrms(orbit) < 0.32
        # Covariance matrix
        Γ = covariance(orbit)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 3
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL_CAR = [-1.272253324405492E-01, -9.466101191080750E-01, -4.526816519305000E-01,
            2.048875630636477E-02, -2.272010701451988E-04, 3.213028469127346E-03]
        JPL_GEO = [-4.582770483046963E-02, -1.328523323350602E-02, -4.831618660681472E-02,
            3.607102715393575E-03, 9.878186144576845E-04, 3.739277819983182E-03]
        JPL_KEP = [2.232655267540509E+00, 5.480018623993814E-01, 8.456323361249829E+00,
            1.322293109123257E+01, 2.778787983622283E+02, 3.530120833962170E+02]
        JPL_EQN = keplerian2equinoctial(JPL_KEP, epoch(orbit) + MJD2000; μ = μ_S)
        JPL_ATTR = cartesian2attributable(JPL_GEO)
        jpl_compatibility_tests(orbit, params, (1.6E-01, 2.4E0, 2.5E-11, 4.3E-12, 1.5E-12),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 21.7 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit, params)
        @test 120.2 < Da < Dc < Db < 269.2
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit, params)
        @test 2.2E9 < Ma < Mc < Mb < 2.8E10
        # MPEC
        @test isnothing(print_mpec(orbit, params))
        println()
    end

    @testset "Admissible region" begin
        using NEOs: AdmissibleRegion, arenergydis, rangerate, rangerates,
            argoldensearch, arboundary, _helmaxrange, R_SI, k_gauss, μ_ES,
            boundary_projection, topo2bary, bary2topo, arW, ardW, ard2W,
            arS, ardS, ard2S, arG

        # Read optical astrometry
        optical = read_optical_mpc80(joinpath(TEST_DATA, "2024BX1.txt"))
        # Parameters
        params = Parameters()
        # First tracklet
        optical = optical[1:3]
        tracklet = reduce_tracklets(optical)[1]
        # Admissible region
        A = AdmissibleRegion(tracklet, params)

        # Values by August 9, 2026

        # Zero AdmissibleRegion
        @test iszero(zero(AdmissibleRegion{Float64}))
        # Custom print
        @test string(A) == "AE: [116.61547, 45.39840, -3.21667, 5.76667] \
            t: 2024-01-20T21:50:15.360 obs: GINOP-KHK, Piszkesteto"
        # Coefficients
        @test length(A.coeffs) == 6
        @test A.coeffs[3] == A.vra^2 * cos(A.dec)^2 + A.vdec^2  # proper motion squared
        # Boundary functions
        xmin, xmax = A.ρ_domain
        @test arW(A, xmin) * arW(A, xmax) > 0
        @test ardW(A, xmin) * ardW(A, xmax) < 0
        @test ard2W(A, xmin) * ard2W(A, xmax) > 0
        @test arS(A, xmin) * arS(A, xmax) > 0
        @test ardS(A, xmin) * ardS(A, xmax) > 0
        @test ard2S(A, xmin) * ard2S(A, xmax) > 0
        @test ard2S(A, xmin) * ard2S(A, xmax) > 0
        @test arG(A, xmin) * arG(A, xmax) < 0
        # Energy discriminant
        @test arenergydis(A, A.ρ_domain[1], :outer) > 0
        @test arenergydis(A, A.ρ_domain[1], :inner) > 0
        @test arenergydis(A, A.ρ_domain[2], :outer) ≈ 0 atol = 1e-18
        ρ0 = min(R_SI, cbrt(2 * k_gauss^2 * μ_ES / A.coeffs[3]))
        @test arenergydis(A, ρ0, :inner) ≈ 0 atol = 1e-18
        @test arenergydis(A, A.ρ_domain[2] + 1.0, :outer) < 0
        @test arenergydis(A, ρ0 + 1.0, :inner) < 0
        # Range-rate
        @test rangerates(A, A.ρ_domain[1], :outer) == A.v_ρ_domain
        a, b = rangerates(A, A.ρ_domain[1], :inner)
        @test a ≈ -b atol = 1e-18
        @test minimum(rangerates(A, A.ρ_domain[2], :outer)) == A.Fs[3, 2]
        @test !isempty(rangerates(A, ρ0, :inner))
        @test rangerates(A, ρ0, :inner) == [zero(ρ0)]
        @test isempty(rangerates(A, A.ρ_domain[2] + 1.0, :outer))
        @test isempty(rangerates(A, ρ0 + 1.0, :inner))
        @test rangerate(A, A.ρ_domain[1], :min, :outer) == A.v_ρ_domain[1]
        @test rangerate(A, A.ρ_domain[1], :max, :outer) == A.v_ρ_domain[2]
        @test rangerate(A, A.ρ_domain[1], :min, :inner) ≈
            -rangerate(A, A.ρ_domain[1], :max, :inner) atol = 1e-18
        # Golden section search
        ρ, v_ρ = argoldensearch(A, A.ρ_domain..., :min, :outer, 1e-20)
        @test A.ρ_domain[1] ≤ ρ ≤ A.ρ_domain[2]
        @test v_ρ ≤ A.v_ρ_domain[1]
        ρ, v_ρ = argoldensearch(A, A.ρ_domain..., :max, :outer, 1e-20)
        @test A.ρ_domain[1] ≤ ρ ≤ A.ρ_domain[2]
        @test v_ρ ≥ A.v_ρ_domain[2]
        ρ, v_ρ = argoldensearch(A, A.ρ_domain[1], ρ0, :min, :inner, 1e-20)
        @test A.ρ_domain[1] ≤ ρ ≤ A.ρ_domain[2]
        @test A.v_ρ_domain[1] ≤ v_ρ ≤ A.v_ρ_domain[2]
        ρ, v_ρ = argoldensearch(A, A.ρ_domain[1], ρ0, :max, :inner, 1e-20)
        @test A.ρ_domain[1] ≤ ρ ≤ A.ρ_domain[2]
        @test A.v_ρ_domain[1] ≤ v_ρ ≤ A.v_ρ_domain[2]
        # Outer boundary
        O0 = arboundary(A, 0.0, :outer, :linear)
        O1 = arboundary(A, 1.0, :outer, :linear)
        O2 = arboundary(A, 2.0, :outer, :linear)
        O3 = arboundary(A, 3.0, :outer, :linear)
        @test O0[1] == O1[1] == A.ρ_domain[1]
        @test [O0[2], O1[2]] == A.v_ρ_domain
        @test O2[1] == _helmaxrange(A.coeffs, A.a_max) == A.ρ_domain[2]
        @test norm(O0 - O3) < 8e-18
        @test O0 == A.Fs[1, :]
        @test O1 == A.Fs[2, :]
        @test O2 == A.Fs[3, :]
        @test norm(O3 - A.Fs[1, :]) < 8e-18
        L0 = arboundary(A, 0.0, :outer, :log)
        L1 = arboundary(A, 1.0, :outer, :log)
        L2 = arboundary(A, 2.0, :outer, :log)
        L3 = arboundary(A, 3.0, :outer, :log)
        @test L0[1] == log10(O0[1])
        @test L1[1] == log10(O1[1])
        @test L2[1] == log10(O2[1])
        @test L3[1] ≈ log10(O3[1]) atol = 6e-15
        # Inner boundary
        I0 = arboundary(A, 0.0, :inner, :linear)
        I1 = arboundary(A, 1.0, :inner, :linear)
        I2 = arboundary(A, 2.0, :inner, :linear)
        @test I0[1] ≈ I2[1] atol = 1e-18
        @test I0[2] ≈ -I2[2] atol = 1e-18
        @test I1[1] ≈ ρ0 atol = 1e-18
        @test I1[2] ≈ 0.0 atol = 1e-10
        P0 = arboundary(A, 0.0, :inner, :log)
        P1 = arboundary(A, 1.0, :inner, :log)
        P2 = arboundary(A, 2.0, :inner, :log)
        @test P0[1] == P2[1] == log10(I0[1]) == log10(I2[1])
        @test P1[1] == log10(I1[1])
        # In
        @test A.Fs[1, :] in A
        @test A.Fs[2, :] in A
        @test A.Fs[3, :] in A
        @test [sum(A.ρ_domain), sum(A.v_ρ_domain)] / 2 in A
        # Topocentric to barycentric conversion
        @test norm(bary2topo(A, topo2bary(A, A.Fs[3, :]...)) .- A.Fs[3, :]) < 8e-6
        # Curvature
        w8s = Veres17(optical)
        C, Γ_C = curvature(optical, w8s)
        σ_C = sqrt.(diag(Γ_C))
        @test all( abs.(C) ./ σ_C .> 0.02)
        χ2 = C' * inv(Γ_C) * C
        @test χ2 > 2.7
        # Boundary projection
        xmin, xmax = A.ρ_domain
        ymin, ymax = A.v_ρ_domain
        xmid, ymid = (xmin + xmax) / 2, (ymin + ymax) / 2
        function distance(A, x, ρ, v_ρ)
            m = v_ρ > ymid ? :max : :min
            y = rangerate(A, x, m)
            return hypot(x - ρ, y - v_ρ)
        end
        ρs = [10^(-3.5), xmid, 1.0]
        v_ρs = [-0.02, ymid, 0.03]
        for (ρ, v_ρ) in Iterators.product(ρs, v_ρs)
            x, y = boundary_projection(A, ρ, v_ρ)
            if (ρ, v_ρ) in A
                @test x == ρ && y == v_ρ
            elseif ρ ≤ xmin
                @test x == xmin && ymin ≤ y ≤ ymax
            elseif ρ ≥ xmax
                @test x == xmax && y == ymid
            else
                @test distance(A, x, ρ, v_ρ) < distance(A, x - 1E-4, ρ, v_ρ)
                @test distance(A, x, ρ, v_ρ) < distance(A, x + 1E-4, ρ, v_ρ)
            end
        end
    end

    @testset "Too Short Arc" begin
        # Read optical astrometry
        optical = read_optical_mpc80(joinpath(TEST_DATA, "2008EK68.txt"))
        # Subset of optical for IOD
        # suboptical = iodsuboptical(optical, 3)
        # In this case: optical == suboptical

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           tsaorder = 2, adamiter = 500, adamQtol = 1e-5,
           jtlsorder = 2, jtlsiter = 20, lsiter = 10,
           significance = 0.99, outrej = false
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, optical)

        @test isa(string(od), String)
        @test scalartype(od) == Float64
        @test opticaltype(od) == OpticalMPC80{Float64}
        @test radartype(od) == Nothing && !hasradar(od)
        @test dof(od) == 6
        @test NEOs.optical(od) == optical && isnothing(NEOs.radar(od))

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by August 9, 2026

        # Curvature
        C, Γ_C = curvature(optical, od.weights)
        σ_C = sqrt.(diag(Γ_C))
        @test all( abs.(C) ./ σ_C .> 3.4)
        χ2 = C' * inv(Γ_C) * C
        @test χ2 > 1_006
        # Check type
        @test isa(orbit, OpticalOrbit{Float64})
        # Tracklets
        @test designation(orbit) == "2008 EK68"
        @test length(optical) == nobs(od) == nobs(orbit) == 10
        @test numberofdays(optical) == numberofdays(orbit) < 0.05
        @test minmaxdates(orbit) == (date(optical[1]), date(optical[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 1
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].indices[1] == 1
        @test orbit.tracklets[end].indices[end] == length(optical)
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), mean(r -> dtutc2days(date(r)), optical), atol = 7e-5)
        @test dtutc2days(date(optical[1])) > firsttime(orbit)
        @test all( norm.(orbit.bwd.p, Inf) .< 2 )
        @test dtutc2days(date(optical[end])) < lasttime(orbit)
        @test all( norm.(orbit.fwd.p, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 10
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test isa(string(orbit.fit), String)
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 6e-3 )
        @test all( snr(orbit) .> 4.1)
        @test chi2(orbit) < 14.24
        @test nrms(orbit) < 0.85
        # Covariance matrix
        Γ = covariance(orbit)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL_CAR = [-9.698411378106104E-01, 2.403528398096335E-01, 1.028828398440475E-01,
            -9.512276508152102E-03, -1.532548560550362E-02, -8.094637971143012E-03]
        JPL_GEO = [-1.068488198434094E-02, 2.545558390227619E-03, -1.475591869336636E-04,
            -4.812208060869299E-03, -3.675325693678584E-06, -1.451573958100394E-03]
        JPL_KEP = [1.484444069122194E+00, 3.966972238396319E-01, 3.961849525038958E+00,
            1.294946899344538E+02, 3.442349857158359E+02, 2.206177172198979E+01]
        JPL_EQN = keplerian2equinoctial(JPL_KEP, epoch(orbit) + MJD2000; μ = μ_S)
        JPL_ATTR = cartesian2attributable(JPL_GEO)
        jpl_compatibility_tests(orbit, params, (1.3E-02, 1.7E-02, 1.9E-11, 3.8E-13, 3.8E-13),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 29.6 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit, params)
        @test 3.1 < Da < Dc < Db < 7.3
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit, params)
        @test 4.3E4 < Ma < Mc < Mb < 5.2E5
        # MPEC
        @test isnothing(print_mpec(orbit, params))
        println()
    end

    @testset "Outlier Rejection" begin
        # Read optical astrometry
        optical = read_optical_mpc80(joinpath(TEST_DATA, "2007VV7.txt"))
        # Subset of optical for IOD
        suboptical = iodsuboptical(optical, 3)

        # Parameters
        params = Parameters(
            bwdoffset = 0.007, fwdoffset = 0.007,
            gaussorder = 2, jtlsorder = 2, jtlsiter = 20, lsiter = 10,
            outrej = true, χ2_rec = 1.0, χ2_rej = 1.25, fudge = 0.0
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, suboptical)

        @test isa(string(od), String)
        @test scalartype(od) == Float64
        @test opticaltype(od) == OpticalMPC80{Float64}
        @test radartype(od) == Nothing && !hasradar(od)
        @test dof(od) == 6
        @test NEOs.optical(od) == suboptical && isnothing(NEOs.radar(od))

        # Initial Orbit Determination (with outlier rejection)
        orbit = initialorbitdetermination(od, params)

        # Values by August 9, 2026

        # Check type
        @test isa(orbit, OpticalOrbit{Float64})
        # Tracklets
        @test designation(orbit) == "2007 VV7"
        @test length(suboptical) == nobs(od) == nobs(orbit) == 18
        @test numberofdays(suboptical) == numberofdays(orbit) < 2.16
        @test minmaxdates(orbit) == (date(suboptical[1]), date(suboptical[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 3
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].indices[1] == 1
        @test orbit.tracklets[end].indices[end] == length(suboptical)
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), dtutc2days(date(od.tracklets[2])), atol = 2e-3)
        @test dtutc2days(date(suboptical[1])) > firsttime(orbit)
        @test all( norm.(orbit.bwd.p, Inf) .< 2 )
        @test dtutc2days(date(suboptical[end])) < lasttime(orbit)
        @test all( norm.(orbit.fwd.p, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 16
        @test nout(orbit.ores) == 2
        # Least squares fit
        @test isa(string(orbit.fit), String)
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 4e-3 )
        @test all( snr(orbit) .> 50)
        @test chi2(orbit) < 1.55
        @test nrms(orbit) < 0.22
        # Covariance matrix
        Γ = covariance(orbit)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 7
        # @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL_CAR = [7.673359267780500E-01, 6.484889639435911E-01, 2.932326837353227E-01,
            -1.102334378117082E-02, 1.539269707113238E-02, 6.528842011032119E-03]
        JPL_GEO = [3.447012800435960E-02, 3.060747693168282E-02, 2.544634286791332E-02,
            8.618815060243545E-04, 3.787188937089640E-03, 1.496710946113642E-03]
        JPL_KEP = [1.776244846691859E+00, 4.381984418639090E-01, 7.819612775042287E-01,
            9.751283439586027E+01, 2.742918197067644E+02, 1.116208224849003E+01]
        JPL_EQN = keplerian2equinoctial(JPL_KEP, epoch(orbit) + MJD2000; μ = μ_S)
        JPL_ATTR = cartesian2attributable(JPL_GEO)
        jpl_compatibility_tests(orbit, params, (7.6E-02, 1.2E-01, 9.0E-12, 1.9E-14, 2.0E-12),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 26.7 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit, params)
        @test 11.9 < Da < Dc < Db < 27.1
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit, params)
        @test 2.2E6 < Ma < Mc < Mb < 2.8E7
        # MPEC
        @test isnothing(print_mpec(orbit, params))
        println()

        # Add remaining observations
        NEOs.update!(od, optical)
        # Refine orbit (with outlier rejection)
        orbit1 = orbitdetermination(od, orbit, params)

        # Check type
        @test isa(orbit1, OpticalOrbit{Float64})
        # Tracklets
        @test designation(orbit1) == "2007 VV7"
        @test length(optical) == nobs(od) == nobs(orbit1) == 21
        @test numberofdays(optical) == numberofdays(orbit1) < 3.03
        @test minmaxdates(orbit1) == (date(optical[1]), date(optical[end]))
        @test length(od.tracklets) == length(orbit1.tracklets) == 4
        @test od.tracklets == orbit1.tracklets
        @test orbit1.tracklets[1].indices[1] == 1
        @test orbit1.tracklets[end].indices[end] == length(optical)
        @test issorted(orbit1.tracklets)
        # Backward (forward) integration
        @test epoch(orbit1) == epoch(orbit)
        @test dtutc2days(date(optical[1])) > firsttime(orbit1)
        @test all( norm.(orbit1.bwd.p, Inf) .< 2 )
        @test dtutc2days(date(optical[end])) < lasttime(orbit1)
        @test all( norm.(orbit1.fwd.p, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit1.ores) == 19
        @test nout(orbit1.ores) == 2
        # Least squares fit
        @test isa(string(orbit1.fit), String)
        @test orbit1.fit.success
        @test all( sigmas(orbit1) .< 3e-4 )
        @test all( snr(orbit1) .> 574)
        @test chi2(orbit1) < 2.38
        @test nrms(orbit1) < 0.25
        # Covariance matrix
        Γ = covariance(orbit1)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit1.qs, 1) == 6
        # (26/04/2025) There are roundoff differences in the nrms of the two
        # jtls iterations; hence, in some os/julia versions, the first (second)
        # iteration has the lowest nrms.
        # @test size(orbit1.qs, 2) == length(orbit1.Qs) == 1
        @test issorted(orbit1.Qs, rev = true)
        @test orbit1.Qs[end] == nrms(orbit1)
        # Compatibility with JPL
        jpl_compatibility_tests(orbit1, params, (4.0E-03, 1.2E-02, 1.4E-10, 8.2E-11, 1.0E-10),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit1, params)
        @test H - dH ≤ 26.7 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit1, params) == 9
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit1, params)
        @test 11.8 < Da < Dc < Db < 26.8
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit1, params)
        @test 2.2E6 < Ma < Mc < Mb < 2.7E7
        # MPEC
        @test isnothing(print_mpec(orbit1, params))
        println()
    end

    @testset "Interesting NEOs" begin

        # 2014 AA hit the Earth around January 2, 2014, 02:49 UTC

        # Read optical astrometry
        optical = read_optical_mpc80(joinpath(TEST_DATA, "2014AA.txt"))
        # Subset of optical for IOD
        # suboptical = iodsuboptical(optical, 3)
        # In this case: optical == suboptical

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           tsaorder = 2, adamiter = 500, adamQtol = 1e-5,
           jtlsorder = 2, jtlsiter = 20, lsiter = 10,
           significance = 0.99, outrej = false
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, optical)

        @test isa(string(od), String)
        @test scalartype(od) == Float64
        @test opticaltype(od) == OpticalMPC80{Float64}
        @test radartype(od) == Nothing && !hasradar(od)
        @test dof(od) == 6
        @test NEOs.optical(od) == optical && isnothing(NEOs.radar(od))

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by August 9, 2026

        # Curvature
        C, Γ_C = curvature(optical, od.weights)
        σ_C = sqrt.(diag(Γ_C))
        @test all( abs.(C) ./ σ_C .> 5.7)
        χ2 = C' * inv(Γ_C) * C
        @test χ2 > 5.93e5
        # Check type
        @test isa(orbit, OpticalOrbit{Float64})
        # Tracklets
        @test designation(orbit) == "2014 AA"
        @test length(optical) == nobs(od) == nobs(orbit) == 7
        @test numberofdays(optical) == numberofdays(orbit) < 0.05
        @test minmaxdates(orbit) == (date(optical[1]), date(optical[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 1
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].indices[1] == 1
        @test orbit.tracklets[end].indices[end] == length(optical)
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), mean(r -> dtutc2days(date(r)), optical), atol = 2e-5)
        @test dtutc2days(date(optical[1])) > firsttime(orbit)
        @test all( norm.(orbit.bwd.p, Inf) .< 2 )
        @test dtutc2days(date(optical[end])) < lasttime(orbit)
        @test all( norm.(orbit.fwd.p, Inf) .< 1e9 )
        # Vector of residuals
        @test notout(orbit.ores) == 7
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test isa(string(orbit.fit), String)
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 3e-4 )
        @test all( snr(orbit) .> 20.5)
        @test chi2(orbit) < 0.23
        @test nrms(orbit) < 0.13
        # Covariance matrix
        Γ = covariance(orbit)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL_CAR = [-1.793429069495696E-01, 8.874118618831563E-01, 3.841433972166043E-01,
            -1.755785111623553E-02, -5.781634216294369E-03, -2.007510615561216E-03]
        JPL_GEO = [3.100782033379628E-04, 2.595364571202023E-03, 6.648062725412545E-04,
            -3.646375932548784E-04, -2.824074490602325E-03, -7.248940385846410E-04]
        JPL_KEP = [1.163575955666616E+00, 2.128185264087166E-01, 1.423597471953649E+00,
            5.237781301766019E+01, 1.016022028285875E+02, 3.243429036265208E+02]
        JPL_EQN = keplerian2equinoctial(JPL_KEP, epoch(orbit) + MJD2000; μ = μ_S)
        JPL_ATTR = cartesian2attributable(JPL_GEO)
        jpl_compatibility_tests(orbit, params, (3.0E-01, 3.8E-01, 1.5E-10, 7.9E-13, 3.2E-12),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 30.9 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit, params)
        @test 1.6 < Da < Dc < Db < 4.0
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit, params)
        @test 7.0E3 < Ma < Mc < Mb < 8.1E4
        # MPEC
        @test isnothing(print_mpec(orbit, params))
        println()

        # 2008 TC3 entered the Earth's atmosphere around October 7, 2008, 02:46 UTC

        # Read optical astrometry
        optical = read_optical_mpc80(joinpath(TEST_DATA, "2008TC3.txt"))
        # Subset of optical for IOD
        suboptical = iodsuboptical(optical, 3)

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           gaussorder = 2, jtlsorder = 2, jtlsiter = 20, lsiter = 10,
           significance = 0.99, outrej = false, parse_eqs = false
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, suboptical)

        @test isa(string(od), String)
        @test scalartype(od) == Float64
        @test opticaltype(od) == OpticalMPC80{Float64}
        @test radartype(od) == Nothing && !hasradar(od)
        @test dof(od) == 6
        @test NEOs.optical(od) == suboptical && isnothing(NEOs.radar(od))

        @test isa(string(od), String)
        @test scalartype(od) == Float64
        @test opticaltype(od) == OpticalMPC80{Float64}
        @test radartype(od) == Nothing && !hasradar(od)
        @test dof(od) == 6
        @test NEOs.optical(od) == suboptical && isnothing(NEOs.radar(od))

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by August 9, 2026

        # Check type
        @test isa(orbit, OpticalOrbit{Float64})
        # Tracklets
        @test designation(orbit) == "2008 TC3"
        @test length(suboptical) == nobs(od) == nobs(orbit) == 18
        @test numberofdays(suboptical) == numberofdays(orbit) < 0.34
        @test minmaxdates(orbit) == (date(suboptical[1]), date(suboptical[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 3
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].indices[1] == 1
        @test orbit.tracklets[end].indices[end] == length(suboptical)
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), dtutc2days(date(od.tracklets[2])), atol = 4e-3)
        @test dtutc2days(date(suboptical[1])) > firsttime(orbit)
        @test all( norm.(orbit.bwd.p, Inf) .< 2 )
        @test dtutc2days(date(suboptical[end])) < lasttime(orbit)
        @test all( norm.(orbit.fwd.p, Inf) .< 1e4 )
        # Vector of residuals
        @test notout(orbit.ores) == 18
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test isa(string(orbit.fit), String)
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 2e-5 )
        @test all( snr(orbit) .> 644)
        @test chi2(orbit) < 4.35
        @test nrms(orbit) < 0.35
        # Covariance matrix
        Γ = covariance(orbit)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL_CAR = [9.739753595124904E-01, 2.154167298799998E-01, 9.401075962536619E-02,
            -7.896756748173077E-03, 1.606197827206136E-02, 6.135361399316373E-03]
        JPL_GEO = [2.933267343256818E-03, -5.413518151214789E-04, 4.289238206468831E-04,
            -3.634877734669857E-03, 7.643113463866795E-04, -4.963086752649505E-04]
        JPL_KEP = [1.273091758414584E+00, 2.870222798582721E-01, 2.341999526552296E+00,
            2.339645303327229E+02, 1.941265709953888E+02, 3.288450951861228E+02]
        JPL_EQN = keplerian2equinoctial(JPL_KEP, epoch(orbit) + MJD2000; μ = μ_S)
        JPL_ATTR = cartesian2attributable(JPL_GEO)
        jpl_compatibility_tests(orbit, params, (5.0E-01, 3.0E-01, 1.2E-08, 1.5E-11, 1.5E-09),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 30.4 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 8
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit, params)
        @test 2.0 < Da < Dc < Db < 5.0
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit, params)
        @test 1.3E4 < Ma < Mc < Mb < 1.7E5
        # MPEC
        @test isnothing(print_mpec(orbit, params))
        println()

        # Add observations
        suboptical = iodsuboptical(optical, 10)

        # Refine orbit
        NEOs.update!(od, suboptical)
        orbit1 = orbitdetermination(od, orbit, params)

        # Check type
        @test isa(orbit1, OpticalOrbit{Float64})
        # Tracklets
        @test designation(orbit1) == "2008 TC3"
        @test length(suboptical) == nobs(od) == nobs(orbit1) == 97
        @test numberofdays(suboptical) == numberofdays(orbit1) < 0.70
        @test minmaxdates(orbit1) == (date(suboptical[1]), date(suboptical[end]))
        @test length(od.tracklets) == length(orbit1.tracklets) == 10
        @test od.tracklets == orbit1.tracklets
        @test orbit1.tracklets[1].indices[1] == 1
        @test orbit1.tracklets[end].indices[end] == 93
        @test issorted(orbit1.tracklets)
        # Backward (forward) integration
        @test epoch(orbit1) == epoch(orbit)
        @test dtutc2days(date(suboptical[1])) > firsttime(orbit1)
        @test all( norm.(orbit1.bwd.p, Inf) .< 1 )
        @test dtutc2days(date(suboptical[end])) < lasttime(orbit1)
        @test all( norm.(orbit1.fwd.p, Inf) .< 1e15 )
        # Vector of residuals
        @test notout(orbit1.ores) == 97
        @test nout(orbit1.ores) == 0
        # Least squares fit
        @test isa(string(orbit1.fit), String)
        @test orbit1.fit.success
        @test all( sigmas(orbit1) .< 4e-7 )
        @test all( snr(orbit1) .> 21_880)
        @test chi2(orbit1) < 54.85
        @test nrms(orbit1) < 0.53
        # Covariance matrix
        Γ = covariance(orbit1)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit1.qs, 1) == 6
        @test size(orbit1.qs, 2) == length(orbit1.Qs) == 2
        @test issorted(orbit1.Qs, rev = true)
        @test orbit1.Qs[end] == nrms(orbit1)
        # Compatibility with JPL
        jpl_compatibility_tests(orbit1, params, (1.1E+01, 5.0E-01, 2.2E-07, 1.8E-8, 1.4E-09),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit1, params)
        @test H - dH ≤ 30.4 ≤ H + dH
        # Parameters uncertainty
        @test all(sigmas(orbit1) .< sigmas(orbit))
        # TODO: understand better differences wrt JPL solutions
        # @test nrms(orbit1) < nrms(orbit)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit1, params) == 5
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit1, params)
        @test 2.1 < Da < Dc < Db < 5.0
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit1, params)
        @test 1.4E4 < Ma < Mc < Mb < 1.7E5
        # MPEC
        @test isnothing(print_mpec(orbit1, params))
        println()
    end

    @testset "research/2025CMDA/orbitdetermination.jl" begin
        using NEOs: iodinitcond, issatellite, hascoord

        function opticalfilter(optical::AbstractOpticalVector)
            # Eliminate observations before oficial discovery
            firstobs = findfirst(r -> !isempty(r.discovery), optical)
            isnothing(firstobs) && return false, optical
            optical = optical[firstobs:end]
            # Filter out incompatible observations
            filter!(optical) do r
                hascoord(r.observatory) && !issatellite(r.observatory) &&
                date(r) > DateTime(2000, 1, 1, 12)
            end
            length(optical) < 3 && return false, optical
            # Find the first set of 3 tracklets with a < 15 days timespan
            tracklets = reduce_tracklets(optical)
            for i in 1:length(tracklets)-2
                numberofdays(tracklets[i:i+2]) > 15.0 && continue
                tracklets = tracklets[i:i+2]
                optical = optical[indices(tracklets)]
                sort!(optical)
                break
            end
            return numberofdays(optical) <= 15.0, optical
        end

        # Fetch and filter optical astrometry
        optical = read_optical_mpc80(joinpath(TEST_DATA, "2023QR6.txt"))
        flag, optical = opticalfilter(optical)

        # Parameters
        params = Parameters(
            coeffstol = Inf, bwdoffset = 0.042, fwdoffset = 0.042,          # Propagation
            gaussorder = 2, safegauss = true, refscale = :log,              # Gauss method
            tsaorder = 2, adamiter = 500, adamQtol = 1e-5,                  # ADAM
            jtlsorder = 2, jtlsiter = 20, lsiter = 20, significance = 0.99, # Least squares
            outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,                      # Outlier rejection
            fudge = 100.0, max_per = 20.0
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, optical)

        @test isa(string(od), String)
        @test scalartype(od) == Float64
        @test opticaltype(od) == OpticalMPC80{Float64}
        @test radartype(od) == Nothing && !hasradar(od)
        @test dof(od) == 6
        @test NEOs.optical(od) == optical && isnothing(NEOs.radar(od))

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params; initcond = iodinitcond)

        # Values by August 9, 2026

        # Check type
        @test isa(orbit, OpticalOrbit{Float64})
        # Tracklets
        @test flag
        @test designation(orbit) == "2023 QR6"
        @test length(optical) == nobs(od) == nobs(orbit) == 6
        @test numberofdays(optical) == numberofdays(orbit) < 6.22
        @test minmaxdates(orbit) == (date(optical[1]), date(optical[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 3
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].indices[1] == 1
        @test orbit.tracklets[end].indices[end] == length(optical)
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), dtutc2days(date(od.tracklets[2])), atol = 3e-3)
        @test dtutc2days(date(optical[1])) > firsttime(orbit)
        @test all( norm.(orbit.bwd.p, Inf) .< 2 )
        @test dtutc2days(date(optical[end])) < lasttime(orbit)
        @test all( norm.(orbit.fwd.p, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 6
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test isa(string(orbit.fit), String)
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 0.018 )
        @test all( snr(orbit) .> 7.00)
        @test chi2(orbit) < 0.91
        @test nrms(orbit) < 0.28
        # Covariance matrix
        Γ = covariance(orbit)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 5
        # @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL_CAR = [8.272632404899825E-01, -8.060647872503925E-01, -6.506192873607640E-01,
            1.660015292864715E-02, -5.614753953372325E-03, 2.899476523012778E-03]
        JPL_GEO = [-4.175116355294079E-02, -3.431743572536730E-01, -4.501902423837642E-01,
            8.333024069432302E-03, -1.924883145539331E-02, -3.011265449368804E-03]
        JPL_KEP = [2.279881975143154E+00, 7.595854949078547E-01, 3.859999487430863E+01,
            2.293324466043271E+02, 3.254353160044697E+02, 2.033667388008325E+01]
        JPL_EQN = keplerian2equinoctial(JPL_KEP, epoch(orbit) + MJD2000; μ = μ_S)
        JPL_ATTR = cartesian2attributable(JPL_GEO)
        jpl_compatibility_tests(orbit, params, (6.0E-01, 1.4E+00, 2.1E-12, 5.9E-14, 5.9E-14),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 18.5 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit, params)
        @test 558.7 < Da < Dc < Db < 1250.1
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit, params)
        @test 2.2E11 < Ma < Mc < Mb < 2.8E12
        # MPEC
        @test isnothing(print_mpec(orbit, params))
        println()
    end

    @testset "Radar astrometry" begin

        # Load observations
        optical = read_optical_mpc80(joinpath(pkgdir(NEOs), "data",
            "99942_2004_2020.dat"))
        filter!(x -> Date(2005, 1, 27) < date(x) < Date(2005, 1, 31), optical)
        radar = read_radar_jpl(joinpath(pkgdir(NEOs, "test", "data",
            "99942_RADAR_2005_2013.json")))
        filter!(x -> Date(2005, 1, 27) < date(x) < Date(2005, 1, 31), radar)

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           gaussorder = 2, safegauss = false,
           tsaorder = 2, adamiter = 500, adamQtol = 1e-5, jtlsorder = 2,
           jtlsmask = false, jtlsiter = 20, lsiter = 10, significance = 0.99,
           outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,
           fudge = 100.0, max_per = 34.0,
       )
        # Orbit determination problem (only optical astrometry)
        od0 = ODProblem(newtonian!, optical)

        @test isa(string(od0), String)
        @test scalartype(od0) == Float64
        @test opticaltype(od0) == OpticalMPC80{Float64}
        @test radartype(od0) == Nothing && !hasradar(od0)
        @test dof(od0) == 6
        @test NEOs.optical(od0) == optical && isnothing(NEOs.radar(od0))

        # Preliminary orbit (only optical astrometry)
        loadjpleph()
        jd0 = datetime2julian(DateTime(2005, 1, 29))
        q00 = kmsec2auday(apophisposvel199(julian2etsecs(jd0)))
        orbit0 = LeastSquaresOrbit(od0, q00, jd0, params)

        # Shift reference epoch
        _jd0_ = mean(dtutc2days, minmaxdates(orbit0)) + PE.J2000
        _orbit0_ = shiftepoch(orbit0, _jd0_, params)
        t0, _t0_ = epoch(orbit0), epoch(_orbit0_)
        @test mre(orbit0(t0), _orbit0_(t0), sigmas(orbit0)) < eps()
        @test mre(orbit0(_t0_), _orbit0_(_t0_), sigmas(_orbit0_)) < eps()

        # Orbit determination problem (both optical and radar astrometry)
        od1 = ODProblem(newtonian!, optical, radar)

        @test isa(string(od1), String)
        @test scalartype(od1) == Float64
        @test opticaltype(od1) == OpticalMPC80{Float64}
        @test radartype(od1) == RadarJPL{Float64} && hasradar(od1)
        @test dof(od1) == 6
        @test NEOs.optical(od1) == optical && NEOs.radar(od1) == radar

        # Refine orbit (both optical and radar astrometry)
        orbit1 = orbitdetermination(od1, orbit0, params)

        # Values by August 9, 2026

        # Check type
        @test isa(orbit1, RadarOrbit{Float64})
        # Astrometry
        @test designation(orbit1) == "99942"
        @test length(optical) == noptical(od1) == noptical(orbit1) == 24
        @test length(radar) == nradar(od1) == nradar(orbit1) == 5
        @test length(optical) + length(radar) == nobs(od1) == nobs(orbit1) == 24 + 5
        @test numberofdays(radar) == numberofdays(orbit1) < 2.04
        @test minmaxdates(orbit1) == (date(radar[1]), date(radar[end]))
        @test length(od1.tracklets) == length(orbit1.tracklets) == 7
        @test od1.tracklets == orbit1.tracklets
        @test od1.radar == orbit1.radar
        @test orbit1.tracklets[1].indices[1] == 1
        @test orbit1.tracklets[end].indices[end] == 23
        @test issorted(orbit1.tracklets)
        @test issorted(orbit1.radar)
        # Backward (forward) integration
        @test epoch(orbit1) == datetime2julian(date(radar[2])) - PE.J2000
        @test dtutc2days(date(optical[1])) > firsttime(orbit1)
        @test dtutc2days(date(radar[1])) > firsttime(orbit1)
        @test all( norm.(orbit1.bwd.p, Inf) .< 2 )
        @test dtutc2days(date(optical[end])) < lasttime(orbit1)
        @test dtutc2days(date(radar[end])) < lasttime(orbit1)
        @test all( norm.(orbit1.fwd.p, Inf) .< 2 )
        # Vectors of residuals
        @test notout(orbit1.ores) == 24
        @test notout(orbit1.rres) == 5
        @test nout(orbit1.ores) == 0
        @test nout(orbit1.rres) == 0
        # Least squares fit
        @test isa(string(orbit1.fit), String)
        @test orbit1.fit.success
        @test all( sigmas(orbit1) .< 2.9e-7 )
        @test all( snr(orbit1) .> 8_342)
        @test chi2(orbit1) < 19.7
        @test nrms(orbit1) < 0.61
        # Covariance matrix
        Γ = covariance(orbit1)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit1.qs, 1) == 6
        @test size(orbit1.qs, 2) == length(orbit1.Qs) <= 2
        @test issorted(orbit1.Qs, rev = true)
        @test orbit1.Qs[end] == nrms(orbit1)
        # Compatibility with JPL
        JPL_CAR = [-5.229992130937651E-01, 8.689454573480734E-01, 3.096174868699621E-01,
            -1.413639580483663E-02, -5.510379552549767E-03, -2.413003153288419E-03]
        JPL_GEO = [9.403261433284629E-02, 1.678184682535047E-01, 5.767627532379624E-03,
            -5.089377531997436E-04, 4.493620961259394E-03, 1.923654593847778E-03]
        JPL_KEP = [9.223295977030230E-01, 1.911190963976789E-01, 3.330797253820763E+00,
            1.263744754026979E+02, 2.044798558304837E+02, 1.361032871672047E+02]
        JPL_EQN = keplerian2equinoctial(JPL_KEP, epoch(orbit1) + MJD2000; μ = μ_S)
        JPL_ATTR = cartesian2attributable(JPL_GEO)
        jpl_compatibility_tests(orbit1, params, (4.5E+00, 1.9E+02, 4.0E-08, 1.5E-11, 8.0E-10),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit1, params)
        @test H - dH ≤ 18.93 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit1, params) == 5
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit1, params)
        @test 433.3 < Da < Dc < Db < 969.8
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit1, params)
        @test 1.0E11 < Ma < Mc < Mb < 1.4E12
        # MPEC
        @test isnothing(print_mpec(orbit1, params))
        println()
    end

    @testset "Linkage" begin
        # Read optical astrometry
        optical = read_optical_mpc80(joinpath(TEST_DATA, "2020NB1.txt"))

        # Parameters
        params = Parameters(
            maxsteps = 20_000, order = 15, abstol = 1E-12, parse_eqs = true,
            coeffstol = Inf, bwdoffset = 0.05, fwdoffset = 0.05,
            gaussorder = 2, safegauss = false, refscale = :log,
            tsaorder = 2, adamiter = 500, adamQtol = 1E-5,
            jtlsorder = 2, jtlsmask = false, jtlsiter = 20, lsiter = 10,
            jtlsproject = true, significance = 0.99, verbose = true,
            outrej = true, χ2_rec = sqrt(9.21), χ2_rej = sqrt(10), fudge = 100.0,
            max_per = 33.3
        )

        # Orbit determination problem with only the observations from 2020
        idxs = findall(x -> date(x) > Date(2020), optical)
        od = ODProblem(newtonian!, optical[idxs]; weights = Veres17, debias = Eggl20)

        @test isa(string(od), String)
        @test scalartype(od) == Float64
        @test opticaltype(od) == OpticalMPC80{Float64}
        @test radartype(od) == Nothing && !hasradar(od)
        @test dof(od) == 6
        @test NEOs.optical(od) == optical[idxs] && isnothing(NEOs.radar(od))

        # Orbit Determination with only the observations from 2020
        orbit = initialorbitdetermination(od, params)

        # Include the observations from 2010
        NEOs.update!(od, optical)
        # Linkage
        orbit = linkage(od, orbit, params)

        # Values by August 9, 2026

        # Check type
        @test isa(orbit, OpticalOrbit{Float64})
        # Tracklets
        @test designation(orbit) == "2020 NB1"
        @test length(optical) == nobs(od) == nobs(orbit) == 44
        @test numberofdays(optical) == numberofdays(orbit) < 3766
        @test minmaxdates(orbit) == (date(optical[1]), date(optical[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 18
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].indices[1] == 1
        @test orbit.tracklets[end].indices[end] == length(optical)
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), mean(r -> dtutc2days(date(r)), optical[idxs]), atol = 6E+0)
        @test dtutc2days(date(optical[1])) > firsttime(orbit)
        @test all( norm.(orbit.bwd.p, Inf) .< 3.7 )
        @test dtutc2days(date(optical[end])) < lasttime(orbit)
        @test all( norm.(orbit.fwd.p, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 43
        @test nout(orbit.ores) == 1
        # Least squares fit
        @test isa(string(orbit.fit), String)
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 1.4E-6 )
        @test all( snr(orbit) .> 1E+5)
        @test chi2(orbit) < 10.40
        @test nrms(orbit) < 0.35
        # Covariance matrix
        Γ = covariance(orbit)
        @test size(Γ) == (6, 6)
        @test eigmin(Γ) > 0
        @test isapprox(Γ, Γ')
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL_CAR = [4.848426285691420E-01, -1.257122809354289E+00, 7.554427206018359E-02,
            2.838166871048631E-03, 1.669151616500038E-02, 4.630710656278403E-03]
        JPL_GEO = [-1.090443706381998E-01, -5.114340419413202E-01, 3.987237578405338E-01,
            -1.075454656354707E-02, 7.432601214257392E-03, 6.161736954055316E-04]
        JPL_KEP = [2.315552921489858E+00, 8.475227268358573E-01, 3.315750659948262E+01,
            1.778735224777621E+02, 2.484374014196204E+02, 3.416088210735048E+02]
        JPL_EQN = keplerian2equinoctial(JPL_KEP, epoch(orbit) + MJD2000; μ = μ_S)
        JPL_ATTR = cartesian2attributable(JPL_GEO)
        jpl_compatibility_tests(orbit, params, (2.8E+0, 3.8E+0, 1.8E-7, 5.0E-9, 1.7E-8),
                                JPL_CAR, JPL_KEP, JPL_EQN, JPL_ATTR)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 20.59 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 1
        # Diameter
        Da, Db = minmax(diameter(H, 0.05), diameter(H, 0.25))
        Dc = diameter(orbit, params)
        @test 200 < Da < Dc < Db < 448
        # Mass
        Ma, Mb = minmax(mass(2_600, Da), mass(2_600, Db))
        Mc = mass(orbit, params)
        @test 1.0E+10 < Ma < Mc < Mb < 1.3E+11
        # MPEC
        @test isnothing(print_mpec(orbit, params))
        println()
    end

end
