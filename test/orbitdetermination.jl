# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using PlanetaryEphemeris
using LinearAlgebra
using Test

using NEOs: AbstractOpticalVector, μ_S, indices, equatorial2ecliptic
using Statistics: mean

function iodsuboptical(optical::AbstractOpticalVector, N::Int = 3)
    tracklets = reduce_tracklets(optical)
    idxs = indices(tracklets[1:N])
    suboptical = optical[idxs]
    return suboptical
end

@testset "Orbit Determination" begin
    @testset "Straight Gauss Method" begin
        # Load observations
        optical = read_optical_mpc80(joinpath(pkgdir(NEOs), "test", "data",
            "2023DW_OPTICAL.dat"))
        # Subset of optical for IOD
        suboptical = iodsuboptical(optical, 3)

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           gaussorder = 2, jtlsorder = 2, jtlsiter = 200, lsiter = 1,
           significance = 0.99, outrej = false, parse_eqs = false
        )
        params = Parameters(params, parse_eqs = true)
        # Orbit determination problem
        od = ODProblem(newtonian!, suboptical)

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by Sep 28, 2025

        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
            typeof(suboptical)})
        # Tracklets
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
        @test dtutc2days(date(suboptical[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(suboptical[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 9
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 9e-4 )
        @test all( snr(orbit) .> 14.5)
        @test chi2(orbit) < 0.53
        @test nrms(orbit) < 0.18
        # Jacobian
        @test size(orbit.jacobian) == (6, 6)
        @test isdiag(orbit.jacobian)
        @test maximum(orbit.jacobian) < 9e-4
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        q0, σ0 = orbit(), sigmas(orbit)
        JPL_CAR = [-0.9867704701732631, 0.3781890325424674, 0.14094513213009532,
            -0.008773157203087259, -0.00947109649687576, -0.005654229864757284]
        @test all(@. abs(q0 - JPL_CAR) / σ0 < 0.76)
        # Osculating orbital elements
        osc = osculating(orbit, params)
        @test iselliptic(osc)
        @test osc.mu == μ_S
        @test epoch(osc) == epoch(orbit) + MJD2000
        @test osc.frame == :ecliptic
        q0 = equatorial2ecliptic(orbit() - params.eph_su(epoch(orbit)))
        @test norm(q0 - osc(), Inf) < 4.0e-14
        q0, σ0 = elements(osc), sigmas(osc)
        JPL_OSC = [8.198835710815939E-01, 3.962989389356275E-01, 5.807184452352074E+00,
            4.045356567755751E+01, 3.261403954217091E+02, 1.216665112960870E+02]
        @test all(@. abs(q0 - JPL_OSC) / σ0 < 0.83)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 24.3 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(orbit, params, 0.05), diameter(orbit, params, 0.25))
        @test 34.3 < Da < Db < 76.9
        # Mass
        Ma, Mb = minmax(mass(orbit, params, 2_600, 0.05), mass(orbit, params, 2_600, 0.25))
        @test 5.5E7 < Ma < Mb < 6.2E8

        # Add observations
        suboptical = iodsuboptical(optical, 15)

        # Refine orbit
        NEOs.update!(od, suboptical)
        orbit1 = orbitdetermination(od, orbit, params)

        # Check type
        @test isa(orbit1, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
            typeof(suboptical)})
        # Tracklets
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
        @test dtutc2days(date(suboptical[1])) > orbit1.bwd.t0 + orbit1.bwd.t[end]
        @test all( norm.(orbit1.bwd.x, Inf) .< 1.2 )
        @test dtutc2days(date(suboptical[end])) < orbit1.fwd.t0 + orbit1.fwd.t[end]
        @test all( norm.(orbit1.fwd.x, Inf) .< 1.2 )
        # Vector of residuals
        @test notout(orbit1.ores) == 43
        @test nout(orbit1.ores) == 0
        # Least squares fit
        @test orbit1.fit.success
        @test all( sigmas(orbit1) .< 2e-4 )
        @test all( snr(orbit1) .> 866 )
        @test chi2(orbit1) < 11.64
        @test nrms(orbit1) < 0.37
        # Jacobian
        @test size(orbit1.jacobian) == (6, 6)
        @test isdiag(orbit1.jacobian)
        @test maximum(orbit1.jacobian) < 9e-4
        # Convergence history
        @test size(orbit1.qs, 1) == 6
        @test size(orbit1.qs, 2) == length(orbit1.Qs) == 6
        @test issorted(orbit1.Qs, rev = true)
        @test orbit1.Qs[end] == nrms(orbit1)
        # Compatibility with JPL
        q0, σ0 = orbit1(), sigmas(orbit1)
        @test all(@. abs(q0 - JPL_CAR) / σ0 < 0.31)
        # Osculating orbital elements
        osc = osculating(orbit1, params)
        @test iselliptic(osc)
        @test osc.mu == μ_S
        @test epoch(osc) == epoch(orbit1) + MJD2000
        @test osc.frame == :ecliptic
        q0 = equatorial2ecliptic(orbit1() - params.eph_su(epoch(orbit1)))
        @test norm(q0 - osc(), Inf) < 3.8e-14
        q0, σ0 = elements(osc), sigmas(osc)
        @test all(@. abs(q0 - JPL_OSC) / σ0 < 0.31)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit1, params)
        @test H - dH ≤ 24.3 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit1, params) == 7
        # Diameter
        Da, Db = minmax(diameter(orbit1, params, 0.05), diameter(orbit1, params, 0.25))
        @test 36.3 < Da < Db < 81.3
        # Mass
        Ma, Mb = minmax(mass(orbit1, params, 2_600, 0.05), mass(orbit1, params, 2_600, 0.25))
        @test 6.5E7 < Ma < Mb < 7.3E8
    end

    @testset "Unsafe Gauss Method" begin
        # Load observations
        optical = fetch_optical_mpc80("2005 TM173", MPC)
        # Subset of optical for IOD
        # suboptical = iodsuboptical(optical, 3)
        # In this case: optical == suboptical

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           gaussorder = 2, jtlsorder = 2, jtlsiter = 200, lsiter = 1,
           significance = 0.99, outrej = false, safegauss = false
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, optical)

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by Sep 28, 2025

        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
            typeof(optical)})
        # Tracklets
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
        @test dtutc2days(date(optical[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(optical[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 6
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 5e-3 )
        @test all( snr(orbit) .> 21.4)
        @test chi2(orbit) < 2.53
        @test nrms(orbit) < 0.46
        # Jacobian
        @test size(orbit.jacobian) == (6, 6)
        @test isdiag(orbit.jacobian)
        @test maximum(orbit.jacobian) < 6.8e-4
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        q0, σ0 = orbit(), sigmas(orbit)
        JPL_CAR = [1.0042569058151192, 0.2231639040146286, 0.11513854178693468,
            -0.010824212819531798, 0.017428798232689943, 0.0071046780555307385]
        @test all(@. abs(q0 - JPL_CAR) / σ0 < 5.6e-3)
        # Osculating orbital elements
        osc = osculating(orbit, params)
        @test iselliptic(osc)
        @test osc.mu == μ_S
        @test epoch(osc) == epoch(orbit) + MJD2000
        @test osc.frame == :ecliptic
        q0 = equatorial2ecliptic(orbit() - params.eph_su(epoch(orbit)))
        @test norm(q0 - osc(), Inf) < 6.1e-14
        q0, σ0 = elements(osc), sigmas(osc)
        JPL_OSC = [2.872424697642789E+00, 6.749395051551541E-01, 1.282355986214476E+00,
            1.725712172730245E+02, 2.413793589329106E+02, 3.538743602962668E+02 - 360]
        @test all(@. abs(q0 - JPL_OSC) / σ0 < 5.5e-3)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 24.0 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(orbit, params, 0.05), diameter(orbit, params, 0.25))
        @test 41.2 < Da < Db < 92.2
        # Mass
        Ma, Mb = minmax(mass(orbit, params, 2_600, 0.05), mass(orbit, params, 2_600, 0.25))
        @test 9.5E7 < Ma < Mb < 1.1E9
    end

    @testset "Gauss Method with ADAM refinement" begin
        # Load observations
        optical = fetch_optical_mpc80("2024 MK", MPC)
        # Subset of optical for IOD
        suboptical = optical[10:21]

        # Parameters
        params = Parameters(
            coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
            gaussorder = 2, safegauss = true,
            tsaorder = 2, adamiter = 500, adamQtol = 1e-5, jtlsorder = 6,
            jtlsmask = false, jtlsiter = 20, lsiter = 10, significance = 0.99,
            outrej = true, χ2_rec = sqrt(9.21), χ2_rej = sqrt(10),
            fudge = 100.0, max_per = 34.0,
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, suboptical)

        # Initial Orbit Determination
        orbit = gaussiod(od, params)

        # Values by Sep 28, 2025

        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
            typeof(suboptical)})
        # Tracklets
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
        @test dtutc2days(date(suboptical[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(suboptical[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 12
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 6.6e-4 )
        @test all( snr(orbit) .> 38.8)
        @test chi2(orbit) < 2.43
        @test nrms(orbit) < 0.32
        # Jacobian
        @test size(orbit.jacobian) == (6, 6)
        @test isdiag(orbit.jacobian)
        @test maximum(orbit.jacobian) < 9.5e-6
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 3
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        q0, σ0 = orbit(), sigmas(orbit)
        JPL_CAR = [-0.12722461679828806, -0.9466098076903212, -0.4526816007640767,
            0.02048875631534963, -0.00022720097573790754, 0.00321302850930331]
        @test all(@. abs(q0 - JPL_CAR) / σ0 < 0.16)
        # Osculating orbital elements
        osc = osculating(orbit, params)
        @test iselliptic(osc)
        @test osc.mu == μ_S
        @test epoch(osc) == epoch(orbit) + MJD2000
        @test osc.frame == :ecliptic
        q0 = equatorial2ecliptic(orbit() - params.eph_su(epoch(orbit)))
        @test norm(q0 - osc(), Inf) < 5.0e-14
        q0, σ0 = elements(osc), sigmas(osc)
        JPL_OSC = [2.232655272359251E+00, 5.480018648354085E-01, 8.456325272115306E+00,
            1.322293305412568E+01, 2.778787985886902E+02, 3.530120962605258E+02 - 360]
        @test all(@. abs(q0 - JPL_OSC) / σ0 < 0.17)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 21.7 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(orbit, params, 0.05), diameter(orbit, params, 0.25))
        @test 120.3 < Da < Db < 269.1
        # Mass
        Ma, Mb = minmax(mass(orbit, params, 2_600, 0.05), mass(orbit, params, 2_600, 0.25))
        @test 2.3E9 < Ma < Mb < 2.7E10
    end

    @testset "Admissible region" begin
        using NEOs: AdmissibleRegion, arenergydis, rangerate, rangerates,
            argoldensearch, arboundary, _helmaxrange, R_SI, k_gauss, μ_ES,
            topo2bary, bary2topo

        # Fetch optical astrometry
        optical = fetch_optical_mpc80("2024 BX1", MPC)
        # Parameters
        params = Parameters()
        # First tracklet
        optical = optical[1:3]
        tracklet = reduce_tracklets(optical)[1]
        # Admissible region
        A = AdmissibleRegion(tracklet, params)

        # Values by Sep 28, 2025

        # Zero AdmissibleRegion
        @test iszero(zero(AdmissibleRegion{Float64}))
        # Custom print
        @test string(A) == "AE: [116.61547, 45.39840, -3.21667, 5.76667] \
            t: 2024-01-20T21:50:15.360 obs: GINOP-KHK, Piszkesteto"
        # Coefficients
        @test length(A.coeffs) == 6
        @test A.coeffs[3] == A.vra^2 * cos(A.dec)^2 + A.vdec^2  # proper motion squared
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
    end

    @testset "Too Short Arc" begin
        # Fetch optical astrometry
        optical = fetch_optical_mpc80("2008 EK68", MPC)
        # Subset of optical for IOD
        # suboptical = iodsuboptical(optical, 3)
        # In this case: optical == suboptical

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

        # Values by Sep 28, 2025

        # Curvature
        C, Γ_C = curvature(optical, od.weights)
        σ_C = sqrt.(diag(Γ_C))
        @test all( abs.(C) ./ σ_C .> 3.4)
        χ2 = C' * inv(Γ_C) * C
        @test χ2 > 1_006
        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
            typeof(optical)})
        # Tracklets
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
        @test dtutc2days(date(optical[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(optical[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 10
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 6e-3 )
        @test all( snr(orbit) .> 4.1)
        @test chi2(orbit) < 14.24
        @test nrms(orbit) < 0.85
        # Jacobian
        @test size(orbit.jacobian) == (6, 6)
        @test isdiag(orbit.jacobian)
        @test maximum(orbit.jacobian) < 9.2e-3
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        q0, σ0 = orbit(), sigmas(orbit)
        JPL_CAR = [-0.9698405495747651, 0.24035304578776012, 0.10288276585828428,
            -0.009512301266159554, -0.01532548565855646, -0.00809464581680694]
        @test all(@. abs(q0 - JPL_CAR) / σ0 < 0.013)
        # Osculating orbital elements
        osc = osculating(orbit, params)
        @test iselliptic(osc)
        @test osc.mu == μ_S
        @test epoch(osc) == epoch(orbit) + MJD2000
        @test osc.frame == :ecliptic
        q0 = equatorial2ecliptic(orbit() - params.eph_su(epoch(orbit)))
        @test norm(q0 - osc(), Inf) < 5.2e-14
        q0, σ0 = elements(osc), sigmas(osc)
        JPL_OSC = [1.484448954296998E+00, 3.966995063832199E-01, 3.961868860866990E+00,
            1.294946774085543E+02, 3.442349856198504E+02, 2.206164242012555E+01]
        @test all(@. abs(q0 - JPL_OSC) / σ0 < 0.013)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 29.6 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(orbit, params, 0.05), diameter(orbit, params, 0.25))
        @test 3.2 < Da < Db < 7.2
        # Mass
        Ma, Mb = minmax(mass(orbit, params, 2_600, 0.05), mass(orbit, params, 2_600, 0.25))
        @test 4.4E4 < Ma < Mb < 5.0E5
    end

    @testset "Outlier Rejection" begin
        # Fetch optical astrometry
        optical = fetch_optical_mpc80("2007 VV7", MPC)
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

        # Initial Orbit Determination (with outlier rejection)
        orbit = initialorbitdetermination(od, params)

        # Values by Sep 28, 2025

        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
            typeof(suboptical)})
        # Tracklets
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
        @test dtutc2days(date(suboptical[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(suboptical[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 16
        @test nout(orbit.ores) == 2
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 4e-3 )
        @test all( snr(orbit) .> 50)
        @test chi2(orbit) < 1.55
        @test nrms(orbit) < 0.22
        # Jacobian
        @test size(orbit.jacobian) == (6, 6)
        @test isdiag(orbit.jacobian)
        @test maximum(orbit.jacobian) < 5e-4
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 7
        # @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        q0, σ0 = orbit(), sigmas(orbit)
        JPL_CAR = [0.7673366466815864, 0.6484892781853565, 0.29323267343908294,
            -0.011023343781911974, 0.015392697071667377, 0.006528842022004942]
        @test all(@. abs(q0 - JPL_CAR) / σ0 < 0.08)
        # Osculating orbital elements
        osc = osculating(orbit, params)
        @test iselliptic(osc)
        @test osc.mu == μ_S
        @test epoch(osc) == epoch(orbit) + MJD2000
        @test osc.frame == :ecliptic
        q0 = equatorial2ecliptic(orbit() - params.eph_su(epoch(orbit)))
        @test norm(q0 - osc(), Inf) < 3.6e-14
        q0, σ0 = elements(osc), sigmas(osc)
        JPL_OSC = [1.776244846691859E+00, 4.381984418639090E-01, 7.819612775042287E-01,
            9.751283439586027E+01, 2.742918197067644E+02, 1.116208224849003E+01]
        @test all(@. abs(q0 - JPL_OSC) / σ0 < 0.07)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 26.7 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(orbit, params, 0.05), diameter(orbit, params, 0.25))
        @test 12.0 < Da < Db < 27.0
        # Mass
        Ma, Mb = minmax(mass(orbit, params, 2_600, 0.05), mass(orbit, params, 2_600, 0.25))
        @test 2.3E6 < Ma < Mb < 2.7E7

        # Add remaining observations
        NEOs.update!(od, optical)
        # Refine orbit (with outlier rejection)
        orbit1 = orbitdetermination(od, orbit, params)

        # Check type
        @test isa(orbit1, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
            typeof(optical)})
        # Tracklets
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
        @test dtutc2days(date(optical[1])) > orbit1.bwd.t0 + orbit1.bwd.t[end]
        @test all( norm.(orbit1.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(optical[end])) < orbit1.fwd.t0 + orbit1.fwd.t[end]
        @test all( norm.(orbit1.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit1.ores) == 19
        @test nout(orbit1.ores) == 2
        # Least squares fit
        @test orbit1.fit.success
        @test all( sigmas(orbit1) .< 3e-4 )
        @test all( snr(orbit1) .> 574)
        @test chi2(orbit1) < 2.38
        @test nrms(orbit1) < 0.25
        # Jacobian
        @test size(orbit1.jacobian) == (6, 6)
        @test isdiag(orbit1.jacobian)
        @test maximum(orbit1.jacobian) < 4e-3
        # Convergence history
        @test size(orbit1.qs, 1) == 6
        # (26/04/2025) There are roundoff differences in the nrms of the two
        # jtls iterations; hence, in some os/julia versions, the first (second)
        # iteration has the lowest nrms.
        # @test size(orbit1.qs, 2) == length(orbit1.Qs) == 1
        @test issorted(orbit1.Qs, rev = true)
        @test orbit1.Qs[end] == nrms(orbit1)
        # Compatibility with JPL
        q0, σ0 = orbit1(), sigmas(orbit1)
        @test all(@. abs(q0 - JPL_CAR) / σ0 < 6.9e-4)
        # Osculating orbital elements
        osc = osculating(orbit1, params)
        @test iselliptic(osc)
        @test osc.mu == μ_S
        @test epoch(osc) == epoch(orbit1) + MJD2000
        @test osc.frame == :ecliptic
        q0 = equatorial2ecliptic(orbit1() - params.eph_su(epoch(orbit1)))
        @test norm(q0 - osc(), Inf) < 4.0e-14
        q0, σ0 = elements(osc), sigmas(osc)
        @test all(@. abs(q0 - JPL_OSC) / σ0 < 6.9e-4)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 26.7 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit1, params) == 9
        # Diameter
        Da, Db = minmax(diameter(orbit1, params, 0.05), diameter(orbit1, params, 0.25))
        @test 11.9 < Da < Db < 26.7
        # Mass
        Ma, Mb = minmax(mass(orbit1, params, 2_600, 0.05), mass(orbit1, params, 2_600, 0.25))
        @test 2.3E6 < Ma < Mb < 2.6E7
    end

    @testset "Interesting NEOs" begin

        # 2014 AA hit the Earth around January 2, 2014, 02:49 UTC

        # Fetch optical astrometry
        optical = fetch_optical_mpc80("2014 AA", MPC)
        # Subset of optical for IOD
        # suboptical = iodsuboptical(optical, 3)
        # In this case: optical == suboptical

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

        # Values by Sep 28, 2025

        # Curvature
        C, Γ_C = curvature(optical, od.weights)
        σ_C = sqrt.(diag(Γ_C))
        @test all( abs.(C) ./ σ_C .> 5.7)
        χ2 = C' * inv(Γ_C) * C
        @test χ2 > 5.93e5
        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
            typeof(optical)})
        # Tracklets
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
        @test dtutc2days(date(optical[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(optical[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 1e9 )
        # Vector of residuals
        @test notout(orbit.ores) == 7
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 3e-4 )
        @test all( snr(orbit) .> 20.5)
        @test chi2(orbit) < 0.23
        @test nrms(orbit) < 0.13
        # Jacobian
        @test size(orbit.jacobian) == (6, 6)
        @test isdiag(orbit.jacobian)
        @test maximum(orbit.jacobian) < 3e-4
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        q0, σ0 = orbit(), sigmas(orbit)
        JPL_CAR = [-0.1793421909678032, 0.8874121750891107, 0.3841434101167349,
            -0.017557851117612377, -0.005781634223099801, -0.0020075106081869185]
        @test all(@. abs(q0 - JPL_CAR) / σ0 < 0.30)
        # Osculating orbital elements
        osc = osculating(orbit, params)
        @test iselliptic(osc)
        @test osc.mu == μ_S
        @test epoch(osc) == epoch(orbit) + MJD2000
        @test osc.frame == :ecliptic
        q0 = equatorial2ecliptic(orbit() - params.eph_su(epoch(orbit)))
        @test norm(q0 - osc(), Inf) < 6.2e-14
        q0, σ0 = elements(osc), sigmas(osc)
        JPL_OSC = [1.163575955666616E+00, 2.128185264087166E-01, 1.423597471953649E+00,
            5.237781301766019E+01, 1.016022028285875E+02, 3.243429036265208E+02 - 360]
        @test all(@. abs(q0 - JPL_OSC) / σ0 < 0.38)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 30.9 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(orbit, params, 0.05), diameter(orbit, params, 0.25))
        @test 1.7 < Da < Db < 3.9
        # Mass
        Ma, Mb = minmax(mass(orbit, params, 2_600, 0.05), mass(orbit, params, 2_600, 0.25))
        @test 7.1E3 < Ma < Mb < 8.0E4

        # 2008 TC3 entered the Earth's atmosphere around October 7, 2008, 02:46 UTC

        # Fetch optical astrometry
        optical = fetch_optical_mpc80("2008 TC3", MPC)
        # Subset of optical for IOD
        suboptical = iodsuboptical(optical, 3)

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           gaussorder = 2, jtlsorder = 2, jtlsiter = 200, lsiter = 1,
           significance = 0.99, outrej = false, parse_eqs = false
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, suboptical)

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by Sep 28, 2025

        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
            typeof(suboptical)})
        # Tracklets
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
        @test dtutc2days(date(suboptical[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(suboptical[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 1e4 )
        # Vector of residuals
        @test notout(orbit.ores) == 18
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 2e-5 )
        @test all( snr(orbit) .> 644)
        @test chi2(orbit) < 4.35
        @test nrms(orbit) < 0.35
        # Jacobian
        @test size(orbit.jacobian) == (6, 6)
        @test isdiag(orbit.jacobian)
        @test maximum(orbit.jacobian) < 2e-5
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        q0, σ0 = orbit(), sigmas(orbit)
        JPL_CAR = [0.9739760787551061, 0.21541704400792083, 0.09401075290627411,
            -0.00789675674941779, 0.0160619782715116, 0.006135361409943397]
        @test all(@. abs(q0 - JPL_CAR) / σ0 < 0.20)
        # Osculating orbital elements
        osc = osculating(orbit, params)
        @test iselliptic(osc)
        @test osc.mu == μ_S
        @test epoch(osc) == epoch(orbit) + MJD2000
        @test osc.frame == :ecliptic
        q0 = equatorial2ecliptic(orbit() - params.eph_su(epoch(orbit)))
        @test norm(q0 - osc(), Inf) < 5.0e-14
        q0, σ0 = elements(osc), sigmas(osc)
        JPL_OSC = [1.273091758414584E+00, 2.870222798582721E-01, 2.341999526552296E+00,
            2.339645303327229E+02, 1.941265709953888E+02, 3.288450951861228E+02 - 360]
        @test all(@. abs(q0 - JPL_OSC) / σ0 < 0.09)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 30.4 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 8
        # Diameter
        Da, Db = minmax(diameter(orbit, params, 0.05), diameter(orbit, params, 0.25))
        @test 2.1 < Da < Db < 4.9
        # Mass
        Ma, Mb = minmax(mass(orbit, params, 2_600, 0.05), mass(orbit, params, 2_600, 0.25))
        @test 1.4E4 < Ma < Mb < 1.6E5

        # Add observations
        suboptical = iodsuboptical(optical, 10)

        # Refine orbit
        NEOs.update!(od, suboptical)
        orbit1 = orbitdetermination(od, orbit, params)

        # Check type
        @test isa(orbit1, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
            typeof(suboptical)})
        # Tracklets
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
        @test dtutc2days(date(suboptical[1])) > orbit1.bwd.t0 + orbit1.bwd.t[end]
        @test all( norm.(orbit1.bwd.x, Inf) .< 1 )
        @test dtutc2days(date(suboptical[end])) < orbit1.fwd.t0 + orbit1.fwd.t[end]
        @test all( norm.(orbit1.fwd.x, Inf) .< 1e15 )
        # Vector of residuals
        @test notout(orbit1.ores) == 97
        @test nout(orbit1.ores) == 0
        # Least squares fit
        @test orbit1.fit.success
        @test all( sigmas(orbit1) .< 4e-7 )
        @test all( snr(orbit1) .> 21_880)
        @test chi2(orbit1) < 54.85
        @test nrms(orbit1) < 0.53
        # Jacobian
        @test size(orbit1.jacobian) == (6, 6)
        @test isdiag(orbit1.jacobian)
        @test maximum(orbit1.jacobian) < 2e-5
        # Convergence history
        @test size(orbit1.qs, 1) == 6
        @test size(orbit1.qs, 2) == length(orbit1.Qs) == 2
        @test issorted(orbit1.Qs, rev = true)
        @test orbit1.Qs[end] == nrms(orbit1)
        # Compatibility with JPL
        q0, σ0 = orbit1(), sigmas(orbit1)
        @test all(@. abs(q0 - JPL_CAR) / σ0 < 0.17)
        # Osculating orbital elements
        osc = osculating(orbit1, params)
        @test iselliptic(osc)
        @test osc.mu == μ_S
        @test epoch(osc) == epoch(orbit1) + MJD2000
        @test osc.frame == :ecliptic
        q0 = equatorial2ecliptic(orbit1() - params.eph_su(epoch(orbit1)))
        @test norm(q0 - osc(), Inf) < 5.8e-14
        q0, σ0 = elements(osc), sigmas(osc)
        @test all(@. abs(q0 - JPL_OSC) / σ0 < 0.13)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 30.4 ≤ H + dH
        # Parameters uncertainty
        @test all(sigmas(orbit1) .< sigmas(orbit))
        # TODO: understand better differences wrt JPL solutions
        # @test nrms(orbit1) < nrms(orbit)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit1, params) == 5
        # Diameter
        Da, Db = minmax(diameter(orbit1, params, 0.05), diameter(orbit1, params, 0.25))
        @test 2.4 < Da < Db < 5.4
        # Mass
        Ma, Mb = minmax(mass(orbit1, params, 2_600, 0.05), mass(orbit1, params, 2_600, 0.25))
        @test 1.9E4 < Ma < Mb < 2.2E5
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
        optical = fetch_optical_mpc80("2023 QR6", MPC)
        flag, optical = opticalfilter(optical)

        # Parameters
        params = Parameters(
            coeffstol = Inf, bwdoffset = 0.042, fwdoffset = 0.042, # Propagation
            safegauss = true, refscale = :log,                     # Gauss method
            adamiter = 500, adamQtol = 1e-5,                       # ADAM
            jtlsiter = 20, lsiter = 10, significance = 0.99,       # Least squares
            outrej = true, χ2_rec = 7.0, χ2_rej = 8.0,             # Outlier rejection
            fudge = 100.0, max_per = 20.0
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, optical)

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params; initcond = iodinitcond)

        # Values by Sep 28, 2025

        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
            typeof(optical)})
        # Tracklets
        @test flag
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
        @test dtutc2days(date(optical[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(optical[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.ores) == 6
        @test nout(orbit.ores) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 0.018 )
        @test all( snr(orbit) .> 7.14)
        @test chi2(orbit) < 0.91
        @test nrms(orbit) < 0.28
        # Jacobian
        @test size(orbit.jacobian) == (6, 6)
        @test isdiag(orbit.jacobian)
        @test maximum(orbit.jacobian) < 0.007
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) <= 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        q0, σ0 = orbit(), sigmas(orbit)
        JPL_CAR = [0.827266656726981, -0.8060653913101916, -0.6506187674672722,
            0.01660013577219304, -0.005614737443087259, 0.002899489877794496]
        @test all(@. abs(q0 - JPL_CAR) / σ0 < 0.59)
        # Osculating orbital elements
        osc = osculating(orbit, params)
        @test iselliptic(osc)
        @test osc.mu == μ_S
        @test epoch(osc) == epoch(orbit) + MJD2000
        @test osc.frame == :ecliptic
        q0 = equatorial2ecliptic(orbit() - params.eph_su(epoch(orbit)))
        @test norm(q0 - osc(), Inf) < 6.0e-14
        q0, σ0 = elements(osc), sigmas(osc)
        JPL_OSC = [2.279881958167905E+00, 7.595854924208774E-01, 3.859999487009846E+01,
            2.293324466168822E+02, 3.254353160068554E+02, 2.033836565098925E+01]
        @test all(@. abs(q0 - JPL_OSC) / σ0 < 0.59)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 18.5 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
        # Diameter
        Da, Db = minmax(diameter(orbit, params, 0.05), diameter(orbit, params, 0.25))
        @test 558.9 < Da < Db < 1250.0
        # Mass
        Ma, Mb = minmax(mass(orbit, params, 2_600, 0.05), mass(orbit, params, 2_600, 0.25))
        @test 2.3E11 < Ma < Mb < 2.7E12
    end

    @testset "Radar astrometry" begin
        using NEOs: RadarResidual

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

        # Preliminary orbit (only optical astrometry)
        loadjpleph()
        jd0 = datetime2julian(DateTime(2005, 1, 29))
        q00 = kmsec2auday(apophisposvel199(julian2etsecs(jd0)))
        orbit0 = LeastSquaresOrbit(od0, q00, jd0, params)

        # Orbit determination problem (both optical and radar astrometry)
        od1 = ODProblem(newtonian!, optical, radar)

        # Refine orbit (both optical and radar astrometry)
        orbit1 = orbitdetermination(od1, orbit0, params)

        # Values by Sep 28, 2025

        # Check type
        @test isa(orbit1, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64,
            typeof(optical), typeof(radar), Vector{RadarResidual{Float64, Float64}}})
        # Astrometry
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
        @test dtutc2days(date(optical[1])) > orbit1.bwd.t0 + orbit1.bwd.t[end]
        @test dtutc2days(date(radar[1])) > orbit1.bwd.t0 + orbit1.bwd.t[end]
        @test all( norm.(orbit1.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(optical[end])) < orbit1.fwd.t0 + orbit1.fwd.t[end]
        @test dtutc2days(date(radar[end])) < orbit1.fwd.t0 + orbit1.fwd.t[end]
        @test all( norm.(orbit1.fwd.x, Inf) .< 2 )
        # Vectors of residuals
        @test notout(orbit1.ores) == 24
        @test notout(orbit1.rres) == 5
        @test nout(orbit1.ores) == 0
        @test nout(orbit1.rres) == 0
        # Least squares fit
        @test orbit1.fit.success
        @test all( sigmas(orbit1) .< 2.9e-7 )
        @test all( snr(orbit1) .> 8_342)
        @test chi2(orbit1) < 19.7
        @test nrms(orbit1) < 0.61
        # Jacobian
        @test size(orbit1.jacobian) == (6, 6)
        @test isdiag(orbit1.jacobian)
        @test maximum(orbit1.jacobian) < 4.1e-3
        # Convergence history
        @test size(orbit1.qs, 1) == 6
        @test size(orbit1.qs, 2) == length(orbit1.Qs) <= 2
        @test issorted(orbit1.Qs, rev = true)
        @test orbit1.Qs[end] == nrms(orbit1)
        # Compatibility with JPL
        q0, σ0 = orbit1(), sigmas(orbit1)
        JPL_CAR = [-5.229992130937651E-01, 8.689454573480734E-01, 3.096174868699621E-01,
            -1.413639580483663E-02, -5.510379552549767E-03, -2.413003153288419E-03]
        @test all(@. abs(q0 - JPL_CAR) / σ0 < 4.5)
        # Osculating orbital elements
        osc = osculating(orbit1, params)
        @test iselliptic(osc)
        @test osc.mu == μ_S
        @test epoch(osc) == epoch(orbit1) + MJD2000
        @test osc.frame == :ecliptic
        q0 = equatorial2ecliptic(orbit1() - params.eph_su(epoch(orbit1)))
        @test norm(q0 - osc(), Inf) < 6.0e-14
        q0, σ0 = elements(osc), sigmas(osc)
        JPL_OSC = [9.223295977030230E-01, 1.911190963976789E-01, 3.330797253820763E+00,
            1.263744754026979E+02, 2.044798558304837E+02, 1.361032871672047E+02]
        @test all(@. abs(q0 - JPL_OSC) / σ0 < 0.65)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit1, params)
        @test H - dH ≤ 18.93 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit1, params) == 5
        # Diameter
        Da, Db = minmax(diameter(orbit1, params, 0.05), diameter(orbit1, params, 0.25))
        @test 433.6 < Da < Db < 969.7
        # Mass
        Ma, Mb = minmax(mass(orbit1, params, 2_600, 0.05), mass(orbit1, params, 2_600, 0.25))
        @test 1.1E11 < Ma < Mb < 1.3E12
    end

end