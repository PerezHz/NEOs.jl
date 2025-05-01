# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using PlanetaryEphemeris
using LinearAlgebra
using Test

using NEOs: RadecMPC, numberofdays, nobs, minmaxdates, notout, nout, reduce_tracklets,
    indices
using Statistics: mean

function iodsubradec(radec::Vector{RadecMPC{T}}, N::Int = 3) where {T <: Real}
    tracklets = reduce_tracklets(radec)
    idxs = indices(tracklets[1:N])
    subradec = radec[idxs]
    return subradec
end

@testset "Orbit Determination" begin
    @testset "Straight Gauss Method" begin
        # Load observations
        radec = read_radec_mpc(joinpath(pkgdir(NEOs), "test", "data", "RADEC_2023_DW.dat"))
        # Subset of radec for IOD
        subradec = iodsubradec(radec, 3)

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           gaussorder = 2, jtlsorder = 2, jtlsiter = 200, lsiter = 1,
           significance = 0.99, outrej = false, parse_eqs = false
        )
        params = Parameters(params, parse_eqs = true)
        # Orbit determination problem
        od = ODProblem(newtonian!, subradec)

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by Apr 25, 2025

        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64})
        # Tracklets
        @test length(subradec) == nobs(od) == nobs(orbit) == 9
        @test numberofdays(subradec) == numberofdays(orbit) < 0.18
        @test minmaxdates(orbit) == (date(subradec[1]), date(subradec[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 3
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].radec[1] == subradec[1]
        @test orbit.tracklets[end].radec[end] == subradec[end]
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), dtutc2days(date(od.tracklets[2])), atol = 4e-4)
        @test dtutc2days(date(subradec[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(subradec[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.res) == 9
        @test nout(orbit.res) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 9e-4 )
        @test all( snr(orbit) .> 14.5)
        @test chi2(orbit) < 0.53
        @test nrms(orbit) < 0.18
        # Jacobian
        @test size(orbit.J) == (6, 6)
        @test isdiag(orbit.J)
        @test maximum(orbit.J) < 9e-4
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) == 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL = [-0.9867704701732631, 0.3781890325424674, 0.14094513213009532,
            -0.008773157203087259, -0.00947109649687576, -0.005654229864757284]
        @test all(abs.(orbit() - JPL) ./ sigmas(orbit) .< 0.76)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 24.3 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 10

        # Add observations
        subradec = iodsubradec(radec, 15)

        # Refine orbit
        NEOs.update!(od, subradec)
        orbit1 = orbitdetermination(od, orbit, params)

        # Check type
        @test isa(orbit1, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64})
        # Tracklets
        @test length(subradec) == nobs(od) == nobs(orbit1) == 43
        @test numberofdays(subradec) == numberofdays(orbit1) < 2.76
        @test minmaxdates(orbit1) == (date(subradec[1]), date(subradec[end]))
        @test length(od.tracklets) == length(orbit1.tracklets) == 15
        @test od.tracklets == orbit1.tracklets
        @test orbit1.tracklets[1].radec[1] == subradec[1]
        @test orbit1.tracklets[end].radec[end] == subradec[end]
        @test issorted(orbit1.tracklets)
        # Backward (forward) integration
        @test epoch(orbit1) == epoch(orbit)
        @test dtutc2days(date(subradec[1])) > orbit1.bwd.t0 + orbit1.bwd.t[end]
        @test all( norm.(orbit1.bwd.x, Inf) .< 1.2 )
        @test dtutc2days(date(subradec[end])) < orbit1.fwd.t0 + orbit1.fwd.t[end]
        @test all( norm.(orbit1.fwd.x, Inf) .< 1.2 )
        # Vector of residuals
        @test notout(orbit1.res) == 43
        @test nout(orbit1.res) == 0
        # Least squares fit
        @test orbit1.fit.success
        @test all( sigmas(orbit1) .< 2e-4 )
        @test all( snr(orbit1) .> 866 )
        @test chi2(orbit1) < 11.64
        @test nrms(orbit1) < 0.37
        # Jacobian
        @test size(orbit1.J) == (6, 6)
        @test isdiag(orbit1.J)
        @test maximum(orbit1.J) < 9e-4
        # Convergence history
        @test size(orbit1.qs, 1) == 6
        @test size(orbit1.qs, 2) == length(orbit1.Qs) == 6
        @test issorted(orbit1.Qs, rev = true)
        @test orbit1.Qs[end] == nrms(orbit1)
        # Compatibility with JPL
        @test all(abs.(orbit1() - JPL) ./ sigmas(orbit1) .< 0.31)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit1, params)
        @test H - dH ≤ 24.3 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit1, params) == 7
    end

    @testset "Unsafe Gauss Method" begin
        # Load observations
        radec = fetch_radec_mpc("2005 TM173")
        # Subset of radec for IOD
        # subradec = iodsubradec(radec, 3)
        # In this case: radec == subradec

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           gaussorder = 2, jtlsorder = 2, jtlsiter = 200, lsiter = 1,
           significance = 0.99, outrej = false, safegauss = false
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, radec)

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by Apr 25, 2025

        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64})
        # Tracklets
        @test length(radec) == nobs(od) == nobs(orbit) == 6
        @test numberofdays(radec) == numberofdays(orbit) < 1.95
        @test minmaxdates(orbit) == (date(radec[1]), date(radec[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 2
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].radec[1] == radec[1]
        @test orbit.tracklets[end].radec[end] == radec[end]
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), dtutc2days(date(radec[4])), atol = 3e-4)
        @test dtutc2days(date(radec[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(radec[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.res) == 6
        @test nout(orbit.res) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 5e-3 )
        @test all( snr(orbit) .> 21.4)
        @test chi2(orbit) < 2.53
        @test nrms(orbit) < 0.46
        # Jacobian
        @test size(orbit.J) == (6, 6)
        @test isdiag(orbit.J)
        @test maximum(orbit.J) < 6.8e-4
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) == 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL = [1.0042569058151192, 0.2231639040146286, 0.11513854178693468,
            -0.010824212819531798, 0.017428798232689943, 0.0071046780555307385]
        @test all(abs.(orbit() - JPL) ./ sigmas(orbit) .< 6e-3)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 24.0 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 10
    end

    @testset "Gauss Method with ADAM refinement" begin
        # Load observations
        radec = fetch_radec_mpc("2024 MK")
        # Subset of radec for IOD
        subradec = radec[10:21]

        # Parameters
        params = Parameters(bwdoffset = 0.007, fwdoffset = 0.007)
        # Orbit determination problem
        od = ODProblem(newtonian!, subradec)

        # Initial Orbit Determination
        orbit = gaussiod(od, params)

        # Values by Apr 25, 2025

        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64})
        # Tracklets
        @test length(subradec) == nobs(od) == nobs(orbit) == 12
        @test numberofdays(subradec) == numberofdays(orbit) < 42.8
        @test minmaxdates(orbit) == (date(subradec[1]), date(subradec[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 3
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].radec[1] == subradec[1]
        @test orbit.tracklets[end].radec[end] == subradec[end]
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), dtutc2days(date(od.tracklets[2])), atol = 4e-4)
        @test dtutc2days(date(subradec[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(subradec[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.res) == 12
        @test nout(orbit.res) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 6.6e-4 )
        @test all( snr(orbit) .> 38.8)
        @test chi2(orbit) < 2.43
        @test nrms(orbit) < 0.32
        # Jacobian
        @test size(orbit.J) == (6, 6)
        @test isdiag(orbit.J)
        @test maximum(orbit.J) < 9.5e-6
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) == 3
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL = [-0.12722461679828806, -0.9466098076903212, -0.4526816007640767,
            0.02048875631534963, -0.00022720097573790754, 0.00321302850930331]
        @test all(abs.(orbit() - JPL) ./ sigmas(orbit) .< 0.16)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 21.7 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 9
    end

    @testset "Admissible region" begin
        using NEOs: AdmissibleRegion, arenergydis, rangerate, rangerates,
            argoldensearch, arboundary, _helmaxrange, R_SI, k_gauss, μ_ES,
            topo2bary, bary2topo

        # Fetch optical astrometry
        radec = fetch_radec_mpc("2024 BX1")
        # Parameters
        params = Parameters()
        # First tracklet
        radec = radec[1:3]
        tracklet = reduce_tracklets(radec)[1]
        # Admissible region
        A = AdmissibleRegion(tracklet, params)

        # Values by Apr 23, 2024

        # Zero AdmissibleRegion
        @test iszero(zero(AdmissibleRegion{Float64}))
        # Custom print
        @test string(A) == "AE: [116.61547, 45.39840, -3.21667, 5.76667] \
            t: 2024-01-20T21:50:15.360 obs: GINOP-KHK, Piszkesteto"
        # Coefficients
        @test length(A.coeffs) == 6
        @test A.coeffs[3] == A.v_α^2 * cos(A.δ)^2 + A.v_δ^2  # proper motion squared
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
        C, Γ_C = curvature(radec)
        σ_C = sqrt.(diag(Γ_C))
        @test all( abs.(C) ./ σ_C .> 0.02)
        χ2 = C' * inv(Γ_C) * C
        @test χ2 > 2.7
    end

    @testset "Too Short Arc" begin
        # Fetch optical astrometry
        radec = fetch_radec_mpc("2008 EK68")
        # Subset of radec for IOD
        # subradec = iodsubradec(radec, 3)
        # In this case: radec == subradec

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           tsaorder = 2, adamiter = 500, adamQtol = 1e-5,
           jtlsorder = 2, jtlsiter = 200, lsiter = 1,
           significance = 0.99, outrej = false
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, radec)

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by Apr 25, 2025

        # Curvature
        C, Γ_C = curvature(radec)
        σ_C = sqrt.(diag(Γ_C))
        @test all( abs.(C) ./ σ_C .> 5.5)
        χ2 = C' * inv(Γ_C) * C
        @test χ2 > 2_516
        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64})
        # Tracklets
        @test length(radec) == nobs(od) == nobs(orbit) == 10
        @test numberofdays(radec) == numberofdays(orbit) < 0.05
        @test minmaxdates(orbit) == (date(radec[1]), date(radec[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 1
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].radec[1] == radec[1]
        @test orbit.tracklets[end].radec[end] == radec[end]
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), mean(r -> dtutc2days(date(r)), radec), atol = 7e-5)
        @test dtutc2days(date(radec[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(radec[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.res) == 10
        @test nout(orbit.res) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 6e-3 )
        @test all( snr(orbit) .> 4.1)
        @test chi2(orbit) < 14.24
        @test nrms(orbit) < 0.85
        # Jacobian
        @test size(orbit.J) == (6, 6)
        @test isdiag(orbit.J)
        @test maximum(orbit.J) < 9.2e-3
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) == 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL = [-0.9698405495747651, 0.24035304578776012, 0.10288276585828428,
            -0.009512301266159554, -0.01532548565855646, -0.00809464581680694]
        @test all(abs.(orbit() - JPL) ./ sigmas(orbit) .< 0.013)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 29.6 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 10
    end

    @testset "Outlier Rejection" begin
        # Fetch optical astrometry
        radec = fetch_radec_mpc("2007 VV7")
        # Subset of radec for IOD
        subradec = iodsubradec(radec, 3)

        # Parameters
        params = Parameters(bwdoffset = 0.007, fwdoffset = 0.007,
            outrej = true, χ2_rec = 1.0, χ2_rej = 1.25, fudge = 0.0)
        # Orbit determination problem
        od = ODProblem(newtonian!, subradec)

        # Initial Orbit Determination (with outlier rejection)
        orbit = initialorbitdetermination(od, params)

        # Values by Apr 25, 2025

        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64})
        # Tracklets
        @test length(subradec) == nobs(od) == nobs(orbit) == 18
        @test numberofdays(subradec) == numberofdays(orbit) < 2.16
        @test minmaxdates(orbit) == (date(subradec[1]), date(subradec[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 3
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].radec[1] == subradec[1]
        @test orbit.tracklets[end].radec[end] == subradec[end]
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), dtutc2days(date(od.tracklets[2])), atol = 2e-3)
        @test dtutc2days(date(subradec[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(subradec[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.res) == 16
        @test nout(orbit.res) == 2
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 4e-3 )
        @test all( snr(orbit) .> 50)
        @test chi2(orbit) < 1.55
        @test nrms(orbit) < 0.22
        # Jacobian
        @test size(orbit.J) == (6, 6)
        @test isdiag(orbit.J)
        @test maximum(orbit.J) < 5e-4
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) == 5
        # @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL = [0.7673366466815864, 0.6484892781853565, 0.29323267343908294,
            -0.011023343781911974, 0.015392697071667377, 0.006528842022004942]
        @test all(abs.(orbit() - JPL) ./ sigmas(orbit) .< 0.07)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 26.7 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 10

        # Add remaining observations
        NEOs.update!(od, radec)
        # Refine orbit (with outlier rejection)
        orbit1 = orbitdetermination(od, orbit, params)

        # Check type
        @test isa(orbit1, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64})
        # Tracklets
        @test length(radec) == nobs(od) == nobs(orbit1) == 21
        @test numberofdays(radec) == numberofdays(orbit1) < 3.03
        @test minmaxdates(orbit1) == (date(radec[1]), date(radec[end]))
        @test length(od.tracklets) == length(orbit1.tracklets) == 4
        @test od.tracklets == orbit1.tracklets
        @test orbit1.tracklets[1].radec[1] == radec[1]
        @test orbit1.tracklets[end].radec[end] == radec[end]
        @test issorted(orbit1.tracklets)
        # Backward (forward) integration
        @test epoch(orbit1) == epoch(orbit)
        @test dtutc2days(date(radec[1])) > orbit1.bwd.t0 + orbit1.bwd.t[end]
        @test all( norm.(orbit1.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(radec[end])) < orbit1.fwd.t0 + orbit1.fwd.t[end]
        @test all( norm.(orbit1.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit1.res) == 19
        @test nout(orbit1.res) == 2
        # Least squares fit
        @test orbit1.fit.success
        @test all( sigmas(orbit1) .< 3e-4 )
        @test all( snr(orbit1) .> 574)
        @test chi2(orbit1) < 2.38
        @test nrms(orbit1) < 0.25
        # Jacobian
        @test size(orbit1.J) == (6, 6)
        @test isdiag(orbit1.J)
        @test maximum(orbit1.J) < 4e-3
        # Convergence history
        @test size(orbit1.qs, 1) == 6
        # (26/04/2025) There are roundoff differences in the nrms of the two
        # jtls iterations; hence, in some os/julia versions, the first (second)
        # iteration has the lowest nrms.
        # @test size(orbit1.qs, 2) == length(orbit1.Qs) == 1
        @test issorted(orbit1.Qs, rev = true)
        @test orbit1.Qs[end] == nrms(orbit1)
        # Compatibility with JPL
        @test all(abs.(orbit1() - JPL) ./ sigmas(orbit1) .< 7e-4)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 26.7 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit1, params) == 9
    end

    @testset "Interesting NEOs" begin

        # 2014 AA hit the Earth around January 2, 2014, 02:49 UTC

        # Fetch optical astrometry
        radec = fetch_radec_mpc("2014 AA")
        # Subset of radec for IOD
        # subradec = iodsubradec(radec, 3)
        # In this case: radec == subradec

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           tsaorder = 2, adamiter = 500, adamQtol = 1e-5,
           jtlsorder = 2, jtlsiter = 200, lsiter = 1,
           significance = 0.99, outrej = false
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, radec)

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by Apr 25, 2025

        # Curvature
        C, Γ_C = curvature(radec)
        σ_C = sqrt.(diag(Γ_C))
        @test all( abs.(C) ./ σ_C .> 7.5)
        χ2 = C' * inv(Γ_C) * C
        @test χ2 > 1.03e6
        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64})
        # Tracklets
        @test length(radec) == nobs(od) == nobs(orbit) == 7
        @test numberofdays(radec) == numberofdays(orbit) < 0.05
        @test minmaxdates(orbit) == (date(radec[1]), date(radec[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 1
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].radec[1] == radec[1]
        @test orbit.tracklets[end].radec[end] == radec[end]
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), mean(r -> dtutc2days(date(r)), radec), atol = 2e-5)
        @test dtutc2days(date(radec[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(radec[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 1e9 )
        # Vector of residuals
        @test notout(orbit.res) == 7
        @test nout(orbit.res) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 3e-4 )
        @test all( snr(orbit) .> 20.5)
        @test chi2(orbit) < 0.23
        @test nrms(orbit) < 0.13
        # Jacobian
        @test size(orbit.J) == (6, 6)
        @test isdiag(orbit.J)
        @test maximum(orbit.J) < 3e-4
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) == 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL = [-0.1793421909678032, 0.8874121750891107, 0.3841434101167349,
            -0.017557851117612377, -0.005781634223099801, -0.0020075106081869185]
        @test all(abs.(orbit() - JPL) ./ sigmas(orbit) .< 0.3)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 30.9 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 10

        # 2008 TC3 entered the Earth's atmosphere around October 7, 2008, 02:46 UTC

        # Fetch optical astrometry
        radec = fetch_radec_mpc("2008 TC3")
        # Subset of radec for IOD
        subradec = iodsubradec(radec, 3)

        # Parameters
        params = Parameters(
           coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007,
           gaussorder = 2, jtlsorder = 2, jtlsiter = 200, lsiter = 1,
           significance = 0.99, outrej = false, parse_eqs = false
        )
        # Orbit determination problem
        od = ODProblem(newtonian!, subradec)

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params)

        # Values by Apr 25, 2025

        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64})
        # Tracklets
        @test length(subradec) == nobs(od) == nobs(orbit) == 18
        @test numberofdays(subradec) == numberofdays(orbit) < 0.34
        @test minmaxdates(orbit) == (date(subradec[1]), date(subradec[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 3
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].radec[1] == subradec[1]
        @test orbit.tracklets[end].radec[end] == subradec[end]
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), dtutc2days(date(od.tracklets[2])), atol = 4e-3)
        @test dtutc2days(date(subradec[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(subradec[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 1e4 )
        # Vector of residuals
        @test notout(orbit.res) == 18
        @test nout(orbit.res) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 2e-5 )
        @test all( snr(orbit) .> 644)
        @test chi2(orbit) < 4.35
        @test nrms(orbit) < 0.35
        # Jacobian
        @test size(orbit.J) == (6, 6)
        @test isdiag(orbit.J)
        @test maximum(orbit.J) < 2e-5
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) == 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL = [0.9739760787551061, 0.21541704400792083, 0.09401075290627411,
            -0.00789675674941779, 0.0160619782715116, 0.006135361409943397]
        @test all(abs.(orbit() - JPL) ./ sigmas(orbit) .< 0.20)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 30.4 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 8

        # Add observations
        subradec = iodsubradec(radec, 10)

        # Refine orbit
        NEOs.update!(od, subradec)
        orbit1 = orbitdetermination(od, orbit, params)

        # Check type
        @test isa(orbit1, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64})
        # Tracklets
        @test length(subradec) == nobs(od) == nobs(orbit1) == 97
        @test numberofdays(subradec) == numberofdays(orbit1) < 0.70
        @test minmaxdates(orbit1) == (date(subradec[1]), date(subradec[end]))
        @test length(od.tracklets) == length(orbit1.tracklets) == 10
        @test od.tracklets == orbit1.tracklets
        @test orbit1.tracklets[1].radec[1] == subradec[1]
        @test orbit1.tracklets[end].radec[end] == subradec[93]
        @test issorted(orbit1.tracklets)
        # Backward (forward) integration
        @test epoch(orbit1) == epoch(orbit)
        @test dtutc2days(date(subradec[1])) > orbit1.bwd.t0 + orbit1.bwd.t[end]
        @test all( norm.(orbit1.bwd.x, Inf) .< 1 )
        @test dtutc2days(date(subradec[end])) < orbit1.fwd.t0 + orbit1.fwd.t[end]
        @test all( norm.(orbit1.fwd.x, Inf) .< 1e15 )
        # Vector of residuals
        @test notout(orbit1.res) == 97
        @test nout(orbit1.res) == 0
        # Least squares fit
        @test orbit1.fit.success
        @test all( sigmas(orbit1) .< 4e-7 )
        @test all( snr(orbit1) .> 21_880)
        @test chi2(orbit1) < 54.85
        @test nrms(orbit1) < 0.53
        # Jacobian
        @test size(orbit1.J) == (6, 6)
        @test isdiag(orbit1.J)
        @test maximum(orbit1.J) < 2e-5
        # Convergence history
        @test size(orbit1.qs, 1) == 6
        @test size(orbit1.qs, 2) == length(orbit1.Qs) == 2
        @test issorted(orbit1.Qs, rev = true)
        @test orbit1.Qs[end] == nrms(orbit1)
        # Compatibility with JPL
        @test all(abs.(orbit1() - JPL) ./ sigmas(orbit1) .< 0.17)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 30.4 ≤ H + dH
        # Parameters uncertainty
        @test all(sigmas(orbit1) .< sigmas(orbit))
        # TODO: understand better differences wrt JPL solutions
        # @test nrms(orbit1) < nrms(orbit)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit1, params) == 5
    end

    @testset "research/Celestial Mech. Dyn. Astron. 2025/orbitdetermination.jl" begin
        using NEOs: iodinitcond, issatellite

        function radecfilter(radec::Vector{RadecMPC{T}}) where {T <: Real}
            # Eliminate observations before oficial discovery
            firstobs = findfirst(r -> !isempty(r.discovery), radec)
            isnothing(firstobs) && return false, radec
            radec = radec[firstobs:end]
            # Filter out incompatible observations
            filter!(radec) do r
                hascoord(r.observatory) && !issatellite(r.observatory) &&
                date(r) > DateTime(2000, 1, 1, 12)
            end
            length(radec) < 3 && return false, radec
            # Find the first set of 3 tracklets with a < 15 days timespan
            tracklets = reduce_tracklets(radec)
            for i in 1:length(tracklets)-2
                numberofdays(tracklets[i:i+2]) > 15.0 && continue
                tracklets = tracklets[i:i+2]
                radec = reduce(vcat, getfield.(tracklets, :radec))
                sort!(radec)
                break
            end
            return numberofdays(radec) <= 15.0, radec
        end

        # Fetch and filter optical astrometry
        radec = fetch_radec_mpc("2023 QR6")
        flag, radec = radecfilter(radec)

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
        od = ODProblem(newtonian!, radec)

        # Initial Orbit Determination
        orbit = initialorbitdetermination(od, params; initcond = iodinitcond)

        # Values by Apr 25, 2025

        # Check type
        @test isa(orbit, LeastSquaresOrbit{typeof(newtonian!), Float64, Float64})
        # Tracklets
        @test flag
        @test length(radec) == nobs(od) == nobs(orbit) == 6
        @test numberofdays(radec) == numberofdays(orbit) < 6.22
        @test minmaxdates(orbit) == (date(radec[1]), date(radec[end]))
        @test length(od.tracklets) == length(orbit.tracklets) == 3
        @test od.tracklets == orbit.tracklets
        @test orbit.tracklets[1].radec[1] == radec[1]
        @test orbit.tracklets[end].radec[end] == radec[end]
        @test issorted(orbit.tracklets)
        # Backward (forward) integration
        @test isapprox(epoch(orbit), dtutc2days(date(od.tracklets[2])), atol = 3e-3)
        @test dtutc2days(date(radec[1])) > orbit.bwd.t0 + orbit.bwd.t[end]
        @test all( norm.(orbit.bwd.x, Inf) .< 2 )
        @test dtutc2days(date(radec[end])) < orbit.fwd.t0 + orbit.fwd.t[end]
        @test all( norm.(orbit.fwd.x, Inf) .< 2 )
        # Vector of residuals
        @test notout(orbit.res) == 6
        @test nout(orbit.res) == 0
        # Least squares fit
        @test orbit.fit.success
        @test all( sigmas(orbit) .< 0.018 )
        @test all( snr(orbit) .> 7.14)
        @test chi2(orbit) < 0.91
        @test nrms(orbit) < 0.28
        # Jacobian
        @test size(orbit.J) == (6, 6)
        @test isdiag(orbit.J)
        @test maximum(orbit.J) < 0.007
        # Convergence history
        @test size(orbit.qs, 1) == 6
        @test size(orbit.qs, 2) == length(orbit.Qs) == 2
        @test issorted(orbit.Qs, rev = true)
        @test orbit.Qs[end] == nrms(orbit)
        # Compatibility with JPL
        JPL = [ 0.827266656726981, -0.8060653913101916, -0.6506187674672722,
            0.01660013577219304, -0.005614737443087259, 0.002899489877794496]
        @test all(abs.(orbit() - JPL) ./ sigmas(orbit) .< 0.59)
        # Absolute magnitude
        H, dH = absolutemagnitude(orbit, params)
        @test H - dH ≤ 18.5 ≤ H + dH
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(orbit, params) == 10
    end

end