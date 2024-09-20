# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using PlanetaryEphemeris
using LinearAlgebra
using Test

using NEOs: NEOSolution, numberofdays

@testset "Orbit Determination" begin
    @testset "Gauss Method (without ADAM)" begin
        # Load observations
        radec = read_radec_mpc(joinpath(pkgdir(NEOs), "test", "data", "RADEC_2023_DW.dat"))
        # Parameters
        params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007, parse_eqs = false)
        params = NEOParameters(params, parse_eqs=true)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by Sep 20, 2024

        # Vector of observations
        @test length(radec) == 123
        @test numberofdays(radec) < 21.0
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 44
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test datetime2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 123
        @test iszero(count(outlier.(sol.res)))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 4.2e-5 )
        @test all( snr(sol) .> 4_700)
        @test nrms(sol) < 0.38
        # Jacobian
        @test size(sol.jacobian) == (6, 6)
        @test !isdiag(sol.jacobian)
        @test maximum(sol.jacobian) < 8e-4
        # Compatibility with JPL
        JPL = [-1.1003236797145037, 0.20774505704014837, 0.04203643429323372,
               -0.004736048200346307, -0.010626587751050683, -0.006016238714906758]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.29)

        # Refine orbit
        sol1 = orbitdetermination(radec, sol, params)

        # Orbit solution
        @test isa(sol1, NEOSolution{Float64, Float64})
        # Tracklets
        @test sol1.tracklets == sol.tracklets
        # Backward integration
        @test datetime2days(date(radec[1])) > sol1.bwd.t0 + sol1.bwd.t[end]
        @test all( norm.(sol1.bwd.x, Inf) .< 1.2 )
        @test isempty(sol1.t_bwd) && isempty(sol1.x_bwd) && isempty(sol1.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[end])) < sol1.fwd.t0 + sol1.fwd.t[end]
        @test all( norm.(sol1.fwd.x, Inf) .< 1.2 )
        @test isempty(sol1.t_fwd) && isempty(sol1.x_fwd) && isempty(sol1.g_fwd)
        # Vector of residuals
        @test length(sol1.res) == 123
        @test iszero(count(outlier.(sol1.res)))
        # Least squares fit
        @test sol1.fit.success
        @test all( sigmas(sol1) .< 4.2e-5 )
        @test all( snr(sol1) .> 4_700 )
        @test nrms(sol1) < 0.38
        # Jacobian
        @test size(sol1.jacobian) == (6, 6)
        @test isdiag(sol1.jacobian)
        @test maximum(sol1.jacobian) < 2e-6
        # Compatibility with JPL
        @test all(abs.(sol1() - JPL) ./ sigmas(sol1) .< 0.29)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(radec, sol1, params) == 6
    end

    @testset "Gauss Method (with ADAM)" begin
        # Load observations
        radec = fetch_radec_mpc("2016 TU93")

        # Parameters
        params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007, adamhelp = true)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by June 1, 2024

        # Vector of observations
        @test length(radec) == 9
        @test numberofdays(radec) < 13.1
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 3
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test datetime2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 9
        @test iszero(count(outlier.(sol.res)))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 8e-5 )
        @test all( snr(sol) .> 36)
        @test nrms(sol) < 0.46
        # Jacobian
        @test size(sol.jacobian) == (6, 6)
        @test isdiag(sol.jacobian)
        @test maximum(sol.jacobian) < 1e-5
        # Compatibility with JPL
        JPL = [1.0102564188486982, 0.2934743828145318, 0.10467187893161536,
               -0.0002634434601757652, 0.01837381321202214, 0.007208485181422459]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 4.8e-2)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(radec, sol, params) == 8
    end

    @testset "Admissible region" begin
        using NEOs: AdmissibleRegion, reduce_tracklets, arenergydis, rangerate,
                    rangerates, argoldensearch, arboundary, _helmaxrange, R_SI,
                    k_gauss, μ_ES, topo2bary, bary2topo

        # Fetch optical astrometry
        radec = fetch_radec_mpc("2024 BX1")
        # Parameters
        params = NEOParameters()
        # First tracklet
        radec = radec[1:3]
        tracklet = reduce_tracklets(radec)[1]
        # Admissible region
        A = AdmissibleRegion(tracklet, params)

        # Values by August 19, 2024

        # Zero AdmissibleRegion
        @test iszero(zero(AdmissibleRegion{Float64}))
        # Custom print
        @test string(A) == "AE: [116.61547, 45.39840, -3.21667, 5.76667] t: 2024-01-20T21:50:15.360 obs: GINOP-KHK, Piszkesteto"
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
        @test isempty(rangerates(A, ρ0, :inner))
        @test isempty(rangerates(A, A.ρ_domain[2] + 1.0, :outer))
        @test isempty(rangerates(A, ρ0 + 1.0, :inner))
        @test rangerate(A, A.ρ_domain[1], :min, :outer) == A.v_ρ_domain[1]
        @test rangerate(A, A.ρ_domain[1], :max, :outer) == A.v_ρ_domain[2]
        @test rangerate(A, A.ρ_domain[1], :min, :inner) ≈ -rangerate(A, A.ρ_domain[1], :max, :inner) atol = 1e-18
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

        # Parameters
        params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by May 31, 2024

        # Vector of observations
        @test length(radec) == 10
        @test numberofdays(radec) < 0.05
        # Curvature
        C, Γ_C = curvature(radec)
        σ_C = sqrt.(diag(Γ_C))
        @test all( abs.(C) ./ σ_C .> 5.5)
        χ2 = C' * inv(Γ_C) * C
        @test χ2 > 2_517
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 1
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test datetime2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 10
        @test iszero(count(outlier.(sol.res)))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 6e-3 )
        @test all( snr(sol) .> 4.1)
        @test nrms(sol) < 0.85
        # Jacobian
        @test size(sol.jacobian) == (6, 6)
        @test isdiag(sol.jacobian)
        @test maximum(sol.jacobian) < 1e-5
        # Compatibility with JPL
        JPL = [-0.9698333704468041, 0.2403646120906576, 0.10288887497365079,
               -0.009512521364762891, -0.015325432155116774, -0.008094623535119534]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.015)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(radec, sol, params) == 11
    end

    @testset "Outlier Rejection" begin
        # Fetch optical astrometry
        radec = fetch_radec_mpc("2007 VV7")

        # Parameters
        params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007)

        # Orbit Determination
        sol = orbitdetermination(radec, params)
        # Outlier rejection
        sol = outlier_rejection(radec, sol, params)

        # Values by Sep 20, 2024

        # Vector of observations
        @test length(radec) == 21
        @test numberofdays(radec) < 3.03
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 4
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test datetime2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 21
        @test count(outlier.(sol.res)) == 2
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 3e-4 )
        @test all( snr(sol) .> 574)
        @test nrms(sol) < 0.25
        # Jacobian
        @test size(sol.jacobian) == (6, 6)
        @test isdiag(sol.jacobian)
        @test maximum(sol.jacobian) < 8e-7
        # Compatibility with JPL
        JPL = [0.7673449629397204, 0.6484776654615118, 0.2932277478785896,
               -0.011023192686652665, 0.015392823966811551, 0.0065288994881745974]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 2e-3)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(radec, sol, params) == 8
    end

    @testset "Interesting NEOs" begin

        # 2014 AA hit the Earth around January 2, 2014, 02:49 UTC

        # Fetch optical astrometry
        radec = fetch_radec_mpc("2014 AA")

        # Parameters
        params = NEOParameters(coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007)

        # Orbit Determination
        sol = orbitdetermination(radec, params)

        # Values by May 31 2024

        # Vector of observations
        @test length(radec) == 7
        @test numberofdays(radec) < 0.05
        # Curvature
        C, Γ_C = curvature(radec)
        σ_C = sqrt.(diag(Γ_C))
        @test all( abs.(C) ./ σ_C .> 7.5)
        χ2 = C' * inv(Γ_C) * C
        @test χ2 > 1.03e6
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 1
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test datetime2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 1e9 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 7
        @test iszero(count(outlier.(sol.res)))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 3e-4 )
        @test all( snr(sol) .> 20.5)
        @test nrms(sol) < 0.13
        # Jacobian
        @test size(sol.jacobian) == (6, 6)
        @test isdiag(sol.jacobian)
        @test maximum(sol.jacobian) < 9e-6
        # Compatibility with JPL
        JPL = [-0.17932853781716676, 0.8874166708195785, 0.38414497112938667,
               -0.017557883503263098, -0.005781328976995571, -0.0020073946372627465]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.3)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(radec, sol, params) == 10

        # 2008 TC3 entered the Earth's atmosphere around October 7, 2008, 02:46 UTC

        # Fetch optical astrometry
        radec = fetch_radec_mpc("2008 TC3")

        # Parameters
        params = NEOParameters(coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007)

        # Observations with <1" weight
        idxs = findall(x -> w8sveres17(x) < 1, radec)
        # Restricted Orbit Determination
        sol = orbitdetermination(radec[idxs], params)

        # Values by May 31, 2024

        # Vector of observations
        @test length(radec) == 883
        @test numberofdays(radec) < 0.80
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 2
        @test sol.tracklets[1].radec[1] == radec[idxs[1]]
        @test sol.tracklets[end].radec[end] == radec[idxs[end]]
        @test issorted(sol.tracklets)
        # Backward integration
        @test datetime2days(date(radec[idxs[1]])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[idxs[end]])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 1e4 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 20
        @test iszero(count(outlier.(sol.res)))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 2e-5 )
        @test all( snr(sol) .> 732)
        @test nrms(sol) < 0.30
        # Jacobian
        @test size(sol.jacobian) == (6, 6)
        @test isdiag(sol.jacobian)
        @test maximum(sol.jacobian) < 1e-5
        # Compatibility with JPL
        JPL = [0.9741070120439872, 0.2151506132683409, 0.0939089782825125,
               -0.007890445003489774, 0.016062726197895585, 0.006136042044307594]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.2)

        # Update orbit with more observations
        sol1 = orbitdetermination(radec[1:30], sol, params)

        # Orbit solution
        @test isa(sol1, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol1.tracklets) == 5
        @test sol1.tracklets[1].radec[1] == radec[1]
        @test sol1.tracklets[end].radec[end] == radec[30]
        @test issorted(sol1.tracklets)
        # Backward integration
        @test datetime2days(date(radec[idxs[1]])) > sol1.bwd.t0 + sol1.bwd.t[end]
        @test all( norm.(sol1.bwd.x, Inf) .< 1 )
        @test isempty(sol1.t_bwd) && isempty(sol1.x_bwd) && isempty(sol1.g_bwd)
        # Forward integration
        @test datetime2days(date(radec[idxs[end]])) < sol1.fwd.t0 + sol1.fwd.t[end]
        @test all( norm.(sol1.fwd.x, Inf) .< 1e4 )
        @test isempty(sol1.t_fwd) && isempty(sol1.x_fwd) && isempty(sol1.g_fwd)
        # Vector of residuals
        @test length(sol1.res) == 30
        @test iszero(count(outlier.(sol1.res)))
        # Least squares fit
        @test sol1.fit.success
        @test all( sigmas(sol1) .< 2e-6 )
        @test all( snr(sol1) .> 4_281)
        @test nrms(sol1) < 0.37
        # Jacobian
        @test size(sol1.jacobian) == (6, 6)
        @test isdiag(sol1.jacobian)
        @test maximum(sol1.jacobian) < 1e-6
        # Compatibility with JPL
        @test all(abs.(sol1() - JPL) ./ sigmas(sol1) .< 1.21)
        # Parameters uncertainty
        @test all(sigmas(sol1) .< sigmas(sol))
        # TODO: understand better differences wrt JPL solutions
        # @test nrms(sol1) < nrms(sol)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(radec[1:30], sol1, params) == 6
    end

end