# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using PlanetaryEphemeris
using LinearAlgebra
using Test

using NEOs: NEOSolution, RadecMPC, reduce_tracklets,
    indices, numberofdays, nout, scaled_variables

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

        @test length(subradec) == 9
        @test numberofdays(subradec) < 0.18

        # Parameters
        params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007, parse_eqs = false)
        params = NEOParameters(params, parse_eqs = true)
        # Orbit determination problem
        od = ODProblem(newtonian!, subradec)

        # Initial Orbit Determination
        sol = orbitdetermination(od, params)

        # Values by Oct 1, 2024

        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 3
        @test sol.tracklets[1].radec[1] == subradec[1]
        @test sol.tracklets[end].radec[end] == subradec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test dtutc2days(date(subradec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test dtutc2days(date(subradec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 9
        @test iszero(nout(sol.res))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 9e-4 )
        @test all( snr(sol) .> 14.5)
        @test nrms(sol) < 0.18
        # Jacobian
        @test size(sol.jacobian) == (6, 6)
        @test !isdiag(sol.jacobian)
        @test maximum(sol.jacobian) < 3e-4
        # Compatibility with JPL
        JPL = [-0.9867704701732631, 0.3781890325424674, 0.14094513213009532,
            -0.008773157203087259, -0.00947109649687576, -0.005654229864757284]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.76)

        # Add observations
        subradec = iodsubradec(radec, 15)

        @test length(subradec) == 43
        @test numberofdays(subradec) < 2.76

        # Refine orbit
        NEOs.update!(od, subradec)
        sol1 = orbitdetermination(od, sol, params)

        # Orbit solution
        @test isa(sol1, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol1.tracklets) == 15
        @test sol1.tracklets[1].radec[1] == subradec[1]
        @test sol1.tracklets[end].radec[end] == subradec[end]
        @test issorted(sol1.tracklets)
        # Backward integration
        @test dtutc2days(date(subradec[1])) > sol1.bwd.t0 + sol1.bwd.t[end]
        @test all( norm.(sol1.bwd.x, Inf) .< 1.2 )
        @test isempty(sol1.t_bwd) && isempty(sol1.x_bwd) && isempty(sol1.g_bwd)
        # Forward integration
        @test dtutc2days(date(subradec[end])) < sol1.fwd.t0 + sol1.fwd.t[end]
        @test all( norm.(sol1.fwd.x, Inf) .< 1.2 )
        @test isempty(sol1.t_fwd) && isempty(sol1.x_fwd) && isempty(sol1.g_fwd)
        # Vector of residuals
        @test length(sol1.res) == 43
        @test iszero(nout(sol1.res))
        # Least squares fit
        @test sol1.fit.success
        @test all( sigmas(sol1) .< 2e-4 )
        @test all( snr(sol1) .> 866 )
        @test nrms(sol1) < 0.37
        # Jacobian
        @test size(sol1.jacobian) == (6, 6)
        @test isdiag(sol1.jacobian)
        @test maximum(sol1.jacobian) < 1e-6
        # Compatibility with JPL
        @test all(abs.(sol1() - JPL) ./ sigmas(sol1) .< 0.31)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(od, sol1, params) == 7
    end

    @testset "Unsafe Gauss Method" begin
        # Load observations
        radec = fetch_radec_mpc("2005 TM173")
        # Subset of radec for IOD
        # subradec = iodsubradec(radec, 3)
        # In this case: radec == subradec

        @test length(radec) == 6
        @test numberofdays(radec) < 1.95

        # Parameters
        params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007, safegauss = false)
        # Orbit determination problem
        od = ODProblem(newtonian!, radec)

        # Initial Orbit Determination
        sol = orbitdetermination(od, params)

        # Values by Dec 11, 2024

        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 2
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test dtutc2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test dtutc2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 6
        @test iszero(nout(sol.res))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 5e-3 )
        @test all( snr(sol) .> 21.5)
        @test nrms(sol) < 0.46
        # Jacobian
        @test size(sol.jacobian) == (6, 6)
        @test !isdiag(sol.jacobian)
        @test maximum(sol.jacobian) < 6.1e-4
        # Compatibility with JPL
        JPL = [1.0042569058151192, 0.2231639040146286, 0.11513854178693468,
            -0.010824212819531798, 0.017428798232689943, 0.0071046780555307385]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 8.1e-3)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(od, sol, params) == 10
    end

    @testset "Gauss Method with ADAM refinement" begin
        # Load observations
        radec = fetch_radec_mpc("2024 MK")
        # Subset of radec for IOD
        subradec = radec[10:21]

        @test length(subradec) == 12
        @test numberofdays(subradec) < 42.8

        # Parameters
        params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007)
        # Orbit determination problem
        od = ODProblem(newtonian!, subradec)

        # Initial Orbit Determination
        varorder = max(params.tsaorder, params.gaussorder, params.jtlsorder)
        scaled_variables("dx", ones(6); order = varorder)
        sol = gaussiod(od, params)

        # Values by Dec 11, 2024

        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 3
        @test sol.tracklets[1].radec[1] == subradec[1]
        @test sol.tracklets[end].radec[end] == subradec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test dtutc2days(date(subradec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test dtutc2days(date(subradec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 12
        @test iszero(nout(sol.res))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 6.6e-4 )
        @test all( snr(sol) .> 38.8)
        @test nrms(sol) < 0.32
        # Jacobian
        @test size(sol.jacobian) == (6, 6)
        @test isdiag(sol.jacobian)
        @test maximum(sol.jacobian) < 9.5e-6
        # Compatibility with JPL
        JPL = [-0.12722461679828806, -0.9466098076903212, -0.4526816007640767,
            0.02048875631534963, -0.00022720097573790754, 0.00321302850930331]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.16)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(od, sol, params) == 9
    end

    @testset "Admissible region" begin
        using NEOs: AdmissibleRegion, arenergydis, rangerate, rangerates,
            argoldensearch, arboundary, _helmaxrange, R_SI, k_gauss, μ_ES,
            topo2bary, bary2topo

        # Fetch optical astrometry
        radec = fetch_radec_mpc("2024 BX1")
        # Parameters
        params = NEOParameters()
        # First tracklet
        radec = radec[1:3]
        tracklet = reduce_tracklets(radec)[1]
        # Admissible region
        A = AdmissibleRegion(tracklet, params)

        # Values by Oct 1, 2024

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
        @test isempty(rangerates(A, ρ0, :inner))
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

        @test length(radec) == 10
        @test numberofdays(radec) < 0.05

        # Parameters
        params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007)
        # Orbit determination problem
        od = ODProblem(newtonian!, radec)

        # Initial Orbit Determination
        sol = orbitdetermination(od, params)

        # Values by Oct 1, 2024

        # Curvature
        C, Γ_C = curvature(radec)
        σ_C = sqrt.(diag(Γ_C))
        @test all( abs.(C) ./ σ_C .> 5.5)
        χ2 = C' * inv(Γ_C) * C
        @test χ2 > 2_516
        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 1
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test dtutc2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test dtutc2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 10
        @test iszero(nout(sol.res))
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
        JPL = [-0.9698405495747651, 0.24035304578776012, 0.10288276585828428,
            -0.009512301266159554, -0.01532548565855646, -0.00809464581680694]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.012)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(od, sol, params) == 11
    end

    @testset "Outlier Rejection" begin
        # Fetch optical astrometry
        radec = fetch_radec_mpc("2007 VV7")
        # Subset of radec for IOD
        subradec = iodsubradec(radec, 3)

        @test length(subradec) == 18
        @test numberofdays(subradec) < 2.16

        # Parameters
        params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007,
            outrej = true, χ2_rec = 1.0, χ2_rej = 1.25, fudge = 0.0)
        # Orbit determination problem
        od = ODProblem(newtonian!, subradec)

        # Initial Orbit Determination (with outlier rejection)
        sol = orbitdetermination(od, params)

        # Values by Oct 11, 2024

        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 3
        @test sol.tracklets[1].radec[1] == subradec[1]
        @test sol.tracklets[end].radec[end] == subradec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test dtutc2days(date(subradec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test dtutc2days(date(subradec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 18
        @test nout(sol.res) == 2
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 4e-3 )
        @test all( snr(sol) .> 50)
        @test nrms(sol) < 0.22
        # Jacobian
        @test size(sol.jacobian) == (6, 6)
        @test !isdiag(sol.jacobian)
        @test maximum(sol.jacobian) < 5e-4
        # Compatibility with JPL
        JPL = [0.7673366466815864, 0.6484892781853565, 0.29323267343908294,
            -0.011023343781911974, 0.015392697071667377, 0.006528842022004942]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.07)

        # Add remaining observations
        NEOs.update!(od, radec)
        # Refine orbit (with outlier rejection)
        sol1 = orbitdetermination(od, sol, params)

        # Orbit solution
        @test isa(sol1, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol1.tracklets) == 4
        @test sol1.tracklets[1].radec[1] == radec[1]
        @test sol1.tracklets[end].radec[end] == radec[end]
        @test issorted(sol1.tracklets)
        # Backward integration
        @test dtutc2days(date(radec[1])) > sol1.bwd.t0 + sol1.bwd.t[end]
        @test all( norm.(sol1.bwd.x, Inf) .< 2 )
        @test isempty(sol1.t_bwd) && isempty(sol1.x_bwd) && isempty(sol1.g_bwd)
        # Forward integration
        @test dtutc2days(date(radec[end])) < sol1.fwd.t0 + sol1.fwd.t[end]
        @test all( norm.(sol1.fwd.x, Inf) .< 2 )
        @test isempty(sol1.t_fwd) && isempty(sol1.x_fwd) && isempty(sol1.g_fwd)
        # Vector of residuals
        @test length(sol1.res) == 21
        @test nout(sol1.res) == 2
        # Least squares fit
        @test sol1.fit.success
        @test all( sigmas(sol1) .< 3e-4 )
        @test all( snr(sol1) .> 574)
        @test nrms(sol1) < 0.25
        # Jacobian
        @test size(sol1.jacobian) == (6, 6)
        @test isdiag(sol1.jacobian)
        @test maximum(sol1.jacobian) < 8e-7
        # Compatibility with JPL
        @test all(abs.(sol1() - JPL) ./ sigmas(sol1) .< 7e-4)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(od, sol1, params) == 8
    end

    @testset "Interesting NEOs" begin

        # 2014 AA hit the Earth around January 2, 2014, 02:49 UTC

        # Fetch optical astrometry
        radec = fetch_radec_mpc("2014 AA")
        # Subset of radec for IOD
        # subradec = iodsubradec(radec, 3)
        # In this case: radec == subradec

        @test length(radec) == 7
        @test numberofdays(radec) < 0.05

        # Parameters
        params = NEOParameters(coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007)
        # Orbit determination problem
        od = ODProblem(newtonian!, radec)

        # Initial Orbit Determination
        sol = orbitdetermination(od, params)

        # Values by Oct 1, 2024

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
        @test dtutc2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test dtutc2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 1e9 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 7
        @test iszero(nout(sol.res))
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
        JPL = [-0.1793421909678032, 0.8874121750891107, 0.3841434101167349,
            -0.017557851117612377, -0.005781634223099801, -0.0020075106081869185]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.3)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(od, sol, params) == 10

        # 2008 TC3 entered the Earth's atmosphere around October 7, 2008, 02:46 UTC

        # Fetch optical astrometry
        radec = fetch_radec_mpc("2008 TC3")
        # Subset of radec for IOD
        subradec = iodsubradec(radec, 3)

        @test length(subradec) == 18
        @test numberofdays(subradec) < 0.34

        # Parameters
        params = NEOParameters(coeffstol = Inf, bwdoffset = 0.007, fwdoffset = 0.007)
        # Orbit determination problem
        od = ODProblem(newtonian!, subradec)

        # Initial Orbit Determination
        sol = orbitdetermination(od, params)

        # Values by Oct 1, 2024

        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 3
        @test sol.tracklets[1].radec[1] == subradec[1]
        @test sol.tracklets[end].radec[end] == subradec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test dtutc2days(date(subradec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test dtutc2days(date(subradec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 1e4 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 18
        @test iszero(nout(sol.res))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 2e-5 )
        @test all( snr(sol) .> 645)
        @test nrms(sol) < 0.35
        # Jacobian
        @test size(sol.jacobian) == (6, 6)
        @test !isdiag(sol.jacobian)
        @test maximum(sol.jacobian) < 4e-6
        # Compatibility with JPL
        JPL = [0.9739760787551061, 0.21541704400792083, 0.09401075290627411,
            -0.00789675674941779, 0.0160619782715116, 0.006135361409943397]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.20)

        # Add observations
        subradec = iodsubradec(radec, 10)

        @test length(subradec) == 97
        @test numberofdays(subradec) < 0.70

        # Refine orbit
        NEOs.update!(od, subradec)
        sol1 = orbitdetermination(od, sol, params)

        # Orbit solution
        @test isa(sol1, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol1.tracklets) == 10
        @test sol1.tracklets[1].radec[1] == subradec[1]
        @test sol1.tracklets[end].radec[end] == subradec[93]
        @test issorted(sol1.tracklets)
        # Backward integration
        @test dtutc2days(date(subradec[1])) > sol1.bwd.t0 + sol1.bwd.t[end]
        @test all( norm.(sol1.bwd.x, Inf) .< 1 )
        @test isempty(sol1.t_bwd) && isempty(sol1.x_bwd) && isempty(sol1.g_bwd)
        # Forward integration
        @test dtutc2days(date(subradec[end])) < sol1.fwd.t0 + sol1.fwd.t[end]
        @test all( norm.(sol1.fwd.x, Inf) .< 1e15 )
        @test isempty(sol1.t_fwd) && isempty(sol1.x_fwd) && isempty(sol1.g_fwd)
        # Vector of residuals
        @test length(sol1.res) == 97
        @test iszero(nout(sol1.res))
        # Least squares fit
        @test sol1.fit.success
        @test all( sigmas(sol1) .< 4e-7 )
        @test all( snr(sol1) .> 21_880)
        @test nrms(sol1) < 0.53
        # Jacobian
        @test size(sol1.jacobian) == (6, 6)
        @test isdiag(sol1.jacobian)
        @test maximum(sol1.jacobian) < 1e-6
        # Compatibility with JPL
        @test all(abs.(sol1() - JPL) ./ sigmas(sol1) .< 0.17)
        # Parameters uncertainty
        @test all(sigmas(sol1) .< sigmas(sol))
        # TODO: understand better differences wrt JPL solutions
        # @test nrms(sol1) < nrms(sol)
        # MPC Uncertainty Parameter
        @test uncertaintyparameter(od, sol1, params) == 5
    end

    @testset "research/2024/orbitdetermination.jl" begin
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

        @test flag
        @test length(radec) == 6
        @test numberofdays(radec) < 6.22

        # Parameters
        params = NEOParameters(
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
        sol = orbitdetermination(od, params; initcond = iodinitcond)

        # Values by Nov 21, 2024

        # Orbit solution
        @test isa(sol, NEOSolution{Float64, Float64})
        # Tracklets
        @test length(sol.tracklets) == 3
        @test sol.tracklets[1].radec[1] == radec[1]
        @test sol.tracklets[end].radec[end] == radec[end]
        @test issorted(sol.tracklets)
        # Backward integration
        @test dtutc2days(date(radec[1])) > sol.bwd.t0 + sol.bwd.t[end]
        @test all( norm.(sol.bwd.x, Inf) .< 2 )
        @test isempty(sol.t_bwd) && isempty(sol.x_bwd) && isempty(sol.g_bwd)
        # Forward integration
        @test dtutc2days(date(radec[end])) < sol.fwd.t0 + sol.fwd.t[end]
        @test all( norm.(sol.fwd.x, Inf) .< 2 )
        @test isempty(sol.t_fwd) && isempty(sol.x_fwd) && isempty(sol.g_fwd)
        # Vector of residuals
        @test length(sol.res) == 6
        @test iszero(nout(sol.res))
        # Least squares fit
        @test sol.fit.success
        @test all( sigmas(sol) .< 0.018 )
        @test all( snr(sol) .> 7.14)
        @test nrms(sol) < 0.28
        # Jacobian
        @test size(sol.jacobian) == (6, 6)
        # Julia v1.11.1 gives a different result for this test (28/Nov)
        # @test isdiag(sol.jacobian)
        @test maximum(sol.jacobian) < 0.003
        # Compatibility with JPL
        JPL = [ 0.827266656726981, -0.8060653913101916, -0.6506187674672722,
            0.01660013577219304, -0.005614737443087259, 0.002899489877794496]
        @test all(abs.(sol() - JPL) ./ sigmas(sol) .< 0.59)
    end

end