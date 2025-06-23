# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using PlanetaryEphemeris
using Dates
using TaylorIntegration
using JLD2
using Test

using InteractiveUtils: methodswith

@testset "Propagation" begin

    @testset "Integration methods" begin

        using TaylorIntegration: jetcoeffs!, _allocate_jetcoeffs!

        @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_ng_eph_threads!}, jetcoeffs!))
        @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_ng_eph_threads!}, _allocate_jetcoeffs!))

        @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_eph_threads!}, jetcoeffs!))
        @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_eph_threads!}, _allocate_jetcoeffs!))

        @test !isempty(methodswith(Val{newtonian!}, jetcoeffs!))
        @test !isempty(methodswith(Val{newtonian!}, _allocate_jetcoeffs!))

    end

    using PlanetaryEphemeris: PlanetaryEphemerisSerialization, selecteph, ea, su, daysec,
          auday2kmsec
    using Statistics

    @testset "Orbit propagation without nongravs: 2023 DW" begin

        # Dynamical function
        dynamics = RNp1BP_pN_A_J23E_J2S_eph_threads!
        # Initial time [Julian date TDB]
        jd0 = datetime2julian(DateTime(2023, 2, 25, 0, 0, 0))
        # Time of integration [years]
        nyears = 0.2
        # Unperturbed initial condition
        q0 = [-9.759018085743707E-01, 3.896554445697074E-01, 1.478066121706831E-01,
              -9.071450085084557E-03, -9.353197026254517E-03, -5.610023032269034E-03]
        # Propagation parameters
        params = Parameters(maxsteps = 1, order = 25, abstol = 1e-20, parse_eqs = true)

        # Initial time [days since J2000]
        t0 = jd0 - PE.J2000
        # Sun's ephemeris
        eph_su = selecteph(NEOs.sseph, su, t0 - nyears*yr, t0 + nyears*yr)
        # Earth's ephemeris
        eph_ea = selecteph(NEOs.sseph, ea, t0 - nyears*yr, t0 + nyears*yr)

        # Warmup propagation (forward)
        NEOs.propagate(dynamics, q0, jd0, nyears, params)
        # Propagate orbit
        params = Parameters(params; maxsteps = 1_000)
        sol_bwd = NEOs.propagate(dynamics, q0, jd0, -nyears, params)
        sol_fwd = NEOs.propagate(dynamics, q0, jd0, nyears, params)

        # Check that solution saves correctly
        @test JLD2.writeas(typeof(sol_fwd)) == PlanetaryEphemerisSerialization{Float64}
        jldsave("test.jld2"; sol_bwd, sol_fwd)
        recovered_sol_fwd = JLD2.load("test.jld2", "sol_fwd")
        recovered_sol_bwd = JLD2.load("test.jld2", "sol_bwd")
        @test sol_fwd == recovered_sol_fwd
        @test sol_bwd == recovered_sol_bwd
        rm("test.jld2")

        @test sol_bwd.t0 == sol_fwd.t0 == t0
        @test (sol_bwd.t[end] - sol_bwd.t[1]) / yr ≈ -nyears
        @test (sol_fwd.t[end] - sol_fwd.t[1]) / yr ≈ nyears
        @test sol_fwd(sol_fwd.t0) == q0
        q_fwd_end = [-1.0168239304400228, -0.3800432452351079, -0.2685901784950398,
                     0.007623614213394988, -0.00961901551025335, -0.004682171726467166]
        @test norm(sol_fwd(sol_fwd.t0 + sol_fwd.t[end]) - q_fwd_end, Inf) < 1e-12
        @test sol_bwd(sol_bwd.t0) == q0
        q_bwd_end = [0.2689956497466164, 0.4198851302334139, 0.2438053951982368,
                     -0.018875911266050937, 0.0167349306087375, 0.007789382070881366]
        @test norm(sol_bwd(sol_bwd.t0 + sol_bwd.t[end])-q_bwd_end, Inf) < 1e-12

        # Read optical astrometry file
        optical_2023DW = read_optical_mpc80(joinpath(pkgdir(NEOs), "test", "data",
            "2023DW_OPTICAL.dat"))
        # Make weigths and debiasing corrections
        w8s = Veres17(optical_2023DW).w8s
        bias = Eggl20(optical_2023DW).bias

        # Compute normalized residuals
        _res_ = NEOs.residuals(
            optical_2023DW,
            w8s, bias;
            xve = t -> auday2kmsec(eph_ea(t/daysec)),
            xvs = t -> auday2kmsec(eph_su(t/daysec)),
            xva = t -> auday2kmsec(sol_fwd(t/daysec))
        )
        res, _, _ = unfold(_res_)

        mean_optical0 = mean(res)
        std_optical0 = std(res)
        chi2_optical0 = chi2(res)
        nms_optical0 = nms(res)
        nrms_optical0 = nrms(res)

        @test mean_optical0 ≈ -0.900 atol=1e-3
        @test std_optical0 ≈ 1.070 atol=1e-3
        @test chi2_optical0 ≈ 479.83 atol=1e-2
        @test nms_optical0 ≈ 1.951 atol=1e-3
        @test nrms_optical0 ≈ 1.397 atol=1e-3

        # Propagate orbit with perturbed initial conditions
        q1 = q0 + vcat(1e-3randn(3), 1e-5randn(3))
        sol1 = NEOs.propagate(dynamics, q1, jd0, nyears, params)

        # Check that solution saves correctly
        jldsave("test.jld2"; sol1)
        recovered_sol1 = JLD2.load("test.jld2", "sol1")
        @test sol1 == recovered_sol1
        rm("test.jld2")

        # Compute residuals for orbit with perturbed initial conditions
        _res1_ = NEOs.residuals(
            optical_2023DW,
            w8s, bias,
            xve = t -> auday2kmsec(eph_ea(t/daysec)),
            xvs = t -> auday2kmsec(eph_su(t/daysec)),
            xva = t -> auday2kmsec(sol1(t/daysec))
        )
        res1, _, _ = unfold(_res1_)

        mean_optical1 = mean(res1)
        std_optical1 = std(res1)
        chi2_optical1 = chi2(res1)
        nms_optical1 = nms(res1)
        nrms_optical1 = nrms(res1)

        @test abs(mean_optical1) ≥ abs(mean_optical0)
        @test std_optical1 ≥ std_optical0
        @test chi2_optical1 ≥ chi2_optical0
        @test nms_optical1 ≥ nms_optical0
        @test nrms_optical1 ≥ nrms_optical0

    end

    @testset "Orbit propagation with nongravs: (99942) Apophis" begin

        using NEOs: isdelay, isdoppler

        # Dynamical function
        dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
        # Initial time [Julian date TDB]
        jd0 = datetime2julian(DateTime(2004, 6, 1))
        # Time of integration [years]
        nyears = 9.0
        # JPL #199 solution for Apophis at June 1st, 2004
        q0 = [-1.0506628055913627, -0.06064314196134998, -0.04997102228887035,
              0.0029591421121582077, -0.01423233538611057, -0.005218412537773594,
              -5.592839897872e-14, 0.0]
        # Propagation parameters
        params = Parameters(maxsteps = 1, order = 25, abstol = 1e-20, parse_eqs = true)

        # Initial time [days since J2000]
        t0 = jd0 - PE.J2000
        # Sun's ephemeris
        eph_su = selecteph(NEOs.sseph, su, t0, t0 + nyears*yr)
        # Earth's ephemeris
        eph_ea = selecteph(NEOs.sseph, ea, t0, t0 + nyears*yr)

        # Warmup propagation
        sol = NEOs.propagate(dynamics, q0, jd0, nyears, params)
        # Propagate orbit
        params = Parameters(params, maxsteps = 5_000)
        sol = NEOs.propagate(dynamics, q0, jd0, nyears, params)

        # Check that solution saves correctly
        @test JLD2.writeas(typeof(sol)) == PlanetaryEphemerisSerialization{Float64}
        jldsave("test.jld2"; sol)
        recovered_sol = JLD2.load("test.jld2", "sol")
        @test sol == recovered_sol
        rm("test.jld2")

        # Read optical astrometry file
        optical_Apophis = read_optical_mpc80(joinpath(pkgdir(NEOs), "test", "data",
            "99942_Tholen_etal_2013.dat"))
        # Make weights and debiasing corrections
        w8s = Veres17(optical_Apophis).w8s
        bias = Eggl20(optical_Apophis).bias

        # Compute optical astrometry residuals
        res_optical = NEOs.residuals(
            optical_Apophis,
            w8s, bias,
            xve = t -> auday2kmsec(eph_ea(t/daysec)),
            xvs = t -> auday2kmsec(eph_su(t/daysec)),
            xva = t -> auday2kmsec(sol(t/daysec))
        )
        res_ra, res_dec = @. ra(res_optical), dec(res_optical)

        # Compute mean optical astrometric residual (right ascension and declination)
        mean_ra, mean_dec = mean(res_ra), mean(res_dec)
        std_ra, std_dec = std(res_ra), std(res_dec)
        chi2_ra, chi2_dec = chi2(res_ra), chi2(res_dec)
        nms_ra, nms_dec = nms(res_ra), nms(res_dec)
        nrms_ra, nrms_dec = nrms(res_ra), nrms(res_dec)

        @test mean_ra ≈ 0.0096 atol=1e-4
        @test mean_dec ≈ -0.0048 atol=1e-4
        @test std_ra ≈ 0.0852 atol=1e-4
        @test std_dec ≈ 0.0414 atol=1e-4
        @test chi2_ra ≈ 3.1699 atol=1e-4
        @test chi2_dec ≈ 0.7496 atol=1e-4
        @test nms_ra ≈ 0.0073 atol=1e-4
        @test nms_dec ≈ 0.0017 atol=1e-4
        @test nrms_ra ≈ 0.0857 atol=1e-4
        @test nrms_dec ≈ 0.0417 atol=1e-4

        # Read radar astrometry file
        radar_Apophis = read_radar_jpl(joinpath(pkgdir(NEOs), "test", "data",
            "99942_RADAR_2005_2013.json"))

        # Compute mean radar (time-delay and Doppler-shift) residuals
        res_radar = residuals(
            radar_Apophis,
            xve = t -> auday2kmsec(eph_ea(t/daysec)),
            xvs = t -> auday2kmsec(eph_su(t/daysec)),
            xva = t -> auday2kmsec(sol(t/daysec)),
            niter = 4,
            tord = 5
        )
        res_del = residual.(res_radar[isdelay.(radar_Apophis)])
        res_dop = residual.(res_radar[isdoppler.(radar_Apophis)])

        mean_del, mean_dop = mean(res_del), mean(res_dop)
        std_del, std_dop = std(res_del), std(res_dop)
        chi2_del, chi2_dop = chi2(res_del), chi2(res_dop)
        nms_del, nms_dop = nms(res_del), nms(res_dop)
        nrms_del, nrms_dop = nrms(res_del), nrms(res_dop)

        @test mean_del ≈ 0.0768 atol=1e-4
        @test mean_dop ≈ -0.5533 atol=1e-4
        @test std_del ≈ 1.6094 atol=1e-4
        @test std_dop ≈ 1.5501 atol=1e-4
        @test chi2_del ≈ 41.5448 atol=1e-4
        @test chi2_dop ≈ 76.1602 atol=1e-4
        @test nms_del ≈ 2.4438 atol=1e-4
        @test nms_dop ≈ 2.6262 atol=1e-4
        @test nrms_del ≈ 1.5633 atol=1e-4
        @test nrms_dop ≈ 1.6206 atol=1e-4

        res = vcat(res_ra, res_dec, res_del, res_dop)

        # Total statistics
        @test mean(res) ≈ -0.0140 atol=1e-4
        @test std(res) ≈ 0.3655 atol=1e-4
        @test chi2(res) ≈ 121.6245 atol=1e-4
        @test nms(res) ≈ 0.1337 atol=1e-4
        @test nrms(res) ≈ 0.3656 atol=1e-4

    end

    @testset "Jet transport propagation and TaylorN serialization" begin

        using PlanetaryEphemeris: TaylorInterpolantNSerialization

        # Test integration (Apophis)

        # Dynamical function
        dynamics = RNp1BP_pN_A_J23E_J2S_eph_threads!
        # Initial date of integration [Julian date TDB]
        jd0 = dtutc2jdtdb(DateTime(2029, 4, 13, 20))
        # Time of integration [years]
        nyears = 0.02
        # Perturbation to nominal initial condition (Taylor1 jet transport)
        dq = scaled_variables(order = 1)
        # Initial conditions
        q0 = [-0.9170913888342959, -0.37154308794738056, -0.1610606989484252,
              0.009701519087787077, -0.012766026792868212, -0.0043488589639194275] + dq
        # Propagation parameters
        params = Parameters(maxsteps = 10, order = 25, abstol = 1e-20, parse_eqs = true)

        # Test parsed vs non-parsed propagation
        sol = NEOs.propagate(dynamics, q0, jd0, 1.0, params)
        params = Parameters(params, parse_eqs = false)
        solnp = NEOs.propagate(dynamics, q0, jd0, 1.0, params)
        @test sol.t == solnp.t
        # TODO: fix roundoff differences near deep close approach in 2029
        @test norm(sol.x - solnp.x, Inf) / norm(solnp.x, Inf) < 4e-18 # 3.757708512785821e-20
        @test JLD2.writeas(typeof(sol)) == TaylorInterpolantNSerialization{Float64}
        jldsave("test.jld2"; sol)
        recovered_sol = JLD2.load("test.jld2", "sol")
        @test sol == recovered_sol
        rm("test.jld2")

        params = Parameters(params; maxsteps = 1)
        sol, tvS, xvS, gvS = NEOs.propagate_root(dynamics, q0, jd0, nyears, params)

        @test JLD2.writeas(typeof(sol)) == TaylorInterpolantNSerialization{Float64}
        jldsave("test.jld2"; sol, tvS, xvS, gvS)
        recovered_sol = JLD2.load("test.jld2", "sol")
        recovered_tvS = JLD2.load("test.jld2", "tvS")
        recovered_xvS = JLD2.load("test.jld2", "xvS")
        recovered_gvS = JLD2.load("test.jld2", "gvS")
        @test sol == recovered_sol
        @test tvS == recovered_tvS
        @test xvS == recovered_xvS
        @test gvS == recovered_gvS
        rm("test.jld2")

        # It is unlikely that such a short integration generates a non-trivial tvS, xvS and gvS.
        # Therefore, to test TaylorNSerialization I suggest to generate random TaylorN and check
        # it saves correctly...
        random_TaylorN = [cos(sum(dq .* rand(6))), sin(sum(dq .* rand(6))), tan(sum(dq .* rand(6)))]
        jldsave("test.jld2"; random_TaylorN)
        recovered_taylorN = JLD2.load("test.jld2", "random_TaylorN")
        @test recovered_taylorN == random_TaylorN
        rm("test.jld2")

    end

    @testset "Jet transport orbit propagation and astrometric observables: (99942) Apophis" begin

        using NEOs: isdelay, isdoppler

        # Dynamical functions
        dynamicsg  = RNp1BP_pN_A_J23E_J2S_eph_threads!
        dynamicsng = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
        # Integration parameters
        nyears = 10.0
        varorder = 1
        # Julian date (TDB) of integration initial time
        jd0 = datetime2julian(DateTime(2004, 6, 1))
        # 7-DOF nominal solution from pha/apophis.jl script at epoch 2004-06-01T00:00:00.000 (TDB)
        q00 = [-1.0506627988664696, -0.060643124245514164, -0.0499709975200415,
               0.0029591416313078838, -0.014232335581939919, -0.0052184125285361415,
               -2.898870403031058e-14, 0.0]
        scalings = vcat(fill(1e-8, 6), 1e-14)
        dq = scaled_variables("δx", scalings, order = varorder)
        q0 = q00 + vcat(dq, zero(dq[1]))

        # Test parsed vs non-parsed propagation: gravity-only model
        params = Parameters(maxsteps = 10, order = 25, abstol = 1e-20, parse_eqs = true)
        sol   = NEOs.propagate(dynamicsg, q0[1:6], jd0, nyears, params)
        params = Parameters(params, parse_eqs = false)
        solnp = NEOs.propagate(dynamicsg, q0[1:6], jd0, nyears, params)
        @test sol.t == solnp.t
        @test norm(sol.x - solnp.x, Inf) < 1e-16
        @test sol == solnp

        # Test parsed vs non-parsed propagation: nongravitational model
        params = Parameters(params, parse_eqs = true)
        sol   = NEOs.propagate(dynamicsng, q0, jd0, nyears, params)
        params = Parameters(params, parse_eqs = false)
        solnp = NEOs.propagate(dynamicsng, q0, jd0, nyears, params)
        @test sol.t == solnp.t
        @test norm(sol.x - solnp.x, Inf) < 1e-16
        @test sol == solnp

        # Propagate orbit (nongrav model)
        params = Parameters(params, maxsteps = 2_000, parse_eqs = true)
        sol = NEOs.propagate(dynamicsng, q0, jd0, nyears, params)

        # Sun's ephemeris
        eph_su = selecteph(NEOs.sseph, su, sol.t0, sol.t0 + sol.t[end])
        # Earth's ephemeris
        eph_ea = selecteph(NEOs.sseph, ea, sol.t0, sol.t0 + sol.t[end])

        # Apophis
        # Change t, x, v units, resp., from days, au, au/day to sec, km, km/sec
        xva(et) = auday2kmsec(sol(et/daysec))
        # Earth
        # Change x, v units, resp., from au, au/day to km, km/sec
        xve(et) = auday2kmsec(eph_ea(et/daysec))
        # Sun
        # Change x, v units, resp., from au, au/day to km, km/sec
        xvs(et) = auday2kmsec(eph_su(et/daysec))

        # Read optical astrometry file
        optical_Apophis = read_optical_mpc80(joinpath(pkgdir(NEOs), "test", "data",
            "99942_Tholen_etal_2013.dat"))
        # Make weights and debiasing corrections
        w8s = Veres17(optical_Apophis).w8s
        bias = Eggl20(optical_Apophis).bias

        # Compute optical astrometry residuals
        res_optical = NEOs.residuals(optical_Apophis, w8s, bias; xvs, xve, xva)
        res_ra, res_dec = @. ra(res_optical), dec(res_optical)

        # Compute mean optical astrometric residual (right ascension and declination)
        mean_ra, mean_dec = mean(res_ra), mean(res_dec)
        std_ra, std_dec = std(res_ra), std(res_dec)
        chi2_ra, chi2_dec = chi2(res_ra), chi2(res_dec)
        nms_ra, nms_dec = nms(res_ra), nms(res_dec)
        nrms_ra, nrms_dec = nrms(res_ra), nrms(res_dec)

        @test mean_ra ≈ 0.0047 atol=1e-4
        @test mean_dec ≈ -0.0061 atol=1e-4
        @test std_ra ≈ 0.0858 atol=1e-4
        @test std_dec ≈ 0.0413 atol=1e-4
        @test chi2_ra ≈ 3.1824 atol=1e-4
        @test chi2_dec ≈ 0.7517 atol=1e-4
        @test nms_ra ≈ 0.0074 atol=1e-4
        @test nms_dec ≈ 0.0017 atol=1e-4
        @test nrms_ra ≈ 0.0858 atol=1e-4
        @test nrms_dec ≈ 0.0417 atol=1e-4

        # Read radar astrometry file
        radar_Apophis = NEOs.read_radar_jpl(joinpath(pkgdir(NEOs), "test", "data",
            "99942_RADAR_2005_2013.json"))
        filter!(x -> year(date(x)) == 2005, radar_Apophis)
        mask_del, mask_dop = @. isdelay(radar_Apophis), isdoppler(radar_Apophis)
        del, dop = radar_Apophis[mask_del], radar_Apophis[mask_dop]

        # Compute mean radar (time-delay and Doppler-shift) residuals
        res_radar = residuals(radar_Apophis; xvs, xve, xva, niter = 10, tord = 10)
        res_del, w_del = @. residual(res_radar[mask_del]), weight(res_radar[mask_del])
        res_dop, w_dop = @. residual(res_radar[mask_dop]), weight(res_radar[mask_dop])

        # Doppler astrometry normalized residuals (i.e., residual/sigma) at nominal solution
        @test abs(res_dop[1]() / w_dop[1]) ≤ rms(dop[1])
        @test abs(res_dop[2]() / w_dop[2]) ≤ rms(dop[2])
        @test abs(res_dop[3]() / w_dop[3]) ≤ rms(dop[3])
        @test abs(res_dop[4]() / w_dop[4]) ≤ rms(dop[4])
        # delay astrometry normalized residuals (i.e., residual/sigma) at nominal solution
        @test abs(res_del[1]() / w_del[1]) ≤ rms(del[1])
        @test abs(res_del[2]() / w_del[2]) ≤ rms(del[2])

        dq_sample = 2ones(7)
        # Doppler astrometry normalized residuals at non-nominal solution
        @test abs(res_dop[1]()) ≤ abs(res_dop[1](dq_sample))
        @test abs(res_dop[2]()) ≤ abs(res_dop[2](dq_sample))
        @test abs(res_dop[3]()) ≤ abs(res_dop[3](dq_sample))
        @test abs(res_dop[4]()) ≤ abs(res_dop[4](dq_sample))
        # delay astrometry normalized residuals at non-nominal solution
        @test abs(res_del[1]()) ≤ abs(res_del[1](dq_sample))
        @test abs(res_del[2]()) ≤ abs(res_del[2](dq_sample))

    end

end
