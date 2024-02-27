# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using PlanetaryEphemeris
using Dates
using TaylorIntegration
using JLD2
using Test

using InteractiveUtils: methodswith

@testset "Orbit propagation" begin

    @testset "Integration methods" begin

        @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_ng_eph_threads!}, TaylorIntegration.jetcoeffs!))
        @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_ng_eph_threads!}, TaylorIntegration._allocate_jetcoeffs!))

        @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_eph_threads!}, TaylorIntegration.jetcoeffs!))
        @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_eph_threads!}, TaylorIntegration._allocate_jetcoeffs!))

    end

    using PlanetaryEphemeris: selecteph, ea, su, daysec, auday2kmsec
    using Statistics

    @testset "Orbit propagation without nongravs: 2023 DW" begin

        # Dynamical function
        dynamics = RNp1BP_pN_A_J23E_J2S_eph_threads!
        # Initial time [Julian date]
        jd0 = datetime2julian(DateTime(2023, 2, 25, 0, 0, 0))
        # Time of integration [years]
        nyears = 0.2
        # Unperturbed initial condition
        q0 = [-9.759018085743707E-01, 3.896554445697074E-01, 1.478066121706831E-01,
              -9.071450085084557E-03, -9.353197026254517E-03, -5.610023032269034E-03]
        # Propagation parameters
        params = NEOParameters(maxsteps = 1, order = 25, abstol = 1e-20,
                            parse_eqs = true)

        # Initial time [days since J2000]
        t0 = jd0 - PE.J2000
        # Solar System ephemeris
        sseph = loadpeeph(NEOs.sseph, t0 - nyears*yr, t0 + nyears*yr)
        # Sun's ephemeris
        eph_su = selecteph(sseph, su)
        # Earth's ephemeris
        eph_ea = selecteph(sseph, ea)

        # Warmup propagation (forward)
        NEOs.propagate(
            dynamics,
            jd0,
            nyears,
            q0,
            params
        )
        # Propagate orbit
        params = NEOParameters(params; maxsteps = 1_000)
        sol_bwd = NEOs.propagate(
            dynamics,
            jd0,
            -nyears,
            q0,
            params
        )
        sol = NEOs.propagate(
            dynamics,
            jd0,
            nyears,
            q0,
            params
        )

        # Check that solution saves correctly
        jldsave("test.jld2"; sol_bwd, sol)
        recovered_sol = JLD2.load("test.jld2", "sol")
        recovered_sol_bwd = JLD2.load("test.jld2", "sol_bwd")
        @test sol == recovered_sol
        @test sol_bwd == recovered_sol_bwd
        rm("test.jld2")

        @test sol_bwd.t0 == sol.t0
        @test (sol_bwd.t[end]-sol_bwd.t[1])/yr ≈ -nyears
        @test (sol.t[end]-sol.t[1])/yr ≈ nyears
        @test sol(sol.t0) == q0
        q_fwd_end = [-1.0168239304400228, -0.3800432452351079, -0.2685901784950398,
                     0.007623614213394988, -0.00961901551025335, -0.004682171726467166]
        @test norm(sol(sol.t0 + sol.t[end])-q_fwd_end, Inf) < 1e-12
        @test sol_bwd(sol_bwd.t0) == q0
        q_bwd_end = [0.2689956497466164, 0.4198851302334139, 0.2438053951982368,
                     -0.018875911266050937, 0.0167349306087375, 0.007789382070881366]
        @test norm(sol_bwd(sol_bwd.t0 + sol_bwd.t[end])-q_bwd_end, Inf) < 1e-12

        # Read optical astrometry file

        obs_radec_mpc_2023DW = read_radec_mpc(joinpath("data", "RADEC_2023_DW.dat"))

        # Compute residuals
        _res_ = NEOs.residuals(
            obs_radec_mpc_2023DW,
            params,
            xve = t -> auday2kmsec(eph_ea(t/daysec)),
            xvs = t -> auday2kmsec(eph_su(t/daysec)),
            xva = t -> auday2kmsec(sol(t/daysec))
        )
        res, w = NEOs.unfold(_res_)

        mean_radec0 = mean(res)
        std_radec0 = std(res)
        rms_radec0 = nrms(res, ones(length(res))) # un-normalized RMS
        @test mean_radec0 ≈ -0.667 atol=1e-3
        @test std_radec0 ≈ 0.736 atol=1e-3
        @test rms_radec0 ≈ 0.992 atol=1e-3

        # Propagate orbit with perturbed initial conditions
        q1 = q0 + vcat(1e-3randn(3), 1e-5randn(3))
        sol1 = NEOs.propagate(
            dynamics,
            jd0,
            nyears,
            q1,
            params
        )

        # Check that solution saves correctly
        jldsave("test.jld2"; sol1 = sol1)
        recovered_sol1 = JLD2.load("test.jld2", "sol1")
        @test sol1 == recovered_sol1
        rm("test.jld2")

        # Compute residuals for orbit with perturbed initial conditions
        _res1_ = NEOs.residuals(
            obs_radec_mpc_2023DW,
            params,
            xve = t -> auday2kmsec(eph_ea(t/daysec)),
            xvs = t -> auday2kmsec(eph_su(t/daysec)),
            xva = t -> auday2kmsec(sol1(t/daysec))
        )
        res1, _ = NEOs.unfold(_res1_)
        mean_radec1 = mean(res1)
        std_radec1 = std(res1)
        rms_radec1 = nrms(res1, ones(length(res1)))

        @test abs(mean_radec1) ≥ abs(mean_radec0)
        @test std_radec1 ≥ std_radec0
        @test rms_radec1 ≥ rms_radec0
    end

    @testset "Orbit propagation with nongravs: (99942) Apophis" begin

        # Dynamical function
        dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
        # Initial time [Julian date]
        jd0 = datetime2julian(DateTime(2004, 6, 1))
        # Time of integration [years]
        nyears = 9.0
        # JPL #199 solution for Apophis at June 1st, 2004
        q0 = [-1.0506628055913627, -0.06064314196134998, -0.04997102228887035,
              0.0029591421121582077, -0.01423233538611057, -0.005218412537773594,
              -5.592839897872e-14, 0.0]
        # Propagation parameters
        params = NEOParameters(maxsteps = 1, order = 25, abstol = 1e-20,
                            parse_eqs = true)

        # Initial time [days since J2000]
        t0 = jd0 - PE.J2000
        # Solar System ephemeris
        sseph = loadpeeph(NEOs.sseph, t0, t0 + nyears*yr)
        # Sun's ephemeris
        eph_su = selecteph(sseph, su)
        # Earth's ephemeris
        eph_ea = selecteph(sseph, ea)

        # Warmup propagation
        sol = NEOs.propagate(
            dynamics,
            jd0,
            nyears,
            q0,
            params
        )
        # Propagate orbit
        params = NEOParameters(params, maxsteps = 5_000)
        sol = NEOs.propagate(
            dynamics,
            jd0,
            nyears,
            q0,
            params
        )

        # Check that solution saves correctly
        jldsave("test.jld2"; sol = sol)
        recovered_sol = JLD2.load("test.jld2", "sol")
        @test sol == recovered_sol
        rm("test.jld2")

        # Read optical astrometry file
        obs_radec_mpc_apophis = read_radec_mpc(joinpath("data", "99942_Tholen_etal_2013.dat"))

        # Compute optical astrometry residuals
        _res_radec_ = NEOs.residuals(
            obs_radec_mpc_apophis,
            params,
            xve = t -> auday2kmsec(eph_ea(t/daysec)),
            xvs = t -> auday2kmsec(eph_su(t/daysec)),
            xva = t -> auday2kmsec(sol(t/daysec))
        )
        res_radec, w_radec = NEOs.unfold(_res_radec_)
        nobsopt = round(Int, length(res_radec))

        # Compute mean optical astrometric residual (right ascension and declination)
        res_ra = res_radec[1:round(Int,nobsopt/2)]
        res_dec = res_radec[round(Int,nobsopt/2)+1:end]
        mean_ra = mean(res_ra)
        mean_dec = mean(res_dec)
        std_ra = std(res_ra)
        std_dec = std(res_dec)
        rms_ra = nrms(res_ra,ones(length(res_ra)))
        rms_dec = nrms(res_dec,ones(length(res_dec)))
        @test mean_ra ≈ 0.0224 atol=1e-4
        @test std_ra ≈ 0.136 atol=1e-3
        @test rms_ra ≈ std_ra atol=1e-2
        @test mean_dec ≈ -0.0124 atol=1e-2
        @test std_dec ≈ 0.0714 atol=1e-2
        @test rms_dec ≈ std_dec atol=1e-2

        # Read radar astrometry file
        deldop_2005_2013 = NEOs.read_radar_jpl(joinpath("data", "99942_RADAR_2005_2013.dat"))

        # Compute mean radar (time-delay and Doppler-shift) residuals
        @time res_del, w_del, res_dop, w_dop = residuals(
            deldop_2005_2013,
            xve = t -> auday2kmsec(eph_ea(t/daysec)),
            xvs = t -> auday2kmsec(eph_su(t/daysec)),
            xva = t -> auday2kmsec(sol(t/daysec)),
            niter = 4,
            tord = 5
        )

        mean_del = mean(res_del)
        mean_dop = mean(res_dop)
        std_del = std(res_del)
        std_dop = std(res_dop)

        @test mean_del ≈ 0.281 atol=1e-2
        @test mean_dop ≈ -0.084 atol=1e-2
        @test std_del ≈ 1.246 atol=1e-2
        @test std_dop ≈ 0.286 atol=1e-2

        res = vcat(res_radec, res_del, res_dop)
        w = vcat(w_radec, w_del, w_dop)

        # Total normalized RMS
        @test nrms(res, w) ≈  0.366 atol=1e-3
    end

    @testset "Jet transport propagation and TaylorN serialization" begin

        # Test integration (Apophis)

        # Dynamical function
        dynamics = RNp1BP_pN_A_J23E_J2S_eph_threads!
        # Initial date of integration [julian days]
        jd0 = datetime2julian(DateTime(2029, 4, 13, 20))
        # Time of integration [years]
        nyears = 0.02
        # Perturbation to nominal initial condition (Taylor1 jet transport)
        dq = NEOs.scaled_variables(order=1)
        # Initial conditions
        q0 = [-0.9170913888342959, -0.37154308794738056, -0.1610606989484252,
              0.009701519087787077, -0.012766026792868212, -0.0043488589639194275] + dq
        # Propagation parameters
        params = NEOParameters(maxsteps = 10, order = 25, abstol = 1e-20, parse_eqs = true)

        # test parsed vs non-parsed propagation
        sol = NEOs.propagate(dynamics, jd0, 1.0, q0, params)
        params = NEOParameters(params, parse_eqs = false)
        solnp = NEOs.propagate(dynamics, jd0, 1.0, q0, params)
        @test sol.t == solnp.t
        # TODO: fix roundoff differences near deep close approach in 2029
        @test norm(sol.x-solnp.x, Inf)/norm(solnp.x, Inf) < 4e-18 # 3.757708512785821e-20
        jldsave("test.jld2"; sol)
        recovered_sol = JLD2.load("test.jld2", "sol")
        @test sol == recovered_sol
        rm("test.jld2")

        params = NEOParameters(params; maxsteps = 1)
        sol, tvS, xvS, gvS = NEOs.propagate_root(dynamics, jd0, nyears, q0, params)

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
        jldsave("test.jld2"; random_TaylorN = random_TaylorN)
        recovered_taylorN = JLD2.load("test.jld2", "random_TaylorN")
        @test recovered_taylorN == random_TaylorN
        rm("test.jld2")

    end

    @testset "Jet transport orbit propagation and astrometric observables: (99942) Apophis" begin

        # Dynamical functions
        dynamicsg  = RNp1BP_pN_A_J23E_J2S_eph_threads!
        dynamicsng = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
        # integration parameters
        nyears::Float64 = 10.0
        varorder::Int = 1
        jd0::Float64 = datetime2julian(DateTime(2004,6,1)) #Julian date of integration initial time
        # 7-DOF nominal solution from pha/apophis.jl script at epoch 2004-06-01T00:00:00.000 (TDB)
        q00::Vector{Float64} = [-1.0506627988664696, -0.060643124245514164, -0.0499709975200415, 0.0029591416313078838, -0.014232335581939919, -0.0052184125285361415, -2.898870403031058e-14, 0.0]
        scalings::Vector{Float64} = vcat(fill(1e-8, 6), 1e-14)
        dq::Vector{TaylorN{Float64}} = NEOs.scaled_variables("δx", scalings, order = varorder)
        q0::Vector{TaylorN{Float64}} = q00 + vcat(dq, zero(dq[1]))

        # test parsed vs non-parsed propagation: gravity-only model
        params = NEOParameters(maxsteps=10, order=25, abstol=1e-20, parse_eqs=true)
        sol   = NEOs.propagate(dynamicsg, jd0, nyears, q0[1:6], params)
        params = NEOParameters(params, parse_eqs=false)
        solnp = NEOs.propagate(dynamicsg, jd0, nyears, q0[1:6], params)
        @test sol.t == solnp.t
        @test norm(sol.x-solnp.x, Inf) < 1e-16
        @test sol == solnp

        # test parsed vs non-parsed propagation: nongravitational model
        params = NEOParameters(params, parse_eqs=true)
        sol   = NEOs.propagate(dynamicsng, jd0, nyears, q0, params)
        params = NEOParameters(params, parse_eqs=false)
        solnp = NEOs.propagate(dynamicsng, jd0, nyears, q0, params)
        @test sol.t == solnp.t
        @test norm(sol.x-solnp.x, Inf) < 1e-16
        @test sol == solnp

        # propagate orbit (nongrav model)
        params = NEOParameters(params, maxsteps = 2_000, parse_eqs = true)
        sol = NEOs.propagate(
            dynamicsng,
            jd0,
            nyears,
            q0,
            params
        )

        # Solar System ephemeris
        sseph_obs = loadpeeph(NEOs.sseph, sol.t0, sol.t0 + sol.t[end])
        # Sun's ephemeris
        eph_su = selecteph(sseph_obs, su)
        # Earth's ephemeris
        eph_ea = selecteph(sseph_obs, ea)

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
        obs_radec_mpc_apophis = read_radec_mpc(joinpath("data", "99942_Tholen_etal_2013.dat"))

        # Compute optical astrometry residuals
        _res_radec_ = NEOs.residuals(obs_radec_mpc_apophis, params; xvs, xve, xva)
        res_radec, w_radec = NEOs.unfold(_res_radec_)
        nobsopt = round(Int, length(res_radec))

        # Compute mean optical astrometric residual (right ascension and declination)
        res_ra = res_radec[1:round(Int,nobsopt/2)]()
        res_dec = res_radec[round(Int,nobsopt/2)+1:end]()
        mean_ra = mean(res_ra)
        mean_dec = mean(res_dec)
        std_ra = std(res_ra)
        std_dec = std(res_dec)
        rms_ra = nrms(res_ra,ones(length(res_ra)))
        rms_dec = nrms(res_dec,ones(length(res_dec)))
        @test mean_ra ≈  0.0083 atol=1e-4
        @test std_ra ≈ 0.136 atol=1e-3
        @test rms_ra ≈ std_ra atol=1e-2
        @test mean_dec ≈ -0.0124 atol=1e-2
        @test std_dec ≈ 0.0714 atol=1e-2
        @test rms_dec ≈ std_dec atol=1e-2

        # Read radar astrometry file
        deldop_2005_2013 = NEOs.read_radar_jpl(joinpath("data", "99942_RADAR_2005_2013.dat"))

        # Compute mean radar (time-delay and Doppler-shift) residuals
        @time res_del, w_del, res_dop, w_dop = residuals(deldop_2005_2013[1:4]; xvs, xve,
            xva, niter=10, tord=10)

        # Doppler astrometry normalized residuals (i.e., residual/sigma) at nominal solution
        @test abs(res_dop[1]()) ≤ deldop_2005_2013[1].Δν_σ
        @test abs(res_dop[2]()) ≤ deldop_2005_2013[2].Δν_σ
        @test abs(res_dop[3]()) ≤ deldop_2005_2013[3].Δν_σ
        @test abs(res_dop[4]()) ≤ deldop_2005_2013[4].Δν_σ
        # delay astrometry normalized residuals (i.e., residual/sigma) at nominal solution
        @test abs(res_del[1]()) ≤ deldop_2005_2013[2].Δτ_σ
        @test abs(res_del[2]()) ≤ deldop_2005_2013[2].Δτ_σ

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
