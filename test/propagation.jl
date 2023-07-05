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

        objname = "2023DW"
        maxsteps = 1000
        nyears = 0.2
        dynamics = RNp1BP_pN_A_J23E_J2S_eph_threads!
        # Initial time [Julian date]
        jd0 = datetime2julian(DateTime(2023,2,25,0,0,0))
        # Initial time [days since J2000]
        t0 = jd0 - PE.J2000
        # unperturbed initial condition
        q0 = [-9.759018085743707E-01, 3.896554445697074E-01, 1.478066121706831E-01, -9.071450085084557E-03, -9.353197026254517E-03, -5.610023032269034E-03]
        # Solar System ephemeris
        sseph = loadpeeph(NEOs.sseph, t0 - nyears*yr, t0 + nyears*yr)
        # Sun's ephemeris
        eph_su = selecteph(sseph, su)
        # Earth's ephemeris
        eph_ea = selecteph(sseph, ea)

        # warmup propagation (backward and forward)
        NEOs.propagate(
            dynamics,
            1,
            jd0,
            -nyears,
            nyears,
            q0,
            Val(true),
            order = 25,
            abstol = 1e-20,
            parse_eqs = true
        )

        # propagate orbit
        sol_bwd, sol = NEOs.propagate(
            dynamics,
            maxsteps,
            jd0,
            -nyears,
            nyears,
            q0,
            Val(true),
            order = 25,
            abstol = 1e-20,
            parse_eqs = true
        )

        # check that solution saves correctly
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
        q_fwd_end = [-1.0168239304400228, -0.3800432452351079, -0.2685901784950398, 0.007623614213394988, -0.00961901551025335, -0.004682171726467166]
        @test norm(sol(sol.t0 + sol.t[end])-q_fwd_end, Inf) < 1e-12
        @test sol_bwd(sol_bwd.t0) == q0
        q_bwd_end = [0.2689956497466164, 0.4198851302334139, 0.2438053951982368, -0.018875911266050937, 0.0167349306087375, 0.007789382070881366]
        @test norm(sol_bwd(sol_bwd.t0 + sol_bwd.t[end])-q_bwd_end, Inf) < 1e-12

        # Read optical astrometry file

        obs_radec_mpc_2023DW = NEOs.read_radec_mpc(joinpath("data", "RADEC_2023_DW.dat"))

        # Compute residuals
        res, _ = NEOs.residuals(
            obs_radec_mpc_2023DW,
            xve=t->auday2kmsec(eph_ea(t/daysec)),
            xvs=t->auday2kmsec(eph_su(t/daysec)),
            xva=t->auday2kmsec(sol(t/daysec))
        )

        mean_radec0 = mean(res)
        std_radec0 = std(res)
        rms_radec0 = nrms(res,ones(length(res))) # un-normalized RMS
        @test mean_radec0 ≈ -0.667 atol=1e-2
        @test std_radec0 ≈ 0.736 atol=1e-2
        @test rms_radec0 ≈ 0.992 atol=1e-2

        # propagate orbit with perturbed initial conditions
        q1 = q0 + vcat(1e-3randn(3), 1e-5randn(3))
        sol1 = NEOs.propagate(
            dynamics,
            maxsteps,
            jd0,
            nyears,
            q1,
            Val(true),
            order = 25,
            abstol = 1e-20,
            parse_eqs = true
        )

        # check that solution saves correctly
        jldsave("test.jld2"; sol1 = sol1)
        recovered_sol1 = JLD2.load("test.jld2", "sol1")
        @test sol1 == recovered_sol1
        rm("test.jld2")

        # compute residuals for orbit with perturbed initial conditions
        res1, _ = NEOs.residuals(
            obs_radec_mpc_2023DW,
            xve=t->auday2kmsec(eph_ea(t/daysec)),
            xvs=t->auday2kmsec(eph_su(t/daysec)),
            xva=t->auday2kmsec(sol1(t/daysec))
        )
        mean_radec1 = mean(res1)
        std_radec1 = std(res1)
        rms_radec1 = nrms(res1,ones(length(res1)))

        @test abs(mean_radec1) ≥ abs(mean_radec0)
        @test std_radec1 ≥ std_radec0
        @test rms_radec1 ≥ rms_radec0
    end

    @testset "Orbit propagation with nongravs: (99942) Apophis" begin

        # integration parameters
        objname = "Apophis"
        maxsteps = 5000
        nyears = 9.0
        dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
        # Initial time [Julian date]
        jd0 = datetime2julian(DateTime(2004,6,1))
        # Initial time [days since J2000]
        t0 = jd0 - PE.J2000
        # JPL #199 solution for Apophis at June 1st, 2004
        q0 = [-1.0506628055913627, -0.06064314196134998, -0.04997102228887035, 0.0029591421121582077, -0.01423233538611057, -0.005218412537773594, -5.592839897872e-14, 0.0]
        # Solar System ephemeris
        sseph = loadpeeph(NEOs.sseph, t0, t0 + nyears*yr)
        # Sun's ephemeris
        eph_su = selecteph(sseph, su)
        # Earth's ephemeris
        eph_ea = selecteph(sseph, ea)

        # warmup propagation
        sol = NEOs.propagate(
            dynamics,
            1,
            jd0,
            nyears,
            q0,
            Val(true),
            order = 25,
            abstol = 1e-20,
            parse_eqs = true
        )

        # propagate orbit
        sol = NEOs.propagate(
            dynamics,
            maxsteps,
            jd0,
            nyears,
            q0,
            Val(true),
            order = 25,
            abstol = 1e-20,
            parse_eqs = true
        )

        # check that solution saves correctly
        jldsave("test.jld2"; sol = sol)
        recovered_sol = JLD2.load("test.jld2", "sol")
        @test sol == recovered_sol
        rm("test.jld2")

        # Read optical astrometry file
        obs_radec_mpc_apophis = NEOs.read_radec_mpc(joinpath("data", "99942_Tholen_etal_2013.dat"))

        # Compute optical astrometry residuals
        res_radec, w_radec = NEOs.residuals(
            obs_radec_mpc_apophis,
            xve=t->auday2kmsec(eph_ea(t/daysec)),
            xvs=t->auday2kmsec(eph_su(t/daysec)),
            xva=t->auday2kmsec(sol(t/daysec))
        )
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
        @test mean_ra ≈ 0.0224 atol=1e-2
        @test std_ra ≈ 0.136 atol=1e-2
        @test rms_ra ≈ std_ra atol=1e-2
        @test mean_dec ≈ -0.0124 atol=1e-2
        @test std_dec ≈ 0.0714 atol=1e-2
        @test rms_dec ≈ std_dec atol=1e-2

        # Read radar astrometry file
        deldop_2005_2013 = NEOs.read_radar_jpl(joinpath("data", "99942_RADAR_2005_2013.dat"))

        # Compute mean radar (time-delay and Doppler-shift) residuals
        @time res_del, w_del, res_dop, w_dop = residuals(
            deldop_2005_2013,
            xve=t->auday2kmsec(eph_ea(t/daysec)),
            xvs=t->auday2kmsec(eph_su(t/daysec)),
            xva=t->auday2kmsec(sol(t/daysec)),
            niter=4,
            tord=5
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
        @test nrms(res, w) ≈ 0.375 atol=1e-2
    end

    @testset "Jet transport propagation and TaylorN serialization" begin

        # Test integration (Apophis)

        # Dynamical function
        local dynamics = RNp1BP_pN_A_J23E_J2S_eph_threads!
        # Order of Taylor polynomials
        local order = 25
        # Absolute tolerance
        local abstol = 1e-20
        # Whether to use @taylorize
        local parse_eqs = true
        # Perturbation to nominal initial condition (Taylor1 jet transport)
        local dq = NEOs.scaled_variables()
        # Initial date of integration (julian days)
        local jd0 = datetime2julian(DateTime(2029, 4, 13, 20))
        # Initial conditions
        local q0 = [-0.9170913888342959, -0.37154308794738056, -0.1610606989484252,
                    0.009701519087787077, -0.012766026792868212, -0.0043488589639194275] .+ dq

        sol = NEOs.propagate(dynamics, 10, jd0, 0.02, q0, Val(true); order, abstol, parse_eqs)
        jldsave("test.jld2"; sol)
        recovered_sol = JLD2.load("test.jld2", "sol")
        @test sol == recovered_sol
        rm("test.jld2")

        sol, tvS, xvS, gvS = NEOs.propagate_root(dynamics, 1, jd0, 0.02, q0, Val(true); order, abstol, parse_eqs)

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

        # It is unlikely that such a short integration generates a non-trivial tvS, xvS and gvS. Therefore, to test
        # VectorTaylorNSerialization I suggest to generate random TaylorN and check it saves correctly...
        local random_TaylorN = [cos(sum(dq .* rand(6))), sin(sum(dq .* rand(6))), tan(sum(dq .* rand(6)))]
        jldsave("test.jld2"; random_TaylorN = random_TaylorN)
        recovered_taylorN = JLD2.load("test.jld2", "random_TaylorN")
        @test recovered_taylorN == random_TaylorN
        rm("test.jld2")

    end

    @testset "Jet transport orbit propagation and astrometric observables: (99942) Apophis" begin

        # integration parameters
        objname::String = "Apophis"
        maxsteps::Int = 2000
        nyears::Float64 = 10.0
        varorder::Int = 1
        dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
        jd0::Float64 = datetime2julian(DateTime(2004,6,1)) #Julian date of integration initial time
        # JPL #199 solution for Apophis at June 1st, 2004
        q00::Vector{Float64} = [-1.0506627941258015, -0.06064313293987095, -0.049970989369473584, 0.002959141747263133, -0.014232335663044254, -0.005218412470120484, -2.789420270048772e-14, 0.0]
        dq::Vector{TaylorN{Float64}} = NEOs.scaled_variables("δx", vcat(fill(1e-8, 6), 1e-14), order = varorder)
        q0::Vector{TaylorN{Float64}} = q00 .+ vcat(dq, 0dq[1])

        # propagate orbit
        sol = NEOs.propagate(
            dynamics,
            maxsteps,
            jd0,
            nyears,
            q0,
            Val(true),
            order = 25,
            abstol = 1e-20,
            parse_eqs = true
        )

        sseph_obs::TaylorInterpolant{Float64,Float64,2} = loadpeeph(NEOs.sseph, sol.t0, sol.t0 + sol.t[end])
        # Sun's ephemeris
        eph_su::TaylorInterpolant{Float64,Float64,2} = selecteph(sseph_obs, su)
        # Earth's ephemeris
        eph_ea::TaylorInterpolant{Float64,Float64,2} = selecteph(sseph_obs, ea)

        # Read optical astrometry file
        obs_radec_mpc_apophis = NEOs.read_radec_mpc(joinpath("data", "99942_Tholen_etal_2013.dat"))

        # Compute optical astrometry residuals
        res_radec, w_radec = NEOs.residuals(
            obs_radec_mpc_apophis,
            xve=t->auday2kmsec(eph_ea(t/daysec)),
            xvs=t->auday2kmsec(eph_su(t/daysec)),
            xva=t->auday2kmsec(sol(t/daysec))
        )
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
        @test mean_ra ≈ 0.0224 atol=1e-2
        @test std_ra ≈ 0.136 atol=1e-2
        @test rms_ra ≈ std_ra atol=1e-2
        @test mean_dec ≈ -0.0124 atol=1e-2
        @test std_dec ≈ 0.0714 atol=1e-2
        @test rms_dec ≈ std_dec atol=1e-2

        # Read radar astrometry file
        deldop_2005_2013 = NEOs.read_radar_jpl(joinpath("data", "99942_RADAR_2005_2013.dat"))

        # Compute mean radar (time-delay and Doppler-shift) residuals
        @time res_del, w_del, res_dop, w_dop = residuals(
            deldop_2005_2013[1:4],
            xve=t->auday2kmsec(eph_ea(t/daysec)),
            xvs=t->auday2kmsec(eph_su(t/daysec)),
            xva=t->auday2kmsec(sol(t/daysec)),
            niter=10,
            tord=10
        )

        @test abs(res_dop[1]()) ≤ deldop_2005_2013[1].Δν_σ
        @test abs(res_del[1]()) ≤ deldop_2005_2013[2].Δτ_σ
        @test abs(res_dop[2]()) ≤ deldop_2005_2013[2].Δν_σ
        @test abs(res_del[2]()) ≤ deldop_2005_2013[2].Δτ_σ
        @test abs(res_dop[3]()) ≤ deldop_2005_2013[3].Δν_σ
        @test abs(res_dop[4]()) ≤ deldop_2005_2013[4].Δν_σ # TODO: fix this residual ("high" residual artifact due to non-optimal initial condition)

        dq_sample = ones(7)
        @test abs(res_dop[1](dq_sample)) ≤ deldop_2005_2013[1].Δν_σ
        @test abs(res_del[1](dq_sample)) ≤ deldop_2005_2013[2].Δτ_σ
        @test abs(res_dop[2](dq_sample)) ≤ deldop_2005_2013[2].Δν_σ
        @test abs(res_del[2](dq_sample)) ≤ deldop_2005_2013[2].Δτ_σ
        @test abs(res_dop[3](dq_sample)) ≤ deldop_2005_2013[3].Δν_σ
        @test abs(res_dop[4](dq_sample)) ≤ deldop_2005_2013[4].Δν_σ
    end

end
