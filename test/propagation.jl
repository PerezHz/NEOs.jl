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

        @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_ng_eph!}, TaylorIntegration.jetcoeffs!))
        @test !isempty(methodswith(Val{RNp1BP_pN_A_J23E_J2S_ng_eph!}, TaylorIntegration._allocate_jetcoeffs!))

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
        dense = true
        quadmath = false # use quadruple precision
        dynamics = RNp1BP_pN_A_J23E_J2S_eph_threads!
        jd0 = datetime2julian(DateTime(2023,2,25,0,0,0)) #Julian date of integration initial time
        # unperturbed initial condition
        q0 = [-9.759018085743707E-01, 3.896554445697074E-01, 1.478066121706831E-01, -9.071450085084557E-03, -9.353197026254517E-03, -5.610023032269034E-03]
        sseph = NEOs.sseph
        # Sun's ephemeris
        eph_su = selecteph(sseph, su)
        # Earth's ephemeris
        eph_ea = selecteph(sseph, ea)

        # warmup propagation
        NEOs.propagate(
            dynamics,
            1,
            jd0,
            nyears,
            sseph,
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
            sseph,
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

        obs_radec_mpc_2023DW = NEOs.read_radec_mpc(joinpath("data", "RADEC_2023_DW.dat"))

        # Compute residuals
        loadjpleph() # load JPL ephemeris
        res, _ = NEOs.residuals(
            obs_radec_mpc_2023DW,
            xve=t->auday2kmsec(eph_ea(t)),
            xvs=t->auday2kmsec(eph_su(t)),
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
            sseph,
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
            xve=t->auday2kmsec(eph_ea(t)),
            xvs=t->auday2kmsec(eph_su(t)),
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
        dense = true
        quadmath = false
        dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph_threads!
        jd0 = datetime2julian(DateTime(2004,6,1)) #Julian date of integration initial time
        # JPL #199 solution for Apophis at June 1st, 2004
        q0 = [-1.0506628055913627, -0.06064314196134998, -0.04997102228887035, 0.0029591421121582077, -0.01423233538611057, -0.005218412537773594, -5.592839897872e-14, 0.0]
        sseph = NEOs.sseph
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
            sseph,
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
            sseph,
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
            xve=t->auday2kmsec(eph_ea(t)),
            xvs=t->auday2kmsec(eph_su(t)),
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
        println("Computing radar astrometric residuals")
        @time res_del, w_del, res_dop, w_dop = residuals(
            deldop_2005_2013,
            xve=t->auday2kmsec(eph_ea(t)),
            xvs=t->auday2kmsec(eph_su(t)),
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
        # Solar System ephemeris
        local sseph = NEOs.sseph
        # Perturbation to nominal initial condition (Taylor1 jet transport)
        local dq = NEOs.scaled_variables()
        # Initial date of integration (julian days)
        local jd0 = datetime2julian(DateTime(2020,12,17))
        # Initial conditions
        local q0 = [-0.18034747703273316, 0.9406910666200128, 0.3457360259054398,
                    -0.016265942170279046, 4.392889725556651e-5, -0.00039519931615139716] .+ dq

        sol = NEOs.propagate(dynamics, 1, jd0, 0.02, sseph, q0, Val(true), order = order, abstol = abstol, parse_eqs = parse_eqs)
        jldsave("test.jld2"; sol = sol)
        recovered_sol = JLD2.load("test.jld2", "sol")
        @test sol == recovered_sol
        rm("test.jld2")

        sol, tvS, xvS, gvS = NEOs.propagate_root(dynamics, 1, jd0, 0.02, sseph, q0, Val(true), order = order, abstol = abstol,
                                                 parse_eqs = parse_eqs)
        jldsave("test.jld2"; sol = sol, tvS = tvS, xvS = xvS, gvS = gvS)
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

end
