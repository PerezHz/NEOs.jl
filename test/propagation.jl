# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using TaylorIntegration
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
        maxsteps = 2000
        nyears = 4.0
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

        # Read optical astrometry file
        obs_radec_mpc_apophis = NEOs.read_radec_mpc(joinpath("data", "99942_Tholen_etal_2013.dat"))

        # Compute residuals
        res, _ = NEOs.residuals(
            obs_radec_mpc_apophis,
            xve=t->auday2kmsec(eph_ea(t)),
            xvs=t->auday2kmsec(eph_su(t)),
            xva=t->auday2kmsec(sol(t/daysec))
        )

        mean_radec = mean(res)
        std_radec = std(res)
        rms_radec = nrms(res,ones(length(res)))
        @test mean_radec ≈ 0.005 atol=1e-2
        @test std_radec ≈ 0.110 atol=1e-2
        @test rms_radec ≈ std_radec atol=1e-2
    end

end
