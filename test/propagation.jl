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

    @testset "Orbit propagation without nongravs: 2023 DW" begin

        using PlanetaryEphemeris: selecteph, ea, su, daysec, auday2kmsec
        using Statistics

        objname = "2023DW"
        maxsteps = 1000
        nyears = 0.2
        dense = true
        quadmath = false # use quadruple precision
        dynamics = RNp1BP_pN_A_J23E_J2S_eph_threads!
        jd0 = datetime2julian(DateTime(2023,2,25,0,0,0)) #Julian date of integration initial time
        q0 = [-9.759018085743707E-01, 3.896554445697074E-01, 1.478066121706831E-01, -9.071450085084557E-03, -9.353197026254517E-03, -5.610023032269034E-03]
        sseph = NEOs.sseph
        # Sun's ephemeris
        eph_su = selecteph(sseph, su)
        # Earth's ephemeris
        eph_ea = selecteph(sseph, ea)

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
        nobs = length(obs_radec_mpc_2023DW)

        # Compute residuals
        loadjpleph() # load JPL ephemeris
        @time res, w = NEOs.residuals(
            obs_radec_mpc_2023DW,
            xve=t->auday2kmsec(eph_ea(t)),
            xvs=t->auday2kmsec(eph_su(t)),
            xva=t->auday2kmsec(sol(t/daysec))
        )

        @test mean(res) ≈ -0.6667868924169593
        @test std(res) ≈ 0.7360918596772479
        # pre-fit normalized RMS using Veres et al. (2017) weights
        @test nrms(res, w) ≈ 1.656629515975594
        # pre-fit normalized RMS using ESA/NEOCC/NEODyS weights
        rms_NEOCC = [0.600,0.600,0.600,0.600,0.600,0.600,0.600,0.600,0.600,0.600,0.600,0.600,0.600,
            0.600,0.600,0.600,0.600,0.600,0.600,1.000,1.000,1.000,0.600,0.600,0.600,0.600,1.000,
            0.600,0.600,0.600,0.600,0.600,0.600,1.000,1.000,0.600,0.600,0.600,0.600,0.500,0.500,
            0.600,0.600,0.600,0.600,0.600,1.000,1.000,1.000,0.500,0.500,0.500,0.500,0.500,0.500,
            0.130,0.150,0.160,0.500,0.500,0.600,0.600,0.400,0.400,0.400,0.600,0.600,0.600,0.600,
            0.600,0.500,0.500,0.400,0.400,0.400,0.600,0.600,0.600,0.600,0.600,0.400,0.400,0.400,
            1.000,1.000,0.300,1.000,0.300,1.000,0.300,0.500,0.300,0.500,0.500,0.600,0.600,0.600,
            0.600,0.600,0.500,0.500,0.400,0.400,0.600,0.400,0.600,0.060,0.090,0.400,0.600,0.600,
            0.600,0.600,0.600,0.600,0.100,0.100,0.100,0.060,0.050,0.050,0.060,0.050]
        w_NEOCC = repeat(1.0./rms_NEOCC.^2, 2)
        @test nrms(res, w_NEOCC) ≈ 3.5579681705883575

    end

end
