# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Dates
using Test

using NEOs: src_path

@testset "Observations" begin

    @testset "CatalogueMPC" begin

        using NEOs: CATALOGUE_MPC_REGEX, CatalogueMPC, isunknown

        # Check global variable NEOs.CATALOGUES_MPC[]
        @test allunique(NEOs.CATALOGUES_MPC[])
        @test isa(NEOs.CATALOGUES_MPC[], Vector{CatalogueMPC})

        # Parse CatalogueMPC
        gaia_s = "  6    Gaia2016"
        gaia_m = match(CATALOGUE_MPC_REGEX, gaia_s)
        gaia = CatalogueMPC(gaia_m)
        @test gaia.code == "6"
        @test gaia.name == "Gaia2016"

        # Unknown catalogue
        unkcat = unknowncat()
        @test isunknown(unkcat)
        @test !isunknown(gaia)

        # Catalogue equality
        @test unkcat == unkcat
        @test gaia == gaia
        @test gaia != unkcat

        # Read/write catalogues file
        check_file = joinpath(dirname(src_path), "test", "data", "CatalogueCodes.txt")
        write_catalogues_mpc(NEOs.CATALOGUES_MPC[], check_file)
        check_cat = read_catalogues_mpc(check_file)
        rm(check_file)
        @test NEOs.CATALOGUES_MPC[] == check_cat

        # Update catalogues file
        update_catalogues_mpc()
        @test allunique(NEOs.CATALOGUES_MPC[])
        @test isa(NEOs.CATALOGUES_MPC[], Vector{CatalogueMPC})

        # Search catalogue code
        cat = search_cat_code("6")
        @test cat == gaia
    end

    @testset "ObservatoryMPC" begin

        using NEOs: OBSERVATORY_MPC_REGEX, ObservatoryMPC, isunknown

        # Check global variable NEOs.OBSERVATORIES_MPC[]
        @test allunique(NEOs.OBSERVATORIES_MPC[])
        @test isa(NEOs.OBSERVATORIES_MPC[], Vector{ObservatoryMPC{Float64}})

        # Parse ObservatoryMPC
        arecibo_s = "251 293.246920.949577+0.312734Arecibo"
        arecibo_m = match(OBSERVATORY_MPC_REGEX, arecibo_s)
        arecibo = ObservatoryMPC(arecibo_m)
        @test arecibo.code == "251"
        @test arecibo.long == 293.24692
        @test arecibo.cos == 0.949577
        @test arecibo.sin == +0.312734
        @test arecibo.name == "Arecibo"

        hubble_s = "250                           Hubble Space Telescope"
        hubble_m = match(OBSERVATORY_MPC_REGEX, hubble_s)
        hubble = ObservatoryMPC(hubble_m)
        @test hubble.code == "250"
        @test isnan(hubble.long)
        @test isnan(hubble.cos)
        @test isnan(hubble.sin)
        @test hubble.name == "Hubble Space Telescope"

        # Unknown observatory
        unkobs = unknownobs()
        @test isunknown(unkobs)
        @test !isunknown(arecibo)
        @test !isunknown(hubble)

        # hascoord
        @test hascoord(arecibo)
        @test !hascoord(hubble)
        @test !hascoord(unkobs)

        # Catalogue equality
        @test unkobs == unkobs
        @test arecibo == arecibo
        @test hubble == hubble
        @test arecibo != unkobs
        @test hubble != unkobs
        @test arecibo != hubble

        # Read/write observatories file
        check_file = joinpath(dirname(src_path), "test", "data", "ObsCodes.txt")
        write_observatories_mpc(NEOs.OBSERVATORIES_MPC[], check_file)
        check_obs = read_observatories_mpc(check_file)
        rm(check_file)
        @test NEOs.OBSERVATORIES_MPC[] == check_obs

        # Update observatories file
        update_observatories_mpc()
        @test allunique(NEOs.OBSERVATORIES_MPC[])
        @test isa(NEOs.OBSERVATORIES_MPC[], Vector{ObservatoryMPC{Float64}})

        # Search observatory code
        obs = search_obs_code("250")
        @test obs == hubble
        obs = search_obs_code("251")
        @test obs == arecibo
    end

    @testset "RadecMPC" begin

        using NEOs: RADEC_MPC_REGEX, RadecMPC
        using Dates

        # Parse RadecMPC
        apophis_s = "99942K04M04N  C2004 03 15.10789 04 06 08.08 +16 55 04.6                om6394691"
        apophis_m = match(RADEC_MPC_REGEX, apophis_s)
        apophis = RadecMPC(apophis_m)
        @test apophis.num == "99942"
        @test apophis.tmpdesig == "K04M04N"
        @test apophis.discovery == ""
        @test apophis.publishnote == ""
        @test apophis.obstech == "C"
        @test apophis.date == DateTime("2004-03-15T02:35:21.696")
        @test apophis.α == 1.0739650841580173
        @test apophis.δ == 0.2952738332250385
        @test apophis.info1 == ""
        @test isnan(apophis.mag)
        @test apophis.band == ""
        @test apophis.catalogue == search_cat_code("o")
        @test apophis.info2 == "m6394"
        @test apophis.observatory == search_obs_code("691")

        # RadecMPC equality
        @test apophis == apophis

        # Read/write radec file
        source_file = joinpath(pkgdir(NEOs), "test", "data", "RADEC_2023_DW.dat")
        source_radec = read_radec_mpc(source_file)

        @test isa(source_radec, Vector{RadecMPC{Float64}})
        @test issorted(source_radec)
        @test allunique(source_radec)
        @test all( length.(string.(source_radec)) .== 80)

        check_file = joinpath(pkgdir(NEOs), "test", "data", "RADEC_2023_DW_.dat")
        write_radec_mpc(source_radec, check_file)
        check_radec = read_radec_mpc(check_file)
        rm(check_file)

        @test source_radec == check_radec

        # Get RadecMPC
        source_file = joinpath(pkgdir(NEOs), "test", "data", "99942.txt")
        get_radec_mpc("99942", source_file)

        @test isfile(source_file)

        source_radec = read_radec_mpc(source_file)
        rm(source_file)

        @test isa(source_radec, Vector{RadecMPC{Float64}})
        @test issorted(source_radec)
        @test allunique(source_radec)
        @test all( map(x -> length(string(x)) ∈ [80, 161], source_radec))

        check_file = joinpath(pkgdir(NEOs), "test", "data", "99942_.txt")
        write_radec_mpc(source_radec, check_file)
        check_radec = read_radec_mpc(check_file)
        rm(check_file)

        @test source_radec == check_radec
    end

    @testset "RadarJPL" begin

        using NEOs: RADAR_JPL_REGEX, RadarJPL, RADAR_JPL_DATEFORMAT, ismonostatic

        # Parse RadarJPL
        apophis_s = "99942 Apophis (2004 MN4)	2005-01-27 23:31:00	-100849.1434	0.250	Hz	2380	251	251	C"
        apophis_m = match(RADAR_JPL_REGEX, apophis_s)
        apophis = RadarJPL(Val(false), apophis_m)
        @test apophis.id == "99942 Apophis (2004 MN4)"
        @test apophis.date == DateTime("2005-01-27 23:31:00", RADAR_JPL_DATEFORMAT)
        @test isnan(apophis.Δτ)
        @test isnan(apophis.Δτ_σ)
        @test apophis.Δτ_units == ""
        @test apophis.Δν == -100849.1434
        @test apophis.Δν_σ == 0.250
        @test apophis.Δν_units == "Hz"
        @test apophis.freq == 2380.0
        @test apophis.rcvr == search_obs_code("251")
        @test apophis.xmit == search_obs_code("251")
        @test apophis.bouncepoint == "C"
        @test ismonostatic(apophis)
        @test !hasdelay(apophis)
        @test hasdoppler(apophis)

        # RadarJPL equality
        @test apophis == apophis

        # Read/write radar file
        source_file = joinpath(dirname(src_path), "data/99942_RADAR_2005_2013.dat")
        source_radar = read_radar_jpl(source_file)
        check_file = joinpath(dirname(src_path), "data/99942_RADAR_2005_2013_.dat")
        write_radar_jpl(source_radar, check_file)
        check_radar = read_radar_jpl(check_file)
        rm(check_file)

        @test source_radar == check_radar
    end

    @testset "NEOCPObject" begin
        using NEOs: NEOCPObject

        # Fetch current NEOCP list
        objects = fetch_objects_neocp()
        @test isa(objects, Vector{NEOCPObject{Float64}})
        # Pick the first object
        object = objects[1]
        # Check NEOCPObject fields
        @test 0 <= object.score <= 100
        @test object.date >= today() - Month(3)
        @test 0 <= object.α <= 2π
        @test -π <= object.δ <= π
        @test object.nobs > 0
        @test object.arc >= 0
        @test object.notseen > 0
        # Check optical astrometry
        filename = get_radec_neocp(object.tmpdesig)
        @test isfile(filename)
        radec1 = read_radec_mpc(filename)
        radec2 = fetch_radec_neocp(object.tmpdesig)
        rm(filename)
        @test length(radec1) == length(radec2) >= object.nobs
        @test radec1 == radec2
        notseen = (now(UTC) - date(radec1[end])).value / 86_400_000
        @test round(notseen, RoundUp, digits = 3) >= object.notseen
        # Check orbits
        filename = get_orbits_neocp(object.tmpdesig)
        @test isfile(filename)
        rm(filename)
    end

    @testset "Time conversions" begin
        t0 = now()
        @test et2dtutc(dtutc2et(t0)) == t0
        @test jdtdb2dtutc(dtutc2jdtdb(t0)) == t0
    end

    @testset "Topocentric" begin
        using NEOs: TimeOfDay, sunriseset, obsposECEF, obsposvelECI

        # Ground observation
        radec_1 = read_radec_mpc("""
        99942        |C2012 12 12.33230011 28 40.300-26 29 32.10         17.70Vu~0mfl807
        99942        |C2012 12 12.33730011 28 38.970-26 29 34.80         17.60Vu~0mfl807
        99942        |C2012 12 12.34221011 28 37.640-26 29 37.50         17.50Vu~0mfl807
        99942        |C2012 12 12.34712011 28 36.330-26 29 40.00         17.50Vu~0mfl807
        99942        |C2012 12 12.35054011 28 35.400-26 29 41.80         17.50Vu~0mfl807
        """)
        # Sattellite observation
        radec_2 = read_radec_mpc("""
        99942         S2020 12 18.97667011 30 15.530-10 46 20.20         19.00RL~4ROFC51
        99942         s2020 12 18.9766701 - 5634.1734 - 2466.2657 - 3038.3924   ~4ROFC51
        99942         S2020 12 19.10732011 30 22.510-10 48 20.00               L~4ROFC51
        99942         s2020 12 19.1073201 - 5654.1816 - 2501.9465 - 2971.1902   ~4ROFC51
        99942         S2020 12 19.23810011 30 29.500-10 50 19.60               L~4ROFC51
        99942         s2020 12 19.2381001 - 5645.7831 - 2512.1036 - 2978.6411   ~4ROFC51
        99942         S2020 12 19.23822011 30 29.570-10 50 19.20               L~4ROFC51
        99942         s2020 12 19.2382201 - 5617.3465 - 2486.4031 - 3053.2209   ~4ROFC51
        """)

        # Check parsing
        @test length(radec_1) == 5
        @test all( map(x -> x.observatory.code, radec_1) .== "807")
        @test length(radec_2) == 4
        @test all( map(x -> x.observatory.code, radec_2) .== "C51")

        # TimeOfDay
        tod_1 = TimeOfDay.(radec_1)
        tod_2 = TimeOfDay.(radec_2)
        # Check
        @test allequal(tod_1)
        @test tod_1[1].light == :night
        @test tod_1[1].start == Date(2012, 12, 11)
        @test tod_1[1].stop == Date(2012, 12, 12)
        @test tod_1[1].utc == -5
        @test allunique(tod_2)
        @test all( getfield.(tod_2, :light) .== :space )
        @test all( date.(radec_2) .== getfield.(tod_2, :start) .== getfield.(tod_2, :start) )
        @test all( getfield.(tod_2, :utc) .== 0 )

        # Sunrise and sunset
        radec = read_radec_mpc("99942        8C2020 12 08.15001011 20 07.510-08 02 54.20         18.50GV~4ROF094")
        sun = sunriseset(radec[1])
        @test dtutc2jdtdb(sun[1]) ≈ dtutc2jdtdb(DateTime("2020-12-08T05:05:59.384"))
        @test dtutc2jdtdb(sun[2]) ≈ dtutc2jdtdb(DateTime("2020-12-08T14:05:49.386"))

        # obsposECEF
        ecef_2 = obsposECEF.(radec_2)
        @test ecef_2[1] ≈ [-3462.643557087632, 5076.197661798687, -3049.6756672719907]
        @test ecef_2[2] ≈ [1351.315736765706, 6027.937408384214, -2982.5146167937583]
        @test ecef_2[3] ≈ [5332.067839021762, 3112.403799578623, -2989.9547254809945]
        @test ecef_2[4] ≈ [5308.786202404402, 3079.725220466387, -3064.4773721684687]

        # obsposvelECI
        eci_2 = obsposvelECI.(radec_2)
        @test eci_2[1] == [-5634.1734, -2466.2657, -3038.3924, 0.0, 0.0, 0.0]
        @test eci_2[2] == [-5654.1816, -2501.9465, -2971.1902, 0.0, 0.0, 0.0]
        @test eci_2[3] == [-5645.7831, -2512.1036, -2978.6411, 0.0, 0.0, 0.0]
        @test eci_2[4] == [-5617.3465, -2486.4031, -3053.2209, 0.0, 0.0, 0.0]

    end

    @testset "Tracklet" begin
        using NEOs: reduce_tracklets

        # Choose this example because of the discontinuity in α

        # Fetch optical astrometry
        radec = fetch_radec_mpc("2020 TJ6")
        # Reduce tracklets
        tracklets = reduce_tracklets(radec)

        # Values by December 19, 2023
        @test length(tracklets) == 5
        @test getfield.(tracklets, :nobs) == [4, 4, 3, 4, 3]
        @test tracklets[3].α ≈ 6.2831 atol = 1e-4
        @test tracklets[3].observatory == search_obs_code("J04")
        @test tracklets[3].night.utc == -1
    end

end
