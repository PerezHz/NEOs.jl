# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Test

@testset "Observations" begin

    @testset "CatalogueMPC" begin

        using NEOs: mpc_catalogue_regex, CatalogueMPC, CatalogueCodes_path, observations_path

        # Check that local catalogues file exists
        @test isfile(joinpath(observations_path, "CatalogueCodes.txt"))

        # Check global variable NEOs.mpc_catalogues[]
        @test allunique(NEOs.mpc_catalogues[])
        @test isa(NEOs.mpc_catalogues[], Vector{CatalogueMPC})

        # Parse CatalogueMPC
        gaia_s = "  6    Gaia2016"
        gaia_m = match(mpc_catalogue_regex, gaia_s)
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
        source_cat = read_catalogues_mpc(CatalogueCodes_path)
        check_file = joinpath(observations_path, "CatalogueCodes_.txt")
        write_catalogues_mpc(source_cat, check_file)
        check_cat = read_catalogues_mpc(check_file)
        rm(check_file)
        @test source_cat == check_cat

        # Update catalogues file
        update_catalogues_mpc(true)
        @test allunique(NEOs.mpc_catalogues[])
        @test isa(NEOs.mpc_catalogues[], Vector{CatalogueMPC})

        # Search catalogue code
        cat = search_cat_code("6")
        @test cat == gaia
    end

    @testset "ObservatoryMPC" begin

        using NEOs: mpc_observatory_regex, ObservatoryMPC, ObsCodes_path

        # Check that local observatories file exists
        @test isfile(joinpath(observations_path, "ObsCodes.txt"))

        # Check global variable NEOs.mpc_observatories[]
        @test allunique(NEOs.mpc_observatories[])
        @test isa(NEOs.mpc_observatories[], Vector{ObservatoryMPC{Float64}})

        # Parse ObservatoryMPC
        arecibo_s = "251 293.246920.949577+0.312734Arecibo"
        arecibo_m = match(mpc_observatory_regex, arecibo_s)
        arecibo = ObservatoryMPC(arecibo_m)
        @test arecibo.code == "251"
        @test arecibo.long == 293.24692
        @test arecibo.cos == 0.949577
        @test arecibo.sin == +0.312734
        @test arecibo.name == "Arecibo"

        hubble_s = "250                           Hubble Space Telescope"
        hubble_m = match(mpc_observatory_regex, hubble_s)
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
        source_obs = read_observatories_mpc(ObsCodes_path)
        check_file = joinpath(observations_path, "ObsCodes_.txt")
        write_observatories_mpc(source_obs, check_file)
        check_obs = read_observatories_mpc(check_file)
        rm(check_file)
        @test source_obs == check_obs

        # Update observatories file
        update_observatories_mpc(true)
        @test allunique(NEOs.mpc_observatories[])
        @test isa(NEOs.mpc_observatories[], Vector{ObservatoryMPC{Float64}})

        # Search observatory code
        obs = search_obs_code("250")
        @test obs == hubble
        obs = search_obs_code("251")
        @test obs == arecibo
    end

    @testset "RadecMPC" begin

        using NEOs: mpc_radec_regex, RadecMPC, src_path
        using Dates

        # Parse RadecMPC
        apophis_s = "99942K04M04N  C2004 03 15.10789 04 06 08.08 +16 55 04.6                om6394691"
        apophis_m = match(mpc_radec_regex, apophis_s)
        apophis = RadecMPC(apophis_m)
        @test apophis.num == "99942"
        @test apophis.tmpdesig == "K04M04N"
        @test apophis.discovery == " "
        @test apophis.publishnote == " "
        @test apophis.obstech == "C"
        @test apophis.date == DateTime("2004-03-15T02:35:21.696")
        @test apophis.α == 1.0739650841580173
        @test apophis.δ == 0.2952738332250385
        @test apophis.info1 == "         "
        @test apophis.mag == "     "
        @test apophis.band == " "
        @test apophis.catalogue == search_cat_code("o")
        @test apophis.info2 == "m6394"
        @test apophis.observatory == search_obs_code("691")

        # RadecMPC equality
        @test apophis == apophis

        # Read/write radec file
        source_file = joinpath(dirname(src_path), "data/99942.dat")
        source_radec = read_radec_mpc(source_file)
        check_file = joinpath(dirname(src_path), "data/99942_.dat")
        write_radec_mpc(source_radec, check_file)
        check_radec = read_radec_mpc(check_file)
        rm(check_file)

        @test source_radec == check_radec

        # Search MPC circulars
        function search_apophis(m::RegexMatch)
            if (m["num"] == "99942") || (m["tmpdesig"] == "N00hp15")
                return true
            else
                return false
            end
        end

        obs = search_circulars_mpc(
            search_apophis,
            "https://minorplanetcenter.net/mpec/K20/K20YA9.html",
            "https://minorplanetcenter.net/mpec/K21/K21JL0.html";
            max_iter = 10
        )

        @test isa(obs, Vector{RadecMPC{Float64}})
        @test issorted(obs)
        @test allunique(obs)
        @test all( (getfield.(obs, :num) .== "99942") .| (getfield.(obs, :tmpdesig) .== "N00hp15") )
    end

    @testset "RadarJPL" begin

        using NEOs: jpl_radar_regex, RadarJPL, jpl_radar_dateformat

        # Parse RadarJPL
        apophis_s = "99942 Apophis (2004 MN4)	2005-01-27 23:31:00	-100849.1434	0.250	Hz	2380	251	251	C"
        apophis_m = match(jpl_radar_regex, apophis_s)
        apophis = RadarJPL(Val(false), apophis_m)
        @test apophis.id == "99942 Apophis (2004 MN4)"
        @test apophis.date == DateTime("2005-01-27 23:31:00", jpl_radar_dateformat)
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

end
