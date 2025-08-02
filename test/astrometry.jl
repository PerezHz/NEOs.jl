using NEOs
using Test

using NEOs: SCRATCH_PATH

@testset "AbstractAstrometry" begin

    @testset "Utils" begin

        # This tests are based on the examples in:
        # https://minorplanetcenter.net/iau/info/PackedDes.html

        # (Un)packdesig
        packedids = ["J95X00A", "J95X01L", "J95F13B", "J98SA8Q", "J98SC7V", "J98SG2S",
            "K99AJ3Z", "K08Aa0A", "K07Tf8A"]
        unpackedids = ["1995 XA", "1995 XL1", "1995 FB13", "1998 SQ108", "1998 SV127",
            "1998 SS162", "2099 AZ193", "2008 AA360", "2007 TA418"]

        @test all(unpackdesig.(packedids) .== unpackedids)
        @test all(packdesig.(unpackedids) .== packedids)

        packednums = ["03202", "50000", "A0345", "a0017", "K3289", "~0000",
            "~000z", "~AZaz", "~zzzz", "0003I"]
        unpackednums = ["3202", "50000", "100345", "360017", "203289", "620000",
            "620061", "3140113", "15396335", "3I"]

        @test all(unpacknum.(packednums) .== unpackednums)
        @test all(packnum.(unpackednums) .== packednums)

    end

    @testset "CatalogueMPC" begin

        using NEOs: CATALOGUES_MPC, CatalogueMPC, isdeprecated, isunknown, unknowncat,
            parse_catalogues_mpc, read_catalogues_mpc

        # Check global variable CATALOGUES_MPC[]
        @test allunique(CATALOGUES_MPC[])
        @test issorted(CATALOGUES_MPC[])
        @test isa(CATALOGUES_MPC[], Vector{CatalogueMPC})
        @test count(isdeprecated, CATALOGUES_MPC[]) == 42
        @test count(isunknown, CATALOGUES_MPC[]) == 1

        # Parse CatalogueMPC
        gaia_s = """
        [
            {
                "Value": "Gaia2016",
                "Meaning": "Gaia epoch 2016",
                "Deprecated": "No"
            }
        ]
        """
        gaia_p = parse_catalogues_mpc(gaia_s)
        @test isa(gaia_p, Vector{CatalogueMPC})
        @test isone(length(gaia_p))
        gaia = first(gaia_p)
        @test gaia.code == '6'
        @test gaia.value == "Gaia2016"
        @test gaia.meaning == "Gaia epoch 2016"
        @test !isdeprecated(gaia)
        @test gaia in CATALOGUES_MPC[]

        # Unknown catalogue
        unkcat = unknowncat()
        @test isunknown(unkcat)
        @test !isunknown(gaia)

        # CatalogueMPC equality
        @test unkcat == unkcat
        @test gaia == gaia
        @test gaia != unkcat

        # Read catalogues file
        CATALOGUES_PATH = joinpath(SCRATCH_PATH[], "astCat_photCat.json")
        @test isfile(CATALOGUES_PATH)
        cats = read_catalogues_mpc(CATALOGUES_PATH)
        @test CATALOGUES_MPC[] == cats

        # Update catalogues file
        update_catalogues_mpc()
        @test allunique(CATALOGUES_MPC[])
        @test issorted(CATALOGUES_MPC[])
        @test isa(CATALOGUES_MPC[], Vector{CatalogueMPC})
        @test count(isdeprecated, CATALOGUES_MPC[]) == 42
        @test count(isunknown, CATALOGUES_MPC[]) == 1

        # Search catalogue code and value
        cat = search_catalogue_code('6')
        @test cat == gaia
        cat = search_catalogue_value("Gaia2016")
        @test cat == gaia

    end

    @testset "ObservatoryMPC" begin

        using NEOs: OBSERVATORIES_MPC, ObservatoryMPC, isoptical, isoccultation,
            issatellite, isradar, isroving, isgeocentric, hascoord, isunknown,
            unknownobs, parse_observatories_mpc, read_observatories_mpc

        # Check global variable OBSERVATORIES_MPC[]
        @test allunique(OBSERVATORIES_MPC[])
        @test issorted(OBSERVATORIES_MPC[])
        @test isa(OBSERVATORIES_MPC[], Vector{ObservatoryMPC{Float64}})
        @test count(isoptical, OBSERVATORIES_MPC[]) == length(OBSERVATORIES_MPC[]) - 35
        @test count(isoccultation, OBSERVATORIES_MPC[]) == 2
        @test count(issatellite, OBSERVATORIES_MPC[]) == 19
        @test count(isradar, OBSERVATORIES_MPC[]) == 12
        @test count(isroving, OBSERVATORIES_MPC[]) == 2
        @test count(isgeocentric, OBSERVATORIES_MPC[]) == 2
        @test count(hascoord, OBSERVATORIES_MPC[]) == length(OBSERVATORIES_MPC[]) - 22
        @test count(isunknown, OBSERVATORIES_MPC[]) == 0

        # Parse ObservatoryMPC
        lemmon_s = """
        {
            "G96": {
                "created_at": "Sat, 25 May 2019 00:11:26 GMT",
                "firstdate": "2007-09-14",
                "lastdate": null,
                "longitude": "249.21128",
                "name": "Mt. Lemmon Survey",
                "name_latex": "Mt. Lemmon Survey",
                "name_utf8": "Mt. Lemmon Survey",
                "obscode": "G96",
                "observations_type": "optical",
                "old_names": null,
                "rhocosphi": "0.845107",
                "rhosinphi": "0.533611",
                "short_name": "Mt. Lemmon Survey",
                "updated_at": "Tue, 15 Apr 2025 20:52:50 GMT",
                "uses_two_line_observations": false,
                "web_link": "http://www.lpl.arizona.edu/css/"
            }
        }
        """
        lemmon_p = parse_observatories_mpc(lemmon_s)
        @test isa(lemmon_p, Vector{ObservatoryMPC{Float64}})
        @test isone(length(lemmon_p))
        lemmon = first(lemmon_p)
        @test lemmon.code == "G96"
        @test lemmon.coords == [249.21128, 0.845107, 0.533611]
        @test lemmon.name == "Mt. Lemmon Survey"
        @test lemmon.uses_two_line_observations == false
        @test lemmon.observations_type == "optical"
        @test lemmon in OBSERVATORIES_MPC[]

        # Unknown observatory
        unkobs = unknownobs()
        @test !hascoord(unkobs)
        @test isunknown(unkobs)
        @test !isunknown(lemmon)

        # ObservatoryMPC equality
        @test unkobs == unkobs
        @test lemmon == lemmon
        @test lemmon != unkobs

        # Read observatories file
        OBSERVATORIES_PATH = joinpath(SCRATCH_PATH[], "observatoriesmpc.json")
        @test isfile(OBSERVATORIES_PATH)
        obs = read_observatories_mpc(OBSERVATORIES_PATH)
        @test OBSERVATORIES_MPC[] == obs

        # Update observatories file
        update_observatories_mpc()
        @test allunique(OBSERVATORIES_MPC[])
        @test issorted(OBSERVATORIES_MPC[])
        @test isa(OBSERVATORIES_MPC[], Vector{ObservatoryMPC{Float64}})
        @test count(isoptical, OBSERVATORIES_MPC[]) == length(OBSERVATORIES_MPC[]) - 35
        @test count(isoccultation, OBSERVATORIES_MPC[]) == 2
        @test count(issatellite, OBSERVATORIES_MPC[]) == 19
        @test count(isradar, OBSERVATORIES_MPC[]) == 12
        @test count(isroving, OBSERVATORIES_MPC[]) == 2
        @test count(isgeocentric, OBSERVATORIES_MPC[]) == 2
        @test count(hascoord, OBSERVATORIES_MPC[]) == length(OBSERVATORIES_MPC[]) - 22
        @test count(isunknown, OBSERVATORIES_MPC[]) == 0

        # Search observatory code
        obs = search_observatory_code("G96")
        @test obs == lemmon

        # Types of observations
        obs = search_observatory_code("F51")
        @test hascoord(obs)
        @test isoptical(obs)

        obs = search_observatory_code("244")
        @test hascoord(obs)
        @test isoccultation(obs)
        obs = search_observatory_code("275")
        @test !hascoord(obs)
        @test isoccultation(obs)

        obs = search_observatory_code("274")
        @test !hascoord(obs)
        @test issatellite(obs)

        obs = search_observatory_code("251")
        @test hascoord(obs)
        @test isradar(obs)

        obs = search_observatory_code("270")
        @test !hascoord(obs)
        @test isroving(obs)

        obs = search_observatory_code("500")
        @test hascoord(obs)
        @test isgeocentric(obs)

        # Fetch obervatory information
        obs = fetch_observatory_information("G96")
        @test obs["created_at"] == "Sat, 25 May 2019 00:11:26 GMT"
        @test obs["firstdate"] == "2007-09-14"
        @test isnothing(obs["lastdate"])
        @test obs["longitude"] == "249.21128"
        @test obs["name"] == "Mt. Lemmon Survey"
        @test obs["name_latex"] == "Mt. Lemmon Survey"
        @test obs["name_utf8"] == "Mt. Lemmon Survey"
        @test obs["obscode"] == "G96"
        @test obs["observations_type"] == "optical"
        @test isnothing(obs["old_names"])
        @test obs["rhocosphi"] == "0.845107"
        @test obs["rhosinphi"] == "0.533611"
        @test obs["short_name"] == "Mt. Lemmon Survey"
        @test obs["updated_at"] == "Tue, 15 Apr 2025 20:52:50 GMT"
        @test obs["uses_two_line_observations"] == false
        @test obs["web_link"] == "http://www.lpl.arizona.edu/css/"

    end

end