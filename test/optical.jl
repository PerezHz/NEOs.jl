using NEOs
using Dates
using PlanetaryEphemeris
using LinearAlgebra
using DataFrames
using Query
using Test

using NEOs: AbstractOpticalAstrometry

const TEST_DATA = joinpath(pkgdir(NEOs), "test", "data")

isvalidEros(x::AbstractOpticalAstrometry) =
    Date(2000) < date(x) < Date(2025) && observatorycode(x) != "K73"

@testset "AbstractOpticalAstrometry" begin

    @testset "OpticalMPC80" begin

        using NEOs: OpticalMPC80, parse_optical_mpc80, isoccultation

        # Parse OpticalMPC80
        apophis_s = """
        99942K04M04N  C2004 03 15.10789 04 06 08.08 +16 55 04.6                om6394691
        """
        apophis_p = parse_optical_mpc80(apophis_s)
        @test isa(apophis_p, Vector{OpticalMPC80{Float64}})
        @test isone(length(apophis_p))
        apophis = first(apophis_p)
        @test apophis.number == "99942"
        @test apophis.desig == "K04M04N"
        @test apophis.discovery == ' '
        @test apophis.note1 == ' '
        @test apophis.note2 == 'C'
        @test apophis.date == date(apophis) == DateTime("2004-03-15T02:35:21.696")
        @test apophis.ra == ra(apophis) == 1.0739650841580173
        @test apophis.dec == dec(apophis) == 0.2952738332250385
        @test apophis.info1 == ""
        @test isnan(apophis.mag) && isnan(mag(apophis))
        @test apophis.band == band(apophis) == search_magnitude_band("UNK")
        @test apophis.catalogue == catalogue(apophis) == search_catalogue_code('o')
        @test apophis.info2 == "m6394"
        @test apophis.observatory == observatory(apophis) == search_observatory_code("691")
        @test apophis.source == apophis_s

        @test NEOs.numtype(apophis) == Float64
        @test designation(apophis) == "99942"
        @test measure(apophis) == (1.0739650841580173, 0.2952738332250385)
        @test rms(apophis) == (1.0, 1.0)
        @test debias(apophis) == (0.0, 0.0)
        @test corr(apophis) == 0.0
        @test trackletid(apophis) == ""
        @test !isdiscovery(apophis)
        @test !isdeprecated(apophis)
        @test vconversion(apophis) == 0.0
        @test cataloguecode(apophis) == 'o'
        @test observatorycode(apophis) == "691"

        # Custom show
        @test string(apophis) == "99942 α: 61.53367° δ: 16.91794° t: 2004-03-15T02:35:21.696 \
                                  obs: Steward Observatory, Kitt Peak-Spacewatch"

        # OpticalMPC80 equality
        @test apophis == apophis

        # Fetch OpticalMPC80
        optical1 = fetch_optical_mpc80("433", MPC)
        filter!(isvalidEros, optical1)
        @test isa(optical1, Vector{OpticalMPC80{Float64}})
        @test issorted(optical1)
        @test allunique(optical1)

        # Read/write OpticalMPC80 file
        filename = joinpath(TEST_DATA, "433.txt")
        optical2 = read_optical_mpc80(filename)
        filter!(isvalidEros, optical2)
        @test isa(optical2, Vector{OpticalMPC80{Float64}})
        @test issorted(optical2)
        @test allunique(optical2)

        filename = joinpath(TEST_DATA, "433_.txt")
        write_optical_mpc80(optical2, filename)
        optical3 = read_optical_mpc80(filename)
        rm(filename)
        @test isa(optical3, Vector{OpticalMPC80{Float64}})
        @test issorted(optical3)
        @test allunique(optical3)

        @test optical1 == optical2 == optical3

        # DataFrames
        df1 = DataFrame(optical1)
        df2 = DataFrame(optical2)
        df3 = DataFrame(optical3)

        @test nrow(df1) == length(optical1)
        @test nrow(df2) == length(optical2)
        @test nrow(df3) == length(optical3)

        @test all(names(df1) .== String.(fieldnames(OpticalMPC80{Float64})))
        @test all(names(df2) .== String.(fieldnames(OpticalMPC80{Float64})))
        @test all(names(df3) .== String.(fieldnames(OpticalMPC80{Float64})))

        @test all(eltype.(eachcol(df1)) .== fieldtypes(OpticalMPC80{Float64}))
        @test all(eltype.(eachcol(df2)) .== fieldtypes(OpticalMPC80{Float64}))
        @test all(eltype.(eachcol(df3)) .== fieldtypes(OpticalMPC80{Float64}))

        # Query
        optical4 = optical1 |> @filter(year(date(_)) > 2011 &&
            isoccultation(observatory(_))) |> DataFrame

        @test nrow(optical4) == 3
        @test all(@. year(optical4.date) > 2011)
        @test all(@. isoccultation(optical4.observatory))

    end

    @testset "NEOCPObject" begin

        using NEOs: NEOCPObject, parse_neocp_objects, isunknown

        # Parse NEOCPObject
        X85177_s = """
        X85177  100 2025 05 24.6  15.6718  -6.9766 21.5 Updated June 13.68 UT            3   0.03 19.1 20.282
        """
        X85177_p = parse_neocp_objects(X85177_s)
        @test isa(X85177_p, Vector{NEOCPObject{Float64}})
        @test isone(length(X85177_p))
        X85177 = first(X85177_p)
        @test NEOs.numtype(X85177) == Float64
        @test X85177.desig == "X85177"
        @test X85177.score == 100
        @test X85177.date == date(X85177) == DateTime("2025-05-24T14:24:00")
        @test X85177.ra == ra(X85177) == 4.10286764571071
        @test X85177.dec == dec(X85177) == -0.1217646405946364
        @test X85177.V == mag(X85177) == 21.5
        @test band(X85177) == ' '
        @test X85177.updated == "Updated June 13.68 UT"
        @test X85177.nobs == 3
        @test X85177.arc == 0.03
        @test X85177.H == 19.1
        @test X85177.notseen == 20.282
        @test X85177.source == X85177_s[1:end-1]
        @test measure(X85177) == (4.10286764571071, -0.1217646405946364)
        @test isunknown(observatory(X85177))
        @test isunknown(catalogue(X85177))
        @test all(isnan, rms(X85177))
        @test all(isnan, debias(X85177))
        @test isnan(corr(X85177))
        @test trackletid(X85177) == ""

        # Custom show
        @test string(X85177) == "X85177 α: 235.07700° δ: -6.97660° t: 2025-05-24T14:24:00"

        # NEOCPObject equality
        @test X85177 == X85177

        # Fetch NEOCPObject
        # Note: the tests below have been commented as they break consistently
        # in Github because the MPC's website rejects the HTTP request needed
        # to download the NEOCP file (16/12/2025)
        # objects1 = fetch_neocp_objects()
        # @test isa(objects1, Vector{NEOCPObject{Float64}})
        # @test issorted(objects1)
        # @test allunique(objects1)

        # Read/write NEOCPObject file
        filename = joinpath(TEST_DATA, "neocp.txt")
        objects2 = read_neocp_objects(filename)
        @test isa(objects2, Vector{NEOCPObject{Float64}})
        @test issorted(objects2)
        @test allunique(objects2)

        filename = joinpath(TEST_DATA, "neocp_.txt")
        write_neocp_objects(objects2, filename)
        objects3 = read_neocp_objects(filename)
        rm(filename)
        @test isa(objects3, Vector{NEOCPObject{Float64}})
        @test issorted(objects3)
        @test allunique(objects3)

        @test objects2 == objects3

    end

    @testset "OpticalRWO" begin

        using NEOs: NEODyS2_OPTICAL_HEADER, OpticalRWO, parse_optical_rwo, isoccultation

        # Parse OpticalRWO
        apophis_s = """
        version =   3
        errmod  = 'vfcc17'
        RMSast  =   3.02422E-01
        RMSmag  =   4.31723E-01
        END_OF_HEADER
        ! Object   Obser ============= Date ============= ================== Right Ascension =================  ================= Declination ===================== ==== Magnitude ==== Ast Obs  Residual SEL
        ! Design   K T N YYYY MM DD.dddddddddd   Accuracy HH MM SS.sss  Accuracy      RMS  F     Bias    Resid sDD MM SS.ss  Accuracy      RMS  F     Bias    Resid Val  B   RMS  Resid Cat Cod       Chi A M
         99942     O C   2004 03 15.10789       1.000E-05 04 06 08.080  1.435E-01    0.612 F   -0.247    0.251 +16 55 04.60  1.000E-01    0.612 F    0.140   -0.070                       o 691      0.43 1 0
        """
        apophis_p = parse_optical_rwo(apophis_s)
        @test isa(apophis_p, Vector{OpticalRWO{Float64}})
        @test isone(length(apophis_p))
        apophis = first(apophis_p)
        @test apophis.design == "99942"
        @test apophis.K == 'O'
        @test apophis.T == 'C'
        @test apophis.N == ' '
        @test apophis.date == date(apophis) == DateTime("2004-03-15T02:35:21.696")
        @test apophis.date_accuracy == 1e-5
        @test apophis.ra == ra(apophis) == 1.0739650841580173
        @test apophis.ra_accuracy == 0.1435
        @test apophis.ra_rms == 0.612
        @test apophis.ra_flag == false
        @test apophis.ra_bias == -0.247
        @test apophis.ra_resid == 0.251
        @test apophis.dec == dec(apophis) == 0.2952738332250385
        @test apophis.dec_accuracy == 0.1
        @test apophis.dec_rms == 0.612
        @test apophis.dec_flag == false
        @test apophis.dec_bias == 0.14
        @test apophis.dec_resid == -0.07
        @test isnan(apophis.mag) && isnan(mag(apophis))
        @test apophis.mag_band == band(apophis) == search_magnitude_band("UNK")
        @test isnan(apophis.mag_rms)
        @test isnan(apophis.mag_resid)
        @test apophis.catalogue == catalogue(apophis) == search_catalogue_code('o')
        @test apophis.observatory == observatory(apophis) == search_observatory_code("691")
        @test apophis.chi == 0.43
        @test apophis.sel_A == 1
        @test apophis.sel_M == 0
        idxs = findfirst(NEODyS2_OPTICAL_HEADER, apophis_s)
        @test apophis.source == apophis_s[last(idxs)+1:end]
        idxs = findfirst("END_OF_HEADER", apophis_s)
        @test apophis.header == apophis_s[1:first(idxs)-1]

        @test NEOs.numtype(apophis) == Float64
        @test designation(apophis) == "99942"
        @test measure(apophis) == (1.0739650841580173, 0.2952738332250385)
        @test rms(apophis) == (0.612, 0.612)
        @test debias(apophis) == (-0.247, 0.14)
        @test corr(apophis) == 0.0
        @test trackletid(apophis) == ""
        @test !isdeprecated(apophis)
        @test vconversion(apophis) == 0.0
        @test cataloguecode(apophis) == 'o'
        @test observatorycode(apophis) == "691"

        # Custom show
        string(apophis) == "99942 α: 61.53367° δ: 16.91794° t: 2004-03-15T02:35:21.696 \
                            obs: Steward Observatory, Kitt Peak-Spacewatch"

        # OpticalRWO equality
        @test apophis == apophis

        # Fetch OpticalRWO
        optical1 = fetch_optical_rwo("433", NEOCC)
        filter!(isvalidEros, optical1)
        @test isa(optical1, Vector{OpticalRWO{Float64}})
        @test issorted(optical1)
        @test allunique(optical1)

        optical2 = fetch_optical_rwo("433", NEODyS2)
        filter!(isvalidEros, optical2)
        @test isa(optical2, Vector{OpticalRWO{Float64}})
        @test issorted(optical2)
        @test allunique(optical2)

        # Read/write OpticalRWO file
        filename = joinpath(TEST_DATA, "433.rwo")
        optical3 = read_optical_rwo(filename)
        filter!(isvalidEros, optical3)
        @test isa(optical3, Vector{OpticalRWO{Float64}})
        @test issorted(optical3)
        @test allunique(optical3)

        filename = joinpath(TEST_DATA, "433_.rwo")
        write_optical_rwo(optical3, filename)
        optical4 = read_optical_rwo(filename)
        rm(filename)
        @test isa(optical4, Vector{OpticalRWO{Float64}})
        @test issorted(optical4)
        @test allunique(optical4)

        @test optical1 == #= optical2 == =# optical3 == optical4

        # DataFrames
        df1 = DataFrame(optical1)
        df2 = DataFrame(optical2)
        df3 = DataFrame(optical3)
        df4 = DataFrame(optical4)

        @test nrow(df1) == length(optical1)
        @test nrow(df2) == length(optical2)
        @test nrow(df3) == length(optical3)
        @test nrow(df4) == length(optical4)

        @test all(names(df1) .== String.(fieldnames(OpticalRWO{Float64})))
        @test all(names(df2) .== String.(fieldnames(OpticalRWO{Float64})))
        @test all(names(df3) .== String.(fieldnames(OpticalRWO{Float64})))
        @test all(names(df4) .== String.(fieldnames(OpticalRWO{Float64})))

        @test all(eltype.(eachcol(df1)) .== fieldtypes(OpticalRWO{Float64}))
        @test all(eltype.(eachcol(df2)) .== fieldtypes(OpticalRWO{Float64}))
        @test all(eltype.(eachcol(df3)) .== fieldtypes(OpticalRWO{Float64}))
        @test all(eltype.(eachcol(df4)) .== fieldtypes(OpticalRWO{Float64}))

        # Query
        optical5 = optical1 |> @filter(year(date(_)) > 2011 &&
            isoccultation(observatory(_))) |> DataFrame

        @test nrow(optical5) == 3
        @test all(@. year(optical5.date) > 2011)
        @test all(@. isoccultation(optical5.observatory))

    end

    @testset "OpticalADES" begin

        using NEOs: OpticalADES, parse_optical_ades, unknowncat, isoccultation

        # Parse OpticalADES
        apophis_s = """
        <?xml version="1.0" encoding="UTF-8"?>
        <ades version="2022">
          <optical>
            <permID>99942</permID>
            <provID>2004 MN4</provID>
            <trkSub>K04M04N</trkSub>
            <obsID>JqcHxe000000DaKP010000001</obsID>
            <trkID>000002w-NJ</trkID>
            <mode>CCD</mode>
            <stn>691</stn>
            <obsTime>2004-03-15T02:35:21.696Z</obsTime>
            <ra>61.53367</ra>
            <dec>16.91794</dec>
            <astCat>USNOB1</astCat>
            <ref>MPS   126394</ref>
            <subFmt>M92</subFmt>
            <precTime>10</precTime>
            <precRA>0.01</precRA>
            <precDec>0.1</precDec>
          </optical>
        </ades>
        """
        apophis_p = parse_optical_ades(apophis_s)
        @test isa(apophis_p, Vector{OpticalADES{Float64}})
        @test isone(length(apophis_p))
        apophis = first(apophis_p)
        @test apophis.permid == "99942"
        @test apophis.provid == "2004 MN4"
        @test apophis.trksub == "K04M04N"
        @test apophis.obsid == "JqcHxe000000DaKP010000001"
        @test apophis.trkid == "000002w-NJ"
        @test apophis.mode == "CCD"
        @test apophis.stn == observatory(apophis) == search_observatory_code("691")
        @test apophis.sys == ""
        @test apophis.ctr == 0
        @test isnan(apophis.pos1)
        @test isnan(apophis.pos2)
        @test isnan(apophis.pos3)
        @test isnan(apophis.vel1)
        @test isnan(apophis.vel2)
        @test isnan(apophis.vel3)
        @test isnan(apophis.poscov11)
        @test isnan(apophis.poscov12)
        @test isnan(apophis.poscov13)
        @test isnan(apophis.poscov22)
        @test isnan(apophis.poscov23)
        @test isnan(apophis.poscov33)
        @test apophis.prog == ""
        @test apophis.obstime == date(apophis) == DateTime("2004-03-15T02:35:21.696")
        @test isnan(apophis.rmstime)
        @test apophis.ra == ra(apophis) == 1.073965142335659
        @test apophis.dec == dec(apophis) == 0.2952737556548495
        @test isnan(apophis.rastar)
        @test isnan(apophis.decstar)
        @test isnan(apophis.deltara)
        @test isnan(apophis.deltadec)
        @test isnan(apophis.rmsra)
        @test isnan(apophis.rmsdec)
        @test isnan(apophis.rmscorr)
        @test apophis.astcat == catalogue(apophis) == search_catalogue_code('o')
        @test isnan(apophis.mag) && isnan(mag(apophis))
        @test isnan(apophis.rmsmag)
        @test apophis.band == band(apophis) == search_magnitude_band("UNK")
        @test apophis.photcat == unknowncat()
        @test apophis.ref == "MPS   126394"
        @test apophis.disc == ""
        @test apophis.subfmt == "M92"
        @test apophis.prectime == 10
        @test apophis.precra == 4.84813681109536e-8
        @test apophis.precdec == 4.84813681109536e-7
        @test isnan(apophis.unctime)
        @test apophis.notes == ""
        @test apophis.remarks == ""
        @test apophis.deprecated == ""
        @test replace(apophis.source, " " => "") == replace(apophis_s[64:end-9], " " => "")

        @test NEOs.numtype(apophis) == Float64
        @test designation(apophis) == "99942"
        @test measure(apophis) == (1.073965142335659, 0.2952737556548495)
        @test rms(apophis) == (1.0, 1.0)
        @test debias(apophis) == (0.0, 0.0)
        @test corr(apophis) == 0.0
        @test trackletid(apophis) == "000002w-NJ"
        @test !isdiscovery(apophis)
        @test !isdeprecated(apophis)
        @test vconversion(apophis) == 0.0
        @test cataloguecode(apophis) == 'o'
        @test observatorycode(apophis) == "691"

        # Custom show
        string(apophis) == "99942 α: 61.53367° δ: 16.91794° t: 2004-03-15T02:35:21.696 \
                            obs: Steward Observatory, Kitt Peak-Spacewatch"

        # OpticalRWO equality
        @test apophis == apophis

        # Fetch OpticalADES
        optical1 = fetch_optical_ades("433", MPC)
        filter!(isvalidEros, optical1)
        @test isa(optical1, Vector{OpticalADES{Float64}})
        @test issorted(optical1)
        @test allunique(optical1)

        # Read/write OpticalRWO file
        filename = joinpath(TEST_DATA, "433.xml")
        optical2 = read_optical_ades(filename)
        filter!(isvalidEros, optical2)
        @test isa(optical2, Vector{OpticalADES{Float64}})
        @test issorted(optical2)
        @test allunique(optical2)

        filename = joinpath(TEST_DATA, "433_.xml")
        write_optical_ades(optical2, filename)
        optical3 = read_optical_ades(filename)
        rm(filename)
        @test isa(optical3, Vector{OpticalADES{Float64}})
        @test issorted(optical3)
        @test allunique(optical3)

        @test optical1 == optical2 == optical3

        # DataFrames
        df1 = DataFrame(optical1)
        df2 = DataFrame(optical2)
        df3 = DataFrame(optical3)

        @test nrow(df1) == length(optical1)
        @test nrow(df2) == length(optical2)
        @test nrow(df3) == length(optical3)

        @test all(names(df1) .== String.(fieldnames(OpticalADES{Float64})))
        @test all(names(df2) .== String.(fieldnames(OpticalADES{Float64})))
        @test all(names(df3) .== String.(fieldnames(OpticalADES{Float64})))

        @test all(eltype.(eachcol(df1)) .== fieldtypes(OpticalADES{Float64}))
        @test all(eltype.(eachcol(df2)) .== fieldtypes(OpticalADES{Float64}))
        @test all(eltype.(eachcol(df3)) .== fieldtypes(OpticalADES{Float64}))

        # Query
        optical4 = optical1 |> @filter(year(date(_)) > 2011 &&
            isoccultation(observatory(_))) |> DataFrame

        @test nrow(optical4) == 3
        @test all(@. year(optical4.obstime) > 2011)
        @test all(@. isoccultation(optical4.stn))

    end

    # Load optical astrometry
    mpc80_file = joinpath(TEST_DATA, "433.txt")
    rwo_file = joinpath(TEST_DATA, "433.rwo")
    ades_file = joinpath(TEST_DATA, "433.xml")

    optical1 = read_optical_mpc80(mpc80_file)
    optical2 = read_optical_rwo(rwo_file)
    optical3 = read_optical_ades(ades_file)

    @testset "load_optical_astrometry" begin
        mpc80_optical, format = load_optical_astrometry(mpc80_file)
        @test mpc80_optical == optical1 && format == "mpc80"
        mpc80_optical, format = load_optical_astrometry(mpc80_file; format = :obs80)
        @test mpc80_optical == optical1 && format == "mpc80"
        # Note: I've commented out the tests below to limit the number of calls
        # made to the MPC Observations API during the tests (@LuEdRaMo, 16/08/26)
        # mpc80_optical, format = load_optical_astrometry("433"; format = :mpc80)
        # @test issubset(optical1, mpc80_optical) && format == "mpc80"

        ades_optical, format = load_optical_astrometry(ades_file)
        @test ades_optical == optical3 && format == "ades"
        ades_optical, format = load_optical_astrometry(ades_file; format = "xml")
        @test ades_optical == optical3 && format == "ades"
        # Note: I've commented out the tests below to limit the number of calls
        # made to the MPC Observations API during the tests (@LuEdRaMo, 16/08/26)
        # ades_optical, format = load_optical_astrometry("433"; format = :ades)
        # @test issubset(optical3, ades_optical) && format == "ades"

        @test_throws ArgumentError load_optical_astrometry(mpc80_file; format = :unknown)
        @test_throws ArgumentError load_optical_astrometry("433"; format = :unknown)
    end

    # Delete observations before J200 and from observatory K73
    filter!(isvalidEros, optical1)
    filter!(isvalidEros, optical2)
    filter!(isvalidEros, optical3)

    @testset "Topocentric" begin

        using NEOs: isday, isnight, issatellite, isgeocentric, sunriseset,
                    obsposECEF, obsposvelECI, parse_optical_ades

        # TimeOfDay
        light1 = timeofday.(optical1)
        light2 = timeofday.(optical2)
        light3 = timeofday.(optical3)

        mask1, mask2, mask3 = @. isday(light1), isday(light2), isday(light3)
        @test count(mask1) == count(mask2) == count(mask3) == 0
        @test light1[mask1] == light2[mask2] == light3[mask3]
        mask1, mask2, mask3 = @. isnight(light1), isnight(light2), isnight(light3)
        @test count(mask1) == count(mask2) == count(mask3) == 6_692
        @test light1[mask1] == light2[mask2] == light3[mask3]
        mask1, mask2, mask3 = @. issatellite(light1), issatellite(light2), issatellite(light3)
        @test count(mask1) == count(mask2) == count(mask3) == 1_790
        # @test light1[mask1] == light2[mask2] == light3[mask3]
        mask1, mask2, mask3 = @. isgeocentric(light1), isgeocentric(light2), isgeocentric(light3)
        @test count(mask1) == count(mask2) == count(mask3) == 0
        @test light1[mask1] == light2[mask2] == light3[mask3]

        @test extrema(getfield.(light1, :utc)) == extrema(getfield.(light2, :utc)) ==
            extrema(getfield.(light3, :utc)) == (-10, 12)

        # Sunrise and sunset
        sun1, sun2, sun3 = sunriseset(optical1[1]), sunriseset(optical2[1]), sunriseset(optical3[1])

        @test abs(datediff(sun1[1], DateTime("2000-01-07T19:00:19.709"))) < 2
        @test abs(datediff(sun2[1], DateTime("2000-01-07T19:00:19.709"))) < 2
        @test abs(datediff(sun3[1], DateTime("2000-01-07T19:00:19.709"))) < 2

        @test abs(datediff(sun1[2], DateTime("2000-01-08T09:03:20.951"))) < 2
        @test abs(datediff(sun2[2], DateTime("2000-01-08T09:03:20.951"))) < 2
        @test abs(datediff(sun3[2], DateTime("2000-01-08T09:03:20.951"))) < 2

        # obsposECEF
        posECEF1 = obsposECEF.(optical1)
        posECEF2 = obsposECEF.(optical2)
        posECEF3 = obsposECEF.(optical3)

        @test maximum(norm, posECEF1 - posECEF2) < 1.4e-10
        @test maximum(norm, posECEF1 - posECEF3) < 1.28
        @test maximum(norm, posECEF2 - posECEF3) < 1.28

        rECEF1, rECEF2, rECEF3 = @. norm(posECEF1), norm(posECEF2), norm(posECEF3)
        mask1, mask2, mask3 = @. isnight(light1), isnight(light2), isnight(light3)
        @test 6_362 ≤ minimum(rECEF1[mask1]) == minimum(rECEF2[mask2]) == minimum(rECEF3[mask3])
        @test maximum(rECEF1[mask1]) == maximum(rECEF2[mask2]) == maximum(rECEF3[mask3]) ≤ 6_3780

        # obsposvelECI
        posvelECI1 = obsposvelECI.(optical1)
        posvelECI2 = obsposvelECI.(optical2)
        posvelECI3 = obsposvelECI.(optical3)

        @test maximum(norm, posvelECI1 - posvelECI2) < 8e-11
        @test maximum(norm, posvelECI1 - posvelECI3) < 0.17
        @test maximum(norm, posvelECI2 - posvelECI3) < 0.17

        rECI1, rECI2, rECI3 = @. norm(getindex(posvelECI1, $Ref(1:3))),
            norm(getindex(posvelECI2, $Ref(1:3))), norm(getindex(posvelECI3, $Ref(1:3)))
        mask1, mask2, mask3 = @. isnight(light1), isnight(light2), isnight(light3)
        @test 6362 ≤ minimum(rECI1[mask1]) ≈ minimum(rECI2[mask2]) ≈ minimum(rECI3[mask3])
        @test maximum(rECI1[mask1]) ≈ maximum(rECI2[mask2]) ≈ maximum(rECI3[mask3]) ≤ 6380

        @test maximum(abs, rECEF1 - rECI1) < 2.4e-10
        @test maximum(abs, rECEF2 - rECI2) < 2.4e-10
        @test maximum(abs, rECEF3 - rECI3) < 2.4e-10

        # Consistency between ICRF_KM and ICRF_AU
        s = """
        <?xml version='1.0' encoding='UTF-8'?>
        <ades version="2022">
        <optical>
            <trkSub>SYNTH</trkSub>
            <trkID>0000000195</trkID>
            <mode>UNK</mode>
            <stn>247</stn>
            <sys>ICRF_KM</sys>
            <ctr>399</ctr>
            <pos1>44560751.403</pos1>
            <pos2>249847769.733</pos2>
            <pos3>63875407.816</pos3>
            <posCov11>1.0</posCov11>
            <posCov12>0.0</posCov12>
            <posCov13>0.0</posCov13>
            <posCov22>1.0</posCov22>
            <posCov23>0.0</posCov23>
            <posCov33>1.0</posCov33>
            <obsTime>2028-04-12T08:16:03.68Z</obsTime>
            <ra>87.975492</ra>
            <dec>25.292972</dec>
            <rmsRA>1</rmsRA>
            <rmsDec>1</rmsDec>
            <astCat>UNK</astCat>
            <mag>10.8</mag>
            <rmsMag>1</rmsMag>
            <band>V</band>
        </optical>
        <optical>
            <trkSub>SYNTH</trkSub>
            <trkID>0000000195</trkID>
            <mode>UNK</mode>
            <stn>247</stn>
            <sys>ICRF_KM</sys>
            <ctr>399</ctr>
            <pos1>43351538.866</pos1>
            <pos2>250481696.410</pos2>
            <pos3>64215174.021</pos3>
            <posCov11>1.0</posCov11>
            <posCov12>0.0</posCov12>
            <posCov13>0.0</posCov13>
            <posCov22>1.0</posCov22>
            <posCov23>0.0</posCov23>
            <posCov33>1.0</posCov33>
            <obsTime>2028-04-12T20:16:03.68Z</obsTime>
            <ra>212.234596</ra>
            <dec>-8.836389</dec>
            <rmsRA>1</rmsRA>
            <rmsDec>1</rmsDec>
            <astCat>UNK</astCat>
            <mag>-4.2</mag>
            <rmsMag>1</rmsMag>
            <band>V</band>
        </optical>
        </ades>
        """
        optical4 = parse_optical_ades(s)

        re = r"<pos(?<i>\d)>(?<x>[\d\.]+)</pos.>"
        ss = replace(s, "<sys>ICRF_KM</sys>" => "<sys>ICRF_AU</sys>")
        for m in eachmatch(re, s)
           i, x = m
           y = parse(Float64, x) / au
           ss = replace(ss, m.match => "<pos$i>$y</pos$i>")
       end
       _optical4_ = parse_optical_ades(ss)

       @test maximum(norm, @. obsposECEF(optical4) - obsposECEF(_optical4_)) < eps()
       @test maximum(norm, @. obsposvelECI(optical4) - obsposvelECI(_optical4_)) < eps()
    end

    using NEOs: OpticalTracklet, OpticalMPC80, OpticalRWO, OpticalADES,
                opticalindices, reduce_tracklets, isunknown, closest_tracklet,
                daysbetween

    @testset "Tracklet" begin

        # Reduce tracklets
        @test_throws ArgumentError reduce_tracklets(OpticalMPC80{Float64}[])
        @test_throws ArgumentError reduce_tracklets(OpticalRWO{Float64}[])
        @test_throws ArgumentError reduce_tracklets(OpticalADES{Float64}[])

        trks1 = reduce_tracklets(optical1)
        trks2 = reduce_tracklets(optical2)
        trks3 = reduce_tracklets(optical3)

        @test isa(trks1, Vector{OpticalTracklet{Float64}})
        @test isa(trks2, Vector{OpticalTracklet{Float64}})
        @test isa(trks3, Vector{OpticalTracklet{Float64}})
        @test length(trks1) == length(trks2) > length(trks3)

        @test all(Base.Fix2(isa, String), sprint.(show, trks1))
        @test all(Base.Fix2(isa, String), sprint.(show, trks2))
        @test all(Base.Fix2(isa, String), sprint.(show, trks3))

        @test all(Base.Fix2(isa, String), sprint.(Ref(show), Ref(MIME("text/plain")), trks1))
        @test all(Base.Fix2(isa, String), sprint.(Ref(show), Ref(MIME("text/plain")), trks2))
        @test all(Base.Fix2(isa, String), sprint.(Ref(show), Ref(MIME("text/plain")), trks3))

        @test nobs(trks1) == length(optical1)
        @test nobs(trks2) == length(optical2)
        @test nobs(trks3) == length(optical3)

        @test opticalindices(trks1) == collect(eachindex(optical1))
        @test opticalindices(trks2) == collect(eachindex(optical2))
        @test opticalindices(trks3) == collect(eachindex(optical3))

        @test opticalindices(trks1, 1:10) == opticalindices(trks1[1:10])
        @test opticalindices(trks2, 1:10) == opticalindices(trks2[1:10])
        @test opticalindices(trks3, 1:10) == opticalindices(trks3[1:10])

        @test all(isempty, trackletid.(trks1))
        @test all(isempty, trackletid.(trks2))
        @test all(isempty, trackletid.(trks3))

        @test closest_tracklet(dtutc2days(optical1[1]), trks1) == 1
        @test closest_tracklet(dtutc2days(optical2[1]), trks2) == 1
        @test closest_tracklet(dtutc2days(optical3[1]), trks3) == 1

        @test maximum(datediff.(trks1, trks2)) == 0
        @test maximum(abs, ra.(trks1) - ra.(trks2)) == 0.0
        @test maximum(abs, dec.(trks1) - dec.(trks2)) == 0.0
        @test all(x -> isnan(x) || iszero(x), mag.(trks1) - mag.(trks2))
        @test observatory.(trks1) == observatory.(trks2)
        @test string.(trks1) == string.(trks2)

        mask13 = @. !isnothing($indexin(trks1, trks3))
        mask23 = @. !isnothing($indexin(trks2, trks3))
        mask31 = @. !isnothing($indexin(trks3, trks1))
        mask32 = @. !isnothing($indexin(trks3, trks2))
        @test mask13 == mask23
        @test mask31 == mask32
        @test count(mask13) == count(mask23) == count(mask31) == count(mask32)

        @test maximum(datediff.(trks1[mask13], trks3[mask31])) == 0
        @test maximum(datediff.(trks2[mask23], trks3[mask32])) == 0

        @test all(@.(
            daysbetween(trks1[mask13], trks3[mask31]) ==
            daysbetween(date(trks1[mask13]), trks3[mask31]) ==
            daysbetween(trks1[mask13], date(trks3[mask31])) ==
            daysbetween(date(trks1[mask13]), date(trks3[mask31]))
        ))
        @test all(@.(
            daysbetween(trks2[mask23], trks3[mask32]) ==
            daysbetween(date(trks2[mask23]), trks3[mask32]) ==
            daysbetween(trks2[mask23], date(trks3[mask32])) ==
            daysbetween(date(trks2[mask23]), date(trks3[mask32]))
        ))

        @test maximum(abs, ra.(trks1[mask13]) - ra.(trks3[mask31])) < 2.8e-7
        @test maximum(abs, ra.(trks2[mask23]) - ra.(trks3[mask32])) < 2.8e-7

        @test maximum(abs, dec.(trks1[mask13]) - dec.(trks3[mask31])) < 2.1e-7
        @test maximum(abs, dec.(trks2[mask23]) - dec.(trks3[mask32])) < 2.1e-7

        @test all(x -> isnan(x) || iszero(x), mag.(trks1[mask13]) - mag.(trks3[mask31]))
        @test all(x -> isnan(x) || iszero(x), mag.(trks2[mask23]) - mag.(trks3[mask32]))

        @test all(==(' '), band.(trks1))
        @test all(==(' '), band.(trks2))
        @test all(==(' '), band.(trks3))

        @test observatory.(trks1[mask13]) == observatory.(trks3[mask31])
        @test observatory.(trks2[mask23]) == observatory.(trks3[mask32])

        @test all(isunknown, catalogue.(trks1))
        @test all(isunknown, catalogue.(trks2))
        @test all(isunknown, catalogue.(trks3))

        @test all(x -> isnan(x[1]) && isnan(x[2]), rms.(trks1))
        @test all(x -> isnan(x[1]) && isnan(x[2]), rms.(trks2))
        @test all(x -> isnan(x[1]) && isnan(x[2]), rms.(trks3))

        @test all(x -> isnan(x[1]) && isnan(x[2]), debias.(trks1))
        @test all(x -> isnan(x[1]) && isnan(x[2]), debias.(trks2))
        @test all(x -> isnan(x[1]) && isnan(x[2]), debias.(trks3))

        @test all(isnan, corr.(trks1))
        @test all(isnan, corr.(trks2))
        @test all(isnan, corr.(trks3))

        @test string.(trks1[mask13]) == string.(trks3[mask31])
        @test string.(trks2[mask23]) == string.(trks3[mask32])

    end

    @testset "Apparition" begin

        using NEOs: RadarJPL, radartype, opticaltype, scalartype, radarindices

        apps1 = apparitions(optical1)
        apps2 = apparitions(optical2)
        apps3 = apparitions(optical3)

        @test all(Base.Fix2(isa, String), sprint.(show, apps1))
        @test all(Base.Fix2(isa, String), sprint.(show, apps2))
        @test all(Base.Fix2(isa, String), sprint.(show, apps3))

        @test all(Base.Fix2(isa, String), sprint.(Ref(show), Ref(MIME("text/plain")), apps1))
        @test all(Base.Fix2(isa, String), sprint.(Ref(show), Ref(MIME("text/plain")), apps2))
        @test all(Base.Fix2(isa, String), sprint.(Ref(show), Ref(MIME("text/plain")), apps3))

        @test all(==(RadarJPL{Float64}), radartype.(apps1))
        @test all(==(RadarJPL{Float64}), radartype.(apps2))
        @test all(==(RadarJPL{Float64}), radartype.(apps3))

        @test all(==(OpticalMPC80{Float64}), opticaltype.(apps1))
        @test all(==(OpticalRWO{Float64}), opticaltype.(apps2))
        @test all(==(OpticalADES{Float64}), opticaltype.(apps3))

        @test all(==(Float64), scalartype.(apps1))
        @test all(==(Float64), scalartype.(apps2))
        @test all(==(Float64), scalartype.(apps3))

        @test length(apps1) == length(apps2) == length(apps3)
        @test nradar(apps1) == nradar(apps2) == nradar(apps3) == 0
        @test noptical(apps1) == noptical(apps2) == noptical(apps3)
        @test nobs(apps1) == nobs(apps2) == nobs(apps3)

        @test radar(apps1) == radar(apps2) == radar(apps3) == RadarJPL[]
        @test optical(apps1) == optical1
        @test optical(apps2) == optical2
        @test optical(apps3) == optical3

        @test numberofdays(apps1) < numberofdays(optical1)
        @test numberofdays(apps2) < numberofdays(optical2)
        @test numberofdays(apps3) < numberofdays(optical3)

        @test radarindices(apps1) == radarindices(apps2) == radarindices(apps3) == Int[]
        @test opticalindices(apps1) == opticalindices(apps2) == opticalindices(apps3)

    end

    using NEOs: getid, update!

    @testset "AbstractWeightingScheme" begin

        # Weighting schemes
        w11 = UniformWeights(optical1)
        w12 = SourceWeights(optical1)
        w13 = Veres17(optical1)

        @test isa(string(w11), String) && isa(string(w12), String) && isa(string(w13), String)
        @test getid(w11) == "Uniform"
        @test getid(w12) == "Source"
        @test getid(w13) == "Veres et al. (2017)"

        Γ11, Γ12, Γ13 = covariances(w11), covariances(w12), covariances(w13)
        @test length(Γ11) == length(Γ12) == length(Γ13)
        @test all(==(radeccovariance(1.0, 1.0)), Γ11)
        @test all(==(radeccovariance(1.0, 1.0)), Γ12)
        @test all(Γ -> Γ[1, 1] == Γ[2, 2] > 0 && iszero(Γ[1, 2]) && iszero(Γ[2, 1]), Γ13)
        @test Γ11 == covariances.(Ref(w11), eachindex(Γ11)) == covariances(w11, eachindex(Γ11))
        @test Γ12 == covariances.(Ref(w12), eachindex(Γ12)) == covariances(w12, eachindex(Γ12))
        @test Γ13 == covariances.(Ref(w13), eachindex(Γ13)) == covariances(w13, eachindex(Γ13))

        ow11, ow12, ow13 = outliers(w11), outliers(w12), outliers(w13)
        @test length(ow11) == length(ow12) == length(ow13)
        @test all(!, ow11)
        @test all(!, ow12)
        @test all(!, ow13)
        @test ow11 == outliers.(Ref(w11), eachindex(ow11)) == outliers(w11, eachindex(ow11))
        @test ow12 == outliers.(Ref(w12), eachindex(ow12)) == outliers(w12, eachindex(ow12))
        @test ow13 == outliers.(Ref(w13), eachindex(ow13)) == outliers(w13, eachindex(ow13))

        setoutlier!.(Ref(w11), eachindex(ow11), Ref(true))
        setoutlier!.(Ref(w12), eachindex(ow12), Ref(true))
        setoutlier!.(Ref(w13), eachindex(ow13), Ref(true))
        @test all(ow11) && all(ow12) && all(ow13)

        setoutlier!(w11, eachindex(ow11), falses(length(ow11)))
        setoutlier!(w12, eachindex(ow12), falses(length(ow12)))
        setoutlier!(w13, eachindex(ow13), falses(length(ow13)))
        @test all(!, ow11) && all(!, ow12) && all(!, ow13)

        w21 = UniformWeights(optical2)
        w22 = SourceWeights(optical2)
        w23 = Veres17(optical2)

        @test isa(string(w21), String) && isa(string(w22), String) && isa(string(w23), String)
        @test getid(w21) == "Uniform"
        @test getid(w22) == "Source"
        @test getid(w23) == "Veres et al. (2017)"

        Γ21, Γ22, Γ23 = covariances(w21), covariances(w22), covariances(w23)
        @test length(Γ21) == length(Γ22) == length(Γ23)
        @test all(==(radeccovariance(1.0, 1.0)), Γ21)
        @test all( @. Γ22 == radeccovariance(getfield(optical2, :ra_rms),
            getfield(optical2, :dec_rms)) )
        @test all(Γ -> Γ[1, 1] == Γ[2, 2] > 0 && iszero(Γ[1, 2]) && iszero(Γ[2, 1]), Γ23)
        @test Γ21 == covariances.(Ref(w21), eachindex(Γ21)) == covariances(w21, eachindex(Γ21))
        @test Γ22 == covariances.(Ref(w22), eachindex(Γ22)) == covariances(w22, eachindex(Γ22))
        @test Γ23 == covariances.(Ref(w23), eachindex(Γ23)) == covariances(w23, eachindex(Γ23))

        ow21, ow22, ow23 = outliers(w21), outliers(w22), outliers(w23)
        @test length(ow21) == length(ow22) == length(ow23)
        @test all(!, ow21)
        @test count(ow22) == count(isoutlier, optical2)
        @test all(!, ow23)
        @test ow21 == outliers.(Ref(w21), eachindex(ow21)) == outliers(w21, eachindex(ow21))
        @test ow22 == outliers.(Ref(w22), eachindex(ow22)) == outliers(w22, eachindex(ow22))
        @test ow23 == outliers.(Ref(w23), eachindex(ow23)) == outliers(w23, eachindex(ow23))

        setoutlier!.(Ref(w21), eachindex(ow21), Ref(true))
        setoutlier!.(Ref(w22), eachindex(ow22), Ref(true))
        setoutlier!.(Ref(w23), eachindex(ow23), Ref(true))
        @test all(ow21) && all(ow22) && all(ow23)

        setoutlier!(w21, eachindex(ow21), falses(length(ow21)))
        setoutlier!(w22, eachindex(ow22), falses(length(ow22)))
        setoutlier!(w23, eachindex(ow23), falses(length(ow23)))
        @test all(!, ow21) && all(!, ow22) && all(!, ow23)

        w31 = UniformWeights(optical3)
        w32 = SourceWeights(optical3)
        w33 = Veres17(optical3)

        @test isa(string(w31), String) && isa(string(w32), String) && isa(string(w33), String)
        @test getid(w31) == "Uniform"
        @test getid(w32) == "Source"
        @test getid(w33) == "Veres et al. (2017)"

        Γ31, Γ32, Γ33 = covariances(w31), covariances(w32), covariances(w33)
        @test length(Γ31) == length(Γ32) == length(Γ33)
        @test all(==(radeccovariance(1.0, 1.0)), Γ31)
        @test all(map(Γ32, optical3) do Γ, x
            σx = rms(x)
            ρ = isnan(x.rmscorr) ? 0.0 : x.rmscorr
            Γ == radeccovariance(σx[1], σx[2], ρ)
        end)
        @test all(Γ -> Γ[1, 1] == Γ[2, 2] > 0 && iszero(Γ[1, 2]) && iszero(Γ[2, 1]), Γ33)
        @test Γ31 == covariances.(Ref(w31), eachindex(Γ31)) == covariances(w31, eachindex(Γ31))
        @test Γ32 == covariances.(Ref(w32), eachindex(Γ32)) == covariances(w32, eachindex(Γ32))
        @test Γ33 == covariances.(Ref(w33), eachindex(Γ33)) == covariances(w33, eachindex(Γ33))

        ow31, ow32, ow33 = outliers(w31), outliers(w32), outliers(w33)
        @test length(ow31) == length(ow32) == length(ow33)
        @test all(!, ow31)
        @test all(!, ow32)
        @test all(!, ow33)
        @test ow31 == outliers.(Ref(w31), eachindex(ow31)) == outliers(w31, eachindex(ow31))
        @test ow32 == outliers.(Ref(w32), eachindex(ow32)) == outliers(w32, eachindex(ow32))
        @test ow33 == outliers.(Ref(w33), eachindex(ow33)) == outliers(w33, eachindex(ow33))

        setoutlier!.(Ref(w31), eachindex(ow31), Ref(true))
        setoutlier!.(Ref(w32), eachindex(ow32), Ref(true))
        setoutlier!.(Ref(w33), eachindex(ow33), Ref(true))
        @test all(ow31) && all(ow32) && all(ow33)

        setoutlier!(w31, eachindex(ow31), falses(length(ow31)))
        setoutlier!(w32, eachindex(ow32), falses(length(ow32)))
        setoutlier!(w33, eachindex(ow33), falses(length(ow33)))
        @test all(!, ow31) && all(!, ow32) && all(!, ow33)

        @test Γ13 == Γ23 == Γ33

        update!(w11, optical1[1:1])
        update!(w12, optical1[1:1])
        update!(w13, optical1[1:1])

        @test covariances(w11) == [radeccovariance(1.0, 1.0)] && outliers(w11) == [false]
        @test covariances(w12) == [radeccovariance(1.0, 1.0)] && outliers(w12) == [false]
        @test covariances(w13) == [radeccovariance(1.0, 1.0)] && outliers(w13) == [false]

        update!(w21, optical2[1:1])
        update!(w22, optical2[1:1])
        update!(w23, optical2[1:1])

        @test covariances(w21) == [radeccovariance(1.0, 1.0)] && outliers(w21) == [false]
        @test covariances(w22) == [radeccovariance(1.0, 1.0)] && outliers(w22) == [false]
        @test covariances(w23) == [radeccovariance(1.0, 1.0)] && outliers(w23) == [false]

        update!(w31, optical3[1:1])
        update!(w32, optical3[1:1])
        update!(w33, optical3[1:1])

        @test covariances(w31) == [radeccovariance(1.0, 1.0)] && outliers(w31) == [false]
        @test covariances(w32) == [radeccovariance(1.0, 1.0)] && outliers(w32) == [false]
        @test covariances(w33) == [radeccovariance(1.0, 1.0)] && outliers(w33) == [false]

        # Custom weights
        w14 = Veres17(optical1)
        c14 = CustomWeights(w14)
        @test isa(string(c14), String)
        @test getid(c14) == "Custom"
        @test covariances(c14) == covariances(w14) && outliers(c14) == outliers(w14)
        @test covariances(c14) !== covariances(w14) && outliers(c14) !== outliers(w14)
        c14s = CustomWeights(w14, 1:2)
        @test nobs(c14s) == 2 && covariances(c14s) == covariances(w14, 1:2)
        cbad = CustomWeights([radeccovariance(1.0, 1.0, 1.0), radeccovariance(1.0, 1.0),
                              radeccovariance(1.0, 1.0, NaN)])
        @test outliers(cbad) == [true, false, true]
        @test outliers(CustomWeights(cbad)) == outliers(cbad)
        @test isnothing(update!(c14, optical1))
        @test_throws ArgumentError update!(cbad, optical1[1:1])
        od14 = ODProblem(newtonian!, optical1; weights = c14, debias = ZeroDebiasing)
        @test od14.weights === c14
        @test covariances(od14) == covariances(w14)
        @test_throws AssertionError ODProblem(newtonian!, optical1[1:2]; weights = c14,
            debias = ZeroDebiasing)

    end

    @testset "Observer position uncertainty" begin

        using NEOs: observercovariance, isposdefcov

        # Synthetic satellite observations: ICRF_KM with an anisotropic position
        # covariance, ICRF_AU with an isotropic one and ICRF_KM without covariance
        sat = read_optical_ades(joinpath(TEST_DATA, "satellite.xml"))
        @test length(sat) == 3
        @test all(x -> observatory(x).frame in ("ICRF_KM", "ICRF_AU"), sat)

        # Observer's position covariance matrices [km²]
        Σ1, Σ2, Σ3 = poscov.(sat)
        @test Σ1 == [1.0 0.5 0.0; 0.5 2.0 0.0; 0.0 0.0 3.0]
        @test Σ2 ≈ (NEOs.au^2 * 1E-16) * I(3)
        @test isnothing(Σ3)
        @test all(isnothing, poscov.(optical1))

        # Observer to asteroid vector [km] (|ρ| = 5E4 km)
        ρ_vec = [3E4, -4E4, 0.0]
        k = rad2arcsec(1.0)^2 / 5E4^2

        # No covariance: no contribution
        @test iszero(observercovariance(sat[3], ρ_vec))
        @test iszero(observercovariance(optical1[1], ρ_vec))

        # Isotropic covariance: Γ = (σ / ρ)² I
        Γ2 = observercovariance(sat[2], ρ_vec)
        @test Γ2 ≈ k * (NEOs.au^2 * 1E-16) * I(2)

        # Anisotropic covariance: Γ = P Σ Pᵀ / ρ²
        Γ1 = observercovariance(sat[1], ρ_vec)
        α, δ = measure(sat[1])
        P = [-sin(α) cos(α) 0.0; -sin(δ)*cos(α) -sin(δ)*sin(α) cos(δ)]
        @test Γ1 ≈ k * P * Σ1 * P'
        @test Γ1[1, 2] == Γ1[2, 1]
        @test isposdefcov(Γ1)
        # The contribution scales as 1 / ρ²
        @test observercovariance(sat[1], 2 * ρ_vec) ≈ Γ1 / 4

        # Total weight matrix
        Γε = radeccovariance(1.0, 1.0)
        W = weightmatrix(Γε + Γ1)
        @test W * (Γε + Γ1) ≈ I(2)
        @test all(diag(W) .< diag(weightmatrix(Γε)))

    end

    @testset "AbstractDebiasingScheme" begin

        # Debiasing schemes
        d11 = ZeroDebiasing(optical1)
        d12 = SourceDebiasing(optical1)
        d13 = Farnocchia15(optical1)
        d14 = Eggl20(optical1)

        @test isa(string(d11), String) && isa(string(d12), String) && isa(string(d13), String) &&
              isa(string(d14), String)
        @test getid(d11) == "Zero"
        @test getid(d12) == "Source"
        @test getid(d13) == "Farnocchia et al. (2015)"
        @test getid(d14) == "Eggl et al. (2020)"

        @test length(debias(d11)) == length(debias(d12)) == length(debias(d13)) == length(debias(d14))
        @test all(==((0.0, 0.0)), debias(d11))
        @test all(==((0.0, 0.0)), debias(d12))

        d21 = ZeroDebiasing(optical2)
        d22 = SourceDebiasing(optical2)
        d23 = Farnocchia15(optical2)
        d24 = Eggl20(optical2)

        @test isa(string(d21), String) && isa(string(d22), String) && isa(string(d23), String) &&
              isa(string(d24), String)
        @test getid(d21) == "Zero"
        @test getid(d22) == "Source"
        @test getid(d23) == "Farnocchia et al. (2015)"
        @test getid(d24) == "Eggl et al. (2020)"

        @test length(debias(d21)) == length(debias(d22)) == length(debias(d23)) == length(debias(d24))
        @test all(==((0.0, 0.0)), debias(d21))
        @test all(@. $debias(d22) == tuple(getfield(optical2, :ra_bias), getfield(optical2, :dec_bias)))

        d31 = ZeroDebiasing(optical3)
        d32 = SourceDebiasing(optical3)
        d33 = Farnocchia15(optical3)
        d34 = Eggl20(optical3)

        @test isa(string(d31), String) && isa(string(d32), String) && isa(string(d33), String) &&
              isa(string(d34), String)
        @test getid(d31) == "Zero"
        @test getid(d32) == "Source"
        @test getid(d33) == "Farnocchia et al. (2015)"
        @test getid(d34) == "Eggl et al. (2020)"

        @test length(debias(d31)) == length(debias(d32)) == length(debias(d33)) == length(debias(d34))
        @test all(==((0.0, 0.0)), debias(d31))
        @test all(==((0.0, 0.0)), debias(d32))

        @test debias(d13) == debias(d23)
        @test maximum(map((x, y) -> hypot(x[1] - y[1], x[2] - y[2]), debias(d13), debias(d33))) < 0.79
        @test maximum(map((x, y) -> hypot(x[1] - y[1], x[2] - y[2]), debias(d23), debias(d33))) < 0.79
        @test debias(d14) == debias(d24)
        @test maximum(map((x, y) -> hypot(x[1] - y[1], x[2] - y[2]), debias(d14), debias(d34))) < 0.95
        @test maximum(map((x, y) -> hypot(x[1] - y[1], x[2] - y[2]), debias(d24), debias(d34))) < 0.95

        @testset "Debiasing table reuse" begin
            # Schemes share the table, but never the observation corrections.
            @test d13.table === d23.table === d33.table
            @test d14.table === d24.table === d34.table
            @test d13.table !== d14.table
            @test debias(d13) !== debias(d23)
            @test debias(d14) !== debias(d24)

            for (id, scheme) in (("2014", d13), ("2018", d14))
                data = @inferred NEOs.select_debiasing_table(Float64, id)
                reference = NEOs.read_debiasing_table(Float64, id)
                @test data == reference
                @test last(data) === scheme.table
                @test debias(scheme) == NEOs.debiasing(optical1, reference...)
            end

            # Simultaneous requests for another number type share a single table.
            tasks = [Threads.@spawn NEOs.select_debiasing_table(Float32, "2014")
                     for _ in 1:4]
            tables = fetch.(tasks)
            table32 = last(first(tables))
            @test all(x -> last(x) === table32, tables)
            @test table32 isa Matrix{Float32}
            @test table32 == Float32.(d13.table)
            @test last(@inferred NEOs.select_debiasing_table(Float32, "2014")) === table32

            # A failed load must leave subsequent requests usable.
            @test_throws ArgumentError NEOs.select_debiasing_table(Float64, "unknown")
            @test_throws ArgumentError NEOs.select_debiasing_table(Float64, "unknown")
            @test last(NEOs.select_debiasing_table(Float64)) === d14.table
        end

        update!(d11, optical1[1:1])
        update!(d12, optical1[1:1])
        update!(d13, optical1[1:1])
        update!(d14, optical1[1:1])

        @test d11.debias == [(0.0, 0.0)]
        @test d12.debias == [(0.0, 0.0)]
        @test d13.debias != [(0.0, 0.0)]
        @test d14.debias != [(0.0, 0.0)]

        update!(d21, optical2[1:1])
        update!(d22, optical2[1:1])
        update!(d23, optical2[1:1])
        update!(d24, optical2[1:1])

        @test d21.debias == [(0.0, 0.0)]
        @test d22.debias != [(0.0, 0.0)]
        @test d23.debias != [(0.0, 0.0)]
        @test d24.debias != [(0.0, 0.0)]

        update!(d31, optical3[1:1])
        update!(d32, optical3[1:1])
        update!(d33, optical3[1:1])
        update!(d34, optical3[1:1])

        @test d31.debias == [(0.0, 0.0)]
        @test d32.debias == [(0.0, 0.0)]
        @test d33.debias != [(0.0, 0.0)]
        @test d34.debias != [(0.0, 0.0)]

        @test d13.table === d23.table === d33.table
        @test d14.table === d24.table === d34.table
        @test length(debias(d13)) == length(debias(d14)) == 1

    end

end
