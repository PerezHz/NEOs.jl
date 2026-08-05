using NEOs
using Dates
using PlanetaryEphemeris
using LinearAlgebra
using DataFrames
using Query
using Test

const TEST_DATA = joinpath(pkgdir(NEOs), "test", "data")

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
        filter!(x -> Date(2000, 1) < date(x) < Date(2025, 6), optical1)
        @test isa(optical1, Vector{OpticalMPC80{Float64}})
        @test issorted(optical1)
        @test allunique(optical1)

        # Read/write OpticalMPC80 file
        filename = joinpath(TEST_DATA, "433.txt")
        optical2 = read_optical_mpc80(filename)
        filter!(x -> Date(2000, 1) < date(x) < Date(2025, 6), optical2)
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
        filter!(x -> Date(2000, 1) < date(x) < Date(2025, 6), optical1)
        @test isa(optical1, Vector{OpticalRWO{Float64}})
        @test issorted(optical1)
        @test allunique(optical1)

        optical2 = fetch_optical_rwo("433", NEODyS2)
        filter!(x -> Date(2000, 1) < date(x) < Date(2025, 6), optical2)
        @test isa(optical2, Vector{OpticalRWO{Float64}})
        @test issorted(optical2)
        @test allunique(optical2)

        # Read/write OpticalRWO file
        filename = joinpath(TEST_DATA, "433.rwo")
        optical3 = read_optical_rwo(filename)
        filter!(x -> Date(2000, 1) < date(x) < Date(2025, 6), optical3)
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
        filter!(x -> Date(2000, 1) < date(x) < Date(2025, 6), optical1)
        @test isa(optical1, Vector{OpticalADES{Float64}})
        @test issorted(optical1)
        @test allunique(optical1)

        # Read/write OpticalRWO file
        filename = joinpath(TEST_DATA, "433.xml")
        optical2 = read_optical_ades(filename)
        filter!(x -> Date(2000, 1) < date(x) < Date(2025, 6), optical2)
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
    filename = joinpath(TEST_DATA, "433.txt")
    optical1 = read_optical_mpc80(filename)
    filter!(x -> Date(2000) < date(x) < Date(2025), optical1)
    filename = joinpath(TEST_DATA, "433.rwo")
    optical2 = read_optical_rwo(filename)
    filter!(x -> Date(2000) < date(x) < Date(2025), optical2)
    filename = joinpath(TEST_DATA, "433.xml")
    optical3 = read_optical_ades(filename)
    filter!(x -> Date(2000) < date(x) < Date(2025), optical3)

    @testset "read_optical_astrometry" begin
        mpc80_file = joinpath(TEST_DATA, "433.txt")
        ades_file = joinpath(TEST_DATA, "433.xml")

        @test read_optical_astrometry(mpc80_file) == read_optical_mpc80(mpc80_file)
        @test read_optical_astrometry(mpc80_file; format = :obs80) ==
            read_optical_mpc80(mpc80_file)
        @test read_optical_astrometry(ades_file) == read_optical_ades(ades_file)
        @test read_optical_astrometry(ades_file; format = "xml") ==
            read_optical_ades(ades_file)
        @test_throws ArgumentError read_optical_astrometry(mpc80_file; format = :unknown)
    end

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
        @test count(mask1) == count(mask2) == count(mask3) == 6_765
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
                indices, reduce_tracklets, isunknown, closest_tracklet

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

        @test nobs(trks1) == length(optical1)
        @test nobs(trks2) == length(optical2)
        @test nobs(trks3) == length(optical3)

        @test indices(trks1) == collect(eachindex(optical1))
        @test indices(trks2) == collect(eachindex(optical2))
        @test indices(trks3) == collect(eachindex(optical3))

        @test indices(trks1, 1:10) == indices(trks1[1:10])
        @test indices(trks2, 1:10) == indices(trks2[1:10])
        @test indices(trks3, 1:10) == indices(trks3[1:10])

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

        apps1 = apparitions(optical1)
        apps2 = apparitions(optical2)
        apps3 = apparitions(optical3)

        @test all(Base.Fix2(isa, String), string.(apps1))
        @test all(Base.Fix2(isa, String), string.(apps2))
        @test all(Base.Fix2(isa, String), string.(apps3))

        @test length(apps1) == length(apps2) == length(apps3)
        @test indices(apps1) == eachindex(optical1)
        @test indices(apps2) == eachindex(optical2)
        @test indices(apps3) == eachindex(optical3)

        @test optical(apps1) == optical1
        @test optical(apps2) == optical2
        @test optical(apps3) == optical3

        @test numberofdays(apps1) < numberofdays(optical1)
        @test numberofdays(apps2) < numberofdays(optical2)
        @test numberofdays(apps3) < numberofdays(optical3)

        @test noptical(apps1) == length(optical1)
        @test noptical(apps2) == length(optical2)
        @test noptical(apps3) == length(optical3)

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

        ww11, ww12, ww13 = weights(w11), weights(w12), weights(w13)
        @test length(ww11) == length(ww12) == length(ww13)
        @test all(==((1.0, 1.0)), ww11)
        @test all(==((1.0, 1.0)), ww12)
        @test ww11 == weights.(Ref(w11), eachindex(ww11)) == weights(w11, eachindex(ww11))
        @test ww12 == weights.(Ref(w12), eachindex(ww12)) == weights(w12, eachindex(ww12))
        @test ww13 == weights.(Ref(w13), eachindex(ww13)) == weights(w13, eachindex(ww13))

        cw11, cw12, cw13 = corrs(w11), corrs(w12), corrs(w13)
        @test length(cw11) == length(cw12) == length(cw13)
        @test all(iszero, cw11)
        @test all(iszero, cw12)
        @test all(iszero, cw13)
        @test cw11 == corrs.(Ref(w11), eachindex(cw11)) == corrs(w11, eachindex(cw11))
        @test cw12 == corrs.(Ref(w12), eachindex(cw12)) == corrs(w12, eachindex(cw12))
        @test cw13 == corrs.(Ref(w13), eachindex(cw13)) == corrs(w13, eachindex(cw13))

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

        ww21, ww22, ww23 = weights(w21), weights(w22), weights(w23)
        @test length(ww21) == length(ww22) == length(ww23)
        @test all(==((1.0, 1.0)), ww21)
        @test all( @. ww22 == tuple(1 / getfield(optical2, :ra_rms),
            1 / getfield(optical2, :dec_rms)) )
        @test ww21 == weights.(Ref(w21), eachindex(ww21)) == weights(w21, eachindex(ww21))
        @test ww22 == weights.(Ref(w22), eachindex(ww22)) == weights(w22, eachindex(ww22))
        @test ww23 == weights.(Ref(w23), eachindex(ww23)) == weights(w23, eachindex(ww23))

        cw21, cw22, cw23 = corrs(w21), corrs(w22), corrs(w23)
        @test length(cw21) == length(cw22) == length(cw23)
        @test all(iszero, cw21)
        @test all(iszero, cw22)
        @test all(iszero, cw23)
        @test cw21 == corrs.(Ref(w21), eachindex(cw21)) == corrs(w21, eachindex(cw21))
        @test cw22 == corrs.(Ref(w22), eachindex(cw22)) == corrs(w22, eachindex(cw22))
        @test cw23 == corrs.(Ref(w23), eachindex(cw23)) == corrs(w23, eachindex(cw23))

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

        ww31, ww32, ww33 = weights(w31), weights(w32), weights(w33)
        @test length(ww31) == length(ww32) == length(ww33)
        @test all(==((1.0, 1.0)), ww31)
        @test all(map(ww32, optical3) do w, x
            σx = rms(x)
            w == (1 / σx[1], 1 / σx[2])
        end)
        @test ww31 == weights.(Ref(w31), eachindex(ww31)) == weights(w31, eachindex(ww31))
        @test ww32 == weights.(Ref(w32), eachindex(ww32)) == weights(w32, eachindex(ww32))
        @test ww33 == weights.(Ref(w33), eachindex(ww33)) == weights(w33, eachindex(ww33))

        cw31, cw32, cw33 = corrs(w31), corrs(w32), corrs(w33)
        @test length(cw31) == length(cw32) == length(cw33)
        @test all(iszero, cw31)
        @test all(map(cw32, optical3) do ρ, x
            ρ == (isnan(x.rmscorr) ? 0.0 : x.rmscorr)
        end)
        @test all(iszero, cw33)
        @test cw31 == corrs.(Ref(w31), eachindex(cw31)) == corrs(w31, eachindex(cw31))
        @test cw32 == corrs.(Ref(w32), eachindex(cw32)) == corrs(w32, eachindex(cw32))
        @test cw33 == corrs.(Ref(w33), eachindex(cw33)) == corrs(w33, eachindex(cw33))

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

        @test ww13 == ww23 == ww33
        @test cw13 == cw23 == cw33

        update!(w11, optical1[1:1])
        update!(w12, optical1[1:1])
        update!(w13, optical1[1:1])

        @test w11.weights == [(1.0, 1.0)] && w11.corrs == [0.0] && w11.outliers == [false]
        @test w12.weights == [(1.0, 1.0)] && w12.corrs == [0.0] && w12.outliers == [false]
        @test w13.weights == [(1.0, 1.0)] && w13.corrs == [0.0] && w13.outliers == [false]

        update!(w21, optical2[1:1])
        update!(w22, optical2[1:1])
        update!(w23, optical2[1:1])

        @test w21.weights == [(1.0, 1.0)] && w21.corrs == [0.0] && w21.outliers == [false]
        @test w22.weights == [(1.0, 1.0)] && w22.corrs == [0.0] && w22.outliers == [false]
        @test w23.weights == [(1.0, 1.0)] && w23.corrs == [0.0] && w23.outliers == [false]

        update!(w31, optical3[1:1])
        update!(w32, optical3[1:1])
        update!(w33, optical3[1:1])

        @test w31.weights == [(1.0, 1.0)] && w31.corrs == [0.0] && w31.outliers == [false]
        @test w32.weights == [(1.0, 1.0)] && w32.corrs == [0.0] && w32.outliers == [false]
        @test w33.weights == [(1.0, 1.0)] && w33.corrs == [0.0] && w33.outliers == [false]

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

    end

end
