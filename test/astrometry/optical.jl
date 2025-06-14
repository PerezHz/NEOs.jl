using NEOs
using Dates
using Test

@testset "AbstractOpticalAstrometry" begin

    @testset "OpticalMPC80" begin

        using NEOs: OpticalMPC80, parse_optical_mpc80

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
        @test apophis.date == DateTime("2004-03-15T02:35:21.696")
        @test apophis.ra == 1.0739650841580173
        @test apophis.dec == 0.2952738332250385
        @test apophis.info1 == ""
        @test isnan(apophis.mag)
        @test apophis.band == ' '
        @test apophis.catalogue == search_catalogue_code('o')
        @test apophis.info2 == "m6394"
        @test apophis.observatory == search_observatory_code("691")
        @test apophis.source == apophis_s

        # OpticalMPC80 equality
        @test apophis == apophis

        # Fetch OpticalMPC80
        optical1 = fetch_optical_mpc80("433", MPC)
        @test isa(optical1, Vector{OpticalMPC80{Float64}})
        @test issorted(optical1)
        @test allunique(optical1)

        # Read/write OpticalMPC80 file
        filename = joinpath(pkgdir(NEOs), "test", "data", "433.txt")
        optical2 = read_optical_mpc80(filename)
        @test isa(optical2, Vector{OpticalMPC80{Float64}})
        @test issorted(optical2)
        @test allunique(optical2)

        filename = joinpath(pkgdir(NEOs), "test", "data", "433_.txt")
        write_optical_mpc80(optical2, filename)
        optical3 = read_optical_mpc80(filename)
        rm(filename)
        @test isa(optical3, Vector{OpticalMPC80{Float64}})
        @test issorted(optical3)
        @test allunique(optical3)

        @test optical1 == optical2 == optical3

    end

    @testset "NEOCPObject" begin

        using NEOs: NEOCPObject, parse_neocp_objects

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
        @test X85177.date == DateTime("2025-05-24T14:24:00")
        @test X85177.ra == 4.10286764571071
        @test X85177.dec == -0.1217646405946364
        @test X85177.V == 21.5
        @test X85177.updated == "Updated June 13.68 UT"
        @test X85177.nobs == 3
        @test X85177.arc == 0.03
        @test X85177.H == 19.1
        @test X85177.notseen == 20.282
        @test X85177.source == X85177_s[1:end-1]

        # NEOCPObject equality
        @test X85177 == X85177

        # Fetch NEOCPObject
        objects1 = fetch_neocp_objects()
        @test isa(objects1, Vector{NEOCPObject{Float64}})
        @test issorted(objects1)
        @test allunique(objects1)

        # Read/write NEOCPObject file
        filename = joinpath(pkgdir(NEOs), "test", "data", "neocp.txt")
        objects2 = read_neocp_objects(filename)
        @test isa(objects2, Vector{NEOCPObject{Float64}})
        @test issorted(objects2)
        @test allunique(objects2)

        filename = joinpath(pkgdir(NEOs), "test", "data", "neocp_.txt")
        write_neocp_objects(objects2, filename)
        objects3 = read_neocp_objects(filename)
        rm(filename)
        @test isa(objects3, Vector{NEOCPObject{Float64}})
        @test issorted(objects3)
        @test allunique(objects3)

        @test objects2 == objects3

    end

    @testset "OpticalRWO" begin

        using NEOs: OpticalRWO, parse_optical_rwo

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
        @test apophis.date == DateTime("2004-03-15T02:35:21.696")
        @test apophis.date_accuracy == 1e-5
        @test apophis.ra == 1.0739650841580173
        @test apophis.ra_accuracy == 0.1435
        @test apophis.ra_rms == 0.612
        @test apophis.ra_flag == false
        @test apophis.ra_bias == -0.247
        @test apophis.ra_resid == 0.251
        @test apophis.dec == 0.2952738332250385
        @test apophis.dec_accuracy == 0.1
        @test apophis.dec_rms == 0.612
        @test apophis.dec_flag == false
        @test apophis.dec_bias == 0.14
        @test apophis.dec_resid == -0.07
        @test isnan(apophis.mag)
        @test apophis.mag_band == ' '
        @test isnan(apophis.mag_rms)
        @test isnan(apophis.mag_resid)
        @test apophis.catalogue == search_catalogue_code('o')
        @test apophis.observatory == search_observatory_code("691")
        @test apophis.chi == 0.43
        @test apophis.sel_A == true
        @test apophis.sel_M == false
        @test apophis.source == apophis_s[492:end]
        @test apophis.header == apophis_s[1:81]

        # OpticalRWO equality
        @test apophis == apophis

        # Fetch OpticalRWO
        optical1 = fetch_optical_rwo("433", NEOCC)
        @test isa(optical1, Vector{OpticalRWO{Float64}})
        @test issorted(optical1)
        @test allunique(optical1)

        # Read/write OpticalRWO file
        filename = joinpath(pkgdir(NEOs), "test", "data", "433.rwo")
        optical2 = read_optical_rwo(filename)
        @test isa(optical2, Vector{OpticalRWO{Float64}})
        @test issorted(optical2)
        @test allunique(optical2)

        filename = joinpath(pkgdir(NEOs), "test", "data", "433_.rwo")
        write_optical_rwo(optical2, filename)
        optical3 = read_optical_rwo(filename)
        rm(filename)
        @test isa(optical3, Vector{OpticalRWO{Float64}})
        @test issorted(optical3)
        @test allunique(optical3)

        @test optical1 == optical2 == optical3

    end

    @testset "OpticalADES" begin
        using NEOs: OpticalADES, parse_optical_ades, unknowncat

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
        @test apophis.obsid == "JqcHxe000000DaKP010000001"
        @test apophis.trkid == "000002w-NJ"
        @test apophis.mode == "CCD"
        @test apophis.stn == search_observatory_code("691")
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
        @test apophis.obstime == DateTime("2004-03-15T02:35:21.696")
        @test isnan(apophis.rmstime)
        @test apophis.ra == 1.073965142335659
        @test apophis.dec == 0.2952737556548495
        @test isnan(apophis.rastar)
        @test isnan(apophis.decstar)
        @test isnan(apophis.deltara)
        @test isnan(apophis.deltadec)
        @test isnan(apophis.rmsra)
        @test isnan(apophis.rmsdec)
        @test isnan(apophis.rmscorr)
        @test apophis.astcat == search_catalogue_code('U')
        @test isnan(apophis.mag)
        @test isnan(apophis.rmsmag)
        @test apophis.band == ""
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

        # OpticalRWO equality
        @test apophis == apophis

        # Fetch OpticalADES
        optical1 = fetch_optical_ades("433", MPC)
        @test isa(optical1, Vector{OpticalADES{Float64}})
        @test issorted(optical1)
        @test allunique(optical1)

        # Read/write OpticalRWO file
        filename = joinpath(pkgdir(NEOs), "test", "data", "433.xml")
        optical2 = read_optical_ades(filename)
        @test isa(optical2, Vector{OpticalADES{Float64}})
        @test issorted(optical2)
        @test allunique(optical2)

        filename = joinpath(pkgdir(NEOs), "test", "data", "433_.xml")
        write_optical_ades(optical2, filename)
        optical3 = read_optical_ades(filename)
        rm(filename)
        @test isa(optical3, Vector{OpticalADES{Float64}})
        @test issorted(optical3)
        @test allunique(optical3)

        @test optical1 == optical2 == optical3

    end

end