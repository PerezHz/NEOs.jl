using NEOs
using Dates
using Test

@testset "AbstractRadarAstrometry" begin

    @testset "RadarJPL" begin

        using NEOs: RadarJPL, parse_radar_jpl, isdelay, isdoppler, ismonostatic

        # Parse RadarJPL
        apophis_s = "{\"signature\":{\"source\":\"NASA/JPL Small-Body Radar Astrometry API\",\
            \"version\":\"1.1\"},\"count\":\"50\",\"fields\":[\"des\",\"epoch\",\"value\",\
            \"sigma\",\"units\",\"freq\",\"rcvr\",\"xmit\",\"bp\"],\"data\":[[\"99942\",\
            \"2005-01-27 23:31:00\",\"-100849.1434\",\"0.250\",\"Hz\",\"2380\",\"-1\",\"-1\",\
            \"C\"]]}\n"
        apophis_p = parse_radar_jpl(apophis_s)
        @test isa(apophis_p, Vector{RadarJPL{Float64}})
        @test isone(length(apophis_p))
        apophis = first(apophis_p)
        @test apophis.des == "99942"
        @test apophis.epoch == DateTime("2005-01-27 23:31:00", "yyyy-mm-dd HH:MM:SS")
        @test apophis.value == -100849.1434
        @test apophis.sigma == 0.250
        @test apophis.units == "Hz"
        @test apophis.freq == 2380.0
        @test apophis.rcvr == search_observatory_code("251")
        @test apophis.xmit == search_observatory_code("251")
        @test apophis.bp == "C"

        @test !isdelay(apophis)
        @test isdoppler(apophis)
        @test ismonostatic(apophis)

        # RadarJPL equality
        @test apophis == apophis

        # Fetch RadarJPL
        radar1 = fetch_radar_jpl("spk" => "20000433", JPL)
        @test isa(radar1, Vector{RadarJPL{Float64}})
        @test issorted(radar1)
        @test allunique(radar1)

        # Read/write OpticalMPC80 file
        filename = joinpath(pkgdir(NEOs), "test", "data", "433.json")
        radar2 = read_radar_jpl(filename)
        @test isa(radar2, Vector{RadarJPL{Float64}})
        @test issorted(radar2)
        @test allunique(radar2)

        filename = joinpath(pkgdir(NEOs), "test", "data", "433_.json")
        write_radar_jpl(radar2, filename)
        radar3 = read_radar_jpl(filename)
        rm(filename)
        @test isa(radar3, Vector{RadarJPL{Float64}})
        @test issorted(radar3)
        @test allunique(radar3)

        @test radar1 == radar2 == radar3

    end

    @testset "RadarRWO" begin

        using NEOs: RadarRWO, parse_radar_rwo, isdelay, isdoppler, ismonostatic

        # Parse RadarRWO
        apophis_s =
        """
        version =   3
        errmod  = 'vfcc17'
        RMSast  =   3.02422E-01
        RMSmag  =   4.31723E-01
        END_OF_HEADER
        ! Object   Obser ====== Date =======  ============ Radar range/range rate (km or km/d) =============  Station  ====  Residual
        ! Design   K T N YYYY MM DD hh:mm:ss         Measure  Accuracy       rms F        Bias       Resid    TRX  RCX        Chi   S
         99942     V c   2005 01 27 23:31:00    548781.78980   1.36040   1.36040 F     0.00000     0.22763    -1   -1         0.17  1
        """
        apophis_p = parse_radar_rwo(apophis_s)
        @test isa(apophis_p, Vector{RadarRWO{Float64}})
        @test isone(length(apophis_p))
        apophis = first(apophis_p)
        @test apophis.design == "99942"
        @test apophis.K == 'V'
        @test apophis.T == 'c'
        @test apophis.N == ' '
        @test apophis.date == DateTime("2005-01-27T23:31:00")
        @test apophis.measure == 548781.7898
        @test apophis.accuracy == 1.3604
        @test apophis.rms == 1.3604
        @test !apophis.flag
        @test apophis.bias == 0.0
        @test apophis.resid == 0.22763
        @test apophis.trx == search_observatory_code("251")
        @test apophis.rcx == search_observatory_code("251")
        @test apophis.chi == 0.17
        @test apophis.S == true

        @test !isdelay(apophis)
        @test isdoppler(apophis)
        @test ismonostatic(apophis)

        # RadarRWO equality
        @test apophis == apophis

        # Fetch RadarRWO
        radar1 = fetch_radar_rwo("433", NEOCC)
        @test isa(radar1, Vector{RadarRWO{Float64}})
        @test issorted(radar1)
        @test allunique(radar1)

        # Read/write OpticalMPC80 file
        filename = joinpath(pkgdir(NEOs), "test", "data", "433.rwo")
        radar2 = read_radar_rwo(filename)
        @test isa(radar2, Vector{RadarRWO{Float64}})
        @test issorted(radar2)
        @test allunique(radar2)

        filename = joinpath(pkgdir(NEOs), "test", "data", "433_.rwo")
        write_radar_rwo(radar2, filename)
        radar3 = read_radar_rwo(filename)
        rm(filename)
        @test isa(radar3, Vector{RadarRWO{Float64}})
        @test issorted(radar3)
        @test allunique(radar3)

        @test radar1 == radar2 == radar3

    end

end