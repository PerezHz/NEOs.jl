using NEOs
using Dates
using DataFrames
using Query
using Test

using NEOs: read_radec_mpc, RadecMPC, read_radar_jpl, RadarJPL

@testset "DataFrame(::AbstractAstrometry)" begin

    # RadecMPC
    radec_2023DW = read_radec_mpc(joinpath(pkgdir(NEOs), "test", "data", "RADEC_2023_DW.dat"))

    df_2023DW = DataFrame(radec_2023DW)

    @test nrow(df_2023DW) == length(radec_2023DW)
    @test all(names(df_2023DW) .== String.(fieldnames(RadecMPC{Float64})))
    @test all(eltype.(eachcol(df_2023DW)) .== fieldtypes(RadecMPC{Float64}))

    recovered_2023DW = Vector{RadecMPC{Float64}}(df_2023DW)

    @test recovered_2023DW == radec_2023DW

    # Query RadecMPC
    I52_2023_02_28 = radec_2023DW |> @filter(Date(date(_)) == Date(2023, 2, 28) && observatory(_) == search_obs_code("I52")) |> DataFrame

    @test nrow(I52_2023_02_28) == 4
    @test all(Date.(I52_2023_02_28.date) .== Date(2023, 2, 28))
    @test allequal(I52_2023_02_28.observatory) && I52_2023_02_28.observatory[1] == search_obs_code("I52")

    # RadarJPL
    radar_Apophis = read_radar_jpl(joinpath("data", "99942_RADAR_2005_2013.dat"))

    df_Apophis = DataFrame(radar_Apophis)

    @test nrow(df_Apophis) == length(radar_Apophis)
    @test all(names(df_Apophis) .== String.(fieldnames(RadarJPL{Float64})))
    @test all(eltype.(eachcol(df_Apophis)) .== fieldtypes(RadarJPL{Float64}))

    recovered_Apophis = Vector{RadarJPL{Float64}}(df_Apophis)

    @test recovered_Apophis == radar_Apophis

    # Query RadarJPL
    _251_before_2013 = radar_Apophis |> @filter(date(_) < Date(2013) && rcvr(_) == search_obs_code("251")) |> DataFrame

    @test nrow(_251_before_2013) == 5
    @test all(_251_before_2013.date .< Date(2013))
    @test allequal(_251_before_2013.rcvr) && _251_before_2013.rcvr[1] == search_obs_code("251")

end