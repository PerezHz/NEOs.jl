using NEOs
using Tables
using DataFrames
using Test

using NEOs: read_radec_mpc, RadecMPC, read_radar_jpl, RadarJPL

@testset "DataFramesExt" begin

    # RadecMPC
    radec_2023DW = read_radec_mpc(joinpath("data", "RADEC_2023_DW.dat"))

    df_2023DW = DataFrame(radec_2023DW)

    @test nrow(df_2023DW) == length(radec_2023DW)
    @test all(names(df_2023DW) .== String.(fieldnames(RadecMPC{Float64})))

    recovered_2023DW = RadecMPC(df_2023DW)

    @test recovered_2023DW == radec_2023DW

    # RadarJPL
    radar_Apophis = read_radar_jpl(joinpath("data", "99942_RADAR_2005_2013.dat"))

    df_Apophis = DataFrame(radar_Apophis)

    @test nrow(df_Apophis) == length(radar_Apophis)
    @test all(names(df_Apophis) .== String.(fieldnames(RadarJPL{Float64})))

    recovered_Apophis = RadarJPL(df_Apophis)

    @test recovered_Apophis == radar_Apophis
end