# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Test

@testset "Read/write catalogues file" begin

    cats_1 = read_catalogues_mpc(NEOs.CatalogueCodes_path)
    check_file = joinpath(NEOs.src_path, "observations/CatalogueCodes_.txt")
    write_catalogues_mpc(cats_1, check_file)
    cats_2 = read_catalogues_mpc(check_file)
    rm(check_file)

    @test cats_1 == cats_2
end

@testset "Update catalogues file" begin
    
    update_catalogues_mpc()
    @test isa(NEOs.mpc_catalogues[], Vector{NEOs.CatalogueMPC}) 

end

@testset "Read/write observatories file" begin

    obs_1 = read_observatories_mpc(NEOs.ObsCodes_path)
    check_file = joinpath(NEOs.src_path, "observations/ObsCodes_.txt")
    write_observatories_mpc(obs_1, check_file)
    obs_2 = read_observatories_mpc(check_file)
    rm(check_file)

    @test obs_1 == obs_2
end

@testset "Update observatories file" begin
    
    update_observatories_mpc()
    @test isa(NEOs.mpc_observatories[], Vector{NEOs.ObservatoryMPC{Float64}}) 

end

@testset "Read/write radec file" begin

    source_file = joinpath(dirname(NEOs.src_path), "data/99942.dat")
    obs_1 = read_radec_mpc(source_file)
    target_file = joinpath(dirname(NEOs.src_path), "data/99942_.dat")
    write_radec_mpc(obs_1, target_file)
    obs_2 = read_radec_mpc(target_file)
    rm(target_file)

    @test obs_1 == obs_2
end

@testset "Search MPC circulars" begin
    
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
        max_iter = 100
    )

    @test isa(obs, Vector{NEOs.RadecMPC{Float64}})
    @test issorted(obs)
    @test allunique(obs)
end

@testset "Read/write radar file" begin

    source_file = joinpath(dirname(NEOs.src_path), "data/99942_RADAR_2005_2013.dat")
    radar_1 = read_radar_jpl(source_file)
    target_file = joinpath(dirname(NEOs.src_path), "data/99942_RADAR_2005_2013_.dat")
    write_radar_jpl(radar_1, target_file)
    radar_2 = read_radar_jpl(target_file)
    rm(target_file)

    @test radar_1 == radar_2
end
