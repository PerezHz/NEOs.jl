# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Test

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