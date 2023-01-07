# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Test

using NEOs: mpc_catalogue_regex, CatalogueMPC, CatalogueCodes_path, observations_path

@testset "CatalogueMPC" begin  

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
    update_catalogues_mpc()
    @test isa(NEOs.mpc_catalogues[], Vector{NEOs.CatalogueMPC}) 

    # Search catalogue code 
    cat = search_cat_code("6")
    @test cat == gaia
end 

using NEOs: mpc_observatory_regex, ObservatoryMPC, ObsCodes_path

@testset "ObservatoryMPC" begin
    
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
    update_observatories_mpc()
    @test isa(NEOs.mpc_observatories[], Vector{NEOs.ObservatoryMPC{Float64}}) 

    # Search observatory code 
    obs = search_obs_code("250")
    @test obs == hubble
    obs = search_obs_code("251")
    @test obs == arecibo
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
