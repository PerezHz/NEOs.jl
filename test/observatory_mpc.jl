# This file is part of the NEOs.jl package; MIT licensed

using NEOs
using Test

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