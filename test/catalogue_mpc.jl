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