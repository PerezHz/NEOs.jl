function __init__()
    # Initialize scratch space
    global scratch_path[] = @get_scratch!("NEOsScratch")
    # Load catalogues 
    CatalogueCodes_path = joinpath(scratch_path[], "CatalogueCodes.txt")
    if isfile(CatalogueCodes_path)
        global mpc_catalogues[] = read_catalogues_mpc(CatalogueCodes_path)
    else 
        update_catalogues_mpc()
    end 
    # Load observatories 
    ObsCodes_path = joinpath(scratch_path[], "ObsCodes.txt")
    if isfile(ObsCodes_path)
        global mpc_observatories[] = read_observatories_mpc(ObsCodes_path)
    else 
        update_observatories_mpc()
    end 
    # Extensions 
    @static if !isdefined(Base, :get_extension)
        @require Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c" begin
            @require DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0" include("../ext/DataFramesExt.jl")
        end
        @require Query = "1a8c2f83-1ff3-5112-b086-8aa67b057ba1" include("../ext/QueryExt.jl")
    end
end