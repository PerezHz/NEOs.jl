function __init__()
    # Initialize scratch space
    global scratch_path[] = @get_scratch!("NEOsScratch")
    # Load catalogues 
    CatalogueCodes_path = joinpath(scratch_path[], "CatalogueCodes.txt")
    if isfile(CatalogueCodes_path)
        global CATALOGUES_MPC[] = read_catalogues_mpc(CatalogueCodes_path)
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
end