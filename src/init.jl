# List of star catalogues recognized by the MPC
const CATALOGUES_MPC = Ref(CatalogueMPC[])

# List of observatories recognized by the MPC
const OBSERVATORIES_MPC = Ref(ObservatoryMPC{Float64}[])

function __init__()
    # Initialize scratch space
    global SCRATCH_PATH[] = @get_scratch!("NEOs")
    # Load catalogues
    path = joinpath(SCRATCH_PATH[], "astCat_photCat.json")
    if !isfile(path)
        write(path, read(CATALOGUES_PATH))
    end
    global CATALOGUES_MPC[] = read_catalogues_mpc(path)
    # Load observatories
    path = joinpath(SCRATCH_PATH[], "observatoriesmpc.json")
    if isfile(path)
        global OBSERVATORIES_MPC[] = read_observatories_mpc(path)
    else
        update_observatories_mpc()
    end
end