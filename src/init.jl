# List of star catalogues recognized by the MPC
const CATALOGUES_MPC = Ref(CatalogueMPC[])

# List of observatories recognized by the MPC
const OBSERVATORIES_MPC = Ref(ObservatoryMPC{Float64}[])

# List of visual magnitude bands recognized by the MPC
const MAGNITUDE_BANDS_MPC = Ref(MagnitudeBandMPC{Float64}[])

# Earth orientation parameters (EOP) 2000
const EOP_IAU2000A::EopIau2000A = download_iers_eop()

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
    if !isfile(path)
        write(path, read(OBSERVATORIES_PATH))
    end
    global OBSERVATORIES_MPC[] = read_observatories_mpc(path)
    # Load magnitude bands
    path = joinpath(SCRATCH_PATH[], "magnitudebandsmpc.json")
    if !isfile(path)
        write(path, read(MAGNITUDE_BANDS_PATH))
    end
    global MAGNITUDE_BANDS_MPC[] = read_magnitude_bands_mpc(path)
end