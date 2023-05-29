# This script generates Julia artifacts for NEOs.jl, including
# - debiasing tables from Chesley et al. (2010), Farnocchia et al. (2015) and Eggl et al. (2020)
# - JPL lsk and spk
# - PlanetaryEphemeris.jl Solar System ephemeris

# TO DO: adapt script for julia1.6+
# TO DO: adapt script for 2010 debiasing tables (requires extra unpacking)

using Pkg.Artifacts
using Pkg.PlatformEngines: sha256, unpack
using Downloads

# This is the path to the Artifacts.toml we will manipulate
const artifact_toml = joinpath(dirname(@__DIR__), "Artifacts.toml")

# JPL FTP URL
const ftp_jpl = "https://ssd.jpl.nasa.gov/ftp/"

# Debiasing tables names and URLs
const names_debias = ["debias", "debias_2014", "debias_2018", "debias_hires2018"]
const urls_debias = ftp_jpl .* "ssd/debias/" .* names_debias .* ".tgz"

# JPL lsk and spk names and URLs
const names_lsk_spk = ["naif0012", "de430", "TTmTDBde430", "a99942"]
const urls_lsk_spk = [
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
    ftp_jpl*"eph/planets/bsp/de430_1850-2150.bsp",
    ftp_jpl*"eph/planets/bsp/TTmTDB.de430.19feb2015.bsp",
    "https://raw.githubusercontent.com/PerezHz/jpleph/main/a99942.tar.gz",
]

# PlanetaryEphemeris.jl Solar System ephemeris name and URL
const name_sseph = "sseph_p100"
const url_sseph = "https://github.com/LuEdRaMo/sseph/raw/main/sseph343ast016_p100y_et.tar.gz"

# Collect names and URLs
const urls = vcat(urls_lsk_spk, urls_debias, url_sseph)
const names = vcat(names_lsk_spk, names_debias, name_sseph)

for (url, name) in zip(urls, names)
    # Query the `Artifacts.toml` file for the hash bound to the name
    # (returns `nothing` if no such binding exists)
    artfct_hash = artifact_hash(name, artifact_toml)
    # If the name was not bound, or the hash it was bound to does not exist, create it!
    if isnothing(artfct_hash) || !artifact_exists(artfct_hash)
        tb = Tuple[]
        # create_artifact() returns the content-hash of the artifact directory once we're finished creating it
        artfct_hash = create_artifact() do artifact_dir
            # We create the artifact by simply downloading a few files into the new artifact directory
            @show url name
            tmp_dir = mktempdir()
            basename_url = basename(url)
            @show basename_url
            dot_idx = findlast(isequal('.'), basename_url)
            ext = basename_url[dot_idx+1:end]
            if ext in ["gz", "tgz"]  # If artifact is a tarball with many files
                tarball = Downloads.download(url, joinpath(tmp_dir, basename(url)))
                try
                    tarball_hash = open(tarball) do file
                        bytes2hex(sha256(file))
                    end
                    unpack(tarball, artifact_dir)
                    tb = [(url, tarball_hash)]
                finally
                    rm(tarball)
                end
            else # url containing a single non-tarball file
                fpath = Downloads.download(url, joinpath(artifact_dir, basename(url)))
                file_hash = open(fpath) do file
                    bytes2hex(sha256(file))
                end
                tb = [(url, file_hash)]
            end
        end
        @show  artfct_hash typeof(artfct_hash)
        # Now bind that hash within our `Artifacts.toml`. `force = true` means that if it already exists,
        # just overwrite the corresponding entry in `Artifacts.toml` with the new content-hash.
        # Unless the source files change, we do not expect the content hash to
        # change, so this should not cause unnecessary version control churn.
        bind_artifact!(artifact_toml, name, artfct_hash, lazy=true, force=true, download_info=tb)
    end

    # Get the path of the artifact, either newly created or previously generated.
    # this should be something like `~/.julia/artifacts/dbd04e28be047a54fbe9bf67e934be5b5e0d357a`
    artfct_path = artifact_path(artfct_hash)
    println(name, " artifact saved to ", artfct_path)
end
