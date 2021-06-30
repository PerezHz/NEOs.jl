# this script generates Julia artifacts for debiasing tables from
# Chesley et al. (2010), Farnocchia et al. (2015) and Eggl et al. (2020)
# TODO: adapt script for julia1.6+
# TODO: adapt script for 2010 debiasing tables (requires extra unpacking)

using Pkg.Artifacts
using Pkg.PlatformEngines
using Pkg.PlatformEngines: sha256
using Pkg
Pkg.PlatformEngines.probe_platform_engines!()

# This is the path to the Artifacts.toml we will manipulate
artifact_toml = joinpath(dirname(@__DIR__), "Artifacts.toml")

ftp_jpl = "ftp://ssd.jpl.nasa.gov/pub/ssd/debias/"

names_debias = ["debias", "debias_2014", "debias_2018", "debias_hires2018"]
urls_debias = ftp_jpl .* names_debias .* ".tgz"

names_lsk_spk = ["naif0012", "de430", "TTmTDBde430", "a99942", "ttmtdb_DE430_2003_2030"]
urls_lsk_spk = [
    "https://raw.githubusercontent.com/PerezHz/jpleph/main/naif0012.tar.gz",
    "https://raw.githubusercontent.com/PerezHz/jpleph/main/de430.tar.gz",
    "https://raw.githubusercontent.com/PerezHz/jpleph/main/TTmTDBde430.tar.gz",
    "https://raw.githubusercontent.com/PerezHz/jpleph/main/a99942.tar.gz",
    "https://raw.githubusercontent.com/PerezHz/jpleph/main/ttmtdb_DE430_2003_2030.tar.gz"
]

urls_ = vcat(urls_lsk_spk, urls_debias)
names_ = vcat(names_lsk_spk, names_debias)

for (url, name) in map((x,y)->(x,y), urls_, names_)
    # Query the `Artifacts.toml` file for the hash bound to the name
    # (returns `nothing` if no such binding exists)
    artfct_hash = artifact_hash(name, artifact_toml)

    # If the name was not bound, or the hash it was bound to does not exist, create it!
    if artfct_hash == nothing || !artifact_exists(artfct_hash)
        tb = Tuple[]
        # create_artifact() returns the content-hash of the artifact directory once we're finished creating it
        artfct_hash = create_artifact() do artifact_dir
            # We create the artifact by simply downloading a few files into the new artifact directory
            @show url
            tarball = download(url, joinpath(artifact_dir, basename(url)))
            try
                tarball_hash = open(tarball) do file
                    bytes2hex(sha256(file))
                end
                @show tarball_hash
                tb = [(url, tarball_hash)]
                unpack(tarball, artifact_dir)
            finally
                rm(tarball)
            end
        end
        # Now bind that hash within our `Artifacts.toml`. `force = true` means that if it already exists,
        # just overwrite the corresponding entry in `Artifacts.toml` with the new content-hash.
        # Unless the source files change, we do not expect the content hash to
        # change, so this should not cause unnecessary version control churn.
        bind_artifact!(artifact_toml, name, artfct_hash, lazy=true, force=true, download_info=tb)
    end

    # Get the path of the artifact, either newly created or previously generated.
    # this should be something like `~/.julia/artifacts/dbd04e28be047a54fbe9bf67e934be5b5e0d357a`
    artfct_path = artifact_path(artfct_hash)
    @show artfct_path
end
