# this script generates Julia artifacts for debiasing tables from
# Chesley et al. (2010), Farnocchia et al. (2015) and Eggl et al. (2020)
# TODO: adapt script for 2010 debiasing tables (requires extra unpacking)

using Pkg.Artifacts
using Pkg.PlatformEngines
using Pkg.PlatformEngines: sha256
using Pkg
Pkg.PlatformEngines.probe_platform_engines!()

# This is the path to the Artifacts.toml we will manipulate
artifact_toml = joinpath(pathof(Apophis), "Artifacts.toml")

ftp_jpl = "ftp://ssd.jpl.nasa.gov/pub/ssd/debias/"

names = ["debias", "debias_2014", "debias_2018", "debias_hires2018"]
urls =ftp_jpl .* names .* ".tgz"

for (url, name) in map((x,y)->(x,y), urls, names)
    # Query the `Artifacts.toml` file for the hash bound to the name
    # (returns `nothing` if no such binding exists)
    artfct_hash = artifact_hash(name, artifact_toml)

    # If the name was not bound, or the hash it was bound to does not exist, create it!
    if artfct_hash == nothing || !artifact_exists(artfct_hash)
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
                global tb = [(url, tarball_hash)]
                ### unpack(tarball, joinpath(artifact_dir, name))
                unpack(tarball, artifact_dir)
            finally
                rm(tarball)
            end
        end
        # Now bind that hash within our `Artifacts.toml`.  `force = true` means that if it already exists,
        # just overwrite with the new content-hash.  Unless the source files change, we do not expect
        # the content hash to change, so this should not cause unnecessary version control churn.
        bind_artifact!(artifact_toml, name, artfct_hash, lazy=true, download_info=tb) ####, force=true
    end

    # Get the path of the artifact, either newly created or previously generated.
    # this should be something like `~/.julia/artifacts/dbd04e28be047a54fbe9bf67e934be5b5e0d357a`
    artfct_path = artifact_path(artfct_hash)
    @show artfct_path
end
