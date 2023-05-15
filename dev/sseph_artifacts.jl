using Pkg.Artifacts
# using Pkg.Artifacts: artifact_exists

# URL of Solar System 2000-2100 ephemeris 
const sseph_p100_url = "https://github.com/LuEdRaMo/sseph/raw/main/sseph343ast016_p100y_et.jld2"

# Path to Artifacts.toml
artifact_toml = joinpath(dirname(@__DIR__), "Artifacts.toml")

# Query the `Artifacts.toml` file for the hash bound to the name "sseph_p100"
# (returns `nothing` if no such binding exists)
sseph_p100_hash = artifact_hash("sseph_p100", artifact_toml)

# If the name was not bound, or the hash it was bound to does not exist, create it!

if isnothing(sseph_p100_hash) || !artifact_exists(sseph_p100_hash)
    # create_artifact() returns the content-hash of the artifact directory once we're finished creating it
    sseph_p100_hash = create_artifact() do artifact_dir
        # We create the artifact by simply downloading the file into the new artifact directory
        download(sseph_p100_url, joinpath(artifact_dir, "sseph343ast016_p100y_et.jld2"))
    end
    # Now bind that hash within our `Artifacts.toml`.  `force = true` means that if it already exists,
    # just overwrite with the new content-hash.  Unless the source files change, we do not expect
    # the content hash to change, so this should not cause unnecessary version control churn.
    bind_artifact!(artifact_toml, "sseph_p100", sseph_p100_hash, lazy = true, force = true, 
                   download_info = Tuple[(sseph_p100_url, string(sseph_p100_hash))])
end