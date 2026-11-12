# Solar system ephemeris source file
const SSEPH_ARTIFACT_PATH = joinpath(artifact"sseph_p100", "sseph343ast016_p100y_et.jld2")
const SSEPH_SOURCE = @load_preference("SSEPH_SOURCE", SSEPH_ARTIFACT_PATH)
# Order of the Taylor polynomials of the solar system ephemeris
const SSEPHORDER::Int = @load_preference("SSEPHORDER", 25)

# NAIF ids of the Sun, Mercury, Venus, the Earth, the Moon, the barycenters of Mars,
# Jupiter, Saturn, Uranus and Neptune, and the barycenter of Pluto; i.e. the first 11
# bodies of every solar system ephemeris
const SSEPH_PLANETS_SPKIDS = [10, 199, 299, 399, 301, 4, 5, 6, 7, 8, 9]

"""
    set_sseph_source(filename)

Set the solar system ephemeris source file to `filename`.

# Extended help

`filename` must be a valid `.jld2` file produced by PlanetaryEphemeris.jl (e.g. by its
`scripts/de430.jl` or `scripts/spice_to_pe.jl` scripts), meaning that it contains the
following fields:
- `sseph`: solar system ephemeris. Its first 11 bodies must be the Sun, Mercury,
    Venus, the Earth, the Moon, and the barycenters of Mars, Jupiter, Saturn, Uranus,
    Neptune and Pluto (in that order). Files produced by PlanetaryEphemeris.jl versions
    older than 0.14 (e.g. the default artifact) store it as `ss16ast_eph`, which is
    also accepted.
- `spkids` (optional): NAIF ids of the bodies in `sseph`.
- `μ` (optional): gravitational parameters [au³/day²] of the bodies in `sseph`.
- `model` (optional): name of the ephemeris model (e.g. `"DE440"`).
- `params` (optional): model parameters, with the names and units of JPL's technical
    comments: `"RE"` [km], `"J2E"` and `"J3E"`.
- `acc_eph` (optional): acceleration ephemeris.
- `pot_eph` (optional): Newtonian potential ephemeris.

`spkids` and `μ` must be given together; if both are missing, the bodies and
gravitational parameters of the JPL DE430 model (as integrated by PlanetaryEphemeris.jl)
are assumed. `model` and `params` must also be given together; if both are missing,
the file must come from the JPL DE430 model, whose parameters are then assumed.
If `acc_eph` and `pot_eph` are missing, or their order is smaller than
`SSEPHORDER` (see `set_sseph_order`), they are computed when NEOs is loaded and
then saved to `filename`.

See also [`set_sseph_order`](@ref).
"""
function set_sseph_source(filename::AbstractString)
    @set_preferences!("SSEPH_SOURCE" => filename)
    @info("Solar system ephemeris source file set to $filename; \
          restart your Julia session for this change to take effect!")
end

"""
    set_sseph_order(order)

Set the order of the Taylor polynomials of the solar system
ephemeris to `order` (default: `25`).

# Extended help

When NEOs is loaded, the polynomials in the source file
(see [`set_sseph_source`](@ref)) are truncated, if their
order is larger than `order`, or padded with zeros otherwise.
"""
function set_sseph_order(order::Integer)
    order ≥ 1 || throw(ArgumentError("The solar system ephemeris order must be positive"))
    @set_preferences!("SSEPHORDER" => Int(order))
    @info("Solar system ephemeris order set to $order; \
          restart your Julia session for this change to take effect!")
end

# Model parameters required in the solar system ephemeris source file
const SSEPH_PARAMS_NAMES = ["RE", "J2E", "J3E"]
# Units of the model parameters (dimensionless parameters are omitted)
const SSEPH_PARAMS_UNITS = Dict("RE" => "km")

# Ephemeris model and model parameters of the JPL DE430 model
const DE430_MODEL = "DE430"
const DE430_PARAMS = Dict{String, Float64}("RE" => PE.RE, "J2E" => PE.J2E,
    "J3E" => PE.J3E)

# NAIF ids of the first `N` bodies of the JPL DE430 model integrated by
# PlanetaryEphemeris.jl (e.g. the ephemeris in SSEPH_ARTIFACT_PATH)
de430_spkids(N::Int) = vcat(SSEPH_PLANETS_SPKIDS,
    2_000_000 .+ PE.de430_343ast_ids[1:N-length(SSEPH_PLANETS_SPKIDS)])

# Read `key` from the open .jld2 `file` if it is a dense TaylorSolution whose
# polynomials have order ≥ SSEPHORDER, and return it with polynomials of order
# SSEPHORDER; otherwise, return `nothing`
function readsseph(file, key::AbstractString)
    haskey(file, key) || return nothing
    eph::DensePropagation2{Float64, Float64} = file[key]
    all(x -> TaylorSeries.order(x) ≥ SSEPHORDER, eph.p) || return nothing
    return changeorder(eph, SSEPHORDER)
end

# Read the solar system ephemeris (with polynomials of order SSEPHORDER), the NAIF ids
# and gravitational parameters of its bodies and, if available, the accelerations and
# Newtonian potentials ephemerides (with polynomials of order SSEPHORDER) from the open
# .jld2 file `file`
function readssephsource(file)
    # Solar system ephemeris key: `sseph` (PlanetaryEphemeris.jl ≥ 0.14) or
    # `ss16ast_eph` (older files, e.g. the default artifact)
    key = haskey(file, "sseph") ? "sseph" : "ss16ast_eph"
    @assert haskey(file, key) "Solar system ephemeris source file has no key \
        \"sseph\" (or \"ss16ast_eph\")"
    # Solar system ephemeris, with polynomials of order SSEPHORDER
    sseph::DensePropagation2{Float64, Float64} = changeorder(file[key], SSEPHORDER)
    N = numberofbodies(sseph)
    # NAIF ids and gravitational parameters [au³/day²] of the bodies in sseph
    @assert haskey(file, "spkids") == haskey(file, "μ") "Solar system ephemeris \
        source file must contain both or none of the keys \"spkids\" and \"μ\""
    if haskey(file, "spkids")
        spkids = Vector{Int}(file["spkids"])
        gms = Vector{Float64}(file["μ"])
    else
        # Bodies and gravitational parameters of the JPL DE430 model
        spkids = de430_spkids(N)
        gms = PE.μ[1:N]
    end
    @assert length(spkids) == length(gms) == N "Solar system ephemeris source file \
        has $N bodies, but $(length(spkids)) NAIF ids and $(length(gms)) \
        gravitational parameters"
    @assert spkids[1:length(SSEPH_PLANETS_SPKIDS)] == SSEPH_PLANETS_SPKIDS "The \
        first $(length(SSEPH_PLANETS_SPKIDS)) bodies of the solar system ephemeris \
        must have NAIF ids $SSEPH_PLANETS_SPKIDS"
    # Ephemeris model and model parameters
    @assert haskey(file, "model") == haskey(file, "params") "Solar system ephemeris \
        source file must contain both or none of the keys \"model\" and \"params\""
    if haskey(file, "model")
        model = String(file["model"])
        params = Dict{String, Float64}(file["params"])
        missingparams = setdiff(SSEPH_PARAMS_NAMES, keys(params))
        @assert isempty(missingparams) "Solar system ephemeris source file is missing \
            the model parameters $missingparams"
    else
        # Files without model parameters (PlanetaryEphemeris.jl < 0.14) are only
        # accepted if they come from the JPL DE430 model; i.e. if the gravitational
        # parameters of the Sun, the planets, the Moon and Pluto match
        n = length(SSEPH_PLANETS_SPKIDS)
        @assert isapprox(gms[1:n], PE.μ[1:n]; rtol = 1E-14) "Solar system ephemeris \
            source file has no model parameters and does not come from the JPL DE430 \
            model; convert it again with PlanetaryEphemeris.jl ≥ 0.14"
        model, params = DE430_MODEL, copy(DE430_PARAMS)
    end
    # Accelerations and Newtonian potentials, with polynomials of order SSEPHORDER
    # (only if they are stored with polynomials of order ≥ SSEPHORDER)
    acceph = readsseph(file, "acc_eph")
    poteph = readsseph(file, "pot_eph")
    return sseph, spkids, gms, model, params, acceph, poteph
end

# Replace the .jld2 file `filename` by a new file with the same name and contents,
# except for the keys `acc_eph` and `pot_eph`, which are set to `acceph` and `poteph`.
# A new file is written (instead of modifying `filename` in place) because JLD2 does
# not free the space of deleted keys
function savesseph(filename::AbstractString, acceph::DensePropagation2,
                   poteph::DensePropagation2)
    # Temporary file in the same directory, so the final move is a rename
    tmp = tempname(dirname(abspath(filename)); cleanup = false)
    try
        jldopen(filename, "r") do src
            jldopen(tmp, "w") do dst
                for key in keys(src)
                    key in ("acc_eph", "pot_eph") && continue
                    dst[key] = src[key]
                end
                dst["acc_eph"] = acceph
                dst["pot_eph"] = poteph
            end
        end
        mv(tmp, filename; force = true)
    finally
        isfile(tmp) && rm(tmp)
    end
    return nothing
end

function load_sseph()
    sseph, spkids, gms, model, params, acceph, poteph = jldopen(readssephsource,
        SSEPH_SOURCE, "r")
    if isnothing(acceph) || isnothing(poteph)
        @info "Computing acceph and poteph"
        acceph, poteph = PE.loadeph(sseph, gms)
        # Save them to SSEPH_SOURCE, so they are not computed again
        try
            savesseph(SSEPH_SOURCE, acceph, poteph)
            @info "Saved acceph and poteph to $SSEPH_SOURCE"
        catch e
            @warn "Could not save acceph and poteph to $SSEPH_SOURCE" exception = e
        end
    end
    ttmtdb::DensePropagation1{Float64, Float64} = TaylorSolution(sseph.t, sseph.p[:, end])
    return sseph, acceph, poteph, ttmtdb, spkids, gms, model, params
end

# Load Solar System, accelerations, newtonian potentials and TT-TDB ephemeris, the NAIF
# ids and gravitational parameters of its bodies, and the ephemeris model and its
# parameters
const sseph, acceph, poteph, ttmtdb, SSEPHSPKIDS, SSEPHGMS, SSEPHMODEL, SSEPHPARAMS =
    load_sseph()

const SSEPHNBODIES::Int = numberofbodies(sseph)
const MINDTTDB::DateTime = julian2datetime(firsttime(sseph) + JD_J2000)
const MAXDTTDB::DateTime = julian2datetime(lasttime(sseph) + JD_J2000)

"""
    print_sseph_summary()

Print a summary of the solar system ephemeris currently loaded by NEOs.
"""
function print_sseph_summary()
    t = repeat(' ', 4)
    print(
        "Solar system ephemeris summary\n",
        t, rpad("Expansion order:", 21), SSEPHORDER, "\n",
        t, rpad("Number of bodies:", 21), SSEPHNBODIES, "\n",
        t, rpad("Initial time:", 21), MINDTTDB, " TDB", "\n",
        t, rpad("Final time:", 21), MAXDTTDB, " TDB", "\n",
        t, rpad("Model:", 21), SSEPHMODEL, "\n",
        t, "Model parameters:\n",
    )
    # Required parameters first, then any other parameter in the file
    others = sort!(collect(setdiff(keys(SSEPHPARAMS), SSEPH_PARAMS_NAMES)))
    for name in vcat(SSEPH_PARAMS_NAMES, others)
        unit = haskey(SSEPH_PARAMS_UNITS, name) ? " " * SSEPH_PARAMS_UNITS[name] : ""
        print(t, t, rpad(name * ":", 17), SSEPHPARAMS[name], unit, "\n")
    end
end
