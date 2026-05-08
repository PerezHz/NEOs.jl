"""
    MagnitudeBandMPC{T} <: AbstractAstrometryMagnitudeBand{T}

A visual magnitude band recognized by the Minor Planet Center.

# Fields

- `allowed::Bool`: true if currently allowed in submissions.
- `band::String`: abbreviated name of band.
- `notes::String`: brief description of band.
- `v_conversion::T`: V-conversion parameter.

!!! reference
    The Minor Planet Center magnitude bands API is described at:
    - https://docs.minorplanetcenter.net/mpc-ops-docs/apis/mag-band/
"""
@auto_hash_equals fields = (band,) struct MagnitudeBandMPC{T} <: AbstractAstrometryMagnitudeBand{T}
    allowed::Bool
    band::String
    notes::String
    v_conversion::T
end

# Order in MagnitudeBandMPC is given by band
isless(a::MagnitudeBandMPC, b::MagnitudeBandMPC) = a.band < b.band

isallowed(x::MagnitudeBandMPC) = x.allowed
isunknown(x::MagnitudeBandMPC) = x.band == "UNK"
unknownband(::Type{T} = Float64) where {T <: Real} =
    MagnitudeBandMPC{T}( true, "UNK", "Unknown", zero(T))

# Print method for MagnitudeBandMPC
show(io::IO, x::MagnitudeBandMPC) = print(io, x.band, " magnitude band ",
    isallowed(x) ? "(allowed)" : "(not allowed)")

# Constructor
function MagnitudeBandMPC(dict)
    allowed = astrometryparse(Bool, dict["allowed"])
    band = astrometryparse(String, dict["band"])
    notes = astrometryparse(String, dict["notes"])
    v_conversion = astrometryparse(Float64, dict["v_conversion"])

    return MagnitudeBandMPC{Float64}(allowed, band, notes, v_conversion)
end

# Parsing
function parse_magnitude_bands_mpc(text::AbstractString)
    # Parse JSON
    dict = JSON.parse(text)
    # Parse magnitude bands
    bands = Vector{MagnitudeBandMPC{Float64}}(undef, length(dict["response"]))
    for (i, d) in enumerate(dict["response"])
        bands[i] = MagnitudeBandMPC(d)
    end
    # Eliminate repeated entries
    unique!(bands)
    # Sort by three-character code
    sort!(bands)

    return bands
end

function read_magnitude_bands_mpc(filename::AbstractString)
    # Read file
    text = read(filename, String)
    # Parse magnitude bands
    bands = parse_magnitude_bands_mpc(text)

    return bands
end

"""
    update_magnitude_bands_mpc()

Update the local visual magnitude bands list.
"""
function update_magnitude_bands_mpc()
    # Download and parse magnitude bands file
    text = fetch_http_text(MPC; mode = -3)
    # Parse magnitude bands
    bands = parse_magnitude_bands_mpc(text)
    # Write magnitude bands to local file
    path = joinpath(SCRATCH_PATH[], "magnitudebandsmpc.json")
    open(path, "w") do file
        write(file, text)
    end
    # Update global variable
    global MAGNITUDE_BANDS_MPC[] = bands

    return nothing
end

"""
    search_magnitude_band(band::AbstractString)

Return the visual magnitude band in `NEOs.MAGNITUDE_BANDS_MPC` that matches `band`.
"""
function search_magnitude_band(band::AbstractString)
    i = findfirst(x -> x.band == band, MAGNITUDE_BANDS_MPC[])
    return isnothing(i) ? unknownband() : MAGNITUDE_BANDS_MPC[][i]
end

"""
    fetch_magnitude_band_information(band::AbstractString)

Return the information given by the Minor Planet Center magnitude
bands API about the visual magnitude band that matches `band`.

!!! reference
    The Minor Planet Center magnitude bands API is described at:
    - https://docs.minorplanetcenter.net/mpc-ops-docs/apis/mag-band/
"""
function fetch_magnitude_band_information(band::AbstractString)
    # Parse HTTP response as String
    text = fetch_http_text(MPC; mode = -3, band)
    # Parse JSON
    dict = JSON.parse(text)

    return dict
end