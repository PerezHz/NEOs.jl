"""
    RadarJPL{T} <: AbstractRadarAstrometry{T}

A radar astrometric observation in the JPL format.

# Fields

- `des::String`: primary designation of the object.
- `epoch::DateTime`: epoch of the measurement (UTC).
- `value::T`: radar delay [us] or Doppler measurement [Hz].
- `sigma::T`: estimated 1-sigma uncertainty of the measurement.
- `units::String`: units associated with the measurement; `us`
    (microseconds) for delay and `Hz` (hertz) for Doppler.
- `freq::T`: frequency of the transmitter [MHz].
- `rcvr::ObservatoryMPC{T}`: receiving station.
- `xmit::ObservatoryMPC{T}`: transmitting station.
- `bp::String`: measurement bounce point with possible values
    of `P` (peak power) and `C` (center of mass).

!!! reference
    The JPL format is described at:
    - https://ssd-api.jpl.nasa.gov/doc/sb_radar.html
"""
@auto_hash_equals fields = (epoch, value, units, rcvr, xmit) struct RadarJPL{T} <: AbstractRadarAstrometry{T}
    des::String
    epoch::DateTime
    value::T
    sigma::T
    units::String
    freq::T
    rcvr::ObservatoryMPC{T}
    xmit::ObservatoryMPC{T}
    bp::String
end

# AbstractAstrometryObservation interface
date(x::RadarJPL) = x.epoch
measure(x::RadarJPL) = x.value
frequency(x::RadarJPL) = x.freq
observatory(x::RadarJPL) = x.rcvr
rms(x::RadarJPL) = x.sigma
debias(::RadarJPL{T}) where {T <: Real} = zero(T)

isdelay(r::RadarJPL) = r.units == "us"
isdoppler(r::RadarJPL) = r.units == "Hz"
ismonostatic(r::RadarJPL) = r.rcvr == r.xmit
issband(freq::Real) = 2000.0 ≤ freq ≤ 4000.0
isxband(freq::Real) = 8000.0 ≤ freq ≤ 12000.0

# Print method for RadarJPL
show(io::IO, r::RadarJPL) = print(io, r.des, " value: ", @sprintf("%+e", r.value),
    " ", r.units, " sigma: ", @sprintf("%.3e", r.sigma), " ", r.units, " t: ",
    r.epoch, " rcvr: ", r.rcvr.name)

# Parsing

jplparse(_, ::Type{String}, x) = String(x)
jplparse(_, ::Type{T}, x) where {T <: Real} = parse(T, x)
jplparse(_, ::Type{DateTime}, x) = DateTime(x, RADAR_JPL_DATEFORMAT)
jplparse(_, ::Type{ObservatoryMPC{T}}, x) where {T <: Real} =
    search_observatory_code(JPL_TO_MPC_OBSCODES[x])

RadarJPL(r::DataFrameRow) = RadarJPL{Float64}(
    r.des, r.epoch, r.value, r.sigma, r.units, r.freq,
    r.rcvr, r.xmit, r.bp
)

function parse_radar_jpl(text::AbstractString)
    # Parse JSON
    dict = JSON.parse(text)
    # NOTE: As of Sep 18th 2025, Ceduna (30-m, UTAS) [-74] has no MPC code
    filter!(x -> x[7] != "-74" && x[8] != "-74", dict["data"])
    L = length(dict["data"])
    iszero(L) && return RadarJPL{Float64}[]
    # Construct DataFrame
    R = RadarJPL{Float64}
    names, types = fieldnames(R), fieldtypes(R)
    df = DataFrame([fill(astrometrydefault(fieldtype(R, name)), L) for name in names],
        collect(names))
    for (i, line) in enumerate(dict["data"])
        for (name, type, x) in zip(names, types, line)
            df[i, name] = jplparse(name, type, x)
        end
    end
    # Parse observations
    radar = RadarJPL.(eachrow(df))
    # Eliminate repeated entries
    unique!(radar)
    # Sort by date
    sort!(radar)

    return radar
end

"""
    fetch_radar_jpl([id, ] source)

Return the radar astrometry of minor body `id` in the JPL format.
The `source` of the observations can only be `JPL`. If `id` is
omitted, return all the radar astrometry available.

!!! reference
    The JPL Small-Body Radar Astrometry API is described at:
    - https://ssd-api.jpl.nasa.gov/doc/sb_radar.html
"""
function fetch_radar_jpl(::Type{JPL})
    # Get and parse HTTP response
    text = fetch_http_text(JPL; id = "" => "")
    # Parse observations
    radar = parse_radar_jpl(text)

    return radar
end

function fetch_radar_jpl(id::Pair{String, String}, ::Type{JPL})
    # Get and parse HTTP response
    text = fetch_http_text(JPL; id)
    # Parse observations
    radar = parse_radar_jpl(text)

    return radar
end

# Read / write

jplstr(x::String) = x
jplstr(x::Real) = string(x)
jplstr(x::ObservatoryMPC) = MPC_TO_JPL_OBSCODES[x.code]
jplstr(x::DateTime) = Dates.format(x, RADAR_JPL_DATEFORMAT)

"""
    read_radar_jpl(filename)

Read from `filename` a vector of radar astrometry in the JPL format.
"""
function read_radar_jpl(filename::AbstractString)
    # Read file
    text = read(filename, String)
    # Parse observations
    radar = parse_radar_jpl(text)

    return radar
end

"""
    write_radar_jpl(obs, filename)

Write `obs` to `filename` in the JPL format.
"""
function write_radar_jpl(obs::AbstractVector{RadarJPL{T}},
                         filename::AbstractString) where {T <: Real}
    names = fieldnames(RadarJPL{T})
    fields = [string(name) for name in names]
    signature = Dict("source"  => "NASA/JPL Small-Body Radar Astrometry API",
                     "version" => "1.1")
    count = string(length(obs))
    data = Vector{Vector{String}}(undef, length(obs))
    for i in eachindex(data)
        data[i] = [jplstr(getfield(obs[i], name)) for name in names]
    end
    dict = Dict("fields" => fields, "signature" => signature, "count" => count,
                "data" => data)
    open(filename, "w") do file
        write(file, JSON.json(dict))
    end
end