"""
    ObservatoryMPC{T} <: AbstractAstrometryObservatory{T}

An observatory recognized by the Minor Planet Center.

# Fields

- `code::String`: three-character observatory code.
- `frame::String`: coordinate frame for the observatory position.
- `coords::SVector{3, T}`: position of the observatory in the aforementioned frame.
- `name::String`: name of the observatory, ASCII formatting.
- `uses_two_line_observations::Bool`: flag of whether the observatory uses "two line"
    observations (most do not).
- `observations_type::String`: the observation format that the observatory provides;
    one of 'optical', 'occultation', 'satellite', 'radar' or 'roving'.

!!! reference
    The list of observatories is available at:
    - https://minorplanetcenter.net/iau/lists/ObsCodesF.html
    - https://minorplanetcenter.net/mpcops/documentation/valid-ades-values/#stn
    The Minor Planet Center observatories API is described at:
    - https://www.minorplanetcenter.net/mpcops/documentation/obscodes-api/

# Extended help

The Minor Planet Center recognizes six different reference frames to report
the position of an observatory, all which are described below.

For ground-based static observatories:

- `MPC`: the default geocentric reference frame used by the MPC. Coordinates
    represent longitude East of the prime meridian [deg] and parallax constants
    `rho cos(phi')` and `rho sin(phi')`, where `phi'` is the geocentric latitude
    and `rho` is the geocentric distance [earth radii].

For ground-based roving observatories:

- `WGS84`: geodetic reference ellipsoid, GPS coordinates are normally obtained in
    this frame. Coordinates represent east longitude [deg], latitude [deg] and
    altitude [m].
- `ITRF`: cylindrical International Terrestrial Reference Frame. Coordinates
    represent east longitude [deg], Rxy [km] and Rz [km].
- `IAU`: IAU planetary cartographic model for bodies other than Earth. Coordinates
    represent longitude [deg], latitude [deg] and altitude [m], as defined by
    corresponding IAU cartography standard.

For space-based observatories:

- `ICRF_AU` or `ICRF_KM`: cartesian IAU International Celestial Reference Frame.
    Coordinates represent equatorial rectangular coordinates [au or km].
"""
@auto_hash_equals fields = (code,) struct ObservatoryMPC{T} <: AbstractAstrometryObservatory{T}
    code::String
    frame::String
    coords::SVector{3, T}
    name::String
    uses_two_line_observations::Bool
    observations_type::String
end

# Order in ObservatoryMPC is given by code
isless(a::ObservatoryMPC, b::ObservatoryMPC) = a.code < b.code

istwoliner(o::ObservatoryMPC) = o.uses_two_line_observations

isoptical(o::ObservatoryMPC) = o.observations_type == "optical"
isoccultation(o::ObservatoryMPC)= o.observations_type == "occultation"
issatellite(o::ObservatoryMPC) = o.observations_type == "satellite"
isradar(o::ObservatoryMPC) = o.observations_type == "radar"
isroving(o::ObservatoryMPC) = o.observations_type == "roving"

isgeocentric(o::ObservatoryMPC) = o.code == "500" || o.code == "244"
hascoord(o::ObservatoryMPC) = all(!isnan, o.coords)

isunknown(o::ObservatoryMPC) = isempty(o.code)
unknownobs(::Type{T} = Float64) where {T <: Real} = ObservatoryMPC{T}("", "",
    SVector{3, T}(NaN, NaN, NaN), "", false, "")

# Print method for ObservatoryMPC
function show(io::IO, o::ObservatoryMPC)
    if isunknown(o)
        print(io, "Unknown observatory")
    elseif hascoord(o)
        print(io, o.name, " [", o.code, "] coords: ", o.coords)
    else
        print(io, o.name, " [", o.code, "]")
    end
end

# Constructors
function ObservatoryMPC(obs::ObservatoryMPC{T};
    code::String = obs.code,
    frame::String = obs.frame,
    coords::SVector{3, T} = obs.coords,
    name::String = obs.name,
    uses_two_line_observations::Bool = obs.uses_two_line_observations,
    observations_type::String = obs.observations_type
    ) where {T <: Real}
    return ObservatoryMPC{T}(code, frame, coords, name,
        uses_two_line_observations, observations_type)
end

function ObservatoryMPC(pair)
    code = astrometryparse(String, first(pair))
    dict = last(pair)
    name = astrometryparse(String, dict["name"])
    uses_two_line_observations = astrometryparse(Bool, dict["uses_two_line_observations"])
    observations_type = astrometryparse(String, dict["observations_type"])
    if observations_type == "optical" || observations_type == "radar"
        frame = "MPC"
    elseif observations_type == "occultation" || observations_type == "satellite"
        frame = "ICRF"
    else # observations_type == "roving"
        frame = "WGS84"
    end
    x = astrometryparse(Float64, dict["longitude"])
    y = astrometryparse(Float64, dict["rhocosphi"])
    z = astrometryparse(Float64, dict["rhosinphi"])
    coords = SVector{3, Float64}(x, y, z)
    if uses_two_line_observations && all(iszero, coords)
        coords = SVector{3, Float64}(NaN, NaN, NaN)
    end

    return ObservatoryMPC{Float64}(code, frame, coords, name,
        uses_two_line_observations, observations_type)
end

# Parsing
function parse_observatories_mpc(text::AbstractString)
    # Parse JSON
    dict = JSON.parse(text)
    # Parse observatories
    obs = Vector{ObservatoryMPC{Float64}}(undef, length(dict))
    for (i, d) in enumerate(dict)
        obs[i] = ObservatoryMPC(d)
    end
    # Eliminate repeated entries
    unique!(obs)
    # Sort by three-character code
    sort!(obs)

    return obs
end

function read_observatories_mpc(filename::AbstractString)
    # Read file
    text = read(filename, String)
    # Parse observatories
    obs = parse_observatories_mpc(text)

    return obs
end

"""
    update_observatories_mpc()

Update the local observatories list.
"""
function update_observatories_mpc()
    # Download and parse observatories file
    text = fetch_http_text(MPC; mode = 1)
    # Parse observatories
    obs = parse_observatories_mpc(text)
    # Write observatories to local file
    path = joinpath(SCRATCH_PATH[], "observatoriesmpc.json")
    open(path, "w") do file
        write(file, text)
    end
    # Update global variable
    global OBSERVATORIES_MPC[] = obs

    return nothing
end

"""
    search_observatory_code(code::AbstractString)

Return the observatory in `NEOs.OBSERVATORIES_MPC` that matches `code`.
"""
function search_observatory_code(code::AbstractString)
    i = findfirst(x -> x.code == code, OBSERVATORIES_MPC[])
    return isnothing(i) ? unknownobs() : OBSERVATORIES_MPC[][i]
end

"""
    fetch_observatory_information(code::AbstractString)

Return the information given by the Minor Planet Center observatories
API about the observatory that matches `code`.

!!! reference
    The Minor Planet Center observatories API is described at:
    - https://www.minorplanetcenter.net/mpcops/documentation/obscodes-api/
"""
function fetch_observatory_information(code::AbstractString)
    # Parse HTTP response as String
    text = fetch_http_text(MPC; mode = 1, id = code)
    # Parse JSON
    dict = JSON.parse(text)

    return dict
end