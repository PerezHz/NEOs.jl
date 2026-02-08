"""
    MPC <: AbstractAstrometrySource

Type representing the Minor Planet Center's APIs, except for those
related to the NEO Confirmation Page.

See also [`NEOCP`](@ref).
"""
struct MPC <: AbstractAstrometrySource end

"""
    NEOCP <: AbstractAstrometrySource

Type representing the Minor Planet Center's APIs related to the
NEO Confirmation Page.

See also [`MPC`](@ref).
"""
struct NEOCP <: AbstractAstrometrySource end

"""
    NEOCC <: AbstractAstrometrySource

Type representing the Near-Earth Objects Coordination Centre
automated data access APIs.

See also [`NEODyS2`](@ref).
"""
struct NEOCC <: AbstractAstrometrySource end

"""
    NEODyS2 <: AbstractAstrometrySource

Type representing the Near-Earth Objects Dynamic Site 2 APIs.

See also [`NEOCC`](@ref).
"""
struct NEODyS2 <: AbstractAstrometrySource end

"""
    JPL <: AbstractAstrometrySource

Type representing the Jet Propulsion Laboratory Solar System
Dynamics APIs.
"""
struct JPL <: AbstractAstrometrySource end

function get_http_response(::Type{MPC}; mode::Int = 0, id::AbstractString = "",
                           ids::AbstractVector{<:AbstractString} = [""],
                           format::AbstractString = "", connect_timeout = 180,
                           readtimeout = 180, kwargs...)
    @assert -1 <= mode <= 2 "`mode`  must be an integer between -1 and 2"
    # Designations
    if mode == -1
        params = JSON.json(Dict("ids" => ids))
        url = DESIGNATIONS_MPC_API
    # CatalogueMPC
    elseif mode == 0
        params = JSON.json(Dict())
        url = CATALOGUES_MPC_FILE_URL
    # ObservatoryMPC
    elseif mode == 1
        params = isempty(id) ? JSON.json(Dict()) : JSON.json(Dict("obscode" => id))
        url = OBSERVATORIES_MPC_API
    # OpticalMPC80 / OpticalADES
    else
        params = JSON.json(Dict("desigs" => [id], "output_format" => [format]))
        url = OBSERVATIONS_MPC_API
    end
    # HTTP response (HTTP.get retries four times by default)
    resp = HTTP.get(url, ("Content-Type" => "application/json",), params;
        connect_timeout, readtimeout, kwargs...)

    return resp
end

# NEOCPObject
function get_http_response(::Type{NEOCP}; mode::Int = 0, id::AbstractString = "",
                           format::AbstractString = "", connect_timeout = 180,
                           readtimeout = 180, kwargs...)
    # All objects listed in the NEOCP
    if mode == 0
        resp = HTTP.get(NEOCP_OBJECTS_FILE_URL; connect_timeout, readtimeout, kwargs...)
    # Specific object from the NEOCP
    else
        params = JSON.json(Dict("trksubs" => [id], "output_format" => [format]))
        resp = HTTP.get(OBSERVATIONS_NEOCP_API, ("Content-Type" => "application/json",),
            params; connect_timeout, readtimeout, kwargs...)
    end

    return resp
end

# OpticalMPC80 / OpticalADES / RadarRWO
function get_http_response(::Type{NEOCC}; id::AbstractString, connect_timeout = 180,
                           readtimeout = 180, kwargs...)
    url = OBSERVATIONS_NEOCC_API * id * ".rwo"
    resp = HTTP.get(url; connect_timeout, readtimeout, kwargs...)

    return resp
end

# OpticalMPC80 / OpticalADES / RadarRWO
function get_http_response(::Type{NEODyS2}; id::AbstractString, connect_timeout = 180,
                           readtimeout = 180, kwargs...)
    url = OBSERVATIONS_NEODyS2_API * id * ".rwo"
    resp = HTTP.get(url; connect_timeout, readtimeout, kwargs...)

    return resp
end

# RadarJPL
function get_http_response(::Type{JPL}; id::Pair{String, String},
                           connect_timeout = 180, readtimeout = 180, kwargs...)
    query = isempty(first(id)) ? Dict{String, String}() : Dict{String, String}(id)
    resp = HTTP.get(RADAR_JPL_API; query, connect_timeout, readtimeout,
        status_exception = false, require_ssl_verification = false, kwargs...)

    return resp
end

function fetch_http_text(::Type{S}; kwargs...) where {S <: AbstractAstrometrySource}
    # HTTP response
    resp = get_http_response(S; kwargs...)
    # Convert to String
    text = String(resp.body)

    return text
end

# Parsing

astrometrydefault(::Type{Int}) = 0
astrometrydefault(::Type{Char}) = ' '
astrometrydefault(::Type{String}) = ""
astrometrydefault(::Type{Bool}) = false
astrometrydefault(::Type{DateTime}) = DateTime(2000, 1, 1)
astrometrydefault(::Type{T}) where {T <: Real} = T(NaN)
astrometrydefault(::Type{CatalogueMPC}) = unknowncat()
astrometrydefault(::Type{ObservatoryMPC{T}}) where {T <: Real} = unknownobs(T)

astrometryparse(::Type{T}, x::T) where {T} = x
astrometryparse(::Type{Bool}, x::AbstractString) = x == "Yes"
astrometryparse(::Type{T}, x::AbstractString) where {T <: Real} = parse(T, x)