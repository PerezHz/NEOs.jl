# Astrometry sources
struct MPC <: AbstractAstrometrySource end
struct NEOCP <: AbstractAstrometrySource end
struct NEOCC <: AbstractAstrometrySource end
struct NEODyS2 <: AbstractAstrometrySource end
struct JPL <: AbstractAstrometrySource end

function get_http_response(::Type{MPC}; mode::Int = 0, id::AbstractString = "",
                           format::AbstractString = "", connect_timeout = 180,
                           readtimeout = 180, kwargs...)
    @assert 0 <= mode <= 2 "`mode`  must be an integer between 0 and 2"
    # CatalogueMPC
    if mode == 0
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

function get_http_response(::Type{NEOCP}; mode::Int = 0, id::AbstractString = "",
                           format::AbstractString = "", connect_timeout = 180,
                           readtimeout = 180, kwargs...)
    # NEOCPObject
    if mode == 0
        resp = HTTP.get(NEOCP_OBJECTS_FILE_URL; connect_timeout, readtimeout, kwargs...)
    else
        params = JSON.json(Dict("trksubs" => [id], "output_format" => [format]))
        resp = HTTP.get(OBSERVATIONS_NEOCP_API, ("Content-Type" => "application/json",),
            params; connect_timeout, readtimeout, kwargs...)
    end

    return resp
end

function get_http_response(::Type{NEOCC}; id::AbstractString, connect_timeout = 180,
                           readtimeout = 180, kwargs...)
    url = OBSERVATIONS_NEOCC_API * id * ".rwo"
    resp = HTTP.get(url; connect_timeout, readtimeout, kwargs...)

    return resp
end

function get_http_response(::Type{NEODyS2}; id::AbstractString, connect_timeout = 180,
                           readtimeout = 180, kwargs...)
    url = OBSERVATIONS_NEODyS2_API * id * ".rwo"
    resp = HTTP.get(url; connect_timeout, readtimeout, kwargs...)

    return resp
end

function get_http_response(::Type{JPL}; id::Pair{String, String},
                           connect_timeout = 180, readtimeout = 180, kwargs...)
    query = Dict(id)
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