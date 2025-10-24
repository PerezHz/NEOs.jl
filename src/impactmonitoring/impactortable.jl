"""
    AbstractImpactMonitoringSource <: AbstractImpactMonitoring

Supertye for the impact monitoring sources interface.
"""
abstract type AbstractImpactMonitoringSource <: AbstractImpactMonitoring end

const SENTRY_JPL_API = "https://ssd-api.jpl.nasa.gov/sentry.api"

"""
    SENTRY <: AbstractImpactMonitoringSource

Type representing the Jet Propulsion Laboratory Sentry APIs.
"""
struct SENTRY <: AbstractImpactMonitoringSource end


function get_http_response(::Type{SENTRY}; id::Pair{String, String},
                           connect_timeout = 180, readtimeout = 180, kwargs...)
    query = isempty(first(id)) ? Dict{String, String}() : Dict{String, String}(id)
    resp = HTTP.get(SENTRY_JPL_API; query, connect_timeout, readtimeout,
        status_exception = false, require_ssl_verification = false, kwargs...)

    return resp
end

"""
    SentryVirtualImpactor{T} <: AbstractVirtualImpactor{T}


"""
struct SentryVirtualImpactor{T} <: AbstractVirtualImpactor{T}
    date::DateTime
    sigma::T
    ip::T
    ps::T
    ts::T
end

impactorsdefault(::Type{DateTime}) = DateTime(2000, 1, 1)
impactorsdefault(::Type{T}) where {T <: Real} = T(NaN)

sentryparse(_, ::Type{T}, x) where {T <: Real} = parse(Float64, x)
function sentryparse(_, ::Type{DateTime}, x)
    date = DateTime(view(x, 1:10), "yyyy-mm-dd")
    fraction = parse(Float64, view(x, 11:length(x)))
    return date + Microsecond(round(Int, fraction * 8.64e10))
end

function parse_impactors_sentry(text::AbstractString)
    # Parse JSON
    dict = JSON.parse(text)
    L = length(dict["data"])
    iszero(L) && return SentryVirtualImpactor{Float64}[]
    # Construct DataFrame
    R = SentryVirtualImpactor{Float64}
    names, types = fieldnames(R), fieldtypes(R)
    df = DataFrame([fill(impactorsdefault(fieldtype(R, name)), L) for name in names],
        collect(names))
    # Find the name of the sigma field
    snames = [string(name) for name in names]
    sigma_names = ["sigma_imp", "sigma_lov", "sigma_mc", "sigma_vi"]
    i = findfirst(Base.Fix1(haskey, first(dict["data"])), sigma_names)
    snames[2] = sigma_names[i]
    for (i, line) in enumerate(dict["data"])
        for (name, type, key) in zip(names, types, snames)
            x = line[key]
            df[i, name] = sentryparse(name, type, x)
        end
    end
    # Parse virtual impactors
    VIs = SentryVirtualImpactor.(eachrow(df))
    # Eliminate repeated entries
    unique!(VIs)
    # Sort by date
    sort!(VIs, by = date)

    return VIs
end

"""
    fetch_impactors_sentry(id, source)

Return the virtual impactors of minor body `id` in the SENTRY format.
The `source` of the observations can only be `SENTRY`.

!!! reference
    The JPL Sentry API is described at:
    - https://ssd-api.jpl.nasa.gov/doc/sentry.html
"""
function fetch_impactors_sentry(id::Pair{String, String}, ::Type{SENTRY})
    # Get and parse HTTP response
    text = fetch_http_text(SENTRY; id)
    # Parse observations
    VIs = parse_impactors_sentry(text)

    return VIs
end






using DataFrames: select!, transform!, ByRow
export SENTRY, impactor_table

function impactor_table(id::Pair{String, String}, ::Type{SENTRY})
    # Get and parse HTTP response
    text = fetch_http_text(SENTRY; id)
    # Parse JSON
    dict = JSON.parse(text)
    if !haskey(dict, "data")
        throw(ArgumentError(dict["error"]))
    end
    # Construct DataFrame
    df = DataFrame(dict["data"])
    # Reorder columns
    cols = [:date, :sigma, :ip, :energy, :ps, :ts]
    colnames = names(df)
    i = findfirst(startswith("sigma"), names(df))
    cols[2] = Symbol(colnames[i])
    select!(df, cols)
    # Sort rows
    sort!(df, :date)
    # Transform columns
    f = ByRow(sentryparse)
    transform!(df, cols[2:6] .=> f, renamecols = false)

    return df
end