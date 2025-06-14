"""
    NEOCPObject{T} <: AbstractOpticalAstrometry{T}

An object listed in the Minor Planet Center NEO Confirmation Page.

# Fields

- `desig::String`: temporary designation.
- `score::T`: NEO desirability score.
- `date::DateTime`: discovery date.
- `α::T`: right ascension [rad].
- `δ::T`: declination [rad].
- `V::T`: rough current magnitude.
- `updated::String`: last update.
- `nobs::Int`: number of observations.
- `arc::T`: arc length [days].
- `H::T`: nominal absolute magnitude.
- `notseen::T`: days since last observation.
- `source::String`: original record.

!!! reference
    The Minor Planet Center Confirmation Page is available at:
    - https://www.minorplanetcenter.net/iau/NEO/toconfirm_tabular.html
"""
@auto_hash_equals fields = (date, ra, dec) struct NEOCPObject{T} <: AbstractOpticalAstrometry{T}
    desig::String
    score::T
    date::DateTime
    ra::T
    dec::T
    V::T
    updated::String
    nobs::Int
    arc::T
    H::T
    notseen::T
    source::String
end

date(x::NEOCPObject) = x.date
observatory(::NEOCPObject{T}) where {T} = unknownobs(T)

# Print method for NEOCPObject
function show(io::IO, x::NEOCPObject)
    print(io, x.desig, " α: ", @sprintf("%.5f", rad2deg(x.ra)), "° δ: ",
        @sprintf("%.5f", rad2deg(x.dec)), "° t: ", x.date)
end

# Parsing

neocpparse(_, ::Type{String}, x) = String(x)

neocpparse(_, ::Type{Int}, x) = parse(Int, x)

function neocpparse(_, ::Type{DateTime}, x)
    date = DateTime(view(x, 1:10), "yyyy mm dd")
    fraction = parse(Float64, view(x, 11:length(x)))
    return date + Microsecond(round(Int, fraction * 8.64e10))
end

function neocpparse(name, ::Type{T}, x) where {T <: Real}
    y = isempty(x) ? T(NaN) : parse(T, x)
    if name == :ra
        y = deg2rad(y * 15)
    elseif name == :dec
        y = deg2rad(y)
    end
    return y
end

NEOCPObject(r::DataFrameRow) = NEOCPObject{Float64}(
    r.desig, r.score, r.date, r.ra, r.dec, r.V, r.updated,
    r.nobs, r.arc, r.H, r.notseen, r.source
)

function parse_neocp_objects(text::AbstractString)
    # Parse lines
    lines = split(text, '\n', keepempty = false)
    L = length(lines)
    # Construct DataFrame
    R = NEOCPObject{Float64}
    names, types = fieldnames(R), fieldtypes(R)
    df = DataFrame([fill(astrometrydefault(fieldtype(R, name)), L) for name in names],
        collect(names))
    for (i, line) in enumerate(lines)
        for (name, type, idxs) in zip(names, types, NEOCP_OBJECT_COLUMNS)
            x = strip(view(line, idxs))
            df[i, name] = neocpparse(name, type, x)
        end
    end
    # Source string
    df.source = lines
    # Parse objects
    objects = NEOCPObject.(eachrow(df))
    # Eliminate repeated entries
    unique!(objects)
    # Sort by date
    sort!(objects)

    return objects
end

"""
    fetch_neocp_objects()

Return the current objects listed in the Minor Planet Center
NEO Confirmation Page.

!!! reference
    The Minor Planet Center Confirmation Page is available at:
    - https://www.minorplanetcenter.net/iau/NEO/toconfirm_tabular.html
"""
function fetch_neocp_objects()
    # Get and parse HTTP response
    text = fetch_http_text(NEOCP; mode = 0)
    # Parse objects
    objects = parse_neocp_objects(text)

    return objects
end

# Read / write

"""
    read_neocp_objects(filename)

Read from `filename` a vector of objects listed in the Minor
Planet Center NEO Confirmation Page.
"""
function read_neocp_objects(filename::AbstractString)
    # Read file
    text = read(filename, String)
    # Parse observations
    objects = parse_neocp_objects(text)

    return objects
end

"""
    write_neocp_objects(objs, filename)

Write a vector of objects listed in the Minor Planet Center
NEO Confirmation Page `objs` to `filename`.
"""
function write_neocp_objects(obs::AbstractVector{NEOCPObject{T}},
                             filename::AbstractString) where {T <: Real}
    open(filename, "w") do file
        for i in eachindex(obs)
            write(file, obs[i].source, '\n')
        end
    end
end

#=
"""
    get_orbits_neocp(id::AbstractString [, filename::AbstractString])

Download sample orbits of NEOCP [1] object `id` and save the output to `filename`
(default: `id.txt`).

!!! reference
    [1] https://www.minorplanetcenter.net/iau/NEO/toconfirm_tabular.html.
"""
function get_orbits_neocp(id::AbstractString, filename::AbstractString =
    replace(id, " " => "_") * ".txt")
    # Assemble URL
    url = string(NEOCP_SHOWORBS_URL, "?Obj=", id, "&orb=y")
    # HTTP response
    resp = HTTP.get(url)
    # Parse response body as String
    text = String(resp.body)
    # Break lines
    lines = split(text, "\n", keepempty = false)
    # Save result to filename
    open(filename, "w") do file
        for i in 2:length(lines)-1
            write(file, lines[i], "\n")
        end
    end

    return filename
end
=#