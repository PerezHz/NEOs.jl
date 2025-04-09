@doc raw"""
    NEOCPObject{T <: Real}

An object listed in the Minor Planet Center NEO Confirmation Page [1].

## Fields
- `tmpdesig::String`: temporary designation.
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

!!! reference
    [1] https://www.minorplanetcenter.net/iau/NEO/toconfirm_tabular.html.
"""
@auto_hash_equals struct NEOCPObject{T <: Real}
    tmpdesig::String
    score::T
    date::DateTime
    α::T
    δ::T
    V::T
    updated::String
    nobs::Int
    arc::T
    H::T
    notseen::T
    # Inner constructor
    function NEOCPObject{T}(tmpdesig::String, score::T, date::DateTime, α::T,
                            δ::T, V::T, updated::String, nobs::Int, arc::T,
                            H::T, notseen::T) where {T <: Real}
        return new{T}(tmpdesig, score, date, α, δ, V, updated, nobs,
                      arc, H, notseen)
    end
end

function NEOCPObject(tmpdesig::String, score::T, date::DateTime, α::T,
                     δ::T, V::T, updated::String, nobs::Int, arc::T,
                     H::T, notseen::T) where {T <: Real}
    return NEOCPObject{T}(tmpdesig, score, date, α, δ, V, updated, nobs,
                          arc, H, notseen)
end

# Print method for NEOCPObject
# Examples:
# P11RPKq
# P21RPBC
function show(io::IO, x::NEOCPObject{T}) where {T <: Real}
    print(io, x.tmpdesig, " α: ", @sprintf("%.5f", rad2deg(x.α)), "° δ: ",
          @sprintf("%.5f", rad2deg(x.δ)), "° t: ", x.date)
end

@doc raw"""
    NEOCPObject(dict::Dict{String, Any})

Convert a dictionary to `NEOCPObject`.
"""
function NEOCPObject(dict::Dict{String, Any})
    # Temporary designation
    tmpdesig = String(dict["Temp_Desig"])
    # NEO desirability score
    score = Float64(dict["Score"])
    # Discovery date
    year = Int(dict["Discovery_year"])
    month = Int(dict["Discovery_month"])
    day = Float64(dict["Discovery_day"])
    day_i = clamp(floor(Int, day), 1, daysinmonth(Date(year, month)))
    utc = day - day_i
    date = DateTime(year, month, day_i) + Microsecond( round(1e6*daysec*utc) )
    # Right ascension [rad]
    α = deg2rad(Float64(dict["R.A."]) * 15)
    # Declination [rad]
    δ = deg2rad(Float64(dict["Decl."]))
    # Rough current magnitude
    V = Float64(dict["V"])
    # Last update
    updated = String(dict["Updated"])
    # Number of observations
    nobs = Int(dict["NObs"])
    # Arc length [days]
    arc = Float64(dict["Arc"])
    # Nominal absolute magnitude
    H = Float64(dict["H"])
    # Days since last observation
    notseen = Float64(dict["Not_Seen_dys"])

    return NEOCPObject(tmpdesig, score, date, α, δ, V, updated,
                       nobs, arc, H, notseen)
end

@doc raw"""
    fetch_objects_neocp()

Download the current Minor Planet Center NEO Confirmation Page [1]
list and return the output as `Vector{NEOCPObject{Float64}}`.

!!! reference
    [1] https://www.minorplanetcenter.net/iau/NEO/toconfirm_tabular.html.
"""
function fetch_objects_neocp()
    # HTTP response
    resp = HTTP.get(NEOCP_FILE_URL)
    # Parse response body as String
    text = String(resp.body)
    # Parse JSON
    dict = JSON.parse(text)
    # Parse NEOCPObject
    objs = NEOCPObject.(dict)

    return objs
end

@doc raw"""
    fetch_radec_neocp(id::AbstractString)

Download MPC optical astrometry of NEOCP object `id` and
return the output as `Vector{RadecMPC{Float64}}`.

!!! reference
    MPC NEOCP Observations API documentation:

    https://www.minorplanetcenter.net/mpcops/documentation/neocp-observations-api/.
"""
function fetch_radec_neocp(id::AbstractString)
    # HTTP parameters
    params = JSON.json(Dict("trksubs" => [id], "output_format" => ["OBS80"]))
    resp = HTTP.get(MPC_NEOCP_OBS_API_URL, ("Content-Type" => "application/json",), params)
    # Convert to String
    text = String(resp.body)
    # Parse JSON
    dict = JSON.parse(text)
    # Parse observations
    radec = RadecMPC.(eachmatch(RADEC_MPC_REGEX, dict[1]["OBS80"]))
    # Eliminate unsuccessful matches and repeated entries
    filter!(r -> isa(r, RadecMPC{Float64}), radec)
    unique!(radec)

    return radec
end

@doc raw"""
    get_radec_neocp(id::AbstractString [, filename::AbstractString])

Download MPC optical astrometry of NEOCP object `id` and save the output to `filename`
(default: `id.txt`).

!!! reference
    MPC NEOCP Observations API documentation:

    https://www.minorplanetcenter.net/mpcops/documentation/neocp-observations-api/.
"""
function get_radec_neocp(id::AbstractString, filename::AbstractString =
    replace(id, " " => "_") * ".txt")
    # Fetch optical astrometry
    radec = fetch_radec_neocp(id)
    # Write observations to file
    write_radec_mpc(radec, filename)

    return filename
end

@doc raw"""
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