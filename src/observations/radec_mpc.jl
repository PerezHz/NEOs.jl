@doc raw"""
    RadecMPC{T <: AbstractFloat} 

An optical (α, δ) measurement in MPC format. The format is described in https://minorplanetcenter.net/iau/info/OpticalObs.html
and discussed thoroughly in pages 158-181 of https://doi.org/10.1016/j.icarus.2010.06.003.

# Fields 

- `num::String`: object's number. 
- `tmpdesig::String`: provisional / temporary designation.
- `discovery::String`: discovery asterisk.
- `publishnote::String`: publishable note.
- `obstech::String`: observation technique.
- `date::DateTime`: date of observation.
- `α::T`: right ascension [rad].
- `δ::T`: declination [rad].
- `info1::String`: additional information. 
- `mag::String`: observed magnitude. 
- `band::String`: magnitude band. 
- `catalog::String`: catalog. 
- `info2::String`: additional information. 
- `observatory::ObservatoryMPC{T}`: observatory. 
"""
struct RadecMPC{T <: AbstractFloat} 
    num::String 
    tmpdesig::String
    discovery::String
    publishnote::String
    obstech::String
    date::DateTime
    α::T
    δ::T
    info1::String
    mag::String
    band::String
    catalog::String
    info2::String
    observatory::ObservatoryMPC{T}
    # Inner constructor 
    function RadecMPC{T}(num::String, tmpdesig::String, discovery::String, publishnote::String, 
                         obstech::String, date::DateTime, α::T, δ::T, info1::String, mag::String, 
                         band::String, catalog::String, info2::String, observatory::ObservatoryMPC{T}) where {T <: AbstractFloat}
        new{T}(num, tmpdesig, discovery, publishnote, obstech, date, α, δ, info1, mag, band, 
               catalog, info2, observatory)
    end
end

# Outer constructor
function RadecMPC(num::String, tmpdesig::String, discovery::String, publishnote::String, 
                  obstech::String, date::DateTime, α::T, δ::T, info1::String, mag::String, 
                  band::String, catalog::String, info2::String, observatory::ObservatoryMPC{T}) where {T <: AbstractFloat} 
    RadecMPC{T}(num, tmpdesig, discovery, publishnote, obstech, date, α, δ, info1, mag, band,
                catalog, info2, observatory)
end

# Two RadecMPC are equal if ther date, α, δ and observatory are equal
function hash(a::RadecMPC{T}, h::UInt) where {T <: AbstractFloat}
    return hash((a.date, a.α, a.δ, a.observatory), h)
end

function ==(a::RadecMPC{T}, b::RadecMPC{T}) where {T <: AbstractFloat}
    return hash(a) == hash(b)
end

# Print method for RadecMPC
# Examples: 
# N00hp15 α: 608995.65 δ: -25653.3 t: 2020-12-04T10:41:43.209 obs: 703
# 99942 α: 422475.3 δ: 97289.49 t: 2021-05-12T06:28:35.904 obs: F51
function show(io::IO, m::RadecMPC{T}) where {T <: AbstractFloat} 
    # If there is no number, use temporary designation
    id_str = filter(!isspace, m.num) == "" ? m.tmpdesig : m.num

    print(io, id_str, " α: ", m.α, " δ: ", m.δ, " t: ", m.date,
              " obs: ", m.observatory.name)
end

# Order in RadecMPC is given by date 
isless(a::RadecMPC{T}, b::RadecMPC{T}) where {T <: AbstractFloat} = a.date < b.date

# Regular expression to parse an optical measurement in MPC format
const mpc_radec_regex = Regex(join(
    [
        # Number regex (columns 1-5)
        raw"(?P<num>.{5})",
        # Temporary designation regex (columns 6-12)
        raw"(?P<tmpdesig>.{7})",
        # Discovery asterisk regex (column 13)
        raw"(?P<discovery>.{1})",
        # Publishable note regex (column 14)
        raw"(?P<publishnote>.{1})",
        # Observation technique regex (column 15)
        raw"(?P<obstech>.{1})",
        # Year regex + space (columns 16-20)
        raw"(?P<year>\d{4}) ",
        # Month regex + space (columns 21-23)
        raw"(?P<month>\d{2}) ",
        # Day regex (columns 24-25)
        raw"(?P<day>\d{2})",
        # Fraction of days regex (columns 26-32) 
        raw"(?P<utc>\.[\d\s]{6})",
        # α hours regex + space (columns 33-35)
        raw"(?P<α_hrs>\d{2}) ",
        # α minutes regex + space (columns 36-38)
        raw"(?P<α_min>\d{2}) ",
        # α seconds regex (columns 39-44)
        raw"(?P<α_sec>\d{2}\.[\d\s]{3})",
        # δ sign regex (column 45)
        raw"(?P<δ_sgn>\+|\-)",
        # δ degrees regex + space (columns 46-48)
        raw"(?P<δ_deg>\d{2}) ",
        # δ minutes regex + space (columns 49-51)
        raw"(?P<δ_min>\d{2}) ",
        # δ seconds regex (columns 52-56)
        raw"(?P<δ_sec>\d{2}\.[\d\s]{2})",
        # Info 1 regex (columns 57-65)
        raw"(?P<info1>.{9})",
        # Magnitude regex (columns 66-70)
        raw"(?P<mag>.{5})",
        # Band regex (column 71)
        raw"(?P<band>.{1})",
        # Catalog regex (column 72)
        raw"(?P<catalog>.{1})",
        # Info 2 regex (columns 73-77)
        raw"(?P<info2>.{5})",
        # Observatory code regex (columns 78-80)
        raw"(?P<obscode>.{3})"
    ]
))

@doc raw"""
    DateTime(year::Int, month::Int, day::Int, utc::T) where {T <: Real}

Construct a `DateTime` type by parts. `utc` is the fraction of day. 
"""
function DateTime(year::Int, month::Int, day::Int, utc::T) where {T <: Real}
    return DateTime(year, month, day) + Microsecond( round(1e6*86_400*utc) )
end

@doc raw"""
    ra(hrs::Int, min::Int, sec::T) where {T <: Real}
    ra(obs::RadecMPC{T}) where {T <: AbstractFloat}

Returns the right ascension in rad. 
"""
function ra(hrs::Int, min::Int, sec::T) where {T <: Real}
    # Convert hours minutes seconds to deg
    α_deg =  15*(hrs + min/60 + sec/3_600)
    # Convert deg to rad
    α_rad = deg2rad(α_deg)
    
    return α_rad

end

ra(obs::RadecMPC{T}) where {T <: AbstractFloat} = getfield(obs, :α)

@doc raw"""
    dec(sgn::String, deg::Int, min::Int, sec::T) where {T <: Real}
    dec(obs::RadecMPC{T}) where {T <: AbstractFloat}

Returns the declination in rad. 
"""
function dec(sgn::String, deg::Int, min::Int, sec::T) where {T <: Real}
    # Convert degrees minutes seconds to degrees
    if sgn == "+"
        δ_deg = +(deg + min/60 + sec/3_600)
    elseif sgn == "-"
        δ_deg = -(deg + min/60 + sec/3_600)
    end

    δ_rad = deg2rad(δ_deg)

    return δ_rad
end

dec(obs::RadecMPC{T}) where {T <: AbstractFloat} = getfield(obs, :δ)

@doc raw"""
    RadecMPC(m::RegexMatch)

Converts a match of `NEOs.mpc_radec_regex` to `RadecMPC`.
"""
function RadecMPC(m::RegexMatch)
    date = DateTime(
        Meta.parse(m["year"]), 
        Meta.parse(m["month"]), 
        Meta.parse(m["day"]),
        Meta.parse(m["utc"])
    )
    α = ra(
        Meta.parse(m["α_hrs"]), 
        Meta.parse(m["α_min"]), 
        Meta.parse(m["α_sec"])
    )
    δ = dec(
        string(m["δ_sgn"]), 
        Meta.parse(m["δ_deg"]), 
        Meta.parse(m["δ_min"]), 
        Meta.parse(m["δ_sec"])
    )
    i = findfirst(x -> x.code == string(m["obscode"]), mpc_observatories[])
    if isnothing(i)
        observatory = unknownobs()
    else
        observatory = mpc_observatories[][i]
    end
    
    return RadecMPC(
        string(m["num"]),
        string(m["tmpdesig"]),
        string(m["discovery"]),
        string(m["publishnote"]),
        string(m["obstech"]),
        date, 
        α, 
        δ, 
        string(m["info1"]),
        string(m["mag"]),
        string(m["band"]),
        string(m["catalog"]),
        string(m["info2"]),
        observatory
    )
end

@doc raw"""
    read_radec_mpc(filename::String)

Returns the matches of `NEOs.mpc_radec_regex` in `filename` as `RadecMPC`.
"""
function read_radec_mpc(filename::String)
    # Read lines of mpc formatted file 
    lines = readlines(filename)
    # Apply regular expressions
    matches = match.(mpc_radec_regex, lines)
    # Eliminate nothings
    filter!(!isnothing, matches)
    # Convert matches to RadecMPC
    obs = RadecMPC.(matches)
    # Sort observations by date
    sort!(obs)
    # Eliminate repeated observations
    unique!(obs)
    
    return obs
end

@doc raw"""
    parse_radec_mpc(text::String)
    parse_radec_mpc(f::Function, text::String)

Returns the matches of `NEOs.mpc_radec_regex` in `text`. A function `f(m::RegexMatch) -> Bool` 
can be passed to filter the matches. 
"""
function parse_radec_mpc(f::Function, text::String)

    # Vector of observations 
    radecs = Vector{RadecMPC{Float64}}(undef, 0)
    # Iterate over the matches 
    for m in eachmatch(mpc_radec_regex, text)
        # Filter by f
        if f(m)
            push!(radecs, RadecMPC(m))
        end
    end
    # If there is at least one observation
    if length(radecs) > 0
        # Sort observations by date
        sort!(radecs)
        # Eliminate repeated observations
        unique!(radecs)
    end
    
    return radecs
end
parse_radec_mpc(text::String) = parse_radec_mpc(t -> true, text)

# MPC main page url 
const mpc_url = "https://minorplanetcenter.net"

# Regex for next circular url 
const next_circular_regex = r"<a href=\"(?P<next>.*)\"><img src=\"/iau/figs/RArrow.gif\""

@doc raw"""
    search_circulars_mpc(url1::String, url2::String; max_iter::Int = 10_000)
    search_circulars_mpc(f::Function, url1::String, url2::String; max_iter::Int = 10_000)

Iterates MPC circulars from `url1` to `url2` and returns the matches of `NEOs.mpc_radec_regex`. 
A function `f(m::RegexMatch) -> Bool` can be passed to filter the observations. If `url2` is not
reached before `max_iter` iterations, the function will print a warning and return the 
matches found so far. 
"""
function search_circulars_mpc(f::Function, url1::String, url2::String; max_iter::Int = 10_000)

    # Vector of observations
    obs = Vector{RadecMPC{Float64}}(undef, 0)
    
    # Number of urls checked 
    n = 0
    # First url 
    u = url1

    while true 
        n += 1
        if n > max_iter
            @warn("$n pages checked before getting to $url2")
            break 
        end
        # Raw html text of webpage u 
        text = get_raw_html(u)
        # Observations found in text 
        obs_ = parse_radec_mpc(f, text)
        # Add new observations 
        obs = vcat(obs, obs_) 
        # Final url 
        if u == url2
            break
        end
        # Next circular url 
        next = match(next_circular_regex, text)["next"]
        u = mpc_url * next
    end
    # If there is at least one observation
    if length(obs) > 0
        # Sort observations by date 
        sort!(obs)
        # Eliminate repeated observations
        unique!(obs)
    end

    return obs
end

search_circulars_mpc(url1::String, url2::String; max_iter::Int = 10_000) = search_circulars_mpc(t -> true, url1, url2; max_iter = max_iter)

@doc raw"""
    mpc_date_str(date::DateTime)

Returns the date in MPC format. 
"""
function mpc_date_str(date::DateTime)
    
    # Year string 
    year_s = lpad(Dates.year(date), 4)
    # Month string 
    month_s = lpad(Dates.month(date), 2, "0")
    # Hours [days] 
    hrs = Dates.hour(date) / 24
    # Minutes [days] 
    min = Dates.minute(date) / 24 / 60
    # Seconds [days] 
    sec = Dates.second(date) / 24 / 60 / 60
    # Milliseconds [days] 
    mls = Dates.millisecond(date) / 24 / 60 / 60 / 1_000
    # Days 
    day_val = Dates.day(date) + hrs + min + sec + mls
    # Days string 
    day_s = @sprintf("%09.6f", day_val)
    # Join everything 
    date_s = join([
        year_s,
        " ",
        month_s,
        " ",
        day_s,
    ])

    return date_s
end

@doc raw"""
    mpc_α_str(α::T) where {T <: Number}

Returns the right ascension [rad] in MPC format. 
"""
function mpc_α_str(α::T) where {T <: Number}
    # Convert rad to deg 
    α_deg = rad2deg(α)
    # Hours 
    hrs_, hrs = modf(α_deg / 15)
    hrs = Int(hrs)
    # Hours string 
    hrs_s = lpad(hrs, 2, "0")
    # Minutes 
    min_, min = modf(60 * hrs_)
    min = Int(min)
    # Minutes string 
    min_s = lpad(min, 2, "0")
    # Seconds 
    sec = 60 * min_
    # Seconds string 
    sec_s = @sprintf("%06.3f", sec)
    # Join everything
    α_s = join([
        hrs_s,
        " ",
        min_s,
        " ",
        sec_s,
    ])

    return α_s 
end

@doc raw"""
    mpc_δ_str(δ::T) where {T <: Number}

Returns the declination [rad] in MPC format. 
"""
function mpc_δ_str(δ::T) where {T <: Number}
    # Sign string 
    sgn_s = δ >= 0 ? "+" : "-"
    # Convert rad to deg 
    δ_deg = abs(rad2deg(δ))
    # Degrees 
    deg_, deg = modf(δ_deg)
    deg = Int(deg)
    # Degrees string 
    deg_s = lpad(deg, 2, "0")
    # Minutes 
    min_, min = modf(60 * deg_)
    min = Int(min)
    # Minutes string 
    min_s = lpad(min, 2, "0")
    # Seconds
    sec = 60 * min_
    # Seconds string 
    sec_s = @sprintf("%05.2f", sec)
    # Join everything
    δ_s = join([
        sgn_s,
        deg_s,
        " ",
        min_s,
        " ",
        sec_s,
    ])

    return δ_s 
end

@doc raw"""
    mpc_radec_str(obs::RadecMPC{T}) where {T <: AbstractFloat}

Returns an observation in MPC format. 
"""
function mpc_radec_str(obs::RadecMPC{T}) where {T <: AbstractFloat}
    # Date string 
    date_s = mpc_date_str(obs.date)
    # Right ascension string 
    α_s = mpc_α_str(obs.α)
    # Declination string 
    δ_s = mpc_δ_str(obs.δ)
    # Join everything
    obs_s = join([
        obs.num,
        obs.tmpdesig,
        obs.discovery,
        obs.publishnote,
        obs.obstech,
        date_s,
        α_s,
        δ_s,
        obs.info1,
        obs.mag,
        obs.band,
        obs.catalog,
        obs.info2,
        obs.observatory.code,
        "\n"
    ])

    return obs_s
end

@doc raw"""
    write_radec_mpc(obs::Vector{RadecMPC{T}}, filename::String) where {T <: AbstractFloat}

Writes `obs` to `filename` in MPC format. 
"""
function write_radec_mpc(obs::Vector{RadecMPC{T}}, filename::String) where {T <: AbstractFloat}
    open(filename, "w") do file
        for i in eachindex(obs)
            line = mpc_radec_str(obs[i])
            write(file, line)
        end 
    end
end