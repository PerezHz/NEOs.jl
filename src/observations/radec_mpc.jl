abstract type AbstractAstrometry end

@doc raw"""
    RadecMPC{T <: AbstractFloat} <: AbstractAstrometry

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
- `catalogue::CatalogueMPC`: catalogue.
- `info2::String`: additional information.
- `observatory::ObservatoryMPC{T}`: observatory.
"""
@auto_hash_equals struct RadecMPC{T <: AbstractFloat} <: AbstractAstrometry
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
    catalogue::CatalogueMPC
    info2::String
    observatory::ObservatoryMPC{T}
    # Inner constructor
    function RadecMPC{T}(num::String, tmpdesig::String, discovery::String, publishnote::String,
                         obstech::String, date::DateTime, α::T, δ::T, info1::String, mag::String,
                         band::String, catalogue::CatalogueMPC, info2::String, observatory::ObservatoryMPC{T}) where {T <: AbstractFloat}
        new{T}(num, tmpdesig, discovery, publishnote, obstech, date, α, δ, info1, mag, band,
               catalogue, info2, observatory)
    end
end

# Outer constructor
function RadecMPC(num::String, tmpdesig::String, discovery::String, publishnote::String,
                  obstech::String, date::DateTime, α::T, δ::T, info1::String, mag::String,
                  band::String, catalogue::CatalogueMPC, info2::String, observatory::ObservatoryMPC{T}) where {T <: AbstractFloat}
    RadecMPC{T}(num, tmpdesig, discovery, publishnote, obstech, date, α, δ, info1, mag, band,
                catalogue, info2, observatory)
end

function RadecMPC(date::DateTime, α::T, δ::T, observatory::ObservatoryMPC{T}) where {T <: AbstractFloat}
    RadecMPC{T}("", "", "", "", "", date, α, δ, "", "", "", unknowncat(), "", observatory)
end

# Print method for RadecMPC
# Examples:
# N00hp15 α: 608995.65 δ: -25653.3 t: 2020-12-04T10:41:43.209 obs: 703
# 99942 α: 422475.3 δ: 97289.49 t: 2021-05-12T06:28:35.904 obs: F51
function show(io::IO, m::RadecMPC{T}) where {T <: AbstractFloat}
    # If there is no number, use temporary designation
    id_str = filter(!isspace, m.num) == "" ? m.tmpdesig : m.num

    print(io, id_str, " α: ", @sprintf("%.5f", rad2deg(m.α)), "° δ: ", @sprintf("%.5f", rad2deg(m.δ)), "° t: ", m.date,
              " obs: ", m.observatory.name)
end

# Functions to get specific fields of a RadecMPC object
num(r::RadecMPC{T}) where {T <: AbstractFloat} = r.num
tmpdesig(r::RadecMPC{T}) where {T <: AbstractFloat} = r.tmpdesig
discovery(r::RadecMPC{T}) where {T <: AbstractFloat} = r.discovery
publishnote(r::RadecMPC{T}) where {T <: AbstractFloat} = r.publishnote
obstech(r::RadecMPC{T}) where {T <: AbstractFloat} = r.obstech
ra(r::RadecMPC{T}) where {T <: AbstractFloat} = r.α
dec(r::RadecMPC{T}) where {T <: AbstractFloat} = r.δ
info1(r::RadecMPC{T}) where {T <: AbstractFloat} = r.info1
mag(r::RadecMPC{T}) where {T <: AbstractFloat} = r.mag
band(r::RadecMPC{T}) where {T <: AbstractFloat} = r.band
catalogue(r::RadecMPC{T}) where {T <: AbstractFloat} = r.catalogue
info2(r::RadecMPC{T}) where {T <: AbstractFloat} = r.info2
observatory(r::RadecMPC{T}) where {T <: AbstractFloat} = r.observatory

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
        # Catalogue regex (column 72)
        raw"(?P<catalogue>.{1})",
        # Info 2 regex (columns 73-77)
        raw"(?P<info2>.{5})",
        # Observatory code regex (columns 78-80)
        raw"(?P<obscode>.{3})"
    ]
))

@doc raw"""
    datetime(year::Int, month::Int, day::Int, utc::T) where {T <: Real}

Construct a `DateTime` type by parts. `utc` is the fraction of day.
"""
function datetime(year::Int, month::Int, day::Int, utc::T) where {T <: Real}
    return DateTime(year, month, day) + Microsecond( round(1e6*86_400*utc) )
end

@doc raw"""
    ra(hrs::Int, min::Int, sec::T) where {T <: Real}
    ra(obs::RadecMPC{T}) where {T <: AbstractFloat}

Return the right ascension in rad.
"""
function ra(hrs::Int, min::Int, sec::T) where {T <: Real}
    # Convert hours minutes seconds to deg
    α_deg =  15*(hrs + min/60 + sec/3_600)
    # Convert deg to rad
    α_rad = deg2rad(α_deg)

    return α_rad

end

@doc raw"""
    dec(sgn::String, deg::Int, min::Int, sec::T) where {T <: Real}
    dec(obs::RadecMPC{T}) where {T <: AbstractFloat}

Return the declination in rad.
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

@doc raw"""
    RadecMPC(m::RegexMatch)

Convert a match of `NEOs.mpc_radec_regex` to `RadecMPC`.
"""
function RadecMPC(m::RegexMatch)
    date = datetime(
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
    # Find catalogue in mpc_catalogues[] that matches catalogue
    catalogue = search_cat_code(string(m["catalogue"]))
    # Find observatory in mpc_observatories[] that matches obscode
    observatory = search_obs_code(string(m["obscode"]))

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
        catalogue,
        string(m["info2"]),
        observatory
    )
end

@doc raw"""
    read_radec_mpc(filename::String)

Return the matches of `NEOs.mpc_radec_regex` in `filename` as `RadecMPC`.
"""
function read_radec_mpc(filename::String)
    # Check that filename is a file
    @assert isfile(filename) "Cannot open file: $filename"
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
    parse_radec_mpc(f::F, text::String) where F

Return the matches of `NEOs.mpc_radec_regex` in `text`. A function `f(m::RegexMatch) -> Bool`
can be passed to filter the matches.
"""
function parse_radec_mpc(f::F, text::String) where F

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
    search_circulars_mpc(f::F, url1::String, url2::String; max_iter::Int = 10_000) where F

Iterate MPC circulars from `url1` to `url2` and return the matches of `NEOs.mpc_radec_regex`.
A function `f(m::RegexMatch) -> Bool` can be passed to filter the observations. If `url2` is not
reached before `max_iter` iterations, the function will print a warning and return the
matches found so far.
"""
function search_circulars_mpc(f::F, url1::String, url2::String; max_iter::Int = 10_000) where F

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

Return the date in MPC format.
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

Return the right ascension [rad] in MPC format.
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

Return the declination [rad] in MPC format.
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

Return an observation in MPC format.
"""
function mpc_radec_str(obs::RadecMPC{T}) where {T <: AbstractFloat}
    # Date string
    date_s = mpc_date_str(obs.date)
    # Right ascension string
    α_s = mpc_α_str(obs.α)
    # Declination string
    δ_s = mpc_δ_str(obs.δ)
    # Catalogue string
    catalogue_s = isunknown(obs.catalogue) ? " " : obs.catalogue.code
    # Observatory string
    obscode_s = isunknown(obs.observatory) ? "   " : obs.observatory.code
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
        catalogue_s,
        obs.info2,
        obscode_s,
        "\n"
    ])

    return obs_s
end

@doc raw"""
    write_radec_mpc(obs::Vector{RadecMPC{T}}, filename::String) where {T <: AbstractFloat}

Write `obs` to `filename` in MPC format.
"""
function write_radec_mpc(obs::Vector{RadecMPC{T}}, filename::String) where {T <: AbstractFloat}
    open(filename, "w") do file
        for i in eachindex(obs)
            line = mpc_radec_str(obs[i])
            write(file, line)
        end
    end
end

@doc raw"""
    get_radec_mpc(id::AbstractString, filename::AbstractString)

Download MPC optical astrometry of NEO `id` and save the output to `filename`. 
"""
function get_radec_mpc(id::AbstractString, filename::AbstractString)
    # MPC search url 
    search_url = search_mpc_url *  replace(id, " " => "+")
    # MPC observations file url 
    obs_url = obs_mpc_url * replace(id, " " => "_") * ".txt"
    # Download database search 
    download(search_url, filename)
    # Download observations file 
    download(obs_url, filename)
    return nothing
end 