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
- `obscode::String`: observatory code. 
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
    obscode::String
    # Inner constructor 
    function RadecMPC{T}(num::String, tmpdesig::String, discovery::String, publishnote::String, 
                         obstech::String, date::DateTime, α::T, δ::T, info1::String, mag::String, 
                         band::String, catalog::String, info2::String, obscode::String) where {T <: AbstractFloat}
        new{T}(num, tmpdesig, discovery, publishnote, obstech, date, α, δ, info1, mag, band, 
               catalog, info2, obscode)
    end
end

# Outer constructor
function RadecMPC(num::String, tmpdesig::String, discovery::String, publishnote::String, 
                  obstech::String, date::DateTime, α::T, δ::T, info1::String, mag::String, 
                  band::String, catalog::String, info2::String, obscode::String) where {T <: AbstractFloat} 
    RadecMPC{T}(num, tmpdesig, discovery, publishnote, obstech, date, α, δ, info1, mag, band,
                catalog, info2, obscode)
end

# Two RadecMPC are equal if ther date, α, δ and observatory code are equal
function hash(a::RadecMPC{T}, h::UInt) where {T <: AbstractFloat}
    return hash((a.date, a.α, a.δ, a.obscode), h)
end

function ==(a::RadecMPC{T}, b::RadecMPC{T}) where {T <: AbstractFloat}
    return hash(a) == hash(b)
end

# Print method for RadecMPC
# Examples: 
# id: N00hp15 α: 608995.65 δ: -25653.3 t: 2020-12-04T10:41:43.209 obs: 703
# id: 99942 α: 422475.3 δ: 97289.49 t: 2021-05-12T06:28:35.904 obs: F51
function show(io::IO, m::RadecMPC{T}) where {T <: AbstractFloat} 
    # If there is no number, use temporary designation
    id_str = filter(!isspace, m.num) == "" ? m.tmpdesig : m.num

    print(io, "id: ", id_str, " α: ", m.α, " δ: ", m.δ, " t: ", m.date,
              " obs: ", m.obscode)
end

sort(obs::Vector{RadecMPC{T}}) where {T <: AbstractFloat} = sort(obs, by = x -> x.date)
sort!(obs::Vector{RadecMPC{T}}) where {T <: AbstractFloat} = sort!(obs, by = x -> x.date)

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
        string(m["obscode"])
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
    get_raw_html(url::String = "https://minorplanetcenter.net/mpec/K20/K20YA9.html")

Returns the raw html text of webpage `url`.
"""
function get_raw_html(url::String = "https://minorplanetcenter.net/mpec/K20/K20YA9.html")
    # Get raw html 
    resp = get(url)
    # Convert to string 
    text = String(resp.body)
    
    return text
end

