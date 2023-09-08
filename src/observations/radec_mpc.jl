abstract type AbstractAstrometry end
struct RightAscension <: AbstractAstrometry end
struct Declination <: AbstractAstrometry end

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
- `mag::T`: observed magnitude.
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
    mag::T
    band::String
    catalogue::CatalogueMPC
    info2::String
    observatory::ObservatoryMPC{T}
    # Inner constructor
    function RadecMPC{T}(num::String, tmpdesig::String, discovery::String, publishnote::String,
                         obstech::String, date::DateTime, α::T, δ::T, info1::String, mag::T,
                         band::String, catalogue::CatalogueMPC, info2::String, observatory::ObservatoryMPC{T}) where {T <: AbstractFloat}
        new{T}(num, tmpdesig, discovery, publishnote, obstech, date, α, δ, info1, mag, band,
               catalogue, info2, observatory)
    end
end

# Outer constructor
function RadecMPC(num::String, tmpdesig::String, discovery::String, publishnote::String,
                  obstech::String, date::DateTime, α::T, δ::T, info1::String, mag::T,
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

function neoparse(x::RegexMatch, i::Int, ::Type{Int64})
    y = tryparse(Int64, x[i])
    if isnothing(y)
        return 0
    else
        return y
    end
end

function neoparse(x::RegexMatch, i::Int, ::Type{DateTime})
    date = DateTime(x[i][1:10], "yyyy mm dd")
    utc = parse(Float64, x[i][11:end])
    return date + Microsecond( round(1e6*86_400*utc) )
end

function neoparse(x::RegexMatch, i::Int, ::Type{RightAscension})
    # Unfold 
    hrs = parse(Int, x[i][1:2])
    min = parse(Int, x[i][4:5])
    sec = parse(Float64, x[i][7:end])
    # Convert hours minutes seconds to deg
    α_deg =  15*(hrs + min/60 + sec/3_600)
    # Convert deg to rad
    α_rad = deg2rad(α_deg)

    return α_rad

end

function neoparse(x::RegexMatch, i::Int, ::Type{Declination})
    # Unfold 
    sgn = Char(x[i][1])
    deg = parse(Int, x[i][2:3])
    min = parse(Int, x[i][5:6])
    sec = parse(Float64, x[i][8:end])
    # Convert degrees minutes seconds to degrees
    if sgn == '+'
        δ_deg = +(deg + min/60 + sec/3_600)
    elseif sgn == '-'
        δ_deg = -(deg + min/60 + sec/3_600)
    end

    δ_rad = deg2rad(δ_deg)

    return δ_rad
end

neoparse(x::RegexMatch, i::Int, ::Type{CatalogueMPC}) = search_cat_code(String(x[i]))

function neoparse(m::RegexMatch, i::Int, ::Type{ObservatoryMPC{Float64}})
    _observatory_ = search_obs_code(String(m[i]))
    if isnothing(m["optional"])
        return _observatory_
    else 
        units = neoparse(m, 22, Int64)
        x = neoparse(m, 23, Float64)
        y = neoparse(m, 24, Float64)
        z = neoparse(m, 25, Float64)
        date = neoparse(m, 6, DateTime)
        if units == 2
            x *= au
            y *= au
            z *= au
        end
        return ObservatoryMPC(_observatory_.code, x, y, z, _observatory_.name, date, :satellite, units)
    end
end

# Regular expression to parse an optical measurement in MPC format
const RADEC_MPC_REGEX = Regex(string(
    # Number regex (columns 1-5)
    raw"(?P<num>.{5})",
    # Temporary designation regex (columns 6-12)
    raw"(?P<tmpdesig>.{7})",
    # Discovery asterisk regex (column 13)
    raw"(?P<discovery>[\*\s]{1})",
    # Publishable note regex (column 14)
    raw"(?P<publishnote>.{1})",
    # Observation technique regex (column 15)
    raw"(?P<obstech>.{1})",
    # Date of observation regex (columns 16-32)
    raw"(?P<date>\d{4}\s\d{2}\s\d{2}\.[\d\s]{6})",
    # Right ascension regex (columns 33-44)
    raw"(?P<α>\d{2}\s\d{2}\s\d{2}\.[\d\s]{3})",
    # Declination regex (columns 45-56)
    raw"(?P<δ>[\+|\-]{1}\d{2}\s\d{2}\s\d{2}\.[\d\s]{2})",
    # Info 1 regex (columns 57-65)
    raw"(?P<info1>.{9})",
    # Magnitude regex (columns 66-70)
    raw"(?P<mag>[\.\s\d]{5})",
    # Band regex (column 71)
    raw"(?P<band>[\w\s]{1})",
    # Catalogue regex (column 72)
    raw"(?P<catalogue>[\w\s]{1})",
    # Info 2 regex (columns 73-77)
    raw"(?P<info2>.{5})",
    # Observatory code regex (columns 78-80)
    raw"(?P<obscode>\w{3})",
    # Optional fields (in case of satellite observations)
    # Breakline regex
    raw"(?:\n)?",
    # Number regex (columns 1-5)
    raw"(?<optional>(?P<_num_>.{5})?",
    # Temporary designation regex (columns 6-12)
    raw"(?P<_tmpdesig_>.{7})?",
    # Blank space regex (column 13)
    raw"(?P<_discovery_>\s)?",
    # Publishable note regex (column 14)
    raw"(?P<_publishnote_>.{1})?",
    # s regex (column 15)
    raw"(?P<_obstech_>s)?",
    # Date of observation regex (columns 16-32)
    raw"(?P<_date_>\d{4}\s\d{2}\s\d{2}\.[\d\s]{6})",
    # Units + space regex (columns 33-34)
    raw"(?P<_units_>\d\s)",
    # X component of geocentric vector (columns 35-46)
    raw"(?P<_x_>[\-\+]{1}[\.\d\s]{11})",
    # Y component of geocentric vector (columns 47-58)
    raw"(?P<_y_>[\-\+]{1}[\.\d\s]{11})",
    # Z component of geocentric vector (columns 59-70)
    raw"(?P<_z_>[\-\+]{1}[\.\d\s]{11})",
    # Band regex (column 71)
    raw"(?P<_band_>[\w\s]{1})?",
    # Catalogue regex (column 72)
    raw"(?P<_catalogue_>[\w\s]{1})?",
    # Info 2 regex (columns 73-77)
    raw"(?P<_info2_>.{5})?",
    # Observatory code regex (columns 78-80)
    raw"(?P<_obscode_>\w{3})?)?",
))

@doc raw"""
    RadecMPC(m::RegexMatch)

Convert a match of `NEOs.RADEC_MPC_REGEX` to `RadecMPC`.
"""
function RadecMPC(m::RegexMatch)
    # Check that matched regex is correct
    @assert m.regex == RADEC_MPC_REGEX "Only matches of `NEOs.RADEC_MPC_REGEX` can be converted to `RadecMPC`."
    # Field types
    types = fieldtypes(RadecMPC{Float64})
    # RadecMPC{Float64} fields
    args = map(i -> begin
        if i == 7
            neoparse(m, i, RightAscension)
        elseif i == 8
            neoparse(m, i, Declination)
        else 
            neoparse(m, i, types[i])
        end
    end, 1:length(types))

    return RadecMPC(args...)
end

@doc raw"""
    read_radec_mpc(s::String)

Return the matches of `NEOs.RADEC_MPC_REGEX` in `s` as `RadecMPC`.
"""
function read_radec_mpc(s::String)
    if !contains(s, "\n") && isfile(s)
        # Read MPC formatted file
        s = String(read(s))
    end
    # Vector of observations
    radec = Vector{RadecMPC{Float64}}(undef, 0)
    # Iterate over the matches
    for m in eachmatch(RADEC_MPC_REGEX, s)
        push!(radec, RadecMPC(m))
    end
    # Eliminate repeated entries
    unique!(radec)
    # Sort observations by date
    sort!(radec)

    return radec
end

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
    date_s = string(year_s, " ", month_s, " ", day_s)

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
    α_s = string(hrs_s, " ", min_s, " ", sec_s)

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
    δ_s = string(sgn_s, deg_s, " ", min_s, " ", sec_s)

    return δ_s
end

function mpc_x_str(x::T) where {T <: AbstractFloat}
    sgn = x > 0 ? "+" : "-"
    y = @sprintf("%.5f", abs(x))
    y = lpad(y, 11)
    return string(sgn, y)
end

# Convert `obs` to a string according to MPC format.
function string(obs::RadecMPC{T}) where {T <: AbstractFloat}
    # Number string 
    num_s = rpad(obs.num, 5)
    # Temporary designation string
    tmpdesig_s = rpad(obs.tmpdesig, 7)
    # Discovery asterisk string
    discovery_s = rpad(obs.discovery, 1)
    # Publishable note string
    publishnote_s = rpad(obs.publishnote, 1)
    # Observation technique string
    obstech_s = rpad(obs.obstech, 1)
    # Date string
    date_s = mpc_date_str(obs.date)
    # Right ascension string
    α_s = mpc_α_str(obs.α)
    # Declination string
    δ_s = mpc_δ_str(obs.δ)
    # Info 1 string
    info1_s = rpad(obs.info1, 9)    
    # Magnitude string
    mag_s = isnan(obs.mag) ? repeat(" ", 5) : @sprintf("%.2f", obs.mag)
    # Band string 
    band_s = rpad(obs.band, 1)
    # Info 2 string 
    info2_s = rpad(obs.info2, 5)    
    # Catalogue string
    catalogue_s = isunknown(obs.catalogue) ? " " : obs.catalogue.code
    # Observatory string
    obscode_s = isunknown(obs.observatory) ? "   " : obs.observatory.code
    # Join everything
    obs_s = string(num_s, tmpdesig_s, discovery_s, publishnote_s, obstech_s, date_s, α_s, δ_s,
                   info1_s, mag_s, band_s, catalogue_s, info2_s, obscode_s)

    if issatellite(obs.observatory)
        # Units string 
        units = obs.observatory.units
        units_s = rpad(units, 1)
        
        x = obs.observatory.long
        y = obs.observatory.cos
        z = obs.observatory.sin
        if units == 2
            x /= au
            y /= au
            z /= au
        end

        # X component string 
        x_s = mpc_x_str(x)
        # Y component string 
        y_s = mpc_x_str(y)
        # Z component string 
        z_s = mpc_x_str(z)
        obs_s = string(obs_s, "\n", num_s, tmpdesig_s, " ", publishnote_s, "s", date_s,
                       units_s, " ", x_s, y_s, z_s, "  ", info2_s, obscode_s)
    end

    return obs_s
end

@doc raw"""
    write_radec_mpc(obs::Vector{RadecMPC{T}}, filename::String) where {T <: AbstractFloat}

Write `obs` to `filename` in MPC format.
"""
function write_radec_mpc(obs::Vector{RadecMPC{T}}, filename::String) where {T <: AbstractFloat}
    open(filename, "w") do file
        for i in eachindex(obs)
            line = string(obs[i])
            write(file, line, "\n")
        end
    end
end

@doc raw"""
    get_radec_mpc(id::Pair{String, String}, filename::String = replace(id[2], " " => "_") * ".txt")

Download MPC optical astrometry of NEO `id` and save the output to `filename`. 
"""
function get_radec_mpc(id::Pair{String, String}, filename::String = replace(id[2], " " => "_") * ".txt")
    # HTTP query
    query = ["table" => "observations", id]
    resp = get("http://minorplanetcenter.net/search_db"; query = query)
    # Converty to String
    text = String(resp.body)
    # Parse JSON
    obs = JSON.parse(text)
    # Find matches
    matches = Vector{RegexMatch}(undef, length(obs))
    for i in eachindex(obs)
        s = obs[i]["original_record"]
        matches[i] = match(RADEC_MPC_REGEX, s)
    end
    filter!(!isnothing, matches)
    # Parse RadecMPC
    radec = RadecMPC.(matches)
    # Write observations to file
    write_radec_mpc(radec, filename)

    return filename
end 

# Methods to convert a Vector{<:AbstractAstrometry} to a DataFrame
istable(::Type{Vector{<:AbstractAstrometry}}) = true
rowaccess(::Type{Vector{<:AbstractAstrometry}}) = true
rows(x::Vector{<:AbstractAstrometry}) = x
schema(::Vector{T}) where {T <: AbstractAstrometry} = Schema(fieldnames(T), Tuple{fieldtypes(T)...})

# Methods to convert a DataFrame to a Vector{<:AbstractAstrometry}
function Vector{T}(df::DataFrame) where {T <: AbstractAstrometry}
    @assert all(String.(fieldnames(T)) .== names(df)) "`DataFrame` column names don't match `$T` fieldnames"
    @assert all(fieldtypes(T) .== eltype.(eachcol(df))) "`DataFrame` column types don't match `$T` fieldtypes"
    obs = Vector{T}(undef, nrow(df)) 
    for (i, row) in zip(eachindex(obs), eachrow(df))
        obs[i] = T(values(row)...)
    end 
    return obs 
end 

convert(::Type{Vector{T}}, df::DataFrame) where {T <: AbstractAstrometry} = Vector{T}(df)
