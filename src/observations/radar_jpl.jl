@doc raw"""
    RadarJPL{T <: AbstractFloat}

A radar measurement in JPL format.

# Fields

- `id::String`: object's name and number.
- `date::DateTime`: date of observation. 
- `τ::T`: time-delay.
- `τ_σ`: time-delay uncertainty.
- `τ_units::String`: units of time delay. 
- `Δν`: Doppler shift.
- `Δν_σ`: Doppler shift uncertainty.
- `Δν_units`: units of Doppler shift.
- `freq`: frequency of the measurement.
- `rcvr`: ID of reciever antenna.
- `xmit`: ID of emission antenna. 
- `bouncepoint`: bounce point. 

---

    RadarJPL(delay::RegexMatch, doppler::RegexMatch)
    RadarJPL(delay::RegexMatch, doppler::Val{false})
    RadarJPL(delay::Val{false}, doppler::RegexMatch)
    RadarJPL(delay::Val{false}, doppler::Val{false})

Converts a match of `NEOs.jpl_radar_regex` to `RadarJPL`. A `Val{false}` indicates that one or both of the measurements 
(time delay or Doppler shift) are missing. 
"""
struct RadarJPL{T <: AbstractFloat}
    id::String
    date::DateTime
    Δτ::T
    Δτ_σ::T
    Δτ_units::String
    Δν::T
    Δν_σ::T
    Δν_units::String
    freq::T
    rcvr::Int
    xmit::Int
    bouncepoint::String
    # Inner constructor
    function RadarJPL{T}(id::String, date::DateTime, Δτ::T, Δτ_σ::T, Δτ_units::String, Δν::T, Δν_σ::T, Δν_units::String, 
                             freq::T, rcvr::Int, xmit::Int, bouncepoint::String) where {T <: AbstractFloat}
        return new{T}(id, date, Δτ, Δτ_σ, Δτ_units, Δν, Δν_σ, Δν_units, freq, rcvr, xmit, bouncepoint)
    end
end

# Outer constructor
function RadarJPL(id::String, date::DateTime, Δτ::T, Δτ_σ::T, Δτ_units::String, Δν::T, Δν_σ::T, Δν_units::String, freq::T, 
                      rcvr::Int, xmit::Int, bouncepoint::String) where {T <: AbstractFloat}
    return RadarJPL{T}(id, date, Δτ, Δτ_σ, Δτ_units, Δν, Δν_σ, Δν_units, freq, rcvr, xmit, bouncepoint)
end


# Two RadarJPL are equal if ther date, Δτ, Δτ_σ, Δτ_units, Δν, Δν_σ, Δν_units, freq, rcvr, xmit and bouncepoint are equal
function hash(a::RadarJPL{T}, h::UInt) where {T <: AbstractFloat}
    return hash((a.date, a.Δτ, a.Δτ_σ, a.Δτ_units, a.Δν, a.Δν_σ, a.Δν_units, a.freq, a.rcvr, a.xmit, a.bouncepoint), h)
end

function ==(a::RadarJPL{T}, b::RadarJPL{T}) where {T <: AbstractFloat}
    return hash(a) == hash(b)
end

hasdelay(r::RadarJPL{T}) where {T <: AbstractFloat} = !isnan(r.Δτ)
hasdoppler(r::RadarJPL{T}) where {T <: AbstractFloat} = !isnan(r.Δν)

# Print method for RadarJPL
# Examples: 
# 99942 Apophis (2004 MN4) Δτ: 1.9202850713e8 us Δν: -102512.9059 Hz t: 2005-01-29T00:00:00
# 99942 Apophis (2004 MN4) Δν: -100849.1434 Hz t: 2005-01-27T23:31:00
# 99942 Apophis (2004 MN4) Δτ: 9.743930871e7 us t: 2013-01-03T10:00:00
function show(io::IO, r::RadarJPL{T}) where {T <: AbstractFloat} 

    delay = hasdelay(r)
    doppler = hasdoppler(r)

    if delay && doppler
        measurement_s = join([" Δτ: ", string(r.Δτ), " ", r.Δτ_units, " Δν: ", string(r.Δν), " ", r.Δν_units])
    elseif delay
        measurement_s = join([" Δτ: ", string(r.Δτ), " ", r.Δτ_units])
    elseif doppler
        measurement_s = join([" Δν: ", string(r.Δν), " ", r.Δν_units])
    else 
        measurement_s = " No measurements "
    end 

    print(io, r.id, measurement_s, " t: ", r.date)

end

# Functions to get specific fields of a RadarJPL object 
data(r::RadarJPL{T}) where {T <: AbstractFloat} = r.date
delay(r::RadarJPL{T}) where {T <: AbstractFloat} = r.Δτ
delay_sigma(r::RadarJPL{T}) where {T <: AbstractFloat} = r.Δτ_σ
delay_units(r::RadarJPL{T}) where {T <: AbstractFloat} = r.Δτ_units
doppler(r::RadarJPL{T}) where {T <: AbstractFloat} = r.Δν
doppler_sigma(r::RadarJPL{T}) where {T <: AbstractFloat} = r.Δν_units
doppler_units(r::RadarJPL{T}) where {T <: AbstractFloat} = r.Δν_units
freq(r::RadarJPL{T}) where {T <: AbstractFloat} = r.freq
rcvr(r::RadarJPL{T}) where {T <: AbstractFloat} = r.rcvr
xmit(r::RadarJPL{T}) where {T <: AbstractFloat} = r.xmit
bouncepoint(r::RadarJPL{T}) where {T <: AbstractFloat} = r.bouncepoint

@doc raw"""
    ismonostatic(rdata::RadarJPL)

Checks whether the setup is monostatic, i.e., receiver and transmitter are the same. 
"""
ismonostatic(r::RadarJPL{T}) where {T <: AbstractFloat} = r.rcvr == r.xmit

@doc raw"""
    issband(f_MHz::T) where {T<:Real}

Checks whether the transmission frequency `f_MHz` (in MHz) belongs to the S band 
(IEEE nomenclature). 
"""
issband(f_MHz::T) where {T<:Real} = 2000.0 ≤ f_MHz ≤ 4000.0

@doc raw"""
    isxband(f_MHz::T) where {T<:Real}

Checks whether the transmission frequency `f_MHz` (in MHz) belongs to the X band 
(IEEE nomenclature). 
"""
isxband(f_MHz::T) where {T<:Real} = 8000.0 ≤ f_MHz ≤ 12000.0

# Regular expression to parse an optical measurement in MPC format
const jpl_radar_regex = Regex(join(
    [
        # ID regex + tab 
        raw"(?P<id>.*)\t",
        # Date regex + tab
        raw"(?P<date>.*)\t",
        # Measurement regex + tab
        raw"(?P<measurement>.*)\t",
        # Uncertainty regex + tab
        raw"(?P<uncertainty>.*)\t",
        # Units regex + tab
        raw"(?P<units>.*)\t",
        # Frequency regex + tab
        raw"(?P<freq>.*)\t",
        # Reciever regex + tab
        raw"(?P<rcvr>.*)\t",
        # Emitter regex + tab
        raw"(?P<xmit>.*)\t",
        # Bouncepoint regex + end of line 
        raw"(?P<bouncepoint>.*)"
    ]
))
# Format of date in JPL radar data files 
const jpl_radar_dateformat = "yyyy-mm-dd HH:MM:SS"

for (X, Y) in Iterators.product( (:(Val{false}), :(RegexMatch)), (:(Val{false}), :(RegexMatch)) )
    @eval begin 
        function RadarJPL(delay::$X, doppler::$Y)

            if $X == RegexMatch
                id = string(delay["id"])
            elseif $Y == RegexMatch
                id = string(doppler["id"])
            else 
                id = ""
            end 

            if $X == RegexMatch
                date = DateTime(delay["date"], jpl_radar_dateformat)
            elseif $Y == RegexMatch
                date = DateTime(doppler["date"], jpl_radar_dateformat)
            else 
                date = DateTime(2000, 1, 1)
            end 

            if $X == RegexMatch
                Δτ = Meta.parse(delay["measurement"])
                Δτ_σ = Meta.parse(delay["uncertainty"])
                Δτ_units = string(delay["units"])
            else 
                Δτ = NaN
                Δτ_σ = NaN
                Δτ_units = ""
            end  

            if $Y == RegexMatch
                Δν = Meta.parse(doppler["measurement"])
                Δν_σ = Meta.parse(doppler["uncertainty"])
                Δν_units = string(doppler["units"])
            else 
                Δν = NaN
                Δν_σ = NaN
                Δν_units = ""
            end  

            if $X == RegexMatch
                freq = string(delay["freq"])
            elseif $Y == RegexMatch
                freq = string(doppler["freq"])
            else 
                freq = ""
            end 

            if $X == RegexMatch
                freq = Float64(Meta.parse(delay["freq"]))
            elseif $Y == RegexMatch
                freq = Float64(Meta.parse(doppler["freq"]))
            else 
                freq = NaN
            end

            if $X == RegexMatch
                rcvr = Meta.parse(delay["rcvr"])
            elseif $Y == RegexMatch
                rcvr = Meta.parse(doppler["rcvr"])
            else 
                rcvr = 0
            end

            if $X == RegexMatch
                xmit = Meta.parse(delay["xmit"])
            elseif $Y == RegexMatch
                xmit = Meta.parse(doppler["xmit"])
            else 
                xmit = 0
            end

            if $X == RegexMatch
                bouncepoint = string(delay["bouncepoint"])
            elseif $Y == RegexMatch
                bouncepoint = string(doppler["bouncepoint"])
            else 
                bouncepoint = ""
            end
            
            return RadarJPL(id, date, Δτ, Δτ_σ, Δτ_units, Δν, Δν_σ, Δν_units, freq, rcvr, xmit, bouncepoint)
        end
    end 
end 

@doc raw"""
    read_radar_jpl(filename::String)

Returns the matches of `NEOs.jpl_radar_regex` in `filename` as `RadarJPL`.
"""
function read_radar_jpl(filename::String)

    # Read lines of mpc formatted file 
    lines = readlines(filename)
    # Apply regular expressions
    matches = match.(jpl_radar_regex, lines)
    # Eliminate nothings
    filter!(!isnothing, matches)
    # Number of matches 
    m = length(matches)

    # Dates of observation (as string)
    dates = getindex.(matches, "date")
    # Eliminate repeated dates 
    unique!(dates)
    # Number of dates 
    m_dates = length(dates)

    # Time delays 
    delay = filter(x -> x["units"] == "us", matches)
    # Number of time delays 
    m_Δτ = length(delay)
    # Doppler shifts
    doppler = filter(x -> x["units"] == "Hz", matches)
    # Number of Doppler shifts
    m_Δν = length(doppler)

    # Check if there were observations with unknown units 
    if m > m_Δτ + m_Δν
        @warn("""$(m - m_Δτ - m_Δν) observations have units other than  us (time delay) or Hz (doppler shift)""")
    end 

    # Vector of radar observations
    radar = Vector{RadarJPL{Float64}}(undef, m_dates)
    # Time delay index
    i_τ = 1
    # Doppler shift index 
    i_ν = 1

    # Iterate over the dates 
    for i in eachindex(dates)
        # The i_τ-th time delay has the i-th date 
        if (i_τ <= m_Δτ) && (delay[i_τ]["date"] == dates[i])
            Δτ = delay[i_τ]
            i_τ += 1 
            flag_τ = true 
        else 
            flag_τ = false 
        end  

        # The i_ν-th Doppler shift has the i-th date 
        if (i_ν <= m_Δν) && (doppler[i_ν]["date"] == dates[i])
            Δν = doppler[i_ν]
            i_ν += 1 
            flag_ν = true 
        else 
            flag_ν = false 
        end  

        if flag_τ && flag_ν
            radar[i] = RadarJPL(Δτ, Δν)
        elseif flag_τ 
            radar[i] = RadarJPL(Δτ, Val(false))
        elseif flag_ν
            radar[i] = RadarJPL(Val(false), Δν)
        else 
            @warn("""Date $(dates[i]) has no time delay nor Doppler shift measurements""")
            radar[i] = RadarJPL(Val(false), Val(false))
        end

    end 

    return radar 
end 

@doc raw"""
    jpl_radar_str(radar::RadarJPL{T}) where {T <: AbstractFloat}

Returns an observation in JPL format. 
"""
function jpl_radar_str(radar::RadarJPL{T}) where {T <: AbstractFloat}

    if hasdelay(radar)
        delay_s = join([
            radar.id, 
            Dates.format(radar.date, jpl_radar_dateformat), 
            @sprintf("%.2f", radar.Δτ), 
            @sprintf("%1.3f", radar.Δτ_σ),
            radar.Δτ_units,
            @sprintf("%.0f", radar.freq),
            radar.rcvr, 
            radar.xmit, 
            radar.bouncepoint, 
            ""
        ], "\t", "\n")
    else 
        delay_s = ""
    end 

    if hasdoppler(radar)
        doppler_s = join([
            radar.id, 
            Dates.format(radar.date, jpl_radar_dateformat), 
            @sprintf("%.4f", radar.Δν),
            @sprintf("%1.3f", radar.Δν_σ), 
            radar.Δν_units,
            @sprintf("%.0f", radar.freq),
            radar.rcvr, 
            radar.xmit, 
            radar.bouncepoint, 
            ""
        ], "\t", "\n")
    else 
        doppler_s = ""
    end 

    # Join everything
    radar_s = join([delay_s, doppler_s])

    return radar_s
end

@doc raw"""
    write_radar_jpl(radar::Vector{RadarJPL{T}}, filename::String) where {T <: AbstractFloat}

Writes `radar` to `filename` in JPL format. 
"""
function write_radar_jpl(radar::Vector{RadarJPL{T}}, filename::String) where {T <: AbstractFloat}
    open(filename, "w") do file
        for i in eachindex(radar)
            line = jpl_radar_str(radar[i])
            write(file, line)
        end 
    end
end