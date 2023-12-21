@doc raw"""
    TimeOfDay

Day/night at a particular timezone.

# Fields 
- `light::Symbol`: 
    - for ground observatories: `:day` or `:night`,
    - for space observatories: `:space`.
- `start::DateTime`.
- `stop::DateTime`.
- `utc::Int`: hours from UTC.
"""
@auto_hash_equals struct TimeOfDay
    light::Symbol
    start::DateTime
    stop::DateTime
    utc::Int 
    function TimeOfDay(date::DateTime, observatory::ObservatoryMPC{T}) where {T <: AbstractFloat}
        if issatellite(observatory)
            return new(:space, date, date, 0)
        end
        # Hours from UTC
        utc = hours_from_UTC(observatory)
        # Today's sunrise / sunset
        today = sunriseset(date, observatory)
        # Yesterday's sunrise / sunset 
        yesterday = sunriseset(date - Day(1), observatory)
        # Tomorrow's sunrise / sunset
        tomorrow = sunriseset(date + Day(1), observatory)
        # Selection
        if yesterday[2] <= date <= today[1]
            return new(:night, Date(yesterday[2]), Date(today[1]), utc)
        elseif today[1] <= date <= today[2]
            return new(:day, Date(today[1]), Date(today[2]), utc)
        elseif today[2] <= date <= tomorrow[1]
            return new(:night, Date(today[2]), Date(tomorrow[1]), utc)
        end

    end
end

TimeOfDay(radec::RadecMPC{T}) where {T <: AbstractFloat} = TimeOfDay(date(radec), observatory(radec))

isday(x::TimeOfDay) = x.light == :day
isnight(x::TimeOfDay) = x.light == :night

# Print method for TimeOfDay
# Examples:
# Night from 2023-06-29 to 2023-06-29 at UTC-7
# Night from 2023-06-29 to 2023-06-30 at UTC+3
function show(io::IO, m::TimeOfDay)
    print(io, uppercasefirst(string(m.light)), " from ", Date(m.start), " to ", Date(m.stop),
          " at UTC", @sprintf("%+d", m.utc))
end

@doc raw"""
    hours_from_UTC(lon::T) where {T <: Real}

Return the naive hour difference between longitude `lon` [rad] and UTC.
"""
hours_from_UTC(lon::T) where {T <: Real} = ceil(Int, 12*lon/π - 0.5)
function hours_from_UTC(observatory::ObservatoryMPC{T}) where  {T <: AbstractFloat} 
    lon, _ = lonlat(observatory)
    return hours_from_UTC(lon)
end

@doc raw"""
    lonlat(observatory::ObservatoryMPC{T}) where {T <: AbstractFloat}

Return longitude and latitude (both in rad) of an observatory.
"""
function lonlat(observatory::ObservatoryMPC{T}) where {T <: AbstractFloat}
    # ECEF [km]
    p_ECEF = obsposECEF(observatory)
    # ECEF [m] -> Geodetic [m]
    lat_geodetic, lon, altitude = ecef_to_geodetic(1_000 * p_ECEF)
    # Geodetic [m] -> Geocentric [m]
    lat_geocentric, _ = geodetic_to_geocentric(lat_geodetic, altitude)

    return lon, lat_geocentric
end 

@doc raw"""
    sunriseset(date::DateTime, observatory::ObservatoryMPC{T}) where {T <: AbstractFloat}
    sunriseset(radec::RadecMPC{T}) where {T <: AbstractFloat}

Return `DateTime` of sunrise and sunset at a particular date and location.

!!! reference
    See "General Solar Position Calculations" by NOAA at
    https://gml.noaa.gov/grad/solcalc/solareqns.PDF.
"""
function sunriseset(date::DateTime, observatory::ObservatoryMPC{T}) where {T <: AbstractFloat}
    # Fractional year [rad]
    γ = 2π * (dayofyear(date) - 1 + (hour(date)-12)/24) / daysinyear(date)
    # Equation of time [min]
    eqtime = 229.18 * (0.000075 + 0.001868*cos(γ) - 0.032077*sin(γ) - 0.014615*cos(2γ)
             - 0.040849*sin(2γ) )
    # Solar declination angle [rad]
    δ = 0.006918 - 0.399912*cos(γ) + 0.070257*sin(γ) - 0.006758*cos(2γ) + 0.000907*sin(2γ)
        - 0.002697*cos(3γ) + 0.00148*sin(3γ)
    # Longitude and latitude [rad]
    lon, lat = lonlat(observatory)
    # Solar hour angle [deg]
    ha_sunrise = acosd(cosd(90.833)/cos(lat)/cos(δ) - tan(lat)*tan(δ))
    ha_sunset = -ha_sunrise
    # Day [DateTime]
    day = DateTime(Date(date))
    # UTC time of sunrise [min]
    sunrise_min = 720 - 4*(rad2deg(lon) + ha_sunrise) - eqtime
    # UTC time of sunrise [DateTime]
    sunrise = day + Microsecond(round(Int, sunrise_min*6e7))
    # UTC time of sunset [min]
    sunset_min = 720 - 4*(rad2deg(lon) + ha_sunset) - eqtime 
    # UTC time of sunset [DateTime]
    sunset = day + Microsecond(round(Int, sunset_min*6e7))

    return sunrise, sunset
end

sunriseset(radec::RadecMPC{T}) where {T <: AbstractFloat} = sunriseset(date(radec), observatory(radec))

@doc raw"""
    obsposECEF(observatory::ObservatoryMPC{T}; kwarg) where {T <: AbstractFloat}
    obsposECEF(x::RadecMPC{T}; kwarg) where {T <: AbstractFloat}
    obsposECEF(x::RadarJPL{T}; kwarg) where {T <: AbstractFloat}

Return the observer's geocentric `[x, y, z]` position vector in Earth-Centered Earth-Fixed
(ECEF) reference frame.

# Keyword argument
- `eop::Union{EopIau1980, EopIau2000A}`: Earth Orientation Parameters (eop).
"""
function obsposECEF(observatory::ObservatoryMPC{T}; eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat}

    # Make sure observatory has coordinates
    @assert hascoord(observatory) "Cannot compute position for observatory [$(observatory.code)] without coordinates"

    if issatellite(observatory) || isoccultation(observatory)
        # Ephemeris seconds since J2000
        et = datetime2et(observatory.date)
        # Earth-Centered Inertial position position of observer
        posvel_ECI = obsposvelECI(observatory, et; eop)
        # UTC seconds
        utc_secs = et - tdb_utc(et)
        # Julian days UTC
        jd_utc = JD_J2000 + utc_secs/daysec
        # State vector
        pv_ECI = OrbitStateVector(jd_utc, posvel_ECI[1:3], posvel_ECI[4:6], zeros(3))

        # Transform position/velocity from Earth-Centered Inertial (ECI) fraom to Earth-Centered Earth-fixed (ECEF) frame
        # ITRF: International Terrestrial Reference Frame
        # GCRF: Geocentric Celestial Reference Frame
        pv_ECEF = sv_eci_to_ecef(pv_ECI, Val(:GCRF), Val(:ITRF), jd_utc, eop)

        # ECEF position
        pos_ECEF = convert(Vector{eltype(pv_ECEF.r)}, pv_ECEF.r)
        
        return pos_ECEF
    else 
        # λ_deg: longitude [degrees east of Greenwich]
        # u: distance from spin axis [km], u = ρ*cos(ϕ')
        # v: height above equatorial plane [km], v = ρ*sin(ϕ'),
        # where ϕ' is the geocentric latitude and ρ is the geocentric distance in km

        # Cilindrical components of Earth-Centered Earth-Fixed position of observer
        λ_deg = observatory.long     # deg
        u = observatory.cos * RE     # km
        v = observatory.sin * RE     # km

        # Cartesian components of Earth-Centered Earth-Fixed position of observer
        λ_rad = deg2rad(λ_deg)       # rad
        x_gc = u * cos(λ_rad)        # km
        y_gc = u * sin(λ_rad)        # km
        z_gc = v                     # km

        # Earth-Centered Earth-Fixed position position of observer
        pos_ECEF = [x_gc, y_gc, z_gc] # km

        return pos_ECEF
    end
end

obsposECEF(x::RadecMPC{T}; eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat} = obsposECEF(x.observatory; eop)
obsposECEF(x::RadarJPL{T}; eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat} = obsposECEF(x.rcvr; eop)

# TODO: avoid sv_ecef_to_ecef overload by defining proper product between DCMs and Taylor1/TaylorN
# method below has been adapted from SatelliteToolboxTransformations.jl, MIT-licensed
#   https://github.com/JuliaSpace/SatelliteToolboxTransformations.jl
function sv_ecef_to_ecef(
    sv::OrbitStateVector,
    T_ECEF1::Val{:ITRF},
    T_ECEF2::Val{:TIRS},
    jd_utc::Taylor1{TaylorN{Float64}},
    eop_data::Union{Nothing, EopIau1980, EopIau2000A} = nothing
)
    D = r_ecef_to_ecef(DCM, T_ECEF1, T_ECEF2, jd_utc, eop_data)

    # Since both frames does not have a significant angular velocity between
    # them, then we just need to convert the representations.
    r_ecef::Vector{Taylor1{TaylorN{Float64}}} = D * sv.r
    v_ecef::Vector{Taylor1{TaylorN{Float64}}} = D * sv.v
    a_ecef::Vector{Taylor1{TaylorN{Float64}}} = D * sv.a
    return OrbitStateVector(sv.t, r_ecef, v_ecef, a_ecef)
end

# TODO: avoid sv_ecef_to_eci overload by defining proper product between DCMs and Taylor1/TaylorN
# method below has been adapted from SatelliteToolboxTransformations.jl, MIT-licensed
#   https://github.com/JuliaSpace/SatelliteToolboxTransformations.jl
function sv_ecef_to_eci(
    sv::OrbitStateVector,
    T_ECEF::Union{Val{:PEF}, Val{:TIRS}},
    T_ECI::Union{T_ECIs, T_ECIs_IAU_2006},
    jd_utc::Taylor1{TaylorN{Float64}},
    eop_data::Union{Nothing, EopIau1980, EopIau2000A} = nothing
)
    # Get the matrix that converts the ECEF to the ECI.
    if eop_data === nothing
        D = r_ecef_to_eci(DCM, T_ECEF, T_ECI, jd_utc)
    else
        D = r_ecef_to_eci(DCM, T_ECEF, T_ECI, jd_utc, eop_data)
    end

    # Since the ECI and ECEF frames have a relative velocity between them, then
    # we must account from it when converting the velocity and acceleration. The
    # angular velocity between those frames is computed using `we` and corrected
    # by the length of day (LOD) parameter of the EOP data, if available.
    ω  = EARTH_ANGULAR_SPEED * (1 - (eop_data !== nothing ? eop_data.lod(jd_utc) / 86400000 : 0))
    vω = [0, 0, ω]

    # Compute the position in the ECI frame.
    r_eci::Vector{Taylor1{TaylorN{Float64}}} = D * sv.r

    # Compute the velocity in the ECI frame.
    vω_x_r = vω × sv.r
    v_eci::Vector{Taylor1{TaylorN{Float64}}} = D * (sv.v + vω_x_r )

    # Compute the acceleration in the ECI frame.
    a_eci::Vector{Taylor1{TaylorN{Float64}}} = D * (sv.a + vω × vω_x_r + 2vω × sv.v)

    return OrbitStateVector(sv.t, r_eci, v_eci, a_eci)
end

@doc raw"""
    obsposvelECI(observatory::ObservatoryMPC{T}, et::U; kwarg) where {T <: AbstractFloat, U <: Number}
    obsposvelECI(x::RadecMPC{T}; kwarg) where {T <: AbstractFloat}
    obsposvelECI(x::RadarJPL{T}; kwarg) where {T <: AbstractFloat}

Return the observer's geocentric `[x, y, z, v_x, v_y, v_z]` "state" vector in Earth-Centered
Inertial (ECI) reference frame. By default, the IAU200A Earth orientation model is used to
transform from Earth-centered, Earth-fixed (ECEF) frame to ECI frame. Other Earth orientation
models, such as the IAU1976/80 model, can be used by importing the
`SatelliteToolboxTransformations.EopIau1980` type and passing it to the `eop` keyword
argument in the function call.

See also [`SatelliteToolboxBase.OrbitStateVector`](@ref) and 
[`SatelliteToolboxTransformations.sv_ecef_to_eci`](@ref).

# Arguments

- `observatory::ObservatoryMPC{T}`: observation site.
- `et::U`: ephemeris time (TDB seconds since J2000.0 epoch).
- `x::RadecMPC{T}/RadarJPL{T}`: astrometric observation.

# Keyword argument

- `eop::Union{EopIau1980, EopIau2000A}`: Earth Orientation Parameters (eop).
"""
function obsposvelECI(observatory::ObservatoryMPC{T}, et::U;
        eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat, U <: Number}
    
    # Make sure observatory has coordinates
    @assert hascoord(observatory) "Cannot compute position for observatory [$(observatory.code)] without coordinates"
    # One with correct type
    oneU = one(et)

    if issatellite(observatory) || isoccultation(observatory)
        # @assert datetime2et(observatory.date) == cte(et)
        return [observatory.long * oneU, observatory.cos * oneU, observatory.sin * oneU,
                zero(et), zero(et), zero(et)]
    else 
        # Earth-Centered Earth-Fixed position position of observer
        pos_ECEF = obsposECEF(observatory)

        # UTC seconds
        utc_secs = et - tdb_utc(et)
        # Julian days UTC
        jd_utc = JD_J2000 + utc_secs/daysec
        # State vector
        pv_ECEF = OrbitStateVector(jd_utc, pos_ECEF, zeros(3), zeros(3))

        # Transform position/velocity from Earth-Centered Earth-fixed (ECEF) frame to Earth-Centered Inertial (ECI) frame
        # ITRF: International Terrestrial Reference Frame
        # GCRF: Geocentric Celestial Reference Frame
        pv_ECI = sv_ecef_to_eci(pv_ECEF, Val(:ITRF), Val(:GCRF), jd_utc, eop)

        # ECI state vector (of correct type)
        r_ECI = Vector{U}(undef, 6)
        # Inertial position
        for i in eachindex(pv_ECI.r)
            r_ECI[i] = pv_ECI.r[i] * oneU
        end
        # Inertial velocity
        for i in eachindex(pv_ECI.v)
            r_ECI[i+3] = pv_ECI.v[i] * oneU
        end
        
        return r_ECI
    end
end

function obsposvelECI(x::RadecMPC{T}; eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat}
    return obsposvelECI(x.observatory, datetime2et(x.date); eop)
end
function obsposvelECI(x::RadarJPL{T}; eop::Union{EopIau1980, EopIau2000A} = eop_IAU2000A) where {T <: AbstractFloat}
    return obsposvelECI(x.rcvr, datetime2et(x.date); eop)
end
